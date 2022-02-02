
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!
!!
!!***
#include "Simulation.h"
subroutine Stencils_lsRedistance3d(phi,phi_orig,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)

        implicit none
        !-----Imported variables
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, dimension(:,:,:), intent(inout):: phi
        real, intent(in) :: dx,dy,dz,dt
        real, dimension(:,:,:), intent(in):: phi_orig

        !------Local variables
        real ::    phio(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
                   sgn(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD), &
                   eps,t,agf,dtL,err1,err2, &
                   ap,an,bp,bn,cp,cn,dp,dn, &
                   sxl,sxr,syl,syr,sm, &
                   ep,en,fp,fn,szl,szr,Dij   

        integer :: i,j,k,n,incrm,m

        eps = 1E-14
        incrm = 1
        m = 1
        err1 = 0.

        !--Set the sign using the ORIGINAL distance function
        do k = kz1,kz2
           do j = jy1,jy2
              do i = ix1,ix2
                 sgn(i,j,k) = phi_orig(i,j,k)/abs(phi_orig(i,j,k)+eps)
              end do
           end do
        end do

        phio = phi

        do k = kz1,kz2
           do j = jy1,jy2
              do i = ix1,ix2

                 phi(i,j,k) = 0.
 
                 !- kpd - If cell is near interface...
                 if ( (phio(i,j,k)*phio(i-1,j,k) .LT. 0) .OR. &
                      (phio(i,j,k)*phio(i+1,j,k) .LT. 0) .OR. &
                      (phio(i,j,k)*phio(i,j-1,k) .LT. 0) .OR. &
                      (phio(i,j,k)*phio(i,j+1,k) .LT. 0) .OR. &
                      (phio(i,j,k)*phio(i,j,k-1) .LT. 0) .OR. &
                      (phio(i,j,k)*phio(i,j,k+1) .LT. 0) ) then


                    !------------------------------------------------------
                    !---Solve Level Set Re-distance equation ---------
                    !------------------------------------------------------
                    phi(i,j,k) = phi_orig(i,j,k) 
                    !------------------------------------------------------


                 !- kpd - If cell is NOT near interface...
                 else

                   !- kpd - Setup distance function stencil.
                   sxl = phio(i-1,j,k)
                   sxr = phio(i+1,j,k)
                   syl = phio(i,j-1,k)
                   syr = phio(i,j+1,k)
                   szl = phio(i,j,k-1)
                   szr = phio(i,j,k+1)

                   sm = phio(i,j,k)

                   ap = max((sm-sxl),0.)/dx
                   an = min((sm-sxl),0.)/dx

                   bp = max((sxr-sm),0.)/dx
                   bn = min((sxr-sm),0.)/dx

                   cp = max((sm-syl),0.)/dy
                   cn = min((sm-syl),0.)/dy

                   dp = max((syr-sm),0.)/dy
                   dn = min((syr-sm),0.)/dy

                   ep = max((sm-szl),0.)/dz
                   en = min((sm-szl),0.)/dz

                   fp = max((szr-sm),0.)/dz
                   fn = min((szr-sm),0.)/dz

                   !---------------------------------------------
                   !- kpd - Compute the magnitude of the gradient
                   !---------------------------------------------
                   if(phi_orig(i,j,k).gt.0.) then
                      agf = sqrt( max(ap**2,bn**2) &
                                + max(cp**2,dn**2) &
                                + max(ep**2,fn**2) ) - 1.0
                   elseif(phi_orig(i,j,k).lt.0.) then
                      agf = sqrt( max(an**2,bp**2) &
                                + max(cn**2,dp**2) &
                                + max(en**2,fp**2) ) - 1.0
                   else
                      agf = 0.
                   end if

                   !------------------------------------------------------
                   !- kpd - Solve Level Set Re-distance equation ---------
                   !------------------------------------------------------
                   phi(i,j,k) = phio(i,j,k) - dt*(sgn(i,j,k)*(agf))
                   !------------------------------------------------------

                 end if

                 if (SIGN(1.0,phi(i,j,k)) .NE. SIGN(1.0,phi_orig(i,j,k))) then
                    print*,"WARNING: LS Dist Function Changed Signs - ",i,j,k
                 end if

                 end do    !do i = ix1,ix2
              end do       !do j = jx1,jx2
           end do          !do k = kx1,kx2

end subroutine Stencils_lsRedistance3d
