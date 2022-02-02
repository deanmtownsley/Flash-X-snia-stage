!!
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
!!***
#include "Simulation.h"
subroutine Stencils_lsRedistance2d(phi,phi_orig,dt,dx,dy,ix1,ix2,jy1,jy2)

        implicit none
        !-----Imported variables
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, dimension(:,:,:), intent(inout):: phi
        real, intent(in) :: dx,dy,dt
        real, dimension(:,:,:), intent(in):: phi_orig

        !------Local variables
        real ::    phio(NXB+2*NGUARD,NYB+2*NGUARD,1), &
                   sgn(NXB+2*NGUARD,NYB+2*NGUARD,1), &
                   eps,t,agf,dtL,err1,err2, &
                   ap,an,bp,bn,cp,cn,dp,dn, &
                   sxl,sxr,syl,syr,sm,Dij   

        real :: phiXP,phiXM,phiYP,phiYM,xP,xM,yP,yM,eX
        real :: dd, dx1, dx2, dx3, dy1, dy2, dy3
        integer :: i,j,k,n,m,itr

        !---For 2-D simulations
        k=1

        eps = 1E-15

        !- kpd- Set the sign using the ORIGINAL distance function
        do j = jy1,jy2
           do i = ix1,ix2
              sgn(i,j,k) = phi_orig(i,j,k)/abs(phi_orig(i,j,k)+eps)
           end do
        end do

        phio  = phi

        !-Loop through interior nodes
        do j = jy1,jy2
           do i = ix1,ix2

            phi(i,j,k) = 0.
 
            sm  = phio(i,j,k)

            sxl = phio(i-1,j,k)
            sxr = phio(i+1,j,k)
            syl = phio(i,j-1,k)
            syr = phio(i,j+1,k)

            if ( (phio(i,j,k)*phio(i-1,j,k) .LT. 0) .OR. &
                 (phio(i,j,k)*phio(i+1,j,k) .LT. 0) .OR. &
                 (phio(i,j,k)*phio(i,j-1,k) .LT. 0) .OR. &
                 (phio(i,j,k)*phio(i,j+1,k) .LT. 0) ) then
               !------------------------------------------------------
               !------------------------------------------------------
               phi(i,j,k) = phi_orig(i,j,k) 
               !------------------------------------------------------

              !------------------------------------------------------
            else

               !-----First Order Upwind... +++++++++++++++++++++++
               ap = max((sm-sxl),0.)/dx
               an = min((sm-sxl),0.)/dx

               bp = max((sxr-sm),0.)/dx
               bn = min((sxr-sm),0.)/dx

               cp = max((sm-syl),0.)/dy
               cn = min((sm-syl),0.)/dy

               dp = max((syr-sm),0.)/dy
               dn = min((syr-sm),0.)/dy
               !+++++++++++++++++++++++++++++++++++++++++++++++++++++

               if(phio(i,j,k).gt.0.) then
                  agf = ( sqrt(max(ap**2,bn**2) + max(cp**2,dn**2)) ) -1.0
               elseif(phio(i,j,k).lt.0.) then
                  agf = ( sqrt(max(an**2,bp**2) + max(cn**2,dp**2)) ) -1.0
               else
                  agf = 0.
               end if
              
               !------------------------------------------------------
               !-------Solve Level Set Re-distance equation ---------
               !------------------------------------------------------
               phi(i,j,k) = phio(i,j,k) - dt*(sgn(i,j,k)*(agf))
               !------------------------------------------------------

            end if  !if interface within one cell

            if (SIGN(1.0,phi(i,j,k)) .NE. SIGN(1.0,phi_orig(i,j,k))) then
               print*,"[Stencils_lsRedistance2d]: WARNING: LS Dist Function Changed Signs - ",i,j,k
            end if

           end do !do j = jy1,jy2
        end do    !do i = iy1,iy2

        return
end subroutine Stencils_lsRedistance2d 
