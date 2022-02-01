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
!!**
#include "Simulation.h"
#include "constants.h"

subroutine ins_setupPoissonRhs_constdens(divu,dt)
        implicit none
        real, dimension(:,:,:), intent(inout) :: divu
        real, intent(in) :: dt
        divu = divu/dt                                               
end subroutine ins_setupPoissonRhs_constdens

subroutine ins_setupPoissonRhs_vardens(divu,&
                                  sigx,sigy,sigz,&
                                  pxn1,pyn1,pzn1,&
                                  pxn2,pyn2,pzn2,&
                                  rhox,rhoy,rhoz,&
                                  rhoGas,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)

        implicit none
        real, dimension(:,:,:), intent(inout) :: divu
        real, dimension(:,:,:), intent(in) :: sigx,sigy,sigz
        real, dimension(:,:,:), intent(in) :: rhox,rhoy,rhoz,pxn1,pyn1,pzn1,&
                                              pxn2,pyn2,pzn2        
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dt,dx,dy,dz,rhoGas
                                                
        integer :: i,j,k
        real :: d1x,d1y,p1x,p1y1,p1y2,p1y,&
                d2x,d2y,p2x,p2y1,p2y2,p2y,&
                d1z,d2z,p1z,p2z1,p2z2,p2z,&
                term1,term2,term3        

        do k=kz1,kz2
         do i=ix1,ix2
           do j=jy1,jy2

              term1 = (sigx(i+1,j,k)-sigx(i,j,k))/dx + (sigy(i,j+1,k)-sigy(i,j,k))/dy
    
              d1x = 1.0-rhoGas*rhox(i+1,j,k)
              p1x = 2*pxn1(i+1,j,k)-pxn2(i+1,j,k)
              d2x = 1.0-rhoGas*rhox(i,j,k)
              p2x = 2*pxn1(i,j,k)-pxn2(i,j,k)

              d1y = 1.0-rhoGas*rhoy(i,j+1,k)
              p1y = 2*pyn1(i,j+1,k)-pyn2(i,j+1,k)
              d2y = 1.0-rhoGas*rhoy(i,j,k)
              p2y = 2*pyn1(i,j,k)-pyn2(i,j,k)

              term2 = (d1x*p1x-d2x*p2x)/dx + (d1y*p1y-d2y*p2y)/dy

#if NDIM == 3
              term1 = term1 + (sigz(i,j,k+1)-sigz(i,j,k))/dz

              d1z = 1.0-rhoGas*rhoz(i,j,k+1)
              p1z = 2*pzn1(i,j,k+1)-pzn2(i,j,k+1)
              d2z = 1.0-rhoGas*rhoz(i,j,k)
              p2z = 2*pzn1(i,j,k)-pzn2(i,j,k)

              term2 = term2 + (d1z*p1z-d2z*p2z)/dz
#endif

              term3 = (rhoGas/dt)*divu(i,j,k)

              divu(i,j,k) = term1+term2+term3
 
          enddo
         enddo
        enddo

end subroutine ins_setupPoissonRhs_vardens

