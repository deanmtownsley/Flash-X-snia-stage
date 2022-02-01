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
subroutine Stencils_lsFaceProps3dFixedSmear(phi,iPropX,iPropY,iPropZ,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2)

    implicit none

    !------Argument List-----
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY,iPropZ

    !------Local variables-----
    integer :: i,j,k
    real :: a1,a2,aa,th
    real, parameter :: eps = 1E-13
    real, dimension(ix1:ix2,jy1:jy2,kz1:kz2) :: pf
    real :: sp, pi

    pi = acos(-1.0)

    pf(ix1:ix2,jy1:jy2,kz1:kz2)   = (sign(1.0,phi(ix1:ix2,jy1:jy2,kz1:kz2))+1.0)/2.0

    !--Property on x-face
    !---Loop through boundary and interior cell faces
    do k = kz1,kz2-1
     do j = jy1,jy2-1
      do i = ix1,ix2-1

         a1 = (pf(i+1,j,k) + pf(i,j,k)) / 2.                       
         a2 = pf(i+1,j,k)  /abs(pf(i+1,j,k)  +eps) * &
              pf(i,j,k)/abs(pf(i,j,k)+eps)

         if(pf(i,j,k).eq.0..and.pf(i+1,j,k).eq.1.) then
         
            th = abs(phi(i+1,j,k))/(abs(phi(i+1,j,k))+abs(phi(i,j,k)))
            aa = th*(1./lsScalarProp) + (1.-th)*(1./iPropX(i+1,j,k))
            iPropX(i+1,j,k) = a1*a2/aa + (1. - a1*a2)/aa
         
         else if(pf(i,j,k).eq.1..and.pf(i+1,j,k).eq.0.) then 
         
            th = abs(phi(i,j,k))/(abs(phi(i,j,k))+abs(phi(i+1,j,k)))
            aa = th*(1./lsScalarProp) + (1.-th)*(1./iPropX(i+1,j,k))
            iPropX(i+1,j,k) = a1*a2/aa + (1. - a1*a2)/aa
         
         else
         
            iPropX(i+1,j,k) = a1*a2*lsScalarProp + (1. - a1*a2)*iPropX(i+1,j,k)            

         end if
       end do
      end do
     end do

     !----Property on y-face
     !----Loop through boundary and interior cell faces
     do k = kz1,kz2-1
      do i = ix1,ix2-1
       do j = jy1,jy2-1

          a1 = (pf(i,j+1,k) + pf(i,j,k)) / 2.           
          a2 = pf(i,j+1,k)  /abs(pf(i,j+1,k)  +eps) * &
               pf(i,j,k)/abs(pf(i,j,k)+eps)

          if(pf(i,j,k).eq.0..and.pf(i,j+1,k).eq.1.) then
          
            th = abs(phi(i,j+1,k))/(abs(phi(i,j+1,k))+abs(phi(i,j,k)))
            aa = th*(1./lsScalarProp) + (1.-th)*(1./iPropY(i,j+1,k))
            iPropY(i,j+1,k) = a1*a2/aa + (1. - a1*a2)/aa

          else if(pf(i,j,k).eq.1..and.pf(i,j+1,k).eq.0.) then 
          
            th = abs(phi(i,j,k))/(abs(phi(i,j,k))+abs(phi(i,j+1,k)))
            aa = th*(1./lsScalarProp) + (1.-th)*(1./iPropY(i,j+1,k))
            iPropY(i,j+1,k) = a1*a2/aa + (1. - a1*a2)/aa

          else
    
            iPropY(i,j+1,k) = a1*a2*lsScalarProp + (1. - a1*a2)*iPropY(i,j+1,k)    

         end if
       end do
      end do
     end do

     !----Property on z-face
     !----Loop through boundary and interior cell faces
     do k = kz1,kz2-1
      do i = ix1,ix2-1
       do j = jy1,jy2-1

          a1 = (pf(i,j,k+1) + pf(i,j,k)) / 2.           
          a2 = pf(i,j,k+1)  /abs(pf(i,j,k+1)  +eps) * &
               pf(i,j,k)/abs(pf(i,j,k)+eps)

          if(pf(i,j,k).eq.0..and.pf(i,j,k+1).eq.1.) then
          
            th = abs(phi(i,j,k+1))/(abs(phi(i,j,k+1))+abs(phi(i,j,k)))
            aa = th*(1./lsScalarProp) + (1.-th)*(1./iPropZ(i,j,k+1))
            iPropZ(i,j,k+1) = a1*a2/aa + (1. - a1*a2)/aa

          else if(pf(i,j,k).eq.1..and.pf(i,j,k+1).eq.0.) then 
          
            th = abs(phi(i,j,k))/(abs(phi(i,j,k))+abs(phi(i,j,k+1)))
            aa = th*(1./lsScalarProp) + (1.-th)*(1./iPropZ(i,j,k+1))
            iPropZ(i,j,k+1) = a1*a2/aa + (1. - a1*a2)/aa

          else
    
            iPropZ(i,j,k+1) = a1*a2*lsScalarProp + (1. - a1*a2)*iPropZ(i,j,k+1)    

         end if
       end do
      end do
     end do


end subroutine Stencils_lsFaceProps3dFixedSmear

subroutine Stencils_lsFaceProps3dUserSmear(phi,iPropX,iPropY,iPropZ,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2,iSmear)

    implicit none

    !------Argument List-----
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY,iPropZ
    real, intent(in) :: iSmear

    !------Local variables-----
    integer :: i,j,k
    real :: a1,a2,aa,th
    real, parameter :: eps = 1E-13
    real, dimension(ix1:ix2,jy1:jy2,kz1:kz2) :: pf
    real :: sp, pi, phiF, smhv

    pi = acos(-1.0)

    pf(ix1:ix2,jy1:jy2,kz1:kz2)   = (sign(1.0,phi(ix1:ix2,jy1:jy2,kz1:kz2))+1.0)/2.0

    !--Property on x-face
    !---Loop through boundary and interior cell faces
    do k = kz1,kz2-1
     do j = jy1,jy2-1
      do i = ix1,ix2-1

         a1 = (pf(i+1,j,k) + pf(i,j,k)) / 2.                       
         a2 = pf(i+1,j,k)  /abs(pf(i+1,j,k)  +eps) * &
              pf(i,j,k)/abs(pf(i,j,k)+eps)

         phiF = 0.5*(phi(i,j,k)+phi(i+1,j,k))

         if(abs(phiF) .le. iSmear) then ! Symmetric smearing - A.Dhruv
         !if(abs(phiF) .le. iSmear .and. phiF .lt. 0.0) then ! Asymmetric smearing - A.Dhruv
            smhv = 0.5 + phiF/(2*iSmear) + sin(2*pi*phiF/(2*iSmear))/(2*pi)
            iPropX(i+1,j,k) = smhv*lsScalarProp + (1-smhv)*iPropX(i+1,j,k)                 
         else         
            iPropX(i+1,j,k) = a1*a2*lsScalarProp + (1. - a1*a2)*iPropX(i+1,j,k)            

         end if

       end do
      end do
     end do

     !----Property on y-face
     !----Loop through boundary and interior cell faces
     do k = kz1,kz2-1
      do i = ix1,ix2-1
       do j = jy1,jy2-1

          a1 = (pf(i,j+1,k) + pf(i,j,k)) / 2.           
          a2 = pf(i,j+1,k)  /abs(pf(i,j+1,k)  +eps) * &
               pf(i,j,k)/abs(pf(i,j,k)+eps)

          phiF = 0.5*(phi(i,j,k)+phi(i,j+1,k))

          if(abs(phiF) .le. iSmear) then ! Symmetric smearing - A.Dhruv
          !if(abs(phiF) .le. iSmear .and. phiF .lt. 0.0) then ! Asymmetric smearing - A.Dhruv
            smhv = 0.5 + phiF/(2*iSmear) + sin(2*pi*phiF/(2*iSmear))/(2*pi)
            iPropY(i,j+1,k) = smhv*lsScalarProp + (1-smhv)*iPropY(i,j+1,k)
          else    
            iPropY(i,j+1,k) = a1*a2*lsScalarProp + (1. - a1*a2)*iPropY(i,j+1,k)    

          end if

       end do
      end do
     end do

     !----Property on z-face
     !----Loop through boundary and interior cell faces
     do k = kz1,kz2-1
      do i = ix1,ix2-1
       do j = jy1,jy2-1

          a1 = (pf(i,j,k+1) + pf(i,j,k)) / 2.           
          a2 = pf(i,j,k+1)  /abs(pf(i,j,k+1)  +eps) * &
               pf(i,j,k)/abs(pf(i,j,k)+eps)
 
          phiF = 0.5*(phi(i,j,k)+phi(i,j,k+1))

          if(abs(phiF) .le. iSmear) then ! Symmetric smearing - A.Dhruv
          !if(abs(phiF) .le. iSmear .and. phiF .lt. 0.0) then ! Asymmetric smearing - A.Dhruv
            smhv = 0.5 + phiF/(2*iSmear) + sin(2*pi*phiF/(2*iSmear))/(2*pi)
            iPropZ(i,j,k+1) = smhv*lsScalarProp + (1-smhv)*iPropZ(i,j,k+1)
 
          else    
            iPropZ(i,j,k+1) = a1*a2*lsScalarProp + (1. - a1*a2)*iPropZ(i,j,k+1)    

          end if

       end do
      end do
     end do

end subroutine Stencils_lsFaceProps3dUserSmear
