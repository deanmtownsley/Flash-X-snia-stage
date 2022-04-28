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
subroutine Stencils_lsNormals2d(phi,iNrmx,iNrmy,dx,dy,ix1,ix2,jy1,jy2)

    implicit none
    !------Argument List-----
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iNrmx, iNrmy
    real, intent(in) :: dx,dy

    !-----Local variables---
    integer :: i,j,k
    real :: gradX,gradY
    real :: tol = 1e-13

    k=1

    do j=jy1+1,jy2-1
      do i=ix1+1,ix2-1
        gradX = (phi(i+1,j,k)-phi(i-1,j,k))/(2.*dx)
        gradY = (phi(i,j+1,k)-phi(i,j-1,k))/(2.*dy)

        iNrmx(i,j,k) = gradX/sqrt(gradX**2+gradY**2+tol)
        iNrmy(i,j,k) = gradY/sqrt(gradX**2+gradY**2+tol)
      end do
    end do

end subroutine Stencils_lsNormals2d
