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
subroutine Stencils_lsCenterPropsSharp(phi,iPropC,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2)

    implicit none

    !------Argument List-----
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropC

    !------Local variables-----
    real, dimension(ix1:ix2,jy1:jy2,kz1:kz2) :: pf

    pf(ix1:ix2,jy1:jy2,kz1:kz2)  = (sign(1.0,phi(ix1:ix2,jy1:jy2,kz1:kz2))+1.0)/2.0

    iPropC(ix1:ix2,jy1:jy2,kz1:kz2) = pf(ix1:ix2,jy1:jy2,kz1:kz2)*lsScalarProp + &
                                     (1.0-pf(ix1:ix2,jy1:jy2,kz1:kz2))*iPropC(ix1:ix2,jy1:jy2,kz1:kz2)

end subroutine Stencils_lsCenterPropsSharp

subroutine Stencils_lsCenterPropsSmeared(phi,iPropC,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2,iSmear)

    implicit none

    !------Argument List-----
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropC
    real, intent(in) :: iSmear

    !------Local variables-----
    integer :: i,j,k
    real, dimension(ix1:ix2,jy1:jy2,kz1:kz2) :: smhv
    real, parameter :: pi=acos(-1.0)

    do k=kz1,kz2
     do j=jy1,jy2
       do i=ix1,ix2

        if(abs(phi(i,j,k)) .le. iSmear) then ! Symmetric smearing - A.Dhruv
        !if(abs(phi(i,j,k)) .le. iSmear .and. phi(i,j,k) .lt. 0.0) then ! Asymmetric smearing - A.Dhruv
              smhv(i,j,k) = 0.5 + phi(i,j,k)/(2*iSmear) + sin(2*pi*phi(i,j,k)/(2*iSmear))/(2*pi)

        else
          if(phi(i,j,k) .ge. 0.0) then
             smhv(i,j,k) = 1.0

          else
             smhv(i,j,k) = 0.0

          end if
        end if

      end do
     end do
    end do

    iPropC(ix1:ix2,jy1:jy2,kz1:kz2) = smhv(ix1:ix2,jy1:jy2,kz1:kz2)*lsScalarProp + &
                                     (1.0-smhv(ix1:ix2,jy1:jy2,kz1:kz2))*iPropC(ix1:ix2,jy1:jy2,kz1:kz2)

end subroutine Stencils_lsCenterPropsSmeared
