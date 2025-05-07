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
#include "constants.h"
#include "Simulation.h"
subroutine Stencils_diffusionConst(rhs, phi, del, Coeff, lo, hi)

  implicit none
  !---Argument List -----
  real, dimension(:,:,:), intent(inout) :: rhs
  real, dimension(:,:,:), intent(in)  :: phi
  real, dimension(MDIM), intent(in) :: del
  real, intent(in) :: Coeff
  integer,dimension(MDIM), intent(in) :: lo, hi

  !---Local Variables
  integer :: i,j,k

  do k=lo(KAXIS),hi(KAXIS)
     do j=lo(JAXIS),hi(JAXIS)
        do i=lo(IAXIS),hi(IAXIS)
           rhs(i,j,k) = rhs(i,j,k) + (Coeff/(del(IAXIS)**2))*(phi(i+1,j,k)+phi(i-1,j,k)-2.*phi(i,j,k))&
                + (Coeff/(del(JAXIS)**2))*(phi(i,j+1,k)+phi(i,j-1,k)-2.*phi(i,j,k))
#if NDIM == 3           
           rhs(i,j,k)=rhs(i,j,k)+(Coeff/(del(KAXIS)**2))*(phi(i,j,k+1)+phi(i,j,k-1)-2.*phi(i,j,k))
#endif           
        end do
     end do
  end do
  return

end subroutine Stencils_diffusionConst

subroutine Stencils_diffusionVar(rhs, phi, del, Coeff, lo, hi)


  implicit none
  !---Argument List -----
  real, dimension(:,:,:), intent(inout) :: rhs
  real, dimension(:,:,:), intent(in)  :: phi
  real, dimension(MDIM), intent(in) :: del
  real, dimension(:,:,:), intent(in) :: Coeff
  integer,dimension(MDIM), intent(in) :: lo, hi

  !---Local Variables
  integer :: i,j,k

  do k=lo(KAXIS),hi(KAXIS)
     do j=lo(JAXIS),hi(JAXIS)
        do i=lo(IAXIS),hi(IAXIS)
           rhs(i,j,k) = rhs(i,j,k) + (Coeff(i,j,k)/(del(IAXIS)**2))*(phi(i+1,j,k)+phi(i-1,j,k)-2.*phi(i,j,k))&
                + (Coeff(i,j,k)/(del(JAXIS)**2))*(phi(i,j+1,k)+phi(i,j-1,k)-2.*phi(i,j,k))
#if NDIM == 3           
           rhs(i,j,k) = rhs(i,j,k) + (Coeff(i,j,k)/(del(KAXIS)**2))*(phi(i,j,k+1)+phi(i,j,k-1)-2.*phi(i,j,k))
#endif           
        end do
     end do
  end do
  return

end subroutine Stencils_diffusionVar
