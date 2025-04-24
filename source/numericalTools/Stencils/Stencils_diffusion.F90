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
subroutine Stencils_diffusionConst(rhs, phi, del, Coeff, lo, hi)
  implicit none
  real, dimension(:,:,:), intent(inout) :: rhs
  real, dimension(:,:,:), intent(in)  :: phi
  real, dimension(3), intent(in) :: del
  real, intent(in) :: Coeff
  integer,dimension(3), intent(in) :: lo, hi
end subroutine Stencils_diffusionConst

subroutine Stencils_diffusionVar(rhs, phi, del, Coeff, lo, hi)
  implicit none
  real, dimension(:,:,:), intent(inout) :: rhs
  real, dimension(:,:,:), intent(in)  :: phi
  real, dimension(3), intent(in) :: del
  real, dimension(:,:,:), intent(in) :: Coeff
  integer, dimension(3), intent(in) :: lo, hi
end subroutine Stencils_diffusionVar
