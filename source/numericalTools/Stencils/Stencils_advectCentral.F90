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
!!**

#include "constants.h"

subroutine Stencils_advectCentral(rhs,phi,u,v,w,delta,lo,hi,face)
  implicit none
  real, dimension(:,:,:), intent(inout):: rhs
  real, dimension(:,:,:), intent(in) :: phi,u,v,w
  real, dimension(MDIM), intent(in) :: delta
  integer, dimension(MDIM), intent(in) :: lo,hi
  integer, dimension(MDIM+1), intent(in) :: face
end subroutine Stencils_advectCentral
