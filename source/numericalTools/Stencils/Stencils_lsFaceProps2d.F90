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
subroutine Stencils_lsFaceProps2dFixedSmear(phi,iPropX,iPropY,lsScalarProp,ix1,ix2,jy1,jy2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY
end subroutine Stencils_lsFaceProps2dFixedSmear

subroutine Stencils_lsFaceProps2dUserSmear(phi,iPropX,iPropY,lsScalarProp,ix1,ix2,jy1,jy2,iSmear)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY
    real, intent(in) :: iSmear
end subroutine Stencils_lsFaceProps2dUserSmear
