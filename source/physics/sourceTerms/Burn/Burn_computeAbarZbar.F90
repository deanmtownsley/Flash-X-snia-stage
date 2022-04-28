!!****f* source/physics/sourceTerms/Burn/Burn_computeAbarZbar
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
!! NAME
!!
!!  Burn_computeAbarZbar
!!
!! SYNOPSIS
!!
!!  call Burn_computeAbarZbar(real, dimension(:,:)(in) :: solnscalars,
!!                            real, dimension(:)(inout) :: abardata,
!!                            real, dimension(:)(inout) :: zbardata)
!!
!! DESCRIPTION
!!
!! Stub
!!
!! ARGUMENTS
!!
!!   solnscalars : solution scalars
!!
!!   abardata : abar info 
!!
!!   zbardata : zbar info
!!
!!
!!
!!***

subroutine Burn_computeAbarZbar(solnScalars, abarData, zbarData)

  implicit none

  real, dimension(:,:), intent(in)  :: solnScalars
  real, dimension(:), intent(inout) :: abarData, zbarData

end subroutine Burn_computeAbarZbar
