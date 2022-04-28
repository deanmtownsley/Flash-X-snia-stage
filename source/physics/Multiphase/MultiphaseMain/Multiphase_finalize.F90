!!****f* source/physics/Multiphase/MultiphaseMain/Multiphase_finalize
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
!!  Multiphase_finalize
!!
!!
!! SYNOPSIS
!!
!!  call Multiphase_finalize
!!
!! ARGUMENTS
!!
!!  none
!!
!! DESCRIPTION
!!
!!
!!***

subroutine Multiphase_finalize()

   use mph_evapInterface, ONLY: mph_evapFinalize

   implicit none

   call mph_evapFinalize()

end subroutine Multiphase_finalize

