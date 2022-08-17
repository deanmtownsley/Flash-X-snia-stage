!!****if* source/numericalTools/MoL/MoLMemory/ml_memZero
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
!!  NAME
!!
!!      ml_memZero
!!
!!  SYNOPSIS
!!
!!      call ml_memZero(integer, intent(in) :: dataStruct)
!!
!!  DESCRIPTION
!!
!!      Zero a specified MoL scratch-memory location
!!
!!  ARGUMENTS
!!
!!      dataStruct : Index of MoL scratch-memory to zero
!!
!!***
subroutine ml_memZero(dataStruct)
   use ml_memData

#include "MoL.h"

   implicit none

   integer, intent(in) :: dataStruct

   if (dataStruct .ge. MOL_INITIAL) scratch_data(:, :, :, :, :, dataStruct) = 0d0
end subroutine ml_memZero
