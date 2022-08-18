!!****if* source/numericalTools/MoL/MoLMemory/ml_memReleaseActiveRHS
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
!!      ml_memReleaseActiveRHS
!!
!!  SYNOPSIS
!!
!!      call ml_memReleaseActiveRHS()
!!
!!  DESCRIPTION
!!
!!      Release the active RHS
!!
!!  ARGUMENTS
!!
!!***
subroutine ml_memReleaseActiveRHS(irhs)
   use ml_memData, only: ml_activeRHS

#include "MoL.h"

   implicit none

   integer, intent(in) :: irhs

   ml_activeRHS = MOL_INVALID
end subroutine ml_memReleaseActiveRHS
