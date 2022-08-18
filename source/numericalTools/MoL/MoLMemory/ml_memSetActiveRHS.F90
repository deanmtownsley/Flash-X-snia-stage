!!****if* source/numericalTools/MoL/MoLMemory/ml_memSetActiveRHS
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
!!      ml_memSetActiveRHS
!!
!!  SYNOPSIS
!!
!!      call ml_memSetActiveRHS(integer, intent(in) :: irhs)
!!
!!  DESCRIPTION
!!
!!      Set the active RHS to be pointed to when MOL_RHS is requested
!!
!!  ARGUMENTS
!!
!!      irhs : the RHS to set as active
!!
!!***
subroutine ml_memSetActiveRHS(irhs)
   use ml_memData, only: ml_activeRHS

   implicit none

   integer, intent(in) :: irhs

   ml_activeRHS = irhs
end subroutine ml_memSetActiveRHS
