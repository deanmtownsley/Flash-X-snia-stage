!!****if* source/numericalTools/MoL/MoLMain/ml_warn
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
!!      ml_warn
!!
!!  SYNOPSIS
!!
!!      call ml_warn(character, intent(in) :: msg(:))
!!
!!  DESCRIPTION
!!
!!      Print a warning message and optionally abort
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine ml_warn(msg)
   use MoL_data

   use Driver_interface, only: Driver_abort

#include "MoL.h"
#include "constants.h"

   implicit none

   character(len=*), intent(in) :: msg

   if (MoL_abortOnWarn) then
      ! Always print warning if abort-on-warning is turned on
      call Driver_abort("[MoL] WARNING: "//msg)
   end if

   if ((MoL_verbosity .ge. MOL_VERBOSE_WARN) .and. (MoL_mpiRank .eq. MASTER_PE)) then
      print *, "[MoL] WARNING: "//msg
   end if
end subroutine ml_warn
