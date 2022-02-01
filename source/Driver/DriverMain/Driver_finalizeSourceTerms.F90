!!****f* source/Driver/DriverMain/Driver_finalizeSourceTerms
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
!!  Driver_finalizeSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_finalizeSourceTerms(logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Finalizes all source terms Units by
!!   calling their respective termination  routines,
!!   viz. Stir_finalize, Burn_finalize, Heat_finalize, Cool_finalize, etc.
!!
!!
!! ARGUMENTS
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_finalizeSourceTerms( restart )

  use Burn_interface, ONLY:  Burn_finalize
  use RadTrans_interface, ONLY : RadTrans_finalize

  implicit none

  logical, intent(in) :: restart

  call Burn_finalize()
  call RadTrans_finalize()

end subroutine Driver_finalizeSourceTerms
