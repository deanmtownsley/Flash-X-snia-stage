!!****if* source/Driver/DriverMain/Driver_initSourceTerms
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
!!  Driver_initSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_initSourceTerms(logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes all source terms Units by 
!!   calling their respective initialization routines
!!  
!! ARGUMENTS
!!   myPE - current processor number
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_initSourceTerms( restart )

  use Burn_interface, ONLY:  Burn_init
  use RadTrans_interface, ONLY : RadTrans_init

  implicit none
  
  logical, intent(in) :: restart

  call Burn_init()
  call RadTrans_init()

end subroutine Driver_initSourceTerms
