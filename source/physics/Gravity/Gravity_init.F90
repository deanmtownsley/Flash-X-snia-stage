!!****f* source/physics/Gravity/Gravity_init
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
!!  Gravity_init
!!  
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables in the Gravity unit, which are typically the 
!!  runtime parameters.  This routine must be called once by Driver_initAll.F90. 
!!  Calling multiple times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!! 
!!  useGravity  BOOLEAN true  Controls turning on/off the compiled gravity unit
!!
!! NOTES
!!   
!!  Each implementation of Gravity has its own runtime parameters.  Be sure to check
!!  the documentation or Config files to see them.
!!
!!***

subroutine Gravity_init ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abort

  implicit none
   

  logical, save :: testUseGravity

!==============================================================================

  !! It is a failure to invoke the stub when useGravity is set TRUE.

  call RuntimeParameters_get ("useGravity", testUseGravity)
  if (testUseGravity) then
     call Driver_abort("Gravity unit seems not to be compiled in, and the Gravity_init stub does not &
          &allow the value of useGravity to be TRUE.")
  end if

end subroutine Gravity_init
