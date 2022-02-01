!!****f* source/physics/Hydro/Hydro_sendOutputData
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
!!  Hydro_sendOutputData
!!
!! SYNOPSIS
!!  
!!  Hydro_sendOutputData()
!!  
!! DESCRIPTION 
!!   
!!  Sends any data from the Hydro unit (stored in Hydro_data) to the
!!  IO unit, through the IO_setScalar interface, for checkpointing.
!!  Can retrieve the data upon restart in Hydro_init through the 
!!  IO_getScalar interface.
!!  
!!  This is the API stub implementation. 
!!
!! ARGUMENTS 
!!
!! USES
!!  
!!  IO_setScalar, which puts a scalar into the IO unit so it can make 
!!  it into a checkpoint. 
!!
!! USED BY
!!
!!  IO_updateScalars, which calls this function to let Hydro know 
!!  it's time to checkpoint.
!!
!! EXAMPLE
!!
!!  Here's what an implementation of this function might do:
!!
!!  USE IO_interface, ONLY : IO_setScalar
!!  call IO_setScalar("time", simTime)
!!
!! 
!!***

subroutine Hydro_sendOutputData()

implicit none
end subroutine Hydro_sendOutputData

