!!****f* source/IO/IO_outputInitial
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
!!  IO_outputInitial
!!
!!
!! SYNOPSIS
!!
!!  IO_outputInitial(integer(in) :: nbegin,
!!                   real(in) :: initialSimTime  
!!                  
!!
!!
!! DESCRIPTION
!!
!!  This routine is called before the main timestep loop.  It outputs the 
!!  initial data to a checkpoint file and plotfile, and particle plotfiles
!!
!!  If particles are not included a stub (empty) routine will be called.
!!
!!
!! ARGUMENTS
!!
!!  nbegin - initial step of simulation
!!  initialSimTime - initial simulation time
!!
!!
!!***

subroutine IO_outputInitial( nbegin, initialSimTime)

  use IO_interface, ONLY : IO_writeIntegralQuantities

  implicit none
  integer, intent(in) ::  nbegin
  real, intent(in)    :: initialSimTime

  call IO_writeIntegralQuantities( 1, initialSimTime)

end subroutine IO_outputInitial
