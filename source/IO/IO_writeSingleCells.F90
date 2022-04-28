!!****f* source/IO/IO_writeSingleCells
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
!!  IO_writeSingleCells
!!
!!  SYNOPSIS
!!
!!  call IO_writeSingleCells( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk))
!!
!!  DESCRIPTION 
!!      Write data for one or several single cells to a file.
!!
!! ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The current list of relevant blocks
!!
!! NOTES
!!
!!  The default implementation is a stub that does nothing.
!!  Simulations need to overwrite this with a real implementation
!!  if such output is desired.
!!
!!  An example implementation can be found in
!!  Simulation/SimulationMain/radflaHD/SupernovaRad1D .
!!***

subroutine IO_writeSingleCells(nblk, blklst)

  implicit none

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)

  return

end subroutine IO_writeSingleCells
