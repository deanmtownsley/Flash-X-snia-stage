!!****if* source/IO/localAPI/io_rescaleCellBoxes
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
!!  io_rescaleCellBoxes
!!
!!
!! SYNOPSIS
!!
!!  call io_rescaleCellBoxes() 
!!
!!
!! DESCRIPTION
!!
!!  Apply linear transformations (stretching / shifting) to the cell
!!  coordinates.
!!
!!  To be called immediately after io_readData has filled in various
!!  Paramesh metadata arrays with information from a checkpoint file
!!  upon restart.
!!
!!  This is a hack.
!!
!!  It is the user's responsibility to adjust the xmin,xmax,ymin,...,zmax
!!  Runtime parameters appropriately when coordinate rescaling is used.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  This manipulates data structures owned by PARAMESH. Therefore it does not
!!  do anything if a different Grid implementation is used.
!!***


subroutine io_rescaleCellBoxes()

  implicit none


end subroutine io_rescaleCellBoxes
