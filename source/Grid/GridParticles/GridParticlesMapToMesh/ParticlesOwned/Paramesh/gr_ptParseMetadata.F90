!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptParseMetadata
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
!!  gr_ptParseMetadata
!!
!! SYNOPSIS
!!
!!  gr_ptParseMetadata(integer(IN) :: bufferSize
!!                     real(IN), dimension(bufferSize) :: dataBuffer
!!                     integer(IN) :: headerPtr
!!                     integer(OUT), dimension(MDIM) :: negh
!!                     integer(OUT), dimension(MDIM) :: neghCornerID
!!                     integer(OUT), dimension(LOW:HIGH,MDIM) :: sectionCoords
!!                     integer(OUT) :: numbElements
!!
!! DESCRIPTION
!!
!! This routine extracts the metadata infomation from the communication
!! buffer.
!! 
!!
!! ARGUMENTS
!!               bufferSize: Size of the communication buffer.
!!               dataBuffer: The actual communication buffer.
!!               headerPtr: The position of the metadata for a particular section 
!!                          in the communication buffer.
!!               negh: The (block,proc,refinementLevel) array of the destination block.
!!               neghCornerID: The corner ID of the destination block.
!!               sectionCoords: The coordinates of the grid points in the destination
!!                              block.
!!               numbElements: The number of grid points in this section.
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptParseMetadata(bufferSize, dataBuffer, headerPtr, &
     negh, neghCornerID, sectionCoords, numbElements)

  use Driver_interface, ONLY : Driver_abort
  use Grid_data, ONLY : gr_meshMe

#include "Simulation.h"
#include "constants.h"
#include "gr_ptMapToMesh.h"

  implicit none
  integer, intent(IN) :: bufferSize
  real, dimension(bufferSize), intent(IN) :: dataBuffer
  integer, intent(IN) :: headerPtr
  integer, dimension(MDIM), intent(OUT) :: negh, neghCornerID
  integer, dimension(LOW:HIGH,MDIM), intent(OUT) :: sectionCoords
  integer, intent(OUT) :: numbElements

  if ((headerPtr + SIZE_HEADER - 1) > bufferSize) then
     call Driver_abort &
          ("[gr_ptParseMetadata]: Metadata extends beyond buffer")
  end if

  negh(BLKID)= int(dataBuffer(headerPtr+BLKID-1))
  negh(BLKPROC) = int(dataBuffer(headerPtr+BLKPROC-1))
  negh(REFLEVELDIF) = int(dataBuffer(headerPtr+REFLEVELDIF-1))
  neghCornerID(IAXIS) = int(dataBuffer(headerPtr+CORNERID-1))
  neghCornerID(JAXIS) = int(dataBuffer(headerPtr+CORNERID))
  neghCornerID(KAXIS) = int(dataBuffer(headerPtr+CORNERID+1))

  sectionCoords(LOW,IAXIS) = int(dataBuffer(headerPtr+COORDSID-1))
  sectionCoords(HIGH,IAXIS) = int(dataBuffer(headerPtr+COORDSID))
  sectionCoords(LOW,JAXIS) = int(dataBuffer(headerPtr+COORDSID+1))
  sectionCoords(HIGH,JAXIS) = int(dataBuffer(headerPtr+COORDSID+2))
  sectionCoords(LOW,KAXIS) = int(dataBuffer(headerPtr+COORDSID+3))
  sectionCoords(HIGH,KAXIS) = int(dataBuffer(headerPtr+COORDSID+4))

  numbElements = (sectionCoords(HIGH,IAXIS)-sectionCoords(LOW,IAXIS)+1) * &
       (sectionCoords(HIGH,JAXIS)-sectionCoords(LOW,JAXIS)+1) * &
       (sectionCoords(HIGH,KAXIS)-sectionCoords(LOW,KAXIS)+1)

  if ((headerPtr + SIZE_HEADER + numbElements - 1) > bufferSize) then
     print *, "Value:", headerPtr + SIZE_HEADER + numbElements, "size:", bufferSize
     call Driver_abort &
          ("[gr_ptParseMetadata]: Data extends beyond buffer")
  end if

end subroutine gr_ptParseMetadata
