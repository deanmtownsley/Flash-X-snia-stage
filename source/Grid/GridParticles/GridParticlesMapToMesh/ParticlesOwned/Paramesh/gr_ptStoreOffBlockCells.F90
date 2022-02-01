!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptStoreOffBlockCells
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
!!  gr_ptStoreOffBlockCells
!!
!! SYNOPSIS
!!
!!   gr_ptStoreOffBlockCells (integer,dimension(MAXBLOCKS), intent(IN) :: particlesPerBlk
!!                            integer,dimension(MAXBLOCKS), intent(IN) :: blockList
!!                            integer,intent(IN) :: blockCount
!!                            integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
!!                            integer,dimension(MDIM), intent(IN) :: blkSize
!!                            integer,dimension(MDIM), intent(IN) :: guard
!!                            integer, intent(OUT) :: BufferSize)
!!
!! DESCRIPTION
!!
!! This subroutine has two purposes.  The first is to calculate the maximum 
!! size of the send / receive buffer on each processor, and the second is to 
!! store useful information in a user defined type (UDT).  
!!
!! This is a stub!
!!
!! ARGUMENTS
!!
!!                          particlesPerBlk: Number of particles residing on each block.
!!                          blockList: List of all leaf blocks existing on this processor.
!!                          blockCount: Number of leaf blocks existing on this processor.
!!                          blkLimitsGC: Upper and lower indicies of source block including guard cells.
!!                          blkSize: Size of the source block (same for each block).
!!                          guard: Number of guard cells for the source block.
!!                          BufferSize: The size of the send / receive buffer.
!! 
!!***

subroutine gr_ptStoreOffBlockCells(particlesPerBlk, blockList, blockCount, blkLimitsGC, blkSize, guard, BufferSize)

  use gr_ptMapData, ONLY : gr_ptDomain
  use Grid_interface, ONLY : Grid_getBlkCornerID, Grid_getBlkBoundBox, Grid_getDeltas
  use gr_ptInterface, ONLY : gr_ptFindNegh, gr_ptGetSrcDestCoords
  use Grid_data, ONLY : gr_meshMe

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Flashx_mpi.h"
#include "gr_ptMapToMesh.h"

  integer,dimension(MAXBLOCKS), intent(IN) :: particlesPerBlk, blockList
  integer,intent(IN) :: blockCount
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
  integer,dimension(MDIM), intent(IN) :: blkSize, guard
  integer, intent(OUT) :: BufferSize

  BufferSize = 0

end subroutine gr_ptStoreOffBlockCells
