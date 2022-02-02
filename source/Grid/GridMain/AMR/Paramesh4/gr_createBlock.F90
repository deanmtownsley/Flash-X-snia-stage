!!****if* source/Grid/GridMain/paramesh/gr_createBlock
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
!!  gr_createBlock
!!
!! SYNOPSIS
!!
!!  call gr_createBlock(real(IN) :: xBlkMin,
!!                      real(IN) :: xBlkMax,
!!                      real(IN) :: yBlkMin,
!!                      real(IN) :: yBlkMax,
!!                      real(IN) :: zBlkMin,
!!                      real(IN) :: zBlkMax,
!!                      integer(IN) :: blockID)
!!  
!! DESCRIPTION 
!!   create a block and initialize its place in the paramesh
!!   tree data structure. The block is bound by the coordinates
!!   (xBlkmin,yBlkMin,zBlkMin) and (xBlkmax,yBlkMax,zBlkMax) on the
!!   physical domain.
!!  
!! ARGUMENTS 
!!
!!   xBlkMin - lower bound coordinate along IAXIS
!!   xBlkMax - upper bound coordinate along IAXIS
!!   yBlkMin - lower bound coordinate along JAXIS
!!   yBlkMax - upper bound coordinate along JAXIS
!!   zBlkMin - lower bound coordinate along KAXIS
!!   zBlkMax - upper bound coordinate along KAXIS
!!   blockID - the identifier for the current block 
!!
!!***

subroutine gr_createBlock(xBlkmin, xBlkmax, &
     yBlkmin, yBlkmax, zBlkmin, zBlkmax,blockID)

!==============================================================================
  use tree, ONLY : lnblocks,lrefine,nodetype,parent,child,bsize,bnd_box,coord,&
                   work_block,gr_btWorkDefaultLeaf
  use Driver_interface, ONLY : Driver_abort
  implicit none
# include "Simulation.h"

  real,intent(IN) :: xBlkmin, xBlkmax, yBlkmin, yBlkmax, zBlkmin, zBlkmax
  integer,intent(IN) :: blockID
  
  
  !=====================================================================
  
  
  if (lnblocks .eq. MAXBLOCKS) then
     call Driver_abort ("createBlock:  cannot add another block")
  endif
  
  !            Increment the local block # counter and set tree information.

  lnblocks = lnblocks + 1
  lrefine(blockID)   = 1
  nodetype(blockID)  = 1
  parent(:,blockID)  = -1
  child(:,:,blockID) = -1
  work_block(blockID) = gr_btWorkDefaultLeaf

!            Set the block's dimensions, bounding box, and center coords.

  bsize(1,blockID) = xBlkmax - xBlkmin
  bsize(2,blockID) = yBlkmax - yBlkmin
  bsize(3,blockID) = zBlkmax - zBlkmin

  bnd_box(1,1,blockID) = xBlkmin
  bnd_box(2,1,blockID) = xBlkmax

  bnd_box(1,2,blockID) = yBlkmin
  bnd_box(2,2,blockID) = yBlkmax

  bnd_box(1,3,blockID) = zBlkmin
  bnd_box(2,3,blockID) = zBlkmax

  coord(:,blockID) = bnd_box(1,:,blockID) + 0.5*bsize(:,blockID)

!====================================================================

  return
end subroutine gr_createBlock
