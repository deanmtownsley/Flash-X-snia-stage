!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptSearchBlk
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
!!  gr_ptSearchBlk
!!
!! SYNOPSIS
!!
!!  gr_ptSearchBlk(integer,dimension(MDIM), intent(IN) :: cornerID, &
!!                 integer,dimension(BLKID:REFLEVELDIF), intent(INOUT) :: negh)
!!
!! DESCRIPTION
!!
!!  Routine to check whether a paricular block on this processor matches the 
!!  corner ID argument.  If there is a match, the negh array is updated.
!!
!! ARGUMENTS
!!               cornerID:   The corner ID we will try to match.
!!               negh:   Array containing information about the destination block(s)
!!
!! PARAMETERS
!! 
!!***


subroutine gr_ptSearchBlk(cornerID,negh)

#include "Simulation.h"
#include "constants.h"
#include "gr_ptMapToMesh.h"

  use tree, ONLY : lnblocks,nodetype
  use Grid_data, ONLY : gr_oneBlock, gr_meshMe

  implicit none  

  integer,dimension(MDIM), intent(IN) :: cornerID
  integer,dimension(BLKID:REFLEVELDIF), intent(INOUT) :: negh
  integer :: lb
  logical :: notFound
  integer, dimension(NDIM) :: diff
  
  lb=0
  
  notFound=.true.
  do while ((lb<lnblocks).and.notFound)
     lb=lb+1
     if(nodetype(lb)==LEAF) then
        diff=cornerID(1:NDIM)-gr_oneBlock(lb)%cornerID(1:NDIM)
        if(maxval(abs(diff))==0)then
           negh(BLKID)=lb
           negh(BLKPROC)=gr_meshMe
           notFound=.false.
        end if
     end if
  end do

  return
end subroutine gr_ptSearchBlk
