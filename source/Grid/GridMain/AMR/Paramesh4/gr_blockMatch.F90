!!****if* source/Grid/GridMain/paramesh/gr_blockMatch
!! NOTICE
!!  Copyright 2020 The University of Chicago
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
!!  gr_blockMatch
!!
!! SYNOPSIS
!!
!!  ans = gr_blockMatch(integer(IN)  :: blkid,
!!                      integer(IN)  :: ntype,
!!                      integer(IN),OPTIONAL  :: refinementlevel)
!!
!! DESCRIPTION
!!
!!  Test whether a block matches a criterion.  Paramesh only.
!!
!! ARGUMENTS
!!
!!   blkid : block ID
!!
!!   ntype : block type, or type of requested match.
!!
!!              For the paramesh implementation, valid values are :
!!
!!              ALL_BLKS    all local blocks on a processor
!!
!!              IBDRY_BLKS  blocks that are on physical boundary along IAXIS
!!              JBDRY_BLKS  blocks that are on physical boundary along JAXIS
!!              KBDRY_BLKS  blocks that are on physical boundary along KAXIS
!!              ANY_BDRY_BLKS  blocks that have any of their faces
!!                          on a physical boundary.
!!              ACTIVE_BLKS all currently active blocks; in the paramesh
!!                          context this means parent and leaf blocks.
!!
!!              LEAF, PARENT_BLK or ANCESTOR  representing
!!                          the type of node in the oct-tree managing the blocks.
!!              REFINEMENT refinement level, the optional refinementlevel
!!                                           argument must be present and valid.
!!
!!              All of these constants are defined in constants.h
!!
!!   refinementlevel : If valid, refinement level to be used as match criterion.
!!                     Values greater that 0 are considered valid.
!!
!!                     If ntype is REFINEMENT, the caller must provide a present
!!                     and valid refinementlevel, and all blocks at that level
!!                     are considered matching.
!!                     For other choices of ntype, the refinementlevel is used
!!                     as additional criterion if present and valid.
!!
!! NOTES
!!
!!   Only implemented for the PARAMESH AMR Grid implementation
!!
!! RETURN TYPE
!!
!!  LOGICAL ans  - indicates whether the block matches the criterion
!!
!! SEE ALSO
!!
!!   gr_pmGetListOfBlocks
!!
!!***

function gr_blockMatch(blkID,ntype,refinementLevel) result(match)

  use tree, ONLY : nodetype,lnblocks,neigh,lrefine

  implicit none

#include "constants.h"
  logical               :: match
  integer,intent(in   ) :: blkID
  integer,intent(in   ) :: ntype
  integer,intent(in   ),OPTIONAL :: refinementLevel

  integer :: i
  logical :: isBnd

  i = blkID
  if (i > lnblocks) then
     match = .FALSE.
     return                   !block number is not valid, return immediately!
  end if

  select case (ntype)
  case (LEAF,PARENT_BLK,ANCESTOR)
     match = (nodetype(i)==ntype)
  case(ALL_BLKS)
     match = .TRUE.
  case (IBDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        isBnd = (neigh(1,1,i) .LE. PARAMESH_PHYSICAL_BOUNDARY)
        isBnd = isBnd.or.(neigh(1,2,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
        match = isBnd
     end if
  case (JBDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        isBnd = (neigh(1,3,i) .LE. PARAMESH_PHYSICAL_BOUNDARY)
        isBnd = isBnd.or.(neigh(1,4,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
        match = isBnd
     end if
  case(KBDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        isBnd = (neigh(1,5,i) .LE. PARAMESH_PHYSICAL_BOUNDARY)
        isBnd = isBnd.or.(neigh(1,6,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
        match = isBnd
     end if
  case(ANY_BDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        match = ANY(neigh(1,1:6,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
     end if
  case(ACTIVE_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
  case(REFINEMENT)
     if (present(refinementLevel)) then
        if (refinementLevel > 0) then
           match = .TRUE.          ! test for matching level is done below
        else
           call Driver_abort("[gr_blockMatch] With ntype=REFINEMENT, argument refinementlevel must be present and > 0.")
        end if
     else
        call Driver_abort("[gr_blockMatch] With ntype=REFINEMENT, optional argument refinementlevel must be present.")
     end if
  case default
     match = .FALSE.
     call Driver_abort("[gr_blockMatch] ntype argument not recognized")
  end select

  if (match .AND. present(refinementLevel)) then
     if (refinementLevel > 0) then
        match = (lrefine(i) == refinementLevel)
     end if
  end if

end function gr_blockMatch
