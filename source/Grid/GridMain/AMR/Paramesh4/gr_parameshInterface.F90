!!****ih* source/Grid/localAPI/gr_parameshInterface
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
!!  gr_parameshInterface
!!
!! SYNOPSIS
!!
!!  use gr_parameshInterface
!!
!! DESCRIPTION
!!
!!  Interfaces for some subprograms private to the paramesh subunit
!!  implementation.
!!
!!***

module gr_parameshInterface
#include "constants.h"
#include "Simulation.h"  
  implicit none
  
  interface
     function gr_blockMatch(blkID, ntype, refinementLevel) result(match)
       implicit none
       integer, intent(IN)           :: blkID
       integer, intent(IN)           :: ntype
       integer, intent(IN), OPTIONAL :: refinementLevel
       logical                       :: match
     end function gr_blockMatch
  end interface
  
  interface 
     subroutine gr_pmGetListOfBlocks(blockType, listOfBlocks,count,refinementLevel,&
          region_bndBox, includePartialBlocks)
       
       integer, intent(in) :: blockType
       integer,dimension(MAXBLOCKS),intent(out) :: listOfBlocks
       integer,intent(out) :: count
       integer,intent(IN), optional :: refinementLevel
       real, dimension(LOW:HIGH,MDIM), intent(IN), optional :: region_bndBox
       logical, intent(IN), optional :: includePartialBlocks
       
     end subroutine gr_pmGetListOfBlocks
  end interface
end module gr_parameshInterface

