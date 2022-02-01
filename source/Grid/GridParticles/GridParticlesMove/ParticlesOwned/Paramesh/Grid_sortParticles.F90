!!****if* source/Grid/GridParticles/Grid_sortParticles
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
!!  Grid_sortParticles
!!
!! SYNOPSIS
!!
!!  call Grid_sortParticles(real(INOUT)    :: dataBuf(props,maxCount),
!!                          integer(IN)    :: props,
!!                          integer(INOUT) :: localCount,
!!                          integer(IN)    :: elementTypes,
!!                          integer(IN)    :: maxCount,
!!                          integer(OUT)   :: elementsPerBlk(MAXBLOCKS,elementTypes),
!!                          integer(IN)    :: attrib1,
!!                 optional,integer(IN)    :: attrib2)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Sorts the elements by block number. There are three types of blocknumber associated with
!!  elements data structures which have valid values, but are not valid blocknumber in the 
!!  mesh. These are "UNKNOWN", "LOST" and "NONEXISTENT". The sorter finds out the number of blocks 
!!  on the current processor in the mesh, and puts all elements associated with blocknumbers 
!!  in the range of 1 to localNumBlocks in the processor into the appropriate bins. In the bin
!!  for localNumBlocks+1 it puts all elements with block number = UNKNOWN, in the bin for
!!  localNumBlocks+2 it puts all elements with block number = LOST, and in localNumBlocks+3
!!  it puts elements with block number NONEXISTENT. If any other value
!!  is encountered in the particle attribute (given by attrib1) that is supposed to
!!  represent a block number, then
!!  * if preprocessor symbol DEBUG_PARTICLES is defined, the routine causes an abort
!!    with a message "undefined block number";
!!  * otherwise the behavior is undefined, and invalid memory access, segfaults,
!!    and/or program aborts should be expected.
!! 
!!
!! ARGUMENTS 
!!
!!  dataBuf : List of elements. It is two dimensional real array, the first dimension
!!              represents each element's properties, and second dimension is index to
!!              elements.
!!
!!  props : number of properties of each element in the dataBuf datastructure
!!
!!  localCount : While coming in it contains the current number of elements mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local elements might change, and the new value is put back into it.
!!               Must be less than or equal to the value of 'maxCount'.
!!  elementTypes  : Count of different types of elements in the simulation,
!!                  or more generally the largest possible bin number
!!                  for the second sort criterium (if 'attrib2' is used).
!!                  Must be 1 or greater in any non-trivial call,
!!                  i.e., when localCount > 0, even when attrib2 is
!!                  not present.
!!                  If attrib2 is present, 'elementTypes' must be
!!                  large enough to accomodate all the values of the
!!                  given attribute that are encountered in any elements.
!!  maxCount : This is parameter determined at runtime, and is the maximum number of local
!!                        elements that a simulation expects to have. All the arrays that hold
!!                        particles in the Particles unit are allocated based on this number.
!!             This is used for sizing the second dimension of the
!!             'dataBuf' argument, and thus must be greater than or
!!             equal to the value of 'localCount'.
!! elementsPerBlk : integer array containing number of elements on each blk.
!!                  If attrib2 is present, will contain number of
!!                  elements per (blk,type) combination, where 'type'
!!                  is the value of the property given by attrib2
!!                  ranging from 1 to elementTypes.
!!                  NOTE that even though the first dimension of elementsPerBlk
!!                  is declared with a range 1:MAXBLOCKS, where MAXBLOCKS is the
!!                  constant taken from Simulation.h, only the first localNumBlocks
!!                  elements within this range can be expected to be nonzero
!!                  on return.
!!  attrib1       : the primary property of the elements on which to sort them.
!!                  THIS MUST BE A BLOCKID-LIKE PROPERTY, in the following sense:
!!                     The values of this property must all be
!!                     integers values (of real type) from the set
!!                     {1, ..., localNumBlocks, UNKNOWN, LOST, NONEXISTENT},
!!                     where localNumBlocks is the value returned by a
!!                     Grid_getLocalNumBlks call and UNKNOWN, LOST,
!!                     NONEXISTENT are symbols for special (negative) values
!!                     defined in constants.h.
!!  attrib2       : if present, then the elements are first sorted on attrib2, and then 
!!                  within each group with a given attrib2 value, on the attrib1.
!!                  THIS MUST BE AN ELEMENTTYPE-LIKE PROPERTY, in the following sense:
!!                     The values of this property must all be
!!                     integers values (of real type) from the set
!!                     {1, ..., elementTypes}, where 'elementTypes' is
!!                     the value passed in as an argument (and which also
!!                     determines the size of the 'elementsPerBlk'
!!                     output array).
!!                     where localNumBlocks is the value returned by a
!!                     Grid_getLocalNumBlks call and UNKNOWN, LOST,
!!                     NONEXISTENT are symbols for special (negative) values
!!                     defined in constants.h.
!!                  This will be the case when Grid_sortParticle is
!!                  called to sort particles maintained and defined by the
!!                  Particles unit (with its usual definitions of
!!                  properties and types), with arguments
!!                  elementTypes=NPART_TYPES and attrib2=TYPE_PART_PROP,
!!                  since then all values of the "type" property should
!!                  have values in the range PART_TYPES_BEGIN (aka 1) ..
!!                  .. PART_TYPES_END (aka NPART_TYPES).
!!                  However, this routines makes no assummption that
!!                  attrib2 actually represents the "type" property and
!!                  will work propery as long as the requirements given
!!                  for elementTypes, attrib2, and the sizing of
!!                  elementsPerBlk are satisfied.
!!
!! SIDE EFFECTS
!!
!!   Elements of gr_ptSourceBuf (in the gr_ptData module) in the
!!   following range are overwritten, since this implementation used
!!   that buffer as temporary workspace: gr_ptSourceBuf(1:props,1:localCount).
!!
!! NOTES
!!   Currently this routine is called (only?) by io_writeParticles in permanent guardcell mode,
!!   by Particles_advance, by Particles_updateGridVar, and by the routines that map mesh to particles
!!   or particles to mesh.
!!   CAUTION : This routine must not be called when "gr_ptSourceBuf" needs to have valid values
!!             since it is used as scratch space by this routine.
!!
!!   The algorithm used for sorting requires (and may silently assume) that
!!     o   the attribute values corresponding to attrib1 in all given
!!         elements are all valid, i.e., they are either in the valid range
!!         1..localNumBlocks or have one of the special values NONEXISTENT,
!!         UNKNOWN, LOST;
!!     o   if attrib1 or attrib2 is BLK_PART_PROP then for each element it
!!         is either in the valid range 1..MAXBLOCKS,
!!         or has one of the special values NONEXISTENT, UNKNOWN, LOST; and
!!     o   if attrib2 is present, then the corresponding attribute values in
!!         all given elements are in the valid range 1..elementTypes.
!!
!!   When this routine is called by code in the Particles unit for particles
!!   owned and defined by it in the usual way, the following apply:
!!     o   attrib1 is BLK_PART_PROP;
!!     o   for each element, the BLK_PART_PROP attribute value is either
!!         in the valid range 1..MAXBLOCKS (actually, 1..localNumBlocks)
!!         or has one of the special values NONEXISTENT, UNKNOWN, LOST;
!!     o   attrib2 is either not present or is TYPE_PART_PROP;
!!     o   if attrib2 is present and is TYPE_PART_PROP, then
!!         elementTypes=NPART_TYPES and the corresponding attribute values
!!         in all given elements are in the valid range 1..elementTypes;
!!     o   if attrib2 is present and isBLK_PART_PROP, then for each element
!!         the attribute value is in the valid range 1..MAXBLOCKS,
!!         and 'elementTypes' has the value MAXBLOCKS (or larger).
!!
!!   Elements with a block value of UNKNOWN or LOST will be sorted to come
!!   after all those elements with block values in the range 1..MAXBLOCKS,
!!   first all those with UNKNOWN followed by all those with UNKNOWN block
!!   values (without any further reordering within those groups). These
!!   elements will *not* be included in the counts in the 'elementsPerBlk'
!!   output array, but will remain included in the count that is returned in
!!   the 'localCount' argument.
!!   Elements with a block value of NONEXISTING will be sorted to come after
!!   all those elements with block values in the range 1..MAXBLOCKS, after
!!   those elements with block values of UNKNOWN and LOST ((without any further
!!   reordering). These elements will *not* be included in the counts in the
!!   'elementsPerBlk' output array, and will also be excluded from the count
!!   that is returned in the 'localCount' argument.
!!   The caller may have to update its view of the number of elements (probably
!!   from the information in 'elementsPerBlk') after Grid_sortElements returns,
!!   if the block association of some of the elements can be UNKNOWN, LOST, or
!!   NONEXISTENT and if it agrees with this routine's understanding of which
!!   elements should count and which shouldn't. Because of the sorting that
!!   will have occurred, the caller can thus effectively drop some or all of
!!   those elements.
!!
!!   When this routine is used as envisioned, the value of attrib1 will usually be the same
!!   as that of the variable gr_ptBlk (in the gr_ptData module), but this is currently
!!   not used or checked.
!!
!! DEV:Missing error handling for various cases of invalid inputs.
!!***

#ifdef DEBUG_ALL
#define DEBUG_SORT_PARTICLES
#endif

subroutine Grid_sortParticles(dataBuf,props, localCount,elementTypes,maxCount,&
                              elementsPerBlk,attrib1, attrib2)
  use gr_ptData, ONLY : gr_ptSourceBuf, gr_ptLogLevel, gr_ptBlk, gr_ptKeepLostParticles
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getLocalNumBlks

  implicit none
#include "Simulation.h"
#include "constants.h"
#include "Particles.h"
  integer,intent(IN) ::  maxCount,props,elementTypes
  integer, intent(INOUT) :: localCount
  integer, parameter :: UNKNOWNBLK=1, LOSTBLK=2, NONEXISTENTBLK=3, &
       ALLBLOCKS=MAXBLOCKS+3
  real,dimension(props,maxCount),intent(INOUT) :: dataBuf
  integer,dimension(MAXBLOCKS,elementTypes),intent(OUT) :: elementsPerBlk
  integer,dimension(ALLBLOCKS,elementTypes)::localElementsPerBlk
  integer, intent(IN) :: attrib1
  integer,optional, intent(IN) :: attrib2
  integer, dimension((MAXBLOCKS)*elementTypes+4) :: pntr
  integer :: i, j,k,m,n,attrib,localNumBlocks,validParticles,lastBlkID
  logical :: pType
  
  integer :: validPointer, unknownID, lostID, nonexistentID

  localElementsPerBlk(:,:)=0

  validPointer=0

#ifdef DEBUG_SORT_PARTICLES
  print*,'Grid_sortParticles local/maxCount=',localCount,maxCount,',attrib1=',attrib1
#endif
  
  pType = present(attrib2)
  attrib = 1
  if(pType) attrib = attrib2 

  !! If no sorting is to be done on type, and there is only one block
  !! then all elements belong to the same block and no more sorting
  !! needs to be done

!  if((.not.pType).and.(localNumBlocks==1)) then
!     localElementsPerBlk(1,1)=localCount
!  else

     !! either there is more than one type, or more than one block, so sort
     !!  This is not an inplace sort, and so is order N in complexity
     !!  There are two passes, one to find out the number of elements
     !!  in each block, and the second one to put elements in order
     
     !! First copy all elements to a scrach, and count number 
     !! of elements per block  
     
     
  k=1 
  call Grid_getLocalNumBlks(localNumBlocks)
  unknownID=localNumBlocks+UNKNOWNBLK
  lostID=localNumBlocks+LOSTBLK
  nonexistentID=localNumBlocks+NONEXISTENTBLK
  
  do i = 1,localCount
     if (pType) k=int(dataBuf(attrib,i)) !DEV: Should protect against attrib > elementTypes here! - KW
     j = int(dataBuf(attrib1,i))
#ifdef DEBUG_PARTICLES
     if((j>localNumBlocks).or.(j==0).or.(j<LOST)) then
        call Driver_abort("Grid_sortParticles : undefined block number")
     end if
#endif
     validPointer=validPointer+1
     gr_ptSourceBuf(1:props,validPointer)=dataBuf(1:props,i)
     if(j==UNKNOWN)j=unknownID
     if(j==LOST)j=lostID
     if(j==NONEXISTENT)j=nonexistentID
     localElementsPerBlk(j,k)=localElementsPerBlk(j,k)+1
  end do

  localCount=validPointer
  !! Setup starting point for each block as it would
  !! once they are sorted
  
  pntr(1)=1
  k=2

  do j=1,elementTypes
     do i = 1,localNumBlocks
        pntr(k)=pntr(k-1)+localElementsPerBlk(i,j)
        k=k+1
     end do
  end do

  do i = localNumBlocks+1,nonexistentID
     pntr(k)=pntr(k-1)+sum(localElementsPerBlk(i,1:elementTypes))
     k=k+1
  end do
  
  !! Drop elements in their right place
  k=0

  do i = 1,localCount
     j=int(gr_ptSourceBuf(attrib1,i))

     if(j>0) then
        if(pType)k=int(gr_ptSourceBuf(attrib,i))-1 !DEV: Should protect against attrib > elementTypes here! - KW
        n=k*localNumBlocks+j
     else
        n=localNumBlocks*elementTypes
        if(j==UNKNOWN)n=n+UNKNOWNBLK
        if(j==LOST)n=n+LOSTBLK
        if(j==NONEXISTENT)n=n+NONEXISTENTBLK
     end if
     m=pntr(n)
     dataBuf(:,m)=gr_ptSourceBuf(:,i)
     pntr(n)=pntr(n)+1
  end do
  elementsPerBlk(1:MAXBLOCKS,:)=localElementsPerBlk(1:MAXBLOCKS,:)
  localCount=localCount-sum(localElementsPerBlk(nonexistentID,1:elementTypes))
  return
end subroutine Grid_sortParticles
