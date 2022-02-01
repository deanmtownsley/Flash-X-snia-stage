!!****if* source/Grid/GridMain/AMR/Paramesh4/Grid_dump
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
!!  Grid_dump
!!
!! SYNOPSIS
!!
!!  call Grid_dump(integer(IN) :: var(num),
!!                 integer(IN) :: num,
!!                 real,dimension(:,:,:,:),POINTER :: solnData,
!!                 Type(Grid_tile_t)(IN) :: blockDesc,
!!                 logical(IN) :: gcell)
!!
!! DESCRIPTION 
!!  
!! Dumps the variables specified in "var" to a file. Can be done from 
!! anywhere in the code, and is useful for diagnostic purposes.
!! With paramesh this function doesn not work in parallel, but works
!! only with a single block
!!  
!! ARGUMENTS 
!!
!! This is not an essential interface, only provided as a convenience
!! for those who need to dump for debugging; test applications may
!! want to override with a customized implementation.
!!
!! ARGUMENTS
!!
!!  var :: 1D integer array containing the indices of the variables
!!         to be dumped (can be conveniently given using the variable
!!         names defined in Simulation.h)
!!  num :: number of variables being dumped.
!!  solnData :: an associated pointer that points to the block's
!!              solution data.
!!              This argument may be unused in some Grid implementations.
!!  blockDesc :: Describes the  block to dump; holds the blockID
!!               (local number of the block) in some Grid
!!               implementations.
!!  gcell :: indicates whether to include guardcells in the dump.
!!             
!! EXAMPLE
!!  
!!  num = 3  !dumping 3 variables
!!  var(1) = DENS_VAR
!!  var(2) = PRES_VAR
!!  var(3) = TEMP_VAR
!!  blockID = 1  ! local block number
!!  gcell = .false.
!!
!!  call Grid_dump(var, num, blockID, gcell)
!!  
!!  will dump the interior cell values of density, pressure and temperature
!!  for local block number 1.
!!
!! NOTES
!!  DENS_VAR, PRES_VAR, TEMP_VAR etc are #defined values in Simulation.h
!!  indicating the index in the physical data array.
!!  The routine calling Grid_dump will need to include Simulation.h .
!! 
!!  This routine works with Paramesh only under very special circumstances
!!***

subroutine Grid_dump(var,num, solnData,blockDesc, gcell)

  use gr_specificData, ONLY : gr_ilo, gr_ihi, gr_jlo, gr_jhi, &
       gr_klo, gr_khi, gr_iloGC, gr_ihiGC, gr_jloGC, gr_jhiGC, &
       gr_kloGC, gr_khiGC

  use physicaldata, ONLY :unk
  use Grid_tile, ONLY : Grid_tile_t

#include "Simulation.h"

  implicit none

  integer, intent(IN) :: num
  integer, dimension(num), intent(IN) :: var
  real,dimension(:,:,:,:),pointer     :: solnData
  type(Grid_tile_t), intent(in)  :: blockDesc
  logical, intent(IN) :: gcell
  

  character(len=80) :: ff1
  integer,dimension(4), save :: filecount = 0
  integer :: blockID
  integer :: i,count

  integer,parameter :: bxn=GRID_IHI_GC-GRID_ILO_GC+1
  integer,parameter :: byn=GRID_JHI_GC-GRID_JLO_GC+1 
  integer,parameter :: bzn=GRID_KHI_GC-GRID_KLO_GC+1 
  real,dimension(1,bxn,byn,bzn,1)::scratch

  blockID = blockDesc % id

  if(.not. gcell) then
     print '(8F11.6)', unk(var(1), gr_ilo:gr_ihi, gr_jlo:gr_jhi, gr_klo:gr_khi, blockID)

  else
     print '(16F7.3)', unk(var(1), gr_iloGc:gr_ihiGc, gr_jloGc:gr_jhiGc, gr_kloGc:gr_khiGc, blockID)
  end if

 

  !print '(8F11.6)', unk(var(1), gr_ilo:gr_ihi, gr_jlo:gr_jhi, gr_klo:gr_khi, blockID) 
 
  ff1 = "FL2"//char(48+filecount(4))//char(48+filecount(3))//&
       char(48+filecount(2))//char(48+filecount(1))
  !print*,'filecount',filecount
  filecount(1) = filecount(1) + 1
  do i = 1,3
     if(filecount(i)==10)then
        filecount(i) = 0
        filecount(i+1)=filecount(i+1)+1
     end if
  end do
  count = 8*bxn*byn*bzn
  if((filecount(3) == 0).and.filecount(2) == 0) then
  
!!! create the data structure for outputting data

     open(41,file=ff1,access='direct',form='unformatted',recl=count)
     do i = 1,num
        !scratch(1,:,:,:,1)=unk(var(i),:,:,:,blockID)
        !write(41,rec=i)scratch
     end do
     close(41)
  end if
  return
end subroutine Grid_dump
