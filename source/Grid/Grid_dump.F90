!!****f* source/Grid/Grid_dump
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
!! Dumps the variables specified in integer array "var" to a file.
!! This routine doesn't require special resources, so it can be done from   
!! anywhere in the code, and is useful for diagnostic purposes.
!! This function can only be used on a single block per processor,
!! and mostly works with Uniform Grid. 
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
!!              solution data
!!  blockDesc :: Describes the  block to dump; holds the blockID
!!               in some Grid implementations. In the UG Grid,
!!               the blockID should always be 1.
!!               This argument may be unused, especially in UG.
!!  gcell :: indicates whether to include guardcells in the dump.
!!             
!! EXAMPLE
!!  
!!  num = 3  !dumping 3 variables
!!  var(1) = DENS_VAR
!!  var(2) = PRES_VAR
!!  var(3) = TEMP_VAR
!!  blockID = 1  ! local block number
!!  blockDesk % id = blockID           ! This is something of a hack
!!  solnData => NULL()
!!  gcell = .false.
!!
!!  call Grid_dump(var, num, solnData, blockDesc, gcell)
!!  
!!  will dump the interior cell values of density, pressure and temperature
!!  for local block number 1.
!!
!!  To explain the use of gcells, consider a global domain with 8x8 points
!!  mapped on 2x2 processors. each processor has blocks of size
!!  4x4.If there are 2 guard cells along each dimension, then the
!!  block size including guardcells is 8x8 and the distribution on
!!  four processors is as shown below ("*"
!!  are the interior points and the )"o" are guard cells.
!!
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!
!!  If gcell is true then dump is of size 16x16 and looks like
!!
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!
!!  and if gcell is false the dump is of size 8x8 and looks like
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!    
!!
!! NOTES
!!  DENS_VAR, PRES_VAR, TEMP_VAR etc are #defined values in Simulation.h
!!  indicating the index in the physical data array.
!!  The routine calling Grid_dump will need to include Simulation.h .
!!
!!***

subroutine Grid_dump(var,num, solnData,blockDesc, gcell)
  use Grid_tile, ONLY : Grid_tile_t

implicit none
  integer, intent(IN) :: num
  integer, dimension(num), intent(IN) :: var
  real,dimension(:,:,:,:),pointer     :: solnData
  type(Grid_tile_t), intent(in)  :: blockDesc
  logical, intent(IN) :: gcell

  return
end subroutine Grid_dump
