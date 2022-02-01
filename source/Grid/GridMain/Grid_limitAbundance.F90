!!****if* source/Grid/GridMain/Grid_limitAbundance
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
!!  Grid_limitAbundance
!!
!!
!! SYNOPSIS
!!
!!  Grid_limitAbundance(integer(IN)   :: blkLimits(2,MDIM),
!!                       real, pointer :: solnData(:,:,:,:))
!!
!!
!! DESCRIPTION
!!
!!  Limit each abundance so it falls between smallx and 1.
!!
!!  This routine is called automatically by Hydro and MHD in FLASH 
!!  if the irenorm runtime parameter is set to 0.
!!
!!
!! ARGUMENTS
!!
!!  blkLimits - the index limits for internal zones of the block to limited
!!  solnData -  Pointer to the block to be limited
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Grid_limitAbundance(blkLimits,solnData)


  use Grid_data, ONLY : gr_smallx

  implicit none
#include "constants.h"
#include "Simulation.h"

  integer, dimension(2,MDIM), INTENT(in) :: blkLimits
  real, POINTER :: solnData(:,:,:,:)
      

  integer :: i, j, k, n
  
  
! we will access the solution data through a pointer return by the 
! accessor method

! loop over the zones in the block
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           do n = 0,NSPECIES-1
                  
! make sure the abundances are between gr_smallx and 1
              solnData(SPECIES_BEGIN+n,i,j,k) =  & 
                   max(gr_smallx, min(1.e0,solnData(SPECIES_BEGIN+n,i,j,k)))
                  
           enddo

        enddo
     enddo
  enddo

  return
end subroutine Grid_limitAbundance
 
