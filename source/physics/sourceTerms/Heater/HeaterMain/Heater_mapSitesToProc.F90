!!***if* source/physics/sourceTerms/Heater/HeaterMain/Heater_mapSitesToProc
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  Heater_mapSitesToProc
!!
!!
!! SYNOPSIS
!!  Heater_mapSitesToProc(logical(IN) :: initial, logical(IN) :: gridChanged)
!!
!! ARGUMENTS
!!  initial - If initial
!!  gridChanged - If grid changed
!!
!!***
#include "constants.h"
#include "Simulation.h"

subroutine Heater_mapSitesToProc(initial, gridChanged)

   use Grid_interface, ONLY: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t
   use Heater_type, ONLY: Heater_type_t
   use Heater_data, ONLY: htr_meshMe, htr_heaterInfo, htr_numHeaters
   use Driver_interface, ONLY: Driver_abort
   use Timers_interface, ONLY: Timers_start, Timers_stop

   implicit none
   include "Flashx_mpi.h"
   logical, intent(in), optional :: gridChanged
   logical, intent(in), optional :: initial

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc
   type(Heater_type_t), pointer :: heater
   real, dimension(:), allocatable :: xSite, ySite, zSite
   integer :: numSitesProc, htr, isite, blockCount
   real :: boundBox(LOW:HIGH, 1:MDIM)
   integer :: numSitesBlk(MAXBLOCKS)

   call Timers_start("Heater_mapSitesToProc")

   ! Check if this procedure is being called after a grid change has occured.
   ! Currently this feature is no implemented and will abort the simulation
   ! if gridChanged == .TRUE.
   !
   ! when this procedure is called with gridChanged == .FALSE. It will print
   ! out a message indicating that it will not do anything
   if (present(gridChanged)) then
      if (gridChanged) then
         call Driver_abort("[Heater_mapSitesToProc] Not Implemented for gridChanged == .TRUE.")
      else
         if (htr_meshMe .eq. MASTER_PE) print *, "[Heater_mapSitesToProc] Doing nothing, gridChanged == .FALSE."
      end if

      ! Check if this procedure is being called at the begining of the time loop after
      ! the heater file has been read.
      ! if initial == .TRUE. then distribute the sites read from the file across processors
      ! and blocks that intersect with them. See the comments in the procedure below
      ! for details
      !
      ! if initial == .FALSE. it will print out a message indicating that it will not do anything
   else if (present(initial)) then
      if (initial) then

         ! First loop over heaters
         do htr = 1, htr_numHeaters

            ! Get pointer to the current heater
            ! and set the loop variable for numSites
            heater => htr_heaterInfo(htr)
            numSitesProc = 0
            numSitesBlk = 0

            ! Second loop over all sites that
            ! have been read from the heater file
            do isite = 1, heater%numSitesAll

               ! Now loop over all blocks and set blockCount to zero
               blockCount = 0
               call Grid_getTileIterator(itor, nodetype=LEAF, tiling=.FALSE.)
               do while (itor%isValid())

                  ! Increase blockCount and get current tile descriptor
                  blockCount = blockCount+1
                  call itor%currentTile(tileDesc)

                  ! Calculate boundBox for the tile descriptor and then check
                  ! If the site lies within the block an perform the logic
                  ! of mapping initial site distribution to block and processor
                  call tileDesc%boundBox(boundBox)

                  if (boundBox(HIGH, IAXIS) .ge. heater%xSiteInit(isite) .and. &
                      boundBox(LOW, IAXIS) .le. heater%xSiteInit(isite) .and. &
                      boundBox(HIGH, JAXIS) .ge. heater%ySiteInit(isite) .and. &
                      boundBox(LOW, JAXIS) .le. heater%ySiteInit(isite) &
#if NDIM == MDIM
                      .and. &
                      boundBox(HIGH, KAXIS) .ge. heater%zSiteInit(isite) .and. &
                      boundBox(LOW, KAXIS) .le. heater%zSiteInit(isite) &
#endif
                      ) then

                     ! Increment numSitesProc that will also act as current index
                     ! and abort with numSitesProc > HTR_MAX_NUMSITES
                     numSitesProc = numSitesProc+1

                     if (numSitesProc > HTR_MAX_NUMSITES) then
                        call Driver_abort("[Heater_mapSitesToProc] Increase HTR_MAX_NUMSITES")
                     end if

                     ! Set local xyz-locations of the site
                     heater%xSiteProc(numSitesProc) = heater%xSiteInit(isite)
                     heater%ySiteProc(numSitesProc) = heater%ySiteInit(isite)
                     heater%zSiteProc(numSitesProc) = heater%zSiteInit(isite)

                     ! Increment numSitesBlk and update the site map
                     ! on the current processor
                     numSitesBlk(blockCount) = numSitesBlk(blockCount)+1
                     heater%siteMapOnProc(blockCount, numSitesBlk(blockCount)) = numSitesProc

                  end if
                  ! Increment the block iterator
                  call itor%next()
               end do ! End loop over blocks
               call Grid_releaseTileIterator(itor)
            end do ! End loop over sites

            ! Set total number of sites present on currrent process
            heater%numSitesProc = numSitesProc
            heater%numSitesBlk = numSitesBlk

            ! Deallocate initial site information and print sites located on current process
            deallocate (heater%xSiteInit, heater%ySiteInit, heater%zSiteInit, heater%radiusInit)
#ifdef DEBUG_ALL
            print *, "Number of sites and block on process: ", heater%numSitesProc, blockCount
#endif
         end do
      else
         if (htr_meshMe .eq. MASTER_PE) print *, "[Heater_mapSitesToProc] Doing nothing, initial == .FALSE."
      end if

   else
      call Driver_abort("[Heater_mapSitesToProc] Unknown option. Supported options gridChanged and intial")
   end if

   call Timers_stop("Heater_mapSitesToProc")

end subroutine Heater_mapSitesToProc
