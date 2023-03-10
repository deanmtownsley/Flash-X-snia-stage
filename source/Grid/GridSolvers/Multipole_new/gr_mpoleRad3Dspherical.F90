!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleRad3Dspherical
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
!!
!!  gr_mpoleRad3Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleRad3Dspherical ()
!!
!! DESCRIPTION
!!
!!  This routine determines the radial sampling for accumulating the moments
!!  for the three-dimensional (3D) spherical case.
!!
!!***

subroutine gr_mpoleRad3Dspherical ()

  use Driver_interface,  ONLY : Driver_abort

  use Grid_data,         ONLY : gr_meshMe,   &
                                gr_meshComm

  use Grid_interface,    ONLY : Grid_getMinCellSizes,   &
                                Grid_getLocalNumBlks,   &
                                Grid_getTileIterator,   &
                                Grid_releaseTileIterator

  use gr_mpoleInterface, ONLY : gr_mpoleSetInnerZoneGrid, &
                                gr_mpoleSetOuterZoneGrid

  use gr_mpoleData,      ONLY : gr_mpolePi,              &
                                gr_mpoleDr,              &
                                gr_mpoleDrInv,           &
                                gr_mpoleDrInnerZone,     &
                                gr_mpoleDrInnerZoneInv,  &
                                gr_mpoleMaxR,            &
                                gr_mpoleIgnoreInnerZone, &
                                gr_mpoleInnerZoneExists, &
                                gr_mpoleInnerZoneMaxR,   &
                                gr_mpoleInnerZoneQmax,   &
                                gr_mpoleInnerZoneSize,   &
                                gr_mpoleOuterZoneExists, &
                                gr_mpoleDomainRmax,      &
                                gr_mpoleDomainZmin,      &
                                gr_mpoleDomainZmax,      &
                                gr_mpoleXcenter,         &
                                gr_mpoleYcenter,         &
                                gr_mpoleZcenter
  
  use Grid_tile,         ONLY : Grid_tile_t
  use Grid_iterator,     ONLY : Grid_iterator_t


!  implicit none

#include "Flashx_mpi_implicitNone.fh"
#include "Simulation.h"
#include "constants.h"
#include "gr_mpole.h"


  integer :: nblks
  integer :: error
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: k,kmin,kmax
  integer :: nPinnerZone
  integer :: nRinnerZone
  integer :: nBlocal
  integer :: nPlocal
  integer :: nRlocal
  integer :: nRlocalPrev

  integer :: localData   (1:2)
  integer :: globalData  (1:2)
  integer :: tileLimits   (LOW:HIGH,1:MDIM)

  real    :: alpha_t, beta_t, alpha_p, beta_p
  real    :: bndBoxILow
  real    :: bndBoxJLow
  real    :: bndBoxKLow
  real    :: cmRsph, cmZ
  real    :: DeltaI, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf, DeltaJSine, DeltaJHalfSine
  real    :: DeltaK, DeltaKHalf, DeltaKSine, DeltaKHalfSine
  real    :: Dxmin
  real    :: maxR, maxRsqr
  real    :: phi, phiCosine, phiSine, phiSineSave
  real    :: theta, thetaCosine, thetaSine, thetaSineSave
  real    :: Rsph
  real    :: r, rsqr
  real    :: x, y, z

  real    :: delta        (1:MDIM)
  real    :: minCellSizes (1:MDIM)
  real    :: bndBox       (LOW:HIGH,1:MDIM)

  logical, allocatable :: blockListInnerZone (:)
  real,    allocatable :: RinnerZone         (:)

  type(Grid_tile_t)     :: tileDesc
  type(Grid_iterator_t) :: itor

!
!
!       ...Get the minimum cell sizes for the linear parts of the domain, which
!          determine the radial distance in 3D space. 
!
!
  call Grid_getMinCellSizes (minCellSizes)

  Dxmin = minCellSizes (IAXIS)
!
!
!       ...Determine the maximum distance from the center of multipole expansion to the
!          computational domain boundaries.

  cmRsph = sqrt (gr_mpoleXcenter * gr_mpoleXcenter + gr_mpoleYcenter * gr_mpoleYcenter + gr_mpoleZcenter * gr_mpoleZcenter)
  gr_mpoleMaxR  = gr_mpoleDomainRmax + cmRsph
  gr_mpoleDr    = HALF * Dxmin
  gr_mpoleDrInv = ONE / gr_mpoleDr
  
  
!
!
!     ...Set initial indicators for inner and outer zone.
!
!
  gr_mpoleInnerZoneExists = .not. gr_mpoleIgnoreInnerZone
  gr_mpoleOuterZoneExists = .not. gr_mpoleInnerZoneExists

  if (gr_mpoleInnerZoneExists) then
!
!
!     ...Proceed with establishing the inner zone (if any). From the determined
!        inner zone atomic length and the previously found maximal radial domain
!        distance, readjust the inner zone size variable. Two cases can happen:
!        1) the size of the inner zone fits into the complete radial doamin (no
!        adjustment needed) or 2) the size of the inner zone exceeds the complete
!        radial domain (adjustment needed). Also override existence criterion for
!        the outer zone, if the largest domain radius exceeds the inner zone radius.
!
!
      gr_mpoleOuterZoneExists = (gr_mpoleInnerZoneSize * gr_mpoleDrInnerZone < gr_mpoleMaxR)

      if ( gr_mpoleInnerZoneSize * gr_mpoleDrInnerZone > gr_mpoleMaxR ) then
           gr_mpoleInnerZoneSize = int (ceiling (gr_mpoleMaxR * gr_mpoleDrInnerZoneInv))
      end if
!
!
!     ...Determine the number of radii to be expected in the inner zone.
!        For each processor, store those local blockID's that actually
!        have radii in the inner zone.
!
!
      call Grid_getLocalNumBlks(nblks)

      allocate (blockListInnerZone (1:nblks))

      gr_mpoleInnerZoneMaxR = real (gr_mpoleInnerZoneSize) * gr_mpoleDrInnerZone
      maxRsqr               = gr_mpoleInnerZoneMaxR * gr_mpoleInnerZoneMaxR

      nBlocal = 0
      nRlocal = 0
      nRlocalPrev = 0

      call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
      do while(itor%isValid())
         call itor%currentTile(tileDesc)
         tileLimits=tileDesc%limits

         call tileDesc%boundBox(bndBox)
         call tileDesc%deltas(delta)
         
         imin = tileLimits (LOW, IAXIS)
         jmin = tileLimits (LOW, JAXIS)
         kmin = tileLimits (LOW, KAXIS)  
         imax = tileLimits (HIGH,IAXIS)
         jmax = tileLimits (HIGH,JAXIS)
         kmax = tileLimits (HIGH,KAXIS)

         DeltaI           = delta (IAXIS)
         DeltaIHalf       = DeltaI * HALF  
        
         DeltaJ           = delta (JAXIS)
         DeltaJHalf       = DeltaJ * HALF
         DeltaJSine       = sin (DeltaJ)
         DeltaJHalfSine   = sin (DeltaJ * HALF)
     
         DeltaK           = delta (KAXIS)
         DeltaKHalf       = DeltaJ * HALF
         DeltaKSine       = sin (DeltaK)
         DeltaKHalfSine   = sin (DeltaK * HALF)

         bndBoxILow = bndBox (LOW,IAXIS)
         bndBoxJLow = bndBox (LOW,JAXIS)
         bndBoxKLow = bndBox (LOW,KAXIS)

         alpha_p       = TWO * DeltaKHalfSine * DeltaKHalfSine
         beta_p        = DeltaKSine
     
         alpha_t       = TWO * DeltaJHalfSine * DeltaJHalfSine       
         beta_t        = DeltaJSine  

         phi       = bndBoxKLow + DeltaKHalf
         phiSine   = sin (phi)
         phiCosine = cos (phi)

         do k = kmin, kmax
            theta       = bndBoxJLow + DeltaJHalf             
            thetaSine   = sin (theta)
            thetaCosine = cos (theta) 
            do j = jmin, jmax
               Rsph = bndBoxILow + DeltaIHalf
               do i = imin, imax

                  x = Rsph * thetaSine * phiCosine - gr_mpoleXcenter
                  y = Rsph * phiSine * thetaSine - gr_mpoleYcenter
                  z = Rsph * thetaCosine - gr_mpoleZcenter

                  rsqr = x * x + y * y + z * z

                  if (rsqr <= maxRsqr) then
                      nRlocal = nRlocal + 1
                  end if

                  Rsph = Rsph + DeltaI
               end do
               thetaSineSave = thetaSine
               thetaSine     = thetaSine   - (alpha_t * thetaSine   - beta_t * thetaCosine  )
               thetaCosine   = thetaCosine - (alpha_t * thetaCosine + beta_t * thetaSineSave)
            end do
            phiSineSave = phiSine
            phiSine     = phiSine   - (alpha_p * phiSine   - beta_p * phiCosine  )
            phiCosine   = phiCosine - (alpha_p * phiCosine + beta_p * phiSineSave)
         end do

         nBlocal = nBlocal + 1
         blockListInnerZone (nBlocal) = (nRlocal > nRlocalPrev)
         nRlocalPrev = nRlocal

         call itor%next()

      end do
      call Grid_releaseTileIterator(itor)
!
!
!     ...Calculate the total number of processors contributing to the inner
!        zone radii and the overall total number of inner zone radii to be
!        expected. Allocate the array that will contain all inner zone radii
!        on all processors. If no inner zone radii are found globally, there
!        is something wrong and the program has to stop.
!
!
      nPlocal = min (nRlocal,1)  ! current processor adds +1, if inner zone radii found

      localData (1) = nPlocal
      localData (2) = nRlocal

      call MPI_AllReduce (localData,     &
                          globalData,    &
                          2,             &
                          FLASH_INTEGER, &
                          MPI_SUM,       &
                          gr_meshComm,   &
                          error          )

      nPinnerZone = globalData (1)
      nRinnerZone = globalData (2)

      if (nRinnerZone == 0) then
          call Driver_abort ('[gr_mpoleRad3Dspherical] ERROR: no inner zone radii found')
      end if

      allocate (RinnerZone (1:nRinnerZone))
!
!
!     ...Calculate and store now all inner zone radii on each processor.
!        Loop only over those local blocks which actually contribute to the
!        inner zone (skip, if no blocks).
!
!
      nRlocal = 0
      nBlocal = 0

      call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
      do while(itor%isValid())
         nBlocal=nBlocal+1
         if(blockListInnerZone(nBlocal)) then

            call itor%currentTile(tileDesc)
            tileLimits=tileDesc%limits

            call tileDesc%boundBox(bndBox)
            call tileDesc%deltas(delta)

             imin = tileLimits (LOW, IAXIS)
             jmin = tileLimits (LOW, JAXIS)
             kmin = tileLimits (LOW, KAXIS)  
             imax = tileLimits (HIGH,IAXIS)
             jmax = tileLimits (HIGH,JAXIS)
             kmax = tileLimits (HIGH,KAXIS)

             DeltaI           = delta (IAXIS)
             DeltaIHalf       = DeltaI * HALF  
        
             DeltaJ           = delta (JAXIS)
             DeltaJHalf       = DeltaJ * HALF
             DeltaJSine       = sin (DeltaJ)
             DeltaJHalfSine   = sin (DeltaJ * HALF)
     
             DeltaK           = delta (KAXIS)
             DeltaKHalf       = DeltaK * HALF
             DeltaKSine       = sin (DeltaK)
             DeltaKHalfSine   = sin (DeltaK * HALF)

             bndBoxILow = bndBox (LOW,IAXIS)
             bndBoxJLow = bndBox (LOW,JAXIS)
             bndBoxKLow = bndBox (LOW,KAXIS)

             alpha_p       = TWO * DeltaKHalfSine * DeltaKHalfSine
             beta_p        = DeltaKSine
     
             alpha_t       = TWO * DeltaJHalfSine * DeltaJHalfSine       
             beta_t        = DeltaJSine  

             phi       = bndBoxKLow + DeltaKHalf
             phiSine   = sin (phi)
             phiCosine = cos (phi)

             do k = kmin, kmax
                theta       = bndBoxJLow + DeltaJHalf             
                thetaSine   = sin (theta)
                thetaCosine = cos (theta)
                do j = jmin, jmax
                   Rsph = bndBoxILow + DeltaIHalf
                   do i = imin, imax

                      x = Rsph * thetaSine * phiCosine - gr_mpoleXcenter
                      y = Rsph * phiSine * thetaSine - gr_mpoleYcenter
                      z = Rsph * thetaCosine - gr_mpoleZcenter

                      r = sqrt (x * x + y * y + z * z) 

                      if (r <= gr_mpoleInnerZoneMaxR) then
                          nRlocal = nRlocal + 1
                          RinnerZone (nRlocal) = r
                      end if

                      Rsph = Rsph + DeltaI
                   end do
                   thetaSineSave = thetaSine
                   thetaSine     = thetaSine   - (alpha_t * thetaSine   - beta_t * thetaCosine  )
                   thetaCosine   = thetaCosine - (alpha_t * thetaCosine + beta_t * thetaSineSave)
                end do
                phiSineSave = phiSine
                phiSine     = phiSine   - (alpha_p * phiSine   - beta_p * phiCosine  )
                phiCosine   = phiCosine - (alpha_p * phiCosine + beta_p * phiSineSave)
             end do

          end if
          call itor%next()

      end do
      call Grid_releaseTileIterator(itor)

      deallocate (blockListInnerZone)
!
!
!       ...Set up the inner zone radial grid.
!
!
      call gr_mpoleSetInnerZoneGrid (nRlocal,     &
                                     nRinnerZone, &
                                     nPinnerZone, &
                                     RinnerZone   )

      deallocate (RinnerZone)

  else
!
!
!       ...No inner zone! Set the inner zone variables to nonexistent.
!
!
      gr_mpoleDrInnerZone   = ZERO
      gr_mpoleInnerZoneMaxR = ZERO
      gr_mpoleInnerZoneQmax = 0

  end if  ! inner zone condition
!
!
!       ...Complete the radial grid picture by setting up the outer (statistical)
!          zone radial grid.
!
!
  call gr_mpoleSetOuterZoneGrid ()
!
!
!       ...Ready!
!
!
  return
end subroutine gr_mpoleRad3Dspherical
