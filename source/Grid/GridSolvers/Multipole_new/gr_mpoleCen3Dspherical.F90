!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCen3Dspherical
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
!!  gr_mpoleCen3Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen3Dspherical (integer, intent (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 3D spherical
!!  geometry. It computes the center of expansion for the multipoles for 3D spherical
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                            integral (r * rho * rho  dr)
!!              Cen (x,y,z) = ----------------------------
!!                              integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                     sum cells (cell center r * cell mass * cell rho)
!!       Cen (x,y,z) = ------------------------------------------------
!!                             sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (x,y,z) has been determined, it is placed on the
!!  the nearest cell corner. The following is computed here:
!!
!!                  1) multipole expansion center (placed on nearest cell corner)
!!                  2) total mass (aborts, if <= 0)
!!                  3) the 'atomic' inner zone length (and its inverse)
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!! NOTES
!!
!!  gr_mpoleXcenter, gr_mpoleYcenter and gr_mpoleZcenter denote the location of
!!  the center of multipole expansion in the 3D cartesian framework (going from
!!  R,theta,phi -> x,y,z).
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen3Dspherical (idensvar)

  use Grid_data,         ONLY : gr_meshMe,   &
                                gr_meshComm

  use Driver_interface,  ONLY : Driver_abort

  use Grid_interface,    ONLY : Grid_getTileIterator,   &
                                Grid_releaseTileIterator,&
                                Grid_getCellCoords


  use gr_mpoleData,      ONLY : gr_mpoleDomainRmax,     &
                                gr_mpoleDomainRmin,     &
                                gr_mpoleDomainThetaMax, &
                                gr_mpoleDomainPhiMax,   &
                                gr_mpoleDrInnerZone,    &
                                gr_mpoleDrInnerZoneInv, &
                                gr_mpolePi,             & 
                                gr_mpoleHalfPi,         & 
                                gr_mpoleTwoPi,          & 
                                gr_mpoleXcenter,        &
                                gr_mpoleYcenter,        &
                                gr_mpoleZcenter,        &
                                gr_mpoleRcenter,        &
                                gr_mpoleThetaCenter,    &
                                gr_mpolePhiCenter,      &
                                gr_mpoleTotalMass
 
 use Grid_tile,         ONLY : Grid_tile_t
 use Grid_iterator,     ONLY : Grid_iterator_t
  
  !implicit none
  
#include "Flashx_mpi_implicitNone.fh"
#include "Simulation.h"
#include "constants.h"
#include "gr_mpole.h"

  
  integer, intent (in) :: idensvar

  logical :: domainRmax, domainThetaMax, domainPhiMax
  logical :: guardCells
  logical :: insideBlock
  logical :: invokeRecv
  logical :: positiveYaxis, negativeYaxis
  logical :: quadrantI, quadrantII, quadrantIII, quadrantIV

  integer :: blockNr
  integer :: blockID
  integer :: error
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: k,kmin,kmax
  integer :: maxEdges
  integer :: messageTag
  integer :: nCellsI, nCellsJ, nCellsK
  integer :: nEdgesI, nEdgesJ, nEdgesK

  integer :: locate      (1:3)
  integer :: status      (MPI_STATUS_SIZE)
  integer :: tileLimits   (LOW:HIGH,1:MDIM)

  real    :: alpha_t, beta_t, alpha_p, beta_p
  real    :: bndBoxILow
  real    :: bndBoxJLow
  real    :: bndBoxKLow
  real    :: cellDelta, cellDensity, cellMass, cellMassDensity, cellVolume
  real    :: cmRsph, cmPhi, cmTheta
  real    :: DeltaI, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf, DeltaJSine, DeltaJHalfSine
  real    :: DeltaK, DeltaKHalf, DeltaKSine, DeltaKHalfSine
  real    :: localMsum, localMDsum, localMDXsum, localMDYsum, localMDZsum
  real    :: maxRsph, maxTheta, maxPhi
  real    :: minRsph, minTheta, minPhi
  real    :: phi, phiCosine, phiSine, phiSineSave
  real    :: theta, thetaCosine, thetaSine, thetaSineSave
  real    :: Rsph
  real    :: totalMassDensityInv
  real    :: x,y,z
  real    :: angularVolumePart
  real, parameter :: sixth = 1./6.

  real    :: delta     (1:MDIM)
  real    :: localData (1:7)
  real    :: totalData (1:5)
  real    :: bndBox    (LOW:HIGH,1:MDIM)

  real, allocatable :: shifts   (:,:)
  real, pointer     :: solnData (:,:,:,:)
  
  integer :: lev
  type(Grid_tile_t)     :: tileDesc
  type(Grid_iterator_t) :: itor
  
  NULLIFY(solnData)
!
!
!     ...Sum quantities over all locally held leaf blocks.
!
!
  localMsum   = ZERO
  localMDsum  = ZERO
  localMDXsum = ZERO
  localMDYsum = ZERO
  localMDZsum = ZERO

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     tileLimits=tileDesc%limits

     call tileDesc%boundBox(bndBox)
     call tileDesc%deltas(delta)
     call tileDesc%getDataPtr(solnData, CENTER)

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
!
!
!             Since location of the center of multipole expansion will be determined in 3D cartesian (x,y,z) coordinates, conversion 
!             from 3D spherical (r, theta, phi) coordinates is necessary.
!
!             The cell volume is:
!
!                      (1/3) * (R^3 - r^3) * (cos[t] - cos[T]) * Dp
!
!             where
!
!                       r  =  left-most (smaller) radial cell distance
!                       R  =  right-most (larger) radial cell distance
!                       t  =  smaller cell angle
!                       T  =  larger cell angle
!                      Dp  =  phi delta value
!
!
!                       r  =  Rsph - Dr/2
!                       R  =  Rsph + Dr/2
!                       t  =  theta - Dt/2
!                       T  =  theta + Dt/2
!
!             with Dt and Dr being the cell's angular and radial delta values.
!             Hence the cell volume becomes, after using some trig identities:
!
!                (1/6) * Dp * Dr * sin[Dt/2] * (12 * (Rsph)^2 + Dr^2) * sin[theta]
!
!
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
           
           angularVolumePart = DeltaJHalfSine * thetaSine
        
           Rsph = bndBoxILow + DeltaIHalf
           do i = imin,imax

              cellVolume  = sixth * angularVolumePart * DeltaI * DeltaK * ((12. * Rsph * Rsph) + (DeltaI * DeltaI))
              cellDensity = solnData (idensvar,i,j,k)
              cellMass    = cellDensity * cellVolume
              cellMassDensity = cellMass * cellDensity

              x = Rsph * thetaSine * phiCosine
              y = Rsph * phiSine * thetaSine
              z = Rsph * thetaCosine

              localMsum   = localMsum   + cellMass
              localMDsum  = localMDsum  + cellMassDensity
              localMDXsum = localMDXsum + cellMassDensity * x
              localMDYsum = localMDYsum + cellMassDensity * y
              localMDZsum = localMDZsum + cellMassDensity * z

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
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()

  end do
call Grid_releaseTileIterator(itor)

!
!
!     ...Prepare for a one-time all reduce call.
!
!
  localData (1) = localMsum
  localData (2) = localMDsum
  localData (3) = localMDXsum
  localData (4) = localMDYsum
  localData (5) = localMDZsum
!
!
!     ...Calculate the total sums and give a copy to each processor.
!
!
  call  MPI_AllReduce (localData,   &
                       totalData,   &
                       5,           &
                       FLASH_REAL,  & 
                       MPI_Sum,     &
                       gr_meshComm, &
                       error        )
!
!
!     ...Analyze total mass obtained. If nonsense, abort.
!
!
  gr_mpoleTotalMass = totalData (1)
     
  if (abs (gr_mpoleTotalMass) < tiny (gr_mpoleTotalMass)) then
      call Driver_abort ('[gr_mpoleCen3Dspherical] ERROR:  gr_mpoleTotalMass <= 0')
  end if
!
!
!     ...Calculate center of multipole expansion cartesian coordinates.
!
!
  totalMassDensityInv = ONE / totalData (2)

  gr_mpoleXcenter = totalData (3) * totalMassDensityInv
  gr_mpoleYcenter = totalData (4) * totalMassDensityInv
  gr_mpoleZcenter = totalData (5) * totalMassDensityInv
!
!
!     ...We need to calculate the 3D spherical coordinates of the center
!        of multipole expansion in order to determine to which block it belongs.
!        FLASH convention is such that the phi angle of the 3D spherical
!        geometry ranges from 0 to 2pi. For accurate angle determination
!        we use the atan function, which returns angles in the range of
!        -pi/2 < phi < pi/2. The following inverse angle calculations are
!        thus performed, depending on which quadrant in the xy-plane the center
!        of multipole expansion lays:
!
!
!                                          y
!                                          |
!                         Quadrant II      |     Quadrant I
!                                          |
!                   phi = pi + atan (y/x)  |  phi = atan (y/x)
!                                          |  
!                                          | /
!                                          |/phi
!                       ----------------------------------------> x
!                                          |
!                                          |
!                         Quadrant III     |     Quadrant I
!                                          |
!                   phi = pi + atan (y/x)  |  phi = 2pi + atan (y/x)
!                                          |
!                                          |
!                                          |
!
!
!        If the center of multipole expansion lays exactly on the y-axis (x = 0),
!        we handle those cases separately.
!
!
 
  cmRsph = sqrt (gr_mpoleXcenter * gr_mpoleXcenter + gr_mpoleYcenter * gr_mpoleYcenter + gr_mpoleZcenter * gr_mpoleZcenter)
  cmTheta = acos (gr_mpoleZcenter / cmRsph) !Theta ranges from 0 to pi in FLASH, so this is OK
  
  quadrantI     = (gr_mpoleXcenter >  ZERO) .and. (gr_mpoleYcenter >= ZERO)
  quadrantII    = (gr_mpoleXcenter <  ZERO) .and. (gr_mpoleYcenter >= ZERO)
  quadrantIII   = (gr_mpoleXcenter <  ZERO) .and. (gr_mpoleYcenter <  ZERO)
  quadrantIV    = (gr_mpoleXcenter >  ZERO) .and. (gr_mpoleYcenter <  ZERO)
  positiveYaxis = (gr_mpoleXcenter == ZERO) .and. (gr_mpoleYcenter >  ZERO)
  negativeYaxis = (gr_mpoleXcenter == ZERO) .and. (gr_mpoleYcenter <  ZERO)

  if (quadrantI) then
      cmPhi = atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (quadrantII) then
      cmPhi = gr_mpolePi + atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (quadrantIII) then
      cmPhi = gr_mpolePi + atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (quadrantIV) then
      cmPhi = gr_mpoleTwoPi + atan (gr_mpoleYcenter / gr_mpoleXcenter)
  else if (positiveYaxis) then
      cmPhi = gr_mpoleHalfPi
  else if (negativeYaxis) then
      cmPhi = gr_mpolePi + gr_mpoleHalfPi
  else
      cmPhi = ZERO
  end if
  
  gr_mpoleRcenter   = cmRsph
  gr_mpoleThetaCenter = cmTheta
  gr_mpolePhiCenter = cmPhi
!
!
!     ...Find the local blockID to which the center of multipole expansion
!        belongs and place the center on the nearest cell corner. Also at
!        this point we determine the inner zone atomic length, since the
!        inner zone is defined around the center of multipole expansion.
!        Whatever processor is doing the relevant calculation sends its
!        final data (updated center of multipole expansion and inner zone
!        atomic length) to the master, which then broadcasts the info.
!
!
  messageTag = 1
  invokeRecv = .true.

!This loop hangs if r_min > 0, regardless if ignoring inner zone or not (probably due to missing region?)
if (gr_mpoleDomainRmin == ZERO) then
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call tileDesc%boundBox (bndBox)

     minRsph  = bndBox (LOW ,IAXIS)
     maxRsph  = bndBox (HIGH,IAXIS)
     minTheta = bndBox (LOW ,JAXIS)
     maxTheta = bndBox (HIGH,JAXIS)
     minPhi   = bndBox (LOW ,KAXIS)
     maxPhi   = bndBox (HIGH,KAXIS)

     insideBlock  =       (cmRsph  >= minRsph ) &
                    .and. (cmTheta >= minTheta) &
                    .and. (cmPhi   >= minPhi  ) &
                    .and. (cmRsph  <  maxRsph ) &     ! the < instead of <= is
                    .and. (cmTheta <  maxTheta) &     ! necessary for finding
                    .and. (cmPhi   <  maxPhi  )       ! the unique block

     domainRmax   =       (cmRsph  == maxRsph                 ) &    ! include (however unlikely) the
                    .and. (cmRsph  == gr_mpoleDomainRmax      )      ! missing R part of the domain
     domainThetaMax   =   (cmTheta == maxTheta                ) &    ! include (however unlikely) the
                    .and. (cmTheta == gr_mpoleDomainThetaMax  )      ! missing theta part of the domain
     domainPhiMax =       (cmPhi   == maxPhi                  ) &    ! include (however unlikely)) the
                    .and. (cmPhi   == gr_mpoleDomainPhiMax    )      ! missing phi part of the domain

     insideBlock  = insideBlock .or. domainRmax .or. domainThetaMax .or. domainPhiMax

     if (insideBlock) then

         lev=tileDesc%level
         call tileDesc%deltas(delta)
         tileLimits=tileDesc%limits

         DeltaI = delta (IAXIS)
         DeltaJ = delta (JAXIS)
         DeltaK = delta (KAXIS)

         gr_mpoleDrInnerZone = HALF * DeltaI   ! based on Rsph only

         imin = tileLimits (LOW, IAXIS)
         jmin = tileLimits (LOW, JAXIS)
         kmin = tileLimits (LOW, KAXIS)  
         imax = tileLimits (HIGH,IAXIS)
         jmax = tileLimits (HIGH,JAXIS)
         kmax = tileLimits (HIGH,KAXIS)

         nCellsI = imax - imin + 1
         nCellsJ = jmax - jmin + 1
         nCellsK = kmax - kmin + 1

         nEdgesI = nCellsI + 1
         nEdgesJ = nCellsJ + 1
         nEdgesK = nCellsK + 1

         maxEdges = max (nEdgesI, nEdgesJ, nEdgesK)

         allocate (shifts (1:maxEdges,3))

         !guardCells = .false.

         call Grid_getCellCoords (IAXIS, FACES, lev, tileLimits(LOW,:), tileLimits(HIGH,:), shifts (1:nEdgesI, 1))
         call Grid_getCellCoords (JAXIS, FACES, lev, tileLimits(LOW,:), tileLimits(HIGH,:), shifts (1:nEdgesJ, 2))
         call Grid_getCellCoords (KAXIS, FACES, lev, tileLimits(LOW,:), tileLimits(HIGH,:), shifts (1:nEdgesK, 3))

         shifts (1:nEdgesI,1) = shifts (1:nEdgesI,1) - cmRsph
         shifts (1:nEdgesJ,2) = shifts (1:nEdgesJ,2) - cmTheta
         shifts (1:nEdgesK,3) = shifts (1:nEdgesK,3) - cmPhi

         locate (1) = minloc (abs (shifts (1:nEdgesI,1)), dim = 1)
         locate (2) = minloc (abs (shifts (1:nEdgesJ,2)), dim = 1)
         locate (3) = minloc (abs (shifts (1:nEdgesK,3)), dim = 1)

         cmRsph  = cmRsph  + shifts (locate (1),1)  ! move multipole center to nearest R edge
         cmTheta = cmTheta + shifts (locate (2),2)  ! move multipole center to nearest Z edge
         cmPhi   = cmPhi   + shifts (locate (3),3)  ! move multipole center to nearest Phi edge

         deallocate (shifts)

         gr_mpoleXcenter   = cmRsph * sin (cmTheta) * cos (cmPhi)    !
         gr_mpoleYcenter   = cmRsph * sin (cmPhi) * sin (cmTheta)    ! convert back to 3D cartesian
         gr_mpoleZcenter   = cmRsph * cos (cmTheta)                  !
         gr_mpoleRcenter   = cmRsph
         gr_mpoleThetaCenter = cmTheta
         gr_mpolePhiCenter = cmPhi
         

         localData (1) = gr_mpoleDrInnerZone
         localData (2) = gr_mpoleXcenter
         localData (3) = gr_mpoleYcenter
         localData (4) = gr_mpoleZcenter
         localData (5) = gr_mpoleRcenter
         localData (6) = gr_mpoleThetaCenter
         localData (7) = gr_mpolePhiCenter

         if (gr_meshMe /= MASTER_PE) then

             call MPI_Send (localData,    &
                            7,            &
                            FLASH_REAL,   &
                            MASTER_PE,    &
                            messageTag,   &
                            gr_meshComm,  &
                            error         )
         else
             invokeRecv = .false.
         end if

         exit

     end if
      call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  if ((gr_meshMe == MASTER_PE) .and. invokeRecv) then

       call MPI_Recv (localData,      &
                      7,              &
                      FLASH_REAL,     &
                      MPI_ANY_SOURCE, &
                      messageTag,     &
                      gr_meshComm,    &
                      status,         &
                      error           )
  end if
else

  gr_mpoleDrInnerZone = ZERO
  localData (1)       = gr_mpoleDrInnerZone
  localData (2)       = gr_mpoleXcenter
  localData (3)       = gr_mpoleYcenter
  localData (4)       = gr_mpoleZcenter
  localData (5)       = gr_mpoleRcenter
  localData (6)       = gr_mpoleThetaCenter
  localData (7)       = gr_mpolePhiCenter


end if
!
!
!     ...At this point, the master has all the info. Broadcast and update all
!        other processors.
!
!
  call MPI_Bcast (localData,   &
                  7,           &
                  FLASH_REAL,  &
                  MASTER_PE,   &
                  gr_meshComm, &
                  error        )

  gr_mpoleDrInnerZone    = localData (1)
  gr_mpoleXcenter        = localData (2)
  gr_mpoleYcenter        = localData (3)
  gr_mpoleZcenter        = localData (4)
  gr_mpoleRcenter        = localData (5)
  gr_mpoleThetaCenter    = localData (6)
  gr_mpolePhiCenter      = localData (7)

  if (gr_mpoleDrInnerZone == ZERO) then
    gr_mpoleDrInnerZoneInv = ZERO
  else
    gr_mpoleDrInnerZoneInv = ONE / gr_mpoleDrInnerZone
  endif

!
!
!     ...Ready!
!
!
  return
end subroutine gr_mpoleCen3Dspherical
