!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMom3Dspherical
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
!!  gr_mpoleMom3Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom3Dspherical (integer (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Prepares for evaluation of the moments in 3D spherical geometry. In this
!!  routine, all the necessary arrays are prepared to enable evaluation of
!!  the moments in radial bin order. Each of the moments are grouped together
!!  according to their radial bins. This will ensure optimum unit stride values
!!  when accessing the big moment arrays and makes threading trivial.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

!!REORDER(4): Uin

subroutine gr_mpoleMom3Dspherical (idensvar)

  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_tile,         ONLY : Grid_tile_t
  use Grid_iterator,     ONLY : Grid_iterator_t

  use gr_mpoleInterface, ONLY : gr_mpoleMomBins3Dspherical

  use gr_mpoleData,      ONLY : gr_mpoleDrInv,                  &
                                gr_mpoleDrInnerZoneInv,         &
                                gr_mpoleMaxQ,                   &
                                gr_mpoleMaxRadialZones,         &
                                gr_mpoleMinRadialZone,          &
                                gr_mpoleZoneRmax,               &
                                gr_mpoleZoneQmax,               &
                                gr_mpoleZoneType,               &
                                gr_mpoleZoneScalarInv,          &
                                gr_mpoleZoneLogNormInv,         &
                                gr_mpoleZoneExponentInv,        &
                                gr_mpoleInnerZoneMaxR,          &
                                gr_mpoleInnerZoneDrRadii,       &
                                gr_mpoleInnerZoneQlower,        &
                                gr_mpoleInnerZoneQupper,        &
                                gr_mpoleInnerZoneResolution,    &
                                gr_mpoleInnerZoneResolutionInv, &
                                gr_mpoleOuterZoneQshift,        &
                                gr_mpoleQ,                      &
                                gr_mpoleQused,                  &
                                gr_mpoleQnumberOfCells,         &
                                gr_mpoleQdataCells3D

  use gr_mpoleData,      ONLY : gr_mpoleXcenter,                &
                                gr_mpoleYcenter,                &
                                gr_mpoleZcenter

  implicit none
  
#include "Simulation.h"
#include "constants.h"
#include "gr_mpole.h"

  integer, intent (in) :: idensvar

  logical :: innerZonePotential

  integer :: blockNr, blockID
  integer :: DrUnit
  integer :: i, imin, imax
  integer :: j, jmin, jmax
  integer :: k, kmin, kmax
  integer :: maxCells
  integer :: nC, nQ
  integer :: Q, Qlocal, Qlower, Qupper
  integer :: type
  integer :: used
  integer :: zone

  integer, save :: maxQtype                ! for multithreading needs to be on stack (save)

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: alpha_t, beta_t, alpha_p, beta_p
  real    :: bndBoxILow, bndBoxJLow, bndBoxKLow
  real    :: cellDelta, cellDensity, cellMass, cellVolume
  real    :: DeltaI, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf, DeltaJSine, DeltaJHalfSine
  real    :: DeltaK, DeltaKHalf, DeltaKSine, DeltaKHalfSine
  real    :: theta, thetaCosine, thetaSine, thetaSineSave
  real    :: phi, phiCosine, phiSine, phiSineSave
  real    :: r, rlocal, rinDrs
  real    :: Rsph
  real    :: sclInv, lgnInv, expInv
  real    :: x,y,z
  real    :: angularVolumePart

  real    :: delta           (1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)
  
  real, parameter :: sixth = 1./6.

  type(Grid_iterator_t) :: itor
  real,              pointer    :: Uin(:,:,:,:)
  type(Grid_tile_t)     :: tileDesc
  integer :: level
  integer, dimension(LOW:HIGH,MDIM) :: tileLimits   ,tileLimitsGC, grownLimits

!
!
!     ...The first pass over all blocks on the current processor will get us information
!        about how many different radial bin indices will be addressed and for each such
!        radial bin index, how many cells it will contain.
!
!
!$omp single
  gr_mpoleQused (:) = 0 
  nullify(Uin)

   call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
   do while(itor%isValid())
      call itor%currentTile(tileDesc)
      tileLimits(:,:)=tileDesc%limits
      tileLimitsGC(:,:)=tileDesc%blkLimitsGC
      grownLimits(:,:)=tileDesc%grownLimits
      call tileDesc%deltas(delta)
      level=tileDesc%level
      call tileDesc%getDataPtr(Uin, CENTER)
      call tileDesc%boundBox(bndBox)

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
!             Since the location of the center of mass has been determined
!             in 3D cartesian (x,y,z) coordinates, converting from 3D spherical (r, theta, phi) coordinates
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
       Rsph = bndBoxILow + DeltaIHalf
       do i = imin, imax

        x = Rsph * thetaSine * phiCosine - gr_mpoleXcenter
        y = Rsph * phiSine * thetaSine - gr_mpoleYcenter
        z = Rsph * thetaCosine - gr_mpoleZcenter

        r = sqrt (x * x + y * y + z * z)  
!
!
!        ...Find the radial bin and increment radial bin counter.
!
!
        innerZonePotential = r <= gr_mpoleInnerZoneMaxR

        if (innerZonePotential) then

            rinDrs = r * gr_mpoleDrInnerZoneInv
            DrUnit = int (ceiling (rinDrs))
            Qlower = gr_mpoleInnerZoneQlower (DrUnit)
            Qupper = gr_mpoleInnerZoneQupper (DrUnit)

            do Q = Qlower,Qupper
               if (rinDrs <= gr_mpoleInnerZoneDrRadii (Q)) exit
            end do

        else

            do zone = gr_mpoleMinRadialZone, gr_mpoleMaxRadialZones
               if (r - gr_mpoleZoneRmax (zone) <= ZERO) exit
            end do

            rlocal = r - gr_mpoleZoneRmax    (zone - 1)
            type   = gr_mpoleZoneType        (zone)
            sclInv = gr_mpoleZoneScalarInv   (zone)
            expInv = gr_mpoleZoneExponentInv (zone)

            if (type == ZONE_EXPONENTIAL) then
                Qlocal = ceiling ( (rlocal * sclInv * gr_mpoleDrInv) ** expInv )
            else if (type == ZONE_LOGARITHMIC) then
                lgnInv = gr_mpoleZoneLogNormInv (zone)
                Qlocal = ceiling ( expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE) )
            end if

            Q = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

        end if

        gr_mpoleQused (Q) = gr_mpoleQused (Q) + 1

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
     call tileDesc%releaseDataPtr(Uin,CENTER)
     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)
!
!
!     ...Create the arrays that will contain the radial info.
!
!
  maxQtype = count  (gr_mpoleQused /= 0)
  maxCells = maxval (gr_mpoleQused     )

  allocate (gr_mpoleQ              (1:maxQtype))
  allocate (gr_mpoleQnumberOfCells (1:maxQtype))
  allocate (gr_mpoleQdataCells3D   (1:maxCells , 1:maxQtype))
!
!
!     ...The second pass over all blocks on the current processor will scatter all
!        the radial bin information into the radial bin info array.
!
!
  gr_mpoleQused (:) = 0 

  nQ = 0

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     tileLimits(:,:)=tileDesc%limits
     tileLimitsGC(:,:)=tileDesc%blkLimitsGC
     grownLimits(:,:)=tileDesc%grownLimits
     call tileDesc%deltas(delta)
     level=tileDesc%level
     call tileDesc%getDataPtr(Uin, CENTER)
     call tileDesc%boundBox(bndBox)

     imin = tileLimits (LOW, IAXIS)
     jmin = tileLimits (LOW, JAXIS)
     kmin = tileLimits (LOW, KAXIS)  
     imax = tileLimits (HIGH,IAXIS)
     jmax = tileLimits (HIGH,JAXIS)
     kmax = tileLimits (HIGH,KAXIS)

     DeltaI         = delta (IAXIS)
     DeltaIHalf     = DeltaI * HALF
     
     DeltaJ         = delta (JAXIS)
     DeltaJHalf     = DeltaJ * HALF
     DeltaJSine     = sin (DeltaJ)
     DeltaJHalfSine = sin (DeltaJHalf)
     
     DeltaK           = delta (KAXIS)
     DeltaKHalf       = DeltaK * HALF
     DeltaKSine       = sin (DeltaK)
     DeltaKHalfSine   = sin (DeltaK * HALF)
     

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
     bndBoxKLow = bndBox (LOW,KAXIS)
!
!
!          ...Create all the cell info needed and place into proper radial bin array places.
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
       do i = imin, imax

       cellVolume  = sixth * angularVolumePart * DeltaI * DeltaK * ((12. * Rsph * Rsph) + (DeltaI * DeltaI))
       cellDensity = Uin (idensvar,i,j,k)
       cellMass    = cellDensity * cellVolume

       x = Rsph * thetaSine * phiCosine - gr_mpoleXcenter
       y = Rsph * phiSine * thetaSine - gr_mpoleYcenter
       z = Rsph * thetaCosine - gr_mpoleZcenter

       r = sqrt (x * x + y * y + z * z) 
!
!
!        ...Find the radial bin.
!
!
        innerZonePotential = r <= gr_mpoleInnerZoneMaxR

        if (innerZonePotential) then

            rinDrs = r * gr_mpoleDrInnerZoneInv
            DrUnit = int (ceiling (rinDrs))
            Qlower = gr_mpoleInnerZoneQlower (DrUnit)
            Qupper = gr_mpoleInnerZoneQupper (DrUnit)

            do Q = Qlower,Qupper
               if (rinDrs <= gr_mpoleInnerZoneDrRadii (Q)) exit
            end do

        else

            do zone = gr_mpoleMinRadialZone, gr_mpoleMaxRadialZones
               if (r - gr_mpoleZoneRmax (zone) <= ZERO) exit
            end do

            rlocal = r - gr_mpoleZoneRmax    (zone - 1)
            type   = gr_mpoleZoneType        (zone)
            sclInv = gr_mpoleZoneScalarInv   (zone)
            expInv = gr_mpoleZoneExponentInv (zone)

            if (type == ZONE_EXPONENTIAL) then
                Qlocal = ceiling ( (rlocal * sclInv * gr_mpoleDrInv) ** expInv )
            else if (type == ZONE_LOGARITHMIC) then
                lgnInv = gr_mpoleZoneLogNormInv (zone)
                Qlocal = ceiling ( expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE) )
            end if

            Q = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

        end if

        used = gr_mpoleQused (Q)

        if (used == 0) then

            nQ = nQ + 1

            gr_mpoleQused                (Q)             = nQ
            gr_mpoleQ                   (nQ)             = Q
            gr_mpoleQnumberOfCells      (nQ)             = 1
            gr_mpoleQdataCells3D      (1,nQ) % coord1    = x
            gr_mpoleQdataCells3D      (1,nQ) % coord2    = y
            gr_mpoleQdataCells3D      (1,nQ) % coord3    = z
            gr_mpoleQdataCells3D      (1,nQ) % cellMass  = cellMass
            gr_mpoleQdataCells3D      (1,nQ) % radius    = r

        else

            nC = gr_mpoleQnumberOfCells (used) + 1

            gr_mpoleQnumberOfCells    (used)             = nC
            gr_mpoleQdataCells3D   (nC,used) % coord1    = x
            gr_mpoleQdataCells3D   (nC,used) % coord2    = y
            gr_mpoleQdataCells3D   (nC,used) % coord3    = z
            gr_mpoleQdataCells3D   (nC,used) % cellMass  = cellMass
            gr_mpoleQdataCells3D   (nC,used) % radius    = r

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

     call tileDesc%releaseDataPtr(Uin,CENTER)
     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)
!$omp end single
!
!
!    ...Call the radial bin clustered moment evaluation routine (all threads).
!
!
  call gr_mpoleMomBins3Dspherical (maxQtype)
!
!
!    ...Deallocate used arrays.
!
!
!$omp single
  deallocate (gr_mpoleQ             )
  deallocate (gr_mpoleQnumberOfCells)
  deallocate (gr_mpoleQdataCells3D  )
!$omp end single
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleMom3Dspherical
