!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMom2Dspherical
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
!!  gr_mpoleMom2Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom2Dspherical (integer (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Prepares for evaluation of the moments in 2D spherical geometry. In this
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

!!REORDER(4): solnData

subroutine gr_mpoleMom2Dspherical (idensvar)

  use Grid_interface,    ONLY : Grid_getTileIterator,           &
                                Grid_releaseTileIterator

  use gr_mpoleInterface, ONLY : gr_mpoleMomBins2Dspherical

  use gr_mpoleData,      ONLY : gr_mpoleThirdPi,                &
                                gr_mpoleFourPi,                 &
                                gr_mpoleDrInv,                  &
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
                                gr_mpoleQdataCells2D

  use gr_mpoleData,      ONLY : gr_mpoleZcenter

  use Grid_tile,         ONLY : Grid_tile_t
  use Grid_iterator,     ONLY : Grid_iterator_t

  implicit none
  
#include "Simulation.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: idensvar

  logical :: innerZonePotential

  integer :: DrUnit
  integer :: i,imin,imax
  integer :: j,jmin,jmax
  integer :: maxCells
  integer :: nC, nQ
  integer :: Q, Qlocal, Qlower, Qupper
  integer :: type
  integer :: used
  integer :: zone

  integer, save :: maxQtype                ! for multithreading needs to be on stack (save)

  integer :: tileLimits   (LOW:HIGH,1:MDIM)

  real    :: alpha, beta
  real    :: angularVolumePart
  real    :: bndBoxILow, bndBoxJLow
  real    :: cellDensity, cellMass, cellVolume
  real    :: DeltaI, DeltaIcube, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf, DeltaJSine, DeltaJHalfSine
  real    :: r, rlocal, rinDrs
  real    :: Rsph
  real    :: sclInv, lgnInv, expInv
  real    :: theta, thetaSine, thetaSineSave, thetaCosine
  real    :: x,z

  real    :: delta  (1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real, pointer :: solnData (:,:,:,:)

  type(Grid_tile_t)     :: tileDesc
  type(Grid_iterator_t) :: itor

  nullify(solnData)

!
!     ...The first pass over all tiles on the current processor will get us information
!        about how many different radial bin indices will be addressed and for each such
!        radial bin index, how many cells it will contain.
!
!
!$omp single
  gr_mpoleQused (:) = 0 

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     tileLimits=tileDesc%limits

     call tileDesc%boundBox(bndBox)
     call tileDesc%deltas(delta)

     imin  = tileLimits (LOW, IAXIS)
     jmin  = tileLimits (LOW, JAXIS)
     imax  = tileLimits (HIGH,IAXIS)
     jmax  = tileLimits (HIGH,JAXIS)

     DeltaI         = delta (IAXIS)
     DeltaIHalf     = DeltaI * HALF
     DeltaJ         = delta (JAXIS)
     DeltaJHalf     = DeltaJ * HALF
     DeltaJSine     = sin (DeltaJ)
     DeltaJHalfSine = sin (DeltaJHalf)

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
!
!
!          ...The 2D spherical case:
!
!
!
!                         |*
!                         |     *
!                         |theta/  *
!                         |    /     *
!                         |   /       *
!                         |  /         *              Rsph --> stored in i-index
!                         | / Rsph      *            theta --> stored in j-index
!                         |/            *
!                         |             *
!                         |            *
!                         |           *
!                         |          *
!                         |        *
!                         |     *
!                         |*
!
!
!             If we would calculate the M=0 Moments around the domain origin, the
!             radial part of the M=0 recursion is equal to Rsph and the z-axis part is
!             Rsph * sin [theta]. However, our Moments will be evaluated around the
!             center of mass, whose coordinates were evaluated in 'cartesian' fashion.
!             Hence we must first convert the original (Rsph,theta) pair into cartesian
!             form and determine the new cartesian pair (x,z) based upon the location of
!             the center of mass. The new radial part for the Moment recursion is then
!             sqrt(x^2+z^2) and the z-axis component is the 'z' from the new cartesian
!             pair (x,z).
!
!
     alpha       = TWO * DeltaJHalfSine * DeltaJHalfSine
     beta        = DeltaJSine

     theta       = bndBoxJLow + DeltaJHalf
     thetaSine   = sin (theta)
     thetaCosine = cos (theta)

     do j = jmin,jmax
      Rsph = bndBoxILow + DeltaIHalf
      do i = imin,imax

       x = Rsph * thetaSine
       z = Rsph * thetaCosine - gr_mpoleZcenter

       r = sqrt (x * x + z * z)
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
      thetaSine     = thetaSine   - (alpha * thetaSine   - beta * thetaCosine  )
      thetaCosine   = thetaCosine - (alpha * thetaCosine + beta * thetaSineSave)
     end do
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
  allocate (gr_mpoleQdataCells2D   (1:maxCells , 1:maxQtype))
!
!
!     ...The second pass over all tiles on the current processor will scatter all
!        the radial bin information into the radial bin info array.
!
!
  gr_mpoleQused (:) = 0 

  nQ = 0

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     tileLimits=tileDesc%limits

     call tileDesc%boundBox(bndBox)
     call tileDesc%deltas(delta)
     call tileDesc%getDataPtr(solnData, CENTER)

     imin  = tileLimits (LOW, IAXIS)
     jmin  = tileLimits (LOW, JAXIS)
     imax  = tileLimits (HIGH,IAXIS)
     jmax  = tileLimits (HIGH,JAXIS)

     DeltaI         = delta (IAXIS)
     DeltaIHalf     = DeltaI * HALF
     DeltaIcube     = DeltaI * DeltaI * DeltaI
     DeltaJ         = delta (JAXIS)
     DeltaJHalf     = DeltaJ * HALF
     DeltaJSine     = sin (DeltaJ)
     DeltaJHalfSine = sin (DeltaJHalf)

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
!
!
!          ...Create all the cell info needed and place into proper radial bin array places.
!             The cell volume is:
!
!                      (2/3) * pi * (R^3 - r^3) * (cos[t] - cos[T])
!
!             where
!
!                       r  =  left-most (smaller) radial cell distance
!                       R  =  right-most (larger) radial cell distance
!                       t  =  smaller cell angle
!                       T  =  larger cell angle
!
!             The above formula can be obtained from the simple formula for a spherical
!             cone and noting that the cell can be obtained by adding and subtracting
!             the appropriate spherical cones. The cell looks like a torus with a 2D polar
!             cell cross-sectional shape. Since our angular and radial measures are based
!             on the 2D polar cell's center, we have:
!
!                       r  =  Rsph - Dr/2
!                       R  =  Rsph + Dr/2
!                       t  =  theta - Dt/2
!                       T  =  theta + Dt/2
!
!             with Dt and Dr being the cell's angular and radial delta values.
!             Hence the cell volume becomes, after using some trig identities:
!
!               (4pi * (Rsph)^2 * Dr + (pi/3) * Dr^3) * sin[Dt/2] * sin[theta]
!
!
     alpha       = TWO * DeltaJHalfSine * DeltaJHalfSine
     beta        = DeltaJSine

     theta       = bndBoxJLow + DeltaJHalf
     thetaSine   = sin (theta)
     thetaCosine = cos (theta)

     do j = jmin,jmax

      angularVolumePart = DeltaJHalfSine * thetaSine

      Rsph = bndBoxILow + DeltaIHalf

      do i = imin,imax

       cellVolume  = angularVolumePart * (gr_mpoleFourPi * Rsph * Rsph * DeltaI + gr_mpoleThirdPi * DeltaIcube)
       cellDensity = solnData (idensvar,i,j,1)
       cellMass    = cellDensity * cellVolume

       x = Rsph * thetaSine
       z = Rsph * thetaCosine - gr_mpoleZcenter

       r = sqrt (x * x + z * z)
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
           gr_mpoleQdataCells2D      (1,nQ) % coord1    = z
           gr_mpoleQdataCells2D      (1,nQ) % cellMass  = cellMass
           gr_mpoleQdataCells2D      (1,nQ) % radius    = r

       else

           nC = gr_mpoleQnumberOfCells (used) + 1

           gr_mpoleQnumberOfCells    (used)             = nC
           gr_mpoleQdataCells2D   (nC,used) % coord1    = z
           gr_mpoleQdataCells2D   (nC,used) % cellMass  = cellMass
           gr_mpoleQdataCells2D   (nC,used) % radius    = r

       end if

       Rsph = Rsph + DeltaI
      end do

      thetaSineSave = thetaSine
      thetaSine     = thetaSine   - (alpha * thetaSine   - beta * thetaCosine  )
      thetaCosine   = thetaCosine - (alpha * thetaCosine + beta * thetaSineSave)
     end do

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
!$omp end single
!
!
!    ...Call the radial bin clustered moment evaluation routine (all threads).
!
!
  call gr_mpoleMomBins2Dspherical (maxQtype)
!
!
!    ...Deallocate used arrays.
!
!
!$omp single
  deallocate (gr_mpoleQ             )
  deallocate (gr_mpoleQnumberOfCells)
  deallocate (gr_mpoleQdataCells2D  )
!$omp end single


!    ...Ready!


  return
end subroutine gr_mpoleMom2Dspherical
