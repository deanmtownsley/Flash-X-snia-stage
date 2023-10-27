!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpolePot3Dspherical
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
!!  gr_mpolePot3Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpolePot3Dspherical  (integer, intent(in) :: ipotvar)
!!
!! DESCRIPTION
!!
!!  Computes the potential field for a three-dimensional spherical geometry
!!  using the mass moments already calculated. On output the variable
!!  indexed by ipotvar contains the potential. The calculations are
!!  entirely local to each processor, since each processor has a local
!!  copy of the moments.
!!
!! ARGUMENTS
!!
!!  ipotvar : index to variable containing the potential
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePot3Dspherical (ipotvar)

    use Grid_interface,    ONLY : Grid_getTileIterator,   &
                                  Grid_releaseTileIterator
  
    use gr_mpoleData,      ONLY : gr_mpoleGravityConstant,        &
                                gr_mpoleNumberInv,              &
                                gr_mpoleTotalNrCosineMoments,   &
                                gr_mpoleDrInv,                  &
                                gr_mpoleDrInnerZoneInv,         &
                                gr_mpoleMaxL,                   &
                                gr_mpoleMax2L,                  &
                                gr_mpoleMaxM,                   &
                                gr_mpoleMaxLM,                  &
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
                                gr_mpoleOuterZoneQshift

  use gr_mpoleData,      ONLY : gr_mpoleXcenter,                &
                                gr_mpoleYcenter,                &
                                gr_mpoleZcenter,                &
                                gr_mpoleQDampingR,              &
                                gr_mpoleQDampingI,              &
                                gr_mpoleMomentR,                &
                                gr_mpoleMomentI

  use Grid_tile,         ONLY : Grid_tile_t
  use Grid_iterator,     ONLY : Grid_iterator_t

  implicit none
  
#include "Simulation.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: ipotvar

  logical :: i2, j2, k2
  logical :: innerZonePotential

  integer :: blockNr, blockID
  integer :: c,s
  integer :: DrUnit
  integer :: imax, imin, iC, iCmax, iF, iFmax, iB
  integer :: jmax, jmin, jC, jCmax, jF, jFmax, jB, jS
  integer :: kmax, kmin, kC, kCmax, kF, kFmax
  integer :: M,MM,L
  integer :: Q, Qlocal, Qlower, Qupper
  integer :: type
  integer :: zone

  integer :: tileLimits   (LOW:HIGH,1:MDIM)

  real    :: alpha_t, beta_t,alpha_p,beta_p
  real    :: bndBoxILow, bndBoxJLow, bndBoxKLow
  real    :: dampI, Idamp
  real    :: dampR, Rdamp
  real    :: DeltaI, DeltaIHalf, DeltaJHalf, DeltaJFourthSine
  real    :: DeltaJ, DeltaJSine, DeltaJHalfSine, DeltaKSine
  real    :: DeltaK, DeltaKHalfSine, DeltaKFourthSine
  real    :: f,g,h
  real    :: facePotential
  real    :: Ic0, Is0, Ic1, Is1, Ic2, Is2, IcL, IsL
  real    :: phi, phiCosine, phiSine, phiSineSave
  real    :: theta, thetaCosine, thetaSine, thetaSineSave
  real    :: Qfloat, QfracI, QfracR
  real    :: r, rR, rI, rlocal, rinDrs, rsqR, rinvI, rsqinvI
  real    :: Rc0, Rs0, Rc1, Rs1, Rc2, Rs2, RcL, RsL
  real    :: Rdamping, Idamping
  real    :: RdampingQuotient, IdampingQuotient
  real    :: RdotI, IdotR
  real    :: sclInv, lgnInv, expInv
  real    :: x, y, z
  real    :: xI, yI, zI
  real    :: xR, yR, zR
  real    :: Rsph

  real    :: delta           (1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real, parameter :: sixth = 1./6.
  real, parameter :: fifth = 1./5.

  real, pointer   :: solnData (:,:,:,:)

  type(Grid_tile_t)     :: tileDesc
  type(Grid_iterator_t) :: itor
!  
!
!     ...Sum quantities over all locally held leaf blocks.
!
!
  nullify(solnData)
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     tileLimits=tileDesc%limits

     call tileDesc%boundBox(bndBox)
     call tileDesc%deltas(delta)
     call tileDesc%getDataPtr(solnData, CENTER)

     imin  = tileLimits (LOW, IAXIS)
     jmin  = tileLimits (LOW, JAXIS)
     kmin  = tileLimits (LOW, KAXIS)  
     imax  = tileLimits (HIGH,IAXIS)
     jmax  = tileLimits (HIGH,JAXIS)
     kmax  = tileLimits (HIGH,KAXIS)


     jCmax = jmax - jmin + 1      ! # of local (inner) cells in j direction
     kCmax = kmax - kmin + 1      ! # of local (inner) cells in k direction
     jFmax = jCmax + jCmax        !  max local (inner) face index in j direction (1st face -> index 0)
     kFmax = kCmax + kCmax        !  max local (inner) face index in k direction (1st face -> index 0)

     DeltaI           = delta (IAXIS)
     DeltaIHalf       = DeltaI * HALF
     
     DeltaJ           = delta (JAXIS)
     DeltaJHalf       = DeltaJ * HALF
     DeltaJSine       = sin (DeltaJ)
     DeltaJHalfSine   = sin (DeltaJ * HALF)
     DeltaJFourthSine = sin (DeltaJ * FOURTH)
     
     DeltaK           = delta (KAXIS)
     DeltaKSine       = sin (DeltaK)
     DeltaKHalfSine   = sin (DeltaK * HALF)
     DeltaKFourthSine = sin (DeltaK * FOURTH)

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)
     bndBoxKLow = bndBox (LOW,KAXIS)

     solnData (ipotvar , imin:imax , jmin:jmax , kmin:kmax) = ZERO

!            The full set of regular and irregular solid harmonics:
!
!                          0 =< L =< gr_mpoleMaxL
!                          0 =< M =< L       (cosine part)
!                          1 =< M =< L       (sine part)
!
!            are calculated and the scalar products are formed with the corresponding
!            moments.
!
!            Explanation of 'QfracR' and 'QfracI':
!            ------------------------------------
!
!            For outer zone (statistical) face potentials, we have to extract the
!            proper fractions of the Q-bin moments. The following picture, showing
!            the relevant radial bins, clarifies what is being done:
!
!
!
!                                            Q-th bin
!
!              ----- MR(Q-1) ---- | ---- MR(Q)----- MI(Q) -----| ------ MI(Q+1)----
!                                                 |
!                                          QfracR | QfracI
!                                                 |
!                                              Qfloat
!
!
!            Here 'Qfloat' is the value obtained by the cell's face radial position
!            and from 'Qfloat' we determined the bin number 'Q'. The values of
!            'QfracR' and 'QfracI' denote the fractional location of 'Qfloat'
!            within the Q-th bin. Hence we always have: QfracR + QfracI = 1.
!            When forming the dot products between the current Q-th bin solid
!            harmonics and the Moments, we take only the relevant fraction of
!            the Moments in Q. This sees the Moments in Q as a kind of statistical
!            average in terms of the radii comprised within the Q-th bin. As an
!            example, we take the dot product between the regular solid harmonics
!            in the Q-th bin R(Q) with all the necessary irregular Moments MI.
!            We have:
!
!
!              RdotI = R(Q) dot { MI(Q+1) + QfracI * [MI(Q) - MI(Q+1)] }
!
!                    = R(Q) dot { [QfracR + QfracI] * MI(Q+1) + QfracI * [MI(Q) - MI(Q+1)] }
!
!                    = R(Q) dot { QfracR * I(Q+1) + QfracI * I(Q) }
!
!
!            The same applies for the dot product between the irregular solid
!            harmonics in the Q-th bin I(Q) and the regular Moments MR. We get:
!
!
!              IdotR = I(Q) dot { QfracI * MR(Q-1) + QfracR * MR(Q) }
!
!
!            For the inner zone face potentials, we combine the regular face moments
!            with the irregular Q-bin moments and the irregular face moments with the
!            regular (Q-1)-bin moments. This can readily be achieved by setting
!            QfracR = 0 and QfracI = 1.
!
!
!


     
     alpha_p     = TWO * DeltaKHalfSine * DeltaKHalfSine      ! for calculating sin (phi + n * DeltaK)
     beta_p      = DeltaKSine   


     phi       = bndBoxKLow                                     ! initial phi axis location of k face
     phiSine   = sin (phi)
     phiCosine = cos (phi)

     do kF = 0, kFmax, 2                                           ! loop over local k face indices

      kC = int (kF/2) + 1                                       ! local (inner) largest k cell index for k face
      k2 = mod (kF,2) == 0 .and. kC > 1 .and. kC < kCmax + 1    ! 2 cells k and k-1 share k face?
      kC = kmin - 1 + min (kC , kCmax)                          ! change to global (inner + guard) k cell index
      jB = 0
      jS = 2

      if (jS == 2) then
       alpha_t       = TWO * DeltaJHalfSine * DeltaJHalfSine          ! for calculating sin (theta + n * DeltaJ)
       beta_t        = DeltaJSine                                     ! and cos (theta + n * DeltaJ)
      else                                                          ! if step size is 1
       alpha_t       = TWO * DeltaJFourthSine * DeltaJFourthSine          ! for calculating sin (theta + n * DeltaJHalf)
       beta_t        = DeltaJHalfSine                                     ! and cos (theta + n * DeltaJHalf)
      endif   

      theta       = bndBoxJLow + jB * DeltaJHalf                 ! initial theta axis location of j face
      thetaSine   = sin (theta)
      thetaCosine = cos (theta)
    
      do jF = jB, jFmax, jS                                     ! loop over local j face indices

      jC = int (jF/2) + 1                                         ! local (inner) largest j cell index for j face
      j2 = mod (jF,2) == 0 .and. jC > 1 .and. jC < jCmax + 1                            ! 2 cells j and j-1 share j face?
      jC = jmin - 1 + min (jC , jCmax)                            ! change to global (inner + guard) j cell index

      Rsph = bndBoxILow + DeltaIHalf                              ! initial Rsph-axis location of i face

       do iC = imin, imax                                          ! loop over all radial cells

        x = Rsph * thetaSine * phiCosine - gr_mpoleXcenter
        y = Rsph * phiSine * thetaSine - gr_mpoleYcenter
        z = Rsph * thetaCosine - gr_mpoleZcenter
        
        r = sqrt (x * x + y * y + z * z)                                   ! radial position of face location
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
            QfracR = ZERO
            QfracI = ONE

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
                Qfloat = (rlocal * sclInv * gr_mpoleDrInv) ** expInv
            else if (type == ZONE_LOGARITHMIC) then
                lgnInv = gr_mpoleZoneLogNormInv (zone)
                Qfloat = expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE)
            end if

            Qlocal = ceiling (Qfloat)
            QfracI = real (Qlocal) - Qfloat
            QfracR = ONE - QfracI
            Q      = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

        end if
!
!
!        ...Calculate the (damped) moments.
!
!
        Rdamping       = gr_mpoleQDampingR (Q)
        Idamping       = gr_mpoleQDampingI (Q)
        IdampingQuotient = gr_mpoleQDampingI (Q+1) / gr_mpoleQDampingI (Q)

        rI = r * Rdamping                  ! irregular solid harmonics will cancel Rdamping
        rinvI = ONE / rI

        RcL = Idamping                     ! start undo 1/(Idamping)^(L+1) in irregular Q moments
        IcL = rinvI * Rdamping             ! start undo (Rdamping)^L in regular Q moments

        Idamp = IdampingQuotient           ! damping for irregular moments for L = 0
        Rdamp = ONE                        ! damping for   regular moments for L = 0

        RdotI = RcL * (QfracI * gr_mpoleMomentI (1,Q) + QfracR * Idamp * gr_mpoleMomentI (1,Q+1))
        IdotR = IcL * (QfracR * gr_mpoleMomentR (1,Q) + QfracI         * gr_mpoleMomentR (1,Q-1))

        facePotential = RdotI + IdotR      ! initialize face potential variable
!
!
!        ...Proceed, if we go beyond monopoles.
!
!
        if (gr_mpoleMaxL > 0) then

         RdampingQuotient = gr_mpoleQDampingR (Q) / gr_mpoleQDampingR (Q-1)

         zR = z * Idamping                 !  regular  solid harmonics will cancel Idamping
         rR = r * Idamping                 !  regular  solid harmonics will cancel Idamping
         zI = z * Rdamping                 ! irregular solid harmonics will cancel Rdamping

         rsqR = rR * rR
         rsqinvI = rinvI * rinvI

         Rc0 = RcL
         Ic0 = IcL
         Rc1 = RcL * zR                    !  regular  solid harmonic for cos 10
         Ic1 = IcL * zI * rsqinvI          ! irregular solid harmonic for cos 10

         Idamp = Idamp * IdampingQuotient  ! store damping for irregular moments for L = 1
         Rdamp = Rdamp * RdampingQuotient  ! store damping for   regular moments for L = 1
         dampI = Idamp                     ! initialize for L accumulation damping
         dampR = Rdamp                     ! initialize for L accumulation damping

         RdotI = Rc1 * (QfracI * gr_mpoleMomentI (2,Q) + QfracR * dampI * gr_mpoleMomentI (2,Q+1))
         IdotR = Ic1 * (QfracR * gr_mpoleMomentR (2,Q) + QfracI * dampR * gr_mpoleMomentR (2,Q-1))

         facePotential = facePotential + RdotI + IdotR

         do L = 2,gr_mpoleMaxL

            h = real (L + L - 1)
            g = real ((L - 1) * (L - 1))
            f = gr_mpoleNumberInv (L) * gr_mpoleNumberInv (L)

            Rc2 = (h * zR * Rc1 - rsqR * Rc0) * f        !  regular  solid harmonic for cos L0
            Ic2 = (h * zI * Ic1 -    g * Ic0) * rsqinvI  ! irregular solid harmonic for cos L0

            dampI = dampI * IdampingQuotient             ! next L damping
            dampR = dampR * RdampingQuotient             ! next L damping

            RdotI = Rc2 * (QfracI * gr_mpoleMomentI (L+1,Q) + QfracR * dampI * gr_mpoleMomentI (L+1,Q+1))
            IdotR = Ic2 * (QfracR * gr_mpoleMomentR (L+1,Q) + QfracI * dampR * gr_mpoleMomentR (L+1,Q-1))

            facePotential = facePotential + RdotI + IdotR

            Rc0 = Rc1                                    !
            Rc1 = Rc2                                    ! swap, to accumulate next L regular
            Ic0 = Ic1                                    ! and irregular solid harmonic
            Ic1 = Ic2                                    !

         end do

         xR = x * Idamping                       !  regular  solid harmonics will cancel Idamping
         yR = y * Idamping                       !  regular  solid harmonics will cancel Idamping
         xI = x * Rdamping                       ! irregular solid harmonics will cancel Rdamping
         yI = y * Rdamping                       ! irregular solid harmonics will cancel Rdamping

         c = gr_mpoleMaxL + 2
         s = gr_mpoleTotalNrCosineMoments + 1

         Rc0 = - RcL * xR * HALF                 !  regular  solid harmonic for cos 11
         Rs0 = - RcL * yR * HALF                 !  regular  solid harmonic for sin 11
         Ic0 = - IcL * xI * rsqinvI              ! irregular solid harmonic for cos 11
         Is0 = - IcL * yI * rsqinvI              ! irregular solid harmonic for sin 11

         RdotI =   Rc0 * (QfracI * gr_mpoleMomentI (c,Q) + QfracR * Idamp * gr_mpoleMomentI (c,Q+1)) &
                 + Rs0 * (QfracI * gr_mpoleMomentI (s,Q) + QfracR * Idamp * gr_mpoleMomentI (s,Q+1))
         IdotR =   Ic0 * (QfracR * gr_mpoleMomentR (c,Q) + QfracI * Rdamp * gr_mpoleMomentR (c,Q-1)) &
                 + Is0 * (QfracR * gr_mpoleMomentR (s,Q) + QfracI * Rdamp * gr_mpoleMomentR (s,Q-1))

         facePotential = facePotential + RdotI + RdotI + IdotR + IdotR    ! include x2 for M /= 0 cases
!
!
!        ...Proceed, if we go beyond dipoles.
!
!
         if (gr_mpoleMaxL > 1) then

          RcL = Rc0                               ! store cos 11 for later (22 case)
          RsL = Rs0                               ! store sin 11 for later (22 case)
          IcL = Ic0                               ! store cos 11 for later (22 case)
          IsL = Is0                               ! store sin 11 for later (22 case)

          Rc1 =         zR * RcL                  !  regular  solid harmonic for cos 21
          Rs1 =         zR * RsL                  !  regular  solid harmonic for sin 21
          Ic1 = THREE * zI * IcL * rsqinvI        ! irregular solid harmonic for cos 21
          Is1 = THREE * zI * IsL * rsqinvI        ! irregular solid harmonic for sin 21

          Idamp = Idamp * IdampingQuotient        ! store damping for irregular moments for L = 2
          Rdamp = Rdamp * RdampingQuotient        ! store damping for   regular moments for L = 2
          dampI = Idamp                           ! initialize for L accumulation damping
          dampR = Rdamp                           ! initialize for L accumulation damping

          RdotI =   Rc1 * (QfracI * gr_mpoleMomentI (c+1,Q) + QfracR * dampI * gr_mpoleMomentI (c+1,Q+1)) &
                  + Rs1 * (QfracI * gr_mpoleMomentI (s+1,Q) + QfracR * dampI * gr_mpoleMomentI (s+1,Q+1))
          IdotR =   Ic1 * (QfracR * gr_mpoleMomentR (c+1,Q) + QfracI * dampR * gr_mpoleMomentR (c+1,Q-1)) &
                  + Is1 * (QfracR * gr_mpoleMomentR (s+1,Q) + QfracI * dampR * gr_mpoleMomentR (s+1,Q-1))

          facePotential = facePotential + RdotI + RdotI + IdotR + IdotR

          do L = 2,gr_mpoleMaxL-1                 ! M=1 and L=3,gr_mpoleMaxL (shifted L-loop!)

             h = real (L + L + 1)
             g = real ((L + 1) * (L - 1))
             f = gr_mpoleNumberInv (L + 2) * gr_mpoleNumberInv (L)

             Rc2 = (h * zR * Rc1 - rsqR * Rc0) * f        !  regular  solid harmonic for cos L1
             Rs2 = (h * zR * Rs1 - rsqR * Rs0) * f        !  regular  solid harmonic for sin L1
             Ic2 = (h * zI * Ic1 -    g * Ic0) * rsqinvI  ! irregular solid harmonic for cos L1
             Is2 = (h * zI * Is1 -    g * Is0) * rsqinvI  ! irregular solid harmonic for sin L1

             dampI = dampI * IdampingQuotient             ! next L damping
             dampR = dampR * RdampingQuotient             ! next L damping

             RdotI =   Rc2 * (QfracI * gr_mpoleMomentI (c+L,Q) + QfracR * dampI * gr_mpoleMomentI (c+L,Q+1)) &
                     + Rs2 * (QfracI * gr_mpoleMomentI (s+L,Q) + QfracR * dampI * gr_mpoleMomentI (s+L,Q+1))
             IdotR =   Ic2 * (QfracR * gr_mpoleMomentR (c+L,Q) + QfracI * dampR * gr_mpoleMomentR (c+L,Q-1)) &
                     + Is2 * (QfracR * gr_mpoleMomentR (s+L,Q) + QfracI * dampR * gr_mpoleMomentR (s+L,Q-1))

             facePotential = facePotential + RdotI + RdotI + IdotR + IdotR

             Rc0 = Rc1                                    !
             Rs0 = Rs1                                    !
             Rc1 = Rc2                                    !
             Rs1 = Rs2                                    ! swap, to accumulate next regular
             Ic0 = Ic1                                    ! and irregular solid harmonic
             Is0 = Is1                                    !
             Ic1 = Ic2                                    !
             Is1 = Is2                                    !

          end do

          c = c + gr_mpoleMaxL
          s = s + gr_mpoleMaxL

          do M = 2, gr_mpoleMaxL-1

             MM = M + M
             h = real (MM + 1)
             g = real (MM - 1)
             f = gr_mpoleNumberInv (MM)

             Rc0 =       (yR * RsL - xR * RcL) * f            !  regular  solid harmonic for cos MM
             Rs0 =     - (yR * RcL + xR * RsL) * f            !  regular  solid harmonic for sin MM
             Ic0 =   g * (yI * IsL - xI * IcL) * rsqinvI      ! irregular solid harmonic for cos MM
             Is0 = - g * (yI * IcL + xI * IsL) * rsqinvI      ! irregular solid harmonic for sin MM

             RdotI =   Rc0 * (QfracI * gr_mpoleMomentI (c,Q) + QfracR * Idamp * gr_mpoleMomentI (c,Q+1)) &
                     + Rs0 * (QfracI * gr_mpoleMomentI (s,Q) + QfracR * Idamp * gr_mpoleMomentI (s,Q+1))
             IdotR =   Ic0 * (QfracR * gr_mpoleMomentR (c,Q) + QfracI * Rdamp * gr_mpoleMomentR (c,Q-1)) &
                     + Is0 * (QfracR * gr_mpoleMomentR (s,Q) + QfracI * Rdamp * gr_mpoleMomentR (s,Q-1))

             facePotential = facePotential + RdotI + RdotI + IdotR + IdotR

             RcL = Rc0                                        ! store cos MM for M+1,M+1
             RsL = Rs0                                        ! store sin MM for M+1,M+1
             IcL = Ic0                                        ! store cos MM for M+1,M+1
             IsL = Is0                                        ! store sin MM for M+1,M+1

             Rc1 =     zR * RcL                               !  regular  solid harmonic for cos L=M+1,M
             Rs1 =     zR * RsL                               !  regular  solid harmonic for sin L=M+1,M
             Ic1 = h * zI * IcL * rsqinvI                     ! irregular solid harmonic for cos L=M+1,M
             Is1 = h * zI * IsL * rsqinvI                     ! irregular solid harmonic for sin L=M+1,M

             Idamp = Idamp * IdampingQuotient                 ! store damping for irregular moments for L = M+1
             Rdamp = Rdamp * RdampingQuotient                 ! store damping for   regular moments for L = M+1
             dampI = Idamp                                    ! initialize for L accumulation damping
             dampR = Rdamp                                    ! initialize for L accumulation damping

             RdotI =   Rc1 * (QfracI * gr_mpoleMomentI (c+1,Q) + QfracR * dampI * gr_mpoleMomentI (c+1,Q+1)) &
                     + Rs1 * (QfracI * gr_mpoleMomentI (s+1,Q) + QfracR * dampI * gr_mpoleMomentI (s+1,Q+1))
             IdotR =   Ic1 * (QfracR * gr_mpoleMomentR (c+1,Q) + QfracI * dampR * gr_mpoleMomentR (c+1,Q-1)) &
                     + Is1 * (QfracR * gr_mpoleMomentR (s+1,Q) + QfracI * dampR * gr_mpoleMomentR (s+1,Q-1))

             facePotential = facePotential + RdotI + RdotI + IdotR + IdotR

             do L = 2, gr_mpoleMaxL-M                         ! L=M+2,gr_mpoleMaxL (shifted L-loop!)

                h = real (L + L + MM - 1)
                g = real (L + MM - 1) * (L - 1)
                f = gr_mpoleNumberInv (L + MM) * gr_mpoleNumberInv (L)

                Rc2 = (h * zR * Rc1 - rsqR * Rc0) * f         !  regular  solid harmonic for cos LM
                Rs2 = (h * zR * Rs1 - rsqR * Rs0) * f         !  regular  solid harmonic for sin LM
                Ic2 = (h * zI * Ic1 -    g * Ic0) * rsqinvI   ! irregular solid harmonic for cos LM
                Is2 = (h * zI * Is1 -    g * Is0) * rsqinvI   ! irregular solid harmonic for sin LM

                dampI = dampI * IdampingQuotient              ! next L damping
                dampR = dampR * RdampingQuotient              ! next L damping

                RdotI =   Rc2 * (QfracI * gr_mpoleMomentI (c+L,Q) + QfracR * dampI * gr_mpoleMomentI (c+L,Q+1)) &
                        + Rs2 * (QfracI * gr_mpoleMomentI (s+L,Q) + QfracR * dampI * gr_mpoleMomentI (s+L,Q+1))
                IdotR =   Ic2 * (QfracR * gr_mpoleMomentR (c+L,Q) + QfracI * dampR * gr_mpoleMomentR (c+L,Q-1)) &
                        + Is2 * (QfracR * gr_mpoleMomentR (s+L,Q) + QfracI * dampR * gr_mpoleMomentR (s+L,Q-1))

                facePotential = facePotential + RdotI + RdotI + IdotR + IdotR

                Rc0 = Rc1                                     !
                Rs0 = Rs1                                     !
                Rc1 = Rc2                                     !
                Rs1 = Rs2                                     ! swap, to accumulate next regular
                Ic0 = Ic1                                     ! and irregular solid harmonic
                Is0 = Is1                                     !
                Ic1 = Ic2                                     !
                Is1 = Is2                                     !

             end do

             c = c + gr_mpoleMaxL - M + 1
             s = s + gr_mpoleMaxL - M + 1

          end do

          g = real (gr_mpoleMax2L - 1)
          f = gr_mpoleNumberInv (gr_mpoleMax2L)

          Rc0 =       (yR * RsL - xR * RcL) * f               !  regular  solid harmonic for cos MaxL,MaxL
          Rs0 =     - (yR * RcL + xR * RsL) * f               !  regular  solid harmonic for sin MaxL,MaxL
          Ic0 =   g * (yI * IsL - xI * IcL) * rsqinvI         ! irregular solid harmonic for cos MaxL,MaxL
          Is0 = - g * (yI * IcL + xI * IsL) * rsqinvI         ! irregular solid harmonic for sin MaxL,MaxL

          RdotI =   Rc0 * (QfracI * gr_mpoleMomentI (c,Q) + QfracR * Idamp * gr_mpoleMomentI (c,Q+1)) &
                  + Rs0 * (QfracI * gr_mpoleMomentI (s,Q) + QfracR * Idamp * gr_mpoleMomentI (s,Q+1))
          IdotR =   Ic0 * (QfracR * gr_mpoleMomentR (c,Q) + QfracI * Rdamp * gr_mpoleMomentR (c,Q-1)) &
                  + Is0 * (QfracR * gr_mpoleMomentR (s,Q) + QfracI * Rdamp * gr_mpoleMomentR (s,Q-1))

          facePotential = facePotential + RdotI + RdotI + IdotR + IdotR

         end if   ! gr_mpoleMaxL > 1       condition
        end if    ! gr_mpoleMaxL > 0       condition

        facePotential = - gr_mpoleGravityConstant * facePotential
!
!
!        ...Add the current face potential to the relevant cell(s) of the potential block.
!           Note, that a face can only be shared by up to 2 cells and no more. The most
!           common case (shared on the j-index) is tested first.
!
!

        if (j2 .and. k2) then
            solnData (ipotvar,iC,jC,kC-1) = solnData (ipotvar,iC,jC,kC-1) + facePotential 
            solnData (ipotvar,iC,jC,kC) = solnData (ipotvar,iC,jC,kC) + facePotential 
            solnData (ipotvar,iC,jC-1,kC) = solnData (ipotvar,iC,jC-1,kC) + facePotential 
            solnData (ipotvar,iC,jC-1 ,kC-1) = solnData (ipotvar,iC,jC-1 ,kC-1) + facePotential 
        else if (k2) then

            solnData (ipotvar,iC,jC,kC-1) = solnData (ipotvar,iC,jC,kC-1) + facePotential 
            solnData (ipotvar,iC,jC,kC  ) = solnData (ipotvar,iC,jC,kC  ) + facePotential 

        else if (j2) then
            solnData (ipotvar,iC,jC,kC  ) = solnData (ipotvar,iC,jC,kC  ) + facePotential 
            solnData (ipotvar,iC,jC-1,kC  ) = solnData (ipotvar,iC,jC-1,kC  ) + facePotential 

        else
            solnData (ipotvar,iC,jC,kC) = solnData (ipotvar,iC,jC,kC) + facePotential 

        end if

        Rsph = Rsph + DeltaI                  ! face increment radial axis
       end do
       thetaSineSave = thetaSine                                                  !
       thetaSine     = thetaSine   - (alpha_t * thetaSine   - beta_t * thetaCosine  ) ! face increment angular theta axis
       thetaCosine   = thetaCosine - (alpha_t * thetaCosine + beta_t * thetaSineSave) !
      end do

      phiSineSave = phiSine                                                      !
      phiSine     = phiSine   - (alpha_p * phiSine   - beta_p * phiCosine  )     ! face increment angular phi axis
      phiCosine   = phiCosine - (alpha_p * phiCosine + beta_p * phiSineSave)     !

     end do
!
!
!    ...Form the potential average in each cell.
!
!
     solnData (ipotvar,imin:imax,jmin:jmax,kmin:kmax) = .25 * solnData (ipotvar,imin:imax,jmin:jmax,kmin:kmax)
!
!
!    ...Get ready for retrieving next LEAF block for the current processor.
!
!
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)

!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpolePot3Dspherical
