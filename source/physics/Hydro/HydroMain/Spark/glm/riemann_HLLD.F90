!!****if* source/physics/Hydro/HydroMain/Spark/glm/riemann
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
!!  riemann
!!
!! SYNOPSIS
!!
!!  riemann( integer(IN) :: dir,
!!               real(IN)    :: VL(HY_NUM_VARS),
!!               real(IN)    :: VR(HY_NUM_VARS),
!!               real(OUT)   :: Fstar(HY_NUM_FLUX),
!!               real(OUT)   :: speed,
!!               integer(OUT):: ierr)
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  VL     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  VR     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!            (includes face pressure at the end)
!!  speed  - fastest signal velocity to compute dt
!!  ierr   - a flag to check unphysical negative state (0 is ok; 1 is bad)
!!
!! DESCRIPTION
!!
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!!   The HLLD Riemann fan:
!!
!!            SL*       SM       SR*
!!   SL        \        |        /        SR
!!     *        \       |       /        *
!!       *   UL* \ UL** | UR** / UR*   *
!!         *      \     |     /      *
!!           *     \    |    /     *
!!             *    \   |   /    *
!!           UL  *   \  |  /   *   UR
!!                 *  \ | /  *
!!                   * \|/ *
!!   --------------------------------------
!!
!! REFERENCE
!!
!!  * Miyoshi & Kusano, JCP, 208:315-344, 2005
!!
!!***

subroutine riemann(dir,VL,VR,inShock,Fstar,speed,ierr)

  use Driver_interface, ONLY : Driver_abort
  use Hydro_data,       ONLY : hy_tiny, hy_hybridRiemann, hy_C_hyp

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_NUM_VARS), intent(IN)  :: VL, VR
  real, dimension(HY_NUM_FLUX),   intent(OUT) :: Fstar
  logical, intent(IN) :: inShock
  real,    intent(OUT) :: speed
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real, parameter :: epsilon = 1.e-4
  real :: cfL,cfR,aL2,aR2,magBL2,magBR2
  real :: magnL,magnR,magtL,magtR,velnL,velnR
  real :: SM,SL,SR,SL2,SR2
  real :: dStarL,dStarR,prestL,prestR,pres,Bn_hll,denom
  real :: scrch1L,scrch1R,scrch2L,scrch2R
  real :: scrch3L,scrch3R,scrch4L,scrch4R,scrch5L,scrch5R
  real :: velxStarL,velyStarL,velzStarL,velxStarR,velyStarR,velzStarR
  real :: magxStarL,magyStarL,magzStarL,magxStarR,magyStarR,magzStarR
  real :: velxStar2,velyStar2,velzStar2,magxStar2,magyStar2,magzStar2
  real :: signumBn
  real, dimension(HY_NUM_FLUX)  :: UL,UR,Uhll,UCstarL,UCstarR,UCstar2L,UCstar2R
  real, dimension(HY_NUM_FLUX) :: FL,FR
  logical :: degeneracyHLLD
  real :: Bn_glm, Psi_glm


  ! Set no error
  ierr = 0

  ! Set no degeneracy for MHD
  degeneracyHLLD=.false.

  ! Normal velocity
  velnL = VL(HY_VELX+dir-1)
  velnR = VR(HY_VELX+dir-1)

  ! Normal and transverse magnetic components
  magnL = VL(HY_MAGX+dir-1)
  magnR = VR(HY_MAGX+dir-1)

  ! Magnitude of B fields
  magBL2 = dot_product(VL(HY_MAGX:HY_MAGZ),VL(HY_MAGX:HY_MAGZ))
  magBR2 = dot_product(VR(HY_MAGX:HY_MAGZ),VR(HY_MAGX:HY_MAGZ))

  ! Total pressure
  prestL = VL(HY_PRES) + 0.5*magBL2
  prestR = VR(HY_PRES) + 0.5*magBR2

  ! Magnitude of transverse fields
  magtL = magBL2-magnL*magnL
  magtR = magBR2-magnR*magnR

  ! Normalize by density
  magBL2= magBL2/VL(HY_DENS)
  magBR2= magBR2/VR(HY_DENS)

  ! Set sound speed
  aL2 = VL(HY_GAMC)*VL(HY_PRES)/VL(HY_DENS)
  aR2 = VR(HY_GAMC)*VR(HY_PRES)/VR(HY_DENS)

  ! Fastest magnetoacoustic waves
  cfL = 0.5*(aL2+magBL2+sqrt((aL2+magBL2)**2-4.*aL2*magnL*magnL/VL(HY_DENS)))
  cfR = 0.5*(aR2+magBR2+sqrt((aR2+magBR2)**2-4.*aR2*magnR*magnR/VR(HY_DENS)))

  ! Check unphysical negativity
  if ((VL(HY_DENS) < hy_tiny .and. VL(HY_DENS) > 0.) .or. &
       (VR(HY_DENS) < hy_tiny .and. VR(HY_DENS) > 0.) .or. &
       (VL(HY_PRES) < hy_tiny .and. VL(HY_PRES) > 0.) .or. &
       (VR(HY_PRES) < hy_tiny .and. VR(HY_PRES) > 0.)) then
     ! This could be vacuum limit. We return with zero flux.
     Fstar = 0.
     return
  elseif (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  end if

  cfL   = sqrt(cfL)
  cfR   = sqrt(cfR)

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velnL - cfL, velnR - cfR)
  SR = max(velnL + cfL, velnR + cfR)

  ! Output maximum local wave speed for dt calculation
  speed = max(abs(SL),abs(SR))

  ! Convert primitive variables to conservative variables
  call prim2con(VL,UL)
  call prim2con(VR,UR)
  call prim2flx(dir,VL,FL)
  call prim2flx(dir,VR,FR)

  ! Get HLL states for later use
  if (SL > 0.) then
     Uhll = UL
  elseif ((SL <= 0.) .and. (SR >= 0.)) then
     Uhll = (SR*UR - SL*UL - FR + FL)/(SR - SL)
  else
     Uhll = UR
  endif


  ! Calculate intermediate states ---------------------------------------------
  Bn_hll = Uhll(HY_FMGX+dir-1) !=(SR*magnR-SL*magnL)/(SR-SL)

  !!***************************************
  !! (I)    UL* and UR* regions           *
  !!***************************************
  ! Normal velocity component and the middle wave SM
  ! SM = u*L = u**L = u*R = u**R
  SM =(VR(HY_DENS)*velnR*(SR-velnR)-VL(HY_DENS)*velnL*(SL-velnL)&
       -prestR+prestL-magnL*magnL+magnR*magnR)/&
       (VR(HY_DENS)*(SR-velnR)-VL(HY_DENS)*(SL-velnL))


  ! Convenient parameters
  scrch1L = SL-velnL
  scrch2L = SL-SM
  scrch3L = SM-velnL

  scrch1R = SR-velnR
  scrch2R = SR-SM
  scrch3R = SM-velnR


  ! Total pressure in the whole Riemann fan
  ! pres*L = pres*R = pres**L = pres**R = pres
  pres = scrch1R*VR(HY_DENS)*prestL &
       -scrch1L*VL(HY_DENS)*prestR &
       +VL(HY_DENS)*VR(HY_DENS)*scrch1R*scrch1L*(velnR-velnL)

  pres = pres/(scrch1R*VR(HY_DENS)-scrch1L*VL(HY_DENS))

  ! Densities in UL* and UR*
  dStarL = UL(HY_MASS)*scrch1L/scrch2L
  dStarR = UR(HY_MASS)*scrch1R/scrch2R

  SL2 = SM - abs(Bn_hll)/sqrt(dStarL) ! = SL*
  SR2 = SM + abs(Bn_hll)/sqrt(dStarR) ! = SR*

  ! Check if degeneracy happens:
  ! This is the case of cf=ca, ca>a
  ! (1) transverse B=0 (note: in general, when trans B=0, cf=ca when ca>a; cs=a when ca<a), and
  ! (2) normal B .ne. 0, and strong to give ca>a.
  ! If this happens, we use HLLC (see PLUTO implementation as well)
  ! In spark, we shall just resort to HLL if this happens.
  if ( (SL2-SL) < epsilon*(SM-SL) ) then
     degeneracyHLLD=.true.
  endif
  if ( (SR-SR2) < epsilon*(SR-SM) ) then
     degeneracyHLLD=.true.
  endif
  if (degeneracyHLLD) then
     if (SL > 0.) then
        Fstar = FL
     elseif (SR < 0.) then
        Fstar = FR
     else !if ((SL <= 0.) .and. (SR >= 0.)) then
        Fstar = (SR*FL - SL*FR + SR*SL*(UR - UL))/(SR - SL)
     endif
     return
  endif




  denom   = VL(HY_DENS)*scrch1L*scrch2L-magnL*magnL
  scrch4L = scrch3L/denom
  scrch5L = (VL(HY_DENS)*scrch1L*scrch1L-magnL*magnL)/denom

  denom   = VR(HY_DENS)*scrch1R*scrch2R-magnR*magnR
  scrch4R = scrch3R/denom
  scrch5R = (VR(HY_DENS)*scrch1R*scrch1R-magnR*magnR)/denom


  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !!    Proceed calculating left star regions     !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  ! Left primitive variables
  select case (dir)
  case (IAXIS)
     magxStarL = Bn_hll
     velxStarL = SM
     magyStarL = VL(HY_MAGY)*scrch5L
     magzStarL = VL(HY_MAGZ)*scrch5L
     velyStarL = VL(HY_VELY) - VL(HY_MAGX)*VL(HY_MAGY)*scrch4L
     velzStarL = VL(HY_VELZ) - VL(HY_MAGX)*VL(HY_MAGZ)*scrch4L
  case (JAXIS)
     magyStarL = Bn_hll
     velyStarL = SM
     magxStarL = VL(HY_MAGX)*scrch5L
     magzStarL = VL(HY_MAGZ)*scrch5L
     velxStarL = VL(HY_VELX) - VL(HY_MAGY)*VL(HY_MAGX)*scrch4L
     velzStarL = VL(HY_VELZ) - VL(HY_MAGY)*VL(HY_MAGZ)*scrch4L
  case (KAXIS)
     magzStarL = Bn_hll
     velzStarL = SM
     magxStarL = VL(HY_MAGX)*scrch5L
     magyStarL = VL(HY_MAGY)*scrch5L
     velxStarL = VL(HY_VELX) - VL(HY_MAGZ)*VL(HY_MAGX)*scrch4L
     velyStarL = VL(HY_VELY) - VL(HY_MAGZ)*VL(HY_MAGY)*scrch4L
  end select

  ! Left conserved variables
  UCstarL(HY_MASS)= dStarL
  UCstarL(HY_XMOM)= UCstarL(HY_MASS)*velxStarL
  UCstarL(HY_YMOM)= UCstarL(HY_MASS)*velyStarL
  UCstarL(HY_ZMOM)= UCstarL(HY_MASS)*velzStarL

  UCstarL(HY_FMGX)= magxStarL
  UCstarL(HY_FMGY)= magyStarL
  UCstarL(HY_FMGZ)= magzStarL
  UCstarL(HY_ENER)= scrch1L*UL(HY_ENER)-prestL*velnL+pres*SM+&
       Bn_hll*(dot_product(VL(HY_VELX:HY_VELZ),VL(HY_MAGX:HY_MAGZ))&
       -velxStarL*UCstarL(HY_FMGX)&
       -velyStarL*UCstarL(HY_FMGY)&
       -velzStarL*UCstarL(HY_FMGZ))
  UCstarL(HY_ENER) =  UCstarL(HY_ENER)/scrch2L



  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! (Id) Proceed calculating right star regions  !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  ! Right primitive variables
  select case (dir)
  case (IAXIS)
     magxStarR = Bn_hll
     velxStarR = SM
     magyStarR = VR(HY_MAGY)*scrch5R
     magzStarR = VR(HY_MAGZ)*scrch5R
     velyStarR = VR(HY_VELY) - VR(HY_MAGX)*VR(HY_MAGY)*scrch4R
     velzStarR = VR(HY_VELZ) - VR(HY_MAGX)*VR(HY_MAGZ)*scrch4R
  case (JAXIS)
     magyStarR = Bn_hll
     velyStarR = SM
     magxStarR = VR(HY_MAGX)*scrch5R
     magzStarR = VR(HY_MAGZ)*scrch5R
     velxStarR = VR(HY_VELX) - VR(HY_MAGY)*VR(HY_MAGX)*scrch4R
     velzStarR = VR(HY_VELZ) - VR(HY_MAGY)*VR(HY_MAGZ)*scrch4R
  case (KAXIS)
     magzStarR = Bn_hll
     velzStarR = SM
     magxStarR = VR(HY_MAGX)*scrch5R
     magyStarR = VR(HY_MAGY)*scrch5R
     velxStarR = VR(HY_VELX) - VR(HY_MAGZ)*VR(HY_MAGX)*scrch4R
     velyStarR = VR(HY_VELY) - VR(HY_MAGZ)*VR(HY_MAGY)*scrch4R
  end select

  ! Right conserved variables
  UCstarR(HY_MASS)= dStarR
  UCstarR(HY_XMOM)= UCstarR(HY_MASS)*velxStarR
  UCstarR(HY_YMOM)= UCstarR(HY_MASS)*velyStarR
  UCstarR(HY_ZMOM)= UCstarR(HY_MASS)*velzStarR

  UCstarR(HY_FMGX)= magxStarR
  UCstarR(HY_FMGY)= magyStarR
  UCstarR(HY_FMGZ)= magzStarR
  UCstarR(HY_ENER)= scrch1R*UR(HY_ENER)-prestR*velnR+pres*SM+&
       Bn_hll*(dot_product(VR(HY_VELX:HY_VELZ),VR(HY_MAGX:HY_MAGZ))&
       -velxStarR*UCstarR(HY_FMGX)&
       -velyStarR*UCstarR(HY_FMGY)&
       -velzStarR*UCstarR(HY_FMGZ))
  UCstarR(HY_ENER) = UCstarR(HY_ENER)/scrch2R
  !! Done with calculating UL* and UR* regions !!



  !!***************************************
  !! (II)    UL** and UR** regions        *
  !!***************************************
  ! First calculate SL* and SR*
  SL2 = SM - abs(Bn_hll)/sqrt(UCstarL(HY_MASS)) ! = SL*
  SR2 = SM + abs(Bn_hll)/sqrt(UCstarR(HY_MASS)) ! = SR*

  ! Densities
  UCstar2L(HY_MASS) = UCstarL(HY_MASS)
  UCstar2R(HY_MASS) = UCstarR(HY_MASS)

  scrch1L = sqrt(UCstarL(HY_MASS))
  scrch1R = sqrt(UCstarR(HY_MASS))
  scrch2L = 1./(scrch1L + scrch1R)
  scrch2R = scrch2L

  signumBn = signum(Bn_hll)

  select case (dir)
  case (IAXIS)
     ! Left primitive variables
     velxStar2 = SM
     velyStar2 = (scrch1L*velyStarL+scrch1R*velyStarR&
          +(UCstarR(HY_FMGY)-UCstarL(HY_FMGY))*signumBn)*scrch2L
     velzStar2 = (scrch1L*velzStarL+scrch1R*velzStarR&
          +(UCstarR(HY_FMGZ)-UCstarL(HY_FMGZ))*signumBn)*scrch2L

     magxStar2 = Bn_hll
     magyStar2 = (scrch1L*magyStarR+scrch1R*magyStarL&
          +scrch1L*scrch1R*(velyStarR-velyStarL)*signumBn)&
          *scrch2L
     magzStar2 = (scrch1L*magzStarR+scrch1R*magzStarL&
          +scrch1L*scrch1R*(velzStarR-velzStarL)*signumBn)&
          *scrch2L

  case (JAXIS)
     ! Left primitive variables
     velxStar2 = (scrch1L*velxStarL+scrch1R*velxStarR&
          +(UCstarR(HY_FMGX)-UCstarL(HY_FMGX))*signumBn)*scrch2L
     velyStar2 = SM
     velzStar2 = (scrch1L*velzStarL+scrch1R*velzStarR&
          +(UCstarR(HY_FMGZ)-UCstarL(HY_FMGZ))*signumBn)*scrch2L

     magxStar2 = (scrch1L*magxStarR+scrch1R*magxStarL&
          +scrch1L*scrch1R*(velxStarR-velxStarL)*signumBn)&
          *scrch2L

     magyStar2 = Bn_hll
     magzStar2 = (scrch1L*magzStarR+scrch1R*magzStarL&
          +scrch1L*scrch1R*(velzStarR-velzStarL)*signumBn)&
          *scrch2L

  case (KAXIS)
     ! Left primitive variables
     velxStar2 = (scrch1L*velxStarL+scrch1R*velxStarR&
          +(UCstarR(HY_FMGX)-UCstarL(HY_FMGX))*signumBn)*scrch2L
     velyStar2 = (scrch1L*velyStarL+scrch1R*velyStarR&
          +(UCstarR(HY_FMGY)-UCstarL(HY_FMGY))*signumBn)*scrch2L
     velzStar2 = SM

     magxStar2 = (scrch1L*magxStarR+scrch1R*magxStarL&
          +scrch1L*scrch1R*(velxStarR-velxStarL)*signumBn)&
          *scrch2L
     magyStar2 = (scrch1L*magyStarR+scrch1R*magyStarL&
          +scrch1L*scrch1R*(velyStarR-velyStarL)*signumBn)&
          *scrch2L
     magzStar2 = Bn_hll

  end select

  ! Left conservative variables
  UCstar2L(HY_XMOM) = UCstar2L(HY_MASS)*velxStar2
  UCstar2L(HY_YMOM) = UCstar2L(HY_MASS)*velyStar2
  UCstar2L(HY_ZMOM) = UCstar2L(HY_MASS)*velzStar2

  UCstar2L(HY_FMGX) = magxStar2
  UCstar2L(HY_FMGY) = magyStar2
  UCstar2L(HY_FMGZ) = magzStar2
  UCstar2L(HY_ENER) = UCstarL(HY_ENER)-sqrt(UCstarL(HY_MASS))*signumBn*&
       (velxStarL*UCstarL (HY_FMGX)&
       +velyStarL*UCstarL (HY_FMGY)&
       +velzStarL*UCstarL (HY_FMGZ)&
       -velxStar2*UCstar2L(HY_FMGX)&
       -velyStar2*UCstar2L(HY_FMGY)&
       -velzStar2*UCstar2L(HY_FMGZ))

  ! Right conservative variables
  UCstar2R(HY_XMOM) = UCstar2R(HY_MASS)*velxStar2
  UCstar2R(HY_YMOM) = UCstar2R(HY_MASS)*velyStar2
  UCstar2R(HY_ZMOM) = UCstar2R(HY_MASS)*velzStar2

  UCstar2R(HY_FMGX) = magxStar2
  UCstar2R(HY_FMGY) = magyStar2
  UCstar2R(HY_FMGZ) = magzStar2
  UCstar2R(HY_ENER) = UCstarR(HY_ENER)+sqrt(UCstarR(HY_MASS))*signumBn*&
       (velxStarR*UCstarR (HY_FMGX)&
       +velyStarR*UCstarR (HY_FMGY)&
       +velzStarR*UCstarR (HY_FMGZ)&
       -velxStar2*UCstar2R(HY_FMGX)&
       -velyStar2*UCstar2R(HY_FMGY)&
       -velzStar2*UCstar2R(HY_FMGZ))
  !! END of calculating all HLLD states !!

  !!***************************************
  !! (III) HLLD fluxes                    *
  !!***************************************
  if (SL >= 0.) then
     Fstar = FL
  elseif ((SL < 0.) .and. (SL2 >= 0.)) then
     Fstar = FL + SL*(UCstarL - UL)
  elseif ((SL2 < 0.) .and. (SM >= 0.)) then
     Fstar = FL + SL2*(UCstar2L - UCstarL) + SL*( UCstarL - UL)
  elseif ((SM < 0.) .and. (SR2 >= 0.)) then
     Fstar = FR + SR2*(UCstar2R - UCstarR) + SR*( UCstarR - UR)
  elseif ((SR2 < 0.) .and. (SR >= 0.)) then
     Fstar = FR + SR*(UCstarR - UR)
  else
     Fstar = FR
  endif

  ! See Mignone & Tzeferacos 2010, sec. 3.2
  Bn_glm  = 0.5*(VL(HY_MAGX+dir-1)+VR(HY_MAGX+dir-1)) - &
       0.5/hy_C_hyp*(VR(HY_PSIB) - VL(HY_PSIB))
  Psi_glm = 0.5*(VL(HY_PSIB)+VR(HY_PSIB)) - &
                   0.5*hy_C_hyp*(VR(HY_MAGX+dir-1)-VL(HY_MAGX+dir-1))
  Fstar(HY_FMGX+dir-1) = Psi_glm
  Fstar(HY_FPSI) = hy_C_hyp*hy_C_hyp*Bn_glm

end subroutine riemann

function signum(x)
  implicit none
  real :: x,signum
  if (x == 0.) then
     signum = 0.
  else
     signum=sign(.5,x)-sign(.5,-x)
  endif
end function signum

#include "prim2con.F90"
#include "prim2flx.F90"
