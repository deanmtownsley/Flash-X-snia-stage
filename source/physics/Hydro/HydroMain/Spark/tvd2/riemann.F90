!!****if* source/physics/Hydro/HydroMain/Spark/tvd2/riemann
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
!!  hy_rk_hlle
!!
!! SYNOPSIS
!!
!!  hy_rk_hlle(integer(IN) :: dir,
!!             real(IN)    :: VL(HY_NUM_VARS),
!!             real(IN)    :: VR(HY_NUM_VARS),
!!             real(OUT)   :: Fstar(HY_NUM_FLUX),
!!             real(OUT)   :: speed,
!!             integer(OUT):: ierr)
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  VL     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  VR     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!  speed  - fastest signal velocity to compute dt
!!  ierr   - a flag to check unphysical negative state (0 is ok; 1 is bad)
!!
!! DESCRIPTION
!!
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!!   The HLLE Riemann fan:
!!
!!            SL                  SR
!!             \                 /
!!              \               /
!!               \      U*     /
!!                \           /
!!                 \         /
!!                  \       /
!!           UL      \     /       UR
!!                    \   /
!!                     \ /
!!   --------------------------------------
!!
!! REFERENCES
!!
!!  * Harten, Lax and van Leer, SIAM  Review, 25(1):35--61, 1983
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!
!!***

subroutine riemann(dir,VL,VR,inShock,Fstar,speed,ierr)

    use Driver_interface, ONLY : Driver_abort

  implicit none

#include "Simulation.h"
#include "Spark.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_NUM_VARS), intent(IN)  :: VL, VR
  logical, intent(IN) :: inShock
  real, dimension(HY_NUM_FLUX), intent(OUT) :: Fstar
  real,    intent(OUT) :: speed
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real :: SL,SR,cfL,cfR,aL2,aR2,velNL,velNR
  real :: magBL2,magBR2,magNL,magNR
  real, dimension(HY_NUM_FLUX) :: UL,UR,FL,FR
  real, parameter :: tiny=1.e-32

  ! Set no error to begin with
  ierr = 0
  Fstar = 0.0

  ! Normal velocity
  velNL = VL(HY_VELX+dir-1)
  velNR = VR(HY_VELX+dir-1)

  ! Set sound speed
  aL2   = VL(HY_GAMC)*VL(HY_PRES)/VL(HY_DENS)
  aR2   = VR(HY_GAMC)*VR(HY_PRES)/VR(HY_DENS)

  ! Check unphysical negativity
  if ((VL(HY_DENS) < tiny .and. VL(HY_DENS) > 0.) .or. &
       (VR(HY_DENS) < tiny .and. VR(HY_DENS) > 0.) .or. &
       (VL(HY_PRES) < tiny .and. VL(HY_PRES) > 0.) .or. &
       (VR(HY_PRES) < tiny .and. VR(HY_PRES) > 0.)) then
     ! This could be vacuum limit. We return with zero flux.
     Fstar = 0.
     return
  elseif (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  else
     cfL = sqrt(aL2)
     cfR = sqrt(aR2)
  endif

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velNL - cfL, velNR - cfR)
  SR = max(velNL + cfL, velNR + cfR)

  ! Output maximum local wave speed for dt calculation
  !speed = abs(velNL)+0.5*(cfL+cfR)
  speed = max(abs(SL),abs(SR))

  ! Convert primitive variables to conservative variables
  call prim2con(VL,UL)
  call prim2con(VR,UR)
  call prim2flx(dir,VL,FL)
  call prim2flx(dir,VR,FR)

  if (SL > 0.) then
     !Ustar(HY_DENS:HY_PRES) = UL(HY_DENS:HY_PRES)
     Fstar = FL
  elseif (SR < 0.) then
     !Ustar(HY_DENS:HY_PRES) = UR(HY_DENS:HY_PRES)
     Fstar = FR
  else !if ((SL <= 0.) .and. (SR >= 0.)) then
     !Ustar(HY_DENS:HY_PRES) = (SR*UR - SL*UL - FR + FL)/(SR - SL)
     Fstar = (SR*FL - SL*FR + SR*SL*(UR - UL))/(SR - SL)
  endif

end subroutine riemann

#include "prim2con.F90"
#include "prim2flx.F90"
