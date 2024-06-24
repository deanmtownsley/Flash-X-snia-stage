!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!> @ingroup physics_Eos
!!
!! @brief Routine to access arbitrary EOS data from WeakLib for one zone only
!!
!! AUTHOR & DATE 
!!   S.M. Couch
!!   J.A. Harris
!!   April 2022
!!


subroutine eos_wlOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,xdEdT,xCs2,xXp,xXn,xXa,xXh,xAbar,xVar,varID,mode)

#include "Simulation.h"
#include "constants.h"

  use Driver_interface, ONLY : Driver_abort
  USE eos_wlData, ONLY: EosNewTable, nVariables

  USE wlInterpolationModule, ONLY: &
      LogInterpolateSingleVariable_3D_Custom_Point, &
      LogInterpolateDifferentiateSingleVariable_3D_Custom_Point
  USE wlEOSInversionModule, ONLY : &
      ComputeTemperatureWith_DEY_Single_Guess, &
      ComputeTemperatureWith_DSY_Single_Guess, &
      ComputeTemperatureWith_DPY_Single_Guess, &
      DescribeEOSInversionError

  implicit none

  real, intent(IN) :: xDens, xYe
  real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
  integer, intent(IN) :: mode, varID
  real, intent(OUT) :: xdEdT,xCs2,xXp,xXn,xXa,xXh,xAbar,xVar

  real :: xTemp_Guess
  real :: xMu_e, xMu_n, xMu_p, xGamc, xAh, xdEdX(3)

  integer :: iD_T, iT_T, iY_T
  integer :: iP_T, iS_T, iE_T
  integer :: iMe_T, iMp_T, iMn_T, iXp_T, iXn_T, iXa_T, iXh_T, iAh_T, iGm_T

  real :: OS_P, OS_E, OS_S
  real :: OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Ah, OS_Gm, OS_Var

  integer :: error_flag

  ! index var mapping:
  !  0 -> full eos call
  !  1 -> pressure
  !  2 -> energy
  !  3 -> entropy
  !  4 -> gamc
  !  5 -> mu_e
  !  6 -> mu_p
  !  7 -> mu_n
  !  8 -> xp
  !  9 -> xn
  ! 10 -> xa
  ! 11 -> xh
  ! 12 -> zh
  ! 13 -> ah
  ! 14 -> bh
  ! 15 -> thermal energy
  ! 16 -> 
  ! 17 -> 
  ! 18 -> 
  ! 19 -> 
  ! 20 -> mu_nue

  iD_T = EosNewTable % TS % Indices % iRho
  iT_T = EosNewTable % TS % Indices % iT
  iY_T = EosNewTable % TS % Indices % iYe

  iP_T  = EosNewTable % DV % Indices % iPressure
  iS_T  = EosNewTable % DV % Indices % iEntropyPerBaryon
  iE_T  = EosNewTable % DV % Indices % iInternalEnergyDensity
  iMe_T = EosNewTable % DV % Indices % iElectronChemicalPotential
  iMp_T = EosNewTable % DV % Indices % iProtonChemicalPotential
  iMn_T = EosNewTable % DV % Indices % iNeutronChemicalPotential
  iXp_T = EosNewTable % DV % Indices % iProtonMassFraction
  iXn_T = EosNewTable % DV % Indices % iNeutronMassFraction
  iXa_T = EosNewTable % DV % Indices % iAlphaMassFraction
  iXh_T = EosNewTable % DV % Indices % iHeavyMassFraction
  iAh_T = EosNewTable % DV % Indices % iHeavyMassNumber
  iGm_T = EosNewTable % DV % Indices % iGamma1

  OS_P  = EosNewTable % DV % Offsets(iP_T)
  OS_S  = EosNewTable % DV % Offsets(iS_T)
  OS_E  = EosNewTable % DV % Offsets(iE_T)
  OS_Me = EosNewTable % DV % Offsets(iMe_T)
  OS_Mp = EosNewTable % DV % Offsets(iMp_T)
  OS_Mn = EosNewTable % DV % Offsets(iMn_T)
  OS_Xp = EosNewTable % DV % Offsets(iXp_T)
  OS_Xn = EosNewTable % DV % Offsets(iXn_T)
  OS_Xa = EosNewTable % DV % Offsets(iXa_T)
  OS_Xh = EosNewTable % DV % Offsets(iXh_T)
  OS_Ah = EosNewTable % DV % Offsets(iAh_T)
  OS_Gm = EosNewTable % DV % Offsets(iGm_T)

  ASSOCIATE( D_T  => EosNewTable % TS % States(iD_T) % Values, &
             T_T  => EosNewTable % TS % States(iT_T) % Values, &
             Y_T  => EosNewTable % TS % States(iY_T) % Values, &
             P_T  => EosNewTable % DV % Variables(iP_T ) % Values, &
             S_T  => EosNewTable % DV % Variables(iS_T ) % Values, &
             E_T  => EosNewTable % DV % Variables(iE_T ) % Values, &
             Me_T => EosNewTable % DV % Variables(iMe_T) % Values, &
             Mn_T => EosNewTable % DV % Variables(iMn_T) % Values, &
             Mp_T => EosNewTable % DV % Variables(iMp_T) % Values, &
             Xp_T => EosNewTable % DV % Variables(iXp_T) % Values, &
             Xn_T => EosNewTable % DV % Variables(iXn_T) % Values, &
             Xa_T => EosNewTable % DV % Variables(iXa_T) % Values, &
             Xh_T => EosNewTable % DV % Variables(iXh_T) % Values, &
             Ah_T => EosNewTable % DV % Variables(iAh_T) % Values, &
             Gm_T => EosNewTable % DV % Variables(iGm_T) % Values )

  SELECT CASE( mode )
  CASE( MODE_DENS_EI )

     xTemp_Guess = xTemp

     CALL ComputeTemperatureWith_DEY_Single_Guess &
        ( xDens, xEner, xYe, D_T, T_T, Y_T, E_T, OS_E, xTemp, xTemp_Guess, &
          error_flag )

     IF( error_flag /= 0 ) THEN
        WRITE(*,*) 'ComputeTemperatureWith_DEY ERROR:'
        WRITE(*,'(5A15,A10)') 'xDens','xEner','xYe','xTemp','Guess','Error'
        WRITE(*,'(5ES15.4,I10)') xDens, xEner, xYe, xTemp, xTemp_Guess, error_flag
        CALL DescribeEOSInversionError( error_flag )
        CALL Driver_abort('[eos_wlOneZone] EOS Inversion Error')
     END IF

  CASE( MODE_DENS_TEMP )


  CASE( MODE_DENS_ENTR )

     xTemp_Guess = xTemp

     CALL ComputeTemperatureWith_DSY_Single_Guess &
        ( xDens, xEntr, xYe, D_T, T_T, Y_T, S_T, OS_S, xTemp, xTemp_Guess, &
          error_flag )

     IF( error_flag /= 0 ) THEN
        WRITE(*,*) 'ComputeTemperatureWith_DSY ERROR:'
        WRITE(*,'(5A15,A10)') 'xDens','xEntr','xYe','xTemp','Guess','Error'
        WRITE(*,'(5ES15.4,I10)') xDens, xEntr, xYe, xTemp, xTemp_Guess, error_flag
        CALL DescribeEOSInversionError( error_flag )
        CALL Driver_abort('[eos_wlOneZone] EOS Inversion Error')
     END IF

  CASE( MODE_DENS_PRES )

     xTemp_Guess = xTemp

     CALL ComputeTemperatureWith_DPY_Single_Guess &
        ( xDens, xPres, xYe, D_T, T_T, Y_T, P_T, OS_P, xTemp, xTemp_Guess, &
          error_flag )

     IF( error_flag /= 0 ) THEN
        WRITE(*,*) 'ComputeTemperatureWith_DPY ERROR:'
        WRITE(*,'(5A15,A10)') 'xDens','xPres','xYe','xTemp','Guess','Error'
        WRITE(*,'(5ES15.4,I10)') xDens, xPres, xYe, xTemp, xTemp_Guess, error_flag
        CALL DescribeEOSInversionError( error_flag )
        CALL Driver_abort('[eos_wlOneZone] EOS Inversion Error')
     END IF

  CASE DEFAULT

     CALL Driver_abort('[eos_wlOneZone] Error: unsupported mode for Nuclear Eos')

  END SELECT

!   xEner = xEner - e_zeroPoint

  IF( varID == 0 ) THEN

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_P , P_T , xPres )

     CALL LogInterpolateDifferentiateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_E , E_T , xEner, xdEdX )

     xdEdT = xdEdX(iT_T)

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_S , S_T , xEntr )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Gm, Gm_T, xGamc )

     xCs2 = xGamc * xPres / xDens

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Xp, Xp_T, xXp )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Xn, Xn_T, xXn )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Xa, Xa_T, xXa )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Xh, Xh_T, xXh )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Ah, Ah_T, xAh )

   !   xAbar = 1.0 / xXp + 1.0 / xXn + 4.0 / xXa + xAh / xXh
     xAbar = 1.0 / (xXp + xXn +  xXa / 4.0 + xXh / xAh)

     xVar = 0.0

  ELSE IF( varID == 20 ) THEN

     ! this mode should only be used when grabbing variables for an already
     ! consistent thermodynamics.  Don't call this if thermo state has changed!

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Me, Me_T, xMu_e )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Mp, Mp_T, xMu_p )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Mn, Mn_T, xMu_n )

     xVar = xMu_e - xMu_n + xMu_p

  ELSE IF( varID > 0 .AND. varID <= nVariables ) THEN

     ! Is this the varID from Eos.h?
     ! Does this need conversion of varID to WL indices?
     OS_Var = EosNewTable % DV % Offsets(varID)

     ASSOCIATE( Var_T => EosNewTable % DV % Variables(varID) % Values )

     CALL LogInterpolateSingleVariable_3D_Custom_Point &
        ( xDens, xTemp, xYe, D_T, T_T, Y_T, OS_Var, Var_T, xVar )

     END ASSOCIATE

  ELSE

     CALL Driver_abort('[eos_wlOneZone] Error: Invalid varID')

  END IF

  END ASSOCIATE

  return
end subroutine eos_wlOneZone
