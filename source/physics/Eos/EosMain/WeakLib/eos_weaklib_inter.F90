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
!! @brief   Interpolation module.
!!   eos_weaklib_short: calls WeakLib interpolate routine
!!   and interpolate for given density, temperature, and
!!   electron fraction.
!!
!! AUTHOR & DATE
!!   R. Chu, Dept. Phys. & Astronomy
!!   U. Tennesee, Knoxville
!!   10/17/2018


MODULE eos_weaklib_inter

  USE eos_wlData, ONLY: EosNewTable, nVariables
  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
      LogInterpolateSingleVariable
  USE wlEOSInversionModule, ONLY : &
      ComputeTemperatureWith_DEY, &
      ComputeTemperatureWith_DSY, &
      ComputeTemperatureWith_DPY

  IMPLICIT NONE

  PUBLIC :: eos_weaklib_short

CONTAINS

  SUBROUTINE eos_weaklib_short &
             ( xDens,xTemp,xYe,xEner,xPres,xEntr,xGamc, &
               mode, err )

    REAL(dp), INTENT(inout), DIMENSION(:) :: xDens,xTemp,xYe,&
                                             xEner,xPres,xEntr, &
                                             xGamc
    INTEGER, INTENT(in)                   :: mode
    INTEGER, INTENT(out)                  :: err

    ! local variables
    INTEGER, DIMENSION(:), ALLOCATABLE  :: error_flag
    LOGICAL,  DIMENSION(nVariables) :: MaskVar
    INTEGER :: vecLen
    INTEGER :: ii

    err = 0
    vecLen = SIZE(xDens)

    ALLOCATE( error_flag(vecLen) )

    ASSOCIATE(   TS => EosNewTable % TS, &
                 DV => EosNewTable % DV, &
            Indices => EosNewTable % DV % Indices )

    SELECT CASE( mode )

      CASE( 0 ) ! MODE_DENS_EI

        CALL ComputeTemperatureWith_DEY &
             ( xDens, xEner, xYe, &
               TS % States(1) % Values, &
               TS % States(2) % Values, &
               TS % States(3) % Values, &
               DV % Variables(Indices % iInternalEnergyDensity) &
               % Values(:,:,:), &
               DV % Offsets(Indices % iInternalEnergyDensity), &
               xTemp, UseInitialGuess_Option = .TRUE., &
               Error_Option = error_flag )

        IF( ANY( error_flag /= 0 ) ) THEN
          WRITE(*,*) 'ComputeTemperatureWith_DEY ERROR:'
          WRITE(*,'(4A15,A10)') 'xDens','xEner','xYe','xTemp','Error'
          DO ii = 1,vecLen
            IF( error_flag(ii) /= 0 ) THEN
              WRITE(*,'(4ES15.4,I10)') &
                xDens(ii), xEner(ii), xYe(ii), xTemp(ii), error_flag(ii)
            END IF
          END DO
          CALL Driver_abort('ComputeTemperatureWith_DEY ERROR')
        END IF

      CASE( 2 ) ! MODE_DENS_ENTR

        CALL ComputeTemperatureWith_DSY &
             ( xDens, xEntr, xYe, &
               TS % States(1) % Values, &
               TS % States(2) % Values, &
               TS % States(3) % Values, &
               DV % Variables(Indices % iEntropyPerBaryon) &
               % Values(:,:,:), &
               DV % Offsets(Indices % iEntropyPerBaryon), &
               xTemp, UseInitialGuess_Option = .TRUE., &
               Error_Option = error_flag )

        IF( ANY( error_flag /= 0 ) ) THEN
          WRITE(*,*) 'ComputeTemperatureWith_DSY ERROR:'
          WRITE(*,'(4A15,A10)') 'xDens','xEntr','xYe','xTemp','Error'
          DO ii = 1,vecLen
            IF( error_flag(ii) /= 0 ) THEN
              WRITE(*,'(4ES15.4,I10)') &
                xDens(ii), xEntr(ii), xYe(ii), xTemp(ii), error_flag(ii)
            END IF
          END DO
          CALL Driver_abort('ComputeTemperatureWith_DSY ERROR')
        END IF

      CASE( 4 ) ! MODE_DENS_PRES

        CALL ComputeTemperatureWith_DPY &
             ( xDens, xPres, xYe, &
               TS % States(1) % Values, &
               TS % States(2) % Values, &
               TS % States(3) % Values, &
               DV % Variables(Indices % iPressure) &
               % Values(:,:,:), &
               DV % Offsets(Indices % iPressure), &
               xTemp, UseInitialGuess_Option = .TRUE., &
               Error_Option = error_flag )

        IF( ANY( error_flag /= 0 ) ) THEN
          WRITE(*,*) 'ComputeTemperatureWith_DPY ERROR:'
          WRITE(*,'(4A15,A10)') 'xDens','xPres','xYe','xTemp','Error'
          DO ii = 1,vecLen
            IF( error_flag(ii) /= 0 ) THEN
              WRITE(*,'(4ES15.4,I10)') &
                xDens(ii), xPres(ii), xYe(ii), xTemp(ii), error_flag(ii)
            END IF
          END DO
          CALL Driver_abort('ComputeTemperatureWith_DPY ERROR')
        END IF

     END SELECT

     CALL LogInterpolateSingleVariable &
          ( xDens, xTemp, xYe, TS % States(1) % Values, &
            TS % States(2) % Values, TS % States(3) % Values, &
            DV % Offsets(Indices % iInternalEnergyDensity), &
            DV % Variables(Indices % iInternalEnergyDensity) &
               % Values(:,:,:), &
            xEner )

     CALL LogInterpolateSingleVariable &
          ( xDens, xTemp, xYe, TS % States(1) % Values, &
            TS % States(2) % Values, TS % States(3) % Values, &
            DV % Offsets(Indices % iPressure), &
            DV % Variables(Indices % iPressure) &
               % Values(:,:,:), &
            xPres )

     CALL LogInterpolateSingleVariable &
          ( xDens, xTemp, xYe, TS % States(1) % Values, &
            TS % States(2) % Values, TS % States(3) % Values, &
            DV % Offsets(Indices % iEntropyPerBaryon), &
            DV % Variables(Indices % iEntropyPerBaryon) &
               % Values(:,:,:), &
            xEntr )

     CALL LogInterpolateSingleVariable &
          ( xDens, xTemp, xYe, TS % States(1) % Values, &
            TS % States(2) % Values, TS % States(3) % Values, &
            DV % Offsets(Indices % iGamma1), &
            DV % Variables(Indices % iGamma1) &
               % Values(:,:,:), &
            xGamc )

    END ASSOCIATE

    DEALLOCATE( error_flag )

    RETURN

  END SUBROUTINE eos_weaklib_short

END MODULE eos_weaklib_inter
