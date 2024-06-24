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
!! @brief Initialization for the WeakLib EOS appropriate for 
!!        core-collapse supernova simulations. Read-in WeakLib Eos 
!!        table fully from a hdf5 file into EosOldTable and compress it 
!!        into a new Eos table EosNewTable for FLASH.
!!

SUBROUTINE eos_weaklibInit()

  USE eos_data, ONLY:                     &
      eos_type,                           &
      eos_meshMe,                         &
      eos_meshComm
  USE eos_wlData
  USE wlKindModule, ONLY: dp
  USE wlDependentVariablesModule, ONLY:   &
      DVIDType
  USE wlEquationOfStateTableModule, ONLY: &
      EquationOfStateTableType,           &
      DeAllocateEquationOfStateTable
  USE wlEOSIOModuleHDF, ONLY:             &
      ReadEquationOfStateTableHDF,        &
      MatchTableStructure,                &
      BroadcastEquationOfStateTableParallel
  USE wlEOSInversionModule, ONLY :        &
      InitializeEOSInversion
  USE RuntimeParameters_interface, ONLY:  &
      RuntimeParameters_get

  IMPLICIT NONE

#include "Eos.h"
#include "constants.h"

  INTEGER                        :: NewnVariables
  TYPE(DVIDType)                 :: NewDVID
  TYPE(EquationOfStateTableType) :: EosOldTable

  INTEGER                        :: ierr ! MPI error output
  LOGICAL                        :: Verbose

  eos_type = EOS_WL

  CALL RuntimeParameters_get('eos_file', eos_file)
  CALL RuntimeParameters_get('eos_wl_muShift', eos_wl_muShift)

  IF ( eos_meshMe == MASTER_PE ) THEN
 
    PRINT*, 'Reading in weaklib Eos table  ', eos_file
  
    CALL ReadEquationOfStateTableHDF( EosOldTable, eos_file )
  
    NewnVariables = 15
  
    NewDVID % iPressure                  = 1
    NewDVID % iEntropyPerBaryon          = 3
    NewDVID % iInternalEnergyDensity     = 2
    NewDVID % iElectronChemicalPotential = 5
    NewDVID % iProtonChemicalPotential   = 6
    NewDVID % iNeutronChemicalPotential  = 7
    NewDVID % iProtonMassFraction        = 8
    NewDVID % iNeutronMassFraction       = 9
    NewDVID % iAlphaMassFraction         = 10
    NewDVID % iHeavyMassFraction         = 11
    NewDVID % iHeavyChargeNumber         = 12
    NewDVID % iHeavyMassNumber           = 13
    NewDVID % iHeavyBindingEnergy        = 14
    NewDVID % iThermalEnergy             = 15
    NewDVID % iGamma1                    = 4 
  
    CALL MatchTableStructure( EosOldTable, EosNewTable, NewDVID, NewnVariables )
   
    CALL DeAllocateEquationOfStateTable( EosOldTable )
  
  END IF

  CALL BroadcastEquationOfStateTableParallel( EOSNewTable, MASTER_PE, eos_meshMe, ierr, eos_meshComm )
  eos_pointer => EosNewTable

  nVariables = EosNewTable % nVariables

  ASSOCIATE &
      ( iDtab  => EosNewTable % TS % Indices % iRho, &
        iTtab  => EosNewTable % TS % Indices % iT, &
        iYtab  => EosNewTable % TS % Indices % iYe, &
        iPtab  => EosNewTable % DV % Indices % iPressure, &
        iStab  => EosNewTable % DV % Indices % iEntropyPerBaryon, &
        iEtab  => EosNewTable % DV % Indices % iInternalEnergyDensity, &
        iMetab => EosNewTable % DV % Indices % iElectronChemicalPotential, &
        iMptab => EosNewTable % DV % Indices % iProtonChemicalPotential, &
        iMntab => EosNewTable % DV % Indices % iNeutronChemicalPotential, &
        iXptab => EosNewTable % DV % Indices % iProtonMassFraction, &
        iXntab => EosNewTable % DV % Indices % iNeutronMassFraction, &
        iXatab => EosNewTable % DV % Indices % iAlphaMassFraction, &
        iXhtab => EosNewTable % DV % Indices % iHeavyMassFraction, &
        iAhtab => EosNewTable % DV % Indices % iHeavyMassNumber, &
        iGmtab => EosNewTable % DV % Indices % iGamma1 )

  ASSOCIATE &
      ( Dtab  => EosNewTable % TS % States(iDtab) % Values, &
        Ttab  => EosNewTable % TS % States(iTtab) % Values, &
        Ytab  => EosNewTable % TS % States(iYtab) % Values, &
        Ptab  => EosNewTable % DV % Variables(iPtab ) % Values, &
        Stab  => EosNewTable % DV % Variables(iStab ) % Values, &
        Etab  => EosNewTable % DV % Variables(iEtab ) % Values, &
        Metab => EosNewTable % DV % Variables(iMetab) % Values, &
        Mptab => EosNewTable % DV % Variables(iMptab) % Values, &
        Mntab => EosNewTable % DV % Variables(iMntab) % Values, &
        Xptab => EosNewTable % DV % Variables(iXptab) % Values, &
        Xntab => EosNewTable % DV % Variables(iXntab) % Values, &
        Xatab => EosNewTable % DV % Variables(iXatab) % Values, &
        Xhtab => EosNewTable % DV % Variables(iXhtab) % Values, &
        Ahtab => EosNewTable % DV % Variables(iAhtab) % Values, &
        Gmtab => EosNewTable % DV % Variables(iGmtab) % Values, &
        OS_P  => EosNewTable % DV % Offsets(iPtab), &
        OS_S  => EosNewTable % DV % Offsets(iStab), &
        OS_E  => EosNewTable % DV % Offsets(iEtab), &
        OS_Me => EosNewTable % DV % Offsets(iMetab), &
        OS_Mp => EosNewTable % DV % Offsets(iMptab), &
        OS_Mn => EosNewTable % DV % Offsets(iMntab), &
        OS_Xp => EosNewTable % DV % Offsets(iXptab), &
        OS_Xn => EosNewTable % DV % Offsets(iXntab), &
        OS_Xa => EosNewTable % DV % Offsets(iXatab), &
        OS_Xh => EosNewTable % DV % Offsets(iXhtab), &
        OS_Ah => EosNewTable % DV % Offsets(iAhtab), &
        OS_Gm => EosNewTable % DV % Offsets(iGmtab) )

  IF ( eos_wl_muShift ) CALL eos_wlApplyMuShift( Mptab, Mntab, OS_Mp, OS_Mn )

  Verbose = ( eos_meshMe == MASTER_PE )
  CALL InitializeEOSInversion &
       ( Dtab, Ttab, Ytab, &
         10.0e0_dp**( Etab ) - OS_E, &
         10.0e0_dp**( Ptab ) - OS_P, &
         10.0e0_dp**( Stab ) - OS_S, &
         Verbose_Option = Verbose )

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: Dtab, Ttab, Ytab, &
    !$OMP          OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Ah, OS_Gm, &
    !$OMP          Ptab, Stab, Etab, Metab, Mptab, Mntab, Xptab, Xntab, Xatab, Xhtab, Ahtab, Gmtab )
#elif defined(WEAKLIB_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN( Dtab, Ttab, Ytab, &
    !$ACC         OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Ah, OS_Gm, &
    !$ACC         Ptab, Stab, Etab, Metab, Mptab, Mntab, Xptab, Xntab, Xatab, Xhtab, Ahtab, Gmtab )
#endif

  END ASSOCIATE

  END ASSOCIATE

  RETURN

  CONTAINS

     ! This should probably exist as a subroutine in a separate file so
     ! that any simulation could easily override it.
     SUBROUTINE eos_wlApplyMuShift( Mp_T, Mn_T, OS_Mp, OS_Mn )
        !For SFHo tables from
        !https://code.ornl.gov/astro/weaklib-tables/-/tree/master/SFHo/LowRes
        !up until git commit hash a36240ed
        !the neutron and proton chemical potentials
        !are tabulated subtracting the neutron-proton-mass difference
        !in order to use the same conventions as used in Chimera
        !to account for this, and get detailed balance factors
        !correct, we add the mass difference back in and then
        !add the SFHo reference masses for the chemical potential
        !(mn for Mn, mp for Mp)
        !For this renomalisation to the original SFHo tables,
        !we need to recalculate the offsets first

        REAL(dp), INTENT(inout), DIMENSION(:,:,:) :: Mp_T, Mn_T
        REAL(dp), INTENT(inout) :: OS_Mp, OS_Mn

        REAL(DP), PARAMETER :: neutron_mass = 939.56542052d0
        REAL(DP), PARAMETER :: proton_mass = 938.2720813d0
        REAL(DP), PARAMETER :: dmnp = 1.29333922d0
        REAL(DP) :: min_M, OS_M_new

        ! Apply the shift for proton chemical potential
        IF ( OS_Mp > 0.0d0 ) THEN
           min_M = -0.5d0 * OS_Mp
        ELSE
           min_M = MINVAL( 10.0d0**Mp_T )
        ENDIF
        OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
        Mp_T     = LOG10( 10.0d0**Mp_T - OS_Mp + proton_mass + dmnp + OS_M_new)
        OS_Mp    = OS_M_new

        ! Apply the shift for neutron chemical potential
        IF ( OS_Mn > 0.0d0 ) THEN
           min_M = -0.5d0 * OS_Mn
        ELSE
           min_M = MINVAL( 10.0d0**Mn_T )
        ENDIF
        OS_M_new = -2.0d0 * MIN( 0.0d0, min_M + proton_mass + dmnp )
        Mn_T     = LOG10( 10.0d0**Mn_T - OS_Mn + proton_mass + dmnp + OS_M_new)
        OS_Mn    = OS_M_new

     END SUBROUTINE eos_wlApplyMuShift

END SUBROUTINE eos_weaklibInit
