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
!! @brief  General parameters (non-array) for EOS Helmholtz
!!

module eos_helmData

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"  
#include "Eos_map.h"

  ! maximum number of iterations for the Newton loop to find T from e
  integer, save :: eos_maxNewton
 
  ! how accurately to do the Newton iteration to get T from e
  real, save :: eos_tol
  
  real, save :: eos_larget
  

  real, save :: eos_fluffDens

  integer,save :: eos_hfetInit

  ! force the iterative solver to leave inputs alone (always true in MODE_DENS_TEMP)
  logical, save :: eos_forceConstantInput

  ! Coulomb multiplier 
  real, save :: eos_coulombMult
  ! abort if pressures become negative
  logical, save :: eos_coulombAbort
 
  logical, save :: eos_useMultiSpecies
  integer,parameter :: EOSIMAX=541,EOSJMAX=201

  ! Minimum vecLen value to use OpenACC implementation of starkiller
  integer, save :: eos_vecLenACC

  real, save :: eos_tlo, eos_thi, eos_tstpi
  real, save :: eos_dlo, eos_dhi, eos_dstpi
  integer, parameter :: eos_t=1,eos_dt=2,eos_dtSqr=3,eos_dtInv=4,eos_dtSqrInv=5
  integer, parameter :: eos_d=1, eos_dd=2,eos_ddSqr=3,eos_ddInv=4,eos_ddSqrInv=5
  real,dimension(EOSJMAX,5),save :: eos_temps
  real,dimension(EOSIMAX,5),save :: eos_rhos
!..for the helmholtz free energy tables
!..for the pressure derivative with density tables
!..for the chemical potential tables
  !..for the number density tables
  integer, parameter :: EOST = 22, EOST_END=27
  real,save,dimension(EOSIMAX,EOSJMAX,EOST_END) :: eos_table
  integer, parameter:: eos_f=1,eos_fd=2,eos_ft=3,eos_fdd=4,&
             eos_ftt=5,eos_fdt=6,eos_fddt=7,eos_fdtt=8,eos_fddtt=9, & 
             eos_dpdf=10,eos_dpdfd=11,eos_dpdft=12,eos_dpdfdt=13, & 
             eos_ef=14,eos_efd=15,eos_eft=16,eos_efdt=17, & 
             eos_xf=18,eos_xfd=19,eos_xft=20,eos_xfdt=21
  integer, parameter :: eos_dpdfdd=22,eos_dpdftt=23,eos_efdd=24,eos_eftt=25,&
             eos_xfdd=26,eos_xftt=27

  real,  save ::  tempRow,denRow, abarRow,zbarRow

  !..totals and their derivatives
  real, save ::  ptotRow,dptRow,             &
       &                             etotRow,detRow, stotRow            
  real, save ::  dedRow, dstRow,dsdRow,dpdRow            
  real,  save ::  deaRow, dezRow  !Calhoun            


  !..electron-positron contributions -- most UNUSED and REMOVED
  real,  save :: pelRow, neRow, etaRow, detatRow
  !..derivative based quantities
  real,save ::    gamcRow
  real,  save ::    cpRow,cvRow
  
  logical, save :: eos_baprox13 = .false. 

!These variables must be threadprivate!!!
!$omp threadprivate(tempRow, denRow, etotRow, abarRow, zbarRow, gamcRow, ptotRow, &
!$omp deaRow, dezRow, detRow, dptRow, dpdRow, dedRow, pelRow, neRow, etaRow, detatRow, &
!$omp cvRow, cpRow, dstRow, dsdRow, stotRow)

end module eos_helmData
