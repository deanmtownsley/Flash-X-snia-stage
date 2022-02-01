!!****if* source/physics/Eos/EosMain/Helmholtz/eos_helmData
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
!!  eos_helmData
!!
!! 
!! SYNOPSIS
!!
!!  use eos_helmData
!!
!! DESCRIPTION
!!
!!  General parameters (non-array) for EOS Helmholtz
!!
!! ARGUMENTS
!!
!!
!!*** 

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
 
  ! Flag for whether to use starkiller helmholtz or old flash version
  logical, save :: eos_useStarkiller
  integer,parameter :: EOSIMAX=541,EOSJMAX=201

  ! Minimum vecLen value to use OpenACC implementation of starkiller
  integer, save :: eos_vecLenACC

  real, save :: eos_tlo, eos_thi, eos_tstpi
  real, save :: eos_dlo, eos_dhi, eos_dstpi
  real,dimension(EOSJMAX),save :: eos_dt,eos_dtSqr,eos_dtInv,eos_dtSqrInv,eos_t
  real,dimension(EOSIMAX),save :: eos_dd,eos_ddSqr,eos_ddInv,eos_ddSqrInv,eos_d

!..for the helmholtz free energy tables
!..for the pressure derivative with density tables
!..for the chemical potential tables
!..for the number density tables
  real,save,dimension(EOSIMAX,EOSJMAX) :: eos_f,eos_fd, eos_ft,eos_fdd,&
                                          eos_ftt,eos_fdt,eos_fddt,&
                                          eos_fdtt, eos_fddtt, & 
                                          eos_dpdf,eos_dpdfd,eos_dpdft,&
                                          eos_dpdfdd,eos_dpdftt,eos_dpdfdt,&
                                          eos_ef,eos_efd,eos_eft,eos_efdd,&
                                          eos_eftt,eos_efdt, & 
                                          eos_xf,eos_xfd,eos_xft,eos_xfdd,&
                                          eos_xftt,eos_xfdt
end module eos_helmData
