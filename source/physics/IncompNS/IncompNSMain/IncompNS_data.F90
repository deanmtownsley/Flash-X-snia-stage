!!****if* source/physics/IncompNS/IncompNSMain/IncompNS_data
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
!!  IncompNS_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE IncompNS_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the Incmp_Navier_Stokes module.
!!
!!***
 
 
module IncompNS_data

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"


  logical, save :: ins_useIncompNS

  integer, save :: ins_isgs
  integer, save :: ins_cflflg
  integer, save :: ins_intschm
  integer, save :: ins_nstep
  integer, save :: ins_meshMe
  integer, save :: ins_meshNumProcs
  integer, save :: ins_meshComm

  real, save :: ins_cfl
  real, save :: ins_sigma
  real, save :: ins_invReynolds
  real, save :: ins_dtspec

  integer, parameter :: ins_rkstep = 3

  real, save :: ins_alf(ins_rkstep)
  real, save :: ins_gam(ins_rkstep)
  real, save :: ins_rho(ins_rkstep)
  real, save :: ins_alfa, ins_gama, ins_rhoa

  integer, save :: ins_intschm_type

  logical, save :: ins_prescorr
  real, save :: ins_prescoeff

  logical, save :: ins_restart

  real, save, dimension(LOW:HIGH,MDIM) :: ins_globalDomain
  integer,save,dimension(2,MDIM) :: ins_domainBC

  logical, save :: ins_predcorrflg = .FALSE.

  integer, save :: ins_prol_method

  ! Container for old timesteps
  real, dimension(ins_rkstep), save :: ins_vardt(-ins_rkstep:0)

  ! Gravitational acceleration:
  real, save    :: ins_gravX, ins_gravY, ins_gravZ 

  ! Fixed pressure gradient, constant mass vars for channel simulations: 
  ! Constant mass Flow in the Z direction.
  real, save    :: ins_dpdx,ins_dpdy,ins_dpdz
  real, save    :: ins_rhoGas, ins_muGas

  integer, dimension(6), save :: ins_pressureBC_types
  real, dimension(2,6), save  :: ins_pressureBC_values
  real, save :: ins_poisfact

  real, save :: ins_mindiv,ins_maxdiv

  real, save :: ins_outflowVel(LOW:HIGH,MDIM)

end module IncompNS_data
