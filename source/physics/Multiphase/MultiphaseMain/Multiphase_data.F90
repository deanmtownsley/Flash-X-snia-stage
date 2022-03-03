!!****if* source/physics/Multiphase/MultiphaseMain/Multiphase_data
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
!!  Multiphase_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE Multiphase_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the Multiphase module.
!!
!!***

module Multiphase_data

#include "Simulation.h"
#include "constants.h"

   logical, save :: mph_useMultiphase

   real, save :: mph_rhoGas
   real, save :: mph_muGas
   real, save :: mph_invWeber
   real, save :: mph_Tsat
   real, save :: mph_thcoGas
   real, save :: mph_CpGas

   real, save :: mph_crmx, mph_crmn

   integer, save :: mph_lsIt, mph_extpIt

   integer, save :: mph_meshMe
   integer, save :: mph_meshNumProcs
   integer, save :: mph_meshComm
   integer, save :: mph_nstep

   real, save :: mph_alfa
   real, save :: mph_rhoa
   real, save :: mph_gama

   integer, save :: mph_iVelFVar, mph_iMuCVar, mph_iRhoFVar, &
                    mph_iJumpVar, mph_iAlphaCVar

   integer, save :: mph_iTempVar, mph_iGfmVar, mph_iMdotVar, &
                    mph_iDivCVar

   real, save :: mph_Stefan, mph_invReynolds, mph_Prandtl

   real, save :: mph_iPropSmear

end module Multiphase_data
