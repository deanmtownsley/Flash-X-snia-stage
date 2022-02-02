!!****if* source/physics/Eos/EosMain/Multigamma/eos_initMgamma
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
!!  eos_initMgamma
!!
!! 
!! SYNOPSIS
!!
!!  call eos_initMgamma()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars needed
!!  by the multigamma EOS from the runtime parameters and physical
!!  constants facilities. This version is for use when multiple species
!!  are present. The gamma values for different species are obtained from
!!  the Multispecies unit, initialized in Simulation_initSpecies.F90 .
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Gamma law Eos for multiple
!!   species with different abundances.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might over write these values with the 
!!   flash.par values for your specific run.  
!!
!!  NOTES
!!
!!  Gamma law Eos defines two mesh-based parameters GAMC_VAR and GAME_VAR in Simulation.h
!!
!!
!!***

subroutine eos_initMgamma()

  use Eos_data, ONLY : eos_gasConstant, eos_smalle, eos_eintSwitch, eos_type
  use eos_mgammaData, ONLY : eos_gc, eos_gammam1j, eos_ggprodj, eos_ggprodinvj, eos_gam1invj

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getProperty

 
  implicit none

#include "Simulation.h"
#include "Eos.h"
#include "Multispecies.h"


  integer :: numCells, istat
  integer :: ispecies


  do ispecies = 1, NSPECIES
     call Multispecies_getProperty(SPECIES_BEGIN + ispecies - 1, GAMMA,  eos_gc(ispecies))
  end do

  ! Note that these are all ARRAYS of size NSPECIES
  eos_gammam1j   = 1. / (eos_gc - 1.)
  eos_ggprodj    = eos_gammam1j * eos_gasConstant
  eos_ggprodinvj = 1. / eos_ggprodj
  eos_gam1invj   = 1. / eos_gammam1j
  eos_type = EOS_MGAM

  return
end subroutine eos_initMgamma
