!!****ih* source/physics/Eos/localAPI/eos_localInterface
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
!!     eos_localInterface
!!
!! SYNOPSIS
!!     use eos_localInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! the Eos unit.
!!
!!***
!!NOVARIANTS

module eos_localInterface
  implicit none
#include "Eos.h"
#include "Simulation.h"


  interface
     subroutine eos_weaklib(mode, vecLen, eosData, massFrac, mask)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
     end subroutine eos_weaklib
  end interface

    interface
     subroutine eos_idealGamma(mode, pres, temp, dens, gamc, eint, entr, abar, zbar , massFrac, mask)
       implicit none
       integer, INTENT(in) :: mode
       real,INTENT(inout) :: pres, temp, dens, gamc, eint, entr, abar, zbar 
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES)    :: massFrac
     end subroutine eos_idealGamma
  end interface
 

  interface
     subroutine eos_helmSpecies(mode, pres, temp, dens, gamc, eint, entr, abar, zbar , massFrac, mask)
       integer, INTENT(in) :: mode
       real,INTENT(inout) :: pres, temp, dens, gamc, eint, entr, abar, zbar 
       real, optional, INTENT(in),dimension(NSPECIES)    :: massFrac
       logical, optional,target, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask
     end subroutine eos_helmSpecies
  end interface

  
  interface 
     subroutine eos_idealGammaInit()
     end subroutine eos_idealGammaInit
  end interface

  interface 
     subroutine eos_helmSpeciesInit()
     end subroutine eos_helmSpeciesInit
  end interface


  interface
     subroutine eos_weaklibInit()
     end subroutine eos_weaklibInit
  end interface

  interface
     subroutine eos_hybridInit()
     end subroutine eos_hybridInit
  end interface

  interface 
     subroutine eos_externalComputeAbarZbar(solnScalars, abarData, zbarData)
       implicit none
       real, dimension(:,:), intent(in) :: solnScalars
       real, dimension(:), intent(out)  :: abarData, zbarData
     end subroutine eos_externalComputeAbarZbar
  end interface

end module eos_localInterface
