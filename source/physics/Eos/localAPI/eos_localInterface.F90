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
     subroutine eos_initGamma()
     end subroutine eos_initGamma
  end interface

  interface 
     subroutine eos_initMgamma()
     end subroutine eos_initMgamma
  end interface

  interface 
     subroutine eos_initHelmholtz()
     end subroutine eos_initHelmholtz
  end interface


  interface
     subroutine eos_initWeaklib()
     end subroutine eos_initWeaklib
  end interface

  interface 
     subroutine eos_externalComputeAbarZbar(solnScalars, abarData, zbarData)
       implicit none
       real, dimension(:,:), intent(in) :: solnScalars
       real, dimension(:), intent(out)  :: abarData, zbarData
     end subroutine eos_externalComputeAbarZbar
  end interface

end module eos_localInterface
