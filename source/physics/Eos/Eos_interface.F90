
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



Module Eos_interface

  implicit none

# include "Eos.h"
# include "Simulation.h"
# include "constants.h"
#include "FortranLangFeatures.fh"
  
  interface
     subroutine Eos_guardCells(eosMode, solnData,corners,layers,skipSrl,blockDesc)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       integer,intent(IN) :: eosMode
       real,dimension(:,:,:,:),pointer::solnData
       logical,intent(IN) :: corners
       integer,dimension(MDIM),optional,intent(IN) :: layers
       logical,optional, intent(IN) :: skipSrl
       type(Grid_tile_t),optional,intent(IN) :: blockDesc
     end subroutine Eos_guardCells
  end interface
  
  interface Eos_wrapped
     subroutine Eos_wrapped(mode,range,solnData)
       integer, intent(in) :: mode
       integer, dimension(2,MDIM), intent(in) :: range
       real, POINTER_INTENT_IN :: solnData(:,:,:,:)
     end subroutine Eos_wrapped
  end interface

  interface
     subroutine Eos(mode, vecLen, eosData,  massFrac, mask, vecBegin,vecEnd,diagFlag)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask     
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer, optional, INTENT(out)    :: diagFlag
     end subroutine Eos
  end interface
  
  interface Eos_init
     subroutine Eos_init()
     end subroutine Eos_init
  end interface
  
  interface Eos_finalize
     subroutine Eos_finalize()
     end subroutine Eos_finalize
  end interface

  interface
     subroutine Eos_getParameters(eintSwitch,inputsAreUnchanged,inputTempIsGuess,constantGammaC,&
          inputMassFracNeeded,smalle,smallE1,smallE2,smallE3)
       real,OPTIONAL,intent(OUT) :: eintSwitch
       logical,OPTIONAL,intent(OUT) :: inputsAreUnchanged
       logical,OPTIONAL,intent(OUT) :: inputTempIsGuess
       logical,OPTIONAL,intent(OUT) :: constantGammaC
       logical,OPTIONAL,intent(OUT) :: inputMassFracNeeded
       real,OPTIONAL,intent(OUT) :: smalle
       real,OPTIONAL,intent(OUT) :: smallE1,smallE2,smallE3
     end subroutine Eos_getParameters
  end interface

  interface Eos_unitTest
     subroutine Eos_unitTest(fileUnit, perfect)
       integer, intent(in) :: fileUnit
       logical, intent(out) :: perfect
     end subroutine Eos_unitTest
  end interface

  interface Eos_getAbarZbar
     subroutine Eos_getAbarZbar(solnVec,abar,zbar,sumY,Ye,massFrac)
       implicit none
       real, OPTIONAL,dimension(NUNK_VARS),intent(IN) :: solnVec
       real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
       real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac
     end subroutine Eos_getAbarZbar
     subroutine Eos_getAbarZbarArraySection(ifirstVar,solnVec,abar,zbar,sumY,Ye,massFrac)
       implicit none
       integer,                            intent(IN) :: ifirstVar
       real, OPTIONAL,                     intent(IN) :: solnVec(ifirstVar:NUNK_VARS)
       real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
       real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac
     end subroutine Eos_getAbarZbarArraySection
  end interface

  interface
     subroutine Eos_logDiagnostics(force)
       implicit none
       logical, intent(IN) :: force
     end subroutine Eos_logDiagnostics
  end interface

  interface
     subroutine Eos_sendOutputData()
       implicit none
     end subroutine Eos_sendOutputData
  end interface

end Module Eos_interface
  
