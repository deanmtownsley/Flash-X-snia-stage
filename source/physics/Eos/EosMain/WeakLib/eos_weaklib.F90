!!****if* source/physics/Eos/EosMain/WeakLib/eos_weaklib
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
!! AUTHOR & DATE 
!!   R. Chu, Dept. Phys. & Astronomy
!!   U. Tennesee, Knoxville
!!   10/17/2018
!!
!! DESCRIPTION
!!
!!  Main driver routine for the weaklib EOS.
!!  It is called by Eos.F90. 
!!  Call eos_weaklib_short() to do the interpolation.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!
!!***
SUBROUTINE eos_weaklib(mode,vecLen,eosData,massFrac,mask)

  USE Driver_interface, ONLY : Driver_abort
  USE eos_weaklib_inter, ONLY: eos_weaklib_short

  IMPLICIT NONE

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"

  !     Arguments
  INTEGER, INTENT(in) :: mode, vecLen
  ! EOS_NUM defined in Eos.h
  REAL, INTENT(inout), DIMENSION(vecLen*EOS_NUM) :: eosData
  REAL, OPTIONAL,INTENT(in), DIMENSION(vecLen*NSPECIES) :: massFrac
  ! must correspond to dimensions of Eos_wrapped
  LOGICAL,OPTIONAL, DIMENSION(EOS_VARS+1:EOS_NUM),INTENT(in)::mask

  INTEGER :: pres, temp, dens, gamc, eint, game, abar, zbar, entr, &
             elef
  REAL, DIMENSION(vecLen) :: xDens,xTemp,xYe,xEner,xPres,xEntr,&
                             xCs2, xGamc, xA, xZ
  INTEGER :: xMode, err

  err = 0

  ! These integers are indexes into the lowest location in UNK that contain the appropriate variable
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  eint = (EOS_EINT-1)*vecLen   
  temp = (EOS_TEMP-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen   
  abar = (EOS_ABAR-1)*vecLen   
  zbar = (EOS_ZBAR-1)*vecLen   
  entr = (EOS_ENTR-1)*vecLen

  elef = (EOS_YE-1)*vecLen

  SELECT CASE(mode)
    CASE(MODE_DENS_EI)
       xMode = 0
    CASE(MODE_DENS_TEMP)
       xMode = 1
    CASE(MODE_DENS_ENTR)
       xMode = 2
    CASE(MODE_DENS_PRES)
       xMode = 4
    CASE default
       CALL Driver_abort&
               ('[Eos] Error: unsupported mode for Nuclear Eos')
  END SELECT

  !      Crank the EOS on the pipes filled above, then fill the FLASH arrays
  !      with the thermodynamic quantities returned by the WeakLib EOS.


  xDens = eosData(dens+1:dens+vecLen)
  xTemp = eosData(temp+1:temp+vecLen)
  xYe   = eosData(elef+1:elef+vecLen)
!!!  xYe   = eosData(zbar+1:zbar+vecLen)/eosData(abar+1:abar+vecLen)
  xEner = eosData(eint+1:eint+vecLen)
  xPres = eosData(pres+1:pres+vecLen)
  xEntr = eosData(entr+1:entr+vecLen)

!!$  PRINT*, 'eos_weaklib  before interpolation  Z', &
!!$          MAXVAL(eosData(zbar+1:zbar+vecLen) ), &
!!$          'A', &
!!$          MAXVAL(eosData(abar+1:abar+vecLen) ), &
!!$          'Ye', MAXVAL(xYe)

  IF( MAXVAL(xDens) < TINY(1.d0) ) THEN
    PRINT*, ' eos_weaklib.F90 line 90 : xDens = zero '
    PRINT*, ' MAXVAL(xDens) ', MAXVAL(xDens), 'MINVAL(xDens) ',MINVAL(xDens) 
    CALL Driver_abort("[EOS] problem with weaklib EOS")
  END IF

  CALL eos_weaklib_short&
       ( xDens,xTemp,xYe,xEner,xPres,xEntr,xGamc,&
         xMode,err )

      IF (err /= 0) THEN
        PRINT*,"ERROR: Printing from eos_weaklib.f90 line 119, inside routine eos_weaklib"
        CALL Driver_abort("[EOS] problem with weaklib EOS")
      ENDIF

!!$  IF( MINVAL(xEner) < TINY(1.d0) ) THEN
!!$    PRINT*, ' eos_weaklib.F90 line 118 : xEner = zero '
!!$    PRINT*, ' MAXVAL(xEner) ', MAXVAL(xEner), 'MINVAL(xEner) ',MINVAL(xEner)
!!$    CALL Driver_abort("[EOS] problem with weaklib EOS")
!!$  END IF
 
      eosData(dens+1:dens+vecLen) = xDens
      eosData(temp+1:temp+vecLen) = xTemp
      eosData(pres+1:pres+vecLen) = xPres
      eosData(eint+1:eint+vecLen) = xEner
      eosData(gamc+1:gamc+vecLen) = xGamc
      eosData(entr+1:entr+vecLen) = xEntr
      eosData(elef+1:elef+vecLen) = xYe 

!!$  PRINT*, 'eos_weaklib  after interpolation  Z', &
!!$          MAXVAL(eosData(zbar+1:zbar+vecLen) ), &
!!$          'A', &
!!$          MAXVAL(eosData(abar+1:abar+vecLen) ), &
!!$          'Ye', MAXVAL(xZ/xA)
  
  RETURN

END SUBROUTINE eos_weaklib
