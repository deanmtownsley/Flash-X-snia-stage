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
SUBROUTINE eos_weaklib(mode,vecLen,xPres, xTemp, xDens, xGamc, xEner, xEntr,xAbar,xZbar,xYe, massFrac,derivs)

  USE Driver_interface, ONLY : Driver_abort
  USE eos_weaklib_inter, ONLY: eos_weaklib_short

  IMPLICIT NONE

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"

  !     Arguments
  INTEGER, INTENT(in) :: mode, vecLen
  REAL, INTENT(inout), DIMENSION(vecLen) :: xDens,xTemp,xPres, xGamc,xEntr,xAbar,xZbar,xEner, xYe
  REAL, OPTIONAL,INTENT(in), DIMENSION(vecLen*NSPECIES) :: massFrac
  real,OPTIONAL, DIMENSION(EOS_VARS+1:EOS_NUM),INTENT(in)::derivs

  REAL, DIMENSION(vecLen) :: xCs2, xA, xZ
  INTEGER :: xMode, err

  err = 0

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

  RETURN

END SUBROUTINE eos_weaklib
