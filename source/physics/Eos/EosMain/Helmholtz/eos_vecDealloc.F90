!!****if* source/physics/Eos/EosMain/Helmholtz/eos_vecDealloc
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
!!  eos_vecDealloc
!! 
!! SYNOPSIS
!!
!!  eos_vecDealloc()
!!
!! DESCRIPTION
!!
!!  Deallocate vectors for Helmholtz which were previously set
!!     up with eos_vecAlloc.  Used only for NONFIXEDBLOCKSIZE routines
!!
!! ARGUMENTS
!!
!!
!!*** 

subroutine eos_vecDealloc()

#include "constants.h"
#include "Simulation.h"

  use eos_vecData, ONLY: tempRow, denRow, abarRow, zbarRow, &
       ptotRow, dptRow, dpdRow, etotRow,&
       detRow, dedRow,  deaRow, dezRow, &
       dstRow, dsdRow,  stotRow, &
       pelRow, neRow, etaRow, detatRow, cpRow, cvRow, gamcRow

  implicit none

#ifndef FIXEDBLOCKSIZE

  deallocate(tempRow)
  deallocate(denRow)
  deallocate(abarRow)
  deallocate(zbarRow)

!..totals and their derivatives
  deallocate(ptotRow)
  deallocate(dptRow)
  deallocate(dpdRow)
! UNUSED  deallocate(dpaRow)
! UNUSED  deallocate(dpzRow)
  deallocate(etotRow)
  deallocate(detRow)
  deallocate(dedRow)
  deallocate(deaRow)
  deallocate(dezRow)
! UNUSED  deallocate(deaRow)
! UNUSED  deallocate(dezRow)
  deallocate(stotRow)
  deallocate(dstRow)
  deallocate(dsdRow)
! UNUSED  deallocate(dsaRow)
! UNUSED  deallocate(dszRow)

  
!..radiation contributions 
!UNUSED  deallocate(pradRow)
!UNUSED  deallocate(eradRow)
!UNUSED  deallocate(sradRow)

!..ion contributions 

!UNUSED  deallocate(pionRow)
!UNUSED  deallocate(eionRow)
!UNUSED  deallocate(sionRow)
!UNUSED  deallocate(xniRow)

!..electron-positron contributions -- most  UNUSED and REMOVED
!  see eos_fixedVecData if you want to see the original names
  deallocate(pelRow)
  deallocate(neRow)
  deallocate(etaRow)
  deallocate(detatRow)

!..coulomb contributions - all UNUSED and REMOVED

!..thermodynamic consistency checks; maxwell relations -- all UNUSED and REMOVED

!..derivative based quantities
  deallocate(cpRow)
  deallocate(cvRow)
  deallocate(gamcRow)

#endif

end subroutine eos_vecDealloc
