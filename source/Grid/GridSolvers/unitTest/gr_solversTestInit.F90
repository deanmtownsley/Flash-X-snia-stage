!!****if* source/Grid/GridSolvers/unitTest/gr_solversTestInit
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
!!  gr_solversInit
!!
!! 
!! SYNOPSIS
!!
!!  call gr_solversTestInit()
!!
!!
!! DESCRIPTION
!!
!!  This routine initializes data used by tests
!!  of grid solvers.
!!
!!***

subroutine gr_solversTestInit()
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use gr_solversTestData, ONLY: gr_testTolL2, gr_testTolLinf
  implicit none 

  call RuntimeParameters_get ("gr_testTolL2",   gr_testTolL2)
  call RuntimeParameters_get ("gr_testTolLinf", gr_testTolLinf)

end subroutine gr_solversTestInit     
