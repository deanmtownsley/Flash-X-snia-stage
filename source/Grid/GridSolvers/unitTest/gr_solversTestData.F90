!!****if* source/Grid/GridSolvers/unitTest/gr_solversTestData
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
!!  gr_solversTestData
!!
!! SYNOPSIS
!!  use gr_solversTestData
!!
!! DESCRIPTION
!!
!!  Defines some data items that are private to
!!  the GridSolver unit test implementation.
!!
!!  
!!***

Module gr_solversTestData 
  
  implicit none

  ! store error norm tolerances
  real, save :: gr_testTolL2, gr_testTolLinf

end Module gr_solversTestData
