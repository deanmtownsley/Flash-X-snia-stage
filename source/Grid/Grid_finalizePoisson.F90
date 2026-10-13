!!****if* source/Grid/GridSolvers/Multipole_new/Grid_solvePoisson
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
!!  Grid_finalizePoisson
!!
!! SYNOPSIS
!!
!!  Grid_finalizePoisson (integer(in)    :: iSoln,
!!                     integer(in)    :: iSrc, 
!!                     integer(6)(in) :: bcTypes,
!!                     real(2,6)(in)  :: bcValues,
!!                     real(inout)    :: poisFact)
!!
!! DESCRIPTION
!!
!!  Poisson solver routine.  This module implements the multipole
!!  summation method for isolated boundary problems.  Periodic
!!  problems are not supported; the multipole method should be
!!  used only for approximately spherical matter distributions.
!!
!! ARGUMENTS
!!
!!  iSoln          : index to variable containing potential
!!  iSrc           : index to variable containing density
!!  bcTypes        : boundary types along various faces,
!!                   only used in verifying that they are isolated
!!  bcValues       : the values to boundary conditions, currently not used
!!  poisFact       : factor to be used in calculation
!!
!!***

subroutine Grid_finalizePoisson (iSoln,                   &
                              iSrc,                    &
                              bcTypes,                 &
                              bcValues,                &
                              poisFact)

  implicit none


  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes (6)
  real,    intent(in)    :: bcValues (2,6)
  real,    intent(inout) :: poisFact
!
  return
end subroutine Grid_finalizePoisson
