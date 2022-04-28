!!****f* source/physics/Eos/Eos_getAbarZbar
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
!! NAME
!!
!!  Eos_getAbarZbar
!! 
!! SYNOPSIS
!!
!!  call Eos_getAbarZbar(
!!             OPTIONAL,real(IN)  :: solnVec(NUNK_VARS),
!!             OPTIONAL,real(OUT) :: abar,
!!             OPTIONAL,real(OUT) :: zbar,
!!             OPTIONAL,real(OUT) :: sumY,
!!             OPTIONAL,real(OUT) :: Ye,
!!             OPTIONAL,real(IN)  :: massFrac(:) )
!!
!!
!!
!! DESCRIPTION
!!
!! Eos_getAbarZbar gets Abar, Zbar, SumY, and/or Ye for one cell.
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   solnVec : optional - the solution vector for one cell
!!
!!   abar    : optional - if present, will be filled with the average atomic mass
!!   zbar    : optional - if present, will be filled with the average ionization level
!!
!!   sumY    : optional - if present, will be filled with the inverse of the
!!                        average atomic mass
!!   Ye      : optional - if present, will be filled with Ye, which is defined
!!                        according to   Ye = zbar / abar
!!              
!!   massFrac : this is an optional argument which may be used when there is more 
!!              than one species in the simulation, as an alternative to providing
!!              the complete solution vector in solnVec
!!   
!!
!!  EXAMPLE 
!!
!!  NOTES
!!
!!      All arguments are optional since there are different ways to
!!      provide input to the routine and different responses a caller
!!      may want.
!!
!!
!!      This interface is defined in Fortran Module 
!!      Eos_interface. All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_getAbarZbar
!!
!!
!!  SEE ALSO
!!
!!     Eos
!!
!!
!!***

#include "Simulation.h"


subroutine Eos_getAbarZbar(solnVec,abar,zbar,sumY,Ye,massFrac)

  implicit none
  
!#include "Eos.h"
!#include "Eos_map.h"
!#include "constants.h"
  
  real, OPTIONAL,dimension(NUNK_VARS),intent(IN) :: solnVec
  real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
  real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac


  if (present(abar)) abar = 0.0
  if (present(zbar)) zbar = 0.0
  if (present(sumY)) sumY = 0.0
  if (present(Ye))   Ye   = 0.0

  return
end subroutine Eos_getAbarZbar



