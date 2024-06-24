!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief Eos_getAbarZbar stub

!> @ingroup physics_Eos
!!
!! @brief return information of interest to other units
!!
!! @details
!! @anchor Eos_getAbarZbar_stub
!!
!! This routine returns values of Abar, Zbar, SumY, and/or Ye for one cell.
!! All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_getAbarZbar
!!
!!   @param solnVec : optional - the solution vector for one cell
!!
!!   @param abar    : optional - if present, will be filled with the average atomic mass
!!   @param zbar    : optional - if present, will be filled with the average ionization level
!!
!!   @param sumY    : optional - if present, will be filled with the inverse of the
!!                        average atomic mass
!!   @param Ye      : optional - if present, will be filled with Ye, which is defined
!!                        according to   Ye = zbar / abar
!!              
!!   @param massFrac : this is an optional argument which may be used when there is more 
!!              than one species in the simulation, as an alternative to providing
!!              the complete solution vector in solnVec
!!   
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



