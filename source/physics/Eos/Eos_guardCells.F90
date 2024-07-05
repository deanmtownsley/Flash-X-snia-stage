!> @copyright Copyright 2024 UChicago Argonne, LLC and contributors
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
!! @brief Eos_guardCells stub

!> @ingroup physics_Eos
!!
!! @brief applies Eos to guardcells that are at fine-coarse boundaries
!!
!! @details
!! @anchor Eos_guardCells_stub
!!
!!
!!  @param eosmode : determines which variables are used as Eos input variables.
!!                   The valid values are MODE_DENS_EI (where density and internal
!!                   energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!                   MODE_DENS_TEMP (density and temperature are inputs).
!!                   These quantities are defined in constants.h.
!!                   The argument is passed unchanged and unexamined to Eos_multiDim calls.
!!
!!  @param Uin : pointer to data structure for one block worth for state variables
!!
!!  @param corners : indicates whether Eos should be called on corner
!!                   guard cells (i.e., diagonal guard cells)
!!  @param layers  : the number of guard cells to be included along each dimension
!!  @param skipSrl : whether to skip guard cell regions that are coming from
!!                   neighboring blocks at the same refinement (or from boundary
!!                   conditions) and thus have not undergone interpolation or
!!                   restrictions.
!!
!!
!!***

subroutine Eos_guardCells(eosMode, Uin,corners,layers,skipSrl)

#include "constants.h"
  implicit none

  integer,intent(IN) :: eosMode
  real,dimension(:,:,:,:),pointer::Uin
  
  logical,intent(IN) :: corners
  integer,dimension(MDIM),optional, intent(IN) :: layers
  logical,optional, intent(IN) :: skipSrl

end subroutine Eos_guardCells
