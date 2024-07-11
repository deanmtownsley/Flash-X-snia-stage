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
!> @ingroup physics_Eos
!!
!! @brief general implementation of the API call for a computing
!!     equation of state on cells in a multidimensional array
!!
!! @subref{Eos_multiDim}
!!
!!

! solnData depends on the ordering on unk
!!REORDER(4): solnData
!!NOVARIANTS

#include "Eos.h"
#include "constants.h"
#include "Simulation.h"

subroutine Eos_multiDim(mode,range,solnData)

  use Driver_interface, ONLY : Driver_abort
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Eos_interface, ONLY : Eos_fillEosData, Eos_getFromEosData, Eos_vector

#include "FortranLangFeatures.fh"

  integer, intent(in) :: mode
  integer, dimension(2,MDIM), intent(in) :: range
  real, POINTER_INTENT_IN :: solnData(:,:,:,:)

  integer :: ierr
  integer :: i,j,k, vecLen

  integer :: n,m
  real :: pres, temp, dens, gamc, eint, entr, abar, zbar, ye

  real, allocatable, dimension(:,:) :: eosData, massFrac

  write(0,*) "Entering Eos_multiDim" !CF debug

  ierr = 1
  select case (mode)
  case (MODE_DENS_PRES)
  ierr = 0
  case (MODE_DENS_TEMP)
  ierr = 0
  case (MODE_DENS_EI)
  ierr = 0
  case (MODE_EOS_NOP,MODE_EOS_WRAPPERONLY)
  ierr = 0
  case (MODE_DENS_TEMP_ALL,MODE_DENS_TEMP_EQUI)
  ierr = 0
  case (MODE_DENS_EI_ALL,MODE_DENS_EI_SCATTER,MODE_DENS_EI_GATHER)
  ierr = 0
  case (MODE_DENS_EI_SELE_GATHER)
  ierr = 0
  case (MODE_DENS_ENTR)
  ierr = 0
  end select
  if(ierr /= 0) then
  call Driver_abort("[Eos_multiDim] "//&
  "invalid mode: must be MODE_DENS_PRES, MODE_DENS_TEMP, MODE_DENS_EI, or variants thereof, or MODE_EOS_NOP")
  end if
  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *

  vecLen = (range(HIGH,IAXIS)-range(LOW,IAXIS)+1)*&
           (range(HIGH,JAXIS)-range(LOW,JAXIS)+1)*&
           (range(HIGH,KAXIS)-range(LOW,KAXIS)+1)

  if (vecLen==0) return ! * Return immediately for empty IAXIS range! (for efficiency and avoiding index range errors)

  ! solnData points to solution data in UNK (or other data structure).
  ! The length of the data being operated upon is determined from the range input array.

  allocate(massFrac(NSPECIES,vecLen))
  allocate(eosData(vecLen, EOS_VARS))

  call Eos_fillEosData(range, solnData,1, vecLen, eosData)
  write(0,*) "Filled EosData" !CF debug
  n = 0
  do k=range(LOW,KAXIS),range(HIGH,KAXIS)
  do j=range(LOW,JAXIS),range(HIGH,JAXIS)
  do i=range(LOW,IAXIS),range(HIGH,IAXIS)
             n=n+1
             massFrac(1:NSPECIES,n) = solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)
  end do
  end do
  end do
  write(0,*) "Off to call Eos_vector. eosData contains ", eosData !CF debug
  call Eos_vector(mode,vecLen,eosData, massFrac)
  call Eos_getFromEosData( 1, vecLen, eosData, range, solnData)

  deallocate(massFrac)
  deallocate(eosData)

end subroutine Eos_multiDim
