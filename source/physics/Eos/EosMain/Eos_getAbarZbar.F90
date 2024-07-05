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
!! @brief Main implementation of the API call for fetching the values of abar and zbar
!! 
!! @subref{Eos_getAbarZbar}
!!
!!
!!
!!***

#include "Simulation.h"


subroutine Eos_getAbarZbar(solnVec,abar,zbar,sumY,Ye,massFrac)
  use Eos_interface, ONLY: Eos_getAbarZbarArraySection
  implicit none
  
  real, OPTIONAL,dimension(NUNK_VARS),intent(IN) :: solnVec
  real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
  real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac

  call Eos_getAbarZbarArraySection(1,solnVec,abar,zbar,sumY,Ye,massFrac)

end subroutine Eos_getAbarZbar


subroutine Eos_getAbarZbarArraySection(ifirstVar,solnVec,abar,zbar,sumY,Ye,massFrac)

#ifdef FLASH_MULTISPECIES
  use Multispecies_interface, ONLY: Multispecies_getSumInv, Multispecies_getSumFrac
#include "Multispecies.h"
#else
#if defined (SUMY_MSCALAR) && defined (YE_MSCALAR)
  use Driver_interface, ONLY : Driver_abort
#else
  use Eos_data, ONLY: eos_singleSpeciesZ, eos_singleSpeciesA
#endif
#endif

  implicit none
  
  integer,                            intent(IN) :: ifirstVar
  real, OPTIONAL,                     intent(IN) :: solnVec(ifirstVar:NUNK_VARS)
  real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
  real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac

  logical :: tryMscalarsFirst

#ifdef FLASH_MULTISPECIES
  integer :: specieStart, specieEnd
  real :: abarValue, zbarValue, abarInv, zbarFrac
#endif


  tryMscalarsFirst = .FALSE.
#if defined (SUMY_MSCALAR) && defined (YE_MSCALAR)
  ! Try deriving information from SumY and Ye mass scalars, first, if possible.
  ! It is assumed that Eos (or some other code) is keeping the contents of
  ! these variables up to date, whether they are derived from Multispecies
  ! or not. - KW
  if (present(solnVec)) then
     if (solnVec(SUMY_MSCALAR) > 0.0) then
        tryMscalarsFirst = .TRUE.
     end if
  end if

  if (tryMscalarsFirst) then
     if (present(abar)) abar = 1.0/solnVec(SUMY_MSCALAR)
     if (present(zbar)) zbar = solnVec(YE_MSCALAR)/solnVec(SUMY_MSCALAR)
     if (present(sumY)) sumY = solnVec(SUMY_MSCALAR)
     if (present(Ye))   Ye   = solnVec(YE_MSCALAR)
     return                     ! We are done, RETURN !
  end if
#endif


#ifdef FLASH_MULTISPECIES
  specieStart = 1
  specieEnd = NSPECIES

  if (present(massFrac)) then
     call Multispecies_getSumInv(A, abarInv ,massFrac(specieStart:specieEnd))
  else
     call Multispecies_getSumInv(A, abarInv ,solnVec(SPECIES_BEGIN:SPECIES_END))
  end if
  abarValue = 1.e0 / abarInv

  if (present(massFrac)) then
     call Multispecies_getSumFrac(Z,zbarFrac,massFrac(specieStart:specieEnd))
  else
     call Multispecies_getSumFrac(Z,zbarFrac,solnVec(SPECIES_BEGIN:SPECIES_END))
  end if
  zbarValue = abarValue*zbarFrac

  if (present(abar)) abar = abarValue
  if (present(zbar)) zbar = zbarValue
  if (present(sumY)) sumY = 1.0/abarValue
  if (present(Ye))   Ye   = zbarValue/abarValue

#else
#if defined (SUMY_MSCALAR) && defined (YE_MSCALAR)
  if (.not. present(solnVec)) then
     ! If SUMY/YE used, the solution vector must be provided,
     ! otherwise there is no way to compute abar or zbar.
     call Driver_abort("[Eos_getAbarZbar] Cannot compute abar or zbar without solution data")
  end if
  
  if (present(abar)) abar = 1.0/solnVec(SUMY_MSCALAR)
  if (present(zbar)) zbar = solnVec(YE_MSCALAR)/solnVec(SUMY_MSCALAR)
  if (present(sumY)) sumY = solnVec(SUMY_MSCALAR)
  if (present(Ye))   Ye   = solnVec(YE_MSCALAR)
#else
  if (present(abar)) abar = eos_singleSpeciesA
  if (present(zbar)) zbar = eos_singleSpeciesZ
  if (present(sumY)) sumY = 1.0/eos_singleSpeciesA
  if (present(Ye))   Ye   = eos_singleSpeciesZ/eos_singleSpeciesA
#endif
#endif

  return
end subroutine Eos_getAbarZbarArraySection



