!!****if* source/physics/Eos/EosMain/Eos_everywhere
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
!!  Eos_everywhere
!!
!!  For more details see the documentation of the NULL implementation
!! 
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData

#include "Eos.h"
#include "constants.h"
#include "Simulation.h"

subroutine Eos_everywhere(mode,gridDataStruct)

  use Driver_interface, ONLY : Driver_abort
  use Grid_interface,ONLY : Grid_getTileIterator, &
                            Grid_releaseTileIterator
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,     ONLY : Grid_tile_t
  use Eos_interface, ONLY : Eos, Eos_putData, Eos_getData
  use Eos_data, ONLY : eos_threadWithinBlock, eos_meshMe
  implicit none

  integer, intent(in) :: mode
  integer, optional, intent(IN) :: gridDataStruct

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  real, pointer :: solnData(:,:,:,:)

  real, allocatable :: eosData(:), massFraction(:)
  real, allocatable :: eosData3D(:), massFraction3D(:)

  logical, target, dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: lo, hi, iSize, jSize, kSize, vecLen, vLen3D
  integer :: iSizeMax, jSizeMax, kSizeMax
  integer :: ipres, idens, itemp, ieint, igamc, iabar, izbar, ientr, iekin, ielef
  integer :: ipres3D, idens3D, itemp3D, ieint3D, igamc3D, iabar3D, izbar3D, ientr3D, iekin3D, ielef3D

  integer :: ierr, istat, dataStruct
  integer :: kBlock


!! ---------------------------------------------------------------------------------
  ! Test calling arguments
#define DEBUG
#ifdef DEBUG
  ierr = 1
  select case (mode)
  case (MODE_DENS_PRES)
     ierr = 0
  case (MODE_DENS_TEMP)
     ierr = 0
  case (MODE_DENS_EI)
     ierr = 0
  case (MODE_EOS_NOP)
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
     call Driver_abort("[Eos_everywhere] "//&
          "invalid mode: must be MODE_DENS_PRES, MODE_DENS_TEMP, MODE_DENS_EI, or variants thereof, or MODE_EOS_NOP")
  end if
#endif

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *

  nullify(solnData)

  vecLen = 0
  iSizeMax = 1
  jSizeMax = 1
  kSizeMax = 1
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits
     iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
     iSizeMax = max(iSizeMax,iSize)
     jSizeMax = max(jSizeMax,jSize)
     kSizeMax = max(kSizeMax,kSize)
     vecLen = vecLen + iSize * jSize * kSize
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
  if (vecLen==0) return ! * Return immediately for empty range! (for efficiency and avoiding index range errors)

  vLen3D = iSizeMax * jSizeMax * kSizeMax
  allocate(massFraction3D(NSPECIES*vLen3D))
  allocate(eosData3D(EOS_NUM*vLen3D))
  allocate(massFraction(NSPECIES*vecLen))
  allocate(eosData(EOS_NUM*vecLen))

  ! solnData points to solution data in UNK (or other data structure).
  ! The length of the data being operated upon is determined from the range input array.

  if(present(gridDataStruct))then
     dataStruct=gridDataStruct
  else
     dataStruct=CENTER
  end if

  eosMask = .FALSE.

  ipres = (EOS_PRES-1)*vecLen
  idens = (EOS_DENS-1)*vecLen
  itemp = (EOS_TEMP-1)*vecLen
  ieint = (EOS_EINT-1)*vecLen
  igamc = (EOS_GAMC-1)*vecLen
  iabar = (EOS_ABAR-1)*vecLen
  izbar = (EOS_ZBAR-1)*vecLen
  ientr = (EOS_ENTR-1)*vecLen
  iekin = (EOS_EKIN-1)*vecLen
#ifdef EOS_YE
  ielef = (EOS_YE  -1)*vecLen
#endif

  ipres3D = (EOS_PRES-1)*vLen3D
  idens3D = (EOS_DENS-1)*vLen3D
  itemp3D = (EOS_TEMP-1)*vLen3D
  ieint3D = (EOS_EINT-1)*vLen3D
  igamc3D = (EOS_GAMC-1)*vLen3D
  iabar3D = (EOS_ABAR-1)*vLen3D
  izbar3D = (EOS_ZBAR-1)*vLen3D
  ientr3D = (EOS_ENTR-1)*vLen3D
  iekin3D = (EOS_EKIN-1)*vLen3D
#ifdef EOS_YE
  ielef3D = (EOS_YE  -1)*vLen3D
#endif

  kBlock = 0
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits
     kBlock = kBlock + 1

     iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
     vLen3D = iSize * jSize * kSize

     call Eos_getData(blkLimits,vLen3D,solnData,dataStruct,eosData3D,massFraction3D)

     lo = (kBlock-1)*vLen3D + 1
     hi =  kBlock   *vLen3D
     eosData(ipres+lo:ipres+hi) = eosData3D(ipres3D+1:ipres3D+vLen3D)
     eosData(idens+lo:idens+hi) = eosData3D(idens3D+1:idens3D+vLen3D)
     eosData(ieint+lo:ieint+hi) = eosData3D(ieint3D+1:ieint3D+vLen3D)
     eosData(itemp+lo:itemp+hi) = eosData3D(itemp3D+1:itemp3D+vLen3D)
     eosData(igamc+lo:igamc+hi) = eosData3D(igamc3D+1:igamc3D+vLen3D)
     eosData(iabar+lo:iabar+hi) = eosData3D(iabar3D+1:iabar3D+vLen3D)
     eosData(izbar+lo:izbar+hi) = eosData3D(izbar3D+1:izbar3D+vLen3D)
     eosData(ientr+lo:ientr+hi) = eosData3D(ientr3D+1:ientr3D+vLen3D)
     eosData(iekin+lo:iekin+hi) = eosData3D(iekin3D+1:iekin3D+vLen3D)
#ifdef EOS_YE
     eosData(ielef+lo:ielef+hi) = eosData3D(ielef3D+1:ielef3D+vLen3D)
#endif

     lo = (kBlock-1)*vLen3D*NSPECIES + 1
     hi =  kBlock   *vLen3D*NSPECIES
     massFraction(lo:hi) = massFraction3D(1:vLen3D*NSPECIES)

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Eos(mode,vecLen,eosData,massFraction,mask=eosMask)

  kBlock = 0
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits
     kBlock = kBlock + 1

     iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
     vLen3D = iSize * jSize * kSize

     lo = (kBlock-1)*vLen3D + 1
     hi =  kBlock   *vLen3D
     eosData3D(ipres3D+1:ipres3D+vLen3D) = eosData(ipres+lo:ipres+hi)
     eosData3D(idens3D+1:idens3D+vLen3D) = eosData(idens+lo:idens+hi)
     eosData3D(ieint3D+1:ieint3D+vLen3D) = eosData(ieint+lo:ieint+hi)
     eosData3D(itemp3D+1:itemp3D+vLen3D) = eosData(itemp+lo:itemp+hi)
     eosData3D(igamc3D+1:igamc3D+vLen3D) = eosData(igamc+lo:igamc+hi)
     eosData3D(iabar3D+1:iabar3D+vLen3D) = eosData(iabar+lo:iabar+hi)
     eosData3D(izbar3D+1:izbar3D+vLen3D) = eosData(izbar+lo:izbar+hi)
     eosData3D(ientr3D+1:ientr3D+vLen3D) = eosData(ientr+lo:ientr+hi)
     eosData3D(iekin3D+1:iekin3D+vLen3D) = eosData(iekin+lo:iekin+hi)
#ifdef EOS_YE
     eosData3D(ielef3D+1:ielef3D+vLen3D) = eosData(ielef+lo:ielef+hi)
#endif

     call Eos_putData(blkLimits,vLen3D,solnData,dataStruct,eosData3D)

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  deallocate(eosData)
  deallocate(massFraction)
  deallocate(eosData3D)
  deallocate(massFraction3D)

  return
end subroutine Eos_everywhere
