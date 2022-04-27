!!****if* source/physics/sourceTerms/Deleptonize/DeleptonizeMain/Deleptonize_init
!!
!! NAME
!!  
!!  Deleptonize_init
!!
!!
!! SYNOPSIS
!! 
!!  call Deleptonize_init()
!!
!!  
!! DESCRIPTION
!!
!!  Perform various initializations (apart from the problem-dependent ones)
!!  for the heat module.
!!
!!
!! ARGUMENTS
!!
!!  
!!
!!***
subroutine Deleptonize_init()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Deleptonize_data
  use Driver_interface, ONLY : Driver_getComm, Driver_getMype, Driver_abort
  use Logfile_interface, ONLY : Logfile_stampMessage, Logfile_stamp

#include "constants.h"
#include "Simulation.h"

  implicit none
  logical :: threadBlockListBuild, threadWithinBlockBuild
  integer :: ioStatus = 0

  call Driver_getComm(MESH_COMM,delep_meshComm)
  call Driver_getMype(MESH_COMM,delep_meshMe)

  call RuntimeParameters_get("useDeleptonize", useDeleptonize)
#ifdef FLASH_LEAKAGE
  call RuntimeParameters_get("useRadTrans", delep_useRadTrans)
#endif

  call RuntimeParameters_get("delep_Enu", delep_Enu)
  call RuntimeParameters_get("delep_rhoOne", delep_rhoOne)
  call RuntimeParameters_get("delep_rhoTwo", delep_rhoTwo)
  call RuntimeParameters_get("delep_yOne", delep_yOne)
  call RuntimeParameters_get("delep_yTwo", delep_yTwo)
  call RuntimeParameters_get("delep_yc", delep_yc)

!   call RuntimeParameters_get("useHeat", delep_useCool)
  call RuntimeParameters_get("useEntr", delep_useEntr)

  delep_minDens = 1.0e6

  delep_maxDens = 0.

  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadDelepBlockList", delep_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadDelepWithinBlock", delep_threadWithinBlock)

  if (delep_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Delep_init]')
     delep_threadBlockList = .false.
  end if
  if (delep_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Delep_init]')
     delep_threadWithinBlock = .false.
  end if

  ! for doing it from data in a file
  call RuntimeParameters_get("delep_fromFile", delep_fromFile)
  if (delep_fromFile) then
     call RuntimeParameters_get("delep_file", delep_file)
     open(unit=512, file=delep_file, status="old")
     delep_nrho = 0
     do while (ioStatus==0)
        delep_nrho = delep_nrho + 1
        read(512,*,IOSTAT=ioStatus) delep_rho(delep_nrho), delep_ye(delep_nrho)
     end do
     delep_nrho = delep_nrho-1
  end if
  
  return
end subroutine Deleptonize_init
