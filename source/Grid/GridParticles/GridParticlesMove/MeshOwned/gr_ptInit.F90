!!****if* source/Grid/GridParticles/gr_ptInit
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
!!  gr_ptInit
!!
!! SYNOPSIS
!!
!!  call gr_ptInit()
!!
!! DESCRIPTION
!!
!!  Initialize values for data in the module gr_ptData.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

#include "Simulation.h"
#ifdef FLASH_PIMG
#include "ProtonImaging.h"
#endif
#ifdef FLASH_PEMI
#include "ProtonEmission.h"
#endif
#ifdef FLASH_EDEP
#include "EnergyDeposition.h"
#endif

subroutine gr_ptInit()
  use gr_ptData, ONLY :   gr_ptMaxPerProc,&
                          gr_ptRemove,gr_ptRemoveAlgo,&
                          gr_ptNumToReduce,gr_ptSieveCheckFreq,&
                          gr_ptLogLevel,gr_ptKeepLostParticles,gr_ptMaxVirtualCount
#ifndef FLASH_GRID_UG
  use gr_ptData, ONLY : gr_ptMaxPerProcUpperThresh, gr_ptMaxPerProcLowerThresh, &
       gr_ptMaxPerProcBlockFactor,gr_ptMaxPerProcBlockNoFuzz, gr_ptRefineOnPtMaxPerProc
  use Grid_data, ONLY : gr_meshMe, gr_refineOnParticleCount
#endif  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage, Logfile_stamp
  use Driver_interface, ONLY : Driver_abort

  implicit none 
  integer, save :: maxPerProc, propCount
 

  call RuntimeParameters_get("gr_ptRemove",gr_ptRemove)
  call RuntimeParameters_get("gr_ptRemoveAlgo",gr_ptRemoveAlgo)
  call RuntimeParameters_get("gr_ptNumToReduce",gr_ptNumToReduce)
  call RuntimeParameters_get("gr_ptSieveCheckFreq",gr_ptSieveCheckFreq)
  call RuntimeParameters_get("keepLostParticles",gr_ptKeepLostParticles)

  call RuntimeParameters_get("pt_logLevel",gr_ptLogLevel)


  gr_ptMaxPerProc=1
  propCount=1

#ifdef NPART_PROPS
  if (NPART_PROPS > 1) then
     call RuntimeParameters_get("pt_maxPerProc",maxPerProc)
     propCount=max(NPART_PROPS,propCount)                       ! DANGER
     gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)          ! OK -> 2nd buffer dimension
  endif
#endif
  gr_ptMaxVirtualCount=gr_ptMaxPerProc

#ifdef RAY_ATTR_COUNT
  call RuntimeParameters_get("ed_maxRayCount",maxPerProc)
  propCount=max(RAY_ATTR_COUNT,propCount)                       ! DANGER
  gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)             ! OK -> 2nd buffer dimension
#endif

#ifdef PROTON_ATTRCOUNT
  call RuntimeParameters_get("pi_maxProtonCount",maxPerProc)
  propCount=max(PROTON_ATTRCOUNT,propCount)                     ! DANGER
  gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)             ! OK -> 2nd buffer dimension
#endif

#ifdef EMPROTON_ATTRCOUNT
  call RuntimeParameters_get("pem_maxProtonCount",maxPerProc)
  propCount=max(EMPROTON_ATTRCOUNT,propCount)                   ! DANGER
  gr_ptMaxPerProc = max(maxPerProc,gr_ptMaxPerProc)             ! OK -> 2nd buffer dimension
#endif
!
!
!    ...The following check has been added, such that the FLASH application
!       stops smoothly with an informative error massage, rather than an
!       uncontrolled code crash.
!
!
#ifdef NPART_PROPS
  if (NPART_PROPS > 1 .and. NPART_PROPS /= propCount) then
      call Logfile_stamp (NPART_PROPS,"[gr_ptInit]: Value of NPART_PROPS is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abort ("[gr_ptInit]: NPART_PROPS must match propCount (see Logfile).")
  endif
#endif

#ifdef RAY_ATTR_COUNT
  if (RAY_ATTR_COUNT /= propCount) then
      call Logfile_stamp (RAY_ATTR_COUNT,"[gr_ptInit]: Value of RAY_ATTR_COUNT is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abort ("[gr_ptInit]: RAY_ATTR_COUNT must match propCount (see Logfile).")
  endif
#endif

#ifdef PROTON_ATTRCOUNT
  if (PROTON_ATTRCOUNT /= propCount) then
      call Logfile_stamp (PROTON_ATTRCOUNT,"[gr_ptInit]: Value of PROTON_ATTRCOUNT is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abort ("[gr_ptInit]: PROTON_ATTRCOUNT must match propCount (see Logfile).")
  endif
#endif


#ifdef EMPROTON_ATTRCOUNT
  if (EMPROTON_ATTRCOUNT /= propCount) then
      call Logfile_stamp (EMPROTON_ATTRCOUNT,"[gr_ptInit]: Value of EMPROTON_ATTRCOUNT is = ")
      call Logfile_stamp (propCount,"[gr_ptInit]: Value of propCount is = ")
      call Logfile_stampMessage ("Different 1st dimensions for Destination and Source Buffer!")
      call Driver_abort ("[gr_ptInit]: EMPROTON_ATTRCOUNT must match propCount (see Logfile).")
  endif
#endif


#ifndef FLASH_GRID_UG
  ! additional RPs for maxPerProc-based refinement criteria - KW
  call RuntimeParameters_get("gr_ptMaxPerProcUpperThresh",gr_ptMaxPerProcUpperThresh)
  call RuntimeParameters_get("gr_ptMaxPerProcLowerThresh",gr_ptMaxPerProcLowerThresh)
  call RuntimeParameters_get("gr_ptMaxPerProcBlockFactor",gr_ptMaxPerProcBlockFactor)
  call RuntimeParameters_get("gr_ptMaxPerProcBlockNoFuzz",gr_ptMaxPerProcBlockNoFuzz)
  call RuntimeParameters_get("gr_ptRefineOnPtMaxPerProc", gr_ptRefineOnPtMaxPerProc)
  if (gr_ptRefineOnPtMaxPerProc .AND. .NOT. gr_refineOnParticleCount) then
     call Logfile_stampMessage( &
          "WARNING: Ignoring gr_ptRefineOnPtMaxPerProc because RP refine_on_particle_count is .FALSE.")
     gr_ptRefineOnPtMaxPerProc = .FALSE.
  end if
#endif
  return
end subroutine gr_ptInit
