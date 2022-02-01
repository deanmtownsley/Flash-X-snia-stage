!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/PttoPt/gr_ptMapInit
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
!!  gr_ptMapInit
!!
!! SYNOPSIS
!!
!!  gr_ptMapInit()
!!
!! DESCRIPTION
!!
!!  Initialize values for all data in the module gr_ptMapData
!!
!! ARGUMENTS
!!
!!
!!***

subroutine gr_ptMapInit()

  use gr_ptMapData, ONLY : gr_ptSmearLen, gr_ptRecvSpecifier, &
       gr_ptRecvSpecifierTmp, gr_ptRecvTotalTmp, gr_ptRecvTotal
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abort
  use Grid_data, ONLY : gr_meshNumProcs

  implicit none 
  integer :: error

  call RuntimeParameters_get("smearLen", gr_ptSmearLen)

  if (gr_ptSmearLen < 0) then
     call Driver_abort("[gr_ptMapInit]: Variable smearLen must be at least 0")
  end if


  allocate(gr_ptRecvSpecifier(0:gr_meshNumProcs-1), & 
       gr_ptRecvSpecifierTmp(0:gr_meshNumProcs-1), &
       gr_ptRecvTotalTmp(0:gr_meshNumProcs-1), &
       gr_ptRecvTotal(0:gr_meshNumProcs-1), STAT=error)
  if (error /= 0) then
     call Driver_abort("[gr_ptMapInit]: Memory cannot be allocated!")
  end if

end subroutine gr_ptMapInit
