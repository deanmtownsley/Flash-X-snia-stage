!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/gr_ptMapInit
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

  use gr_ptMapData, ONLY : gr_ptSmearLen
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abort
  implicit none 


  call RuntimeParameters_get("smearLen", gr_ptSmearLen)


  if (gr_ptSmearLen < 0) then
     call Driver_abort("Variable smearLen must be at least 0")
  end if

end subroutine gr_ptMapInit
