!!****ih* source/Grid/localAPI/gr_sbInterface
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
!!  gr_sbInterface
!!
!! SYNOPSIS
!!
!!  use gr_sbInterface
!!
!! DESCRIPTION
!!
!!  This is the header file for the solid body sub-unit that defines its
!!  interfaces.
!!***

#include "constants.h"
#include "Simulation.h"

module gr_sbInterface
  implicit none

  interface
     Subroutine gr_sbInit()
       implicit none
     End Subroutine gr_sbInit
  end interface

  interface
     Subroutine gr_sbCreateGroups()
       implicit none
     End Subroutine gr_sbCreateGroups
  end interface

  interface
     Subroutine gr_sbCreateParticles()
       implicit none
     End Subroutine gr_sbCreateParticles
  end interface

  interface
     Subroutine gr_sbGetProcBlock()
       implicit none
     End Subroutine gr_sbGetProcBlock
  end interface

  interface
     Subroutine gr_sbSendPosn()
       implicit none
     End Subroutine gr_sbSendPosn
  end interface

  interface
     Subroutine gr_sbStoreParticlesPerProc()
       implicit none
     End Subroutine gr_sbStoreParticlesPerProc
  end interface

  interface
     Subroutine gr_sbSendParticleCount()
       implicit none
     End Subroutine gr_sbSendParticleCount
  end interface

  interface
     Subroutine gr_sbSendParticles()
       implicit none
     End Subroutine gr_sbSendParticles
  end interface

  interface
     Subroutine gr_sbUpdateForces()
       implicit none
     End Subroutine gr_sbUpdateForces
  end interface

  interface
     Subroutine gr_sbSendForces()
       implicit none
     End Subroutine gr_sbSendForces
  end interface

  interface
     Subroutine gr_sbFinalize()
       implicit none
     End Subroutine gr_sbFinalize
  end interface

  interface
     Subroutine gr_sbDistributedForces()
       implicit none
     End Subroutine gr_sbDistributedForces
  end interface

  interface
     Subroutine gr_sbSendBoundBox()
       implicit none
     End Subroutine gr_sbSendBoundBox
  end interface

end module gr_sbInterface
