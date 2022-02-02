!!****ih* source/Grid/localAPI/gr_interfaceTypeDecl
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
!!  gr_interfaceTypeDecl
!!
!! SYNOPSIS
!!
!!  use gr_interfaceTypeDecl
!!
!! DESCRIPTION
!!
!!  Contains derived data type declarations used by gr_interface.  
!!  We need to place the derived data types in a module in order that 
!!  different subroutines have access to the same type declaration.
!!
!!***

module gr_interfaceTypeDecl

#include "constants.h"
#include "Simulation.h"

  implicit none
  type BlockRegion_t
     integer :: numNegh  !No. neighbors at this guard cell region. 0 = external boundary.
     integer, dimension(BLKNO:TYPENO, 2**(NDIM-1)) :: details
  end type BlockRegion_t

  type AllBlockRegions_t
     type (BlockRegion_t), dimension(1:3, 1:1+2*K2D, 1:1+2*K3D) :: regionInfo
  end type AllBlockRegions_t

  type gr_solversDbgContext_t
     integer :: component
     integer :: group
     integer :: libErrCode
     integer :: flashErrCode    ! 1 for ERROR, 2 for INFO
     integer :: retriable       ! 0 for NO, 1 for YES
  end type gr_solversDbgContext_t

end module gr_interfaceTypeDecl
