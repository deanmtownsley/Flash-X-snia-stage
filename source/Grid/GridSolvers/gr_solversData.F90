!!****if* source/Grid/GridSolvers/gr_solversData
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
!!  gr_solversData
!!
!! SYNOPSIS
!!  use gr_solversData, ONLY: gr_solversDbgContext
!!
!! DESCRIPTION
!!
!!  Defines storage for some data items that are common to
!!  several GridSolvers implementation.
!!
!!  
!!***

Module gr_solversData 
  
  use gr_interfaceTypeDecl, ONLY: gr_solversDbgContext_t
  
  implicit none

  ! Structure that holds context information on the current operation,
  ! for debugging
  type(gr_solversDbgContext_t), save :: gr_solversDbgContext
  
end Module gr_solversData
