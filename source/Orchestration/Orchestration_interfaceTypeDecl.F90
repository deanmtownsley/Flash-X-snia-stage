!> @copyright Copyright 2024 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! @file
!> @ingroup Orchestration
!!
!! @brief Public interface for Orchestration unit type declarations
!!
!! @details
!! A Flash-X Fortran module that declares types that are used in some
!! of the public subroutine interfaces declared in Orchestration_interface.
!! of all routine's in the Orchestration unit that are part of this unit's public
!! interface.  The version of this module at the top unit level (stub
!! level) should only contain enough information to allow other units to
!! have syntactically correct references to public API subroutines of the
!! Orchestration unit, even when a concrete implementation of the
!! Orchestration unit is not included in a build.
!! This top-level version of the module will be overridden by a version
!! that contains the necessary details if a concrete implementation of the
!! unit is included.
Module Orchestration_interfaceTypeDecl
  use,intrinsic :: iso_c_binding, ONLY: C_INT, C_PTR
  implicit none

  integer, parameter, public :: MILHOJA_INT = C_INT
  abstract interface
     !> Fortran interface of the runtime's task function.
     !! C_threadId - unique zero-based index of runtime thread calling this
     !!              routine
     !! C_dataItemPtr - C pointer to Grid DataItem to which the task
     !!                 function should be applied
     subroutine milhoja_runtime_taskFunction(C_threadId, C_dataItemPtr) bind(c)
       import
       implicit none
       integer(MILHOJA_INT), intent(IN), value :: C_threadId
       type(C_PTR),          intent(IN), value :: C_dataItemPtr
     end subroutine milhoja_runtime_taskFunction
  end interface

  type, public :: Orchestration_tileCPtrs_t
  end type Orchestration_tileCPtrs_t
  type, public :: Orchestration_tileCInts_t
  end type Orchestration_tileCInts_t
  type, public :: Orchestration_tileCReals_t
  end type Orchestration_tileCReals_t

  type :: Orchestration_tileCInfo_t
  end type Orchestration_tileCInfo_t
end Module Orchestration_interfaceTypeDecl
