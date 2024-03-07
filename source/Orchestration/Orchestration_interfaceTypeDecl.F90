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
