Module Orchestration_interfaceTypeDecl
   use,intrinsic :: iso_c_binding
   implicit none
#include "constants.h"
   type, public :: Orchestration_tileCPtrs_t
      type(c_ptr) :: ccBlkPtr
      type(c_ptr) :: fluxBlkPtrs(MDIM)
   end type Orchestration_tileCPtrs_t
   type, public :: Orchestration_tileCInts_t
      integer :: nCcComp
      integer :: nFluxComp
      integer :: loGC(MDIM), hiGC(MDIM)
      integer :: lo(MDIM),   hi(MDIM)
      integer :: ndim
      integer :: level
      integer :: gridIdxOrBlkId
      integer :: tileIdx
   end type Orchestration_tileCInts_t
   type, public :: Orchestration_tileCReals_t
      real :: deltas(MDIM)
   end type Orchestration_tileCReals_t

   type :: Orchestration_tileCInfo_t
      type(Orchestration_tileCPtrs_t) :: CPtrs
      type(Orchestration_tileCInts_t) :: CInts
      type(Orchestration_tileCReals_t) :: CReals
   end type Orchestration_tileCInfo_t
end Module Orchestration_interfaceTypeDecl
! Local Variables:
! f90-program-indent: 3
! indent-tabs-mode: nil
! End:
