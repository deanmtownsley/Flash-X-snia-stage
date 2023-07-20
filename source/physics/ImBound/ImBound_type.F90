!!****if* source/physics/ImBound/ImBound_bodyType
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
!!  ImBound_bodyType
!!
!!
!! SYNOPSIS
!!
!!  MODULE ImBound_bodyType()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the ImBound module.
!!
!!***
#include "constants.h"

module ImBound_type
   use iso_c_binding
   implicit none

   type ib_elem
      real :: xA, yA
      real :: xB, yB
      real :: xCenter, yCenter
   end type ib_elem

   type ImBound_type_t
      integer :: numElems
      integer :: dim
      type(ib_elem), dimension(:), allocatable :: elems
      real :: velx, vely
      real :: boundBox(LOW:HIGH, IAXIS:KAXIS)
      type(c_ptr) :: kdTree = c_null_ptr
      type(c_ptr) :: kdTreeRC = c_null_ptr
   contains
      procedure, public :: buildTree
      procedure, public :: searchTree
      procedure, public :: searchTreeRc
   end type ImBound_type_t

contains

   !! This builds ANN tree for a given body data
   subroutine buildTree(this)
      use iso_c_binding
      use ANN_mod
      use ANN_types_mod
      !! ANN tree variables
      integer :: p_i, rows, cols, i, cols_rc
      real, dimension(:, :), allocatable, target :: dset_data, dset_data_rc
      class(ImBound_type_t), intent(INOUT)  :: this
      !! build the ANN tree for given body
      !=======================================================================
      this%kdTree = c_null_ptr
      this%kdTreeRC = c_null_ptr
      cols = this%dim ! need to (un)hard code this
      cols_rc = 1
      rows = this%numElems
      allocate (dset_data(rows, cols))
      allocate (dset_data_rc(rows, cols_rc))
      do p_i = 1, this%numElems
         dset_data(p_i, :) = (/this%elems(p_i)%xA, this%elems(p_i)%yA/)
         dset_data_rc(p_i, :) = this%elems(p_i)%yA
      end do
      !!! build the ann tree
      !! for finding nearest neighbors to get level-set value
      call ann_buildTree(rows, cols, c_loc(dset_data), this%kdTree)
      !! for finding nearest neighbors for ray casting
      call ann_buildTree(rows, cols_rc, c_loc(dset_data_rc), this%kdTreeRC)
   end subroutine

   !! search the ANN tree for finding nearest neighbors (NN) for dist value

   !! Inputs:

   ! this - the body pointer
   ! queryPt - point for which to find NN
   ! dim - dimension
   ! nn - #NN needed
   ! eps - distance tolerance

   !! Outputs:

   ! nnIdx - indices of NN
   ! dists - distance of NN

   subroutine searchTree(this, queryPt, nn, nnIdx, dists, eps)
      use ANN_mod
      use ANN_types_mod
      class(ImBound_type_t), intent(IN)  :: this
      integer, intent(IN) :: nn
      ! query point
      real, dimension(:), target, intent(IN) :: queryPt
      ! distance of queryPt from neighbors
      real, dimension(:), allocatable, target, intent(OUT)  :: dists
      ! indices of nearest neighbors
      integer, dimension(:), allocatable, target, intent(OUT):: nnIdx
      real, intent(in) :: eps
      allocate (dists(nn))
      allocate (nnIdx(nn))
      call ann_kSearch(c_loc(queryPt), this%dim, nn, c_loc(nnIdx), c_loc(dists), eps, this%kdTree)
   end subroutine

   !! search the ANN tree for finding nearest neighbors (NN) for dist sign (ray casting)

   !! Inputs:

   ! this - the body pointer
   ! queryPt - point for which to find NN
   ! dim - dimension
   ! nn - #NN needed
   ! eps - distance tolerance

   !! Outputs:

   ! nnIdx - indices of NN
   ! dists - distance of NN

   subroutine searchTreeRC(this, queryPt, nn, nnIdx, dists, eps)
      use ANN_mod
      use ANN_types_mod
      class(ImBound_type_t), intent(IN)  :: this
      integer, intent(IN) :: nn
      ! query point
      real, dimension(:), target, intent(IN) :: queryPt
      ! distance of queryPt from neighbors
      real, dimension(:), allocatable, target, intent(OUT)  :: dists
      ! indices of nearest neighbors
      integer, dimension(:), allocatable, target, intent(OUT):: nnIdx
      real, intent(in) :: eps
      allocate (dists(nn))
      allocate (nnIdx(nn))
      call ann_kSearch(c_loc(queryPt), this%dim, nn, c_loc(nnIdx), c_loc(dists), eps, this%kdTreeRC)
   end subroutine

   !type ImBound_type_t
   !   integer :: numElems
   !   type(ib_elem), dimension(:), allocatable :: elems
   !end type ImBound_type_t

end module ImBound_type
