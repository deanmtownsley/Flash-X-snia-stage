!!****if* source/physics/ImBound/ImBoundMain/ib_annBuildTree
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use body file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!***

subroutine ib_annBuildTree(body)
   use iso_c_binding
   use ANN_mod
   use ANN_types_mod
   use ImBound_type, ONLY: ImBound_type_t

   implicit none

   class(ImBound_type_t), intent(INOUT)  :: body

      !! ANN tree variables
   integer :: p_i, rows, cols, i, cols_rc
   real, dimension(:, :), allocatable, target :: dset_data, dset_data_rc
      !! build the ANN tree for given body
   !=======================================================================
   body%kdTree = c_null_ptr
   body%kdTreeRC = c_null_ptr
   cols = body%dims ! need to (un)hard code body
   cols_rc = 1
   rows = body%numElems
   allocate (dset_data(rows, cols))
   allocate (dset_data_rc(rows, cols_rc))
   do p_i = 1, body%numElems
      dset_data(p_i, :) = (/body%elems(p_i)%xA, body%elems(p_i)%yA/)
      dset_data_rc(p_i, :) = body%elems(p_i)%yA
   end do
      !!! build the ann tree
      !! for finding nearest neighbors to get level-set value
   call ann_buildTree(rows, cols, c_loc(dset_data), body%kdTree)
      !! for finding nearest neighbors for ray casting
   call ann_buildTree(rows, cols_rc, c_loc(dset_data_rc), body%kdTreeRC)

end subroutine ib_annBuildTree
