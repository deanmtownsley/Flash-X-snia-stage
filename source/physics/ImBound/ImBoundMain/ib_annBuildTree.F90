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

   type(ImBound_type_t), intent(INOUT)  :: body

   ! ANN tree variables
   integer :: panelIndex, rowSize, colSize
   real, dimension(:, :), allocatable, target :: annDset

   ! build the ANN tree for given body
   body%kdTree = c_null_ptr
   colSize = body%dims ! need to (un)hard code body
   rowSize = body%numElems

   allocate (annDset(rowSize, colSize))

   do panelIndex = 1, body%numElems
      annDset(panelIndex, :) = body%elems(panelIndex)%center(1:body%dims)
   end do

   ! build the ann tree
   ! for finding nearest neighbors to get level-set value
   call ann_buildTree(rowSize, colSize, c_loc(annDset), body%kdTree)

   deallocate (annDset)

end subroutine ib_annBuildTree
