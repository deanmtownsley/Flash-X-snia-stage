!!****if* source/physics/ImBound/localAPI/ib_annSearchTree
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

subroutine ib_annSearchTree(body, queryPt, annElems, annIdx)
   !
   use iso_c_binding
   use ANN_mod
   use ANN_types_mod
   use ImBound_type, ONLY: ImBound_type_t
   implicit none
   !
   type(ImBound_type_t), intent(IN)  :: body
   integer, intent(IN) :: annElems
   real, dimension(:), target, intent(IN) :: queryPt ! query point
   integer, dimension(:), target, intent(OUT):: annIdx ! indices of nearest neighbors

   ! local variables
   real :: eps = 0.
   real, dimension(:), allocatable, target :: annDists

   allocate (annDists(annElems))

   call ann_kSearch(c_loc(queryPt), body%dims, annElems, c_loc(annIdx), c_loc(annDists), eps, body%kdTree)

   deallocate(annDists)

end subroutine ib_annSearchTree
