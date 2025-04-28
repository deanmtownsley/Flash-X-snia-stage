!!****if* source/physics/ImBound/ImBoundMain/ib_annDestroyTree
!! NOTICE
!!  Copyright 2024 UChicago Argonne, LLC and contributors
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

subroutine ib_annDestroyTree(body)
   use iso_c_binding
   use ANN_mod
   use ANN_types_mod
   use ImBound_type, ONLY: ImBound_type_t

   implicit none
   type(ImBound_type_t), intent(INOUT)  :: body

   call ann_destroyTree(body%kdTree)
end subroutine ib_annDestroyTree
