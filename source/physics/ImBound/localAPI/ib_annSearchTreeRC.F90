!!****if* source/physics/ImBound/localAPI/ib_annSearchTreeRC
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
!!***

subroutine ib_annSearchTreeRC(body, queryPt, nn, nnIdx, dists, eps)
   !
   use ImBound_type, ONLY: ImBound_type_t
   implicit none
   !
   class(ImBound_type_t), intent(IN)  :: body
   integer, intent(IN) :: nn
   ! query point
   real, dimension(:), target, intent(IN) :: queryPt
   ! distance of queryPt from neighbors
   real, dimension(:), allocatable, target, intent(OUT)  :: dists
   ! indices of nearest neighbors
   integer, dimension(:), allocatable, target, intent(OUT):: nnIdx
   real, intent(in) :: eps
end subroutine ib_annSearchTreeRC
