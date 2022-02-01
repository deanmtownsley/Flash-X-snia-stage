!!****f* source/physics/RadTrans/RadTrans_mgdSetBc
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
!!  RadTrans_mgdSetBc
!!
!! SYNOPSIS
!!
!!  call RadTrans_mgdSetBc(integer(in) :: ig,
!!                         integer, optional(in) :: bctypes,
!!                         real, optional(in) :: bcvalues,
!!                         integer, optional(in) :: f,
!!                         integer, optional(in) :: bctype,
!!                         real, optional(in) :: bcvalue)
!!
!! DESCRIPTION
!!
!! Set the boundary condition for a specific group. Can be invoked in
!! one of two modes:
!! 
!! 1. bctypes/bcvalues present: This allows you to easily set the BC on all
!!    of the faces of the domain
!! 2. f/bctype/bcvalue present: This allows you to easily set the BC on a
!!    specific face of the domain
!!
!!
!! ARGUMENTS
!!
!!   ig : energy group number
!!   bctypes : boundary condition type on each boundary
!!   bcvalues : boundary condition value (for dirichlet) on each boundary
!!   f : boundary number
!!   bctype : type of bondary for face f
!!   bcvalue : dirichlet value for face f
!!
!!***

subroutine RadTrans_mgdSetBc(ig, bcTypes, bcValues, f, bcType, bcValue)
  implicit none
  integer, intent(in) :: ig

  integer, optional, intent(in) :: bcTypes(6)
  real, optional, intent(in) :: bcValues(6)

  integer, optional, intent(in) :: f
  integer, optional, intent(in) :: bcType
  real, optional, intent(in) :: bcValue

  ! Stub Implementation

end subroutine RadTrans_mgdSetBc
