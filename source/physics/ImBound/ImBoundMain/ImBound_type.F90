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
      real :: xA, yA, zA
      real :: xB, yB, zB
      real :: xC, yC, zC
      real :: xCenter, yCenter, zCenter
      real :: xNorm, yNorm, zNorm
   end type ib_elem

   type ImBound_type_t
      integer :: numElems
      integer :: dims
      type(ib_elem), dimension(:), allocatable :: elems
      real :: velx, vely, velz
      real :: thetax, thetay, thetaz
      real :: boundBox(LOW:HIGH, IAXIS:KAXIS)
      type(c_ptr) :: kdTree = c_null_ptr
      type(c_ptr) :: kdTreeRC = c_null_ptr
   end type ImBound_type_t

end module ImBound_type
