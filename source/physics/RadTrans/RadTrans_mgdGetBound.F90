!!****f* source/physics/RadTrans/RadTrans_mgdGetBound
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
!!  NAME 
!!
!!  RadTrans_mgdGetBound
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdGetBound( integer(IN) :: g,
!!                             real(OUT) :: b )
!!
!!  DESCRIPTION 
!!      This subroutine is used to access a particular radiation
!!      energy group boundary for MGD.
!!
!! ARGUMENTS
!!
!!      g : The boundary number, group g is bounded by g and g+1
!!      b : The boundary energy [ergs]
!! 
!!***
subroutine RadTrans_mgdGetBound(g, b)
  implicit none

  integer, intent(in) :: g
  real,    intent(out) :: b

  ! Stub implementation
  b = 0.0
  return

end subroutine RadTrans_mgdGetBound
