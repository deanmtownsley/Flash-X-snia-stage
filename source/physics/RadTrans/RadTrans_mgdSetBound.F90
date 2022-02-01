!!****f* source/physics/RadTrans/RadTrans_mgdSetBound
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
!!  RadTrans_mgdSetBound
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdSetBound( integer(IN) :: g,
!!                             real(IN) :: b )
!!
!!  DESCRIPTION 
!!      This subroutine is used to set a particular radiation
!!      energy group boundary for MGD.
!!
!! ARGUMENTS
!!
!!      g : The boundary number, group g is bounded by g and g+1
!!      b : The boundary energy [ergs]
!! 
!!***
subroutine RadTrans_mgdSetBound(g, b)
  implicit none

  integer, intent(in) :: g
  real,    intent(in) :: b
  ! Stub implementation
  return

end subroutine RadTrans_mgdSetBound
