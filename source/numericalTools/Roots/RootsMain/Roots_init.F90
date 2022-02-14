!!****if* source/numericalTools/Roots/RootsMain/Roots_init
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
!!  Roots_init
!!
!! SYNOPSIS
!! 
!!  call Roots_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the Roots unit.
!!
!! ARGUMENTS
!!
!!***

subroutine Roots_init ()

  use Roots_data

  implicit none
!
!
!     ...Set all machine/precision dependent parameters that control
!        quality/accuracy of the roots.
!
!
  rt_macheps = epsilon (1.0)
  rt_LPN     = huge    (1.0)
  rt_sqrtLPN = sqrt (rt_LPN)
!
!
!    ...Ready!
!
!
  return
end subroutine Roots_init
