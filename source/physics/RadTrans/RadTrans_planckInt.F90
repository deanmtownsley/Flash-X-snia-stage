!!****f* source/physics/RadTrans/RadTrans_planckInt
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
!!  RadTrans_planckInt
!!
!!  SYNOPSIS
!!
!!  call RadTrans_planckInt(real(IN)  :: x,
!!                           real(OUT) :: p)
!!
!!  DESCRIPTION 
!!      This routine must evaluate the Planck integral
!!
!!      The Planck integral is:
!!      
!!             x/|        y^3
!!      P(x) =   |  dy ----------
!!             0 |/    exp(y) - 1
!!
!!
!! ARGUMENTS
!!
!!   x : The argument for the planck integral
!!   p : The value of the integral
!!
!!***
subroutine RadTrans_planckInt(x, p)

  implicit none

  ! Arguments:
  real, intent(in)  :: x
  real, intent(out) :: p

  ! Stub implementation
  p = 0.0
  return
end subroutine RadTrans_planckInt
