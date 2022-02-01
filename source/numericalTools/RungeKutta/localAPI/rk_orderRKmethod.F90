!!****if* source/numericalTools/RungeKutta/localAPI/rk_orderRKmethod
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
!!  rk_orderRKmethod
!!
!! SYNOPSIS
!!
!!  rk_orderRKmethod (character (len=*), intent (in) :: method)
!!
!! DESCRIPTION
!!
!!  This integer function returns the order of the specified Runge Kutta method.
!!
!! ARGUMENTS
!!
!!  method : the Runge Kutta method
!!
!! NOTES
!!
!!***

integer function rk_orderRKmethod (method)

  implicit none

  character (len=*), intent (in) :: method

  rk_orderRKmethod = 0

  return
end function rk_orderRKmethod
