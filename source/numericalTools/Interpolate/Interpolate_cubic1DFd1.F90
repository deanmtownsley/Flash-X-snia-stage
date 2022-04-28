!!****f* source/numericalTools/Interpolate/Interpolate_cubic1DFd1
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
!!  Interpolate_cubic1DFd1
!!
!! SYNOPSIS
!!
!!  Interpolate_cubic1DFd1 (real, intent (in) :: a (1:4),
!!                          real, intent (in) :: x)
!!
!! DESCRIPTION
!!
!!  Calculates the function value and the rescaled 1st derivative value for a single [x]
!!  rescaled [0,1] coordinate and the 4 monocubic expansion coefficients. The monocubic
!!  expansion reads, for one line, in terms of the rescaled [0,1] x coordinate:
!!
!!                                 3          i
!!                        F (x) = sum  a (i) x
!!                                i=0
!!
!!  The location index of the a (i) inside the 4-dimensional vector is:
!!
!!                location index of a (i)  =  1 + i
!!
!!  The rescaled derivatives are given by the general formula:
!!
!!                 r    r    3                  i-r
!!                d / dx  = sum (i)  * a (i) * x
!!                          i=r    r
!!
!!  where the Pochhammer symbols are defined as:
!!
!!                  (i)  = i * (i-1) * (i-2) * ... * (i-r+1)
!!                     r
!!
!!  From the derivative formula we see, that the highest non-zero derivative is of
!!  3-rd order.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th monocubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!
!! NOTES
!!
!!  1) The function is defined as a real array of size 2:
!!
!!           Interpolate_cubic1DFd1 (1) = the function value
!!           Interpolate_cubic1DFd1 (2) = the rescaled d/dx value
!!
!!  2) The code checks, if the supplied coordinate [x] is rescaled.
!!
!!***

function Interpolate_cubic1DFd1 (a,x)

  implicit none

  real, intent (in) :: a (1:4)
  real, intent (in) :: x

  real :: Interpolate_cubic1DFd1 (1:2)  ! declares the function as an array

  Interpolate_cubic1DFd1 (1:2) = 0.0

  return
end function Interpolate_cubic1DFd1
