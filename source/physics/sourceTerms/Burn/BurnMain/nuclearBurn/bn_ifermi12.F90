!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_ifermi12
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
!!   bn_ifermi12
!!
!! SYNOPSIS
!!
!!   real function bn_ifermi12(real(IN) :: f)
!!
!! ARGUMENTS
!!
!!  f -- umm, sorry don't know what f is in ifermi12
!!  
!! DESCRIPTION
!!
!!  routine ifermi12 does an inverse fermi integral for bn_sneutx
!!  this routine applies a rational function expansion to get the inverse
!!  fermi-dirac integral of order 1/2 when it is equal to f.
!!  maximum error is 4.19e-9.   reference: antia apjs 84,101 1993
!!
!!***

real function bn_ifermi12(f)

   implicit none

   !! Declare arguments
   real, intent(IN) :: f

   !! Declare local variables
   integer :: i
   real :: rn, den, ff

   !! Constants (converted from saved variables)
   integer, parameter :: m1 = 4, k1 = 3, m2 = 6, k2 = 5
   real, parameter :: an = 0.5e0

   real, dimension(5), parameter :: a1 = (/ &
      1.999266880833e4, 5.702479099336e3, &
      6.610132843877e2, 3.818838129486e1, 1.0e0 /)

   real, dimension(4), parameter :: b1 = (/ &
      1.771804140488e4, -2.014785161019e3, &
      9.130355392717e1, -1.670718177489e0 /)

   real, dimension(7), parameter :: a2 = (/ &
      -1.277060388085e-2, 7.187946804945e-2, &
      -4.262314235106e-1, 4.997559426872e-1, &
      -1.285579118012e0, -3.930805454272e-1, 1.0e0 /)

   real, dimension(6), parameter :: b2 = (/ &
      -9.745794806288e-3, 5.485432756838e-2, &
      -3.299466243260e-1, 4.077841975923e-1, &
      -1.145531476975e0, -6.067091689181e-2 /)

   !! Function logic
   if (f < 4.0e0) then
      rn = f + a1(m1)
      do i = m1 - 1, 1, -1
         rn = rn * f + a1(i)
      end do
      den = b1(k1 + 1)
      do i = k1, 1, -1
         den = den * f + b1(i)
      end do
      bn_ifermi12 = log(f * rn / den)
   else
      ff = 1.0e0 / f**(1.0e0 / (1.0e0 + an))
      rn = ff + a2(m2)
      do i = m2 - 1, 1, -1
         rn = rn * ff + a2(i)
      end do
      den = b2(k2 + 1)
      do i = k2, 1, -1
         den = den * ff + b2(i)
      end do
      bn_ifermi12 = rn / (den * ff)
   end if

   return
end function bn_ifermi12


