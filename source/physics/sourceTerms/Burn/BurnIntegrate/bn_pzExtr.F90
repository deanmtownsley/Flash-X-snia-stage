!!****if* source/physics/sourceTerms/Burn/BurnIntegrate/bn_pzExtr.F90
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
!!  bn_pzExtr
!!
!! SYNOPSIS
!!
!!  subroutine bn_pzExtr(integer(IN) :: iest,
!!                       real(IN)    :: xest,
!!                       real(IN)    :: yest(:),
!!                       real(INOUT)   :: yz(:),
!!                       real(OUT)   :: dy(:),
!!                       integer(IN) :: nv)
!!
!! DESCRIPTION
!!
!!  use polynomial extrapolation to evaluate nv functions at x=0 by fitting
!!  a polynomial to a sequence of estimates with progressively smaller values
!!  x=xest, and corresponding function vectors yest(1:nv). the call is number
!!  iest in the sequence of calls. extrapolated function values are output as
!!  yz(1:nv), and their estimated error is output as dy(1:nv)
!!
!! ARGUMENTS
!!
!!  iest - integer(IN)      number in the sequence of calls
!!  xest - real(IN)         fitting location
!!  yest - real(IN)(1:nv)   function vectors at location xest
!!  yz -   real(INOUT)(1:nv)  extrapolated function values
!!  dy -   real(OUT)(1:nv)  estimated error of extrapolated function values
!!  nv -   integer(IN)      number of functions to evaluate
!!
!!***
subroutine bn_pzExtr(state, iest, xest, yest, yz, dy, nv)

  use Driver_interface, ONLY: Driver_abort
  use bnIntegrate_interface, ONLY: pzExtr_state_t

  implicit none

  !! Declare arguments
  type(pzExtr_state_t), intent(INOUT) :: state
  integer, intent(IN)  :: nv, iest
  real, intent(IN)     :: xest
  real, intent(IN), dimension(nv) :: yest
  real, intent(OUT), dimension(nv) :: dy, yz

  !! Local variables
  integer, parameter :: nmax=50, imax=13
  integer :: j, k1
  real :: delta, f1, f2, q

  !! Check for array bounds
  if (nv > nmax) then
     call Driver_abort('ERROR in bn_pzExtr: nv exceeds nmax')
  end if
  if (iest > imax) then
     call Driver_abort('ERROR in bn_pzExtr: iest exceeds imax')
  end if

  !! Save current independent variable
  state%x(iest) = xest
  do j = 1, nv
     dy(j) = yest(j)
     yz(j) = yest(j)
  end do

  !! Store first estimate in the first column
  if (iest == 1) then
     do j = 1, nv
        state%qcol(j,1) = yest(j)
     end do

  else
     do j = 1, nv
        state%d(j) = yest(j)
     end do
     do k1 = 1, iest - 1
        delta = 1.0e0 / (state%x(iest - k1) - xest)
        f1    = xest * delta
        f2    = state%x(iest - k1) * delta

        !! Propagate tableau one diagonal more
        do j = 1, nv
           q                 = state%qcol(j, k1)
           state%qcol(j, k1) = dy(j)
           delta             = state%d(j) - q
           dy(j)             = f1 * delta
           state%d(j)        = f2 * delta
           yz(j)             = yz(j) + dy(j)
        end do
     end do
     do j = 1, nv
        state%qcol(j, iest) = dy(j)
     end do
  end if

  return
end subroutine bn_pzExtr
