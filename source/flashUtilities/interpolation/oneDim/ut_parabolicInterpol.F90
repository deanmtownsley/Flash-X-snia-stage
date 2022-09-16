!!****if* source/flashUtilities/interpolation/oneDim/ut_parabolicInterpol
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
!!  ut_parabolicInterpol
!!
!! DESCRIPTION
!!
!!   Given a vector of coordinates, x, the size, n, and associated function
!!   values, var, take the zone edges and center (y_l, y_c, y_r), and return
!!   the parabolic interpolation to this.
!!
!!***
subroutine ut_parabolicInterpol(x,var,n,y_l,y_c,y_r,var_interp)
    ! Given a vector of coordinates, x, the size, n, and associated function
    ! values, var, take the zone edges and center (y_l, y_c, y_r), and return
    ! the parabolic interpolation to this.
    !
    ! x(n)        coordinate values
    ! var(n)      function values at x(n)
    ! n           size of x and var
    !
    ! y_l         coordinate of left edge of the zone
    ! y_c         coordinate of center of the zone
    ! y_r         coordinate of right edge of the zone
    !
    ! var_interp  zone average value of the function in that zone, using
    !             a parabolic structure
    !
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x(n), var(n)
    real,intent(in) :: y_l, y_c, y_r
    real,intent(out) :: var_interp

    real :: var_l, var_c, var_r, err_int
    integer, parameter :: op = 2
    real, parameter :: sixth = 1.d0/6.d0
    integer :: kat

    ! get the function value at the left edge of the zone
    kat = 0

    if (y_l < x(1)) then
        ! the x array is monotonic -- if we are below the minimum in x, then
        ! just set the index to the first value
        kat = 1
    else
        call ut_hunt(x,n,y_l,kat)
    endif

    kat = max(1, min(kat - op/2 + 1, n - op + 1))
    call ut_polint(x(kat),var(kat),op,y_l,var_l,err_int)

    ! get the function value at the center of the zone
    call ut_hunt(x,n,y_c,kat)
    kat = max(1, min(kat - op/2 + 1, n - op + 1))
    call ut_polint(x(kat),var(kat),op,y_c,var_c,err_int)

    ! get the function value at the right edge of the zone
    call ut_hunt(x,n,y_r,kat)
    kat = max(1, min(kat - op/2 + 1, n - op + 1))
    call ut_polint(x(kat),var(kat),op,y_r,var_r,err_int)
    
    ! construct the zone averaged value
    var_interp = (var_l + 4.e0*var_c + var_r)*sixth

    return
end subroutine ut_parabolicInterpol
