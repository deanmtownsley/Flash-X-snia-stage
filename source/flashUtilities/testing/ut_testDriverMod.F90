!!****if* source/flashUtilities/testing/ut_testDriverMod
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
!!  ut_testDriverMod
!!
!! SYNOPSIS
!!
!!  use ut_testDriverMod
!!
!! DESCRIPTION
!!
!!  A module that encapsulates common routines and variables that can be used in
!!  unittests.  Typically,  this module is used in the unittests version of 
!!  Driver_evolveAll as follows
!!
!!      <declare variables and setup for test>
!!      ...
!!      call start_test_run
!!      ...
!!      call assertFalse(didAbort, "Action aborted")
!!      call assertAlmostEqual(x, 1.1, 1.0d-15, "Incorrect x value")
!!      ...
!!      call finish_test_run
!!      ...
!!      <clean-up>
!!
!! NOTES
!!  
!!  The finish_test_run prints results to standard output.
!!
!!***

#include "constants.h"

module ut_testDriverMod
    implicit none
    private

    integer, save :: my_n_tests = 0
    integer, save :: my_n_failed = 0
    logical, save :: is_testing = .FALSE.

    interface assertEqual
        procedure :: assertEqualInt
        procedure :: assertEqualReal
    end interface assertEqual

    interface assertAlmostEqual
        procedure :: assertAlmostEqualAuto
        procedure :: assertAlmostEqual
    end interface assertAlmostEqual

    interface assertSetEqual
        procedure :: assertSetEqual2dIntArray
    end interface assertSetEqual

    public :: start_test_run
    public :: finish_test_run

    public :: assertTrue
    public :: assertFalse
    public :: assertEqual
    public :: assertSetEqual
    public :: assertAlmostEqual

contains

    subroutine start_test_run()
        use Driver_Interface, ONLY : Driver_abort

        if (is_testing) then
            call Driver_abort("[start_test_run] Already testing")
        end if

        is_testing = .TRUE.
    end subroutine start_test_run

    function finish_test_run() result(did_succeed)
        use Driver_data,      ONLY : dr_globalMe
        use Driver_Interface, ONLY : Driver_abort

        logical :: did_succeed

        if (.NOT. is_testing) then
            call Driver_abort("[finish_test_run] Not testing yet")
        end if

        is_testing = .FALSE.

        ! DEV: TODO reduction to collect number of tests/fails/max walltime?
        did_succeed = (my_n_failed == 0)
        if (dr_globalMe == MASTER_PE) then
            ! Print result to standard out
            write(*,*)
            if (did_succeed) then
                write(*,*) "SUCCESS - ", &
                           (my_n_tests - my_n_failed), "/", my_n_tests, ' passed' 
            else 
                write(*,*) "FAILURE - ", &
                           (my_n_tests - my_n_failed), "/", my_n_tests, ' passed'
            end if
        end if
    end function finish_test_run

    subroutine assertTrue(a, msg)
        logical,      intent(IN) :: a
        character(*), intent(IN) :: msg

        character(256) :: buffer
        
        if (.NOT. a) then
            write(buffer,'(A)') msg
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertTrue

    subroutine assertFalse(a, msg)
        logical,      intent(IN) :: a
        character(*), intent(IN) :: msg

        character(256) :: buffer
        
        if (a) then
            write(buffer,'(A)') msg
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertFalse

    subroutine assertEqualInt(a, b, msg)
        integer,      intent(IN) :: a
        integer,      intent(IN) :: b
        character(*), intent(IN) :: msg

        character(256) :: buffer

        if (a /= b) then
            write(buffer,'(2A,I0,A,I0)') msg, " ", a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertEqualInt

    subroutine assertEqualReal(a, b, msg)
        real,         intent(IN) :: a
        real,         intent(IN) :: b
        character(*), intent(IN) :: msg

        character(256) :: buffer

        if (a /= b) then
600         format(A,1P,G24.16,A,G24.16)
            write(buffer,600) msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertEqualReal

    subroutine assertAlmostEqualAuto(a, b, msg)
        real,         intent(IN) :: a
        real,         intent(IN) :: b
        character(*), intent(IN) :: msg

        real :: prec
        character(256) :: buffer

        prec = 2.0 * spacing(min(abs(a),abs(b)))
        if (ABS(b - a) > prec) then
700         format(A,' (with auto tol)',6P,G24.16,A,0P,G24.16)
            write(buffer,700) msg, a, " != ", b
200         format((A))
            write(*,200) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertAlmostEqualAuto

    subroutine assertAlmostEqual(a, b, prec, msg)
        real,         intent(IN) :: a
        real,         intent(IN) :: b
        real,         intent(IN) :: prec
        character(*), intent(IN) :: msg

        character(256) :: buffer

        if (ABS(b - a) > prec) then
            write(buffer,'(A,F15.8,A,F15.8)') msg, a, " != ", b
            write(*,*) TRIM(ADJUSTL(buffer))
            my_n_failed = my_n_failed + 1
        end if
        my_n_tests = my_n_tests + 1
    end subroutine assertAlmostEqual

    subroutine assertSetEqual2dIntArray(A, B, msg)
        integer,      intent(IN) :: A(:, :)
        integer,      intent(IN) :: B(:, :)
        character(*), intent(IN) :: msg

        logical        :: in_set
        logical        :: failed
        integer        :: j, k

        my_n_tests = my_n_tests + 1

        ! Confirm A subset of B
        failed = .FALSE.
        do j = 1, SIZE(A, 1)
            in_set = .FALSE.
            do k = 1, SIZE(B, 1)
                if (ALL(A(j, :) == B(k, :))) then
                    in_set = .TRUE.
                    exit
                end if
            end do

            if (.NOT. in_set) then
                write(*,*) msg, " - ", A(j, :), " of A not in B"
                failed = .TRUE.
            end if
        end do

        ! Confirm B subset of A
        do j = 1, SIZE(B, 1)
            in_set = .FALSE.
            do k = 1, SIZE(A, 1)
                if (ALL(B(j, :) == A(k, :))) then
                    in_set = .TRUE.
                    exit
                end if
            end do

            if (.NOT. in_set) then
                write(*,*) msg, " - ", B(j, :), " of B not in A"
                failed = .TRUE.
            end if
        end do

        if (failed) then
            my_n_failed = my_n_failed + 1
        end if
    end subroutine assertSetEqual2dIntArray

end module ut_testDriverMod

