!!****ih* source/numericalTools/RungeKutta/localAPI/rk_interface
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
!!  rk_interface
!!
!! SYNOPSIS
!!
!!  use rk_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Runge Kutta unit that defines its private interfaces.
!!
!!***

Module rk_interface

  interface
    integer function rk_orderRKmethod (method)
      character (len=*), intent (in) :: method
    end function rk_orderRKmethod
  end interface

  interface
    subroutine rk_setButcherTableauRepository () 
    end subroutine rk_setButcherTableauRepository
  end interface

  interface
    subroutine rk_stepEA (f, n, s,                      &
                          x, y,                         &
                          eMin, ePower, eFrac, eBase,   &
                          htry,                         &
                                          hused, hnext, &
                                          yout, eout    )
      interface
        function f (x,y)
          real, intent (in) :: x
          real, intent (in) :: y (:)
          real              :: f (1:size (y))
        end function f
      end interface

      integer, intent (in)  :: n
      integer, intent (in)  :: s
      real,    intent (in)  :: x
      real,    intent (in)  :: y     (:)
      real,    intent (in)  :: eMin
      real,    intent (in)  :: ePower
      real,    intent (in)  :: eFrac
      real,    intent (in)  :: eBase (:)
      real,    intent (in)  :: htry
      real,    intent (out) :: hused
      real,    intent (out) :: hnext
      real,    intent (out) :: yout  (:)
      real,    intent (out) :: eout  (:)
    end subroutine rk_stepEA
  end interface

  interface
    subroutine rk_stepEAC (f, n, nc, s,                  &
                           x, y, ymin, ymax,             &
                           eMin, ePower, eFrac, eBase,   &
                           htry,                         &
                                           hused, hnext, &
                                           yout, eout    )
      interface
        function f (x,y)
          real, intent (in) :: x
          real, intent (in) :: y (:)
          real              :: f (1:size (y))
        end function f
      end interface

      interface
        function ymin (nc,y)
          integer, intent (in) :: nc
          real,    intent (in) :: y (:)
          real                 :: ymin (1:nc)
        end function ymin
      end interface

      interface
        function ymax (nc,y)
          integer, intent (in) :: nc
          real,    intent (in) :: y (:)
          real                 :: ymax (1:nc)
        end function ymax
      end interface

      integer, intent (in)  :: n
      integer, intent (in)  :: nc
      integer, intent (in)  :: s
      real,    intent (in)  :: x
      real,    intent (in)  :: y     (:)
      real,    intent (in)  :: eMin
      real,    intent (in)  :: ePower
      real,    intent (in)  :: eFrac
      real,    intent (in)  :: eBase (:)
      real,    intent (in)  :: htry
      real,    intent (out) :: hused
      real,    intent (out) :: hnext
      real,    intent (out) :: yout  (:)
      real,    intent (out) :: eout  (:)
    end subroutine rk_stepEAC
  end interface

end Module rk_interface
