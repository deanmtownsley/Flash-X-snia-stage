!!****if* source/Grid/GridMain/AMR/Amrex/gr_cleanDensityData
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
!!  gr_cleanDensityData
!!
!! SYNOPSIS
!!
!!  call gr_cleanDensityData( real  (in)       :: smallrho,
!!                           integer(in)       :: lo(MDIM),
!!                           integer(in)       :: hi(MDIM),
!!                           amrex_real(INOUT) :: d(dlo(IAXIS):dhi(IAXIS),
!!                                                  dlo(JAXIS):dhi(JAXIS),
!!                                                  dlo(KAXIS):dhi(KAXIS),
!!                                                  nd),
!!                           integer(in)       :: dlo(MDIM),
!!                           integer(in)       :: dhi(MDIM),
!!                           integer(in)       :: nd,
!!                           integer(in)       :: scomp,
!!                           integer(in)       :: ncomp)
!!
!! DESCRIPTION
!!
!!  Check some data for basic sanity, and optionally either clean it up
!!  a bit or abort if something is found wrong.
!!
!!  This subroutine checks only whether density (the DENS_VAR component
!!  IN UNK) is greater than or equal to smallRho.
!!
!! ARGUMENTS
!!
!!  smallrho : density floor, probably from runtime parameter smlrho
!!  lo - the index of the lower-left corner of the region of the given
!!       box that (potentially) requires cleaning.
!!  hi - the index of the upper-right corner of the region of the given
!!       box that (potentially) requires cleaning.
!!  d - the data array for the box that (potentially) requires cleaning.
!!  dlo - the lower-left index of the given box
!!  dhi - the upper-right index of the given box
!!  nd -  the number of physical quantities in box
!!  scomp - the index of the first physical quantity in the box
!!  ncomp - the number of physical quantities in the box
!!
!! NOTES
!!
!!  This subroutine does a part of what gr_sanitizeDataAfterInterp does
!!  In the paramesh Grid implementation.
!!
!! SEE ALSO
!!  gr_cleanEnergyData
!!  gr_postinterpolationWork
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine gr_cleanDensityData(smallRho, &
                               lo, hi, &
                               d, dlo, dhi, nd)
  use Driver_interface, ONLY : Driver_abort
  use Grid_data,        ONLY : gr_sanitizeDataMode

  implicit none

  real,    intent(in)    :: smallRho
  integer, intent(in)    :: lo(MDIM), hi(MDIM)
  integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
  integer, intent(in)    :: nd
  real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                              dlo(JAXIS):dhi(JAXIS), &
                              dlo(KAXIS):dhi(KAXIS), &
                              nd)

  integer :: i, j, k

#ifdef DENS_VAR

  if (gr_sanitizeDataMode == 0) RETURN

900 format("         Value at",I3,",",I3,",",I3," set to smlrho=",1PG8.2)

  ! DEV: TODO Determine how to implement all modes and all levels of 
  !           verbosity
  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        if (d(i,j,k,DENS_VAR) < smallRho) then
          if      (gr_sanitizeDataMode > 0) then
            write(*,*) "WARNING: [gr_cleanDensityData]"
            write(*,*) "         Density data less than smlrho:", d(i,j,k,DENS_VAR)
          end if
          if      (gr_sanitizeDataMode == 3) then
            write(*,900)         i, j, k, smallRho
            d(i,j,k,DENS_VAR) = max(d(i,j,k,DENS_VAR), smallRho)
          else if (gr_sanitizeDataMode == 4) then
            call Driver_abort("[gr_cleanDensityData] Density data less than smlrho")
          end if
        end if
      end do
    end do
  end do
#endif

end subroutine gr_cleanDensityData

