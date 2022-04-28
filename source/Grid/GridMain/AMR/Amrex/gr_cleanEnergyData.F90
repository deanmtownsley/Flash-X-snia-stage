!!****if* source/Grid/GridMain/AMR/Amrex/gr_cleanEnergyData
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
!!  gr_cleanEnergyData
!!
!! SYNOPSIS
!!
!!  call gr_cleanEnergyData(real(in) :: smallE,
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
!!  This subroutine checks only whether energy variables (the ENER_VAR and
!! EINT_VAR components IN UNK) are greater than or equal to smallE.
!!
!!
!! ARGUMENTS
!!
!!  smallE : energy floor, probably from runtime parameter smallE
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
!!  gr_cleanDensityData
!!  gr_postinterpolationWork
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine gr_cleanEnergyData(smallE, &
                              lo, hi, &
                              d, dlo, dhi, nd)
  use Driver_interface, ONLY : Driver_abort
  use Grid_data,        ONLY : gr_sanitizeDataMode

  implicit none
  
  real,    intent(in)    :: smallE
  integer, intent(in)    :: lo(MDIM), hi(MDIM)
  integer, intent(in)    :: dlo(MDIM), dhi(MDIM)
  integer, intent(in)    :: nd
  real,    intent(inout) :: d(dlo(IAXIS):dhi(IAXIS), &
                              dlo(JAXIS):dhi(JAXIS), &
                              dlo(KAXIS):dhi(KAXIS), &
                              nd)

  integer :: i, j, k

  if (gr_sanitizeDataMode == 0) RETURN

900 format("         Value at",I3,",",I3,",",I3," set to smallE=",1PG8.2)
#ifdef ENER_VAR
  ! DEV: TODO Determine how to implement all modes and all levels of 
  !           verbosity
  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        if (d(i,j,k,ENER_VAR) < smallE*0.999999999) then
           ! The factor 0.999999999 is used to avoid triggering aborts
           ! when some physics code, in particular in the Eos unit,
           ! has put a floor of smallE under a very low internal
           ! energy and this then becomes even slightly smaller as a
           ! result of rounding. If it is intended that situations
           ! like that produce messages or aborts, it should be the
           ! responsibility of the physics unit that applied the
           ! smallE floor to take care that that happens.
          if      (gr_sanitizeDataMode > 0) then
            write(*,*) "WARNING: [gr_cleanEnergyData]"
            write(*,*) "         Total energy data less than smallE:", d(i,j,k,ENER_VAR)
          end if
          if      (gr_sanitizeDataMode == 3) then
            write(*,900)         i, j, k, smallE
            d(i,j,k,ENER_VAR) = max(d(i,j,k,ENER_VAR), smallE)
          else if (gr_sanitizeDataMode == 4) then
            call Driver_abort("[gr_cleanEnergyData] ENER data less than smallE")
          end if
        end if
      end do
    end do
  end do
#endif
#ifdef EINT_VAR
  do     k = lo(KAXIS), hi(KAXIS) 
    do   j = lo(JAXIS), hi(JAXIS) 
      do i = lo(IAXIS), hi(IAXIS)
        if (d(i,j,k,EINT_VAR) < smallE*0.999999999) then
           ! For the factor 0.999999999 see above.
          if      (gr_sanitizeDataMode > 0) then
            write(*,*) "WARNING: [gr_cleanEnergyData]"
            write(*,*) "         Internal energy data less than smallE:", d(i,j,k,EINT_VAR)
          end if
          if      (gr_sanitizeDataMode == 3) then
            write(*,900)         i, j, k, smallE
            d(i,j,k,EINT_VAR) = max(d(i,j,k,EINT_VAR), smallE)
          else if (gr_sanitizeDataMode == 4) then
            call Driver_abort("[gr_cleanEnergyData] EINT data less than smallE")
          end if
        end if
      end do
    end do
  end do
#endif

end subroutine gr_cleanEnergyData

