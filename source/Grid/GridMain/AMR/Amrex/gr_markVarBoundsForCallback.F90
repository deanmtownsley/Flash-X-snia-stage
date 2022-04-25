!!****if* source/Grid/GridMain/AMR/Amrex/gr_markVarBoundsForCallback
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
!!  gr_markVarBoundsForCallback
!!
!! SYNOPSIS
!!  gr_markVarBoundsForCallback(  integer(in) :: var
!!                                real(in)    :: var_bnd_min,
!!                                real(in)    :: var_bnd_max,
!!                                integer(in) :: lev,
!!                                c_ptr(in)   :: tags,
!!                                c_char(in)  :: tagval)
!!
!!
!! PURPOSE
!!  Refine all blocks for which a given variable (Var) exceeds or falls
!!  between bounds defined by var_bnd_min and var_bnd_max.
!!
!! ARGUMENTS
!!  Var -    the variable of interest
!!
!!  var_bnd_min  -  lower limit for bound
!!
!!  var_bnd_max  -  upper limit for bound
!!
!!  lev - the 0-based level index
!!  tags - C-pointer to an AMReX tagbox array.  The elements of this are tag
!!         boxes.  The cells of these tagboxes are set to communicate a need
!!         to refine the associated block.
!!  tagval - for full, rich AMReX tagging, this values should be assigned to
!!           each cell that has insufficient resolution.
!!
!! NOTES
!!
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!
!!
!!***

subroutine gr_markVarBoundsForCallback(Var, var_bnd_min, var_bnd_max, lev, tags, tagval)

!-------------------------------------------------------------------------------
   use iso_c_binding
   use amrex_fort_module, ONLY: wp => amrex_real
   use amrex_box_module, ONLY: amrex_box
   use amrex_tagbox_module, ONLY: amrex_tagboxarray
   use amrex_multifab_module, ONLY: amrex_mfiter, &
                                    amrex_mfiter_build, &
                                    amrex_mfiter_destroy

   use Driver_interface, ONLY: Driver_abort
   use Grid_data, ONLY: gr_geometry
   use Grid_interface, ONLY: Grid_getBlkCenterCoords
   use gr_physicalMultifabs, ONLY: unk
   use Grid_tile, ONLY: Grid_tile_t

#include "constants.h"
#include "Simulation.h"

   implicit none

! Arguments
   integer, intent(IN) :: var
   real, intent(IN) :: var_bnd_min, var_bnd_max
   integer, intent(IN) :: lev
   type(c_ptr), intent(IN) :: tags
   character(c_char), intent(IN) :: tagval

! Local data

   type(amrex_tagboxarray) :: tag
   type(amrex_mfiter)      :: mfi
   type(amrex_box)         :: bx
   type(Grid_tile_t)       :: blockDesc

   character(c_char), contiguous, pointer :: tagData(:, :, :, :)
   real(wp), contiguous, pointer :: solnData(:, :, :, :)

   real, dimension(MDIM) :: blockCenter, blockSize
   real                  :: xl, xr, yl, yr, zl, zr
   integer               :: b
   logical               :: x_in_rect, y_in_rect, z_in_rect

   integer :: i, j, k
   integer, dimension(MDIM) :: lo, hi

   tag = tags

   call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)
   do while (mfi%next())
      bx = mfi%fabbox()

      blockDesc%level = lev + 1
      blockDesc%grid_index = mfi%grid_index()
      blockDesc%limits(LOW, :) = 1
      blockDesc%limits(HIGH, :) = 1
      blockDesc%limits(LOW, 1:NDIM) = bx%lo(1:NDIM) + 1 + NGUARD
      blockDesc%limits(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1 - NGUARD
      blockDesc%blkLimitsGC(LOW, :) = 1
      blockDesc%blkLimitsGC(HIGH, :) = 1
      blockDesc%blkLimitsGC(LOW, 1:NDIM) = bx%lo(1:NDIM) + 1
      blockDesc%blkLimitsGC(HIGH, 1:NDIM) = bx%hi(1:NDIM) + 1
      blockDesc%grownLimits(:, :) = blockDesc%blkLimitsGC(:, :)

      call Grid_getBlkCenterCoords(blockDesc, blockCenter)
      call blockDesc%physicalSize(blockSize)
      blockSize(:) = 0.5*blockSize(:)

      lo = blockDesc%limits(LOW, :)
      hi = blockDesc%limits(HIGH, :)

      tagData => tag%dataptr(mfi)
      solnData => unk(lev)%dataPtr(mfi)

      associate (lo => blockDesc%limits(LOW, :), &
                 hi => blockDesc%limits(HIGH, :), &
                 lo_tag => lbound(tagData), &
                 hi_tag => ubound(tagData))

         if (maxval(solnData(lo(IAXIS):hi(IAXIS), &
                             lo(JAXIS):hi(JAXIS), &
                             lo(KAXIS):hi(KAXIS), var)) > var_bnd_min .AND. &
             minval(solnData(lo(IAXIS):hi(IAXIS), &
                             lo(JAXIS):hi(JAXIS), &
                             lo(KAXIS):hi(KAXIS), var)) < var_bnd_max) then

            i = INT(0.5d0*DBLE(lo_tag(IAXIS) + hi_tag(IAXIS)))
            j = INT(0.5d0*DBLE(lo_tag(JAXIS) + hi_tag(JAXIS)))
            k = INT(0.5d0*DBLE(lo_tag(KAXIS) + hi_tag(KAXIS)))

            ! Fourth index is 1:1
            tagData(i, j, k, 1) = tagval
         end if

      end associate

      nullify (tagData)
      nullify (solnData)
   end do
   call amrex_mfiter_destroy(mfi)

   return
end subroutine gr_markVarBoundsForCallback
