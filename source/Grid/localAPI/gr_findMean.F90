!!****if* source/Grid/localAPI/gr_findMean
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
!!  gr_findMean
!! 
!! SYNOPSIS
!!  call gr_findMean(integer(in)  :: iSrc,
!!                   integer(in)  :: iType,
!!                   logical(in)  :: bGuardcell,
!!                      real(out) :: mean)
!!
!! DESCRIPTION
!!  Calculates the mean of a function
!!
!!
!! ARGUMENTS
!!  iSrc -- the index (e.g. DENS_VAR) into the unk routine to calculate over
!!  iType -- the type of mean.  Valid values will be  (feel free to implement more)
!!     1 = L1 norm = arithmetic average
!!     2 = arithmetic average
!!  bGuardcell -- logical indicating whether guard cells should be included in the calculation
!!  mean -- the requested mean
!!
!!
!!***

subroutine gr_findMean(iSrc, iType, bGuardcell, mean)
  
  implicit none

  integer, intent(in) :: iSrc, iType
  logical, intent(in) :: bGuardcell
  real, intent(out) :: mean


end subroutine gr_findMean

