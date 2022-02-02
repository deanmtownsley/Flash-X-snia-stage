!!****if* source/Grid/localAPI/gr_mpoleRadialSampling
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
!!  gr_mpoleRadialSampling
!!
!! 
!! SYNOPSIS
!!
!!  gr_mpoleRadialSampling()
!!
!!
!! DESCRIPTION
!!
!!  Determines the radial sampling for accumulating the Moments.
!!  It defines the radial bins into which each Moment falls. If the
!!  radial sampling is chosen unwisely, lots of memory for holding
!!  empty Moment bins might be wasted. This routine calls the appropriate
!!  subroutines according to the geometry specified.
!!
!!
!!***

subroutine gr_mpoleRadialSampling()

  implicit none

  return
end subroutine gr_mpoleRadialSampling
