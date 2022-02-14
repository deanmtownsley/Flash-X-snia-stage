!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleDeallocateRadialArrays
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
!!  gr_mpoleDeallocateRadialArrays
!!
!! SYNOPSIS
!!
!!  gr_mpoleDeallocateRadialArrays ()
!!
!! DESCRIPTION
!!
!!  Deallocates the radial multi-bin arrays and the arrays defining the radial inner zone.
!!
!!***

subroutine gr_mpoleDeallocateRadialArrays ()

  use gr_mpoleData,  ONLY : gr_mpoleInnerZoneDrRadii, &
                            gr_mpoleInnerZoneQupper,  &
                            gr_mpoleInnerZoneQlower,  &
                            gr_mpoleMomentR,          &
                            gr_mpoleMomentI,          &
                            gr_mpoleQDampingR,        &
                            gr_mpoleQDampingI,        &
                            gr_mpoleQRadii,           &
                            gr_mpoleQused,            &
                            gr_mpoleScratch
 
  implicit none
!
!
!       ...Deallocate the arrays.
!
!
  if (allocated (gr_mpoleInnerZoneDrRadii)) deallocate (gr_mpoleInnerZoneDrRadii)
  if (allocated (gr_mpoleInnerZoneQupper) ) deallocate (gr_mpoleInnerZoneQupper)
  if (allocated (gr_mpoleInnerZoneQlower) ) deallocate (gr_mpoleInnerZoneQlower)
  if (allocated (gr_mpoleMomentR)         ) deallocate (gr_mpoleMomentR)
  if (allocated (gr_mpoleMomentI)         ) deallocate (gr_mpoleMomentI)
  if (allocated (gr_mpoleQDampingR)       ) deallocate (gr_mpoleQDampingR)
  if (allocated (gr_mpoleQDampingI)       ) deallocate (gr_mpoleQDampingI)
  if (allocated (gr_mpoleQRadii)          ) deallocate (gr_mpoleQRadii)
  if (allocated (gr_mpoleQused)           ) deallocate (gr_mpoleQused)
  if (allocated (gr_mpoleScratch)         ) deallocate (gr_mpoleScratch)
!
!
!       Done.
!
!
  return
end subroutine gr_mpoleDeallocateRadialArrays
