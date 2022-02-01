!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleAllocateRadialArrays
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
!!  gr_mpoleAllocateRadialArrays
!!
!! SYNOPSIS
!!
!!  gr_mpoleAllocateRadialArrays ()
!!
!! DESCRIPTION
!!
!!  Allocates the radial multi-bin arrays. Before this routine can be called,
!!  the maximum number of radial bins for the current iteration must have
!!  been established.
!!
!!***

subroutine gr_mpoleAllocateRadialArrays ()

  use gr_mpoleData,      ONLY : gr_mpoleMaxQ,      &
                                gr_mpoleMaxLM,     &
                                gr_mpoleMomentR,   &
                                gr_mpoleMomentI,   &
                                gr_mpoleQDampingR, &
                                gr_mpoleQDampingI, &
                                gr_mpoleQRadii,    &
                                gr_mpoleQused,     &
                                gr_mpoleScratch

  use Driver_interface,  ONLY : Driver_abort
 
  implicit none

  integer :: status
!
!
!       ...Allocate the radial multi-bin moment arrays and the corresponding
!          scratch array. Note the extra outer radial bin for the irregular
!          moment array and the extra inner radial bin for the regular moment
!          array!
!
!
  allocate (gr_mpoleMomentR (1:gr_mpoleMaxLM,0:gr_mpoleMaxQ  ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleMomentR allocate failed')
  end if

  allocate (gr_mpoleMomentI (1:gr_mpoleMaxLM,1:gr_mpoleMaxQ+1), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleMomentI allocate failed')
  end if

  allocate (gr_mpoleScratch (1:gr_mpoleMaxLM,1:gr_mpoleMaxQ  ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleScratch allocate failed')
  end if
!
!
!       ...Allocate the radial multi-bin radii array.
!
!
  allocate (gr_mpoleQRadii (0:gr_mpoleMaxQ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQRadii allocate failed')
  end if
!
!
!       ...Allocate the radial multi-bin damping arrays.
!
!
  allocate (gr_mpoleQDampingR (0:gr_mpoleMaxQ  ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQDampingR allocate failed')
  end if

  allocate (gr_mpoleQDampingI (1:gr_mpoleMaxQ+1), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQDampingI allocate failed')
  end if
!
!
!       ...Allocate the radial multi-bin monitoring array.
!
!
  allocate (gr_mpoleQused (1:gr_mpoleMaxQ), stat = status)

  if (status > 0) then
      call Driver_abort ('[gr_mpoleAllocateRadialArrays] ERROR: gr_mpoleQused allocate failed')
  end if
!
!
!       Done.
!
!
  return
end subroutine gr_mpoleAllocateRadialArrays
