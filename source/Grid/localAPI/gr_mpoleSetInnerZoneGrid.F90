!!****if* source/Grid/localAPI/gr_mpoleSetInnerZoneGrid
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
!!  gr_mpoleSetInnerZoneGrid
!!
!! 
!! SYNOPSIS
!!
!!  call gr_mpoleSetInnerZoneGrid (integer (in)    :: nRlocal,
!!                                 integer (in)    :: nRinnerZone,
!!                                 integer (in)    :: nPinnerZone,
!!                                 real    (inout) :: RinnerZone (1:nRinnerZone))
!!
!! DESCRIPTION
!!
!!  This routine sets up the inner zone radial grid from all local inner zone
!!  radii on each processor. When exiting this routine, all processors will have
!!  a copy of the inner zone grid and its associated data.
!!
!! ARGUMENTS
!!
!!  nRlocal     : number of inner zone radii on the local processor
!!  nRinnerZone : total number of inner zone radii
!!  nPinnerZone : total number of processors containing the inner zone
!!  RinnerZone  : the collection of all inner zone radii
!!
!!***

subroutine gr_mpoleSetInnerZoneGrid (nRlocal,     &
                                     nRinnerZone, &
                                     nPinnerZone, &
                                     RinnerZone   )

  implicit none

  integer, intent (in)    :: nRlocal
  integer, intent (in)    :: nRinnerZone
  integer, intent (in)    :: nPinnerZone
  real,    intent (inout) :: RinnerZone (1:nRinnerZone)

  RinnerZone = 0.0

  return
end subroutine gr_mpoleSetInnerZoneGrid
