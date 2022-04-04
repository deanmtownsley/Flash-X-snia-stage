!!****f* source/Grid/Grid_notifySolnDataUpdate
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
!!  Grid_notifySolnDataUpdate
!!
!! SYNOPSIS
!!
!!  call Grid_notifySolnDataUpdate(
!!                        OPTIONAL,integer(in) :: gds,
!!                        OPTIONAL,logical(in) :: mask(*))
!!  call Grid_notifySolnDataUpdate(
!!                                 integer(in) :: varList(:),
!!                        OPTIONAL,integer(in) :: gds)
!!
!! DESCRIPTION
!!
!!  Notify the Grid unit that the calling unit intends to update
!!  (or has updated) some solution variables.
!!
!!  This is only necessary, and only makes sense, when a promise has been
!!  made that solution variables would *not* be changed. Such promises are
!!  made when Grid_fillGuardCells is called with unitReadsMeshDataOnly=.TRUE.,
!!  and are valid until either this interface is called, or
!!  Grid_fillGuardCells is called without unitReadsMeshDataOnly=.TRUE. .
!!
!!  As an effect of this notification, skipping of the next guard fill (the
!!  next time Grid_fillGuardCells is called, usually from a different unit) is
!!  canceled.  (Such skipping may have been enabled by an earlier call to
!!  Grid_fillGuardCells with unitReadsMeshDataOnly=.TRUE..) Currently
!!  (2013-02-13), this is the only effect.
!!
!! ARGUMENTS
!!
!!   gds : "grid data struct", i.e., probably one of CENTER,
!!         CENTER_FACES, FACES.  Use the same here as you would for
!!         the 'gridDataStruct' dummy argument of Grid_fillGuardCells.
!!         (CURRENTLY INGORED)
!!
!!   mask: a logical mask, intended to be of the same structure as for
!!         Grid_fillGuardCells
!!         (CURRENTLY INGORED)
!!
!!   varList: a list of UNK (or other) variables
!!           (CURRENTLY INGORED)
!!
!! SIDE EFFECTS
!!
!!  Resets gr_gcellsUpToDate flag internal to the Grid unit.
!!
!! SEE ALSO
!!  Grid_fillGuardCells
!!
!! NOTES
!!  
!!  All arguments are currently ignored. (2022-03-04)
!!
!!  This call has an effect only with PARAMESH, since only the
!!  Paramesh4 implementation of Grid_fillGuardCells recognizes the internal
!!  gr_gcellsUpToDate flag.
!!
!!  Grid_notifySolnDataUpdate is only relevant when there are code units that
!!  call Grid_fillGuardCells with the dummy argument unitReadsMeshDataOnly=
!!  .TRUE.; currently (2022-03-04), this is controlled with the
!!  "reduceGcellFills" runtime parameter in all unit implementations included
!!  with Flash-X.
!!
!!  A unit that updates solution data should either call Grid_fillGuardCells
!!  (without unitReadsMeshDataOnly=.TRUE., of course) or
!!  Grid_notifySolnDataUpdate; it need not do both, since an intent to update
!!  solution variables is considered implied by a call to Grid_fillGuardCells
!!  without unitReadsMeshDataOnly=.TRUE.
!!
!!  If a unit only updates variables that it considers "private" to itself,
!!  it may skip calling Grid_notifySolnDataUpdate. A variable can only be
!!  considered private in this sense if no unit depends on its guard cells
!!  being filled.
!!***

subroutine Grid_notifySolnDataUpdate(gds,mask)
  implicit none
  integer,OPTIONAL,intent(in) :: gds     !"grid data struct", i.e., CENTER, CENTER_FACES, FACES
  logical,OPTIONAL,intent(in) :: mask(*) !optional mask, as for Grid_fillGuardCells
  return
end subroutine Grid_notifySolnDataUpdate

subroutine Grid_notifySolnDataUpdateVlist(varList,gds)
  implicit none
  integer,intent(in)          :: varList(:)     ! list of UNK (or other?) variables
  integer,OPTIONAL,intent(in) :: gds     !"grid data struct", i.e., CENTER, CENTER_FACES, FACES
  return
end subroutine Grid_notifySolnDataUpdateVlist
