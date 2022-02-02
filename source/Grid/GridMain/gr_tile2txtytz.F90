!!****if* source/Grid/GridMain/Grid_getListOfTiles
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
!!  Grid_getListOfTiles
!!
!! NOTES
!!
!!  Maybe this soubroutine is not needed in this standalone form;
!!  the Grid_iterator can contain a private subroutine tile2txtytz instead.
!!  That probably makes more sense since nxt,nyt,nzt specific to an
!!  iterator build should be used rather than from gr_tilePolicyData.
!!***

subroutine gr_tile2txtytz(tileID,tx,ty,tz)
  use gr_tilePolicyData,ONLY: nxt=>gr_tileNxt,nyt=>gr_tileNyt,nzt=>gr_tileNzt
  implicit none
  integer,intent(IN) :: tileID
  integer,intent(OUT) :: tx,ty,tz

  integer :: tn
!!$  tn = mod(tileID,100000)
  tn = tileID

  tx = mod(tn,nxt)
  ty = mod(tn/nxt,nyt)
  tz = tn/(nxt*nyt)
end subroutine gr_tile2txtytz

