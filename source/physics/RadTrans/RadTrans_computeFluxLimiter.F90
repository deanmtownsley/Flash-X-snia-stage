!!****f* source/physics/RadTrans/RadTrans_computeFluxLimiter
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
!!  RadTrans_computeFluxLimiter
!!
!! SYNOPSIS
!!
!!  call RadTrans_computeFluxLimiter(integer(in) :: ifl,
!!                                   integer(in) :: iflOut,
!!                                   integer(in) :: ieddi3,
!!                                   real(INOUT) :: solnData,
!!                                   integer(IN) :: blockID,
!!                                   integer(IN),OPTIONAL :: gcLayers)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ifl : 
!!
!!   iflOut : 
!!
!!   ieddi3 : 
!!
!!   solnData : 
!!
!!   blockID : ID of block in current processor
!!
!!   gcLayers : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine RadTrans_computeFluxLimiter(ifl, iflOut, ieddi3, solnData, blockID,gcLayers)
  implicit none

  integer, intent(in) :: ifl
  integer, intent(in) :: iflOut
  integer, intent(in) :: ieddi3
  real,    intent(INOUT) :: solnData(:,1:,1:,1:)
  integer, intent(IN) :: blockID
  integer, intent(IN),OPTIONAL :: gcLayers

end subroutine RadTrans_computeFluxLimiter
