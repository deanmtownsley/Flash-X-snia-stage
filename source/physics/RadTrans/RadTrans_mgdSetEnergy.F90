!!****f* source/physics/RadTrans/RadTrans_mgdSetEnergy
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
!!  NAME 
!!
!!  RadTrans_mgdSetEnergy
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdSetEnergy( integer(IN) :: blockId,
!!                              integer(IN) :: axis(MDIM),
!!                              integer(IN) :: grpNum,
!!                              real(IN)    :: eg )
!!
!!  DESCRIPTION 
!!
!!      Set the specific energy for a particular energy group in a
!!      particular cell. This routine has been created to make it easy
!!      for users to specify the energy in a group. This can be a
!!      little complicated because of mesh replication - but all of
!!      the details are handled internally in RadTrans
!!
!! ARGUMENTS
!!
!!    blockId : The blockId of the cell
!!    axis    : An array storing the i,j,k coordinate of the cell
!!    grpNum  : The energy group number
!!    eg      : The specific internal energy to use [ergs/g]
!! 
!!***
subroutine RadTrans_mgdSetEnergy(blockId, axis, grpNum, eg)
  implicit none

#include "constants.h"

  integer, intent(in) :: blockId
  integer, intent(in) :: axis(MDIM)
  integer, intent(in) :: grpNum
  real,    intent(in) :: eg

  return

end subroutine RadTrans_mgdSetEnergy
