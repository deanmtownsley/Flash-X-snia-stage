!!****f* source/physics/RadTrans/RadTrans_sumEnergy
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
!!  RadTrans_sumEnergy
!!
!!  SYNOPSIS
!!
!!  call RadTrans_sumEnergy( integer(IN) :: ivar,
!!                           integer(IN) :: nblk,
!!                           integer(IN) :: blklst(nblk),
!!
!!  DESCRIPTION 
!!
!!    This subroutine is useful when mesh replication is active
!!    (meshCopyCount > 1). It takes an unk variable and adds it across
!!    all the meshes. For example, if each mesh had computed ERAD
!!    separately using its own energy groups, this subroutine could be
!!    used to add all of the ERADs to compute the total radiation
!!    energy.
!!
!! ARGUMENTS
!!
!!   ivar   : the unk variable to add
!!   blklst : the list of blocks to cover
!!   nblk   : the number of blocks
!!
!!***
subroutine RadTrans_sumEnergy(ivar, nblk, blklst)
  implicit none
  integer, intent(in) :: ivar
  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  ! Stub Implementation
end subroutine RadTrans_sumEnergy
     
