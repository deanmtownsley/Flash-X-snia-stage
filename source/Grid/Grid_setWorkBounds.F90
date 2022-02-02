!!****f* source/Grid/Grid_setWorkBounds
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
!!  Grid_setWorkBounds
!!
!! SYNOPSIS
!!
!!  Grid_setWorkBounds(real(IN)   :: pwork_bnd,
!!                     real(IN)   :: lwork_bnd,
!!                     integer(IN):: bnd_type)
!!
!! DESCRIPTION
!!  Sets the upper or lower bounds for work.
!!  Used for load distribution in Paramesh.
!!
!! ARGUMENTS
!!  pwork_bnd - work bound for non-leaf blocks
!!  lwork_bnd - work bound for leaf blocks
!!  bnd_type - bound type. Options = {LOW, HIGH}
!!
!!***

      subroutine Grid_setWorkBounds(pwork_bnd,lwork_bnd,bnd_type)
      implicit none
      real,intent(in)    :: pwork_bnd
      real,intent(in)    :: lwork_bnd
      integer,intent(in) :: bnd_type

      !implementation in Paramesh

      return
      end subroutine Grid_setWorkBounds
