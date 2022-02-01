!!****f* source/Grid/GridMain/AMR/Paramesh4/bittree/Grid_setWork
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
!!  Grid_setWork
!!
!! SYNOPSIS
!!
!!  Grid_setWork(Grid_tile_t(IN) :: tileDesc
!!               real(IN)        :: work   )
!!               integer(IN)     :: mode   )
!!  
!! DESCRIPTION 
!!  Sets the work value for a given block. Can pass an operation name (defined in
!!  constants.h) to apply an operation to the work instead of overwriting the
!!  previous value. Used for uneven load distributions when using Paramesh.
!!
!! ARGUMENTS
!!  tileDesc - block tile.
!!  work - desired weight in the work distribution.
!!  mode - (optional) If supplied, specifies the operation to apply the given work.
!!
!!***
      subroutine Grid_setWork(tileDesc, work, mode)
      use Grid_tile, ONLY : Grid_tile_t

      implicit none
      type(Grid_tile_t),intent(in) :: tileDesc
      real,intent(in)              :: work
      integer,intent(in),optional  :: mode

      !implementation in Paramesh

      return
      end subroutine Grid_setWork
