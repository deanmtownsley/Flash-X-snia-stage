!----------------------------------------------------------------------
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/user_coord_transfm
!! NAME
!!
!!   user_coord_transfm
!!
!! SYNOPSIS
!!
!!   Call user_coord_transfm(x,y)
!!   Call user_coord_transfm(integer,integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(inout) :: x, y
!!
!! INCLUDES
!!
!! USES
!!
!! CALLS
!!
!! RETURNS
!!
!! DESCRIPTION
!!  
!!   A user defined routine to perform some coordinate transformation.
!!
!! AUTHORS
!!
!!   You, the user of PARAMESH.
!!
!!***

         Subroutine user_coord_transfm(x,y)

         Integer,Intent(inout) :: x, y

         Return
         End Subroutine user_coord_transfm
