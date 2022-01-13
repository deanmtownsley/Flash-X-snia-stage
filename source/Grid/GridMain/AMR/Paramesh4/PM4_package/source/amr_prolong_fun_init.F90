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

!!****f* source/amr_prolong_fun_init
!! NAME
!!
!!   amr_prolong_fun_init
!!
!! SYNOPSIS
!!
!!   Call amr_prolong_fun_init()
!!
!! ARGUMENTS
!!
!!   No arguments.
!!
!! INCLUDES
!!
!! USES
!!
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_prolong_cc_fun_init
!!   amr_prolong_face_fun_init
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine calls the routines which compute the values of dx,dy and 
!!   dz and some index vectors used during the interpolation process. 
!!   These are used inside the prolongation routines
!!   saving needless repetitive computation at the cost of minimal storage
!!   space.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          July 1997
!!
!!***

      Subroutine amr_prolong_fun_init

      use paramesh_interfaces, only : amr_prolong_cc_fun_init,         & 
                                      amr_prolong_face_fun_init

      Implicit None

!-----Begin executable code.

      Call amr_prolong_cc_fun_init

      Call amr_prolong_face_fun_init

      Return
    End Subroutine amr_prolong_fun_init
