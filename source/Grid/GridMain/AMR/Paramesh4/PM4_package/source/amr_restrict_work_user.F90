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

!!****f* source/amr_restrict_work_user
!! NAME
!!
!!   amr_restrict_work_user
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_work_user()
!!
!! ARGUMENTS
!!
!!   None or user defined.
!!
!! INCLUDES
!! 
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!! CALLS
!!
!! DESCRIPTION
!!
!!   This is a stub routine for performing the interpolation operation
!!   during restriction of data stored in 'work'. It is meant to be replaced 
!!   by a user written routine if the provided default routine is insufficient
!!   for your application.
!!
!! AUTHORS
!!
!!   YOU the users of Paramesh.
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_work_user()

      Return
      End Subroutine amr_restrict_work_user
