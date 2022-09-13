!!****f* source/numericalTools/MoL/MoL_registerPostUpdate
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
!!      MoL_registerPostUpdate
!!
!!  SYNOPSIS
!!
!!      call MoL_registerPostUpdate(integer, intent(in)         :: postUpdateType,
!!                                  procedure(MoL_postUpdate_t) :: postUpdateFunc)
!!
!!  DESCRIPTION
!!
!!      Register a post-update function with MoL
!!
!! ARGUMENTS
!!
!!      postUpdateType : One of the following function types defined in MoL.h
!!                          MOL_POST_UPDATE       -  Post-update (slow) per-stage
!!                          MOL_POST_UPDATE_FAST  -  Post-update (fast) per-stage
!!      postUpdateFunc : Procedure of type MoL_postUpdate_t to register
!!
!!***
subroutine MoL_registerPostUpdate(postUpdateType, postUpdateFunc)
   use MoL_functionTypes, only: MoL_postUpdate_t

   implicit none

   integer, intent(in) :: postUpdateType
   procedure(MoL_postUpdate_t) :: postUpdateFunc

   return
end subroutine MoL_registerPostUpdate
