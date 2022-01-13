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

      subroutine amr_mg_prolong (nprocs, mype2, level)

      use tree
      use amr_mg_common

      implicit none

      integer, intent(in) :: nprocs, mype2, level
      integer             :: lb

      newchild(:) = .FALSE.
      do lb = 1, lnblocks

         nodetype(lb) = nodetype_old(lb)
         if (lrefine(lb) > level) nodetype(lb) = -1
         if (lrefine(lb) == level) then
           nodetype(lb) = 1
           newchild(lb) = .TRUE.
         endif
         if (lrefine(lb) == level-1 .and. nodetype(lb) /= 1)  & 
     &       nodetype(lb) = 2

      end do

      call amr_get_new_nodetypes (nprocs, mype2, level)

!               Call the PARAMESH prolongation routine.

      call amr_prolong (mype2, 2, 1) 

      return
      end subroutine amr_mg_prolong
