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

      subroutine amr_get_new_nodetypes (nprocs,mype,level)
      
      use mpi_morton
      
      integer, intent(in) :: nprocs, mype, level
      integer :: tag_offset
      logical :: lec, lnc, lfulltree
      
      call mpi_amr_read_guard_comm_mg(nprocs,level)
      call mpi_amr_read_prol_comm_mg(nprocs,level)
      call mpi_amr_read_flux_comm_mg(nprocs,level)
      call mpi_amr_read_restrict_comm_mg(nprocs,level)
      
      return
      end subroutine amr_get_new_nodetypes
