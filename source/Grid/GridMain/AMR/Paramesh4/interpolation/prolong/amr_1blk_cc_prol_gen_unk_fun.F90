#include "paramesh_preprocessor.fh"
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


      subroutine amr_1blk_cc_prol_gen_unk_fun &
     &  (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
     &   mype,lb,pe_p,lb_p)


!------------------------------------------------------------------------
!
! This routine is a wrapper routine which calls the functions
! which prolong data for UNK.
!
!------------------------------------------------------------------------

      use timings, ONLY: timing_mpi, timer_amr_1blk_cc_prol_gen_unk

      implicit none

!------------------------------------

      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)


      include 'mpif.h'
      double precision :: time1


!------------------------------------

! local variables

      integer :: ivar

!------------------------------------

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

! Call the minimally changed subroutine for Paramesh2
      call amr_prolong_gen_unk1_fun &
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
     &     mype,lb)
      

      if (timing_mpi) then
              timer_amr_1blk_cc_prol_gen_unk =  &
     &             timer_amr_1blk_cc_prol_gen_unk &
     &                          + mpi_wtime() - time1
      endif

      return
      end subroutine amr_1blk_cc_prol_gen_unk_fun
