!!****if* source/Grid/GridMain/AMR/gr_sortByWorkTools.F90
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
!!
!!  gr_sortByWorkTools
!!
!! SYNOPSIS
!!
!!  USE gr_sortByWorkTools, ONLY: gr_sortByWorkConsolidated, &
!!                                gr_sortByWorkDistributed,  &
!!                                gr_calcProcLocblk
!!
!! DESCRIPTION
!!
!!  This module contains some tools for assigning subsequences of
!!  blocks ordered by Morton numbers to procs.
!!
!!  The idea is to "sort" blocks into buckets in such a way that work,
!!  represented by a weight factor assigned to each block, gets
!!  distributed approximately evenly across buckets.
!!***

module gr_sortByWorkTools
  implicit none

contains

!!****if* source/Grid/GridMain/AMR/gr_sortByWorkTools.F90
!!
!! NAME
!!
!!  gr_sortByWorkConsolidated
!!
!! SYNOPSIS
!!
!!  call gr_sortByWorkConsolidated(nprocs,mype,work,&
!!                                 totblocks,mort_cutoffs)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  integer,intent(in) :: nprocs, mype, totblocks
!!  real,intent(in) :: work(totblocks)
!!  integer,intent(inout) :: mort_cutoffs(nprocs)
!!
!!***

    subroutine gr_sortByWorkConsolidated(nprocs,mype,work,&
                                           totblocks,mort_cutoffs)

      implicit none

!-----Input/output statements
      integer,intent(in)     :: nprocs, mype, totblocks
      real,intent(in)        :: work(totblocks)
      integer,intent(inout)  :: mort_cutoffs(nprocs)

!-----Local variables and arrays
      real    :: tot_work, work_per_proc, work_so_far
      integer :: current_proc, current_block

!-----Initialize to ensure localMortUB(i) <= localMortUB(i+1)
      mort_cutoffs(:) = totblocks

      tot_work = sum(work(1:totblocks))
      work_per_proc = tot_work/real(nprocs)
      work_so_far = 0.
      current_proc = 0
      current_block = 1

!-----Apply algorithm to update mort_cutoffs
      call gr_sortAlgorithm(nprocs,work_per_proc,&
                          work_so_far,current_proc,&
                          current_block,totblocks,totblocks,&
                          work,mort_cutoffs)

    end subroutine gr_sortByWorkConsolidated

!!****if* source/Grid/GridMain/AMR/gr_sortByWorkTools.F90
!!
!! NAME
!!
!!  gr_sortByWorkDistributed
!!
!! SYNOPSIS
!!
!!  call gr_sortByWorkDistributed(nprocs,mype,worktemp,&
!!                                lnblockst,mort_cutoffs)
!!
!! DESCRIPTION
!!
!!
!!  Unlike the above routine, this routine optimizes for memory.
!!  Instead of storing the work of every single block in existence,
!!  this routine only allocates enough memory to store the work
!!  of blocks on each proc in sequence.
!!
!! ARGUMENTS
!!
!!  integer,intent(in) :: nprocs, mype, lnblockst, buflen
!!  real,intent(in) :: worktemp(buflen)
!!  integer,intent(inout) :: mort_cutoffs(nprocs)
!!
!!***

    subroutine gr_sortByWorkDistributed(nprocs,mype,worktemp,&
                                          lnblockst,mort_cutoffs,buflen)

!-----Use statements
      use Grid_data, ONLY : gr_meshComm

#include "Flashx_mpi_implicitNone.fh"

!-----Input/output statements
      integer,intent(in)    :: nprocs, mype, lnblockst, buflen
      real,intent(in)       :: worktemp(buflen)
      integer,intent(inout) :: mort_cutoffs(nprocs)

!-----Local variables and arrays
      integer :: j,ierr
      integer :: totblocks
      integer :: lnblockst_arr(0:nprocs-1)
      real    :: worktempbuf(buflen)
      real    :: locwork, totwork
      real    :: work_so_far, work_per_proc
      integer :: current_proc, current_block


!-----Gather list of lnblockst for each proc
      call MPI_ALLGATHER(lnblockst,1,MPI_INTEGER,  &
                         lnblockst_arr,1,MPI_INTEGER,    &
                         gr_meshComm,ierr)

      totblocks = sum(lnblockst_arr)

!-----Sum total work across procs
      locwork = sum(worktemp(1:lnblockst))
      call MPI_ALLREDUCE(locwork,totwork,1,FLASH_REAL, &
                         MPI_SUM,gr_meshComm,ierr)


!-----Initialize to ensure mort_cutoffs(i) <= mort_cutoffs(i+1)
      mort_cutoffs(:) = totblocks

      work_per_proc = totwork/real(nprocs)
      work_so_far = 0.
      current_proc = 0
      current_block = 1
      do j=0,nprocs-1
!-------Proc-by-proc broadcast work list from worktempbuf
        if (j.eq.mype) then
          worktempbuf(1:lnblockst_arr(j)) = worktemp(1:lnblockst_arr(j))
        end if

        call MPI_BCAST(worktempbuf,lnblockst_arr(j),FLASH_REAL,  &
                       j,gr_meshComm,ierr)

!-------Apply algorithm to update mort_cutoffs
        call gr_sortAlgorithm(nprocs,work_per_proc,&
                            work_so_far,current_proc,&
                            current_block,totblocks,&
                            lnblockst_arr(j),worktempbuf,&
                            mort_cutoffs)

      end do !0,nprocs-1

    end subroutine gr_sortByWorkDistributed

!!****
!!
!! NAME
!!
!!  gr_sortAlgorithm
!!
!! SYNOPSIS
!!
!!  call gr_sortAlgorithm(nprocs,work_per_proc,work_so_far,
!!                      current_proc,current_block,totblocks,
!!                      nblocks,workbuf,mort_cutoffs)
!!
!! DESCRIPTION
!!
!!  Routine that manages the sort-by-work.
!!
!! ARGUMENTS
!!
!!  real,intent(in)       :: work_per_proc
!!  real,intent(inout)    :: work_so_far
!!  integer,intent(inout) :: current_proc,current_block
!!  integer, intent(in)   :: totblocks,nblocks
!!  real,intent(in)       :: workbuf(totblocks)
!!  integer,intent(out)   :: mort_cutoffs(nprocs)
!!
!!***

    subroutine gr_sortAlgorithm(nprocs,work_per_proc,&
                                work_so_far,current_proc,&
                                current_block,totblocks,&
                                nblocks,workbuf,mort_cutoffs)
      
      implicit none

!-----Input/output statements
      integer,intent(in)    :: nprocs
      real,intent(in)       :: work_per_proc
      real,intent(inout)    :: work_so_far
      integer,intent(inout) :: current_proc,current_block
      integer, intent(in)   :: totblocks,nblocks
      real,intent(in)       :: workbuf(totblocks)
      integer,intent(out)   :: mort_cutoffs(nprocs)

!-----Local variables and arrays
      integer :: i
      integer :: proc_num

#ifdef IMPROVED_SORT
      do i=1,nblocks
        work_so_far = work_so_far + workbuf(i)
        proc_num = int((work_so_far - 1.)/work_per_proc)
        if((proc_num.gt.current_proc.or.&
           totblocks-current_block.le.nprocs-(current_proc+1))&
           .and.(current_proc+2).le.nprocs ) then
          mort_cutoffs(current_proc+1) = current_block-1
          current_proc = current_proc + 1
        end if
        current_block = current_block + 1
      end do
#else
      do i=1,nblocks
        work_so_far = work_so_far + workbuf(i)
        proc_num = int((work_so_far - 1.)/work_per_proc)
        if(proc_num.gt.current_proc) then
          if(proc_num.gt.current_proc+1) then
            mort_cutoffs(current_proc+1:proc_num) = current_block-1
          else
            mort_cutoffs(current_proc+1) = current_block-1
          end if
          current_proc = max(min(proc_num,nprocs-1),0)
        end if
        current_block = current_block + 1
      end do
#endif

    end subroutine gr_sortAlgorithm

!!****if* source/Grid/GridMain/AMR/gr_sortByWorkTools.F90
!!
!! NAME
!!
!!  gr_calcProcLocblk
!!
!! SYNOPSIS
!!
!!  call gr_calcProcLocblk(nprocs, cutoffs, mort, proc, locblk)
!!
!! DESCRIPTION
!!
!!  Routine that calculates proc and local block id from global
!!  morton number (index) and a list of morton index cutoffs by proc.
!!
!! ARGUMENTS
!!
!!  integer, intent(in) :: nprocs
!!  integer, intent(in) :: cutoffs(nprocs)
!!  integer, intent(in) :: mort
!!  integer, intent(out) :: proc
!!  integer, intent(out) :: locblk
!!
!!***

    subroutine gr_calcProcLocblk(nprocs, cutoffs, mort, proc, locblk)
      
      implicit none

      integer, intent(in) :: nprocs
      integer, intent(in) :: cutoffs(nprocs)
      integer, intent(in) :: mort
      integer, intent(out) :: proc
      integer, intent(out) :: locblk

      integer :: j

      proc = -1
      locblk = -1
      procloop: do j=1,nprocs
        if (j.eq.1) then
          if(mort.le.cutoffs(j)) then
            proc = j-1
            locblk = mort
            exit procloop
          end if
        else
          if(mort.gt.cutoffs(j-1).AND.mort.le.cutoffs(j) ) then
            proc = j-1
            locblk = mort - cutoffs(j-1)
            exit procloop
          end if
        end if
      end do procloop

      return
    end subroutine gr_calcProcLocblk
end module gr_sortByWorkTools
