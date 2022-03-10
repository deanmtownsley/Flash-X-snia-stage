!!****if* source/IO/IOMain/hdf5/serial/PM/io_readData
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
!!  io_readData
!!
!!
!! SYNOPSIS
!!
!!  call io_readData()
!!
!!
!! DESCRIPTION
!!
!!  This is the reading counterpart to io_writeData.  It reads an HDF5
!!  file and distributes it to the processors to restart a simulation.
!!
!!  this is the reading counterpart to checkpoint_wr.  It eats the HDF
!!  output and distributes it to the processors.
!!
!!  Subroutine to read checkpoint file using AMR package.
!!  Currently reads are done serially by processor 0 and data sent to
!!  other processors.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***

#include "constants.h"
#include "Simulation.h"
#include "io_flash.h"

subroutine io_readData()

   use io_typeInterface, ONLY: io_xfer_tree_data

   use IO_data, ONLY: io_globalMe, io_globalNumProcs, io_globalComm, &
                      io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
                      io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
                      io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
                      io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
                      io_realScalarNames, io_realScalarValues, io_numRealScalars, &
                      io_intScalarNames, io_intScalarValues, io_numIntScalars, &
                      io_logScalarNames, io_logScalarValues, io_numLogScalars, &
                      io_strScalarNames, io_strScalarValues, io_numStrScalars, &
                      io_logToIntScalarValues, io_logToIntParmValuesPrev, io_unklabels, &
                      io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
                      io_baseName, io_checkpointFileNumber, io_outputSplitNum, io_comm, io_chkptFileID, &
                      io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels, &
                      tree_data_t, io_splitNumBlks, io_meshMe, io_meshNumProcs

   use Driver_interface, ONLY: Driver_abort
   use RuntimeParameters_interface, ONLY: RuntimeParameters_bcast
   use Logfile_interface, ONLY: Logfile_stamp
   use Grid_interface, ONLY: Grid_putLocalNumBlks, Grid_receiveInputData, Grid_sendOutputData
   use IO_interface, ONLY: IO_getScalar
   use Particles_interface, ONLY: Particles_createDataStructs
   use Simulation_interface, ONLY: Simulation_initBlock
   use gr_specificData, ONLY: gr_nToLeft, gr_globalOffset
   use gr_amrexInterface, ONLY: gr_fillPhysicalBC, gr_preinterpolationWork, &
                                gr_postinterpolationWork
   use Grid_data, ONLY: gr_globalNumBlocks, gr_globalDomain, gr_interpolator, &
                        lo_bc_amrex, hi_bc_amrex, gr_eosModeInit, &
                        gr_maxRefine, gr_doFluxCorrection
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t
   use Grid_interface, ONLY: Grid_getTileIterator, &
                             Grid_releaseTileIterator
   use gr_physicalMultifabs, ONLY: unk, &
                                   gr_scratchCtr, &
                                   fluxes, &
                                   flux_registers, facevarx, facevary, facevarz
   use amrex_amrcore_module, ONLY: amrex_set_boxarray, amrex_get_numlevels, amrex_get_boxarray, &
                                   amrex_set_distromap, amrex_set_finest_level, amrex_geom, &
                                   amrex_ref_ratio, amrex_get_distromap
   use amrex_fort_module, ONLY: wp => amrex_real
   use amrex_box_module, ONLY: amrex_box
   use amrex_boxarray_module, ONLY: amrex_boxarray_build, &
                                    amrex_boxarray, amrex_print
   use amrex_distromap_module, ONLY: amrex_distromap_build, amrex_distromap, &
                                     amrex_print
   use amrex_multifab_module, ONLY: amrex_multifab, amrex_multifab_build
   use amrex_multifab_module, ONLY: amrex_mfiter, &
                                    amrex_mfiter_build, &
                                    amrex_mfiter_destroy
   use gr_fluxregister_mod, ONLY: gr_fluxregisterBuild
   use amrex_plotfile_module, ONLY: amrex_write_plotfile

   use iso_c_binding, ONLY: c_associated

#include "Flashx_mpi_implicitNone.fh"

   type(amrex_box) :: domain
   type(amrex_boxarray) :: ba
   type(amrex_distromap) :: dm
   type(amrex_multifab), POINTER :: unk_mf(:)

!
   type(tree_data_t) :: tree_data

! create a character variable to hold the string representation of the block
! number.  Note this is set to be 4 characters long (i.e. max = 9999).
   character(len=4) :: fnumStr
   character(len=MAX_STRING_LENGTH) :: filename
   integer :: localNumBlocks, localNumBlockst

   integer :: blockID, var, dir
   integer :: jproc, i, ii, j, k, lev, lb, numx, numy, numz, lblock
   logical :: has_b1
   real :: time
   integer, allocatable :: id_on_lev(:, :)
   logical, save :: isPtDataStructCreated = .false.

   integer :: alnblocks !avg local num blocks
   integer :: lnblockst !temp local num blocks
   integer :: ierr, numContinuousBlocks
   integer :: procBlocks
   integer :: globalOffsett
   integer :: lo(3), hi(3)
   integer, allocatable ::bnd_box_sub(:, :, :)

   integer, allocatable :: lrefine(:), levBlock(:)
   integer, allocatable :: nodetype(:), localBIds(:), gid(:, :)
   real, allocatable :: coord(:, :), bsize(:, :), bnd_box(:, :, :)
   integer, allocatable :: pmap(:), localBlockList(:), tempList(:), numCB(:), startOffset(:)

   real, allocatable :: unkBuf(:, :, :, :, :), unkTemp(:, :, :, :, :)
   real, allocatable :: faceXBuf(:, :, :, :, :)
   real, allocatable :: faceYBuf(:, :, :, :, :)
   real, allocatable :: faceZBuf(:, :, :, :, :)

   integer :: xx, yy, local, startIndex

   integer          :: myGlobalComm
   integer          :: status(MPI_STATUS_SIZE)

   character(len=4) :: recordLabel

! for logfile output
   character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:, :) :: strBuf
   character(len=16)                                    :: numToStr
   integer :: splitOffset, localOffset
   integer :: doread, rank

   integer, parameter :: presentDims = MDIM
   logical, parameter :: do_gsurr_blks_read = .false.

   integer :: min_lev, max_lev, tmp_lev, max_lev_used, temp, blockLocation
   integer :: maxBlockChunks, totBlockChunks
   logical :: nodal(1:MDIM)
   type(Grid_iterator_t) :: itor
   type(Grid_tile_t)     :: tileDesc

   real(wp), contiguous, pointer, dimension(:, :, :, :) :: dp, facexData, faceyData, facezData

   myGlobalComm = io_globalComm

   call MPI_Comm_Rank(myGlobalComm, rank, ierr)
   call io_getOutputName(io_checkpointFileNumber, "hdf5", "_chk_", filename, .false.)

   if (io_globalMe == MASTER_PE) then
      write (*, *) '[io_readData] Opening ', trim(filename), &
         ' for restart'
      write (*, FMT="(1X,A11)", ADVANCE="NO") 'Progress'
      allocate (strBuf(2, 2))
      write (strBuf(1, 1), "(A)") "type"
      write (strBuf(1, 2), "(A)") "checkpoint"
      write (strBuf(2, 1), "(A)") "name"
      write (strBuf(2, 2), "(A)") trim(filename)
      call Logfile_stamp(strBuf, 2, 2, "[io_readData] file opened")
      deallocate (strBuf)
   end if

   call io_h5open_file_for_read(io_chkptFileID, filename, io_comm, io_outputSplitNum)

   call io_h5read_header(io_globalMe, io_chkptFileID, io_unklabels, io_outputSplitNum)

!allocate space for scalar and parameter lists
   call io_prepareListsRead()

   call io_h5read_runtime_parameters(io_chkptFileID, &
                                     io_numRealParmsPrev, &
                                     io_realParmNamesPrev, &
                                     io_realParmValuesPrev, &
                                     io_numIntParmsPrev, &
                                     io_intParmNamesPrev, &
                                     io_intParmValuesPrev, &
                                     io_numStrParmsPrev, &
                                     io_strParmNamesPrev, &
                                     io_strParmValuesPrev, &
                                     io_numLogParmsPrev, &
                                     io_logParmNamesPrev, &
                                     io_logToIntParmValuesPrev)

   call io_h5read_scalars(io_chkptFileID, &
                          io_numRealScalars, &
                          io_realScalarNames, &
                          io_realScalarValues, &
                          io_numIntScalars, &
                          io_intScalarNames, &
                          io_intScalarValues, &
                          io_numStrScalars, &
                          io_strScalarNames, &
                          io_strScalarValues, &
                          io_numLogScalars, &
                          io_logScalarNames, &
                          io_logToIntScalarValues)

   call io_finalizeListsRead()
   call IO_getScalar("globalNumBlocks", gr_globalNumBlocks)
   call io_checkBlockShape(gr_globalNumBlocks)
   if (io_outputSplitNum /= 1) then
      call IO_getScalar("splitNumBlocks", io_splitNumBlks)
   else
      io_splitNumBlks = gr_globalNumBlocks
   end if

   allocate (bnd_box_sub(2, MDIM, gr_globalNumBlocks))

   allocate (tree_data%lrefine(gr_globalNumBlocks))
   allocate (tree_data%nodetype(gr_globalNumBlocks))
   allocate (tree_data%coord(MDIM, gr_globalNumBlocks))
   allocate (tree_data%bsize(MDIM, gr_globalNumBlocks))
   allocate (tree_data%bnd_box(2, MDIM, gr_globalNumBlocks))
   allocate (tree_data%gid(2*MDIM + 1 + 2**MDIM, gr_globalNumBlocks))

   allocate (levBlock(gr_globalNumBlocks))

   if (io_globalMe == MASTER_PE) then
      call io_xfer_tree_data(tree_data, io_chkptFileID, IO_FILE_HDF5, &
                             IO_READ_XFER_MASTER_PE, &
                             gr_globalNumBlocks, 0, presentDims, gr_globalNumBlocks)
   end if

! call mpi_bcast for all tree_data read of master
   call MPI_Bcast(tree_data%coord, size(tree_data%coord), FLASH_REAL, MASTER_PE, &
                  myGlobalComm, ierr)
   call MPI_Bcast(tree_data%bnd_box, size(tree_data%bnd_box), FLASH_REAL, MASTER_PE, &
                  myGlobalComm, ierr)
   call MPI_Bcast(tree_data%lrefine, size(tree_data%lrefine), FLASH_INTEGER, MASTER_PE, &
                  myGlobalComm, ierr)
   call MPI_Bcast(tree_data%bsize, size(tree_data%bsize), FLASH_REAL, MASTER_PE, &
                  myGlobalComm, ierr)
   call MPI_Bcast(tree_data%gid, size(tree_data%gid), FLASH_INTEGER, MASTER_PE, &
                  myGlobalComm, ierr)

! find min an max level
   min_lev = MINVAL(tree_data%lrefine)
   max_lev = MAXVAL(tree_data%lrefine)

! gr_maxRefine can be lower/higher set from par file
   max_lev_used = MIN(max_lev, gr_maxRefine)

! allocate physical multifabs
   allocate (unk(0:gr_maxRefine - 1))

#if NFACE_VARS > 0
   allocate (facevarx(0:gr_maxRefine - 1))
#if NDIM >= 2
   allocate (facevary(0:gr_maxRefine - 1))
#endif
#if NDIM == 3
   allocate (facevarz(0:gr_maxRefine - 1))
#endif
#endif

#if 1
   allocate (gr_scratchCtr(0:gr_maxRefine - 1))
#if NFLUXES > 0
# ifdef USE_LEVELWIDE_FLUXES
   allocate (fluxes(0:gr_maxRefine - 1, 1:NDIM))
# endif
   if (gr_doFluxCorrection) then
      allocate (flux_registers(1:gr_MaxRefine))
   end if
#endif
#endif

! Store the relative index of a block on its own level
   allocate (id_on_lev(max_lev_used, gr_globalNumBlocks))
   nullify (dp, facexData, faceyData, facezData)

   lblock = 0
   allocate (localBlockList(gr_globalNumBlocks))
   allocate (tempList(gr_globalNumBlocks))

   localBlockList = -1

! i is 1 based
   do i = min_lev, max_lev_used
      k = 0
      bnd_box_sub(:, :, :) = 0

      do j = 1, gr_globalNumBlocks
         if ((tree_data%lrefine(j)) == i) then
            k = k + 1
            id_on_lev(i, k) = j !store global block id of level i's block k

            if (NDIM >= 1) bnd_box_sub(1, 1, k) = NXB*int((tree_data%bnd_box(1, 1, j) &
                                                           - gr_globalDomain(LOW, 1))/tree_data%bsize(1, j))
            if (NDIM >= 2) bnd_box_sub(1, 2, k) = NYB*int((tree_data%bnd_box(1, 2, j) &
                                                           - gr_globalDomain(LOW, 2))/tree_data%bsize(2, j))
            if (NDIM == 3) bnd_box_sub(1, 3, k) = NZB*int((tree_data%bnd_box(1, 3, j) &
                                                           - gr_globalDomain(LOW, 3))/tree_data%bsize(3, j))

            !bnd_box_sub(HI,:,k) needs to be decremented by 1 since it should be corner of top right cell
            if (NDIM >= 1) bnd_box_sub(2, 1, k) = NXB*int((tree_data%bnd_box(2, 1, j) &
                                                           - gr_globalDomain(LOW, 1))/tree_data%bsize(1, j)) - 1
            if (NDIM >= 2) bnd_box_sub(2, 2, k) = NYB*int((tree_data%bnd_box(2, 2, j) &
                                                           - gr_globalDomain(LOW, 2))/tree_data%bsize(2, j)) - 1
            if (NDIM == 3) bnd_box_sub(2, 3, k) = NZB*int((tree_data%bnd_box(2, 3, j) &
                                                           - gr_globalDomain(LOW, 3))/tree_data%bsize(3, j)) - 1
         end if
      end do

      !necessary for construction of iterator
      call amrex_set_finest_level(i - 1)
      ! build and set boxarray
      call amrex_boxarray_build(ba, bnd_box_sub(:, :, 1:k))
      call amrex_set_boxarray(i - 1, ba)
      ! build and set distromap
      call amrex_distromap_build(dm, ba)
      allocate (pmap(k))
      ! get pmap, Note: this is a new procedure in AMReX
      call dm%get_pmap(pmap)

      ! use pmap for this level to populate local block list for this rank
      ! localBlocklist
      do local = 1, k
         if (rank .EQ. pmap(local)) then
            lblock = lblock + 1
            localBlockList(lblock) = id_on_lev(i, local)
         end if
      end do
      deallocate (pmap)

!   call amrex_print(dm)
      call amrex_set_distromap(i - 1, dm)

      ! Build multifabs
      call amrex_multifab_build(unk(i - 1), ba, dm, NUNK_VARS, NGUARD)

#if NFACE_VARS > 0
      ! Face variables
      nodal(:) = .FALSE.
      nodal(IAXIS) = .TRUE.
      call amrex_multifab_build(facevarx(i - 1), ba, dm, NFACE_VARS, NGUARD, nodal)
#if NDIM >= 2
      nodal(:) = .FALSE.
      nodal(JAXIS) = .TRUE.
      call amrex_multifab_build(facevary(i - 1), ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#if NDIM == 3
      nodal(:) = .FALSE.
      nodal(KAXIS) = .TRUE.
      call amrex_multifab_build(facevarz(i - 1), ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#endif

#if 1
#if NSCRATCH_CENTER_VARS > 0
      call amrex_multifab_build(gr_scratchCtr(i - 1), ba, dm, NSCRATCH_CENTER_VARS, 0)
#endif

#if NFLUXES > 0
# ifdef USE_LEVELWIDE_FLUXES
      ! No need to store fluxes for guardcells
      do dir = 1, SIZE(fluxes, 2)
         nodal(:) = .FALSE.
         nodal(dir) = .TRUE.
         call amrex_multifab_build(fluxes(i - 1, dir), ba, dm, NFLUXES, 0, &
                                   nodal=nodal)
      end do
# endif

      if ((i > 1) .AND. (gr_doFluxCorrection)) then
#    ifdef USE_AMREX_FLASHFLUXREGISTER
         call gr_fluxregisterBuild(flux_registers(i - 1), fba=ba, cba=unk(i - 2)%ba, fdm=dm, cdm=unk(i - 2)%dm, &
                                   fine_lev=i - 1, ncomp=NFLUXES)
#    else
         call gr_fluxregisterBuild(flux_registers(i - 1), ba, dm, &
                                   amrex_ref_ratio(i - 2), &
                                   i - 1, NFLUXES)
#    endif
      end if
#endif
#endif

   end do !i=min_lev,max_lev_used

! fix the size of localBlockList from gr_globalNumBlocks to the actual number of local blocks
! Remove -1 values to get the number of total number of local blocks
   tempList = localBlockList
   localNumBlocks = 0

   do i = 1, gr_globalNumBlocks
      if (localBlockList(i) .NE. -1) then
         localNumBlocks = localNumBlocks + 1
      end if
   end do

   deallocate (localBlockList)
   allocate (localBlockList(localNumBlocks))

   do i = 1, localNumBlocks
      localBlockList(i) = tempList(i)
   end do

! can be replaced by setting operations of gr_nToLeft in Grid_sendOutput
! this call is need for gr_globalOffset - used for writing particles
   call Grid_sendOutputData()

! set the localNumblocks
   call Grid_putLocalNumBlks(localNumBlocks)

! TODO: make this an external sort function
! Sorting to get the most numContinuous blocks
   DO J = 1, UBOUND(localBlockList, 1)
      DO K = J + 1, UBOUND(localBlockList, 1)
         IF (localBlockList(J) > localBlockList(K)) THEN
            temp = localBlockList(K)
            localBlockList(K) = localBlockList(J)
            localBlockList(J) = temp
         END IF
      END DO
   END DO

! unkBuf variable is available to all procs
   allocate (unkBuf(NUNK_VARS, NXB, NYB, NZB, localNumBlocks))

#if NFACE_VARS > 0
   allocate (faceXBuf(NFACE_VARS, NXB + 1, NYB, NZB, localNumBlocks))
   faceXBuf = 0.0
#if NDIM >= 2
   allocate (faceYBuf(NFACE_VARS, NXB, NYB + 1, NZB, localNumBlocks))
   faceYBuf = 0.0
#endif
#if NDIM == 3
   allocate (faceZBuf(NFACE_VARS, NXB, NYB, NZB + 1, localNumBlocks))
   faceZBuf = 0.0
#endif
#endif

   allocate (numCB(localNumBlocks))
   allocate (startOffset(localNumBlocks))
   allocate (localBIds(localNumBlocks))
! set to zero
   unkBuf = 0.0
   numCB = 0
   startOffset = 0
   localBIds = 0.0

! init local variables for indexes of reading
   localOffset = 0
   numContinuousBlocks = 1
   startIndex = 1
   doread = 1

! call MPI_BARRIER (io_globalComm, ierr)
! print *, rank, "rank has #blocks:", localNumBlocks, "sorted pmap dist:", localBlockList

   totBlockChunks = 1
! get total block chunks and the corresponding:
! continuous blocks, start offset and local start block ids
! these are very helpful for reading later
   do i = 1, localNumBlocks - 1
      ! set the offset to global block id -1 at the beginning
      if (numContinuousBlocks .EQ. 1) then
         localOffset = localBlockList(i) - 1
      end if
      if (localBlockList(i) .EQ. (localBlockList(i + 1) - 1)) then
         numContinuousBlocks = numContinuousBlocks + 1

         if (i .EQ. (localNumBlocks - 1)) then
            numCB(totBlockChunks) = numContinuousBlocks
            startOffset(totBlockChunks) = localOffset
            localBIds(totBlockChunks) = startIndex
         end if
      else
         ! we have the numContinuousBlocks and offset, store it
         numCB(totBlockChunks) = numContinuousBlocks
         startOffset(totBlockChunks) = localOffset
         localBIds(totBlockChunks) = startIndex

         ! reset
         numContinuousBlocks = 1
         startIndex = i + 1
         ! check if this is the penultimate block
         ! if so, store the last block
         if (i .EQ. (localNumBlocks - 1)) then
            totBlockChunks = totBlockChunks + 1
            localOffset = localBlockList(i + 1) - 1
            numCB(totBlockChunks) = numContinuousBlocks
            startOffset(totBlockChunks) = localOffset
            localBIds(totBlockChunks) = startIndex
         end if
         ! increment totBlockChunks
         totBlockChunks = totBlockChunks + 1
      end if
   end do

! if there is only one local block
   if (1 .EQ. localNumBlocks) then
      numCB(1) = 1
      startOffset(1) = localBlockList(1) - 1
      localBIds(1) = 1
   end if

! get the max of totBlockChunks across all processors - maxBlockChunks
! this is the number of times io_h5read_unknowns will be called
! since it is a collective call; all ranks have to call it
! some ranks will call with no blocks to read
   CALL MPI_AllReduce(totBlockChunks, maxBlockChunks, 1, MPI_INTEGER, MPI_MAX, io_globalComm, ierr)

! now loop thru maxBlockChunks and perform the actual read from restart file
   do k = 1, maxBlockChunks
      if (k <= localNumBlocks) then
         localOffset = startOffset(k)
         numContinuousBlocks = numCB(k)
         startIndex = localBIds(k)
      else ! some processors have fewers blocks than maxBlockChunks
         localOffset = 0
         numContinuousBlocks = 0
         startIndex = 0
      end if

      do ii = UNK_VARS_BEGIN, UNK_VARS_END
         recordLabel = io_unklabels(ii)
         allocate (unkTemp(1, NXB, NYB, NZB, numContinuousBlocks))
         unkTemp = 0.0

         call io_h5read_unknowns(io_chkptFileID, &
                                 NXB, &
                                 NYB, &
                                 NZB, &
                                 unkTemp(1, :, :, :, :), &
                                 recordLabel, &
                                 numContinuousBlocks, &
                                 io_splitNumBlks, &
                                 localOffset, &
                                 doread)
         ! set unkBuf if something was read
         if (numContinuousBlocks > 0) then
            unkBuf(ii, :, :, :, startIndex:startIndex + numContinuousBlocks - 1) = unkTemp(1, :, :, :, 1:numContinuousBlocks)
         end if
         deallocate (unkTemp)
      end do

      do ii = 1, NFACE_VARS
         recordLabel = io_faceXVarLabels(ii)
         allocate (unkTemp(1, NXB + 1, NYB, NZB, numContinuousBlocks))
         unkTemp = 0.0

         call io_h5read_unknowns(io_chkptFileID, &
                                 NXB + 1, &
                                 NYB, &
                                 NZB, &
                                 unkTemp(1, :, :, :, :), &
                                 recordLabel, &
                                 numContinuousBlocks, &
                                 io_splitNumBlks, &
                                 localOffset, &
                                 doread)
         ! set unkBuf if something was read
         if (numContinuousBlocks > 0) then
            faceXBuf(ii, :, :, :, startIndex:startIndex + numContinuousBlocks - 1) = unkTemp(1, :, :, :, 1:numContinuousBlocks)
         end if
         deallocate (unkTemp)

#if NDIM >=2
         recordLabel = io_faceYVarLabels(ii)
         allocate (unkTemp(1, NXB, NYB + 1, NZB, numContinuousBlocks))
         unkTemp = 0.0

         call io_h5read_unknowns(io_chkptFileID, &
                                 NXB, &
                                 NYB + 1, &
                                 NZB, &
                                 unkTemp(1, :, :, :, :), &
                                 recordLabel, &
                                 numContinuousBlocks, &
                                 io_splitNumBlks, &
                                 localOffset, &
                                 doread)
         ! set unkBuf if something was read
         if (numContinuousBlocks > 0) then
            faceYBuf(ii, :, :, :, startIndex:startIndex + numContinuousBlocks - 1) = unkTemp(1, :, :, :, 1:numContinuousBlocks)
         end if
         deallocate (unkTemp)
#endif

#if NDIM == 3
         recordLabel = io_faceZVarLabels(ii)
         allocate (unkTemp(1, NXB, NYB, NZB + 1, numContinuousBlocks))
         unkTemp = 0.0

         call io_h5read_unknowns(io_chkptFileID, &
                                 NXB, &
                                 NYB, &
                                 NZB + 1, &
                                 unkTemp(1, :, :, :, :), &
                                 recordLabel, &
                                 numContinuousBlocks, &
                                 io_splitNumBlks, &
                                 localOffset, &
                                 doread)
         ! set unkBuf if something was read
         if (numContinuousBlocks > 0) then
            faceZBuf(ii, :, :, :, startIndex:startIndex + numContinuousBlocks - 1) = unkTemp(1, :, :, :, 1:numContinuousBlocks)
         end if
         deallocate (unkTemp)
#endif

      end do

   end do

! Loop over levels to set actual unk
   do i = min_lev, max_lev_used
! get tile iterator, for each proc it gets the blocks as per dm
      call Grid_getTileIterator(itor, ALL_BLKS, level=i, tiling=.false.)

      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         ! global block id for level i's (tileDesc%grid_index + 1) block
         j = id_on_lev(i, tileDesc%grid_index + 1)
         blockLocation = -1
         ! unkBuf has localblocks and 1 to NumLocalBlocks
         !  get the blockLocation in unk for global block j
         do local = 1, localNumBlocks
            if (localBlockList(local) .EQ. j) then
               blockLocation = local
            end if
         end do

         ! abort no block location unk isn't available to populate
         if (blockLocation .EQ. -1) then
            call Driver_abort("unk data is not read from restart file!")
         end if

         associate (lo => tileDesc%limits(LOW, :), &
                    hi => tileDesc%limits(HIGH, :))
            call tileDesc%getDataPtr(dp, CENTER)
            do var = UNK_VARS_BEGIN, UNK_VARS_END
               ! set all values at once
               dp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), var) = &
                  unkbuf(var, 1:NXB, 1:NYB, 1:NZB, blockLocation)
            end do
            call tileDesc%releaseDataPtr(dp, CENTER)

            call tileDesc%getDataPtr(facexData, FACEX)
#if NDIM >= 2
            call tileDesc%getDataPtr(faceyData, FACEY)
#endif
#if NDIM == 3
            call tileDesc%getDataPtr(facezData, FACEZ)
#endif

            do var = 1, NFACE_VARS
               facexData(lo(1):hi(1) + 1, lo(2):hi(2), lo(3):hi(3), var) = &
                  faceXBuf(var, 1:NXB + 1, 1:NYB, 1:NZB, blockLocation)
#if NDIM >=2
               faceyData(lo(1):hi(1), lo(2):hi(2) + 1, lo(3):hi(3), var) = &
                  faceYBuf(var, 1:NXB, 1:NYB + 1, 1:NZB, blockLocation)
#endif
#if NDIM == 3
               facezData(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3) + 1, var) = &
                  faceZBuf(var, 1:NXB, 1:NYB, 1:NZB + 1, blockLocation)
#endif
            end do

            call tileDesc%releaseDataPtr(facexData, FACEX)
#if NDIM >= 2
            call tileDesc%releaseDataPtr(faceyData, FACEY)
#endif
#if NDIM == 3
            call tileDesc%releaseDataPtr(facezData, FACEZ)
#endif

         end associate

         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
   end do

   deallocate (unkBuf)
#if NFACE_VARS > 0
   deallocate (faceXBuf)
#if NDIM >=2
   deallocate (faceYBuf)
#endif
#if NDIM == 3
   deallocate (faceZBuf)
#endif
#endif

   if (io_globalMe == MASTER_PE) &
      print *, 'read_data:  read ', gr_globalNumBlocks, ' blocks.'

   if (io_globalMe == MASTER_PE) then
      allocate (strBuf(2, 2))
      write (strBuf(1, 1), "(A)") "type"
      write (strBuf(1, 2), "(A)") "checkpoint"
      write (strBuf(2, 1), "(A)") "name"
      write (strBuf(2, 2), "(A)") trim(filename)
      call Logfile_stamp(strBuf, 2, 2, "[io_readData] file_closed")
      if (allocated(strBuf)) deallocate (strBuf)
   end if

   call MPI_BARRIER(io_globalComm, ierr)
   if (io_globalMe == MASTER_PE) &
      print *, 'io_readData:  finished reading input file.'

   return
end subroutine io_readData
