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
!!  io_writeData
!!
!!
!! SYNOPSIS
!!
!!  call io_writeData(integer(io_fileID_t)(in) :: fileID)
!!
!!
!! DESCRIPTION
!!
!!  This function writes the checkpoint data to an hdf5 file to store the
!!  paramesh data.  IO is done in parallel -- no copying of the data to
!!  a single processor
!!  to do the writing is performed.  HDF5 v. 1.4.0 or later is required
!!
!!  HDF5 uses MPI-IO (via ROMIO) to support parallel IO.  Each processor
!!  must open the file, create the datasets, and dataspaces for each HDF
!!  record.
!!
!!  A single record for each of the PARAMESH data structures is created.  A
!!  processor only writes to a subset of this record.  Each record has a
!!  dimension with length = globalNumBlocks.  The offset of a processor into this
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of io_writeData, each variable is given its own
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.
!!
!!
!! ARGUMENTS
!!
!!  fileID - integer file identifier for hdf5 file
!!
!!
!! NOTES
!!
!!  The KIND type parameter io_fileID_t is defined in Fortran module io_intfTypesModule.
!!  It should ensure that fileID is compatible with the hid_t of the HDF5 library version
!!  used.
!!
!!  variables that start with "io_" belong to the data module IO_data.
!!  the "io_" is meant to indicate that these variables have IO unit
!!  scope.   For performance purposes IO particularly, io_writeData
!!  and read data are allowed to access unk directly.  Additionally,
!!  io_writeData needs access to some variables that are specific to
!!  Grid.  These are stored in Grid_ioData and variables associated
!!  with Grid_ioData start with "gio_"
!!
!!***

!!R E O R D E R (5):unk,facevar[xyz],scratch

#include "constants.h"
#include "Simulation.h"
#include "io_flash.h"

subroutine io_writeData(fileID)

#ifdef USE_IO_C_INTERFACE
   use iso_c_binding, ONLY: c_loc, c_char
#endif

   use IO_data, ONLY: io_globalMe, io_globalNumProcs, io_globalComm, &
                      io_realParmNames, io_realParmValues, io_numRealParms, &
                      io_intParmNames, io_intParmValues, io_numIntParms, &
                      io_logParmNames, io_logParmValues, io_numLogParms, &
                      io_strParmNames, io_strParmValues, io_numStrParms, &
                      io_realScalarNames, io_realScalarValues, io_numRealScalars, &
                      io_intScalarNames, io_intScalarValues, io_numIntScalars, &
                      io_logScalarNames, io_logScalarValues, io_numLogScalars, &
                      io_strScalarNames, io_strScalarValues, io_numStrScalars, &
                      io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
                      io_geometry, io_setupCall, io_buildDir, io_flashRelease, &
                      io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
                      io_setupTimeStamp, io_buildTimeStamp, io_doublePrecision, &
                      io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
                      io_plotVarStr, io_nPlotVars, io_outputSplitNum, &
                      io_chkGuardCellsOutput, &
                      io_plotGridVarStr, io_faceXVarLabels, io_faceYVarLabels, &
                      io_faceZVarLabels, io_comm, io_splitNumBlks, io_plotfileMetadataDP, &
                      io_plotfileGridQuantityDP, io_fileFormatVersion, &
                      io_meshMe, io_acrossMe, io_acrossNumProcs, io_unklabelsGlobal, io_unkNonRep, &
                      tree_data_t
   use io_intfTypesModule, ONLY: io_fileID_t

   use Simulation_interface, ONLY: Simulation_mapStrToInt
   use Driver_interface, ONLY: Driver_abort
   use Logfile_interface, ONLY: Logfile_stamp
   use Grid_interface, ONLY: Grid_getLocalNumBlks, Grid_getBlkIndexLimits

   use Grid_data, ONLY: gr_globalNumBlocks
   !use gr_specificData, ONLY : gr_nToLeft, gr_globaloffset, scratch, gr_gid
   use gr_specificData, ONLY: gr_nToLeft, gr_globaloffset

   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t
   use Grid_interface, ONLY: Grid_getTileIterator, &
                             Grid_releaseTileIterator

!  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
   use gr_physicalMultifabs, ONLY: unk, facevarx, facevary, facevarz
   use gr_physicalMultifabs, ONLY: gr_scratchCtr, &
                                   fluxes, &
                                   flux_registers
   use io_typeInterface, ONLY: io_xfer_tree_data

   use gr_specificData, ONLY: gr_ioBlkNodeType, gr_ioBlkCoords, gr_ioBlkBsize, gr_ioBlkBoundBox, &
                              gr_ioBlkLrefine

   implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
   integer, parameter :: c_char = KIND('A')
#endif

#include "Flashx_mpi.h"

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t)     :: tileDesc

   integer(io_fileID_t), intent(in) :: fileID

   type(tree_data_t) :: tree_data
   character(len=MAX_STRING_LENGTH) :: buff
   character(len=16) :: numToStr

   integer :: localNumBlocks, blockID
   integer :: i, lb, j, u, AllocateStatus
   integer :: dowrite
   integer, allocatable :: procnumber(:)

   ! allocate storage to hold a single variable information
   ! this should only be a small memory overhead
   integer, parameter :: single = SELECTED_REAL_KIND(p=6)

!  real,save :: unkBufGC(1,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MAXBLOCKS)
   real, allocatable :: unkBufGC(:, :, :, :, :)

   ! allocate storage to hold the coordinate information and bounding box
   ! information
   real(kind=single) :: tmpSingle(MDIM, MAXBLOCKS)
   real(kind=single) :: bndSingle(2, MDIM, MAXBLOCKS)
   real(kind=single) :: spMax, spMin

   logical :: isPlotVar
   logical, allocatable :: isPlotVars(:)
   real, allocatable :: globalVarMin(:), globalVarMax(:)
   real, allocatable :: globalFaceXMin(:), globalFaceXMax(:)
   real, allocatable :: globalFaceYMin(:), globalFaceYMax(:)
   real, allocatable :: globalFaceZMin(:), globalFaceZMax(:)

   real, allocatable :: unkBuf(:, :, :, :)
   !real, allocatable :: unkBuf(:,:,:,:,:)
   real, allocatable :: faceXBuf(:, :, :, :, :)
   real, allocatable :: faceYBuf(:, :, :, :, :)
   real, allocatable :: faceZBuf(:, :, :, :, :)

   real(kind=single), allocatable :: unkt(:, :, :, :, :)

   logical :: writeGuardCells = .true.

   integer :: blkLimits(2, MDIM), blkLimitsGC(2, MDIM)
   integer :: blkLimitsFX(2, MDIM), blkLimitsGcFX(2, MDIM)
   integer :: blkLimitsFY(2, MDIM), blkLimitsGcFY(2, MDIM)
   integer :: blkLimitsFZ(2, MDIM), blkLimitsGcFZ(2, MDIM)

   !Adjust our offset so that we can eliminate redundant data in a group
   integer :: localOffset, localRank, splitOffset
   integer :: ierr

   !presentDims needs to depend on FILE_FORMAT_VERSION.
   integer, parameter :: xferType = IO_WRITE_XFER, presentDims = MDIM, &
                         libType = IO_FILE_HDF5
   integer :: numFileBlks

   integer, allocatable :: lrefine(:), levBlock(:)
   integer, allocatable :: nodetype(:)
   real, allocatable :: gid(:, :), coord(:, :), bsize(:, :), bnd_box(:, :, :)

   real, dimension(:, :, :, :), POINTER :: solnData
   real, dimension(:, :, :, :, :), POINTER :: scratch
   integer, dimension(MDIM) :: lo, hi

   nullify (solnData)
   nullify (scratch)

  !! call the generic function prepareLists to allocate and
  !! fill the runtime parameter lists and the scalar lists
   call io_prepareListsWrite()

   call Grid_getLocalNumBlks(localNumBlocks)

   allocate (tree_data%lrefine(localNumBlocks))
   allocate (tree_data%nodetype(localNumBlocks))
   allocate (tree_data%coord(MDIM, localNumBlocks))
   allocate (tree_data%bsize(MDIM, localNumBlocks))
   !allocate(tree_data % gid(MDIM, localNumBlocks))
   allocate (tree_data%bnd_box(2, MDIM, localNumBlocks))
   allocate (tree_data%gid(2*MDIM + 1 + 2**MDIM, localNumBlocks))

   allocate (tree_data%procnumber(max(1, localNumBlocks)))

   ! populate tree_data vars
   tree_data%procnumber(:) = io_meshMe
   tree_data%lrefine(1:localNumBlocks) = gr_ioBlkLrefine(1:localNumBlocks)
   tree_data%nodetype(1:localNumBlocks) = gr_ioBlkNodeType(1:localNumBlocks)
   tree_data%bnd_box(:, :, 1:localNumBlocks) = gr_ioBlkBoundBox(:, :, 1:localNumBlocks)
   tree_data%coord(:, 1:localNumBlocks) = gr_ioBlkCoords(:, 1:localNumBlocks)
   tree_data%bsize(:, 1:localNumBlocks) = gr_ioBlkBsize(:, 1:localNumBlocks)

   tree_data%gid(:, :) = -1

   !Find our local offset for split-file IO
   if (io_outputSplitNum > 1) then

      call MPI_ALLREDUCE(gr_globalOffset, splitOffset, 1, FLASH_INTEGER, &
                         MPI_MIN, io_comm, ierr)
      localOffset = gr_globalOffset - splitOffset
      !find number of blocks for a file
      call MPI_ALLREDUCE(localNumBlocks, io_splitNumBlks, 1, FLASH_INTEGER, &
                         MPI_SUM, io_comm, ierr)
   else
      localOffset = gr_globalOffset
      io_splitNumBlks = gr_globalNumBlocks
   end if

   numFileBlks = io_splitNumBlks

   !If we are writing a checkpoint file (io_doublePrecision == .true.)
  !!write all of the unk vars.
   !If writing a plotfile then only certain variables are written
   if (io_doublePrecision) then

     !! write the header info
      call io_h5write_header(io_meshMe, ubound(io_unklabelsGlobal, 1), fileID, io_geometry, &
                             io_unklabelsGlobal, io_setupCall, io_fileCreationTime, io_flashRelease, &
                             io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
                             io_setupTimeStamp, io_buildTimeStamp, io_fileFormatVersion, &
                             io_outputSplitNum)

   else
      call io_h5write_header(io_meshMe, io_nPlotVars, fileID, io_geometry, &
                             io_plotVarStr, io_setupCall, io_fileCreationTime, io_flashRelease, &
                             io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
                             io_setupTimeStamp, io_buildTimeStamp, io_fileFormatVersion, &
                             io_outputSplitNum)

   end if

  !! write the runtime parameters
   call io_h5write_runtime_parameters(io_globalMe, &
                                      fileID, &
                                      io_numRealParms, &
                                      io_realParmNames, &
                                      io_realParmValues, &
                                      io_numIntParms, &
                                      io_intParmNames, &
                                      io_intParmValues, &
                                      io_numStrParms, &
                                      io_strParmNames, &
                                      io_strParmValues, &
                                      io_numLogParms, &
                                      io_logParmNames, &
                                      io_logToIntParmValues, &
                                      io_outputSplitNum)

  !! write the scalars
   call io_h5write_scalars(io_globalMe, &
                           fileID, &
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
                           io_logToIntScalarValues, &
                           io_outputSplitNum)

   call io_finalizeListsWrite()
   call io_createDatasets(fileID, numFileBlks, presentDims)

   ! write tree data to hdf5 file
   call io_xfer_tree_data(tree_data, fileID, IO_FILE_HDF5, IO_WRITE_XFER, &
                          localNumBlocks, localOffset, presentDims, numFileBlks)

   deallocate (tree_data%procnumber)
   nullify (tree_data%procnumber)
   deallocate (tree_data%coord)
   nullify (tree_data%coord)
   deallocate (tree_data%bsize)
   nullify (tree_data%bsize)
   deallocate (tree_data%bnd_box)
   nullify (tree_data%bnd_box)
   deallocate (tree_data%lrefine)
   nullify (tree_data%lrefine)
   deallocate (tree_data%nodetype)
   nullify (tree_data%nodetype)
   deallocate (tree_data%gid)
   nullify (tree_data%gid)

   allocate (globalVarMin(ubound(io_unklabelsGlobal, 1)))
   allocate (globalVarMax(ubound(io_unklabelsGlobal, 1)))

   !get the max and minimum variables
   call io_getVarExtrema(ubound(io_unklabelsGlobal, 1), globalVarMin, globalVarMax, CENTER)

   if (io_doublePrecision .and. io_chkGuardCellsOutput) then ! blkLimits is not used otherwise
      blockID = 1
      call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
#ifndef FL_NON_PERMANENT_GUARDCELLS
      blkLimits = blkLimitsGC
#endif
   end if

   !do i = UNK_VARS_BEGIN,UNK_VARS_END
   do u = 1, ubound(io_unklabelsGlobal, 1)
      call Simulation_mapStrToInt(io_unklabelsGlobal(u), i, MAPBLOCK_UNK)
      ! only write mesh replicated data from mesh 0
      dowrite = 0
      if (i /= NONEXISTENT) then
         if (localNumBlocks > 0 .and. (io_acrossMe .eq. 0 .or. io_unkNonRep(i) > 0)) dowrite = 1
      end if

      ! populate and unkBuf
      if (localNumBlocks > 0) then
         allocate (unkBuf(NXB, NYB, NZB, localNumBlocks))
         unkBuf = 0.0
         call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.); lb = 1
         do while (itor%isValid())
            call itor%currentTile(tileDesc)
            call tileDesc%getDataPtr(solnData, CENTER, localFlag=.TRUE.)
            unkBuf(:, :, :, lb) = solnData(io_ilo:io_ihi, io_jlo:io_jhi, &
                                           io_klo:io_khi, i)
            call tileDesc%releaseDataPtr(solnData, CENTER)
            call itor%next(); lb = lb + 1
         end do
         call Grid_releaseTileIterator(itor)
         nullify (solnData)

         if (io_doublePrecision) then
            if (io_chkGuardCellsOutput) then
               call io_h5write_unknowns(io_globalMe, &
                                        fileID, &
                                        GRID_IHI_GC, &
                                        GRID_JHI_GC, &
                                        GRID_KHI_GC, &
                                        unkBuf, &
                                        globalVarMin(u), &
                                        globalVarMax(u), &
                                        io_unklabelsGlobal(u), &
                                        localNumBlocks, &
                                        io_splitNumBlks, &
                                        localOffset, &
                                        dowrite)
            else
               call io_h5write_unknowns(io_globalMe, &
                                        fileID, &
                                        NXB, &
                                        NYB, &
                                        NZB, &
                                        unkBuf, &
                                        globalVarMin(u), &
                                        globalVarMax(u), &
                                        io_unklabelsGlobal(u), &
                                        localNumBlocks, &
                                        io_splitNumBlks, &
                                        localOffset, &
                                        dowrite)
            end if
         else
            isPlotVar = any(io_unklabelsGlobal(u) == io_plotVarStr(1:io_nPlotVars))
            if (isPlotVar) then
               if (io_plotfileGridQuantityDP) then
                  call io_h5write_unknowns(io_globalMe, &
                                           fileID, &
                                           NXB, &
                                           NYB, &
                                           NZB, &
                                           unkBuf, &
                                           globalVarMin(u), &
                                           globalVarMax(u), &
                                           io_unklabelsGlobal(u), &
                                           localNumBlocks, &
                                           io_splitNumBlks, &
                                           localOffset, &
                                           dowrite)
               else
                  spMax = real(globalVarMax(u), kind=single)
                  spMin = real(globalVarMin(u), kind=single)
                  call io_h5write_unknowns_sp(io_globalMe, &
                                              fileID, &
                                              NXB, &
                                              NYB, &
                                              NZB, &
                                              spMin, &
                                              spMax, &
                                              unkBuf, &
                                              io_unklabelsGlobal(u), &
                                              localNumBlocks, &
                                              io_splitNumBlks, &
                                              localOffset, &
                                              dowrite)
               end if
            end if
         end if

         call MPI_BARRIER(io_globalComm, ierr)
         deallocate (unkBuf)

      end if ! endif localNumBlocks > 0
   end do

   deallocate (globalVarMin)
   deallocate (globalVarMax)

   allocate (unkBuf(NXB, NYB, NZB, 1:localNumBlocks))
   allocate (globalVarMin(NSCRATCH_GRID_VARS))
   allocate (globalVarMax(NSCRATCH_GRID_VARS))
   allocate (isPlotVars(NSCRATCH_GRID_VARS))

   do i = SCRATCH_GRID_VARS_BEGIN, SCRATCH_GRID_VARS_END
      call io_isPlotVar(i, isPlotVars(i), MAPBLOCK_SCRATCH)
   end do

   !get the max and minimum variables
   if (ANY(isPlotVars)) then
      call io_getVarExtrema(NSCRATCH_GRID_VARS, globalVarMin, globalVarMax, SCRATCH)
   end if

   dowrite = 0
   if (localNumBlocks > 0) dowrite = 1

   !write the scratch grid vars if the user defines any in flash.par
   !we can use the same routine as when writing the unknowns.
   do i = SCRATCH_GRID_VARS_BEGIN, SCRATCH_GRID_VARS_END
      isPlotVar = isPlotVars(i)
      if (isPlotVar) then

         if (io_doublePrecision .or. io_plotfileGridQuantityDP) then

            unkBuf(1:NXB, 1:NYB, 1:NZB, 1:localNumBlocks) = &
               scratch(io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi, i, 1:localNumBlocks)
            call io_h5write_unknowns(io_globalMe, &
                                     fileID, &
                                     NXB, &
                                     NYB, &
                                     NZB, &
                                     unkBuf, &
                                     globalVarMin(i), &
                                     globalVarMax(i), &
                                     io_plotGridVarStr(i), &
                                     localNumBlocks, &
                                     io_splitNumBlks, &
                                     localOffset, &
                                     dowrite)

         else

            unkBuf(1:NXB, 1:NYB, 1:NZB, 1:localNumBlocks) = &
               real(scratch(io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi, i, 1:localNumBlocks), kind=single)

            spMax = real(globalVarMax(i), kind=single)
            spMin = real(globalVarMin(i), kind=single)

            call io_h5write_unknowns_sp(io_globalMe, &
                                        fileID, &
                                        NXB, &
                                        NYB, &
                                        NZB, &
                                        spMin, &
                                        spMax, &
                                        unkBuf, &
                                        io_plotGridVarStr(i), &
                                        localNumBlocks, &
                                        io_splitNumBlks, &
                                        localOffset, &
                                        dowrite)

         end if
      end if
   end do

   deallocate (unkBuf)
   deallocate (isPlotVars)
   deallocate (globalVarMin)
   deallocate (globalVarMax)

#if(NFACE_VARS>0)
   if (io_doublePrecision) then

      allocate (globalFaceXMin(NFACE_VARS))
      allocate (globalFaceXMax(NFACE_VARS))
      call io_getVarExtrema(NFACE_VARS, globalFaceXMin, globalFaceXMax, FACEX)

#if NDIM >= 2
      allocate (globalFaceYMin(NFACE_VARS))
      allocate (globalFaceYMax(NFACE_VARS))
      call io_getVarExtrema(NFACE_VARS, globalFaceYMin, globalFaceYMax, FACEY)
#endif

#if NDIM == 3
      allocate (globalFaceZMin(NFACE_VARS))
      allocate (globalFaceZMax(NFACE_VARS))
      call io_getVarExtrema(NFACE_VARS, globalFaceZMin, globalFaceZMax, FACEZ)
#endif

      if (.NOT. io_chkGuardCellsOutput) then
         allocate (faceXBuf(1, NXB + 1, NYB, NZB, localNumBlocks))
         faceXBuf = 0.0
#if NDIM >= 2
         allocate (faceYBuf(1, NXB, NYB + 1, NZB, localNumBlocks))
         faceYBuf = 0.0
#endif
#if NDIM == 3
         allocate (faceZBuf(1, NXB, NYB, NZB + 1, localNumBlocks))
         faceZBuf = 0.0
#endif

         do i = 1, NFACE_VARS
            call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.); lb = 1
            do while (itor%isValid())
               call itor%currentTile(tileDesc)
               call tileDesc%getDataPtr(solnData, FACEX)
               lo = tileDesc%limits(LOW, :)
               hi = tileDesc%limits(HIGH, :)
               faceXBuf(1, 1:NXB + 1, 1:NYB, 1:NZB, lb) = solnData(lo(IAXIS):hi(IAXIS) + 1, &
                                                                   lo(JAXIS):hi(JAXIS), &
                                                                   lo(KAXIS):hi(KAXIS), i)
               call tileDesc%releaseDataPtr(solnData, FACEX)
               call itor%next(); lb = lb + 1
            end do
            call Grid_releaseTileIterator(itor)
            nullify (solnData)

            call io_h5write_unknowns(io_globalMe, &
                                     fileID, &
                                     NXB + 1, &
                                     NYB, &
                                     NZB, &
                                     faceXBuf, &
                                     globalFaceXMin(i), &
                                     globalFaceXMax(i), &
                                     io_faceXVarLabels(i), &
                                     localNumBlocks, &
                                     io_splitNumBlks, &
                                     localOffset, &
                                     dowrite)

#if NDIM >= 2
            call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.); lb = 1
            do while (itor%isValid())
               call itor%currentTile(tileDesc)
               call tileDesc%getDataPtr(solnData, FACEY)
               lo = tileDesc%limits(LOW, :)
               hi = tileDesc%limits(HIGH, :)
               faceYBuf(1, 1:NXB, 1:NYB + 1, 1:NZB, lb) = solnData(lo(IAXIS):hi(IAXIS), &
                                                                   lo(JAXIS):hi(JAXIS) + 1, &
                                                                   lo(KAXIS):hi(KAXIS), i)
               call tileDesc%releaseDataPtr(solnData, FACEY)
               call itor%next(); lb = lb + 1
            end do
            call Grid_releaseTileIterator(itor)
            nullify (solnData)

            call io_h5write_unknowns(io_globalMe, &
                                     fileID, &
                                     NXB, &
                                     NYB + 1, &
                                     NZB, &
                                     faceYBuf, &
                                     globalFaceYMin(i), &
                                     globalFaceYMax(i), &
                                     io_faceYVarLabels(i), &
                                     localNumBlocks, &
                                     io_splitNumBlks, &
                                     localOffset, &
                                     dowrite)
#endif

#if NDIM == 3
            call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.); lb = 1
            do while (itor%isValid())
               call itor%currentTile(tileDesc)
               call tileDesc%getDataPtr(solnData, FACEZ)
               lo = tileDesc%limits(LOW, :)
               hi = tileDesc%limits(HIGH, :)
               faceZBuf(1, 1:NXB, 1:NYB, 1:NZB + 1, lb) = solnData(lo(IAXIS):hi(IAXIS), &
                                                                   lo(JAXIS):hi(JAXIS), &
                                                                   lo(KAXIS):hi(KAXIS) + 1, i)
               call tileDesc%releaseDataPtr(solnData, FACEZ)
               call itor%next(); lb = lb + 1
            end do
            call Grid_releaseTileIterator(itor)
            nullify (solnData)

            call io_h5write_unknowns(io_globalMe, &
                                     fileID, &
                                     NXB, &
                                     NYB, &
                                     NZB + 1, &
                                     faceZBuf, &
                                     globalFaceZMin(i), &
                                     globalFaceZMax(i), &
                                     io_faceZVarLabels(i), &
                                     localNumBlocks, &
                                     io_splitNumBlks, &
                                     localOffset, &
                                     dowrite)
#endif
         end do

      else ! if (.NOT. io_chkGuardCellsOutput)

         allocate (faceXBuf(1, GRID_IHI_GC + 1, GRID_JHI_GC, GRID_KHI_GC, localNumBlocks))
         faceXBuf = 0.0
#if NDIM >= 2
         allocate (faceYBuf(1, GRID_IHI_GC, GRID_JHI_GC + 1, GRID_KHI_GC, localNumBlocks))
         faceYBuf = 0.0
#endif
#if NDIM == 3
         allocate (faceZBuf(1, GRID_IHI_GC, GRID_JHI_GC, GRID_KHI_GC + 1, localNumBlocks))
         faceZBuf = 0.0
#endif
         do i = 1, NFACE_VARS

            call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.); lb = 1
            do while (itor%isValid())
               call itor%currentTile(tileDesc)
               call tileDesc%getDataPtr(solnData, FACEX)
               faceXBuf(1, :, :, :, lb) = solnData(:, :, :, i)
               call tileDesc%releaseDataPtr(solnData, FACEX)
               call itor%next(); lb = lb + 1
            end do
            call Grid_releaseTileIterator(itor)
            nullify (solnData)

            call io_h5write_unknowns(io_globalMe, &
                                     fileID, &
                                     GRID_IHI_GC + 1, &
                                     GRID_JHI_GC, &
                                     GRID_KHI_GC, &
                                     faceXBuf, &
                                     globalFaceXMin(i), &
                                     globalFaceXMax(i), &
                                     io_faceXVarLabels(i), &
                                     localNumBlocks, &
                                     io_splitNumBlks, &
                                     localOffset, &
                                     dowrite)
#if NDIM >= 2
            call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.); lb = 1
            do while (itor%isValid())
               call itor%currentTile(tileDesc)
               call tileDesc%getDataPtr(solnData, FACEY)
               faceYBuf(1, :, :, :, lb) = solnData(:, :, :, i)
               call tileDesc%releaseDataPtr(solnData, FACEY)
               call itor%next(); lb = lb + 1
            end do
            call Grid_releaseTileIterator(itor)
            nullify (solnData)

            call io_h5write_unknowns(io_globalMe, &
                                     fileID, &
                                     GRID_IHI_GC, &
                                     GRID_JHI_GC + 1, &
                                     GRID_KHI_GC, &
                                     faceYBuf, &
                                     globalFaceYMin(i), &
                                     globalFaceYMax(i), &
                                     io_faceYVarLabels(i), &
                                     localNumBlocks, &
                                     io_splitNumBlks, &
                                     localOffset, &
                                     dowrite)
#endif

#if NDIM == 3
            call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.); lb = 1
            do while (itor%isValid())
               call itor%currentTile(tileDesc)
               call tileDesc%getDataPtr(solnData, FACEZ)
               faceZBuf(1, :, :, :, lb) = solnData(:, :, :, i)
               call tileDesc%releaseDataPtr(solnData, FACEZ)
               call itor%next(); lb = lb + 1
            end do
            call Grid_releaseTileIterator(itor)
            nullify (solnData)

            call io_h5write_unknowns(io_globalMe, &
                                     fileID, &
                                     GRID_IHI_GC, &
                                     GRID_JHI_GC, &
                                     GRID_KHI_GC + 1, &
                                     faceZBuf, &
                                     globalFaceZMin(i), &
                                     globalFaceZMax(i), &
                                     io_faceZVarLabels(i), &
                                     localNumBlocks, &
                                     io_splitNumBlks, &
                                     localOffset, &
                                     dowrite)
#endif
         end do

      end if                   !if (.NOT. io_chkGuardCellsOutput)

      deallocate (faceXBuf)
      deallocate (globalFaceXMin)
      deallocate (globalFaceXMax)

#if NDIM >= 2
      deallocate (faceYBuf)
      deallocate (globalFaceYMin)
      deallocate (globalFaceYMax)
#endif

#if NDIM == 3
      deallocate (faceZBuf)
      deallocate (globalFaceZMin)
      deallocate (globalFaceZMax)
#endif
   end if
#endif

   if (io_globalMe .EQ. MASTER_PE) then
      write (numToStr, "(I7)") gr_globalNumBlocks
      buff = "wrote "//numToStr//" blocks"
      call Logfile_stamp(buff, '[io_writeData]')
   end if

   return

end subroutine io_writeData
