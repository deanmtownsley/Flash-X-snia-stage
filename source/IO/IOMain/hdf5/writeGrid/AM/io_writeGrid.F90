!!****if* source/IO/IOMain/hdf5/writeGrid/AM/io_writeGrid
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
!!
!!
!!
!! NAME
!!
!!  io_writeGrid
!!
!!
!! SYNOPSIS
!!
!!  call io_writeGrid()
!!
!!
!! DESCRIPTION
!!
!! This function writes the grid information to an hdf5 file to store the
!! Paramesh or AMReX Grid cell coordinates (Left, Center, Right)
!! and the cell metrics for later use in post-processing
!!
!! Currently only supports hdf5 IO
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine io_writeGrid()

#include "constants.h"
#include "Simulation.h"
#include "io_flash.h"

   use Driver_interface, ONLY: Driver_abort
   use IO_data, ONLY: io_plotFileNumber, io_meshMe, io_meshComm

   use Grid_data, ONLY: gr_globalNumBlocks, gr_globalNumProcs

   use Grid_interface, ONLY: Grid_getNumBlksFromType, Grid_getTileIterator, &
                             Grid_releaseTileIterator, Grid_getCellCoords

   use Grid_tile, ONLY: Grid_tile_t
   use Grid_iterator, ONLY: Grid_iterator_t

   use Timers_interface, ONLY: Timers_start, Timers_stop

   use HDF5

#include "Flashx_mpi_implicitNone.fh"

   integer :: ierr, fc, lb, a, b, c, d, jproc, offset, localNumBlocks
   integer :: status(MPI_STATUS_SIZE)
   character(len=5) :: CrdLbs(MDIM, 4), MtrLbs(MDIM, 3)
   real, dimension(MDIM, 4) :: CrdMax, CrdMin, MtrMax, MtrMin
   real, allocatable, dimension(:) :: iBuff, jBuff, kBuff, dBuff
   real, allocatable, dimension(:, :) :: fCoord
   real, allocatable, dimension(:, :, :, :) :: Coords
   real, allocatable, dimension(:, :) :: Deltas
   integer :: lnblocks
   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor
   integer, dimension(MDIM) :: lo, hi
   real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
   real :: del(MDIM)

   ! locals necessary to read hdf5 file
   integer :: error
   integer(HID_T) :: file_id, dspc_id, dset_id, aspc_id, attr_id
   integer(HSIZE_T), dimension(2) :: dset_dims
   integer(HSIZE_T), dimension(1) :: dset_sngl
   character(len=32) :: dsetname, attrname
   character(len=MAX_STRING_LENGTH) :: filename

   call Grid_getNumBlksFromType(ALL_BLKS, lnblocks)

#ifdef FLASH_IO_HDF5

   call Timers_start("io_writeGrid")

   ! Open an hdf5 file for writing grid information
   if (io_meshMe == MASTER_PE) then

      call io_getOutputName(io_plotFileNumber, "hdf5_", "grd_", filename, .false.)

      call h5open_f(error)

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      if (file_id == -1) then
         print *, "Error: Unable to initialize GRID OUTPUT file"
         call Driver_abort("Unable to initialize hdf5 file")
      end if

   end if

   if (io_meshMe == MASTER_PE) then

      ! Create dataset labels
      CrdLbs(IAXIS, LEFT_EDGE:RIGHT_EDGE + 1) = (/'xxxl', 'xxxc', 'xxxr', 'xxxf'/)
      CrdLbs(JAXIS, LEFT_EDGE:RIGHT_EDGE + 1) = (/'yyyl', 'yyyc', 'yyyr', 'yyyf'/)
      CrdLbs(KAXIS, LEFT_EDGE:RIGHT_EDGE + 1) = (/'zzzl', 'zzzc', 'zzzr', 'zzzf'/)
      MtrLbs(IAXIS, LEFT_EDGE:RIGHT_EDGE) = (/'  ', 'dx', '  '/)
      MtrLbs(JAXIS, LEFT_EDGE:RIGHT_EDGE) = (/'  ', 'dy', '  '/)
      MtrLbs(KAXIS, LEFT_EDGE:RIGHT_EDGE) = (/'  ', 'dz', '  '/)

      allocate (Coords(max(NXB, NYB, NZB), gr_globalNumBlocks, 3, MDIM))
      allocate (Deltas(gr_globalNumBlocks, MDIM))

   end if

   offset = 0
   do jproc = 0, gr_globalNumProcs - 1

      if (io_meshMe == MASTER_PE .and. jproc == MASTER_PE) then

         if (lnblocks > 0) then
            allocate (iBuff(NXB*lnblocks*3), jBuff(NYB*lnblocks*3), kBuff(NZB*lnblocks*3), dBuff(lnblocks*MDIM))
            do fc = 1, 3
               lb = 1
               call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
               do while (itor%isValid())
                  call itor%currentTile(tileDesc)

                  lo = tileDesc%limits(LOW, :)
                  hi = tileDesc%limits(HIGH, :)

                  allocate (xCoord(lo(IAXIS):hi(IAXIS)))
                  allocate (yCoord(lo(JAXIS):hi(JAXIS)))
                  allocate (zCoord(lo(KAXIS):hi(KAXIS)))

                  xCoord = 0.0
                  yCoord = 0.0
                  zCoord = 0.0

                  call Grid_getCellCoords(IAXIS, fc, tileDesc%level, lo, hi, xCoord)
                  call Grid_getCellCoords(JAXIS, fc, tileDesc%level, lo, hi, yCoord)
                  call Grid_getCellCoords(KAXIS, fc, tileDesc%level, lo, hi, zCoord)

                  call tileDesc%deltas(del)

                  iBuff(NXB*lnblocks*(fc - 1) + NXB*(lb - 1) + 1:NXB*lnblocks*(fc - 1) + NXB*lb) = &
                     xCoord(lo(IAXIS):hi(IAXIS))

                  jBuff(NYB*lnblocks*(fc - 1) + NYB*(lb - 1) + 1:NYB*lnblocks*(fc - 1) + NYB*lb) = &
                     yCoord(lo(JAXIS):hi(JAXIS))

                  kBuff(NZB*lnblocks*(fc - 1) + NZB*(lb - 1) + 1:NZB*lnblocks*(fc - 1) + NZB*lb) = &
                     zCoord(lo(KAXIS):hi(KAXIS))

                  dBuff(lnblocks*(fc - 1) + lb) = del(fc)

                  deallocate (xCoord, yCoord, zCoord)

                  lb = lb + 1
                  call itor%next()
               end do
               call Grid_releaseTileIterator(itor)
            end do
            Coords(1:NXB, offset + 1:offset + lnblocks, 1:3, IAXIS) = reshape(iBuff, (/NXB, lnblocks, 3/))
            Coords(1:NYB, offset + 1:offset + lnblocks, 1:3, JAXIS) = reshape(jBuff, (/NYB, lnblocks, 3/))
            Coords(1:NZB, offset + 1:offset + lnblocks, 1:3, KAXIS) = reshape(kBuff, (/NZB, lnblocks, 3/))
            Deltas(offset + 1:offset + lnblocks, 1:MDIM) = reshape(dBuff, (/lnblocks, MDIM/))
            deallocate (iBuff, jBuff, kBuff, dBuff)
            offset = offset + lnblocks
         end if

      end if

      if (io_meshMe == MASTER_PE .and. jproc /= MASTER_PE) then

         call MPI_Recv(localNumBlocks, 1, FLASH_INTEGER, jproc, 1, io_meshComm, status, ierr)
         if (localNumBlocks > 0) then
            allocate (iBuff(NXB*localNumBlocks*3), &
                      jBuff(NYB*localNumBlocks*3), &
                      kBuff(NZB*localNumBlocks*3), &
                      dBuff(localNumBlocks*MDIM))

            call MPI_Recv(iBuff, NXB*localNumBlocks*3, FLASH_REAL, jproc, 1 + IAXIS, io_meshComm, status, ierr)
            call MPI_Recv(jBuff, NYB*localNumBlocks*3, FLASH_REAL, jproc, 1 + JAXIS, io_meshComm, status, ierr)
            call MPI_Recv(kBuff, NZB*localNumBlocks*3, FLASH_REAL, jproc, 1 + KAXIS, io_meshComm, status, ierr)
            call MPI_Recv(dBuff, localNumBlocks*3, FLASH_REAL, jproc, 2 + MDIM, io_meshComm, status, ierr)

            Coords(1:NXB, offset + 1:offset + localNumBlocks, 1:3, IAXIS) = &
               reshape(iBuff, (/NXB, localNumBlocks, 3/))

            Coords(1:NYB, offset + 1:offset + localNumBlocks, 1:3, JAXIS) = &
               reshape(jBuff, (/NYB, localNumBlocks, 3/))

            Coords(1:NZB, offset + 1:offset + localNumBlocks, 1:3, KAXIS) = &
               reshape(kBuff, (/NZB, localNumBlocks, 3/))

            Deltas(offset + 1:offset + localNumBlocks, 1:MDIM) = reshape(dBuff, (/localNumBlocks, MDIM/))

            deallocate (iBuff, jBuff, kBuff, dBuff)

            offset = offset + localNumBlocks
         end if
      end if

      if (io_meshMe /= MASTER_PE .and. jproc == io_meshMe) then

         call MPI_Send(lnblocks, 1, FLASH_INTEGER, MASTER_PE, 1, io_meshComm, ierr)
         if (lnblocks > 0) then
            allocate (iBuff(NXB*lnblocks*3), jBuff(NYB*lnblocks*3), kBuff(NZB*lnblocks*3), dBuff(lnblocks*MDIM))
            do fc = 1, 3
               lb = 1
               call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
               do while (itor%isValid())
                  call itor%currentTile(tileDesc)

                  lo = tileDesc%limits(LOW, :)
                  hi = tileDesc%limits(HIGH, :)

                  allocate (xCoord(lo(IAXIS):hi(IAXIS)))
                  allocate (yCoord(lo(JAXIS):hi(JAXIS)))
                  allocate (zCoord(lo(KAXIS):hi(KAXIS)))

                  xCoord = 0.0
                  yCoord = 0.0
                  zCoord = 0.0

                  call Grid_getCellCoords(IAXIS, fc, tileDesc%level, lo, hi, xCoord)
                  call Grid_getCellCoords(JAXIS, fc, tileDesc%level, lo, hi, yCoord)
                  call Grid_getCellCoords(KAXIS, fc, tileDesc%level, lo, hi, zCoord)

                  call tileDesc%deltas(del)

                  iBuff(NXB*lnblocks*(fc - 1) + NXB*(lb - 1) + 1:NXB*lnblocks*(fc - 1) + NXB*lb) = &
                     xCoord(lo(IAXIS):hi(IAXIS))

                  jBuff(NYB*lnblocks*(fc - 1) + NYB*(lb - 1) + 1:NYB*lnblocks*(fc - 1) + NYB*lb) = &
                     yCoord(lo(JAXIS):hi(JAXIS))

                  kBuff(NZB*lnblocks*(fc - 1) + NZB*(lb - 1) + 1:NZB*lnblocks*(fc - 1) + NZB*lb) = &
                     zCoord(lo(KAXIS):hi(KAXIS))

                  dBuff(lnblocks*(fc - 1) + lb) = del(fc)

                  deallocate (xCoord, yCoord, zCoord)

                  lb = lb + 1
                  call itor%next()
               end do
               call Grid_releaseTileIterator(itor)
            end do
            call MPI_Send(iBuff, NXB*lnblocks*3, FLASH_REAL, MASTER_PE, 1 + IAXIS, io_meshComm, ierr)
            call MPI_Send(jBuff, NYB*lnblocks*3, FLASH_REAL, MASTER_PE, 1 + JAXIS, io_meshComm, ierr)
            call MPI_Send(kBuff, NZB*lnblocks*3, FLASH_REAL, MASTER_PE, 1 + KAXIS, io_meshComm, ierr)
            call MPI_Send(dBuff, lnblocks*3, FLASH_REAL, MASTER_PE, 2 + MDIM, io_meshComm, ierr)
            deallocate (iBuff, jBuff, kBuff, dBuff)
         end if

      end if

   end do

   if (io_meshMe == MASTER_PE) then

      ! Write coordinates to file
      do a = IAXIS, KAXIS

         ! Determine bounds
         select case (a)
         case (IAXIS)
            d = NXB
         case (JAXIS)
            d = NYB
         case (KAXIS)
            d = NZB
         end select

         ! write each face per axis
         do b = LEFT_EDGE, RIGHT_EDGE

            ! find extreme values
            CrdMax(a, b) = maxval(Coords(1:d, :, b, a))
            CrdMin(a, b) = minval(Coords(1:d, :, b, a))
            MtrMax(a, b) = maxval(Deltas(:, a))
            MtrMin(a, b) = minval(Deltas(:, a))

            ! write dimensions
            dsetname = CrdLbs(a, b)
            dset_dims = (/d, gr_globalNumBlocks/)
            call h5screate_simple_f(2, dset_dims, dspc_id, error)
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Coords(1:d, :, b, a), dset_dims, error)

            attrname = "maximum"
            dset_sngl = (/1/)
            call h5screate_simple_f(1, dset_sngl, aspc_id, error)
            call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error)
            call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, CrdMax(a, b), dset_sngl, error)
            call h5aclose_f(attr_id, error)

            attrname = "minimum"
            call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error)
            call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, CrdMin(a, b), dset_sngl, error)
            call h5aclose_f(attr_id, error)
            call h5sclose_f(aspc_id, error)

            call h5dclose_f(dset_id, error)
            call h5sclose_f(dspc_id, error)

            ! write metrics
            if (b == CENTER) then
               dsetname = MtrLbs(a, b)
               dset_sngl = (/gr_globalNumBlocks/)
               call h5screate_simple_f(1, dset_sngl, dspc_id, error)
               call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)
               call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Deltas(:, a), dset_dims, error)

               attrname = "maximum"
               dset_sngl = (/1/)
               call h5screate_simple_f(1, dset_sngl, aspc_id, error)
               call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error)
               call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, MtrMax(a, b), dset_sngl, error)
               call h5aclose_f(attr_id, error)

               attrname = "minimum"
               call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error)
               call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, MtrMin(a, b), dset_sngl, error)
               call h5aclose_f(attr_id, error)
               call h5sclose_f(aspc_id, error)

               call h5dclose_f(dset_id, error)
               call h5sclose_f(dspc_id, error)
            end if

         end do

         ! Consolodate axis "face" points
         allocate (fCoord(d + 1, gr_globalNumBlocks))
         fCoord(1, :) = Coords(1, :, LEFT_EDGE, a)
         fCoord(2:d + 1, :) = Coords(1:d, :, RIGHT_EDGE, a)
         if (a == KAXIS .AND. d == 1) fCoord(2, :) = fCoord(1, :) + 0.000001

         ! find extreme values
         CrdMax(a, 4) = maxval(fCoord)
         CrdMin(a, 4) = minval(fCoord)

         ! write dimensions
         dsetname = CrdLbs(a, 4)
         dset_dims = (/d + 1, gr_globalNumBlocks/)
         call h5screate_simple_f(2, dset_dims, dspc_id, error)
         call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fCoord, dset_dims, error)

         attrname = "maximum"
         dset_sngl = (/1/)
         call h5screate_simple_f(1, dset_sngl, aspc_id, error)
         call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error)
         call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, CrdMax(a, b), dset_sngl, error)
         call h5aclose_f(attr_id, error)

         attrname = "minimum"
         call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error)
         call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, CrdMin(a, b), dset_sngl, error)
         call h5aclose_f(attr_id, error)
         call h5sclose_f(aspc_id, error)

         call h5dclose_f(dset_id, error)
         call h5sclose_f(dspc_id, error)

         deallocate (fCoord)

      end do

      ! release storage arrays
      deallocate (Coords)
      deallocate (Deltas)
   end if

   ! Close the hdf5 file
   if (io_meshMe == MASTER_PE) then

      print *, "*** Wrote grid data to ", trim(filename), " ****"

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end if

   call Timers_stop("io_writeGrid")

#endif

   return

end subroutine io_writeGrid
