!!****if* source/physics/Eos/unitTest/Hybrid/Eos_unitTest
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
!! NAME
!!
!!  Eos_unitTest
!!
!! SYNOPSIS
!!
!!  call Eos_unitTest(integer(IN) :: fileUnit,
!!                    logical(OUT) :: perfect
!!
!! DESCRIPTION
!!
!! This function is the unit test for the Eos unit. It is invoked in
!! the setup unitTest/Eos. The Config file for Eos unit test setup
!! requests a few extra variables in the main grid data structure for
!! Grid scope temporary storage. The Simulation_initBlock of the Eos
!! unit test initializes density in the right place for the DENS_VAR
!! variable (see Simulation.h for DENS_VAR, TEMP_VAR etc definitions), and
!! temperature and pressure in the extra storage space CTMP_VAR
!! and CPRS_VAR. The physical quantities at this point are not in
!! thermal equilibrium.
!!
!! The Eos_unit test starts by copying the initialized
!! temperature into the TEMP_VAR location and calling the
!! Eos_everywhere function with eosMode = MODE_DENS_TEMP, where
!! density and temperature are given and pressure and energy are
!! calculated. Now PRES_VAR and EINT_VAR contain values of pressure
!! and internal energy that are in thermal equilibrium, and the pressure values
!! are not necessarily what was stored in the extra storage space
!! during intialization.
!!
!! At this point in time three quantities; temperature,
!! pressure and energy are saved in the extra storage requested by
!! the unitTest/Eos setup, say OTMP_VAR, OPRS_VAR and OENT_VAR. Now
!! the Eos_unitTest function calls Eos_everywhere with eosMode =
!! MODE_DENS_PRES, followed by eosMode= MODE_DENS_EI.  If the
!! newly calculated values of temperature, pressure and energy are
!! the same as those saved in OTMP_VAR, OPRS_VAR and OENT_VAR, then
!! we can conclude that the Eos is working in MODE_DENS_PRES and
!! MODE_DENS_EI modes. However, we still can't say anything about the
!! MODE_DENS_TEMP mode. So we repeat the process by copying CPRS_VAR
!! into PRES_VAR and calling Eos_everywhere with MODE_DENS_PRES. We
!! again save the calculated values in the extra storage and make two
!! more Eos_everywhere calls with the remaining two modes. This time if
!! the new and old values of variables compare, we can conclude that
!! MODE_DENS_TEMP works too, and hence the unit test is successful.
!!
!!
!!  ARGUMENTS
!!
!!
!!   fileUnit : unit number for file opened by the unitTest/Eos setup
!!              in which to write results of the test
!!
!!   perfect : indicates test ran without error is true.
!!
!!***

!!REORDER(4): solnData

subroutine Eos_unitTest(fileUnit, perfect)

   use Eos_interface, ONLY: Eos_everywhere
   use Grid_interface, ONLY: Grid_getTileIterator, &
                             Grid_releaseTileIterator, &
                             Grid_getBlkType
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t
   use IO_interface, ONLY: IO_writeCheckpoint
   use Eos_data, ONLY: eos_meshMe, eos_meshNumProcs
   use eos_testData, ONLY: tolerance => eos_testTolerance

   use Timers_interface, only: Timers_start, Timers_stop, Timers_getSummary

   implicit none

# include "Eos.h"
# include "constants.h"
# include "Simulation.h"

   integer, intent(in) :: fileUnit
   logical, intent(out) :: perfect
   integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
   type(Grid_iterator_t) :: itor
   type(Grid_tile_t)     :: tileDesc

   integer :: blockID
   real, pointer, dimension(:, :, :, :):: solnData
   logical:: test1allB, test2allB, test3allB, test4allB !for all blocks

   integer, parameter :: maxPrintPE = 20
   integer, save :: nodeType = LEAF
   integer :: ib, ie, jb, je, kb, ke

   character(len=9), dimension(MODE_DENS_TEMP:MODE_DENS_PRES), parameter :: &
      mode_names = ["DENS_TEMP", "DENS_EINT", "DENS_PRES"]

   character(len=:), allocatable :: block_fmt, mode_fmt, set_fmt, err_fmt, &
                                    status_fmt, setup_fmt, init_fmt, result_fmt, &
                                    comp_fmt

   block_fmt  = "(i6,':  Block: ',i8,'  Type: ',i4)" !&

   setup_fmt  = "(i6,':  Setting up test for mode: ',a9,'    In: ',a4,',',a4,'    Out: ',a4,',',a4)" !&
   mode_fmt   = "(i6,':  Running test for mode:    ',a9,'    In: ',a4,',',a4,'    Out: ',a4,',',a4)" !&

   set_fmt    = "(i6,':  Setting ',a4,' = ',a4)"              !&
   init_fmt   = "(i6,':  Initialized ',a4,1x,a3': ',es24.15)" !&
   result_fmt = "(i6,':  Resulting ',a4,1x,a3': ',es24.15)"   !&

   status_fmt = "(i6,':  Result for test for ',a9,':  ',a4)"   !&
   err_fmt    = "(i6,':  Maximum error in ',a4,': ',es24.15)"  !&
   comp_fmt   = "(i6,':     ',a4,2x,es24.15,4x,a4,2x,es24.15)" !&

   nullify (solnData)

   call setupTest(MODE_DENS_TEMP, DENS_VAR, TEMP_VAR, DENS_VAR, CTMP_VAR, EINT_VAR, PRES_VAR, OTMP_VAR, OENT_VAR, OPRS_VAR)
   call IO_writeCheckpoint()   !! This is checkpoint 001

   call runTest(test1allB, MODE_DENS_EI, DENS_VAR, EINT_VAR, TEMP_VAR, PRES_VAR, OTMP_VAR, OPRS_VAR)
   call IO_writeCheckpoint()   !! This is checkpoint 002

   call runTest(test2allB, MODE_DENS_PRES, DENS_VAR, PRES_VAR, EINT_VAR, TEMP_VAR, OENT_VAR, OTMP_VAR)
   call IO_writeCheckpoint()   !! This is checkpoint 003

   call setupTest(MODE_DENS_PRES, DENS_VAR, PRES_VAR, DENS_VAR, PRES_VAR, EINT_VAR, TEMP_VAR, OPRS_VAR, OENT_VAR, OTMP_VAR)
   call IO_writeCheckpoint()   !! This is checkpoint 004

   call runTest(test3allB, MODE_DENS_EI, DENS_VAR, EINT_VAR, TEMP_VAR, PRES_VAR, OTMP_VAR, OPRS_VAR)
   call IO_writeCheckpoint()   !! This is checkpoint 005

   call runTest(test4allB, MODE_DENS_TEMP, DENS_VAR, TEMP_VAR, EINT_VAR, PRES_VAR, OENT_VAR, OPRS_VAR)
   call IO_writeCheckpoint()   !! This is checkpoint 006

   perfect = test1allB .and. test2allB .and. test3allB .and. test4allB
   if (eos_meshMe < maxPrintPE) then
      if (perfect) then
         print "(i6,':  ',a)", eos_meshMe, 'SUCCESS all tests passed'
      else
         print "(i6,':  ',a)", eos_meshMe, 'FAILURE some tests failed'
      end if
   end if

   call Timers_getSummary(1)

contains

   function getBlockID()
      implicit none
      integer :: getBlockID
#ifdef FLASH_GRID_AMREX
      getBlockID = tileDesc%grid_index
#else
      getBlockID = tileDesc%id
#endif
   end function getblockID

   subroutine setupTest(mode, in_1, in_2, in_1_init, in_2_init, out_1, out_2, in_2_exp, out_1_exp, out_2_exp)
      use Simulation_interface, only: Simulation_mapIntToStr

      implicit none

      integer, intent(in) :: mode
      integer, intent(in) :: in_1, in_2
      integer, intent(in) :: in_1_init, in_2_init
      integer, intent(in) :: out_1, out_2
      integer, intent(in) :: in_2_exp, out_1_exp, out_2_exp

      character(len=9) :: name
      character(len=4) :: in_1_name, in_2_name
      character(len=4) :: out_1_name, out_2_name
      character(len=4) :: out_1_exp_name, out_2_exp_name
      character(len=4) :: in_1_init_name, in_2_init_name

      name = mode_names(mode)

      call Simulation_mapIntToStr(in_1, in_1_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(in_2, in_2_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(in_1_init, in_1_init_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(in_2_init, in_2_init_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_1, out_1_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_2, out_2_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_1_exp, out_1_exp_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_2_exp, out_2_exp_name, MAPBLOCK_UNK)

      if (eos_meshMe < maxPrintPE) then
         write (*, *) repeat("=", 88)
         write (*, setup_fmt) eos_meshMe, name, in_1_name, in_2_name, out_1_name, out_2_name
         write (*, *) repeat("=", 88)
      end if

      call Grid_getTileIterator(itor, LEAF, tiling=.false.)
      TileLoopInit: do
         if (.not. itor%isValid()) exit TileLoopInit

         call itor%currentTile(tileDesc)

         blockID = getBlockID()
         call Grid_getBlkType(blockId, nodeType)

         blkLimits = tileDesc%limits

         call tileDesc%getDataPtr(solnData, CENTER)

         ib = blkLimits(LOW, IAXIS); ie = blkLimits(HIGH, IAXIS)
         jb = blkLimits(LOW, JAXIS); je = blkLimits(HIGH, JAXIS)
         kb = blkLimits(LOW, KAXIS); ke = blkLimits(HIGH, KAXIS)

         if (in_1 .ne. in_1_init) &
            solnData(in_1, ib:ie, jb:je, kb:ke) = solnData(in_1_init, ib:ie, jb:je, kb:ke)
         if (in_2 .ne. in_2_init) &
            solnData(in_2, ib:ie, jb:je, kb:ke) = solnData(in_2_init, ib:ie, jb:je, kb:ke)

         if (eos_meshMe < maxPrintPE) then
            write (*, block_fmt) eos_meshMe, blockID, nodeType
            write (*, init_fmt) eos_meshMe, in_1_name, "min", minval(solnData(in_1, ib:ie, jb:je, kb:ke))
            write (*, init_fmt) eos_meshMe, in_1_name, "max", maxval(solnData(in_1, ib:ie, jb:je, kb:ke))
            write (*, init_fmt) eos_meshMe, in_2_name, "min", minval(solnData(in_2, ib:ie, jb:je, kb:ke))
            write (*, init_fmt) eos_meshMe, in_2_name, "max", maxval(solnData(in_2, ib:ie, jb:je, kb:ke))
         end if

         call tileDesc%releaseDataPtr(solnData, CENTER)

         call itor%next()
      end do TileLoopInit

      call Grid_releaseTileIterator(itor)

      call Timers_start(name)
      call Eos_everywhere(mode)
      call Timers_stop(name)

      write (*, *) repeat("-", 88)

      call Grid_getTileIterator(itor, LEAF, tiling=.false.)
      TileLoopResult: do
         if (.not. itor%isValid()) exit TileLoopResult

         call itor%currentTile(tileDesc)

         blockID = getBlockID()
         call Grid_getBlkType(blockId, nodeType)

         blkLimits = tileDesc%limits

         call tileDesc%getDataPtr(solnData, CENTER)

         ib = blkLimits(LOW, IAXIS); ie = blkLimits(HIGH, IAXIS)
         jb = blkLimits(LOW, JAXIS); je = blkLimits(HIGH, JAXIS)
         kb = blkLimits(LOW, KAXIS); ke = blkLimits(HIGH, KAXIS)

         solnData(in_2_exp, ib:ie, jb:je, kb:ke) = solnData(in_2, ib:ie, jb:je, kb:ke)
         solnData(out_1_exp, ib:ie, jb:je, kb:ke) = solnData(out_1, ib:ie, jb:je, kb:ke)
         solnData(out_2_exp, ib:ie, jb:je, kb:ke) = solnData(out_2, ib:ie, jb:je, kb:ke)

         if (eos_meshMe < maxPrintPE) then
            write (*, block_fmt) eos_meshMe, blockID, nodeType
            write (*, result_fmt) eos_meshMe, out_1_name, "min", minval(solnData(out_1, ib:ie, jb:je, kb:ke))
            write (*, result_fmt) eos_meshMe, out_1_name, "max", maxval(solnData(out_1, ib:ie, jb:je, kb:ke))
            write (*, result_fmt) eos_meshMe, out_2_name, "min", minval(solnData(out_2, ib:ie, jb:je, kb:ke))
            write (*, result_fmt) eos_meshMe, out_2_name, "max", maxval(solnData(out_2, ib:ie, jb:je, kb:ke))
         end if

         call tileDesc%releaseDataPtr(solnData, CENTER)

         call itor%next()
      end do TileLoopResult

      call Grid_releaseTileIterator(itor)
   end subroutine setupTest

   subroutine runTest(success, mode, in_1, in_2, out_1, out_2, out_1_exp, out_2_exp, in_1_init, in_2_init)
      use Simulation_interface, only: Simulation_mapIntToStr

      implicit none

      logical, intent(out) :: success
      integer, intent(in) :: mode
      integer, intent(in) :: in_1, in_2
      integer, intent(in) :: out_1, out_2
      integer, intent(in) :: out_1_exp, out_2_exp
      integer, intent(in), optional :: in_1_init, in_2_init

      character(len=9) :: name
      character(len=4) :: in_1_name, in_2_name
      character(len=4) :: out_1_name, out_2_name
      character(len=4) :: out_1_exp_name, out_2_exp_name
      character(len=4) :: in_1_init_name, in_2_init_name

      integer, dimension(MDIM) :: ierr_1, ierr_2
      real :: err_1, err_2

      success = .true.

      name = mode_names(mode)

      call Simulation_mapIntToStr(in_1, in_1_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(in_2, in_2_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_1, out_1_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_2, out_2_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_1_exp, out_1_exp_name, MAPBLOCK_UNK)
      call Simulation_mapIntToStr(out_2_exp, out_2_exp_name, MAPBLOCK_UNK)

      if (present(in_1_init)) &
         call Simulation_mapIntToStr(in_1_init, in_1_init_name, MAPBLOCK_UNK)
      if (present(in_2_init)) &
         call Simulation_mapIntToStr(in_2_init, in_2_init_name, MAPBLOCK_UNK)

      if (eos_meshMe < maxPrintPE) then
         write (*, *) repeat("=", 88)
         write (*, mode_fmt) eos_meshMe, name, in_1_name, in_2_name, out_1_name, out_2_name
         write (*, *) repeat("=", 88)
      end if

      call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

      TileLoopInit: do
         if (.not. itor%isValid()) exit TileLoopInit

         call itor%currentTile(tileDesc)

         blockID = getBlockID()
         call Grid_getBlkType(blockId, nodeType)

         blkLimits = tileDesc%limits

         call tileDesc%getDataPtr(solnData, CENTER)

         ib = blkLimits(LOW, IAXIS); ie = blkLimits(HIGH, IAXIS)
         jb = blkLimits(LOW, JAXIS); je = blkLimits(HIGH, JAXIS)
         kb = blkLimits(LOW, KAXIS); ke = blkLimits(HIGH, KAXIS)

         if (present(in_1_init) .or. present(in_2_init)) then
            if (eos_meshMe < maxPrintPE) then
               write (*, block_fmt) eos_meshMe, blockID, nodeType
            end if
            if (present(in_1_init)) then
               write (*, set_fmt) eos_meshMe, in_1_name, in_1_init_name
               solnData(in_1, ib:ie, jb:je, kb:ke) = solnData(in_1_init, ib:ie, jb:je, kb:ke)
            end if

            if (present(in_2_init)) then
               write (*, set_fmt) eos_meshMe, in_2_name, in_2_init_name
               solnData(in_2, ib:ie, jb:je, kb:ke) = solnData(in_2_init, ib:ie, jb:je, kb:ke)
            end if
         end if

         if (out_1 .ne. TEMP_VAR) solnData(out_1, ib:ie, jb:je, kb:ke) = 0.0
         if (out_2 .ne. TEMP_VAR) solnData(out_2, ib:ie, jb:je, kb:ke) = 0.0

         call tileDesc%releaseDataPtr(solnData, CENTER)

         call itor%next()
      end do TileLoopInit

      call Grid_releaseTileIterator(itor)

      call Timers_start(name)
      call Eos_everywhere(mode)
      call Timers_stop(name)

      call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

      TileLoopCheck: do
         if (.not. itor%isValid()) exit TileLoopCheck

         call itor%currentTile(tileDesc)

         blockID = getBlockID()
         call Grid_getBlkType(blockId, nodeType)

         blkLimits = tileDesc%limits

         call tileDesc%getDataPtr(solnData, CENTER)

         ib = blkLimits(LOW, IAXIS); ie = blkLimits(HIGH, IAXIS)
         jb = blkLimits(LOW, JAXIS); je = blkLimits(HIGH, JAXIS)
         kb = blkLimits(LOW, KAXIS); ke = blkLimits(HIGH, KAXIS)

         if (eos_meshMe < maxPrintPE) then
            write (*, block_fmt) eos_meshMe, blockID, nodeType
         end if

         ierr_1 = maxloc(abs((solnData(out_1, ib:ie, jb:je, kb:ke) &
                              - solnData(out_1_exp, ib:ie, jb:je, kb:ke)) &
                             /solnData(out_1_exp, ib:ie, jb:je, kb:ke)))

         ierr_2 = maxloc(abs((solnData(out_2, ib:ie, jb:je, kb:ke) &
                              - solnData(out_2_exp, ib:ie, jb:je, kb:ke)) &
                             /solnData(out_2_exp, ib:ie, jb:je, kb:ke)))

         err_1 = abs((solnData(out_1, ierr_1(1), ierr_1(2), ierr_1(3)) &
                      - solnData(out_1_exp, ierr_1(1), ierr_1(2), ierr_1(3))) &
                     /solnData(out_1_exp, ierr_1(1), ierr_1(2), ierr_1(3)))

         err_2 = abs((solnData(out_2, ierr_2(1), ierr_2(2), ierr_2(3)) &
                      - solnData(out_2_exp, ierr_2(1), ierr_2(2), ierr_2(3))) &
                     /solnData(out_2_exp, ierr_2(1), ierr_2(2), ierr_2(3)))

         success = tolerance .gt. err_1
         success = success .and. (tolerance .gt. err_2)

         if (eos_meshMe < maxPrintPE) then
            write (*, err_fmt) eos_meshMe, out_1_name, err_1
            write (*, comp_fmt) eos_meshMe, in_1_name, solnData(in_1, ierr_1(1), ierr_1(2), ierr_1(3)), &
               in_2_name, solnData(in_2, ierr_1(1), ierr_1(2), ierr_1(3))
            write (*, comp_fmt) eos_meshMe, out_1_name, solnData(out_1, ierr_1(1), ierr_1(2), ierr_1(3)), &
               out_1_exp_name, solnData(out_1_exp, ierr_1(1), ierr_1(2), ierr_1(3))
            write (*, err_fmt) eos_meshMe, out_2_name, err_2
            write (*, comp_fmt) eos_meshMe, in_1_name, solnData(in_1, ierr_2(1), ierr_2(2), ierr_2(3)), &
               in_2_name, solnData(in_2, ierr_2(1), ierr_2(2), ierr_2(3))
            write (*, comp_fmt) eos_meshMe, out_2_name, solnData(out_2, ierr_2(1), ierr_2(2), ierr_2(3)), &
               out_2_exp_name, solnData(out_2_exp, ierr_2(1), ierr_2(2), ierr_2(3))

            if (success) then
               write (*, status_fmt) eos_meshMe, name, "PASS"
            else
               write (*, status_fmt) eos_meshMe, name, "FAIL"
            end if

            write (*, *) repeat("=", 88)
         end if

         call tileDesc%releaseDataPtr(solnData, CENTER)

         call itor%next()
      end do TileLoopCheck

      call Grid_releaseTileIterator(itor)

   end subroutine runTest

end subroutine Eos_unitTest

