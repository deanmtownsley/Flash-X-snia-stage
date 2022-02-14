!!****if* source/physics/Eos/unitTest/WeakLib/Eos_unitTest
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
!!  PARAMETERS
!!
!!  eintSwitch  a rarely used switch which ensures that internal energy calculations 
!!        maintain sufficient precision. Important only if energyTotal is dominated 
!!        by energyKinetic.
!!
!!***

!!REORDER(4): solnData

subroutine Eos_unitTest(fileUnit, perfect)

  use Eos_interface, ONLY : Eos_everywhere
  use Grid_interface,ONLY : Grid_getTileIterator, &
                            Grid_releaseTileIterator, &
                            Grid_getBlkType
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,     ONLY : Grid_tile_t
  use IO_interface, ONLY : IO_writeCheckpoint
  use Eos_data, ONLY : eos_meshMe, eos_meshNumProcs
  use eos_testData, ONLY: eos_testPresModeStr, &
                          eos_testEintModeStr, &
                          eos_testTempModeStr, &
                          eos_testPresMode, &
                          eos_testEintMode, &
                          eos_testTempMode
  use eos_testData,  ONLY : tolerance => eos_testTolerance
  implicit none

# include "Eos.h"
# include "constants.h"
# include "Simulation.h"

  integer, intent(in) :: fileUnit
  logical, intent(out) :: perfect
  integer :: blockID
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  real :: presErr, tempErr, eintErr

  real, pointer, dimension(:,:,:,:):: solnData
  logical:: test1,test2,test3,test4 !for a block
  logical:: test1allB,test2allB,test3allB,test4allB !for all blocks

  integer :: vecLen, blockOffset,  pres, dens, temp, e, n, m
  integer :: isize, jsize, ksize, i,j,k, nStartsAtOne
  real, dimension(:), allocatable :: eosData
  real, dimension(:), allocatable :: massFrac
  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
  real, allocatable, dimension(:,:,:,:) :: derivedVariables
  real, dimension(:,:,:), allocatable :: deriv1, deriv2

  character(len=7),pointer:: ap
  character(len=7),target :: a
  integer,parameter :: maxPrintPE = 20
  integer,save :: nodeType = LEAF
  integer :: ib,ie,jb,je,kb,ke
  integer, dimension(3) :: startingPos, dataSize
  real :: presErr1, presErr2
  real :: tempErr1, tempErr2
  real :: eintErr1, eintErr2

  nullify(solnData)

  if (eos_meshNumProcs==1) then
     a = ''
     ap => a
  else
20   format(I6,':')
     write(a,20) eos_meshMe
     a = trim(adjustl(a))
     ap => a
  end if

! info for checkpoint output


  mask = .true.

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     !! In Simulation_initBlock,
     !! temperature is initialized in CTMP_VAR and pressure is
     !! initialized in CPRS_VAR. We don't change these variables
     !! we copy them into the usual variable name as needed.

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     ! Testing density/temperature in; energy/pressure out
     if (eos_meshMe<maxPrintPE) print *,ap,'Block',blockID,' type',nodeType
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_TEMP or similar: Density, temperature in; energy, pressure out; mode=', &
          eos_testTempMode,eos_testTempModeStr
     if (eos_meshMe<maxPrintPE) then
        print*,ap,'The initialized extreme values are '
        print*,ap,'Initialized Density min',minval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Density max',maxval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Temperature min',minval(abs(solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Temperature max',maxval(abs(solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Initialized Pressure min',minval(solnData(CPRS_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'Initialized Pressure max',maxval(solnData(CPRS_VAR,ib:ie,jb:je,kb:ke))
     end if 

     solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Eos_everywhere(eos_testTempMode)

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

    !! Summarize results of MODE_DENS_TEMP (or similar) call
    if (eos_meshMe<maxPrintPE) then
        print*,ap,'The resulting extreme values are '
!!$        print*,ap,'Resulting Density min',minval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
!!$        print*,ap,'Resulting Density max',maxval(abs(solnData(DENS_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting Temperature min',minval(abs(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting Temperature max',maxval(abs(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting Pressure min',minval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'Resulting Pressure max',maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'Resulting internal energy min',minval(abs(solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
        print*,ap,'Resulting internal energy max',maxval(abs(solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     end if

     ! Save the equilibrium values in O variables and (hopefully) write them out
     solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
     solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
     solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call IO_writeCheckpoint()   !! This is checkpoint 001

  test1allB = .TRUE.
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     !  Testing density/energy in, temperature/pressure out
     if (eos_meshMe<maxPrintPE) print *,ap,'Block',blockID,' type',nodeType
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI or similar: Density, energy in; temperature, pressure out; mode=', &
       eos_testEintMode,eos_testEintModeStr

     !  Zero output variables
     !  solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=0  ! don't zero TEMP or eos_helm cannot converge in MODE_DENS_EI
     solnData(PRES_VAR,:,:,:)=0 

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Eos_everywhere(eos_testEintMode)

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     if (eos_meshMe<maxPrintPE) then !! Summarize results of MODE_DENS_EI (or similar) call
        print*,ap,'  Temperature min ',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Temperature max ',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Pressure min ',minval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Pressure max ',maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     end if

     !! Calculate error from MODE_DENS_EI (or similar) call.
     tempErr = maxval(   &
             abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)- solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/ &
                 solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))  )
     if (eos_meshMe<maxPrintPE) then
        print*,ap,'  The calculated error in temperature is ',tempErr
     end if

     presErr = maxval(abs((solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/&
          solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) print*,ap,'  The calculated error in pressure is ',presErr

     test1 = (tolerance > tempErr)
     test1 = test1.and.(tolerance > presErr)
     if(test1) then
         if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_EI or similar is fine'
     else
        if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI or similar is BAD!!!'
        test1allB = .FALSE.
     endif

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call IO_writeCheckpoint()  !! This is checkpoint 002

  test2allB = .TRUE.
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     ! Testing density/pressure in, energy/temp out
     if (eos_meshMe<maxPrintPE) print *,ap,'Block',blockID,' type',nodeType
     if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_PRES or similar: Density, pressure in; energy, temperature out; mode=', &
       eos_testPresMode,eos_testPresModeStr

     !solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)
     !solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)
     solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0.0

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Eos_everywhere(eos_testPresMode)

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     !! Summarize results of MODE_DENS_PRES (or similar) call;
     !! calculate error from MODE_DENS_PRES (or similar) call.
     tempErr = maxval(abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
     eintErr = maxval(abs((solnData(EINT_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OENT_VAR,ib:ie,jb:je,kb:ke))/solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) then
        print*,ap,'  Energy min is',minval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Energy max is',maxval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Temperature min is',minval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  Temperature max is',maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
        print*,ap,'  The calculated error in energy is ',eintErr
        print*,ap,'  The calculated error in temperature is ',tempErr
     end if

     test2 = (tolerance > tempErr)
     test2 = test2.and.(tolerance > eintErr)
     if(test2) then
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_PRES or similar is fine'
     else
        if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_PRES or similar is BAD!!!'
        test2allB = .FALSE.
     endif

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call IO_writeCheckpoint()   !! This is checkpoint 003

  test3allB = .TRUE.
  test4allB = .TRUE.
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     if (eos_meshMe<maxPrintPE) print*,ap,'And now to verify the other solutions'

     !!  Do calculations in reverse order
     ! Density and pressure in, energy and temperature out
     solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)
     solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)
     solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0.0

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Eos_everywhere(MODE_DENS_PRES)

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     ! Now we have a "true"  temperature and internal energy; save them for comparison
     solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
     solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
     solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)

     ! Density and energy in, temperature and pressure out
     !! zero output values to make sure they're being calculated
     solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=0.0
     !solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)=solnData(CTMP_VAR,ib:ie,jb:je,kb:ke)

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Eos_everywhere(MODE_DENS_EI)

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     tempErr1 = maxval(solnData(TEMP_VAR,ib:ie,jb:je,kb:ke))
     tempErr2 = maxval(solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'maxval TEMP_VAR OTMP_VAR',tempErr1,tempErr2

     presErr1 = maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     presErr2 = maxval(solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'maxval PRES_VAR OPRS_VAR',presErr1,presErr2

     presErr = maxval(abs(&
          (solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/&
           solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in pressure from EI is ',presErr

     ! NOTE this ALWAYS comes out to zero. suspect  something wrong here....
     tempErr = maxval(abs((solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OTMP_VAR,ib:ie,jb:je,kb:ke))/solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in temperature from EI is ',tempErr

     test3 = (tolerance > tempErr)
     test3 = test3.and.(tolerance > presErr)

     if(test3) then
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_EI is fine '
     else
        if (eos_meshMe<maxPrintPE) print *,ap,'MODE_DENS_EI is BAD!!!'
        test3allB = .FALSE.
     endif

     solnData(EINT_VAR,ib:ie,jb:je,kb:ke)=0.0
     solnData(PRES_VAR,ib:ie,jb:je,kb:ke)=0.0
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  call Eos_everywhere(MODE_DENS_TEMP)

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id     ! only used for some useful screen output
#else
     blockID = tileDesc % grid_index  ! only for some useful output
#endif
     call Grid_getBlkType(blockId,nodeType)
     call tileDesc%getDataPtr(solnData, CENTER)
     blkLimits = tileDesc%limits

     ib=blkLimits(LOW,IAXIS) ; ie=blkLimits(HIGH,IAXIS)
     jb=blkLimits(LOW,JAXIS) ; je=blkLimits(HIGH,JAXIS)
     kb=blkLimits(LOW,KAXIS) ; ke=blkLimits(HIGH,KAXIS)

     eintErr1 = maxval(solnData(EINT_VAR,ib:ie,jb:je,kb:ke))
     eintErr2 = maxval(solnData(OENT_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'maxval EINT_VAR OENT_VAR',eintErr1,eintErr2

     presErr1 = maxval(solnData(PRES_VAR,ib:ie,jb:je,kb:ke))
     presErr2 = maxval(solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))
     if (eos_meshMe<maxPrintPE) print*,ap,'maxval PRES_VAR OPRS_VAR',presErr1,presErr2

     presErr = maxval(abs((solnData(PRES_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OPRS_VAR,ib:ie,jb:je,kb:ke))/solnData(PRES_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in pressure is ',presErr

     eintErr = maxval(abs((solnData(EINT_VAR,ib:ie,jb:je,kb:ke)-&
          solnData(OENT_VAR,ib:ie,jb:je,kb:ke))/solnData(EINT_VAR,ib:ie,jb:je,kb:ke)))
     if (eos_meshMe<maxPrintPE) print*,ap,'The calculated error in energy is ',eintErr

     test4 = (tolerance > presErr)
     test4 = test4.and.(tolerance > eintErr)
     if(test4) then
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_TEMP is fine'
     else
        if (eos_meshMe<maxPrintPE) print*,ap,'MODE_DENS_TEMP is BAD!!!!'
        test4allB = .FALSE.
     endif
     solnData(OPRS_VAR,ib:ie,jb:je,kb:ke)=solnData(PRES_VAR,ib:ie,jb:je,kb:ke)
     solnData(OENT_VAR,ib:ie,jb:je,kb:ke)=solnData(EINT_VAR,ib:ie,jb:je,kb:ke)
     solnData(OTMP_VAR,ib:ie,jb:je,kb:ke)=solnData(TEMP_VAR,ib:ie,jb:je,kb:ke)
     
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  !! Output to get the derived variables
  call IO_writeCheckpoint()

  if (eos_meshMe.EQ.MASTER_PE) print*,'out of the loop'
  perfect = test1allB.and.test2allB.and.test3allB.and.test4allB
  if(perfect) then
     if (eos_meshMe<maxPrintPE) print*,ap,'SUCCESS all tests were fine'
  else
     if (eos_meshMe<maxPrintPE) print*,ap,'FAILURE some tests failed'
  end if
  return
end subroutine Eos_unitTest




