!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/Driver_initAll
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
!!  Driver_initAll
!!
!! SYNOPSIS
!!
!!   Driver_initAll()
!!
!! DESCRIPTION
!!
!!  Performs Flash initializations, which includes:
!!
!!  Call all 'init' routines in units.  Order does matter,
!!  particularly when restarting from a checkpoint file.
!!
!!  For the most part, Driver_initAll calls other units' init
!!  routines directly, like call IO_init or call Grid_init.  This
!!  routine also makes calls to other Driver initialization routines
!!  like Driver_initMaterialProperties or Driver_initSourceTerms.
!!  These routines then call the unit-specific initialization 
!!  routines.  This level of abstraction was added to simplify
!!  the initialization calls.
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!! HISTORY 
!!
!!  2008-03-14 KW   Moved material properties initialization up.
!!
!!***


subroutine Driver_initAll()
  
  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs, dr_globalComm, dr_nbegin, &
       dr_initialSimTime, dr_elapsedWCTime, &
       dr_initialWCTime, dr_restart, dr_dtInit, dr_redshift,dr_particlesInitialized

  use Driver_interface, ONLY : Driver_init, &
    Driver_verifyInitDt, Driver_abort, Driver_setupParallelEnv
  use RuntimeParameters_interface, ONLY : RuntimeParameters_init, RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_init
  use PhysicalConstants_interface, ONLY : PhysicalConstants_init
  use Gravity_interface, ONLY : Gravity_init, &
    Gravity_potential
  use Timers_interface, ONLY : Timers_init, Timers_start, Timers_stop

  use Grid_interface, ONLY : Grid_init, Grid_initDomain
  use Particles_interface, ONLY : Particles_init,  Particles_initData 
  use Eos_interface, ONLY : Eos_init
  use Hydro_interface, ONLY : Hydro_init
  use Simulation_interface, ONLY : Simulation_init
  use IO_interface, ONLY :IO_init, IO_outputInitial
  use Gravity_interface, ONLY :  Gravity_potential
  implicit none       
  
#include "constants.h"
#include "Simulation.h"

  logical :: updateRefine

  dr_elapsedWCTime = 0.0

  !! Initialize message-passing interface and mesh package.
  !! This step and the timer initialization must not depend
  !! on any runtime parameters, because they (the parameters)
  !! have not yet been initialized.

  !! hand process ID out to C routines to avoid architecture-dependent code
  call driver_abortc_set_mype(dr_globalMe)

  !! make sure our stack (and whatever other rlimits) are big enough.  
  !! this should get around the 2Mb stack limit that pthreads
  !! imposes if linked statically (but not dynamically!)
  call dr_set_rlimits(dr_globalMe)

    
  !! Initialize runtime parameters
  call RuntimeParameters_init(  dr_restart)
  call Driver_setupParallelEnv()

  !! Initialize the code timers.  Ideally should be first thing in
  !! code but currently the timing package
  !! uses MPI_WTime(), so Driver_initParallel() must go first, and
  !! uses RuntimeParameters_get(), so RuntimeParameters_init() must go
  !! first.
  call Timers_init(  dr_initialWCTime)
  call Timers_start("initialization")


  !PhysicalConstants init and Multispecies init must come before Logfile
  !since their values are stamped to the logfile
  call PhysicalConstants_init( )

  call Logfile_init( )

  call Grid_init( )
  
  if(dr_globalMe==MASTER_PE)print*,'MaterialProperties initialized'
  

  call RuntimeParameters_get('dtInit',dr_dtInit)

  call Particles_init(  dr_restart)       ! Particles
  
#ifdef DEBUG_DRIVER
  if(dr_globalMe==MASTER_PE)print*,'Particles initialized'
#endif

  if(.not. dr_restart) then     
     

     call Driver_init()

     !Eos must come before Grid
     call Eos_init()

       !must come before Grid since simulation specific values must go on the Grid

     call Simulation_init()


     call Grid_initDomain(dr_restart,dr_particlesInitialized)
     if (dr_globalMe==MASTER_PE)print *, ' Finished with Grid_initDomain, no restart'

     call IO_init( )

  else if(dr_restart) then
     
     call IO_init( )


     call Driver_init()

     call Eos_init()

     call Simulation_init()
     dr_particlesInitialized=.true.
     call Grid_initDomain( dr_restart,dr_particlesInitialized)
     if (dr_globalMe==MASTER_PE) print *, ' Finished with Grid_initDomain, restart'
     
  end if

  !Hydro_init must go before Driver
  if(dr_globalMe==MASTER_PE) print *, 'Ready to call Hydro_init' 
  call Hydro_init()           ! Hydrodynamics, MHD, RHD
  if(dr_globalMe==MASTER_PE)print*,'Hydro initialized'
  
  
  call Gravity_init()         ! Gravity
  if(dr_globalMe==MASTER_PE)print*,'Gravity initialized'


  call Driver_verifyInitDt()
  if(dr_globalMe==MASTER_PE)print*,'Initial dt verified'
 
  !For active particle simulations we must initialize particle 
  !positions before the call to Gravity_potential.
  call Particles_initData(dr_restart,dr_particlesInitialized)
  
  call IO_outputInitial( dr_nbegin, dr_initialSimTime)
  if(dr_globalMe==MASTER_PE)print*,'Initial plotfile written'

  if(dr_globalMe==MASTER_PE)print*,'Driver init all done'

  !!Done with initialization.
  call Timers_stop ("initialization")


  return
end subroutine Driver_initAll
