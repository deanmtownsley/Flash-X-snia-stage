!!****if* source/Particles/localAPI/pt_advanceMC
!!
!! NAME
!!
!!  pt_advanceMC
!!
!! SYNOPSIS
!!
!!  call pt_advanceCharged(real(in)    :: dtnew,
!!                         integer(in) :: iSpecies)
!!
!! DESCRIPTION
!!
!!   Advances particles in time
!!
!! ARGUMENTS
!!
!!   dtNew : current time interval
!!   iSpecies - which species we are advancing
!!
!!
!!
!!***

subroutine pt_advanceMC(dtNew, iSpecies)

  use Particles_data, only : useParticles, particles,&
                             pt_half_rt_timesteps
  use Timers_interface, only : Timers_start, Timers_stop  
  use emission, only : emit_mcps
  use transport, only : transport_mcps
  use rhd, only : apply_rad_source_terms
  implicit none


#include "Simulation.h"
#include "constants.h"
 
  ! Input/Output 
  real, intent(in)  :: dtNew
  integer, intent(in) :: iSpecies

  ! Other parameters
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: numofblks
  real, pointer :: solnVec(:,:,:,:)
  real :: dt_step

  ! Quit if particle module is not in use
  if (.not. useParticles) return


  ! Start performing RT
  call Timers_start("Radiation Transport")

  ! Adjust time step for split hydro solver
  dt_step = dtNew
  if (pt_half_rt_timesteps) dt_step = 0.50*dtNew

  ! Emission of radiation, including both thermal, point, and face
  call emit_mcps(dt_step, iSpecies)

  ! Radiation transport
  call transport_mcps(dt_step, iSpecies)

  ! Deposit radiation source terms
  call apply_rad_source_terms(dt_step)

  ! End of current RT step
  call Timers_stop("Radiation Transport")
  
end subroutine pt_advanceMC
