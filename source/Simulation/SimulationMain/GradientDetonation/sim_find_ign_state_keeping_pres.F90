
!! Dean M. Townsley 2016
!!
!! do a simple root find to get the density for the given temperature and
!! pressure.
!! on input thermstate() contains background (non-ignited) state
!! on output thermstate() containts ignited state at ign_temp and density
!! such that pressure is same as in non-ignited state

subroutine sim_find_ign_state_keeping_pres(ign_temp, thermstate, abund)

   use Eos_interface, ONLY : Eos_vector
   use eos_helmData, ONLY : eos_tol
   use Driver_interface, ONLY : Driver_abort

   implicit none
#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

   real, intent(in) :: ign_temp
   real, intent(inout), dimension(EOS_NUM) :: thermstate
   real, dimension(1,EOS_VARS +1:EOS_NUM) :: derivs
   real, dimension(1,EOS_VARS) :: thermstate_buffer
   real, intent(in), dimension(SPECIES_BEGIN:SPECIES_END) :: abund
   real, dimension(1,NSPECIES) :: abund_buffer

   !logical :: mask(EOS_VARS+1:EOS_NUM)
   real :: guess_dens, dp_drho_T, target_pres, new_guess_dens
   integer :: iter
   integer, parameter :: max_iter = 100
   real :: tol, err

   tol = eos_tol
   !mask(:) = .false.
   !mask(EOS_DPD) = .true.

   iter = 0
   abund_buffer(1,:) = abund(SPECIES_BEGIN:SPECIES_END)
   write(0,*) "Before first Eos call"

   ! get pressure to match from info in thermstate
   thermstate_buffer(1,:) = thermstate(1:EOS_VARS)
   call Eos_vector( MODE_DENS_TEMP, 1, thermstate_buffer, abund_buffer)
   !abund(SPECIES_BEGIN:SPECIES_END) = abund_buffer(1,:)
   thermstate(1:EOS_VARS) = thermstate_buffer(1,:)

   write(0,*) "After first Eos call"
   target_pres = thermstate(EOS_PRES)

   thermstate(EOS_TEMP) = ign_temp

   guess_dens = thermstate(EOS_DENS)
   err=tol*2.0
   do while ( iter < max_iter .and. err > tol )
      iter = iter + 1

      ! evaluate
      thermstate(EOS_DENS) = guess_dens
      thermstate_buffer(1,:) = thermstate(1:EOS_VARS)
      !abund_buffer(1:NSPECIES) = abund(SPECIES_BEGIN: SPECIES_END)

      call Eos_vector(MODE_DENS_TEMP, 1, thermstate_buffer, massFrac=abund_buffer, derivs=derivs)!, mask=mask)
      thermstate(1:EOS_VARS) = thermstate_buffer(1,:)

      dp_drho_T = derivs(1,EOS_DPD)

      ! newton-raphson improved guess
      new_guess_dens = guess_dens + (target_pres-thermstate(EOS_PRES))/dp_drho_T
      new_guess_dens = max( 0.0, new_guess_dens )
      ! but don't change too much in one step (assumes positive)
      new_guess_dens = max( guess_dens*0.5 , min( guess_dens*1.5, new_guess_dens ) )
      err = abs((new_guess_dens-guess_dens)*2.0/(new_guess_dens+guess_dens))

      guess_dens = new_guess_dens
   enddo


   if ( .not. (err <= tol) ) then
      ! NaN error will also fail test
      print *, 'ign_temp=', ign_temp
      print *, 'targpres=', target_pres
      print *, 'guessden=', guess_dens
      print *, 'iter=', iter
      print *, 'err=', err
      call Driver_abort("Failed to converge in sim_find_ignition_state_keeping_pressure")
   endif

   thermstate(EOS_DENS) = guess_dens

   return

end subroutine
