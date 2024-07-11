
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
   real, intent(in), dimension(SPECIES_BEGIN:SPECIES_END) :: abund

   logical :: mask(EOS_VARS+1:EOS_NUM)
   real :: guess_dens, dp_drho_T, target_pres, new_guess_dens
   integer :: iter
   integer, parameter :: max_iter = 100
   real :: tol, err

   tol = eos_tol
   mask(:) = .false.
   mask(EOS_DPD) = .true.

   iter = 0
       
   ! get pressure to match from info in thermstate
   call Eos_vector( MODE_DENS_TEMP, 1, thermstate, abund)
   target_pres = thermstate(EOS_PRES)

   thermstate(EOS_TEMP) = ign_temp

   guess_dens = thermstate(EOS_DENS)
   err=tol*2.0
   do while ( iter < max_iter .and. err > tol )
      iter = iter + 1

      ! evaluate
      thermstate(EOS_DENS) = guess_dens
      call Eos_vector(MODE_DENS_TEMP, 1, thermstate, massFrac=abund)!, mask=mask)
      dp_drho_T = thermstate( EOS_DPD )

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
