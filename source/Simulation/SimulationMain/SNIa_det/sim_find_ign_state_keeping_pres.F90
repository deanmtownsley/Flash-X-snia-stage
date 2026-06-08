
!! Dean M. Townsley 2016
!!
!! do a simple root find to get the density for the given temperature and
!! pressure.
!! on input thermstate() contains background (non-ignited) state
!! on output thermstate() containts ignited state at ign_temp and density
!! such that pressure is same as in non-ignited state

subroutine sim_find_ign_state_keeping_pres(ign_temp, thermstate, abund)

   use Eos_interface, ONLY : Eos
   use eos_helmData, ONLY : eos_tol
   use Driver_interface, ONLY : Driver_abort

   implicit none
#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

   real, intent(in) :: ign_temp
   real, intent(inout), dimension(EOS_NUM) :: thermstate
   real, intent(in), dimension(SPECIES_BEGIN:SPECIES_END) :: abund

   real :: guess_dens, dp_drho_T, target_pres, new_guess_dens
   integer :: iter
   integer, parameter :: max_iter = 100
   real :: tol, err

   tol = eos_tol

   iter = 0
       
   ! get pressure to match from info in thermstate
   call Eos(MODE_DENS_TEMP, thermstate(EOS_PRES), thermstate(EOS_TEMP), &
            thermstate(EOS_DENS), thermstate(EOS_GAMC), thermstate(EOS_EINT), &
            thermstate(EOS_ENTR), thermstate(EOS_ABAR), thermstate(EOS_ZBAR), &
            thermstate(EOS_YE), massFrac=abund )
   target_pres = thermstate(EOS_PRES)

   thermstate(EOS_TEMP) = ign_temp

   guess_dens = thermstate(EOS_DENS)
   err=tol*2.0
   do while ( iter < max_iter .and. err > tol )
      iter = iter + 1

      ! evaluate
      thermstate(EOS_DENS) = guess_dens
      call Eos(MODE_DENS_TEMP, thermstate(EOS_PRES), thermstate(EOS_TEMP), &
               thermstate(EOS_DENS), thermstate(EOS_GAMC), thermstate(EOS_EINT), &
               thermstate(EOS_ENTR), thermstate(EOS_ABAR), thermstate(EOS_ZBAR), &
               thermstate(EOS_YE), massFrac=abund )
      dp_drho_T = thermstate( EOS_DPD )

      ! newton-raphson improved guess
      new_guess_dens = guess_dens + (target_pres-thermstate(EOS_PRES))/dp_drho_T
      err = abs(new_guess_dens-guess_dens)*2.0/(new_guess_dens+guess_dens)

      guess_dens = new_guess_dens
   enddo


   if ( .not. (err <= tol) ) then
      ! NaN error will also fail test
      call Driver_abort("Failed to converge in sim_find_ignition_state_keeping_pressure")
   endif

   thermstate(EOS_DENS) = guess_dens

   return

end subroutine
