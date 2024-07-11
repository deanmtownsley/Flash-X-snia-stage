!!****if* source/physics/sourceTerms/Flame/FlameEffects/Flame_rhJump
!!
!! NAME
!!
!!  Flame_rhJump
!!
!! SYNOPSIS
!!
!!  call Flame_rhJump (
!!                     real(inout) :: eosData_u(EOS_NUM),
!!                     real(inout) :: eosData_b(EOS_NUM),
!!                     real(in)    :: q,
!!                     real(in)    :: s,
!!                     integer(in) :: mode,
!!           optional, real(in)    :: mfrac_u(NSPECIES),
!!           optional, real(in)    :: mfrac_b(NSPECIES) )
!!
!! DESCRIPTION
!!
!!  Computes burned state from a Rankine-Hugoniot jump condition across the flame front
!!  
!!  computes the thermodynamic properties of the ash (dens,temp) given
!!  those of the fuel, the energy release q, the front speed s (w.r.t. the fuel), and the
!!  composition of the ash material
!!
!! NOTES
!!
!!  The unburned state is determined by calling the Eos in mode "mode"
!!  using the supplied eosData_u and mfrac_u.
!!
!!  This routine works with or without USE_EOS_LITE, but the
!!  "compositon" of the burned state is set in different ways:
!!
!!     with USE_EOS_LITE: set eosData_b(EOS_ABAR) and eosData_b(EOS_ZBAR)
!!     without:           set mfrac_b
!!
!! ARGUMENTS
!!
!!   eosData_u - thermodynamic state info about unburned material
!!   eosData_b - thermodynamic state info about burned material
!!      q      - energy release (erg/g)
!!      s      - flame speed (cm/s)
!!      mode   - EOS mode used to calculate full unburned state
!!    mfrac_u  - composition of unburned state (if not present abar and zbar from eosData_u is used)
!!    mfrac_b  - composition of burned state (if not present abar and zbar from eosData_b is used)
!!
!! SEE ALSO
!!
!!  See Flame_interface.F90 for possible updates
!!
!! HISTORY
!!
!!  Dean Townsley 2007,2008
!!
!!***

#include "Simulation.h"

#ifdef DEBUG_ALL
#define DEBUG_FLAME
#endif

subroutine Flame_rhJump(eosData_u, eosData_b, q, s, mode, mfrac_u, mfrac_b)

  use fl_effData, ONLY : fl_effsmlrho
  use Eos_interface, ONLY : Eos_vector
  
  implicit none

#include "constants.h"
#include "Eos.h"  

  real, dimension(EOS_NUM),  intent(inout) :: eosData_u
  real, dimension(EOS_NUM),  intent(inout) :: eosData_b
  real,    intent(in)     :: q, s
  integer, intent(in)     :: mode  !! This is the Eos mode
  real, optional, dimension(NSPECIES), intent(in) :: mfrac_u
  real, optional, dimension(NSPECIES), intent(in) :: mfrac_b

  real, dimension(1,EOS_NUM) :: eosData_buff
  real, dimension(1,EOS_VARS+1:EOS_NUM) :: derivs_u
  real, dimension(1,EOS_VARS+1:EOS_NUM) :: derivs_b
  real                    :: dfd, dft, dgd, dgt, f, g
  real                    :: sq_s_dens, d_inv_dens, determinant_inv
  real                    :: dens_b_old, temp_b_old
  real                    :: dens_u, temp_u, ener_u, pres_u, dens_b, temp_b, ener_b, pres_b
  real                    ::  error
  integer                 :: niters
1 format(5(2X, E10.4))

  eosData_buff(1,:)=eosData_u(:)
  if (present(mfrac_u)) then
     call Eos_vector(mode,1,eosData_buff,mfrac_u,derivs_u)
  else
     call Eos_vector(mode,1,eosData_buff,derivs=derivs_u)
  endif
  eosData_u(:)=eosData_buff(1,:)


#ifdef DEBUG_FLAME
  write (*,*) '------------------------------------------------------'
  write (*,*) ''
  write (*,*) 'rhjump: flame speed, heat release:'
  write (*,1) s, q
  write (*,*) 'rhjump: unburned temp, dens, ener, pres:'
  write (*,1) eosData_u((/EOS_TEMP, EOS_DENS, EOS_EINT, EOS_PRES/))
  write (*,*) ''
#endif

!.. initial (not very good) guess for burned state:

  eosData_b(EOS_DENS) = eosData_u(EOS_DENS)
  eosData_b(EOS_EINT) = eosData_u(EOS_EINT) + q
  ! need a guess temperature
  eosData_b(EOS_TEMP) = eosData_u(EOS_TEMP)
  if ( (.not.present(mfrac_b)) .and. (eosData_b(EOS_ABAR) == 0.0 ) ) then
     call Driver_abort("Must set Abar of burned material in Interfaces_rhjump with no mfrac")
  endif

  eosData_buff(1,:)=eosData_b(:)
  if (present(mfrac_b)) then
     call Eos_vector(MODE_DENS_EI,1,eosData_buff,mfrac_b,derivs_b)
  else
     call Eos_vector(MODE_DENS_EI,1,eosData_buff,derivs=derivs_b)
  endif
  eosData_b(:)=eosData_buff(1,:)


  ! some renames for readability
  dens_u = eosData_u(EOS_DENS)
  temp_u = eosData_u(EOS_TEMP)
  ener_u = eosData_u(EOS_EINT)
  pres_u = eosData_u(EOS_PRES)
  sq_s_dens = (s*dens_u)**2
  dens_b = eosData_b(EOS_DENS)
  temp_b = eosData_b(EOS_TEMP) 
  ener_b = eosData_b(EOS_EINT)
  pres_b = eosData_b(EOS_PRES)

  ! 2d newton loop to find dens_b and temp_b which solve
  !    pres_b = pres_u - (dens_u*s)**2*(1/dens_b - 1/dens_u)
  !    ener_b = ener_u + q - 0.5 ( pres_b+ pres_u)*(1/dens_b-1/dens_u)

  !  the ordering of evaluations and tests here is a little weird -- should be fixed (DMT 2007/4/14)
  error = 1.
  niters = 0
  do while (error > 1.e-8 .and. niters < 100)


     d_inv_dens = 1./dens_b - 1./dens_u

     f = pres_b - pres_u + sq_s_dens * d_inv_dens 
     g = ener_b - ener_u - q + 0.5*(pres_b + pres_u) * d_inv_dens

     dfd = derivs_b(1,EOS_DPD) - sq_s_dens/dens_b**2
     dft = derivs_b(1,EOS_DPT)
     dgd = derivs_b(1,EOS_DED) + 0.5*d_inv_dens*derivs_b(1,EOS_DPD) - 0.5*(pres_b + pres_u) / dens_b**2
     dgt = derivs_b(1,EOS_DET) + 0.5*d_inv_dens*derivs_b(1,EOS_DPT)

     determinant_inv = 1./(dfd*dgt - dft*dgd)

     dens_b_old = dens_b
     temp_b_old = temp_b

     dens_b = dens_b - (f*dgt - g*dft) * determinant_inv
     temp_b = temp_b + (f*dgd - g*dfd) * determinant_inv

     if (dens_b .lt. fl_effsmlrho) then

        dens_b = 0.5*dens_b_old
        temp_b = temp_b_old
    
     elseif (temp_b .lt. temp_u) then

        temp_b = temp_u
        dens_b = dens_b_old

     endif

     ! un-rename new values
     eosData_b(EOS_DENS)=dens_b
     eosData_b(EOS_TEMP)=temp_b

!     write (6,*) 'in loop call', dens_b, temp_b, mfrac_b
     eosData_buff(1,:)=eosData_b(:)
     if (present(mfrac_b)) then
        call Eos_vector(MODE_DENS_TEMP,1,eosData_buff,mfrac_b,derivs_b)
     else
        call Eos_vector(MODE_DENS_TEMP,1,eosData_buff,derivs=derivs_b)
     endif
     eosData_b(:)=eosData_buff(1,:)     

     ! rename
     pres_b=eosData_b(EOS_PRES)
     ener_b=eosData_b(EOS_EINT)

     error = abs(f/pres_u) + abs(g/ener_u)

#ifdef DEBUG_FLAME
     write(*,1) temp_b, dens_b, ener_b, pres_b, error
#endif

     niters = niters + 1

  enddo

#ifdef DEBUG_FLAME
  write (*,*) ''
  write (*,*) 'rhjump: burned temp, dens, ener, pres:'
  write (*,1) temp_b, dens_b, ener_b, pres_b
  write (*,*) 'rhjump:   niters, error:', niters, error
  write (*,*) ''
  write (*,*) '-------------------------------------------------'
#endif
   if (niters >= 100) write (6,*) 'rhjump did not converge'

!------------------------------------------------------------------

end subroutine Flame_rhJump
