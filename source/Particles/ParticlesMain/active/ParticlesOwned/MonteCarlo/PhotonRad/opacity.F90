module opacity 

  implicit none
  contains

!====================!
! ABSORPTION OPACITY !
!====================!
subroutine calc_abs_opac_vectorized(rho, temp, energy, Ye, iSpecies, ka)
  use Particles_data, only : pt_dens_threshold, pt_grey_abs_opac, pt_is_grey
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, intent(in) :: rho(:), temp(:), energy(:), Ye(:)
  integer, intent(in) :: iSpecies
  real, intent(out) :: ka(:)
  
  if (pt_is_grey) then
    ka = pt_grey_abs_opac
    ka = ka * rho
  else
    call Driver_abortFlash("calc_abs_opac:&
                           non-grey emission not implemented yet!")
  end if
  
  where (rho .lt. pt_dens_threshold)
    ka = 0.0d0
  end where
end subroutine calc_abs_opac_vectorized

subroutine calc_abs_opac(cellID, solnVec, energy, iSpecies, ka)
  use Driver_interface,  ONLY : Driver_abortFlash
  implicit none

#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: energy
  integer, intent(in) :: iSpecies
  real, intent(out) :: ka

  ! aux variables
  real :: rho, temp, Ye
  real :: tmp(1)

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  Ye = 0 ! not used

  call calc_abs_opac_vectorized( (/rho/), (/temp/), (/energy/), (/Ye/), iSpecies, tmp)
  ka = tmp(1)
end subroutine calc_abs_opac

!====================!
! SCATTERING OPACITY !
!====================!
subroutine calc_sca_opac_vectorized(rho, temp, energy, Ye, iSpecies, ks)
  use Particles_data, only : pt_dens_threshold, pt_grey_sca_opac, pt_is_grey
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, intent(in) :: rho(:), temp(:), energy(:), Ye(:)
  integer, intent(in) :: iSpecies
  real, intent(out) :: ks(:)

  if (pt_is_grey) then
    ks = pt_grey_sca_opac
    ks = ks * rho
  else
    call Driver_abortFlash("calc_sca_opac:&
                            non-grey emission not implemented yet!")
  end if

  where (rho .lt. pt_dens_threshold)
    ks = 0.0d0
  end where
end subroutine calc_sca_opac_vectorized


subroutine calc_sca_opac(cellID, solnVec, energy, iSpecies, ks)
  use Driver_interface,  ONLY : Driver_abortFlash
  implicit none

#include "Flash.h"

  ! Input/Output
  integer, dimension(MDIM), intent(in) :: cellID
  real, pointer :: solnVec(:,:,:,:)
  real, intent(in) :: energy
  integer, intent(in) :: iSpecies
  real, intent(out) :: ks

  ! aux variables
  real :: rho, temp, Ye
  real :: tmp(1)

  rho = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  Ye = 0 ! not used

  call calc_sca_opac_vectorized((/rho/), (/temp/), (/energy/), (/Ye/), iSpecies, tmp)
  ks = tmp(1)
end subroutine calc_sca_opac


!===================================================!
! amount of energy to emit in the space-energy cell !
!===================================================!
function calc_emit_dE(cellID, solnVec, g, dt, dV, iSpecies)
  use Driver_interface,  ONLY : Driver_abortFlash
  use Particles_data, only : pt_is_grey
  use Particles_data, only : sigma

  implicit none

#include "Flash.h"

  integer, dimension(MDIM), intent(in) :: cellID
  integer, intent(in) :: g, iSpecies ! energy group index
  real, intent(in) :: dt, dV
  real, pointer :: solnVec(:,:,:,:)
  real, parameter :: eps_dummy = 1.0d-11
  real :: dE, calc_emit_dE, ka, temp

  if (pt_is_grey) then
    temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
    call calc_abs_opac(cellID, solnVec, eps_dummy, iSpecies, ka)
    ! ka here is already in co-moving frame
    dE = 4.0 * ka * sigma * (temp**4.0) * dV * dt
  else
    call Driver_abortFlash("thermal_emission: non-grey thermal emission&
                            not yet implemented.")
  end if
end function calc_emit_dE

!======================================!
! flux coming from an emission surface !
!======================================!
function calc_face_flux(iSpecies)
  use Particles_data, only : pt_is_FacePlanck, pt_is_grey, pt_constFaceFlux, pt_FacePlanckTemp
  use Driver_interface,  ONLY : Driver_abortFlash
  use Particles_data, only : sigma

  implicit none

  integer, intent(in) :: iSpecies
  real :: calc_face_flux, flux

  if(.not. pt_is_grey) then
    call Driver_abortFlash("face emission: not yet implemented.")
  end if

  if (pt_is_FacePlanck) then
    flux = sigma * (pt_FacePlanckTemp**4)
  else
    flux = pt_constFaceFlux
  end if

  calc_face_flux = flux

end function calc_face_flux

!=================!
! Planck function !
!=================!
! written in a way that int f dx = 1
elemental function f_planck(x)
  implicit none
#include "constants.h"
  real, intent(in) :: x
  real :: f_planck

  if (x <= 20.0) then
    f_planck = (15.0 * (x**3)) / ((PI ** 4.0) * (EXP(x) - 1.0))
  else
    f_planck = 0.0d0
  end if

  return

end function f_planck

!===========================================================!
! derivative of equilibrium energy density with temperature !
!===========================================================!
subroutine calc_equilibrium_Edens_derivatives(cellID, solnVec, iSpecies, f0, beta, zeta)
  use Particles_data, only : a_rad, R, kB
  use Eos_data, only : eos_singleSpeciesA
  use Particles_data, only : pt_marshak_eos, pt_is_grey
  use Particles_data, only : rt_Emid

  implicit none

  integer, dimension(MDIM), intent(in) :: cellID
  integer, intent(in) :: iSpecies ! energy group index
  real, pointer :: solnVec(:,:,:,:)
  real :: rho, temp, gamc
  real, intent(out) :: beta, zeta
  real :: dE_dtemp_radiation, dE_dtemp_matter, c_V, x
  real, dimension(size(rt_Emid)) :: f0
  integer :: g

  if(.not. pt_is_grey) then
    call Driver_abortFlash("calc_equilibrium_Edens_derivatives: spectral method not implemented for photons")
  end if

  ! get grid variables
  temp = solnVec(TEMP_VAR  , cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  rho  = solnVec(DENS_VAR  , cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  gamc = solnVec(GAMC_VAR  , cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  
  c_V = R / ((gamc - 1.0d0)*eos_singleSpeciesA)  ! erg/K/mol * 1/amu = erg/K/g
  dE_dtemp_matter    = rho * c_V                 ! erg/K/cm^3
  dE_dtemp_radiation = 4.0d0 * a_rad * (temp**3) ! erg/K/cm^3
  beta =  dE_dtemp_radiation / dE_dtemp_matter

  ! quantities not relevant for photons (relates to rate of change of energy with electron fraction)
  zeta = 0

  ! blackbody function
  do g=1, size(rt_Emid)
    x = rt_Emid(g) / (kB * temp)
    f0(g) = f_planck(x) * PI**4 / (15.0 * x**3)
  end do

  ! Override when Marshak EOS is turned on
  if (pt_marshak_eos) then
    beta = 1.0 !0.25
  end if

end subroutine calc_equilibrium_Edens_derivatives


!=========================!
! sample_energy_blackbody !
!=========================!
subroutine sample_energy_blackbody(solnData, cellID, iSpecies, eps)
  use Particles_data, only : pt_is_grey, pt_grey_eps,&
                             pt_energy_min_eV, pt_energy_max_eV
  use Particles_data, only : kB, ev2erg
  use random, only : rand

  implicit none
#include "constants.h"

  ! Input/output
  real, pointer :: solnData(:,:,:,:)
  integer, dimension(MDIM) :: cellID
  integer, intent(in) :: iSpecies
  real, intent(out) :: eps

  ! aux variables
  real :: temp, x, f_x
  real, parameter :: f_max = 0.218886

  temp = solnData(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  if (pt_is_grey) then
    eps = pt_grey_eps
  else ! Thermal emission
    do 
      x = pt_energy_min_eV + rand() * (pt_energy_max_eV - pt_energy_min_eV)
      x = (x * ev2erg) / (kB * temp)
      f_x = f_planck(x)
      if (f_max * rand() .LT. f_x) exit
    end do

    eps = (x * kB * temp) ! in erg
  end if

end subroutine sample_energy_blackbody


function electron_number(iSpecies)
  implicit none

  integer, intent(in) :: iSpecies
  real :: electron_number

  electron_number = 0.0

end function electron_number

end module opacity
