!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity
!!
!! NAME
!!  MagneticResistivity
!!
!! SYNOPSIS
!!  call MagneticResistivity  (real(in)  :: dens,
!!                            real(in)   :: temp,
!!                            real (in)  :: xn(NSPECIES) ,
!!                            real(out)  :: magResist)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron Magnetic Resistivity for all materials,
!! including those with Z > 1. The expressions used here comes from NRL plasma form.
!!
!!  Returns Magnetic Resistivity resPar in magResist.
!!
!! ARGUMENTS
!!
!! ARGUMENTS
!!  temp      - Plasma temperature
!!  dens      - Plasma density
!!  xn        - Species
!!  magResist - Magnetic resistivity
!!
!!***

#include "Flash.h"
#include "constants.h"  

subroutine MagneticResistivity(temp,dens,xn,magResist)
  use MagneticResistivity_data, ONLY: mag_useMagneticResistivity, &
       res_mele, res_qele, res_navo, res_speedlt, res_boltz, res_hbar
  use MagneticResistivity_data, ONLY: res_mUnit
  use MagneticResistivity_data, ONLY: res_coef
  use MagneticResistivity_data, ONLY: res_maxRes
  use Eos_interface, ONLY: Eos_getAbarZbar

  implicit none

  real, intent(IN) :: temp, dens
  real, intent(IN), dimension(NSPECIES) :: xn
  real, intent(OUT):: magResist

  real :: resPar
  real :: tele, tele_eV
  real :: tion
  real :: nele
  real :: nion
  real :: abar
  real :: zbar
  real :: eqtime
  real :: resPerpLoc
  real :: ll

  call Eos_getAbarZbar(massFrac=xn,abar=abar,zbar=zbar)

  nion = dens * res_navo / abar
  nele = zbar * nion

  tele = temp
  tion = temp

  tele_eV = tele/11604.5221
  call loglambda(tele, zbar*nion, zbar, ll)

  resPerpLoc = 8.196e5*zbar*ll/tele_eV**1.5  !! In CGS -- Here tele has to be in eV

  ! This formula is only valid when the magnetic field is strong (and
  ! only for hydrogen):
  resPar = resPerpLoc/1.96

  !! FOR CGS--> SI unit conversion
  if (res_mUnit == "SI" .or. res_mUnit == "si" ) then
     resPerpLoc = resPerpLoc*1.0e-4
     resPar  = resPar *1.0e-4
  end if

  !! cap resPar and resPerp to avoid unphysically cold regions
  resPar = min(resPar, res_maxRes)
  
  magResist = resPar

contains

  subroutine loglambda(tele, nele, zbar, ll)
    implicit none

    ! This subroutine computes the Coulomb logarithm. The formula used
    ! comes from Atzeni.
    
    real, intent(in)  :: tele ! electron temperature [K]
    real, intent(in)  :: nele ! electron number density [cm^-3]
    real, intent(in)  :: zbar ! the average ionization [unitless]
    real, intent(out) :: ll   ! the coulomb logarithm [unitless]

    real :: bmax, bmin, bmin_classic, bmin_quantum
    real, parameter :: ll_floor = 1.0
    
    bmax = sqrt(res_boltz * tele / (4*PI * res_qele**2 * nele))
    
    bmin_classic = zbar * res_qele**2 / (3*res_boltz*tele)
    bmin_quantum = res_hbar / (2*sqrt(3*res_boltz*tele*res_mele))
    bmin = max(bmin_classic, bmin_quantum)

    ll = max(log(bmax/bmin), ll_floor)

  end subroutine loglambda
  
end subroutine MagneticResistivity
