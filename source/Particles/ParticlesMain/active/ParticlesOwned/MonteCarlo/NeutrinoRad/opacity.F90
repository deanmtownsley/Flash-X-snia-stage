module opacity 


  use UnitsModule, ONLY : Gram, Centimeter, Second, AtomicMassUnit, Kelvin, MeV

  implicit none

  real, parameter :: conv_x    = Centimeter
  real, parameter :: conv_dens = Gram / Centimeter**3
  real, parameter :: conv_mom  = Gram / Centimeter**2 / Second ! momentum density
  real, parameter :: conv_enr  = Gram / Centimeter / Second**2 ! energy density
  real, parameter :: erg_to_MeV = (Gram*Centimeter**2/Second**2) / MeV
  real, parameter :: conv_ne   = Gram / Centimeter**3 / AtomicMassUnit
  real, parameter :: conv_J    = Gram/Second**2/Centimeter
  real, parameter :: conv_H    = Gram/Second**3


  contains

!====================!
! ABSORPTION OPACITY !
!====================!
subroutine calc_abs_opac_vectorized(rho, temp, Ye, energy, iSpecies, ka)
  use Particles_data, only : pt_is_grey 
  use Particles_data, only : pt_dens_threshold, pt_grey_abs_opac
  use NeutrinoOpacitiesComputationModule
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, intent(in) :: rho(:), temp(:), energy(:), Ye(:)
  integer, intent(in) :: iSpecies
  real, intent(out) :: ka(:)
  integer :: ip_B, ip_E

  if (pt_is_grey) then
    ka = pt_grey_abs_opac
    ka = ka * rho
  else
    ip_B = 1
    ip_E = size(rho)
    call ComputeNeutrinoOpacities_EC_Vector(iP_B, iP_E, energy * erg_to_MeV * MeV, rho * conv_dens, temp * Kelvin, Ye, iSpecies, ka)
    ka = ka * conv_x
  end if
  
  where (rho .lt. pt_dens_threshold)
    ka = 0.0d0
  end where
  
end subroutine calc_abs_opac_vectorized

subroutine calc_abs_opac(cellID, solnVec, energy, iSpecies, ka)
  use Driver_interface,  ONLY : Driver_abortFlash
  use NeutrinoOpacitiesComputationModule
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
  Ye = solnVec(YE_MSCALAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  call calc_abs_opac_vectorized( (/rho/), (/temp/), (/Ye/), (/energy/), iSpecies, tmp)
  ka = tmp(1)
end subroutine calc_abs_opac

!====================!
! SCATTERING OPACITY !
!====================!
subroutine calc_sca_opac_vectorized(rho, temp, Ye, energy, iSpecies, ks)
  use Particles_data, only : pt_dens_threshold, pt_grey_sca_opac, pt_is_grey
  use NeutrinoOpacitiesComputationModule
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Input/Output
  real, intent(in) :: rho(:), temp(:), energy(:), Ye(:)
  integer, intent(in) :: iSpecies
  real, intent(out) :: ks(:)
  integer :: ip_B, ip_E, imoment

  if (pt_is_grey) then
    ks = pt_grey_sca_opac
    ks = ks * rho
  else
    iP_B = 1
    iP_E = size(rho)
    imoment = 1
    call ComputeNeutrinoOpacities_ES_Vector(iP_B, iP_E, energy * erg_to_MeV * MeV, rho * conv_dens, temp * Kelvin, Ye, iSpecies, imoment, ks)
    ks = ks * conv_x
  end if

  where (rho .lt. pt_dens_threshold)
    ks = 0.0d0
  end where
end subroutine calc_sca_opac_vectorized


subroutine calc_sca_opac(cellID, solnVec, energy, iSpecies, ks)
  use Driver_interface,  ONLY : Driver_abortFlash
  use NeutrinoOpacitiesComputationModule
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
  Ye = solnVec(YE_MSCALAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  call calc_sca_opac_vectorized( (/rho/), (/temp/), (/Ye/), (/energy/), iSpecies, tmp)
  ks = tmp(1)
end subroutine calc_sca_opac

!===================================================!
! amount of energy to emit in the space-energy cell !
!===================================================!
function calc_emit_dE(cellID, solnVec, g, dt, dV, iSpecies)
  use Particles_data, only : pt_is_grey
  use Particles_data, only : clight, h_planck
  use Particles_data, only : rt_Emid, rt_int4piE3dE
  use NeutrinoOpacitiesComputationModule

  implicit none

#include "Flash.h"

  integer, dimension(MDIM), intent(in) :: cellID
  integer, intent(in) :: g, iSpecies ! energy group index
  real, intent(in) :: dt, dV
  real, pointer :: solnVec(:,:,:,:)
  real :: dE, calc_emit_dE, ka, rho, temp, ye, feq
  integer :: iE_B, iE_E

  if (pt_is_grey) then
    call Driver_abortFlash("thermal_emission: grey thermal emission&
                            not yet implemented.")
  else
    rho  = solnVec(DENS_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
    temp = solnVec(TEMP_VAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
    Ye = solnVec(YE_MSCALAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

    ! compute the fermi-dirac function with thornado
    call ComputeEquilibriumDistributions_Point(rt_Emid(g) * erg_to_MeV * MeV, rho * conv_dens, temp * Kelvin, Ye, feq, iSpecies)

    ! get the absorption opacity
    call calc_abs_opac(cellID, solnVec, rt_Emid(g), iSpecies, ka)

    ! compute the emitted energy (erg)
    dE = feq * ka * clight*dt * dV * rt_int4piE3dE(g) / (h_planck*clight)**3
  end if

  calc_emit_dE = dE
end function calc_emit_dE

!======================================!
! flux coming from an emission surface !
!======================================!
function calc_face_flux(iSpecies)
  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

  integer, intent(in) :: iSpecies
  real :: calc_face_flux

  call Driver_abortFlash("face_emission: not yet implemented.")

end function calc_face_flux

!===========================================================!
! derivative of equilibrium energy density with temperature !
!===========================================================!
subroutine calc_equilibrium_Edens_derivatives(cellID, solnVec, iSpecies, f0, beta, zeta)
  use Particles_data, only : rt_Emid, rt_int4piE3dE, rt_int4piE2dE
  use Particles_data, only : h_planck, clight
  use EquationOfStateModule_TABLE, only : BaryonMass
  use UnitsModule, only : gram
  use NeutrinoOpacitiesComputationModule, only : ComputeEquilibriumDistributionAndDerivatives_Point

  implicit none

  integer, dimension(MDIM), intent(in) :: cellID
  integer, intent(in) :: iSpecies ! energy group index
  real, pointer :: solnVec(:,:,:,:)
  real, intent(out) :: beta, zeta
  real, dimension(size(rt_Emid)) :: f0, df0dY, df0dU

  real :: rho, temp, ye

  ! get grid variables
  temp = solnVec(TEMP_VAR  , cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  rho  = solnVec(DENS_VAR  , cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  ye   = solnVec(YE_MSCALAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  ! returns U in units of erg/gram
  ! (see EquationOfStateModule_TABLE.F90:266,1230 and NeutrinoOpacitiesComputationModule.F90:695
  ! E_thornado = E_MeV*UnitE where UnitE=MeV=1.32e-57
  ! rho_thornado = rho_cgs*UnitD where UnitD=g/ccm=7.42e-25
  ! T_thornado = T_K*UnitT where UnitT=K=1.14e-67
  ! U_thornado = U_cgs*UnitU where UnitU=erg/g=1.11e-21
  ! so 1/U_cgs = 1/U_thornado * UnitU
  call ComputeEquilibriumDistributionAndDerivatives_Point(1, size(rt_Emid), rt_Emid * erg_to_MeV * MeV, rho * conv_dens, temp * Kelvin, Ye, iSpecies, f0, df0dY, df0dU)
  df0dU = df0dU * (conv_enr/conv_dens) / rho                                         ! convert to 1/(erg cm^-3)
  beta = sum(df0dU * rt_int4piE3dE / (h_planck*clight)**3)                           ! dimensionless
  zeta = sum(df0dY * rt_int4piE3dE / (h_planck*clight)**3) / (rho*gram / BaryonMass) ! erg

end subroutine calc_equilibrium_Edens_derivatives

!=========================!
! sample_energy_blackbody !
!=========================!
subroutine sample_energy_blackbody(solnData, cellID, iSpecies, eps)
  use Particles_data, only : pt_is_grey, pt_grey_eps
  use Particles_data, only : rt_Emid, rt_int4piE3dE
  use Particles_data, only : kB, ev2erg
  use random, only : rand
  use NeutrinoOpacitiesComputationModule

  implicit none
#include "constants.h"

  ! Input/output
  real, pointer :: solnData(:,:,:,:)
  integer, dimension(MDIM) :: cellID
  integer, intent(in) :: iSpecies
  real, intent(out) :: eps

  ! aux variables
  real :: rho, temp, Ye, x, f_x
  real, dimension(size(rt_Emid)) :: blackbody
  integer :: iE_B, iE_E, g

  iE_B = 1
  iE_E = size(rt_Emid)

  rho  = solnData(DENS_VAR  , cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  temp = solnData(TEMP_VAR  , cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))
  Ye   = solnData(YE_MSCALAR, cellID(IAXIS), cellID(JAXIS), cellID(KAXIS))

  if (pt_is_grey) then
    eps = pt_grey_eps
  else ! Thermal emission
    ! compute the fermi-dirac function with thornado
    ! results are in weird code units but we don't care - they will be normalized out.
    call ComputeEquilibriumDistributions_Point(iE_B, iE_E, rt_Emid * erg_to_MeV * MeV, rho * conv_dens, temp * Kelvin, Ye, blackbody, iSpecies)
    
    ! integrate over energy to get the integrated emissivity in each bin
    ! assume distribution function and opacity are constant in the bin
    ! again, we don't care about constant factors.
    blackbody = blackbody * rt_int4piE3dE

    ! divide by the max value
    blackbody = blackbody / MAXVAL(blackbody)

    ! get energy bin by rejection sampling
    do 
      g = int(size(rt_Emid)*rand()) + 1;
      if (rand() .LT. blackbody(g)) exit
    end do

    eps = rt_Emid(g)
  end if

end subroutine sample_energy_blackbody

function electron_number(iSpecies)
  implicit none

  integer, intent(in) :: iSpecies
  real :: electron_number

  electron_number = 0.0
  ! set neutrino flavor
  if(iSpecies .eq. 1) then
    electron_number = 1.0
  else if(iSpecies .eq. 2) then
    electron_number = -1.0
  end if
end function electron_number

end module opacity
