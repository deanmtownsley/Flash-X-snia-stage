!!****if* source/physics/materialProperties/Conductivity/ConductivityAniso/Conductivity_anisoFullState
!!
!! NAME
!!  Conductivity_anisoFullState
!!
!! SYNOPSIS
!!  call Conductivity_anisoFullState(real(in)    :: solnVec(NUNK_VARS),
!!                     OPTIONAL,real(out)   :: isochoricCond(3),
!!                     OPTIONAL,real(out)   :: diffCoeff(3),
!!                     OPTIONAL,integer(in) :: component)
!!
!! DESCRIPTION
!!
!! Computes the Epperlein-Haines electron conductivities for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from Epperlein & Haines, Phys. Fluids 29, 1029, 1986.
!!
!!  Returns anisotropic thermal conductivities and/or diffusivity coefficients.
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   isochoricCond  :   isochoric conductivities (parallel, perpendicular, cross)
!!   diffCoeff :   diffusion coefficients ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!***

#include "Simulation.h"
#include "constants.h"  


subroutine Conductivity_anisoFullState(solnVec,isochoricCond,diffCoeff,component)
  use cond_commonData, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele
  use cond_commonData, ONLY: cond_navo, cond_c
  use Eos_interface, ONLY: Eos_getAbarZbar
  use cond_interface,ONLY: cond_getCv
  use Driver_interface,ONLY: Driver_abort

  implicit none

  real,    target,   intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff(3)
  real,    OPTIONAL, intent(OUT)  :: isochoricCond(3)
  integer, OPTIONAL, intent(IN) :: component

  real, pointer :: massfrac(:)

  real :: parCondLoc, perpCondLoc, crossCondLoc, isochoricCondLoc(3), diffCoeffLoc(3)
  integer :: componentLoc, tempToUse

  real :: xtemp, xden, xmag
  real :: nele, nion, ll
  real :: abar, zbar
  real :: cv
  real :: mion
  real :: omega, tau, chi, delta
  real :: max_cond
  logical :: useTion, useTele
  integer :: tionVar, teleVar

  real :: g0, g0p, g1p, c0p, c1p, c2p, g0pp, g1pp, c0pp, c1pp, c2pp

  isochoricCondLoc(:) = 0.0
  diffCoeffLoc(:) = 0.0

#ifndef FLASH_3T
  call Driver_abort("[Conductivity_anisoFullState] Anisotropic conductivity only works in 3T")
#endif

#ifndef MAGP_VAR
     call Driver_abort("[Conductivity_anisoFullState] Anisotropic conductivity only works in MHD")
#endif

#if defined(DENS_VAR) && defined(TEMP_VAR)
  if (present(component)) then
     componentLoc = component
  else
     componentLoc = 0
  end if

  tempToUse = TEMP_VAR

  select case (componentLoc)
  case(1)
#ifdef TION_VAR
     tempToUse = TION_VAR
#endif
  case(2)
#ifdef TELE_VAR
     tempToUse = TELE_VAR
#endif
  case(3)
#ifdef TRAD_VAR
     call Driver_abort("[Conductivity_anisoFullState] Anisotropic conductivity does not work with radiation")
#endif
  end select

#ifdef TION_VAR
  useTion = (tempToUse == TION_VAR)
  tionVar = TION_VAR
#else
  useTion = .FALSE.
  tionVar = TEMP_VAR
#endif
#ifdef TELE_VAR
  useTele = (tempToUse == TELE_VAR)
  teleVar = TELE_VAR
#else
  useTele = .FALSE.
  teleVar = TEMP_VAR
#endif

    if (cond_useConductivity) then
       call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)
       xden = solnVec(DENS_VAR)
       nele = zbar * xden * cond_navo / abar
       xmag = sqrt(solnVec(MAGP_VAR)*8.0*PI)

       if(useTele .or. tempToUse == TEMP_VAR) then
          xtemp = solnVec(teleVar)

          call cond_getEppHainCoeffs(zbar,g0,g0p,g1p,c0p,c1p,c2p,&
                                     g0pp,g1pp,c0pp,c1pp,c2pp)

          ! Compute the electron conductivities:
          call cond_logLambda(xtemp, nele, solnVec(tionVar), zbar, ll)
          tau = 6.0*sqrt(2.0*cond_mele)*(PI*cond_boltz*xtemp)**1.5 / &
                (ll * cond_qele**4 * zbar**2 * nele)
          omega = cond_qele*xmag / (cond_mele*cond_c)
          chi = omega*tau

          parCondLoc = nele*cond_boltz**2*xtemp*tau / cond_mele * g0
          perpCondLoc = nele*cond_boltz**2*xtemp*tau / cond_mele * &
                        (g1p*chi + g0p) / (chi**3 + c2p*chi**2 + c1p*chi + c0p)
          crossCondLoc = nele*cond_boltz**2*xtemp*tau / cond_mele * &
                        chi*(g1pp*chi + g0pp) / (chi**3 + c2pp*chi**2 + c1pp*chi + c0pp)

          isochoricCondLoc(:) = (/ parCondLoc, perpCondLoc, crossCondLoc /)

          !! For testing purposes
          !! isochoricCondLoc = (/ 2.e13, 5.e12, 0.0 /)

          if (present(diffCoeff)) then
             cv = cond_getCv(component=2,solnVec=solnVec)
             diffCoeffLoc(:) = isochoricCondLoc(:)/(xden*cv)
          end if

       elseif(useTion) then
          mion = abar / cond_navo
          nion = nele / zbar
          xtemp = solnVec(tionVar)

          ! Compute the ion conductivities:
          call cond_logLambdaII(xtemp, solnVec(teleVar), nele, zbar, ll)
          tau = 12.0*sqrt(mion)*(PI*cond_boltz*xtemp)**1.5 / &
                (ll * cond_qele**4 * zbar**4 * nion)
          omega = zbar*cond_qele*xmag / (mion*cond_c)
          chi = omega*tau
          delta = chi**4 + 2.7*chi**2 + 0.677

          parCondLoc = nion*cond_boltz**2*xtemp*tau / mion * 3.906
          perpCondLoc = nion*cond_boltz**2*xtemp*tau / mion * &
                        (2*chi**2 + 2.645) / delta
          crossCondLoc = nion*cond_boltz**2*xtemp*tau / mion * &
                        chi*(2.5*chi**2 + 4.65) / delta

          isochoricCondLoc(:) = (/ parCondLoc, perpCondLoc, crossCondLoc /)

          !! For testing purposes
          !! isochoricCondLoc = (/ 2.e13, 5.e12, 0.0 /)

          if (present(diffCoeff)) then
             cv = cond_getCv(component=1,solnVec=solnVec)
             diffCoeffLoc(:) = isochoricCondLoc(:)/(xden*cv)
          end if
       end if
    end if

#endif

  if(present(isochoricCond)) isochoricCond(:) = isochoricCondLoc(:)
  if(present(diffCoeff)) diffCoeff(:) = diffCoeffLoc(:)

end subroutine Conductivity_anisoFullState
