!!****if* source/physics/Eos/EosMain/Multigamma/Eos
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
!!
!!  Eos
!!
!! SYNOPSIS
!!
!!  call Eos(integer(IN) :: mode,
!!                  integer(IN) :: vecLen,
!!                  real(INOUT) :: eosData(vecLen*EOS_NUM),
!!        optional, integer(IN) :: vecBegin,
!!        optional, integer(IN) :: vecEnd,
!!        optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!        optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!
!! DESCRIPTION
!!
!!
!!  This routine applies the equation of state to thermodynamic quantities
!!  at one or more grid points.  The number of points is determined by the
!!  argument veclen.  Data should be packaged for this routine in the 1d
!!  array, eosData.  The data in eosData is organized as: 1:vecLen
!!  points contain the first variable, vecLen+1:2*vecLen points contain
!!  the second variable and so on. The number and 
!!  order of variables in the array is determined by the constants defined
!!  in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and internal energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  internal energy and temperature are generated as output. Note that
!!  internal energy is EINT_VAR, not ENER_VAR.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!  
!!  This is a multigamma version, which means there are multiple species
!!  in the fluid, each in different abundances, and each with
!!  a different gamma.  This EOS takes into account the contribution to
!!  the thermodynamic properties of the gas from each species
!!  appropriately.
!!  
!!  The argument, massFrac, holds the mass fractions in an order determined
!!  by the Multispecies unit.  
!!  
!!  
!! ARGUMENTS 
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points (cells) for which the eosData array is sized.
!!             If vecBegin and vecEnd are not present, this is also the
!!             number of points (cells) for which EOS computation is to be done.
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on. 
!!
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is 1.
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
!!
!!  massFrac : Contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!             An implementation that does not need derivative quantities should
!!             set the mask equal to .false.
!!
!!
!!
!! NOTES
!!
!!  NSPECIES is defined in Simulation.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, and MODE_DENS_PRES are defined in constants.h.
!!
!!  User code should not call this implementation routine directly, but
!!  should call Eos and make sure that the desired Multigamma implementation
!!  is included in the simulation configuration.
!!  All code calling the Eos interface should include a 
!!    use Eos_interface 
!!  statement, preferable with "ONLY" attribute, e.g.,
!!    use Eos_interface, ONLY:  Eos
!!  All routines calling this routine directly should include a 
!!    use eos_localInterface
!!  statement, preferable with "ONLY" attribute, e.g.,
!!    use eos_localInterface, ONLY:  Eos
!!
!!  For Gamma and Multigamma routines, the entropy and entropy derivatives 
!!  calculations have not been confirmed to be correct.  Use with caution.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***

#define ORIGINAL_GAMC_AVERAGE

subroutine Eos(mode, vecLen, eosData,massFrac, mask ,vecBegin,vecEnd, diag)

!==============================================================================
  use eos_mgammaData, ONLY: eos_gammam1j,  eos_ggprodj, eos_gc
  use Eos_data, ONLY : eos_gasConstant
  use Driver_interface, ONLY : Driver_abort
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Simulation.h"
#include "Multispecies.h"

  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  real, optional, INTENT(in),dimension(vecLen*NSPECIES)    :: massFrac
  logical,  optional, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask
  logical, optional, INTENT(in) :: diag
  real,dimension(NSPECIES) :: weight
  real :: rt,abarValue, abarInv, zbarValue, zbarFrac
  real :: gmc ! avg these instead of rt - KW
  integer :: specieStart, specieEnd
  integer :: dens, temp, pres, eint, abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, c_v, c_p, gamc, pel, ne, eta
  integer :: i, ilo,ihi, rowLen

  
!==============================================================================
  if (present(vecBegin)) then
     ilo = vecBegin
  else
     ilo = 1
  end if
  if (present(vecEnd)) then
     ihi = vecEnd
  else
     ihi = vecLen
  end if
  rowLen = ihi - ilo + 1
#ifdef DEBUG_EOS
  if (ilo < 1 .OR. ilo > vecLen) then
     print*,'[Eos] ilo is',ilo
     call Driver_abort("[Eos] invalid ilo")
  end if
  if (ihi < 1 .OR. ihi > vecLen) then
     print*,'[Eos] ihi is',ihi
     call Driver_abort("[Eos] invalid ihi")
  end if
  if (rowLen < 0 .OR. rowLen > vecLen) then
     print*,'[Eos] rowLen is',rowLen
     call Driver_abort("[Eos] invalid rowLen")
  end if
#endif
  if (rowLen == 0) then
     print*,'[Eos] rowLen is 0.'
  end if

  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  entr = (EOS_ENTR-1)*vecLen

  if(.not.present(massFrac)) then
     call Driver_abort("[Eos] Multigamma needs mass fractions")
  end if


  ! First some computations that are common to all cases handled below. - KW
  do i = ilo,ihi
     specieStart = (i-1)*NSPECIES + 1
     specieEnd = i*NSPECIES

     call Multispecies_getSumInv(A, abarInv ,massFrac(specieStart:specieEnd))
     abarValue = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z,zbarFrac,massFrac(specieStart:specieEnd))
     zbarValue = abarValue*zbarFrac

#ifdef ORIGINAL_GAMC_AVERAGE
     weight = massFrac(specieStart:specieEnd)*eos_gammam1j
     call Multispecies_getSumInv(A,rt, weight)
#else
     weight = massFrac(specieStart:specieEnd)*eos_gc
     call Multispecies_getSumInv(A, gmc, weight)
#endif

     eosData(abar+i) = abarValue
     eosData(zbar+i) = zbarValue
#ifdef ORIGINAL_GAMC_AVERAGE
     eosData(gamc+i) = 1.0e0 + 1.0e0/(rt*eosData(abar+i))
#else
     eosData(gamc+i) = gmc*eosData(abar+i)
#endif
  end do


! mode:  temperature and density given
  if (mode == MODE_DENS_TEMP) then
     do i = ilo,ihi
        eosData(pres+i) = eos_gasConstant*eosData(dens+i)*eosData(temp+i)/eosData(abar+i)
        eosData(eint+i) = eos_gasConstant*eosData(temp+i)/((eosData(gamc+i)-1.e0)*eosData(abar+i))
        eosData(entr+i) = (eosData(pres+i)/eosData(dens+i) + eosData(eint+i))/eosData(temp+i)
     end do

! mode:  density and internal energy
  elseif (mode == MODE_DENS_EI) then

     do i = ilo,ihi 
        eosData(pres+i) = eosData(dens+i)*eosData(eint+i)*(eosData(gamc+i)-1.e0)
        eosData(temp+i) = eosData(eint+i)*(eosData(gamc+i)-1.e0) * eosData(abar+i)/eos_gasConstant
        eosData(entr+i) = (eosData(pres+i)/eosData(dens+i) + eosData(eint+i))/eosData(temp+i)
     end do

! mode:  density and pressure
  elseif (mode == MODE_DENS_PRES) then
     do i = ilo,ihi
        eosData(eint+i) = eosData(pres+i) / ( ( eosData(gamc+i) - 1.0 ) * eosData(dens+i) )
        eosData(temp+i) = eosData(eint+i) *(eosData(gamc+i) - 1.0) * eosData(abar+i) / eos_gasConstant
        eosData(entr+i) = (eosData(pres+i)/eosData(dens+i) + eosData(eint+i))/eosData(temp+i)
     end do

  else 
     call Driver_abort("[Eos] Unrecognized input mode given to Eos")     
  endif
  
  if(present(mask)) then
     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        !! flash2 equation is dpt = gasconstant*density/abar
        eosData(dpt+ilo:dpt+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi)/eosData(abar+ilo:abar+ihi)
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
        !! flash2 equation is dpd =gasconstant*temperature/abar 
        eosData(dpd+ilo:dpd+ihi) = eos_gasConstant*eosData(temp+ilo:temp+ihi)/eosData(abar+ilo:abar+ihi)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
         eosData(det+ilo:det+ihi) = eos_gasConstant / &
                                    ((eosData(gamc+ilo:gamc+ihi) - 1.0) * eosData(abar+ilo:abar+ihi))
     end if
     if(mask(EOS_DED))then 
        ded = (EOS_DED-1)*vecLen
        eosData(ded+ilo:ded+ihi) = 0.
     end if

     ! Entropy derivatives
     if (mask(EOS_DST)) then
        if (mask(EOS_DET) .AND. mask(EOS_DPT)) then
           det = (EOS_DET-1)*vecLen
           dpt = (EOS_DPT-1)*vecLen
           dst = (EOS_DST-1)*vecLen
           eosData(dst+ilo:dst+ihi) = &
         &       ((eosData(dpt+ilo:dpt+ihi) / eosData(dens+ilo:dens+ihi) + eosData(det+ilo:det+ihi)) - &
         &       (eosData(pres+ilo:pres+ihi) / eosData(dens+ilo:dens+ihi) + eosData(eint+ilo:eint+ihi)) / &
         &       eosData(temp+ilo:temp+ihi) ) / eosData(temp+ilo:temp+ihi)
        else
           call Driver_abort("[Eos] Cannot calculate EOS_DST without EOS_DET and EOS_DPT")
        end if
     end if
     if (mask(EOS_DSD)) then
        if (mask(EOS_DED) .AND. mask(EOS_DPD)) then
           dsd = (EOS_DSD-1)*vecLen
           ded = (EOS_DED-1)*vecLen
           dpd = (EOS_DPD-1)*vecLen
           eosData(dsd+ilo:dsd+ihi) = &
        &       ( ((eosData(dpd+ilo:dpd+ihi) - eosData(pres+ilo:pres+ihi)/eosData(dens+ilo:dens+ihi)) / &
        &          eosData(dens+ilo:dens+ihi)) + eosData(ded+ilo:ded+ihi)) / eosData(temp+ilo:temp+ihi)
         else
           call Driver_abort("[Eos] Cannot calculate EOS_DSD without EOS_DED and EOS_DPD")
        end if
     end if


     if(mask(EOS_PEL))then 
        pel = (EOS_PEL-1)*vecLen
        eosData(pel+ilo:pel+ihi) = 0.
     end if
     if(mask(EOS_NE))then 
        ne = (EOS_NE-1)*vecLen
        eosData(ne+ilo:ne+ihi) = 0.
     end if
     if(mask(EOS_ETA))then 
        eta = (EOS_ETA-1)*vecLen
        eosData(eta+ilo:eta+ihi) = 0.
     end if
     
     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           eosData(c_v+ilo:c_v+ihi) = eosData(det+ilo:det+ihi)
        else
           call Driver_abort("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if
     
     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           eosData(c_p+ilo:c_p+ihi) = eosData(gamc+ilo:gamc+ihi)*eosData(c_v+ilo:c_v+ihi)
        else
           call Driver_abort("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
     end if
  end if
  return
end subroutine Eos

!!..no matter what the input mode compute the entropy
!!..ignore the -chemical_potential*number_density part for now
!!$  dens_inv = 1.0e0/eosData(dens+ilo:+ihi)
!!$  temp_inv = 1.0e0/eosData(temp+ilo:+ihi)
!!$  stot     = (pres*dens_inv + eosData(eint+ilo:+ihi))*temp_inv 
!!$  dstotdd  = (eosData(EOS_DPD)*dens_inv - pres*dens_inv*dens_inv + eosData(EOS_DED))*temp_inv
!!$  dstotdt  = (eosData(EOS_DPT)*dens_inv + eosData(EOS_DET))*temp_inv  - (pres*dens_inv + eosData(eint+ilo:+ihi)) * temp_inv*temp_inv 
!!$  



