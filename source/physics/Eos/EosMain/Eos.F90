!!****if* source/physics/Eos/EosMain/Eos
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!
!!  Eos
!!
!! SYNOPSIS
!!
!!  call      Eos(integer(IN) :: mode,
!!                integer(IN) :: vecLen,
!!                real(INOUT) :: eosData(vecLen*EOS_NUM),
!!      optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional, logical(IN),target :: mask(EOS_VARS+1:EOS_NUM),
!!      optional, integer(IN) :: vecBegin,
!!      optional, integer(IN) :: vecEnd  ,
!!      optional, integer(IN) :: diagFlag  )
!!
!! DESCRIPTION
!!
!!  This routine applies the equation of state to thermodynamic 
!!  quantities at one or more grid cells.  The number of cells is 
!!  determined by the argument veclen.  The routine expects data packaged
!!  for it in the 1d array, eosData.  The data in eosData is organized as:
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
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
!!
!!  ARGUMENTS
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
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells for optimization,
!!             untested.
!!             Only supported by some implementations: currently, supported in
!!             eos_idealGamma and eos_mgamma.
!!             If not present, the default is 1.
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells for optimization,
!!             untested.
!!             Only supported by some implementations: currently, supported in
!!             eos_idealGamma and eos_mgamma.
!!             If not present, the default is vecLen.
!!
!!  diagFlag : This optional dummy can be used by some implementations to pass error
!!             indications or other diagnostic information back to the caller.
!!             If the implementation has notrhing to report, 0 is returned.
!!             The meanoing of other values has not yet been standardized.
!!             Only supported by some implementations: currently, only some
!!             implementations of eos_idealGamma return any nonzero diagnostics.
!!             
!!
!!
!! EXAMPLE
!!
!!
!! ------------------ Row at a time example, with derivates (based on Eos_unitTest) --------
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Simulation.h"       ! for NSPECIES, EOS_NUM
!!  #include "Eos.h"         ! for EOS_VAR order
!!  integer veclen, isize, jsize, ksize, i,j,k, e
!!  real, dimension(:), allocatable :: eosData
!!  real, dimension(:), allocatable :: massFrac
!!  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
!!  real, allocatable, dimension(:,:,:,:) :: derivedVariables
!!  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!!
!!   ! in the Eos_unitTest, this loops over all blocks.... here is a snippet from inside
!!    allocate(derivedVariables( .....)
!!    allocate(eosData(vecLen*EOS_NUM))
!!    allocate(massFrac(vecLen*NSPECIES))
!!    mask = .true.
!!
!!    ! indices into the first location for these variables
!!    pres = (EOS_PRES-1)*vecLen
!!    dens = (EOS_DENS-1)*vecLen
!!    temp = (EOS_TEMP-1)*vecLen
!!
!!
!!    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
!!           do i = 1,vecLen
!!              massFrac((i-1)*NSPECIES+1:i*NSPECIES) = &
!!                   solnData(SPECIES_BEGIN:SPECIES_END,ib+i-1,j,k)
!!           end do
!!
!!           eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
!!           eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
!!           ! Eos Helmholtz needs a good initial estimate of temperature no matter what the mode
!!           eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
!!
!!           call Eos(MODE_DENS_PRES,vecLen,eosData,massFrac,mask)
!!
!!           do e=EOS_VARS+1,EOS_NUM
!!              m = (e-1)*vecLen
!!              derivedVariables(1:vecLen,j-NGUARD,k-NGUARD,e) =  eosData(m+1:m+vecLen)
!!           end do
!!        end do
!!     end do
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
!!  All routines calling this routine should include a 
!!  use Eos_interface 
!!  statement, preferable with "ONLY" attribute e.g.
!!  use Eos_interface, ONLY:  Eos
!!
!!  When using Helmholtz, strange behaviour occurs.  See the notes in subunit
!!  Helmholtz/Eos.F90 directly.  For Helmholtz, when operating in MODE_DENS_EI,
!!  the INPUT energy is updated.
!!  Similarly, when operating in MODE_DENS_PRES, the INPUT pressure is updated. 
!!  Physicists need to be aware of this.
!!  This behavior can be turned off by setting the runtime parameter
!!  eos_forceConstantInput.
!!
!!  For Gamma and Multigamma routines, the entropy and entropy derivatives 
!!  calculations have not been confirmed to be correct.  Use with caution.
!!
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!  eos_idealGamma  actual implementation routine for the Gamma Eos implementation
!!  eos_mgamma  actual implementation routine for the Multigamma Eos implementation
!!
!!***

subroutine Eos(mode, vecLen, eosData, massFrac, mask, vecBegin,vecEnd, diagFlag)

!==============================================================================
  use Driver_interface, ONLY : Driver_abort
  use Eos_data, ONLY : eos_meshMe, eos_type
  use eos_localInterface, ONLY : eos_idealGamma, eos_mgamma, eos_helmholtz,&
      eos_tabulated, eos_nuclear, eos_weaklib
  implicit none
#include "constants.h"
#include "Eos.h"
#include "Simulation.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  integer, optional, INTENT(out)    :: diagFlag
  logical :: pMassFrac_and_mask, pMassFrac, pMask

  if (present(diagFlag)) diagFlag = 0

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *
  if (mode==MODE_EOS_WRAPPERONLY) return ! * Return immediately for MODE_EOS_WRAPPERONLY! *

  pMassFrac = present(massFrac)
  pMask = present(mask)
  pMassFrac_and_mask = pMassFrac.and.pMask
  
  if(pMassFrac_and_mask) then
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac, mask=mask, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac, mask=mask)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData, massFrac=massFrac, mask=mask)
     case(EOS_TAB)
        call eos_tabulated(mode, vecLen, eosData, massFrac=massFrac, mask=mask)
     case(EOS_NUC)
        call eos_nuclear(mode, vecLen, eosData, massFrac, mask=mask)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData, massFrac, mask)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  elseif (pMassFrac) then
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData, massFrac=massFrac)
     case(EOS_TAB)
        call eos_tabulated(mode, vecLen, eosData, massFrac=massFrac)
     case(EOS_NUC)
        call eos_nuclear(mode, vecLen, eosData, massFrac=massFrac)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData, massFrac)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  elseif (pMask) then
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, mask=mask, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, mask=mask)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData, mask=mask)
     case(EOS_TAB)
        call eos_tabulated(mode, vecLen, eosData, mask=mask)
     case(EOS_NUC)
        call eos_nuclear(mode, vecLen, eosData, mask=mask)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData, mask=mask)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  else
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData)
     case(EOS_TAB)
        call eos_tabulated(mode, vecLen, eosData)
     case(EOS_NUC)
        call eos_nuclear(mode, vecLen, eosData)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  end if
  return
end subroutine Eos
