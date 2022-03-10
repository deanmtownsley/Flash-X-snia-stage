!!****if* source/physics/Eos/EosMain/Helmholtz/SpeciesBased/starkiller/Eos
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
!! Eos
!!
!! SYNOPSIS
!!
!!  call Eos(integer(IN) :: mode,
!!                     integer(IN) :: vecLen,
!!                     real(INOUT) :: eosData(vecLen*EOS_NUM),
!!           optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!           optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!
!! DESCRIPTION
!!
!!   Driver for the Helmholtz equation of state (starkiller version).
!!   See the NOTES section for important information about this implementation.
!!
!!  This routine applies the equation of state to thermodynamic 
!!  quantities at one or more grid cells.  The number of cells is 
!!  determined by the argument veclen.  Data is packaged for this 
!!  routine in the 1d array, eosData.  The data in eosData is organized: 
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  energy and temperature are generated as output.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!
!!  The primary difference in the starkiller implementation is that if 
!!  vecLen > eos_vecLenACC, the calculations will be performed on the GPU
!!  through the use of OpenACC directives. This cannot be done if Eos is
!!  called from within an OpenMP parallel region.
!!
!!  ARGUMENTS
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
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
!!             Although this is declared as an optional dummy argument, an
!!             actual argument MUST be provided when calling THIS implementation
!!             of Eos.
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
!!             The Helmholtz EOS kernel calculation ignores the mask setting and calculates
!!             all derivatives, whether needed or not.  This routine does not return
!!             derivatives if the mask is requested, but the calculation is not speeded up
!!             by setting the mask.
!!
!!
!! PARAMETERS
!!
!!  eos_tol    Controls the accuracy of the Newton Rhapson iterations for MODE_DENS_EI and 
!!             MODE_DENS_PRES.
!!
!!
!!
!!*** 


subroutine Eos(mode,vecLen,eosData,massFrac,mask,vecB,vecE,diagFlag)

  use Driver_interface, ONLY : Driver_abort
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
       Multispecies_getSumFrac
  use Logfile_interface, ONLY:  Logfile_stampMessage

  use eos_helmInterface, ONLY : eos_helm

  use eos_helmData, ONLY: eos_tol, eos_maxNewton,&
       eos_forceConstantInput, eos_useStarkiller, eos_vecLenACC
  use Eos_data, ONLY : eos_smallt, eos_meshMe, eos_singleSpeciesA, eos_singleSpeciesZ
  use eos_vecData, ONLY:  tempRow, denRow, etotRow, abarRow, zbarRow, &
       gamcRow, ptotRow, deaRow, dezRow, stotRow, dsdRow, dstRow, &
       detRow, dptRow, dpdRow, dedRow, pelRow, neRow, etaRow, detatRow, cvRow, cpRow, NROWMAX

  use eos_type_module, ONLY : composition, eos_input_rt, eos_input_re, eos_input_rp, clean_state, &
       eos_state, reduced_eos_state, eos_tp, comp_tp
  use actual_eos_module, ONLY : actual_eos, max_newton

  !$ use omp_lib
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
#endif

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  integer, optional, INTENT(in) :: vecB,vecE
  real, INTENT(inout), dimension(vecLen*EOS_NUM) :: eosData
  real, optional,INTENT(in), dimension(vecLen*NSPECIES) :: massFrac
  ! must correspond to dimensions of Eos_wrapped
  logical,optional,target, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
  integer, optional, INTENT(out) :: diagFlag

  integer :: i, k
  integer :: pres, temp, dens, gamc, eint
  integer :: abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, dea, dez, pel, ne, eta, detat, c_v, c_p
  real    :: abarInv, zbarFrac

  ! declare some local storage for the results of the Newton iteration
  real,dimension(vecLen)::  ewantRow, tnew, error,pwantRow
  !  local storage for forcedConstantInput -- could be allocatable, but might be so slow
  !  that it's not worth the small storage save.
  real,dimension(vecLen)::  psaveRow, esaveRow

  integer :: eos_input
  !$acc declare create(eos_input)

  type(eos_tp) :: local_eos_state(vecLen)
  type(comp_tp) :: comp_state
  integer :: vecBegin, vecEnd
  logical :: in_par, do_acc

  !      Fill the pipe with the initial temperature, density, and composition.
  !      The composition is parametrized by abar and zbar, which are computed
  !      from the mass fractions xn().

  !  if you are working with electron abundance mass scalars, then you don't
  !  necessarily have to have mass fractions.
  if(.not.present(massFrac)) then
     call Driver_abort("[Eos] Helmholtz needs mass fractions")
  end if

  if (present(diagFlag)) then
     diagFlag = 0
  endif

  vecBegin = 1
  vecEnd = vecLen

  ! These integers are indexes into the lowest location in UNK that contain the appropriate variable
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen   ! in flash2 eos_helm, this is etot
  gamc = (EOS_GAMC-1)*vecLen   ! in flash2 eos_helm, this is gamc
  abar = (EOS_ABAR-1)*vecLen   
  zbar = (EOS_ZBAR-1)*vecLen   
  entr = (EOS_ENTR-1)*vecLen

  !! For allocatable arrays, set them up now.
#ifndef FIXEDBLOCKSIZE
  call eos_vecAlloc(vecLen)
#endif

  if ( .not. eos_useStarkiller ) then

     do k = 1, vecLen

        tempRow(k)    = eosData(temp+k)
        denRow(k)     = eosData(dens+k)

        ! Note in Eos.F90, we assume the user knows what he's doing.  Eos_wrapped does not.

#ifdef FLASH_MULTISPECIES
        !Calculate the inverse in a way that allows for zero mass fractions
        call Multispecies_getSumInv(A, abarInv,massFrac((k-1)*NSPECIES+1:k*NSPECIES))
        abarRow(k) = 1.e0 / abarInv

        call Multispecies_getSumFrac(Z, zbarFrac, massFrac((k-1)*NSPECIES+1:k*NSPECIES))
        zbarRow(k) = abarRow(k) * zbarFrac
#else
        ! No multispecies defined, use default values (same as Gamma formulation)
        abarRow(k) = eos_singleSpeciesA
        zbarRow(k) = eos_singleSpeciesZ
#endif

     enddo

     eosData(abar+1:abar+vecLen) = abarRow(1:vecLen) 
     eosData(zbar+1:zbar+vecLen) = zbarRow(1:vecLen)

     !==============================================================================

     !      MODE_DENS_TEMP  temperature and density given

     !      Crank the EOS on the pipes filled above, then fill the FLASH arrays
     !      with the thermodynamic quantities returned by the EOS.

     if (mode==MODE_DENS_TEMP) then

        call eos_helm(1,vecLen)

        eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
        eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
        eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
        eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

        !==============================================================================
        !      MODE_DENS_EI  internal energy and density given

     else if (mode==MODE_DENS_EI) then

        ewantRow(1:vecLen)   = eosData(eint+1:eint+vecLen)   ! store desired internal energy for mode=2 case
        if (eos_forceConstantInput) then
           esaveRow = ewantRow
        end if
        ! Initialize the errors
        error(:) = 0.0e0

        ! Do the first eos call with all the zones in the pipe
        !  NOTE that eos_helm can ONLY operate in the equivalent of
        !  MODE_DENS_TEMP, as it returns pressure, energy and derivatives only
        !  So if you send in a crappy temperature here, you'll get a crappy starting
        !  position and the iteration won't converge.
        !  Initial temperature here is what is stored in the grid, even though we 
        !    SUSPECT this is not in equilibrium (or we wouldn't be calling Eos if it was fine)

        call eos_helm(vecBegin,vecEnd)
        !  Now eos_helm has returned ptotRow, etotRow, detRow, and gamcRow


        !  Create initial condition
        do k = vecBegin, vecEnd
           !  ewantRow is our desired EI input

           tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k))  & 
                &           / detRow(k)

           ! Don't allow the temperature to change by more than an order of magnitude 
           ! in a single iteration
           if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
                &           10.e0*tempRow(k)
           if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
                &           0.1e0*tempRow(k)

           ! Compute the error
           error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

           ! Store the new temperature
           tempRow(k) = tnew(k)

           ! Check if we are freezing, if so set the temperature to smallt, and adjust 
           ! the error so we don't wait for this one
           if (tempRow(k) .LT. eos_smallt) then
              tempRow(k) = eos_smallt
              error(k)    = 0.1*eos_tol
           endif

        enddo

        ! Loop over the zones individually now
        do k = vecBegin, vecEnd
           do i = 2, eos_maxNewton
              if (error(k)< eos_tol) goto 70

              call eos_helm(k,k)

              tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k))  & 
                   &              / detRow(k)

              ! Don't allow the temperature to change by more than an order of magnitude 
              ! in a single iteration
              if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
                   &              10.e0*tempRow(k)
              if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
                   &              0.1e0*tempRow(k)

              ! Compute the error
              error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

              ! Store the new temperature
              tempRow(k) = tnew(k)

              ! Check if we are freezing, if so set the temperature to eos_smallt, and adjust 
              ! the error so we don't wait for this one
              if (tempRow(k) .LT. eos_smallt) then
                 tempRow(k) = eos_smallt
                 error(k)    = .1*eos_tol
              endif

           end do  ! end of Newton iterations loop.  Failure drops below, success goes to 70

           ! Land here if too many iterations are needed -- failure

           print *, ' '
           print *, 'Newton-Raphson failed in subroutine Eos'
           print *, '(e and rho as input):'
           print *, ' '
           print *, 'too many iterations', eos_maxNewton
           print *, ' '
           print *, ' temp = ', tempRow(k)
           print *, ' dens = ', denRow(k)
           print *, ' pres = ', ptotRow(k)
           print *, ' eint = ', etotRow(k)

           if (present(diagFlag)) then
              diagFlag = k
              return
           else
              call Driver_abort('[Eos] Error: too many iterations in Helmholtz Eos.')
           endif


           ! Land here if the Newton iteration converged
           !  jumps out of the iterations, but then continues to the next vector location

70         continue           
        end do

        ! Crank through the entire eos one last time

        call eos_helm(vecBegin,vecEnd)

        ! Fill the FLASH arrays with the results.  

        !  In MODE_DENS_EI, we should be generating temperature and pressure (plus gamma and entropy)
        eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
        eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
        eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
        eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

        !  Update the energy to be the true energy, instead of the energy we were trying to meet
        !  ConstantInput LBR and KW believe this is WRONG -- the input arrays should not be changed
        if (eos_forceConstantInput)  then
           eosData(eint+1:eint+vecLen) = esaveRow(1:vecLen)
        else
           eosData(eint+1:eint+vecLen) = etotRow(1:vecLen)
        end if


        !==============================================================================

        !      MODE_DENS_PRES  pressure and density given

     else if (mode==MODE_DENS_PRES) then

        pwantRow(1:vecLen) = eosData(pres+1:pres+vecLen)   ! store desired pressure for mode=3 case
        if (eos_forceConstantInput) then
           psaveRow = pwantRow
        end if
        ! Initialize the errors
        error(:) = 0.0e0

        ! Do the first eos call with all the zones in the pipe
        call eos_helm(vecBegin,vecEnd)

        do k = vecBegin, vecEnd

           tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k))  & 
                &           / dptRow(k)

           ! Don't allow the temperature to change by more than an order of magnitude 
           ! in a single iteration
           if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
                &           10.e0*tempRow(k)
           if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
                &           0.1e0*tempRow(k)

           ! Compute the error
           error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

           ! Store the new temperature
           tempRow(k) = tnew(k)

           ! Check if we are freezing, if so set the temperature to smallt, and adjust 
           ! the error so we don't wait for this one
           if (tempRow(k) .LT. eos_smallt) then
              tempRow(k) = eos_smallt
              error(k)    = 0.1*eos_tol
           endif

        enddo

        ! Loop over the zones individually now
        do k = vecBegin, vecEnd

           do i = 2, eos_maxNewton

              if (error(k) .LT. eos_tol) goto 170

              ! do eos only over this single item
              call eos_helm(k,k)

              tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k))  & 
                   &              / dptRow(k)

              ! Don't allow the temperature to change by more than an order of magnitude 
              ! in a single iteration
              if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
                   &              10.e0*tempRow(k)
              if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
                   &              0.1e0*tempRow(k)

              ! Compute the error
              error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

              ! Store the new temperature
              tempRow(k) = tnew(k)

              ! Check if we are freezing, if so set the temperature to eos_smallt, and adjust 
              ! the error so we don't wait for this one
              if (tempRow(k) .LT. eos_smallt) then
                 tempRow(k) = eos_smallt
                 error(k)    = .1*eos_tol
              endif

           enddo

           ! Land here if too many iterations are needed

           print *, ' '
           print *, 'Newton-Raphson failed in Helmholtz Eos:'
           print *, '(p and rho as input)'
           print *, ' '
           print *, 'too many iterations'
           print *, ' '
           print *, ' k    = ', k,vecBegin,vecEnd
           print *, ' temp = ', tempRow(k)
           print *, ' dens = ', denRow(k)
           print *, ' pres = ', ptotRow(k)
           print *, ' eint = ', etotRow(k)
           print *, ' abar = ', abarRow(k)
           print *, ' zbar = ', zbarRow(k)

           if (present(diagFlag)) then
              diagFlag = k
              return
           else
              call Driver_abort('[Eos] Error: too many iterations in Helmholtz Eos.')
           endif


           ! Land here if the Newton iteration converged

170        continue

        enddo

        ! Crank through the entire eos one last time

        call eos_helm(vecBegin,vecEnd)

        ! Fill the FLASH arrays with the results.  
        eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
        eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
        eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
        eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

        ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
        !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
        if (eos_forceConstantInput) then
           eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
        else
           eosData(pres+1:pres+vecLen) = ptotRow(1:vecLen)
        end if


        !==============================================================================

        ! Unknown EOS mode selected

     else if (mode .NE. MODE_EOS_NOP) then
        if (eos_meshMe .EQ. MASTER_PE) print*, '[Eos] Error: unknown input mode', mode
        call Driver_abort('[Eos] Error: unknown input mode in subroutine Eos')
     end if

     ! Get the optional values
     if(present(mask)) then
        ! Entropy derivatives
        if(mask(EOS_DST)) then
           dst = (EOS_DST-1)*vecLen
           eosData(dst+1:dst+vecLen) = dstRow(1:vecLen)
        end if
        if(mask(EOS_DSD)) then
           dsd = (EOS_DSD-1)*vecLen
           eosData(dsd+1:dsd+vecLen) = dsdRow(1:vecLen)
        end if
        if(mask(EOS_DPT)) then
           dpt = (EOS_DPT-1)*vecLen
           eosData(dpt+1:dpt+vecLen) = dptRow(1:vecLen)
        end if
        if(mask(EOS_DPD)) then
           dpd = (EOS_DPD-1)*vecLen
           eosData(dpd+1:dpd+vecLen) = dpdRow(1:vecLen)
        end if
        if(mask(EOS_DET))then
           det = (EOS_DET-1)*vecLen
           eosData(det+1:det+vecLen) = detRow(1:vecLen)
        end if
        if(mask(EOS_DED))then 
           ded = (EOS_DED-1)*vecLen
           eosData(ded+1:ded+vecLen) = dedRow(1:vecLen)
        end if
        if(mask(EOS_DEA))then 
           dea = (EOS_DEA-1)*vecLen
           eosData(dea+1:dea+vecLen) = deaRow(1:vecLen)
        end if
        if(mask(EOS_DEZ))then 
           dez = (EOS_DEZ-1)*vecLen
           eosData(dez+1:dez+vecLen) = dezRow(1:vecLen)
        end if
        if(mask(EOS_PEL))then 
           pel = (EOS_PEL-1)*vecLen
           eosData(pel+1:pel+vecLen) = pelRow(1:vecLen)
        end if
        if(mask(EOS_NE))then 
           ne = (EOS_NE-1)*vecLen
           eosData(ne+1:ne+vecLen) = neRow(1:vecLen)
        end if
        if(mask(EOS_ETA))then 
           eta = (EOS_ETA-1)*vecLen
           eosData(eta+1:eta+vecLen) = etaRow(1:vecLen)
        end if
        if(mask(EOS_DETAT))then
           detat = (EOS_DETAT-1)*vecLen
           eosData(detat+1:detat+vecLen) = detatRow(1:vecLen)
        end if

        if(mask(EOS_CV))then
           if(mask(EOS_DET)) then
              c_v = (EOS_CV-1)*vecLen
              eosData(c_v+1:c_v+vecLen) = cvRow(1:vecLen)
           else
              call Driver_abort("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
           end if
        end if

        if(mask(EOS_CP))then
           if(mask(EOS_CV).and.mask(EOS_DET)) then
              c_p = (EOS_CP-1)*vecLen
              eosData(c_p+1:c_p+vecLen) = cpRow(1:vecLen)
           else
              call Driver_abort("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
           end if
        end if
     end if

  else ! if eos_useStarkiller

     select case (mode)
     case (MODE_DENS_TEMP)
        eos_input = eos_input_rt
     case (MODE_DENS_EI)
        eos_input = eos_input_re
     case (MODE_DENS_PRES)
        eos_input = eos_input_rp
     case default
        eos_input = -1
     end select

     if ( eos_input > 0 ) then

        ! Use OpenACC if vecLen is large enough and not inside OpenMP parallel region
        ! Use OpenMP if not using OpenACC and not inside OpenMP parallel region
#ifdef _OPENMP
        in_par = omp_in_parallel()
#else
        in_par = .false.
#endif
#ifdef _OPENACC
        do_acc = ( vecLen >= eos_vecLenACC .and. .not. in_par )
#else
        do_acc = .false.
#endif

        local_eos_state(1:vecLen) % rho = eosData(dens+1:dens+vecLen)
        local_eos_state(1:vecLen) % T   = eosData(temp+1:temp+vecLen)
        local_eos_state(1:vecLen) % p   = eosData(pres+1:pres+vecLen)
        local_eos_state(1:vecLen) % e   = eosData(eint+1:eint+vecLen)
        do k = 1, vecLen
           comp_state % xn(1:NSPECIES) = massFrac((k-1)*NSPECIES+1:k*NSPECIES)
           call composition(comp_state)
           local_eos_state(k) % mu_e = comp_state % mu_e
           local_eos_state(k) % y_e  = comp_state % y_e
           local_eos_state(k) % abar = comp_state % abar
           local_eos_state(k) % zbar = comp_state % zbar
           call clean_state(local_eos_state(k))
        end do
        eosData(abar+1:abar+vecLen) = local_eos_state(1:vecLen) % abar
        eosData(zbar+1:zbar+vecLen) = local_eos_state(1:vecLen) % zbar
        if (eos_forceConstantInput) then
           if (mode==MODE_DENS_EI) then
              esaveRow(1:vecLen) = local_eos_state(1:vecLen) % e
           else if (mode==MODE_DENS_PRES) then
              psaveRow(1:vecLen) = local_eos_state(1:vecLen) % p
           end if
        end if

        if ( do_acc ) then

           reduced_eos_state(1:vecLen) % rho  = local_eos_state(1:vecLen) % rho
           reduced_eos_state(1:vecLen) % T    = local_eos_state(1:vecLen) % T
           reduced_eos_state(1:vecLen) % p    = local_eos_state(1:vecLen) % p
           reduced_eos_state(1:vecLen) % e    = local_eos_state(1:vecLen) % e
           reduced_eos_state(1:vecLen) % y_e  = local_eos_state(1:vecLen) % y_e
           reduced_eos_state(1:vecLen) % abar = local_eos_state(1:vecLen) % abar
           reduced_eos_state(1:vecLen) % zbar = local_eos_state(1:vecLen) % zbar

           !$acc update device(eos_input,reduced_eos_state(1:vecLen))
           !$acc kernels
           eos_state(1:vecLen) % rho  = reduced_eos_state(1:vecLen) % rho
           eos_state(1:vecLen) % T    = reduced_eos_state(1:vecLen) % T
           eos_state(1:vecLen) % p    = reduced_eos_state(1:vecLen) % p
           eos_state(1:vecLen) % e    = reduced_eos_state(1:vecLen) % e
           eos_state(1:vecLen) % y_e  = reduced_eos_state(1:vecLen) % y_e
           eos_state(1:vecLen) % abar = reduced_eos_state(1:vecLen) % abar
           eos_state(1:vecLen) % zbar = reduced_eos_state(1:vecLen) % zbar
           eos_state(1:vecLen) % nr   = 0
           !$acc loop
           do k = 1, vecLen
              call actual_eos(eos_input,eos_state(k))
              if ( eos_input /= eos_input_rt ) then
                 ! add an extra call to match the old starkiller == 0 behavior
                 call actual_eos(eos_input_rt,eos_state(k))
               end if
           end do
           !$acc end kernels
           !$acc update host(eos_state(1:vecLen))
           local_eos_state(1:vecLen) = eos_state(1:vecLen)

        else

           !$omp parallel if (.not. in_par) &
           !$omp default(none) &
           !$omp shared(eos_input,local_eos_state,vecLen) &
           !$omp private(k)

           !$omp do schedule(static)
           do k = 1, vecLen
              call actual_eos(eos_input,local_eos_state(k))
              if ( eos_input /= eos_input_rt ) then
                 ! add an extra call to match the old starkiller == 0 behavior
                 call actual_eos(eos_input_rt,local_eos_state(k))
              end if
           end do
           !$omp end do
           !$omp end parallel

        end if
     end if

     if ( any( local_eos_state(1:vecLen) % nr > max_newton ) ) then
        do k = 1, vecLen
           if ( local_eos_state(k) % nr > max_newton ) then
              print *, ' '
              print *, 'Newton-Raphson failed in Helmholtz Eos:'
              if (mode==MODE_DENS_PRES) then
                 print *, '(p and rho as input)'
              else if (mode==MODE_DENS_EI) then
                 print *, '(e and rho as input)'
              end if
              print *, ' '
              print *, 'too many iterations'
              print *, ' '
              print *, ' k    = ', k
              print *, ' temp = ', local_eos_state(k) % T
              print *, ' dens = ', local_eos_state(k) % rho
              print *, ' pres = ', local_eos_state(k) % p
              print *, ' eint = ', local_eos_state(k) % e
              print *, ' abar = ', local_eos_state(k) % abar
              print *, ' zbar = ', local_eos_state(k) % zbar
              if (present(diagFlag)) then
                 diagFlag = k
              end if
           end if
        end do
        if (present(diagFlag)) then
           return
        else
           call Driver_abort('[Eos] Error: too many iterations in Helmholtz Eos.')
        endif
     end if

     if (mode==MODE_DENS_TEMP) then
        eosData(pres+1:pres+vecLen) = local_eos_state(1:vecLen) % p
        eosData(eint+1:eint+vecLen) = local_eos_state(1:vecLen) % e
        eosData(gamc+1:gamc+vecLen) = local_eos_state(1:vecLen) % gam1
        eosData(entr+1:entr+vecLen) = local_eos_state(1:vecLen) % s
     else if (mode==MODE_DENS_EI) then
        eosData(temp+1:temp+vecLen) = local_eos_state(1:vecLen) % T
        eosData(pres+1:pres+vecLen) = local_eos_state(1:vecLen) % p
        eosData(gamc+1:gamc+vecLen) = local_eos_state(1:vecLen) % gam1
        eosData(entr+1:entr+vecLen) = local_eos_state(1:vecLen) % s
        if (eos_forceConstantInput) then
           eosData(eint+1:eint+vecLen) = esaveRow(1:vecLen)
        else
           eosData(eint+1:eint+vecLen) = local_eos_state(1:vecLen) % e
        end if
     else if (mode==MODE_DENS_PRES) then
        eosData(temp+1:temp+vecLen) = local_eos_state(1:vecLen) % T
        eosData(eint+1:eint+vecLen) = local_eos_state(1:vecLen) % e
        eosData(gamc+1:gamc+vecLen) = local_eos_state(1:vecLen) % gam1
        eosData(entr+1:entr+vecLen) = local_eos_state(1:vecLen) % s
        if (eos_forceConstantInput) then
           eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
        else
           eosData(pres+1:pres+vecLen) = local_eos_state(1:vecLen) % p
        end if
     else if (mode .NE. MODE_EOS_NOP) then
        if (eos_meshMe .EQ. MASTER_PE) print*, '[Eos] Error: unknown input mode', mode
        call Driver_abort('[Eos] Error: unknown input mode in subroutine Eos')
     end if

     ! Get the optional values
     if(present(mask)) then

        ! Entropy derivatives
        if(mask(EOS_DST)) then
           dst = (EOS_DST-1)*vecLen
           eosData(dst+1:dst+vecLen) = local_eos_state(1:vecLen) % dsdT
        end if
        if(mask(EOS_DSD)) then
           dsd = (EOS_DSD-1)*vecLen
           eosData(dsd+1:dsd+vecLen) = local_eos_state(1:vecLen) % dsdr
        end if
        if(mask(EOS_DPT)) then
           dpt = (EOS_DPT-1)*vecLen
           eosData(dpt+1:dpt+vecLen) = local_eos_state(1:vecLen) % dpdT
        end if
        if(mask(EOS_DPD)) then
           dpd = (EOS_DPD-1)*vecLen
           eosData(dpd+1:dpd+vecLen) = local_eos_state(1:vecLen) % dpdr
        end if
        if(mask(EOS_DET))then
           det = (EOS_DET-1)*vecLen
           eosData(det+1:det+vecLen) = local_eos_state(1:vecLen) % dedT
        end if
        if(mask(EOS_DED))then 
           ded = (EOS_DED-1)*vecLen
           eosData(ded+1:ded+vecLen) = local_eos_state(1:vecLen) % dedr
        end if
        if(mask(EOS_DEA))then 
           dea = (EOS_DEA-1)*vecLen
           eosData(dea+1:dea+vecLen) = local_eos_state(1:vecLen) % dedA
        end if
        if(mask(EOS_DEZ))then 
           dez = (EOS_DEZ-1)*vecLen
           eosData(dez+1:dez+vecLen) = local_eos_state(1:vecLen) % dedZ
        end if
        if(mask(EOS_PEL))then 
           pel = (EOS_PEL-1)*vecLen
           eosData(pel+1:pel+vecLen) = local_eos_state(1:vecLen) % pele
        end if
        if(mask(EOS_NE))then 
           ne = (EOS_NE-1)*vecLen
           eosData(ne+1:ne+vecLen) = local_eos_state(1:vecLen) % xne
        end if
        if(mask(EOS_ETA))then 
           eta = (EOS_ETA-1)*vecLen
           eosData(eta+1:eta+vecLen) = local_eos_state(1:vecLen) % eta
        end if
        if(mask(EOS_DETAT))then
           detat = (EOS_DETAT-1)*vecLen
           eosData(detat+1:detat+vecLen) = local_eos_state(1:vecLen) % detadt
        end if

        if(mask(EOS_CV))then
           if(mask(EOS_DET)) then
              c_v = (EOS_CV-1)*vecLen
              eosData(c_v+1:c_v+vecLen) = local_eos_state(1:vecLen) % cv
           else
              call Driver_abort("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
           end if
        end if

        if(mask(EOS_CP))then
           if(mask(EOS_CV).and.mask(EOS_DET)) then
              c_p = (EOS_CP-1)*vecLen
              eosData(c_p+1:c_p+vecLen) = local_eos_state(1:vecLen) % cp
           else
              call Driver_abort("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
           end if
        end if
     end if
  end if

  !! Close up arrays if previously allocated
#ifndef FIXEDBLOCKSIZE  
  call eos_vecDealloc()
#endif

  return

end subroutine Eos


