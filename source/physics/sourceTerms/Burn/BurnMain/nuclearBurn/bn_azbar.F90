!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_azbar
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
!! bn_azbar
!!
!! SYNOPSIS
!!
!! bn_azbar()
!!
!! DESCRIPTION
!!
!!  routine azbar computes composition variables from the mass fractions 
!!  
!!  given the mass fractions in xmass(i), return the molar abundances ymass(i), 
!!  total number of moles per gram ytot1, the mean number of nucleons abar, 
!!  mean nucleon charge zbar, mean nucleon charge squared z2bar, and the
!!  electron mole number bye.
!!
!! NOTES
!!   the output variables are stored in data structure Burn_dataEOS
!!  
!!***

subroutine bn_azbar()

  use Burn_dataEOS, ONLY:  abar,zbar,z2bar,ytot1,bye
  use Burn_data, ONLY: xmass,ymass,aion,zion,bion,aioninv,zionsq

  implicit none

#include "constants.h"
#include "Simulation.h"

  !!  local declarations
  integer          i
  real             zbarxx,z2barxx

  zbarxx  = 0.0e0
  z2barxx = 0.0e0
  ytot1   = 0.0e0

  do i=1,NSPECIES
     ymass(i) = xmass(i) * aioninv(i)
     zbarxx   = zbarxx + zion(i) * ymass(i)
     z2barxx  = z2barxx + zionsq(i) * ymass(i)
     ytot1    = ytot1 + ymass(i)
  enddo

  abar   = 1.0e0/ytot1
  zbar   = zbarxx * abar
  z2bar  = z2barxx * abar
  bye    = zbar * ytot1

  return
end subroutine bn_azbar
