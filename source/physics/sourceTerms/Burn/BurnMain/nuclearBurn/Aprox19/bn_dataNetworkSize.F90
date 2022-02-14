!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_dataNetworkSize
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
!!  bn_dataNetworkSize
!!
!!
!! SYNOPSIS
!! 
!!  use bn_dataNetworkSize
!!
!! DESCRIPTION
!!
!!  Contains variables indicating the nuclear network size.
!!
!! NOTES
!!
!!   In Flash2, this routine was called network_size.fh
!!
!!***

Module bn_dataNetworkSize
  implicit none
  integer, parameter :: nrat = 86
  integer, parameter :: nratp1 = nrat+1

!! For communication between bn_networkSparseJakob and bn_networkSparsePonters
!!  were in the common /elccm1/  eloc,nterms
  integer, parameter ::     neloc=112
  integer          eloc(neloc),nterms



end Module bn_dataNetworkSize
