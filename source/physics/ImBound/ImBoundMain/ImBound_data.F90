!!****if* source/physics/ImBound/ImBoundMain/ImBound_data
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
!!  ImBound_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE ImBound_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the ImBound module.
!!
!!***

module ImBound_data

#include "Simulation.h"
#include "constants.h"

   use ImBound_type, ONLY: ImBound_type_t

   implicit none

   logical, save :: ib_useImBound
   logical, save :: ib_enableSelectiveMapping
   logical, save :: ib_bruteForceMapping

   integer, save :: ib_lsIt
   integer, save :: ib_meshMe
   integer, save :: ib_meshNumProcs
   integer, save :: ib_meshComm
   integer, save :: ib_nstep
   integer, save :: ib_annQueries
   integer, save, dimension(:), allocatable :: ib_annIdx

   real, dimension(2), save :: ib_alfa
   real, save :: ib_rhoa
   real, save :: ib_gama

   type(ImBound_type_t), save, dimension(:), allocatable, target :: ib_bodyInfo

   character(len=20), save :: ib_bodyName
   integer, save :: ib_numBodies

   real, save :: ib_invReynolds
   logical, save :: ib_withIncompNS
   integer, save :: ib_iVelFVar, ib_iPGradVar, ib_iVFrcVar
   integer, save :: ib_iVelXVar, ib_iVelYVar, ib_iVelZVar, ib_iPresVar

   real, save :: ib_forceBuffer(MDIM)

end module ImBound_data
