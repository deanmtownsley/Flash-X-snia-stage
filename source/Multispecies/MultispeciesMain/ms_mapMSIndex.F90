!!****if* source/Multispecies/MultispeciesMain/ms_mapMSIndex
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
!!  ms_mapMSIndex
!!
!! SYNOPSIS
!!
!!  ms_mapMSIndex(integer, intent(in)  :: name,
!!                integer, intent(out)  :: msindex)
!!
!! DESCRIPTION
!!
!!  Maps the integer value of the name of the species to the correct
!!  index in the multispecies array.
!!  
!!  This is done under the hood as so not to make the user learn 2
!!  different naming schemes
!!
!! ARGUMENTS
!!
!!    name - name of species define in Simulation.h ie NI56_SPEC
!!    msindex - the corresponding index in the multispecies array
!!
!!***



subroutine ms_mapMSIndex(name, msindex)

  implicit none

#include "Simulation.h"

  integer, intent(in)   :: name
  integer, intent(out)  :: msindex

  msindex = name - NPROP_VARS

end subroutine ms_mapMSIndex
