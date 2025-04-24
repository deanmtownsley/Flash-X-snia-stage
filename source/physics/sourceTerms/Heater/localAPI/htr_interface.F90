!!****if* source/physics/sourceTerms/Heater/localAPI/htr_interface
!!
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
!!
!! SYNOPSIS
!!  htr_interface()
!!
!! DESCRIPTION
!!  This is an interface specific for heater geometry and specifications
!!
!!***

#include "constants.h"
#include "Simulation.h"

Module htr_interface

   implicit none

   interface
      subroutine htr_read(heaterID, heaterFile)
         integer, intent(in)          :: heaterID
         character(len=*), intent(in) :: heaterFile
      end subroutine htr_read
   end interface

End module htr_interface
