!****h* source/physics/SolidMechanics/Imbound_interface
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
!!  Imbound_interface
!!
!! SYNOPSIS
!!
!!  use Imbound_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Immersed boundary (IB)
!! unit that defines its public interfaces.
!!
!!***

Module SolidMechanics_interface

   implicit none

#include "Simulation.h"

   interface !SolidMechanics_init
      subroutine SolidMechanics_init(restart)
         implicit none
         logical, INTENT(IN) :: restart
      end subroutine SolidMechanics_init
   end interface

   interface  !SolidMechanics_finalize
      subroutine SolidMechanics_finalize()
         implicit none
      end subroutine SolidMechanics_finalize
   end interface

   interface
      subroutine SolidMechanics_updateBodyForce(body, time, dt)
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         type(ImBound_type_t), intent(inout) :: body
         real, intent(in) :: time, dt
      end subroutine SolidMechanics_updateBodyForce
   end interface

end Module SolidMechanics_interface
