!****h* source/physics/ImBound/Imbound_interface
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

Module ImBound_interface

  implicit none

#include "Simulation.h"


  interface !ImBound_init
     subroutine ImBound_init(restart)
     implicit none
     logical, INTENT(IN) :: restart
     end subroutine ImBound_init
  end interface

  interface  !ImBound_finalize
     subroutine ImBound_finalize()
     implicit none
     end subroutine ImBound_finalize
  end interface 

  interface
     subroutine ImBound_getScalarProp(name, value)
     implicit none
     character(len=*), intent(in)  :: name
     real, intent(out)             :: value
     end subroutine ImBound_getScalarProp
  end interface

end Module ImBound_interface
