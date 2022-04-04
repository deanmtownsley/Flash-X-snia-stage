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
subroutine IncompNS_getGridVar(name, value)

#include"Simulation.h"

  use Driver_interface, only: Driver_abort

  implicit none
  character(len=*), intent(in)  :: name
  integer, intent(out)          :: value

  select case(name)
     case("Face_Velocity","face_velocity","FACE_VELOCITY")
       value = VELC_FACE_VAR
     case("Center_Pressure","center_pressure","CENTER_PRESSURE")
       value = PRES_VAR
     case("Center_Divergence","center_divergence","CENTER_DIVERGENCE")
       value = DUST_VAR
     case("Face_Density","face_density","FACE_DENSITY")
       value = RHOF_FACE_VAR
     case("Center_Viscosity","center_viscosity","CENTER_VISCOSITY")
       value = VISC_VAR
     case("Face_Pressure_Jump","face_pressure_jump","FACE_PRESSURE_JUMP")
       value = SIGM_FACE_VAR
     case default
       value = -1
       print *,"Error in setting grid var: ",name 
       call Driver_abort("IncompNS_getGridVar: Unknown IncompNS Grid Variable")
  end select  

  return

end subroutine IncompNS_getGridVar
