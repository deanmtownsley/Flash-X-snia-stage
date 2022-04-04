!!***if* source/physics/HeatAD/HeatADMain/varDiffusion/HeatAD_diffusion
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

subroutine HeatAD_getGridVar(name, value)

#include"Simulation.h"

  use Driver_interface, only: Driver_abort

  implicit none
  character(len=*), intent(in)  :: name
  integer, intent(out)          :: value

  select case(name)
     case("Center_Temperature","center_temperature","CENTER_TEMPERATURE")
       value = TEMP_VAR
     case("Center_Thermal_Diffusivity","center_thermal_diffusivity","CENTER_THERMAL_DIFFUSIVITY")
       value = ALPH_VAR
     case("Center_Thermal_Forcing","center_thermal_forcing","CENTER_THERMAL_FORCING")
       value = TFRC_VAR
     case default
       value = -1
       print *,"Error in setting grid var: ",name 
       call Driver_abort("HeatAD_getGridVar: Unknown HeatAD Grid Variable")
  end select  

  return

end subroutine HeatAD_getGridVar
