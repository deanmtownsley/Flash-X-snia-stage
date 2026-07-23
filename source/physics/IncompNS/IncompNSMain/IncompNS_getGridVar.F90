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

#include"Simulation.h"

subroutine IncompNS_getGridVar(name, value)

   use Driver_interface, only: Driver_abort

   implicit none

   character(len=*), intent(in)  :: name
   integer, intent(out)          :: value

   select case (name)
   case ("Face_Velocity", "face_velocity", "FACE_VELOCITY")
      value = VELC_FACE_VAR
   case ("Center_Pressure", "center_pressure", "CENTER_PRESSURE")
      value = PRES_VAR
   case ("Center_Divergence", "center_divergence", "CENTER_DIVERGENCE")
      value = DUST_VAR
   case ("Face_Pressure_Grad", "face_pressure_grad", "FACE_PRESSURE_GRAD")
      value = PGN0_FACE_VAR
   case ("Face_Vel_Forcing", "face_vel_forcing", "FACE_VEL_FORCING")
      value = VFRC_FACE_VAR
#ifdef INCOMPNS_VARDENS
   case ("Face_Pressure_Jump", "face_pressure_jump", "FACE_PRESSURE_JUMP")
      value = SIGM_FACE_VAR
   case ("Face_Density", "face_density", "FACE_DENSITY")
      value = RHOF_FACE_VAR
   case ("Center_Density", "center_density", "CENTER_DENSITY")
      value = RHOC_VAR
   case ("Center_Viscosity", "center_viscosity", "CENTER_VISCOSITY")
      value = VISC_VAR
#endif
#ifdef INCOMPNS_EXTRAS
   case ("Center_Velx", "center_velx", "CENTER_VELX")
      value = VELX_VAR
   case ("Center_Vely", "center_vely", "CENTER_VELY")
      value = VELY_VAR
   case ("Center_Velz", "center_velz", "CENTER_VELZ")
      value = VELZ_VAR
#endif
   case default
      value = -1
      print *, "Error in setting grid var: ", name
      call Driver_abort("IncompNS_getGridVar: Unknown IncompNS Grid Variable")
   end select

   return

end subroutine IncompNS_getGridVar
