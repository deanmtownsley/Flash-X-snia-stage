subroutine Multiphase_getGridVar(name, value)
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

   use Driver_interface, only: Driver_abort

   implicit none
   character(len=*), intent(in)  :: name
   integer, intent(out)          :: value

   select case (name)
   case ("Center_Levelset", "center_levelset", "CENTER_LEVELSET")
      value = DFUN_VAR
   case ("Center_Levelset_Forcing", "center_Levelset_forcing", "CENTER_LEVELSET_FORCING")
      value = DFRC_VAR
#ifdef MULTIPHASE_EVAPORATION
   case ("Center_Massflux", "center_massflux", "CENTER_MASSFLUX")
      value = MFLX_VAR
   case ("Center_Hflux_Liquid", "center_hflux_liquid", "CENTER_HFLUX_LIQUID")
      value = HFLQ_VAR
   case ("Center_Hflux_Gas", "center_hflux_gas", "CENTER_HFLUX_GAS")
      value = HFGS_VAR
#endif
   case default
      value = -1
      print *, "Error in setting grid var: ", name
      call Driver_abort("Multiphase_getGridVar: Unknown Multiphase Grid Variable")
   end select

   return

end subroutine Multiphase_getGridVar
