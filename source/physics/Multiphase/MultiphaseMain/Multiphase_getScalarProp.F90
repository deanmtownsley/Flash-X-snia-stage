subroutine Multiphase_getScalarProp(name, value)
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

   use Multiphase_data
   use Driver_interface, ONLY: Driver_abort

   implicit none
   character(len=*), intent(in)  :: name
   real, intent(out)             :: value

   character(len=100)            :: errorMessage

   select case (name)
   case ("Gas_Density", "gas_density", "GAS_DENSITY")
      value = mph_rhoGas
   case ("Gas_Viscosity", "gas_viscosity", "GAS_VISCOSITY")
      value = mph_muGas
   case ("Gas_Conductivity", "gas_conductivity", "GAS_CONDUCTIVITY")
      value = mph_thcoGas
   case ("Gas_Specific_Heat", "gas_specific_heat", "GAS_SPECIFIC_HEAT")
      value = mph_CpGas
   case ("Saturation_Temperature", "saturation_temperature", "SATURATION_TEMPERATURE")
      value = mph_Tsat
   case default
      value = 0.
      write (errorMessage, *) '[Multiphase_getScalarProp] Unknown scalar: ', name
      call Driver_abort(errorMessage)
   end select

   return

end subroutine Multiphase_getScalarProp
