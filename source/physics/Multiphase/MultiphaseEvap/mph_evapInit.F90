!!****if* source/physics/Multiphase/MultiphaseEvap/mph_evapInit
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
!!
!!
!!******
subroutine mph_evapInit()

   use Multiphase_data
   use HeatAD_interface, ONLY: HeatAD_getGridVar

   implicit none
   call HeatAD_getGridVar('Center_Thermal_Diffusivity', mph_iAlphaCVar)
   call HeatAD_getGridVar('Center_Temperature', mph_iTempVar)
   call HeatAD_getGridVar('Center_Thermal_Forcing', mph_iGfmVar)

end subroutine mph_evapInit
