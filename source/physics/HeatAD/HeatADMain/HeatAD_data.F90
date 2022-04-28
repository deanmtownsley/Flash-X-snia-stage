module HeatAD_data
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
     
     implicit none
     real, save :: ht_Prandtl ! Prandtl Number
     real, save :: ht_Nusselt ! Nusselt Number
     real, save :: ht_Biot    ! Biot Number
     real, save :: ht_Twall_high, ht_Twall_low ! Non Dimensional wall Temp
     real, save :: ht_Tbulk, ht_Tsat ! Bulk and saturation temperature
     real, save :: ht_invReynolds ! Inverse Re 

     real, save :: ht_thcoGas ! Gas thermal conductivity
     real, save :: ht_CpGas   ! Gas specific heat

     integer, save :: ht_meshMe
     integer, save :: ht_meshNumProcs
     integer, save :: ht_meshComm     

     real, save :: ht_alfa
     real, save :: ht_rhoa, ht_gama

     logical, save :: ht_useHeatAD

     integer, save :: ht_iVelFVar

end module HeatAD_data
