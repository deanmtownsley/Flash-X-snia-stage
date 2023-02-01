!!****h* source/Simulation/SimulationMain/unitTest/MoL/sim_interface
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
!! This is the header file for the MoL unit test module
!! that defines its public interfaces.
!!***
module sim_interface

   implicit none

   interface
      subroutine sim_verifySolution(t, valid, maxError)
         implicit none
         real, intent(in) :: t
         logical, intent(out) :: valid
         real, intent(out) :: maxError
      end subroutine sim_verifySolution
   end interface

end module sim_interface
