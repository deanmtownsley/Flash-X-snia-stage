!!****f* source/Simulation/SimulationMain/CCSN_Chimera/Simulation_finalize
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This function cleans up the Simulation unit, deallocates memory, etc.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Simulation_finalize()

  use chimera_model_module

  implicit none

  call close_chimera_file

end subroutine Simulation_finalize
