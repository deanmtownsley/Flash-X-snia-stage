!!****f* source/Simulation/Simulation_freeUserArrays
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
!!  Simulation_freeUserArrays
!!
!! SYNOPSIS
!!  Simulation_freeUserArrays()
!!
!! DESCRIPTION
!!  This is where the user should place code for a setup that needs to
!!  free memory created in Simulation_init that may have been used in
!!  other aspects of initialization but is no longer needed once
!!  the initialization stage is finished. 
!!
!! ARGUMENTS
!!
!!***
subroutine Simulation_freeUserArrays()

  implicit none

  return

end subroutine Simulation_freeUserArrays
