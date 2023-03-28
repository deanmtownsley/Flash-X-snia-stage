!!****f* source/Simulation/Simulation_molPostRegrid
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
!!  NAME
!!
!!      Simulation_molPostRegrid
!!
!!  SYNOPSIS
!!
!!      call Simulation_molPostRegrid(real, intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Perform any post-regrid work.  If a physics unit stored primitive
!!      variables in UNK in additiona to the evolved conserved variables
!!      used by MoL, a con2prim should be called here to maintain consistency.
!!      Refined primitive variables will not necessarily be consistent with
!!      the results of a con2prim call on the refined conserved variables.
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!
!!***
subroutine Simulation_molPostRegrid(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Simulation_molPostRegrid
