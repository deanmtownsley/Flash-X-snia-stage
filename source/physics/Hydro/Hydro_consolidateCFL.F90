!!****f* source/physics/Hydro/Hydro_consolidateCFL
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
!!
!!  Hydro_consolidateCFL
!!
!! SYNOPSIS
!!
!!  call Hydro_consolidateCFL()
!!
!! DESCRIPTION
!!
!!  Consolidate the CFL factors of different MPI tasks by a global
!!  reduction operation.
!!
!!  Different MPI tasks may have computed different values for an updated Hydro
!!  CFL factor. This routine should be called before Hydro_computeDt is called
!!  to compute an updated Hydro time step.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! SIDE EFFECTS
!!
!!  The CFL factor used by the Hydro unit is updated.
!!
!! NOTES
!!
!!  This routine must be called collectively by all MPI tasks in the global
!!  Hydro communicator.
!!
!!***

subroutine Hydro_consolidateCFL()

  implicit none

end subroutine Hydro_consolidateCFL
