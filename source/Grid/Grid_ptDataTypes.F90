module Grid_ptDataTypes
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

  ! Grid_ptContainer_t is currently not used anywhere.
  type Grid_ptContainer_t
  end type Grid_ptContainer_t

  ! In a non-stub implementation of Grid_ptDataTypes, Grid_particle_t
  ! should basically become an alias for the type of mesh-owned particles.
  type Grid_particle_t
  end type Grid_particle_t

end module Grid_ptDataTypes
