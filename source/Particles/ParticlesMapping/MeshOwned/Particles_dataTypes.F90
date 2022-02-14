module Particles_dataTypes
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
  ! Particles_meshParticle_t is an alias for the type of mesh-owned particles.
  ! It is assumed here that we are configured with AMReX.

  use pt_amrexParticleMod, ONLY : Particles_meshParticle_t => pt_amrexParticle_t

end module Particles_dataTypes
