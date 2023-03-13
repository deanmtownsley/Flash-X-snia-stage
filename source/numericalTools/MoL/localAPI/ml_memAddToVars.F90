!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief ml_memAddToVars stub

!> @ingroup MoLPrivate
!!
!! @brief Perform a linear combination of evolved variable
!!
!! @details
!! @anchor ml_memAddToVars_stub
!!
!! Perform a linear combination of source terms into the specified destination
!! for all variables evolved by MoL.  The destination can either be the evolved
!! variables in UNK or their corresponding locations in one of the MoL-specific
!! scratch memory locations.  The source terms will only-ever be taken from
!! MoL-specific scratch memory locations.  The linear combinations will take one
!! of the following form:
!!
!!    `dst = dstFac*dst + fac(1)*src(1) + ... + fac(nsrcs)*src(nsrcs)`
!!
!! Valid locations include (defined in \ref MoL.h):
!! - `MOL_EVOLVED` : Evolved variables in UNK
!! - `MOL_INITIAL` : Copy of the evolved variables at the start of a timestep
!! - `MOL_RHS`     : The currently-being-calculated RHS terms
!! - other         : Each integrator may specify some additional number of
!!                   of scratch-memory for intermediate stages/RHS terms
!!
!! @param dst     Index of the destination location to store the linear combination
!! @param dstFac  Scaling factor for the destination - set this to zero to overwrite
!!                the existing value
!! @param nsrcs   Number of source terms for the general N-src implementation
!! @param srcs    Array of source index locations in MoL scratch memory
!! @param facs    Array of scaling factors for each source term
subroutine ml_memAddToVars(dst, dstFac, nsrcs, srcs, facs)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: nsrcs
   integer, intent(in) :: srcs(nsrcs)
   real, intent(in) :: facs(nsrcs)

   return
end subroutine ml_memAddToVars
