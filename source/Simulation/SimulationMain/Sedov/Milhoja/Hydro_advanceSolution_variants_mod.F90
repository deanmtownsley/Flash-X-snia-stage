!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! @file

!> A module that encapsulates all compilable variants of the PUD-developed
!! static Fortran routines that are finalized by the CODE GENERATOR for use in
!! hydro advance solution task functions.  They have been "finalized" by the
!! CODE GENERATOR in the sense that macros in the original code have been
!! replaced with the appropriate target-/problem-/platform-specific keys.
!! In addition, the CODE GENERATOR will have effectively removed Flash-X
!! offloading hints and directives in CPU variants; replaced with concrete,
!! compilable offloading directives for GPU variants.
!!
!! For the purposes of this example, we are imagining that the recipe for our
!! simulation has stipulated that the solution will be advanced in terms of
!! hydrodynamics effects using only the CPU.
!!
!! This code should be written entirely by the CODE GENERATOR and, under normal
!! coding circumstances, would be derived from the
!! PUD-developed hy_* code that resides in
!! Hydro/HydroMain/simpleUnsplit.  Note that while the PUD-developed
!! code can be made private to the Hydro unit since it has previously been
!! called only by Hydro.F90, the routines in this module are presently called by
!! task functions owned by the Driver unit.  Therefore, they are in the public
!! interface.
!!
!! Note that the contents of this module could be generally useful in the sense
!! that they could be run by code other than the Orchestration unit.
!!
!! @todo Add in all metadata related to the files, the code generator, and the
!!       Simulation that this was created for.  This could include, for example,
!!         - timestamp of creation
!!         - git commit and state of repo at time of generation
!!         - state of PUD-developed files in repo and diff if altered
!!         - machine on which run and username
!!         - name/version of code generator
!!         - setup call
!!         - name/version of recipe that these are to be used in
!! @todo Is this a good design?  Main driver is likely to keep CODE GENERATORs
!!       simple and to have products to generally-useful.  For instance, we
!!       don't want the generated code to be too tightly-linked to Flash-X.
!! @todo OK to put the variants of more than one PUD-developed routine into this
!!       module?  If so, then a variants module is not related to a single
!!       routine, but rather to a single class of task function (e.g. a class
!!       containing the CPU and GPU variants of the same task function).  Is
!!       that sensible?  Will all variants of a routine only be used in a single
!!       class of task functions or across different classes of task functions?
module Hydro_advanceSolution_variants_mod
    implicit none
    private

    public :: Hydro_computeSoundSpeed_block_cpu
    public :: Hydro_computeFluxes_X_block_cpu
    public :: Hydro_computeFluxes_Y_block_cpu
    public :: Hydro_computeFluxes_Z_block_cpu
    public :: Hydro_updateSolution_block_cpu

contains

#include "Hydro_computeSoundSpeed_block_cpu.F90"
#include "Hydro_computeFluxes_X_block_cpu.F90"
#include "Hydro_computeFluxes_Y_block_cpu.F90"
#include "Hydro_computeFluxes_Z_block_cpu.F90"
#include "Hydro_updateSolution_block_cpu.F90"

end module Hydro_advanceSolution_variants_mod

