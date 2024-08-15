#include "constants.h"
#include "Simulation.h"

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
!! A standard Flash-X Fortran module that encapsulates all Grid-private variables
!! as needed by Milhoja.  Note that Grid_init and auxiliary initialization routines
!! called by Grid_init effectively act as the GridMain constructor.  Therefore, all
!! variables in this module should be explicitly initialized once Grid_init is
!! called.
!!
!! @todo gr_numDataStruct, gr_gridDataStruct, gr_gridDataStructSize, and
!! gr_bcEnableApplyMixedGds never initialized.  Does Grid_init need to
!! initialize them with gr_setDataStructInfo?  Only gr_bcApplyToAllBlks seems to
!! use them and this presently aborts if called with Milhoja.  Should gr_bcEnableApplyMixedGds
!! be set in Grid_init prior to calling gr_setDataStructInfo so that it is
!! always explicitly initialized?  It looks like that routine should just
!! use the value rather than set it.
!! @todo Is gr_allPeriodic necessary?
!! @todo Is gr_numRefineVarsMax necessary?
!! @todo Will Milhoja ever be used with Monotonic?  Do we need
!!       gr_intpolStencilWidth?
module Grid_data
    implicit none
    public
    
    !!!!!----- MPI & OpenMP Information    
    integer, save :: gr_globalComm
    integer, save :: gr_globalMe
    integer, save :: gr_globalNumProcs
    
    integer, save :: gr_meshComm
    integer, save :: gr_meshMe
    integer, save :: gr_meshNumProcs
    
    integer, save :: gr_meshAcrossComm
    integer, save :: gr_meshAcrossMe
    integer, save :: gr_meshAcrossNumProcs

    logical, save :: gr_useOrchestration
    integer, save :: gr_envOmpNumThreads

    real,    save :: gr_smallx
    real,    save :: gr_smalle
    real,    save :: gr_smallrho

    integer, save :: gr_numDataStruct
    integer, save :: gr_gridDataStruct(NDATATYPES)
    integer, save :: gr_gridDataStructSize(NDATATYPES)
    logical, save :: gr_bcEnableApplyMixedGds

    integer, save :: gr_lRefineMax
    integer, save :: gr_maxRefine
    integer, save :: gr_geometry
    real,    save :: gr_minCellSize
    real,    save :: gr_minCellSizes(1:MDIM)
    logical, save :: gr_allPeriodic

    integer, save :: gr_dirGeom(1:MDIM)
    logical, save :: gr_dirIsAngular(1:MDIM)
    integer, save :: gr_domainBC(LOW:HIGH, 1:MDIM)
    real,    save :: gr_globalDomain(LOW:HIGH, 1:MDIM)

    integer, save :: gr_vartypes(UNK_VARS_BEGIN:UNK_VARS_END)

    integer, save :: gr_numRefineVars
    integer, save :: gr_numRefineVarsMax
    integer, save :: gr_refine_var(MAXREFVARS)
    real,    save :: gr_refine_cutoff(MAXREFVARS)
    real,    save :: gr_derefine_cutoff(MAXREFVARS)
    real,    save :: gr_refine_filter(MAXREFVARS)
    logical, save :: gr_enforceMaxRefinement

    logical, save :: gr_enableTiling
    logical, save :: gr_useTiling
    integer, save :: gr_tileSize(1:MDIM)

    logical, save :: gr_geometryOverride
    character(len=MAX_STRING_LENGTH), save :: gr_str_geometry

    integer, save :: gr_eosMode
    integer, save :: gr_eosModeInit

    logical, save :: gr_justExchangedGC
    logical, save :: gr_gcellsUpToDate

    ! Copied from AMReX
    integer, parameter :: gr_intpolStencilWidth = 1
end module Grid_data

