!!****if* source/Grid/GridMain/paramesh/Grid_init
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
!!  Grid_init
!!
!! SYNOPSIS
!!
!!  call Grid_init()
!!
!! DESCRIPTION
!!  Initialize Grid_data
!!
!! ARGUMENTS
!!  none
!!
!! PARAMETERS 
!!
!!  nblockx [INTEGER] 
!!     num initial blocks in x dir
!!  nblocky [INTEGER] 
!!     num initial blocks in y dir   
!!  nblockz [INTEGER] 
!!     num initial blocks in z dir   
!!  lrefine_max [INTEGER] 
!!      maximum AMR refinement level
!!  lrefine_min [INTEGER] 
!!      minimum AMR refinement level
!!  nrefs [INTEGER] 
!!      refine/derefine AMR grid every nrefs timesteps
!!
!!  refine_var_1 [INTEGER] 
!!     indicates first variable on which to refine
!!  refine_cutoff_1 [REAL] 
!!      threshold value to trigger refinement for refine_var_1
!!  derefine_cutoff_1 [REAL]
!!      threshold value to trigger derefinement for refine_var_1
!!  refine_filter_1 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_1
!!
!!  refine_var_2 [INTEGER] 
!!     indicates second variable on which to refine
!!  refine_cutoff_2 [REAL] 
!!      threshold value to trigger refinement for refine_var_2
!!  derefine_cutoff_2 [REAL]
!!      threshold value to trigger derefinement for refine_var_2
!!  refine_filter_2 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_2
!!
!!  refine_var_3 [INTEGER] 
!!     indicates third variable on which to refine (if needed)
!!  refine_cutoff_3 [REAL] 
!!      threshold value to trigger refinement for refine_var_3
!!  derefine_cutoff_3 [REAL]
!!      threshold value to trigger derefinement for refine_var_3
!!  refine_filter_3 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_3
!!
!!  refine_var_4 [INTEGER] 
!!     indicates fourth variable on which to refine (if needed)
!!  refine_cutoff_4 [REAL] 
!!      threshold value to trigger refinement for refine_var_4
!!  derefine_cutoff_4 [REAL]
!!      threshold value to trigger derefinement for refine_var_4
!!  refine_filter_4 [REAL]
!!      prevents error calculations from diverging numerically for refine_var_4
!!
!!  flux_correct [BOOLEAN]
!!     turns on or off flux correction
!! small  [REAL]
!!   Generic small value that can be used as floor where needed
!! smlrho [REAL]  
!!   Cutoff value for density    
!! smallp [REAL]  
!!   Cutoff value for pressure
!! smalle [REAL]  
!!   Cutoff value for energy
!! smallt [REAL]  
!!   Cutoff value for temperature
!! smallu [REAL]  
!!   Cutoff value for velocity
!! smallx [REAL]  
!!   Cutoff value for abundances
!! eosMode[STRING]
!!   determines which variables to calculate from the ones
!!   defined. Possible values are "dens_ie", "dens_pres" and "dens_temp"
!! interpol_order [INTEGER]
!!   Order of interpolation, used in Paramesh2 "monotonic" interpolation
!!   for mesh prolongation
!! grid_monotone_hack [BOOLEAN]
!!   If .true., apply radical monotonicity constraints to interpolants,
!!   i.e., completely flatten them if they violate monotonicity.  Used
!!   in Paramesh2 "quadratic_cartesian" interpolation for mesh prolongation.
!! earlyBlockDistAdjustment [BOOLEAN]
!!   If .true., let Paramesh redistribute blocks
!!   across processors early, so that the block distribution chosen by
!!   Paramesh will be in effect when time evolution begins after restart.
!!   If earlyBlockDistAdjustment is .false., the block distribution enacted
!!   by the IO unit when it read a checkpoint file will normally still be
!!   in effect when time evolution begins after a restart.
!!   This flag is ignored if not restarting from a checkpoint.
!! 
!! lrefine_del [INTEGER]
!! gr_lrefineMaxRedDoByTime [BOOLEAN]
!! gr_lrefineMaxRedTRef [REAL]
!! gr_lrefineMaxRedTimeScale [REAL]
!! gr_lrefineMaxRedLogBase [REAL]
!!
!! gr_lrefineMaxRedDoByLogR [BOOLEAN]
!! gr_lrefineMaxRedRadiusFact [REAL]
!! x_refine_center [REAL]
!! y_refine_center [REAL]
!! z_refine_center [REAL]
!!
!! gr_restrictAllMethod [INTEGER]
!!***

subroutine Grid_init()

  use Grid_data
  use gr_ptInterface,   ONLY : gr_ptInit
  use gr_tilePolicyData,ONLY : gr_initTilePolicy
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Driver_interface, ONLY : Driver_abort, Driver_getMype, &
    Driver_getNumProcs, Driver_getComm
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Simulation_interface, ONLY : Simulation_mapStrToInt, Simulation_getVarnameType
  use Grid_interface, only: Grid_getVarNonRep

!  use gr_sbInterface, ONLY: gr_sbInit
  implicit none

#include "Simulation.h"
#include "constants.h"

!!$  logical :: useProtonEmission
!!$  logical :: useProtonImaging
!!$  logical :: useElectronSpectrometry

  integer :: i, j, k, localNumBlocks, ii, numLeafBlks

  character(len=MAX_STRING_LENGTH),save :: refVarname,refVarString,paramString
  character(len=MAX_STRING_LENGTH),save :: refCutoffName,refCutOffString
  character(len=MAX_STRING_LENGTH),save :: derefCutoffName,derefCutOffString
  character(len=MAX_STRING_LENGTH),save :: refFiltername,refFilterString
  character(len=MAX_STRING_LENGTH) :: xl_bcString,xr_bcString
  character(len=MAX_STRING_LENGTH) :: yl_bcString,yr_bcString
  character(len=MAX_STRING_LENGTH) :: zl_bcString,zr_bcString
  character(len=MAX_STRING_LENGTH) :: eosModeString, grav_boundary_type
  integer,save :: refVar
  integer :: countInComm, color, key, ierr
  integer :: nonrep


  call Driver_getMype(GLOBAL_COMM, gr_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, gr_globalNumProcs)
  call Driver_getComm(GLOBAL_COMM, gr_globalComm)

  call Driver_getMype(MESH_COMM, gr_meshMe)
  call Driver_getNumProcs(MESH_COMM, gr_meshNumProcs)
  call Driver_getComm(MESH_COMM, gr_meshComm)

  call Driver_getMype(MESH_ACROSS_COMM, gr_meshAcrossMe)
  call Driver_getNumProcs(MESH_ACROSS_COMM, gr_meshAcrossNumProcs)
  call Driver_getComm(MESH_ACROSS_COMM, gr_meshAcrossComm)

  ! DO THIS EARLY - must be before gr_initGeometry is called:
  !get the boundary conditions stored as strings in the flash.par file
  call RuntimeParameters_get("xl_boundary_type", xl_bcString)
  call RuntimeParameters_get("xr_boundary_type", xr_bcString)
  call RuntimeParameters_get("yl_boundary_type", yl_bcString)
  call RuntimeParameters_get("yr_boundary_type", yr_bcString)
  call RuntimeParameters_get("zl_boundary_type", zl_bcString)
  call RuntimeParameters_get("zr_boundary_type", zr_bcString)

  !map the string boundary conditions to integer constants defined in constants.h
  call RuntimeParameters_mapStrToInt(xl_bcString,gr_domainBC(LOW,IAXIS))
  call RuntimeParameters_mapStrToInt(xr_bcString,gr_domainBC(HIGH,IAXIS))
  call RuntimeParameters_mapStrToInt(yl_bcString,gr_domainBC(LOW,JAXIS))
  call RuntimeParameters_mapStrToInt(yr_bcString,gr_domainBC(HIGH,JAXIS))
  call RuntimeParameters_mapStrToInt(zl_bcString,gr_domainBC(LOW,KAXIS))
  call RuntimeParameters_mapStrToInt(zr_bcString,gr_domainBC(HIGH,KAXIS))

!----------------------------------------------------------------------------------
! mesh geometry - done early so Paramesh_init can use gr_geometry for some checking
!----------------------------------------------------------------------------------
  ! Initialization of gr_geometry etc is done in gr_initGeometry.

  call gr_initGeometry()

#ifdef GRID_WITH_MONOTONIC
  if (NGUARD < 4) then
     if (gr_meshMe==MASTER_PE) then
        print*,'Grid_init: Monotonic grid interpolation requires at least 4 layers of guard cells.'
        print*,' However, NGUARD is only ', NGUARD
        print*," Maybe you want to setup with '-gridinterpolation=native',"
        print*," or make sure that NGUARD is set correctly in Config file."
        call Driver_abort("Please setup with '-gridinterpolation=native', or change NGUARD.")
     end if
  endif
#endif

  ! Initialization of gr_imin,gr_imax,...,gr_kmax is done in gr_initGeometry, called above.

  call RuntimeParameters_get('lrefine_min', gr_minRefine)
  call RuntimeParameters_get('lrefine_max', gr_lrefineMax)
  call RuntimeParameters_get("nrefs", gr_nrefs)
  call RuntimeParameters_get('lrefine_min_init', gr_lrefineMinInit)

  call RuntimeParameters_get("smalle",gr_smalle)
  call RuntimeParameters_get("smlrho",gr_smallrho)
  call RuntimeParameters_get("smallx",gr_smallx) !
!  call RuntimeParameters_get("grid_monotone_hack", gr_monotone) ! for "quadratic_cartesian" interpolation
  call RuntimeParameters_get("interpol_order",gr_intpol) ! for "monotonic" interpolation

  call RuntimeParameters_get("flux_correct", gr_doFluxCorrection)

  call RuntimeParameters_get("gr_enableTiling", gr_enableTiling)
  call RuntimeParameters_get("gr_useTiling", gr_useTiling)
  if (gr_useTiling .AND. .NOT. gr_enableTiling) then
     if (gr_meshMe == MASTER_PE) then
        write(*,*) "[Grid_init] WARNING: Tiling is disabled by gr_enableTiling=.FALSE."
        write(*,*) "                     Therefore gr_useTiling will not take effect."
     end if
  end if
  call RuntimeParameters_get("gr_tileSizeX", gr_tileSize(IAXIS))
  call RuntimeParameters_get("gr_tileSizeY", gr_tileSize(JAXIS))
  call RuntimeParameters_get("gr_tileSizeZ", gr_tileSize(KAXIS))

#ifdef FLASH_PARTICLES
  call RuntimeParameters_get('useParticles',gr_useParticles)
  call RuntimeParameters_get('pt_maxPerProc',gr_maxParticlesPerProc)
#else
  gr_useParticles=.false.
#endif

!!$#ifdef FLASH_EDEP
!!$  call RuntimeParameters_get('useEnergyDeposition', gr_useEnergyDeposition)
!!$  gr_useParticles = gr_useEnergyDeposition
!!$#else
!!$  gr_useEnergyDeposition = .false.
!!$#endif

!!$#ifdef FLASH_GRID_PARTICLES
!!$  call RuntimeParameters_get('useProtonImaging',useProtonImaging)
!!$  if (useProtonImaging) then
!!$      gr_useParticles=.true.
!!$  end if
!!$  call RuntimeParameters_get('useProtonEmission',useProtonEmission)
!!$  if (useProtonEmission) then
!!$      gr_useParticles=.true.
!!$  end if
!!$!  call RuntimeParameters_get('useElectronSpectrometry',useElectronSpectrometry)
!!$!  if (useElectronSpectrometry) then
!!$!      gr_useParticles=.true.
!!$!  end if
!!$#endif

  call RuntimeParameters_get('useOrchestration',gr_useOrchestration)

  !Check if there are gravitational isolated boundary conditions
  !in order to determine which solvers to intialize.
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  gr_isolatedBoundaries = (grav_boundary_type=="isolated")

  gr_allPeriodic = .TRUE.
  do i = 1, NDIM
     if(gr_domainBC(LOW,  i) /= PERIODIC)     gr_allPeriodic = .FALSE.
     if(gr_domainBC(HIGH, i) /= PERIODIC)     gr_allPeriodic = .FALSE.
  end do

  call RuntimeParameters_get("bndPriorityOne",gr_bndOrder(1))
  call RuntimeParameters_get("bndPriorityTwo",gr_bndOrder(2))
  call RuntimeParameters_get("bndPriorityThree",gr_bndOrder(3))

  ! Initialization of gr_globalDomain, containing the same information as gr_imin,...,gr_kmax,
  ! is done in gr_initGeometry, called above.


  call RuntimeParameters_get("eosMode", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosMode)

  call RuntimeParameters_get("eosModeInit", eosModeString)
  call RuntimeParameters_mapStrToInt(eosModeString, gr_eosModeInit)

  gr_eosModeNow = gr_eosModeInit ! may change after initialization is done

  call RuntimeParameters_get("earlyBlockDistAdjustment", gr_earlyBlockDistAdjustment)
  gr_justExchangedGC = .false.

  call RuntimeParameters_get("refine_on_particle_count",gr_refineOnParticleCount)

  call RuntimeParameters_get("min_particles_per_blk",gr_minParticlesPerBlk)
  call RuntimeParameters_get("max_particles_per_blk",gr_maxParticlesPerBlk)


  call gr_initSpecific()

  call RuntimeParameters_get("convertToConsvdForMeshCalls", gr_convertToConsvdForMeshCalls)
  call RuntimeParameters_get("convertToConsvdInMeshInterp", gr_convertToConsvdInMeshInterp)
  if (gr_convertToConsvdInMeshInterp) then
     if (gr_convertToConsvdForMeshCalls) then
        ! For PARAMESH 4, if both ways of conversion to conserved form are requested,
        ! Let the new mechanism win and try to make sure the old one is not used. - KW
        if(gr_meshMe == MASTER_PE) &
             print*,'WARNING: convertToConsvdForMeshCalls ignored since convertToConsvdInMeshInterp is requested'
        call Logfile_stampMessage( &
             "WARNING: convertToConsvdForMeshCalls ignored since convertToConsvdInMeshInterp is requested")
        gr_convertToConsvdForMeshCalls = .FALSE.
     end if
  end if

#ifdef GRID_WITH_MONOTONIC
  gr_intpolStencilWidth = 2     !Could possibly be less if gr_intpol < 2  - KW
#endif

  !! This section of the code identifies the variables to used in
  !! the refinement criterion. If a variable is a refinement variable
  !! then the corresponding refinement/derefinement cutoff and filter
  !! values also have to be fetched. The config file defines
  !! refinement variables as strings, names as "refine_var_1",
  !! "refine_var_2" etc, with the current maximum being 4. The
  !! general utility routine takes the base "refine_var_" and appends
  !! the index at the end of the string to generate the parameter
  !! name and the routine Simulation_mapStrToInt finds its index into UNK.

  call RuntimeParameters_get("refine_var_count",gr_numRefineVarsMax)
  gr_refine_var = NONEXISTENT
  gr_numRefineVars=0

  refVarName='refine_var_'
  refCutoffName='refine_cutoff_'
  derefCutoffName='derefine_cutoff_'
  refFilterName='refine_filter_'

  do i = 1,gr_numRefineVarsMax
     call concatStringWithInt(refVarName,i,refVarString)
     call RuntimeParameters_get( refVarString, paramString)
     if(paramString /= "none") then
        do ! not a real loop
           call Simulation_mapStrToInt(paramString, refVar, MAPBLOCK_UNK)
           if(refVar <= 0) exit
           call Grid_getVarNonRep(MAPBLOCK_UNK, refVar, nonrep)
           if(nonrep > 0) then; refVar = 0; exit; end if

           gr_numRefineVars=gr_numRefineVars+1
           gr_refine_var(gr_numRefineVars)=refVar
           call concatStringWithInt(refCutoffName,gr_numRefineVars,refCutoffString)
           call concatStringWithInt(derefCutoffName,gr_numRefineVars,derefCutOffString)
           call concatStringWithInt(refFilterName,gr_numRefineVars,refFilterString)
           call RuntimeParameters_get( refCutoffString, gr_refine_cutoff(gr_numRefineVars)  )
           call RuntimeParameters_get( derefCutoffString, gr_derefine_cutoff(gr_numRefineVars) )
           call RuntimeParameters_get( refFilterString,  gr_refine_filter(gr_numRefineVars) )
           exit ! told you it wasnt a real loop
        end do
        if(refVar <= 0) then
           if(gr_globalMe == MASTER_PE) &
              print*, 'WARNING: Unrecognized or non-replicated variable name in refine_var_',i,' treating it as "none"'
           call Logfile_stampMessage( &
              'WARNING: Unrecognized or non-replicatedvariable name in refine_var, treating it as "none"')
        end if
     end if
  end do

  gr_enforceMaxRefinement = .FALSE.

  call RuntimeParameters_get("lrefine_del", gr_lrefineDel)


  call RuntimeParameters_get("gr_lrefineMaxRedDoByLogR", gr_lrefineMaxRedDoByLogR)
  call RuntimeParameters_get("gr_lrefineMaxRedRadiusFact", gr_lrefineMaxRedRadiusSq)
  gr_lrefineMaxRedRadiusSq = gr_lrefineMaxRedRadiusSq * gr_lrefineMaxRedRadiusSq
  call RuntimeParameters_get("x_refine_center", gr_lrefineCenterI)
  call RuntimeParameters_get("y_refine_center", gr_lrefineCenterJ)
  call RuntimeParameters_get("z_refine_center", gr_lrefineCenterK)

  call RuntimeParameters_get("gr_lrefineMaxRedDoByTime", gr_lrefineMaxRedDoByTime)
  if (gr_lrefineMaxRedDoByTime) gr_enforceMaxRefinement = .TRUE.
  call RuntimeParameters_get("gr_lrefineMaxRedTimeScale", gr_lrefineMaxRedTimeScale)
  call RuntimeParameters_get("gr_lrefineMaxRedLogBase", gr_lrefineMaxRedLogBase)
  call RuntimeParameters_get("gr_lrefineMaxRedTRef", gr_lrefineMaxRedTRef)

  call RuntimeParameters_get("gr_lrefineMaxByTime", gr_lrefineMaxByTime)
  if (gr_lrefineMaxByTime) then
     gr_enforceMaxRefinement = .TRUE.

     do i = 1, GR_LREFMAXTIMES
        refVarName = "gr_lrefmaxTime_"
        call concatStringWithInt(refVarName,i,refVarString)
        call RuntimeParameters_get(refVarString, gr_lrefmaxTimes(i))

        refVarName = "gr_lrefmaxTimeValue_"
        call concatStringWithInt(refVarName,i,refVarString)
        call RuntimeParameters_get(refVarString, gr_lrefmaxTimeValues(i))

        if(i > 1) then 
           if(gr_lrefmaxTimes(i) > 0.0 .and. &
                gr_lrefmaxTimes(i) < gr_lrefmaxTimes(i-1)) then
              if(gr_meshMe == MASTER_PE) then
                 call Driver_abort('[Grid_init] Custom lrefine_max times must be in increasing order')
              end if
           end if
        end if

        if(gr_lrefmaxTimes(i) > 0.0 .and. gr_lrefmaxTimeValues(i) < 0) then
           call Driver_abort('[Grid_init] If custom lrefine_max time is set, the value must also be set')
        end if
     end do
  end if

  if(gr_numRefineVars==0)then
     if(gr_meshMe == MASTER_PE) print*,'WARNING : Adaptive Grid did not find any refinement variables'
     call Logfile_stampMessage("WARNING : Adaptive Grid did not find any variable to refine")
  end if

  call RuntimeParameters_get("gr_restrictAllMethod", gr_restrictAllMethod)

  gr_anyVarToConvert = .FALSE.
  do i = UNK_VARS_BEGIN,UNK_VARS_END
     call Simulation_getVarnameType(i, gr_vartypes(i))
     if (gr_vartypes(i) .eq. VARTYPE_PER_MASS) gr_anyVarToConvert = .TRUE.
  end do




  !Only call the particle initialization routines when
  !we are using particles.
  if (gr_useParticles .eqv. .true. ) then
     call gr_ptInit()
     call gr_ptMapInit()
  endif

!  call gr_sbInit()

    ! Reduce guard cell fills
  call RuntimeParameters_get ("reduceGcellFills", gr_reduceGcellFills)

  gr_region=0.0
 
  call gr_initTilePolicy

end subroutine Grid_init
