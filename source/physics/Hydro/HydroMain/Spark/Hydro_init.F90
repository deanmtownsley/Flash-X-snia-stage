!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_init
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
!!  Hydro_init
!!
!!
!!  For more details see the documentation of the NULL implementation
!!
!!
!!***

subroutine Hydro_init()

  use Hydro_data
  use Driver_interface,            ONLY : Driver_abort, Driver_getMype, &
       Driver_getNumProcs,    &
       Driver_getComm
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
       RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface,           ONLY : Logfile_stampMessage, &
       Logfile_stampVarMask, &
       Logfile_stamp
  use Grid_interface,              ONLY : Grid_setFluxHandling
  use IO_interface,                ONLY : IO_getScalar, IO_setScalar

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  character(len=MAX_STRING_LENGTH) :: str_geometry
  integer :: i
  logical :: threadBlockListBuild, threadWithinBlockBuild
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  ! Set allocation flag to false. This will allow the scratch array to only be allocated once.
  scratch_allocated = .False.

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,hy_meshMe)
  call Driver_getNumProcs(MESH_COMM,hy_meshNumProcs)
  call Driver_getComm(GLOBAL_COMM,hy_globalComm)

  call RuntimeParameters_get("useHydro",hy_useHydro) !!DEV: not yet systematically honored

  call RuntimeParameters_get("cfl", hy_cfl)

  call RuntimeParameters_get("restart", hy_restart)
  if (hy_restart) then
     call IO_getScalar("dt", hy_dt)
     hy_dtmin = hy_dt
     call IO_getScalar("C_hyp", hy_lChyp)
  else
     hy_dtmin = HUGE(1.0)
     ! Initialize local maximum hyperbolic speed
     hy_lChyp = TINY(1.0)
     call IO_setScalar("C_hyp", hy_lChyp)
  endif

#ifdef FLASH_GRID_PARAMESH
  if ((NGUARD > 4) .and. (NXB < 2*NGUARD)) then
     call Driver_abort&
          ("[Hydro_init]: Hydro requires larger NXB, etc. for the given number of guardcells.")
  endif
#endif

  !! Need to issue a warning or error if PM4dev is not used

  call RuntimeParameters_get("smlrho",              hy_smalldens)
  call RuntimeParameters_get("smallp",              hy_smallpres)
  call RuntimeParameters_get("smallE",              hy_smallE)
  call RuntimeParameters_get("smallx",              hy_smallX)
  call RuntimeParameters_get("smallu",              hy_smallu)

  call RuntimeParameters_get("updateHydroFluxes",   hy_updateHydrofluxes)
  call RuntimeParameters_get("use_hybridRiemann",   hy_hybridRiemann)
  call RuntimeParameters_get("use_flattening",      hy_flattening)

  call RuntimeParameters_get("cvisc",               hy_cvisc)
  call RuntimeParameters_get("lim_rad",             hy_limRad)

  call RuntimeParameters_get("alpha_glm",           hy_alphaGLM)

  !! Geometry ------------------------------------------------------------------
  call RuntimeParameters_get("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, hy_geometry)
  if (hy_geometry .NE. CARTESIAN .AND. hy_meshME == MASTER_PE )  then
     print *, "[Hydro_init]: Using non-Cartesian Geometry!"
  endif
  call RuntimeParameters_get("flux_correct", hy_fluxCorrect)
  call RuntimeParameters_get("flux_correct_perLevel", hy_fluxCorrectPerLevel)
  ! if (NDIM > 1) then
     if (hy_fluxCorrect) then
        if (hy_geometry == CARTESIAN) then
           call Grid_setFluxHandling('consv_flux_densities')
        else
           call Grid_setFluxHandling('consv_fluxes')
        endif
     end if
  ! end if
  call RuntimeParameters_get("hy_useTiling", hy_useTiling)
  !! For correct flux correction in non-Cartesian geometry----------------------
  do i = 1, NFLUXES
     hy_fluxCorVars(i) = i
  enddo

  !! Allow selective guardcell fill calls ---------------------------------------
  hy_gcMaskSize = NUNK_VARS
  hy_gcMask = .FALSE.

  ! Set guardcell masking for needed variables
  ! First, the standard set for plain Hydro
  hy_gcMask(DENS_VAR) = .TRUE.
  hy_gcMask(EINT_VAR) = .TRUE.
  hy_gcMask(ENER_VAR) = .TRUE.
  hy_gcMask(GAMC_VAR) = .TRUE.
  hy_gcMask(GAME_VAR) = .TRUE.
  hy_gcMask(PRES_VAR) = .TRUE.
  hy_gcMask(TEMP_VAR) = .TRUE.
  hy_gcMask(VELX_VAR) = .TRUE.
  hy_gcMask(VELY_VAR) = .TRUE.
  hy_gcMask(VELZ_VAR) = .TRUE.
  ! Now for the MHD variables, including FACEVARS
  ! Now for the gravity variables
#ifdef GPOL_VAR
  hy_gcMask(GPOL_VAR) = .TRUE.
#endif
#ifdef GPOT_VAR
  hy_gcMask(GPOT_VAR) = .TRUE.
#endif
#ifdef SHOK_VAR
  hy_gcMask(SHOK_VAR) = .TRUE.
#endif
#ifdef SPARK_GLM
  hy_gcMask(MAGX_VAR:MAGZ_VAR) = .TRUE.
  hy_gcMask(PSIB_VAR) = .TRUE.
#endif
  ! Fill for all species/mass scalars
  hy_gcMask(SPECIES_BEGIN:MASS_SCALARS_END) = .TRUE.

  call Logfile_stampVarMask(hy_gcMask, .FALSE., '[Hydro_init]', 'gcNeed')


  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadHydroWithinBlock", hy_threadWithinBlock)

  if (hy_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Hydro_init]')
     hy_threadWithinBlock = .false.
  end if

  if (hy_geometry == POLAR) & !Polar in 3D (that's a no no)
       call Driver_abort("[Hydro_computeDt] ERROR: Polar geometry not supported in 3D")

  call PhysicalConstants_get("Newton", hy_gravConst)
  hy_4piGinv = (4.*PI*hy_gravConst)**(-1)

  ! Deal with units for MHD. We will ALWAYS use Gaussian CGS. Jackson 2E FTW.
  call RuntimeParameters_get("unitsystem", hy_units)
  if (hy_units /= "NONE" .AND. hy_units /= "none") then
    hy_bref = sqrt(4.0*PI)
  end if

end subroutine Hydro_init
