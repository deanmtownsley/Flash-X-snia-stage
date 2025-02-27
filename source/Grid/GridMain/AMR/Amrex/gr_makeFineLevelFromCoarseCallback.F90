!!****if* source/Grid/GridMain/AMR/Amrex/gr_makeFineLevelFromCoarseCallback
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
!!  gr_makeFineLevelFromCoarseCallback
!!
!! SYNOPSIS
!!  gr_makeFineLevelFromCoarseCallback(integer(IN)    :: lev,
!!                                     amrex_real(IN) :: time,
!!                                     c_ptr(IN)      :: pba,
!!                                     c_ptr(IN)      :: pdm)
!!
!! DESCRIPTION
!!  This routine is a callback routine that is registered with the AMReX AMR
!!  core at initialization.  AMReX calls this routine to populate a
!!  newly-created level with data prolongated from parent blocks and possibly
!!  with boundary data.  Prolongation is accomplished with AMReX's conservative
!!  linear interpolation routine.
!!
!!  The AMReX fill patch routines called by this routine might trigger calls
!!  to gr_fillPhysicalBC, which requires that data be in primitive form.
!!  Therefore, this routine assumes that no primitive-to-conservative
!!  transformations have been done.  AMReX will use a FLASH-provided routine to
!!  convert primitive form data to conservative form data before interpolation.
!!  It will use a different FLASH-provided routine to convert interpolated
!!  conservative form data to primitive form data as well.
!!
!!  Upon returning, the remade multifab will have data in all interiors 
!!  as well as guardcells.  EoS is not run on blocks so that it is the 
!!  responsibility of the code that triggered this callback to execute
!!  EoS where needed.
!!
!!  This routine should only be invoked by AMReX.
!!
!! ARGUMENTS
!!  lev - a 0-based number identifying the refinement level to create.  The
!!        zeroth level is the coarsest level to be used in the simulation and a
!!        larger integer indicates a finer refinement.
!!  time - IGNORED
!!  pba - a C pointer to the AMReX box array object to use for constructing the
!!        multifab for the given level
!!  pdm - a C pointer to the AMReX distribution mapping of boxes across
!!        processors to be used for constructing the multifab for the given
!!        level.
!!
!! SEE ALSO
!!  gr_preinterpolationWork
!!  gr_postinterpolationWork
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine gr_makeFineLevelFromCoarseCallback(lev, time, pba, pdm) bind(c)
    use iso_c_binding
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_amrcore_module,      ONLY : amrex_ref_ratio, &
                                          amrex_geom
    use amrex_boxarray_module,     ONLY : amrex_boxarray
    use amrex_distromap_module,    ONLY : amrex_distromap
    use amrex_multifab_module,     ONLY : amrex_multifab_build, &
                                          amrex_mfiter, &
                                          amrex_mfiter_build, &
                                          amrex_mfiter_destroy
    use gr_fluxregister_mod,       ONLY : gr_fluxregisterBuild
    use amrex_fillpatch_module,    ONLY : amrex_fillcoarsepatch

    use Grid_data,                 ONLY : gr_doFluxCorrection, &
                                          gr_amrexDidRefinement, &
                                          gr_interpolator, gr_interpolatorFace, &
                                          lo_bc_amrex, hi_bc_amrex, &
                                          lo_bc_amrexFace, hi_bc_amrexFace, &
                                          gr_meshMe
    use gr_amrexInterface,         ONLY : gr_clearLevelCallback, &
                                          gr_preinterpolationWork, &
                                          gr_postinterpolationWork, &
                                          gr_fillPhysicalBC
    use gr_physicalMultifabs,      ONLY : unk, &
                                          gr_scratchCtr, &
                                          facevars, &
                                          fluxes, &
                                          flux_registers

    implicit none

    integer,     intent(IN), value :: lev
    real(wp),    intent(IN), value :: time
    type(c_ptr), intent(IN), value :: pba
    type(c_ptr), intent(IN), value :: pdm

    type(amrex_boxarray)  :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter)    :: mfi

    integer :: nFab

    integer :: dir
    logical :: nodal(1:MDIM)

    ! Communicate to Grid_updateRefinement that we have created a level
    gr_amrexDidRefinement = .TRUE.

    ba = pba
    dm = pdm

    !!!!! CREATE MULTIFABS FOR STORING PHYSICAL DATA AT GIVEN LEVEL
    call gr_clearLevelCallback(lev)

    ! Cell-centered unknowns
    call amrex_multifab_build(unk(lev), ba, dm, NUNK_VARS, NGUARD)

#if NFACE_VARS > 0
    ! Face variables
    nodal(:)     = .FALSE.
    nodal(IAXIS) = .TRUE.
    call amrex_multifab_build(facevars(IAXIS, lev), ba, dm, NFACE_VARS, NGUARD, nodal)
#if NDIM >= 2
    nodal(:)     = .FALSE.
    nodal(JAXIS) = .TRUE.
    call amrex_multifab_build(facevars(JAXIS, lev), ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#if NDIM == 3
    nodal(:)     = .FALSE.
    nodal(KAXIS) = .TRUE.
    call amrex_multifab_build(facevars(KAXIS, lev), ba, dm, NFACE_VARS, NGUARD, nodal)
#endif
#endif

#if NSCRATCH_CENTER_VARS > 0
    call amrex_multifab_build(gr_scratchCtr(lev), ba, dm, NSCRATCH_CENTER_VARS, 0)
#endif

#if NFLUXES > 0
# ifdef USE_LEVELWIDE_FLUXES
    ! No need to store fluxes for guardcells
    do dir = 1, SIZE(fluxes, 2)
        nodal(:)   = .FALSE.
        nodal(dir) = .TRUE.
        call amrex_multifab_build(fluxes(lev, dir), ba, dm, NFLUXES, 0, &
                                  nodal=nodal)
    end do
# endif

    if ((lev > 0) .AND. (gr_doFluxCorrection)) then
#     ifdef USE_AMREX_FLASHFLUXREGISTER
        call gr_fluxregisterBuild(flux_registers(lev), fba=ba, cba=unk(lev-1)%ba, fdm=dm, cdm=unk(lev-1)%dm, &
                                      fine_lev=lev, ncomp=NFLUXES)
#     else
        call gr_fluxregisterBuild(flux_registers(lev), ba, dm, &
                                      amrex_ref_ratio(lev-1), &
                                      lev, NFLUXES)
#     endif
    end if
#endif

    !!!!!----- Fill new refinement level via interpolation from parent block
    ! This *hopefully* will do the guard cell fill as well
    ! NOTE: FLASH does not use sub-cycling (temporal interpolation)
    !
    ! -1 because Fortran variable index starts with 1
    call amrex_fillcoarsepatch(unk(lev), time,     unk(lev-1),  &
                                         time+0.1, unk(lev-1),  &
                                         amrex_geom(lev-1), gr_fillPhysicalBC,  &
                                         amrex_geom(lev  ), gr_fillPhysicalBC,  &
                                         time, &
                                         UNK_VARS_BEGIN, UNK_VARS_BEGIN, NUNK_VARS, &
                                         amrex_ref_ratio(lev-1), gr_interpolator, &
                                         lo_bc_amrex, hi_bc_amrex, &
                                         gr_preinterpolationWork, &
                                         gr_postinterpolationWork)

#if NFACE_VARS > 0
     call amrex_fillcoarsepatch(facevars(:, lev), time,     facevars(:, lev-1), &
                                                  time+0.1, facevars(:, lev-1), &
                                amrex_geom(lev-1), &
                                gr_fillPhysicalBC, gr_fillPhysicalBC, &
#if NDIM == MDIM
        &                       gr_fillPhysicalBC, &
#endif
                                amrex_geom(lev  ), &
                                gr_fillPhysicalBC, gr_fillPhysicalBC, &
#if NDIM == MDIM
        &                       gr_fillPhysicalBC, &
#endif
                             time, 1, 1, NFACE_VARS, &
                             amrex_ref_ratio(lev-1), gr_interpolatorFace, &
                             lo_bc_amrexFace, hi_bc_amrexFace)
#endif

    nFab = 0
    call amrex_mfiter_build(mfi, unk(lev), tiling=.false.)
    do while(mfi%next())
        nFab = nFab + 1 
    end do
    call amrex_mfiter_destroy(mfi)

    if (gr_meshMe == MASTER_PE) then
       write(*,'(A,I0,A,I0,A)') "Made fine level ", lev + 1, " - ", nFab, " blocks"
    end if

end subroutine gr_makeFineLevelFromCoarseCallback

