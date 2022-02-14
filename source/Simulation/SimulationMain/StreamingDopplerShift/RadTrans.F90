!!****cr* source/Simulation/SimulationMain/StreamingDopplerShift/RadTrans
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
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans( real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!      This subroutine performs the radiatiative transfer calculation
!!      for this step. 
!!
!! ARGUMENTS
!!
!!   dt     : The time step
!!   pass   : Reverses direction of solve
!!
!! HISTORY
!!   2020     Initial RadTrans                       - Ran Chu
!!   2021     Added support for thornado ApplyBC     - Ran Chu
!!
!!***

!!REORDER(4): solnData

subroutine RadTrans( dt, pass )

  use Driver_interface, ONLY : Driver_abort
  use Eos_interface, ONLY : Eos_everywhere
  use FluidFieldsModule, ONLY : uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  use GeometryFieldsModule, ONLY : uGF, iGF_h_1, iGF_h_2, iGF_h_3
  use Grid_interface, ONLY : Grid_fillGuardCells, &
      Grid_getMaxRefinement, Grid_getDeltas, &
      Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_data, ONLY: gr_str_geometry
  use RadiationFieldsModule, ONLY : uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  use RadTrans_data, ONLY : rt_useRadTrans, rt_enableTiling, &
     rt_eosModeGc, rt_gcMask, rt_str_geometry, rt_geometry
  use rt_data, ONLY : rt_doExplicit, rt_doImplicit
  use SubcellReconstructionModule, ONLY : ReconstructionMatrix, ProjectionMatrix
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use TimeSteppingModule_Flash, ONLY : ComputeTimeStep_TwoMoment, Update_IMEX_PDARS
  use UnitsModule, ONLY : Gram, Centimeter, Second, AtomicMassUnit

  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

#include "Simulation.h"
#include "constants.h"

  real,    intent(in) :: dt
  integer, intent(in), optional :: pass

  real, pointer, dimension(:,:,:,:) :: solnData
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  integer, dimension(MDIM) :: lo,hi,u_lo,u_hi
  real, dimension(MDIM) :: dx, dx_thornado
  real :: dt_thornado, dt_sub, t_sub

  integer :: i, j, k, ic, jc, kc, ii, jj, kk, n, ioff, ivar, i_sub, n_sub
  integer :: level, maxLev

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t) :: tileDesc

  integer, parameter :: my_ngrow = 2
  integer :: nX(3), swX(3)
  real :: xL(3), xR(3)

  real :: D, S1, S2, S3, E, Nel, ekin
  integer :: iNode, iNodeX, iNodeE, iCF, iCR, iS, iE

  integer :: faces(LOW:HIGH, 1:MDIM)
  integer :: onBoundary(LOW:HIGH, 1:MDIM)
  integer :: ApplyBC

  real, parameter :: conv_x    = Centimeter
  real, parameter :: conv_dens = 1
  real, parameter :: conv_mom  = conv_dens * Centimeter/Second
  real, parameter :: conv_enr  = Gram / Centimeter / Second**2
  real, parameter :: conv_ne   = Gram / Centimeter**3 / AtomicMassUnit
  real, parameter :: conv_J    = Gram/Second**2/Centimeter
  real, parameter :: conv_H    = Gram/Second**3

  nullify(solnData)

  if (.NOT. rt_useRadTrans) return

  if ( my_ngrow /= 2 ) then
     call Driver_abort("Need two ghost cells in call_to_thornado!")
  else if ( my_ngrow*THORNADO_NNODESX > NGUARD ) then
     call Driver_abort("NGUARD must be at least my_ngrow*THORNADO_NNODESX")
  end if

  call Timers_start("RadTrans")

  call Timers_start("rt_gc")
  call Grid_fillGuardCells(CENTER,ALLDIR, &
    minLayers=my_ngrow*THORNADO_NNODESX,doEos=.true.,eosMode=rt_eosModeGc, &
    maskSize=NUNK_VARS,mask=rt_gcMask,makeMaskConsistent=.true.)
  call Timers_stop("rt_gc")

#ifdef FLASH_GRID_UG
  maxLev = 1
#else
  call Grid_getMaxRefinement(maxLev,mode=1) !mode=1 means lrefine_max, which does not change during sim.
#endif
  do level = 1, maxLev

     ! Get cell widths and timestep
     call Grid_getDeltas(level,dx)
     dx_thornado = THORNADO_NNODESX * dx
     call ComputeTimeStep_TwoMoment( dx_thornado, dt_thornado )
     dt_sub = min(dt_thornado,dt)

     t_sub = 0.0
     do while ( t_sub < dt )
        if ( t_sub + dt_sub > dt ) then
          dt_sub = dt - t_sub
          call Driver_abort("No substep in this problem.") 
        end if

        call Grid_getTileIterator(itor, LEAF, level=level, tiling=rt_enableTiling)
        do while(itor%isValid())
           call itor%currentTile(tileDesc)

           ! Get a pointer to solution data
           call tileDesc%getDataPtr(solnData, CENTER)

           ! Check if the block needs boundary condition
           ! RC: This only works for the reflecting inner boundary in 1D
           call tileDesc%faceBCs(faces)
           if( faces(LOW,IAXIS) == USER_DEFINED .and. &
               rt_geometry == CARTESIAN )then
             ApplyBC = 0
           else
             ApplyBC = 0
           end if

           ! get dimensions/limits and coordinates
           nX = 1
           swX = 0
           xL = 0.0
           if ( rt_geometry == CARTESIAN ) then
              xR = 1.0 / conv_x
           else if ( rt_geometry == CYLINDRICAL ) then
              xR = [ 1.0, 1.0, 2.0*PI ]
           else if ( rt_geometry == SPHERICAL ) then
              xR = [ 1.0, PI, 2.0*PI ]
           end if

           lo = tileDesc%limits(LOW,1:MDIM)
           hi = tileDesc%limits(HIGH,1:MDIM)
           nX(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
           swX(1:NDIM) = my_ngrow
           u_lo = 1 - swX
           u_hi = nX + swX

           call tileDesc%boundBox(boundBox)
           xL(1:NDIM) = boundBox(LOW, 1:NDIM)
           xR(1:NDIM) = boundBox(HIGH,1:NDIM)

           ! convert cm to m for Thornado
           xL(1) = xL(1) * conv_x
           xR(1) = xR(1) * conv_x
           if ( rt_geometry == CARTESIAN .or. rt_geometry == CYLINDRICAL) then
              xL(2) = xL(2) * conv_x
              xR(2) = xR(2) * conv_x
              if ( rt_geometry == CARTESIAN ) then
                 xL(3) = xL(3) * conv_x
                 xR(3) = xR(3) * conv_x
              end if
           end if

           ! Setup thornado data structures
           call Timers_start("rt_init")
           call InitThornado_Patch(nX, swX, xL, xR, THORNADO_NSPECIES, rt_str_geometry)
           call Timers_stop("rt_init")

           ! Thornado uses units where c = G = k = 1, Meter = 1
           call Timers_start("rt_reconstruction")
           do kc = u_lo(KAXIS), u_hi(KAXIS)
              do jc = u_lo(JAXIS), u_hi(JAXIS)
                 do ic = u_lo(IAXIS), u_hi(IAXIS)
                    i = lo(IAXIS) + THORNADO_NNODESX*(ic-1)
                    j = lo(JAXIS) + THORNADO_NNODESX*(jc-1)
                    k = lo(KAXIS) + THORNADO_NNODESX*(kc-1)

                    ! Interpolate fluid vars from Flash arrays into Thornado arrays
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       D   = 0.0
                       S1  = 0.0
                       S2  = 0.0
                       S3  = 0.0
                       E   = 0.0
                       Nel = 0.0
                       do n = 1, THORNADO_FLUID_NDOF
                          ii = mod((n-1)                    ,THORNADO_NNODESX) + i
                          jj = mod((n-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                          kk = mod((n-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k
                          D   = D   + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)
                          S1  = S1  + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(VELX_VAR,  ii,jj,kk)
                          S2  = S2  + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(VELY_VAR,  ii,jj,kk)
                          S3  = S3  + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(VELZ_VAR,  ii,jj,kk)
                          E   = E   + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(ENER_VAR,  ii,jj,kk)
                          Nel = Nel + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(YE_MSCALAR,ii,jj,kk)
                       end do
                       uCF(iNodeX,ic,jc,kc,iCF_D ) &
                       = D   * conv_dens
                       uCF(iNodeX,ic,jc,kc,iCF_S1) &
                       = S1  * conv_mom * uGF(iNodeX,ic,jc,kc,iGF_h_1)
                       uCF(iNodeX,ic,jc,kc,iCF_S2) &
                       = S2  * conv_mom * uGF(iNodeX,ic,jc,kc,iGF_h_2)
                       uCF(iNodeX,ic,jc,kc,iCF_S3) &
                       = S3  * conv_mom * uGF(iNodeX,ic,jc,kc,iGF_h_3)
                       uCF(iNodeX,ic,jc,kc,iCF_E ) &
                       = E   * conv_enr
                       uCF(iNodeX,ic,jc,kc,iCF_Ne) &
                       = Nel * conv_ne
                    end do

                    ! Copy radiation vars from Flash arrays into Thornado arrays
                    do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS
                    do iE = 1-THORNADO_SWE, THORNADO_NE+THORNADO_SWE

                       ioff = THORNADO_BEGIN &
                          + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE) &
                                    *THORNADO_NMOMENTS) &
                          + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                          + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                       do iNode = 1, THORNADO_RAD_NDOF

                          iNodeE = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                          iNodeX = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                          ii     = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
                          jj     = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                          kk     = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

                          ivar = ioff + iNodeE - 1

                          uCR(iNode,iE,ic,jc,kc,iCR,iS) = solnData(ivar,ii,jj,kk)

                       end do
                    end do ; end do ; end do

                 end do
              end do
           end do

           call Timers_stop("rt_reconstruction")

#if defined(THORNADO_OMP_OL)
           !$OMP TARGET ENTER DATA &
           !$OMP MAP( to: uCF, uCR, uGF )
#elif defined(THORNADO_OACC)
           !$ACC ENTER DATA &
           !$ACC COPYIN( uCF, uCR, uGF )
#endif

           ! Call the Fortran interface that lives in the Thornado repo
           call Timers_start("rt_imex")
           call Update_IMEX_PDARS(dt_sub*Second, uCF, uCR, &
              Explicit_Option = rt_doExplicit, Implicit_Option = rt_doImplicit, &
              SingleStage_Option = .false., BoundaryCondition_Option = ApplyBC )
           call Timers_stop("rt_imex")

#if defined(THORNADO_OMP_OL)
           !$OMP TARGET EXIT DATA &
           !$OMP MAP( from: uCF, uCR ) &
           !$OMP MAP( release: uGF )
#elif defined(THORNADO_OACC)
           !$ACC EXIT DATA &
           !$ACC COPYOUT( uCF, uCR ) &
           !$ACC DELETE( uGF )
#endif

           ! Copy back from the Thornado arrays into Flash arrays
           call Timers_start("rt_projection")
           do kc = 1, nX(3)
              do jc = 1, nX(2)
                 do ic = 1, nX(1)
                    i = lo(IAXIS) + THORNADO_NNODESX*(ic-1)
                    j = lo(JAXIS) + THORNADO_NNODESX*(jc-1)
                    k = lo(KAXIS) + THORNADO_NNODESX*(kc-1)

                    ! Interpolate fluid vars from Thornado arrays into Flash arrays
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       ii = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
                       jj = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                       kk = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k
                       D   = 0.0
                       S1  = 0.0
                       S2  = 0.0
                       S3  = 0.0
                       E   = 0.0
                       Nel = 0.0
                       do n = 1, THORNADO_FLUID_NDOF
                          D   = D   + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_D )
                          S1  = S1  + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_S1) &
                                / uGF(n,ic,jc,kc,iGF_h_1)
                          S2  = S2  + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_S2) &
                                / uGF(n,ic,jc,kc,iGF_h_2)
                          S3  = S3  + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_S3) &
                                / uGF(n,ic,jc,kc,iGF_h_3)
                          E   = E   + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_E )
                          Nel = Nel + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_Ne)
                       end do
                       solnData(DENS_VAR,  ii,jj,kk) = D   / conv_dens
                       solnData(VELX_VAR,  ii,jj,kk) = S1  / conv_mom  / solnData(DENS_VAR,ii,jj,kk)
                       solnData(VELY_VAR,  ii,jj,kk) = S2  / conv_mom  / solnData(DENS_VAR,ii,jj,kk)
                       solnData(VELZ_VAR,  ii,jj,kk) = S3  / conv_mom  / solnData(DENS_VAR,ii,jj,kk)
                       solnData(ENER_VAR,  ii,jj,kk) = E   / conv_enr  / solnData(DENS_VAR,ii,jj,kk)
                       solnData(YE_MSCALAR,ii,jj,kk) = Nel / conv_ne   / solnData(DENS_VAR,ii,jj,kk)
#ifdef EINT_VAR
                       ekin = 0.5*SUM(solnData(VELX_VAR:VELZ_VAR,ii,jj,kk)**2)
                       solnData(EINT_VAR,ii,jj,kk) = solnData(ENER_VAR,ii,jj,kk) - ekin
#endif
                    end do

                    ! Copy radiation vars from Thornado arrays into Flash arrays
                    do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS
                    do iE = 1-THORNADO_SWE, THORNADO_NE+THORNADO_SWE

                       ioff = THORNADO_BEGIN &
                          + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE) &
                                    *THORNADO_NMOMENTS) &
                          + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                          + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                       do iNode = 1, THORNADO_RAD_NDOF

                          iNodeE = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                          iNodeX = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                          ii     = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
                          jj     = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                          kk     = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

                          ivar = ioff + iNodeE - 1

                          solnData(ivar,ii,jj,kk) = uCR(iNode,iE,ic,jc,kc,iCR,iS)

                       end do
                    end do ; end do ; end do

                 end do
              end do
           end do
           call Timers_stop("rt_projection")

           call Timers_start("rt_finalize")
           call FreeThornado_Patch()
           call Timers_stop("rt_finalize")

           call tileDesc%releaseDataPtr(solnData, CENTER)

           call itor%next()
        end do
        call Grid_releaseTileIterator(itor)

        t_sub = t_sub + dt_sub


        if ( t_sub < dt ) then
           call Timers_start("rt_gc_sub")
           call Grid_fillGuardCells(CENTER,ALLDIR, &
              minLayers=my_ngrow*THORNADO_NNODESX,doEos=.false.,eosMode=rt_eosModeGc, &
              maskSize=NUNK_VARS,mask=rt_gcMask,makeMaskConsistent=.false.)
           call Timers_stop("rt_gc_sub")
        end if
     end do
  end do

  call Timers_start("eos")

  call Eos_everywhere(MODE_DENS_EI)

  call Timers_stop("eos")

  call Timers_stop("RadTrans")

  return

end subroutine RadTrans
