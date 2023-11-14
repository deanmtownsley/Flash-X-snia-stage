!!****f* source/RadTrans/RadTrans_molExplicitRHS
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
!!      RadTrans_molExplicitRHS
!!
!!  SYNOPSIS
!!
!!      call RadTrans_molExplicitRHS(real,    intent(in) :: t,
!                                    integer, intent(in) :: activeRHS
!!                                   real,    intent(in) :: dtWeight)
!!
!!  DESCRIPTION
!!
!!      Calculate explicit RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      t         : Current time
!!      activeRHS : RHS data struct to fill
!!      dtWeight  : Weight timestep (e.g. for flux corrections)
!!
!!***

!!REORDER(4): Uin, dU

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

subroutine RadTrans_molExplicitRHS(t, activeRHS, dtWeight)

  use Driver_interface, ONLY : Driver_abort
  use Grid_data, ONLY: gr_str_geometry
  use Grid_interface, ONLY : Grid_fillGuardCells, &
      Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  use MoL_interface, ONLY : MoL_getDataPtr, MoL_releaseDataPtr
  use RadTrans_data, ONLY : rt_useRadTrans, rt_enableTiling, &
     rt_eosModeGc, rt_gcMask, rt_str_geometry, rt_geometry
  use rt_data, ONLY : rt_doExplicit, rt_doImplicit, rt_irhs, rt_ivar
  use rt_tm_interface, ONLY : rt_tm_reconstruction, rt_tm_projection
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use FluidFieldsModule, ONLY : uCF
  use GeometryFieldsModuleE, ONLY : uGE
  use GeometryFieldsModule, ONLY : uGF
  use RadiationFieldsModule, ONLY : uCR, nCR, nSpecies
  use ProgramHeaderModule, ONLY : nDOFZ, iZ_B0, iZ_B1, iZ_E0, iZ_E1, &
     iX_B0, iX_B1, iX_E0, iX_E1
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch
  use TimeSteppingModule_Flash, ONLY : ApplyBoundaryConditions_Radiation, &
     ApplyBoundaryConditions_Euler_FLASH
  use UnitsModule, ONLY : Centimeter

#if   defined(THORNADO_ORDER_1)
  USE TwoMoment_DiscretizationModule_Streaming, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
#elif defined(THORNADO_ORDER_V)
  USE TwoMoment_DiscretizationModule_Streaming_OrderV, ONLY: &
    ComputeIncrement_TwoMoment_Explicit
#endif

  use Simulation_data, ONLY : sim_globalMe

   implicit none

   real, intent(in) :: t
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

  real, pointer, dimension(:,:,:,:) :: Uin, dU
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  integer, dimension(MDIM) :: lo,hi,u_lo,u_hi

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t) :: tileDesc

  integer, parameter :: my_ngrow = 1
  integer :: nX(3), swX(3), iZ_SW_P(4)
  real :: xL(3), xR(3)
  integer :: iX1, iX2, iX3, iS, iCR, iE, iNodeZ, iNodeE, iNodeX
  integer :: i, j, k, ii, jj, kk, irhs, ivar

  real, allocatable, dimension(:,:,:,:,:,:,:) :: d_uCR

  integer :: faces(LOW:HIGH, 1:MDIM)
  integer :: onBoundary(LOW:HIGH, 1:MDIM)
  integer :: ApplyBC

  real, parameter :: conv_x = Centimeter

  if ( .not. rt_doExplicit ) return

  if ( my_ngrow /= 1 ) then
     call Driver_abort("Need one ghost cell in call_to_thornado!")
  else if ( my_ngrow*THORNADO_NNODESX > NGUARD ) then
     call Driver_abort("NGUARD must be at least my_ngrow*THORNADO_NNODESX")
  end if

  call Timers_start("RadTrans_molExplicitRHS")

  call Timers_start("rt_gc")
  call Grid_fillGuardCells(CENTER,ALLDIR, &
    minLayers=my_ngrow*THORNADO_NNODESX,doEos=.true.,eosMode=rt_eosModeGc, &
    maskSize=NUNK_VARS,mask=rt_gcMask,makeMaskConsistent=.true.)
  call Timers_stop("rt_gc")

  call Grid_getTileIterator(itor, LEAF, tiling=rt_enableTiling)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     ! Get a pointer to solution data
     call MoL_getDataPtr(tileDesc, Uin, MOL_EVOLVED)
     call MoL_getDataPtr(tileDesc, dU, activeRHS)

     ! Check if the block needs boundary condition
     ! RC: This only works for the reflecting inner boundary in 1D
     call tileDesc%faceBCs(faces)
     if( faces(LOW,IAXIS) == REFLECTING .and. rt_geometry == SPHERICAL )then
       ApplyBC = 3
     else
       ApplyBC = 0
     end if

     ! get dimensions/limits and coordinates
     nX = 1
     swX = 0
     iZ_SW_P = 0
     xL = 0.0
     if ( rt_geometry == CARTESIAN ) then
        xR = 1.0
     else if ( rt_geometry == CYLINDRICAL ) then
        xR = [ 1.0, 1.0, 2.0*PI ]
     else if ( rt_geometry == SPHERICAL ) then
        xR = [ 1.0, PI, 2.0*PI ]
     end if

     lo = tileDesc%limits(LOW,1:MDIM)
     hi = tileDesc%limits(HIGH,1:MDIM)
     nX(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
     swX(1:NDIM) = my_ngrow
     iZ_SW_P(2:NDIM+1) = swX(1:NDIM)
     u_lo = 1 - swX
     u_hi = nX + swX

     call tileDesc%boundBox(boundBox)
     xL(1:NDIM) = boundBox(LOW, 1:NDIM)
     xR(1:NDIM) = boundBox(HIGH,1:NDIM)

     ! Thornado uses units where c = G = k = 1, Meter = 1
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

     !@M dir_enter_data &
     !@M dir_copyin()( @M dir_map_type(to) Uin, lo, hi, u_lo, u_hi, nX ) &
     !@M dir_end

     call InitThornado_Patch(nX, swX, xL, xR, THORNADO_NSPECIES, rt_str_geometry)
     allocate( d_uCR(1:nDOFZ, &
                     iZ_B1(1):iZ_E1(1), &
                     iZ_B1(2):iZ_E1(2), &
                     iZ_B1(3):iZ_E1(3), &
                     iZ_B1(4):iZ_E1(4), &
                     1:nCR,1:nSpecies) )
     call Timers_stop("rt_init")

     call Timers_start("rt_reconstruction")
     call rt_tm_reconstruction(Uin, nX, lo, hi, u_lo, u_hi)
     call Timers_stop("rt_reconstruction")

      CALL ApplyBoundaryConditions_Radiation &
             ( iZ_SW_P, iZ_B0, iZ_E0, iZ_B1, iZ_E1, uCR, ApplyBC )

#if   defined(THORNADO_ORDER_V)
      CALL ApplyBoundaryConditions_Euler_FLASH &
             ( iZ_SW_P, iX_B0, iX_E0, iX_B1, iX_E1, uCF, ApplyBC )
#endif

     ! Call the Fortran interface that lives in the Thornado repo
     call Timers_start("rt_explicit")
#if   defined(THORNADO_ORDER_1)
      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
               uGE, uGF, uCR, d_uCR )
#elif defined(THORNADO_ORDER_V)
      CALL ComputeIncrement_TwoMoment_Explicit &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, &
               uGE, uGF, uCF, uCR, d_uCR )
#endif
     call Timers_stop("rt_explicit")

     ! Copy back from the Thornado arrays into Flash arrays
     call Timers_start("rt_projection")
     !call rt_tm_projection(Uin, nX, lo, hi, u_lo, u_hi)

     !@M dir_loop(7) &
     !@M dir_private()( i, j, k, ii, jj, kk, irhs, iNodeE, iNodeX ) &
     !@M dir_present()( Uin, uCR, nX, lo ) &
     !@M dir_end

     do iX3 = 1, nX(3)
        do iX2 = 1, nX(2)
           do iX1 = 1, nX(1)
              do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS
                 do iE = 1-THORNADO_SWE, THORNADO_NE+THORNADO_SWE
                    do iNodeZ = 1, THORNADO_RAD_NDOF

                       i = lo(IAXIS) + THORNADO_NNODESX*(iX1-1)
                       j = lo(JAXIS) + THORNADO_NNODESX*(iX2-1)
                       k = lo(KAXIS) + THORNADO_NNODESX*(iX3-1)

                       iNodeE = mod((iNodeZ-1)                 ,THORNADO_NNODESE   ) + 1
                       iNodeX = mod((iNodeZ-1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                       ii     = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
                       jj     = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                       kk     = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

                       irhs = rt_irhs(iNodeE,iE,iCR,iS)

                       dU(irhs,ii,jj,kk) = dU(irhs,ii,jj,kk) + d_uCR(iNodeZ,iE,iX1,iX2,iX3,iCR,iS)

                    end do
                 end do
              end do ; end do
           end do
        end do
     end do
     call Timers_stop("rt_projection")

     if ( sim_globalMe == MASTER_PE ) then
       ivar = rt_ivar(1,1,1,1)
       write(*,*) Uin(ivar,lo(IAXIS)  ,lo(JAXIS)  ,lo(KAXIS)  )
       write(*,*) Uin(ivar,lo(IAXIS)-1,lo(JAXIS)  ,lo(KAXIS)  ), &
                  Uin(ivar,lo(IAXIS)-2,lo(JAXIS)  ,lo(KAXIS)  )
       write(*,*) Uin(ivar,lo(IAXIS)  ,lo(JAXIS)-1,lo(KAXIS)  ), &
                  Uin(ivar,lo(IAXIS)  ,lo(JAXIS)-2,lo(KAXIS)  )
       write(*,*) Uin(ivar,lo(IAXIS)  ,lo(JAXIS)  ,lo(KAXIS)-1), &
                  Uin(ivar,lo(IAXIS)  ,lo(JAXIS)  ,lo(KAXIS)-2)
       write(*,*) d_uCR(1,1,1,1,1,1,1)
     end if

     call Timers_start("rt_finalize")

     !@M dir_exit_data &
     !@M dir_copyout()( @M dir_map_type(from) Uin ) &
     !@M dir_delete()( @M dir_map_type(release) lo, hi, u_lo, u_hi, nX ) &
     !@M dir_end

     deallocate( d_uCR )
     call FreeThornado_Patch()
     call Timers_stop("rt_finalize")

     call MoL_releaseDataPtr(tileDesc, dU, activeRHS)
     call MoL_releaseDataPtr(tileDesc, Uin, MOL_EVOLVED)

     call itor%next()
  end do ! iterator loop
  call Grid_releaseTileIterator(itor)

  call Timers_stop("RadTrans_molExplicitRHS")

   return
end subroutine RadTrans_molExplicitRHS
