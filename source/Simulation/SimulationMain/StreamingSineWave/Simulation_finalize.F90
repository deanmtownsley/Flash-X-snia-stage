!!****f* source/Simulation/Simulation_finalize
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
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This dummy function cleans up the Simulation unit, deallocates memory, etc.
!!  However, as nothing needs to be done, only this stub is included.
!!
!! ARGUMENTS
!!
!!
!!
!!***

!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"

subroutine Simulation_finalize()

  use Driver_interface, ONLY : Driver_getSimTime
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  use Simulation_data, ONLY : sim_globalMe, sim_globalComm

  use rt_tm_interface, ONLY : rt_tm_reconstruction, rt_tm_projection

  use FluidFieldsModule, ONLY : uCF, nCF
  USE GeometryFieldsModule, ONLY : uGF
  use ProgramHeaderModule, ONLY : iZ_B0, iZ_E0, iZ_B1, iZ_E1
  use RadiationFieldsModule, ONLY : uCR, uPR, uAR, uGR
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch
  use TwoMoment_UtilitiesModule, ONLY : ComputeFromConserved_TwoMoment
  use UnitsModule, ONLY : Centimeter, Second

#include "Flashx_mpi_implicitNone.fh"
  
  real, dimension(:,:,:,:), pointer :: solnData

  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real :: simTime

  real :: MaxError_local(THORNADO_NMOMENTS,THORNADO_NSPECIES)
  real :: MaxError(THORNADO_NMOMENTS,THORNADO_NSPECIES)

  integer, dimension(1:MDIM) :: lo, hi, u_lo, u_hi, loGC, hiGC

  integer, parameter :: my_ngrow = 0
  integer :: nX(3), swX(3), iS, iPR, ierr
  real :: xL(3), xR(3)

  real, parameter :: conv_x = Centimeter

  nullify(solnData)

  call Driver_getSimTime(simTime)

  MaxError_local = 0.0
  MaxError = 0.0

  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     ! Get a pointer to solution data
     call tileDesc%getDataPtr(solnData, CENTER)

     ! get dimensions/limits and coordinates
     nX = 1
     swX = 0
     xL = 0.0
     xR = 1.0

     lo = tileDesc%limits(LOW,1:MDIM)
     hi = tileDesc%limits(HIGH,1:MDIM)
     loGC = tileDesc % blkLimitsGC(LOW ,1:MDIM)
     hiGC = tileDesc % blkLimitsGC(HIGH,1:MDIM)
     nX(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
     swX(1:NDIM) = my_ngrow
     u_lo = 1 - swX
     u_hi = nX + swX

     call tileDesc%boundBox(boundBox)
     xL(1:NDIM) = boundBox(LOW, 1:NDIM)
     xR(1:NDIM) = boundBox(HIGH,1:NDIM)

     ! convert cm to m for Thornado (cartesian geometry assumed)
     xL = xL * conv_x
     xR = xR * conv_x

     call InitThornado_Patch(nX, swX, xL, xR, THORNADO_NSPECIES, 'cartesian')

     call rt_tm_reconstruction(solnData, nX, lo, hi, loGC, hiGC, u_lo, u_hi, tileDesc % level )

#if   defined(THORNADO_OACC)
     !$ACC UPDATE HOST( uGF, uCF, uCR )
#elif defined(THORNADO_OMP_OL)
     !$OMP TARGET UPDATE FROM( uGF, uCF, uCR )
#endif

     ! get primitive radiation fields
     call ComputeFromConserved_TwoMoment &
        ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR, uAR, uGR )

     call ComputeError_SineWaveStreaming( simTime * Second, MaxError_local )

     call FreeThornado_Patch()

     call tileDesc%releaseDataPtr(solnData, CENTER)

     call itor%next()
  end do ! iterator loop
  call Grid_releaseTileIterator(itor)

  call MPI_AllReduce(MaxError_local(1,1),MaxError(1,1), &
    THORNADO_NMOMENTS*THORNADO_NSPECIES,FLASH_REAL,MPI_MAX,sim_globalComm,ierr)

  if ( sim_globalMe == MASTER_PE ) then
     WRITE(*,*)
     WRITE(*,'(A2,A)') '', 'INFO: SineWaveStreaming Error'
     WRITE(*,*)
     WRITE(*,'(A4,A2,4A12)') '', 'Sp', 'N', 'G1', 'G2', 'G3'
     DO iS = 1, THORNADO_NSPECIES
        WRITE(*,'(A4,I2.2,4ES12.4E2)') '', iS, MaxError(:,iS)
     END DO
  end if

  return

contains

  SUBROUTINE ComputeError_SineWaveStreaming( t, MaxError )

    USE KindModule, ONLY: &
       DP, Zero, TwoPi
    USE MeshModule, ONLY: &
       MeshX, NodeCoordinate
    USE ProgramHeaderModule, ONLY: &
       nDOFZ, iE_B0, iE_E0, iX_B0, iX_E0
    USE RadiationFieldsModule, ONLY: &
       uPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
    USE ReferenceElementModule, ONLY: &
       NodeNumberTable

    REAL(DP), INTENT(in) :: t
    REAL(DP), INTENT(inout) :: MaxError(THORNADO_NMOMENTS,THORNADO_NSPECIES)

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeZ, iNodeX1
    REAL(DP) :: X1, D_A, I1_A, I2_A, I3_A

    DO iS  = 1       , THORNADO_NSPECIES
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0   , iE_E0

      DO iNodeZ = 1, nDOFZ

        iNodeX1 = NodeNumberTable(2,iNodeZ)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        D_A  = 0.50_DP + 0.49_DP * SIN( TwoPi * ( X1 - t ) )
        I1_A = D_A
        I2_A = Zero
        I3_A = Zero

        MaxError(iPR_D, iS) &
          = MAX( ABS( D_A  - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_D ,iS) ), &
                 MaxError(iPR_D ,iS) )

        MaxError(iPR_I1,iS) &
          = MAX( ABS( I1_A - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I1,iS) ), &
                 MaxError(iPR_I1,iS) )

        MaxError(iPR_I2,iS) &
          = MAX( ABS( I2_A - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I2,iS) ), &
                 MaxError(iPR_I2,iS) )

        MaxError(iPR_I3,iS) &
          = MAX( ABS( I3_A - uPR(iNodeZ,iE,iX1,iX2,iX3,iPR_I3,iS) ), &
                 MaxError(iPR_I3,iS) )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE ComputeError_SineWaveStreaming

end subroutine Simulation_finalize
