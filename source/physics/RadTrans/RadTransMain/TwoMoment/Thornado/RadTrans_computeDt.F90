!!****f* source/physics/RadTrans/RadTrans_computeDt
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
!!   RadTrans
!!
!!  SYNOPSIS
!!
!!   RadTrans_computeDt(Grid_tile_t(IN) :: tileDesc,
!!                      real,pointer    :: solnData(:,:,:,:),   
!!                      real,(INOUT)    :: dt_radtrans, 
!!                      integer(INOUT)  :: dtMinLoc(:))
!!
!!  DESCRIPTION 
!!    Compute radiative transfer time step
!!
!!  ARGUMENTS
!!    tileDesc      --  meta-information about the tile/block
!!    solnData      --  the physical, solution data from grid
!!    dt_radtrans   --  variable to hold timestep constraint
!!    dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine RadTrans_computeDt(tileDesc, solnData, dt_radtrans, dtMinLoc)
  
  use Grid_tile                   , only : Grid_tile_t
  use RadTrans_data               , only : rt_geometry, rt_str_geometry, &
                                           rt_meshMe, rt_cfl
  use rt_data                     , only : rt_doExplicit
  use UnitsModule                 , only : Centimeter
  use PhysicalConstantsModule     , only : SpeedOfLightCGS
  use ProgramHeaderModule         , only : nDimsX, nNodes, iX_B0, iX_E0
  use MeshModule                  , only : MeshX
  use GeometryFieldsModule        , only : uGF, iGF_h_1, iGF_h_2, iGF_h_3
  use ThornadoInitializationModule, only : InitThornado_Patch, &
                                           FreeThornado_Patch

  implicit none

  type(Grid_tile_t), intent(IN) :: tileDesc
  real,                 pointer :: solnData(:,:,:,:) 
  real,           intent(INOUT) :: dt_radtrans
  integer,        intent(INOUT) :: dtMinLoc(5)

  integer :: nX(3), swX(3), lo(MDIM), hi(MDIM)
  integer :: iX1, iX2, iX3
  real    :: boundBox(LOW:HIGH,MDIM)
  real    :: xL(3), xR(3)
  real    :: CFL, TimeStep, dt(3)

  if ( .not. rt_doExplicit ) then
     dt_radtrans = HUGE( 1.0 )
     dtMinLoc(1) = tileDesc % limits(LOW,IAXIS)
     dtMinLoc(2) = tileDesc % limits(LOW,JAXIS)
     dtMinLoc(3) = tileDesc % limits(LOW,KAXIS)
     dtMinLoc(4) = tileDesc % level
     dtMinLoc(5) = rt_meshMe
     return
  end if

  CFL = rt_cfl / ( SpeedOfLightCGS * DBLE( NDIM * ( 2 * THORNADO_NNODES - 1 ) ) )

  TimeStep = HUGE( 1.0 )
  dt       = HUGE( 1.0 )

  nX  = 1
  xL  = 0.0
  if     ( rt_geometry == CARTESIAN   )then
    xR = 1.0
  else if( rt_geometry == CYLINDRICAL )then
    xR = [ 1.0, 1.0, 2.0*PI ]
  else if( rt_geometry == SPHERICAL   )then
    xR = [ 1.0, PI , 2.0*PI ]
  end if

  lo = tileDesc % limits(LOW ,1:MDIM)
  hi = tileDesc % limits(HIGH,1:MDIM)

  nX(1:NDIM) = ( hi(1:NDIM) - lo(1:NDIM) + 1 ) / THORNADO_NNODESX

  swX = 0 ! --- Don't need ghost cells for this.

  call tileDesc % boundBox( boundBox )

  xL(1:NDIM) = boundBox(LOW ,1:NDIM)
  xR(1:NDIM) = boundBox(HIGH,1:NDIM)

  ! --- Convert From Flash-X's [ cm ] to thornado units

  xL(1) = xL(1) * Centimeter
  xR(1) = xR(1) * Centimeter
  if( rt_geometry == CARTESIAN .or. rt_geometry == CYLINDRICAL )then
    xL(2) = xL(2) * Centimeter
    xR(2) = xR(2) * Centimeter
    if( rt_geometry == CARTESIAN )then
      xL(3) = xL(3) * Centimeter
      xR(3) = xR(3) * Centimeter
    end if
  end if

  call InitThornado_Patch &
         ( nX, swX, xL, xR, THORNADO_NSPECIES, rt_str_geometry )

  ASSOCIATE &
    ( dX1 => MeshX(1) % Width, &
      dX2 => MeshX(2) % Width, &
      dX3 => MeshX(3) % Width )

  DO iX3 = iX_B0(3), iX_E0(3)
  DO iX2 = iX_B0(2), iX_E0(2)
  DO iX1 = iX_B0(1), iX_E0(1)

    dt(1)   = dX1(iX1) * MINVAL( uGF(:,iX1,iX2,iX3,iGF_h_1) )

    IF( nDimsX > 1 )THEN

      dt(2) = dX2(iX2) * MINVAL( uGF(:,iX1,iX2,iX3,iGF_h_2) )

    END IF

    IF( nDimsX > 2 )THEN

      dt(3) = dX3(iX3) * MINVAL( uGF(:,iX1,iX2,iX3,iGF_h_3) )

    END IF

    TimeStep = MIN( TimeStep, CFL * MINVAL( dt ) / Centimeter )

    IF( TimeStep < dt_radtrans )THEN

      dt_radtrans = TimeStep

      dtMinLoc(1) = THORNADO_NNODESX * (iX1-1) + lo(IAXIS)
      dtMinLoc(2) = THORNADO_NNODESX * (iX2-1) + lo(JAXIS)
      dtMinLoc(3) = THORNADO_NNODESX * (iX3-1) + lo(KAXIS)
      dtMinLoc(4) = tileDesc % level
      dtMinLoc(5) = rt_meshMe

    END IF

  END DO
  END DO
  END DO

  END ASSOCIATE ! dX1, etc.

  call FreeThornado_Patch()

  return
end subroutine RadTrans_computeDt
