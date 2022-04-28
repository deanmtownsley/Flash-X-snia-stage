!!****if* source/Simulation/SimulationMain/StreamingDopplerShift/user_bnd
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
!!    user_bnd
!!
!!  SYNOPSIS
!!
!! DESCRIPTION
!!              Called for a block if the boundary condition is bnd_user.
!!              The boundary condition for this StreamingDopplerShift is
!!              fixed radiation field boundary.
!!              This routine is called by Grid_bcApplyToRegion routine.
!!
!!***

subroutine user_bnd( guard, axis, regionData, regionSize )

  use Simulation_data
  use Driver_interface, ONLY : Driver_abort
  use KindModule, ONLY : DP, Zero, One, Three
  use MeshModule, ONLY : NodeCoordinate, MeshE
  use RadiationFieldsModule, ONLY : iCR_N, iCR_G1, iCR_G2, iCR_G3
  use UnitsModule, only : MeV

  implicit none

#include "constants.h"
#include "Simulation.h"

  integer, intent(IN) :: guard, axis
  integer,dimension(REGION_DIM),intent(IN) :: regionSize
  real,dimension(regionSize(BC_DIR),&
       regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),&
       regionSize(STRUCTSIZE)),intent(INOUT)::regionData

  integer  :: i, ivar, je, ke
  integer  :: iNodeE, iE, iCR, iS
  real     :: E, D

  je = regionSize(SECOND_DIR)
  ke = regionSize(THIRD_DIR)

  do ivar = 1,regionSize(STRUCTSIZE)

    ! radiation field variables

    if( ivar >= THORNADO_BEGIN .and. ivar <= THORNADO_END )then

      iNodeE = mod(ivar - THORNADO_BEGIN, THORNADO_NNODESE) + 1
      iE     = mod((ivar-THORNADO_BEGIN) / &
               (THORNADO_NNODESE), THORNADO_NE+2*THORNADO_SWE) &
               + 1 - THORNADO_SWE
      iCR    = (ivar-THORNADO_BEGIN) / &
               (THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) + 1
      iS     = (ivar-THORNADO_BEGIN) / &
               (THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) + 1

      regionData(1:guard,1:je,1:ke,ivar) = Zero

      if( iE >= 1 .and. iE <= THORNADO_NE )then

        E = NodeCoordinate( MeshE,    iE,  iNodeE  )
        E = E / MeV

        select case( TRIM( sim_rad_spectrum ) )
          case( 'Fermi-Dirac' )
            D = One / ( EXP( E / Three - Three ) + One )
          case( 'Bose-Einstein' )
            D = One / ( EXP( E ) - One )
          case default
            D = One / ( EXP( E / Three - Three ) + One )
        end select

        ! Default = Zero
        regionData(1:guard,1:je,1:ke,ivar) = Zero

        ! Number Densities on Comoving Frame
        if (iCR == iCR_N) then
           regionData(1:guard,1:je,1:ke,ivar) = D / sim_s
        end if
        ! Numer Flux Densities on Comoving Frame
        if (iCR == axis + iCR_N ) then ! iCR_G1 == IXIS
           regionData(1:guard,1:je,1:ke,ivar) = 0.999_DP * D / sim_s
        end if

      end if

    else

     ! copy fluid field like outflowing

     do i = 1,guard
        regionData(i,1:je,1:ke,ivar) = regionData(guard+1,1:je,1:ke,ivar)
     end do

    end if

  end do

end subroutine user_bnd
