!!****f* source/Simulation/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockId)
!!
!!
!!
!! DESCRIPTION
!!  This routine applies initial conditions of a specific simulation
!!  to the specified block.
!!
!! 
!! ARGUMENTS
!!
!!  blockId -         the number of the block to update
!!
!! PARAMETERS
!!
!!  eosModeInit -     after this routine sets up initial conditions,
!!                    the grid package calls Eos to insure the
!!                    values are thermodynamically consistent.  This
!!                    parameter controls the mode of application of
!!                    the Eos.  Its default is "dens_ie", and it can
!!                    be one of "dens_ie", "dens_pres", "dens_temp".
!!                    Setting this value to dens_temp, for instance,
!!                    would make it possible to leave this routine
!!                    having just initialized density and temperature,
!!                    leaving Eos to calculate the rest.
!!
!! SEE ALSO
!!
!!  Eos_wrapped
!!***

subroutine Simulation_initBlock(solnData, tileDesc)
  ! the inclusion of pAmbient and rhoAmbient for hydro/eos usage
  use Simulation_data, ONLY: sim_pAmbient, sim_rhoAmbient, sim_tempAmbient, sim_YeAmbient
  use Particles_data, ONLY : R, rt_gamma
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_putPointData,&
       Grid_getCellCoords, Grid_getMinCellSize
  use Grid_tile, ONLY : Grid_tile_t
  USE Eos_data, ONLY : eos_singleSpeciesA

  implicit none
#include "constants.h"
#include "Flash.h"

  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc
  integer  ::  i, j, k
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  ! for hydro/eos
  REAL :: vx, vy, vz, p, ek, e, rho
  REAL :: temp
  ! for cross-cell advancement checking
  REAL, DIMENSION(LOW:HIGH, MDIM) :: bndBox
  REAL, DIMENSION(MDIM) :: deltaCell
  INTEGER, PARAMETER :: ycoordsize = NYB + 2*NGUARD
  REAL, DIMENSION(ycoordsize) :: ycoords
  INTEGER, PARAMETER :: xcoordsize = NXB + 2*NGUARD
  REAL, DIMENSION(xcoordsize) :: xcoords
  INTEGER, PARAMETER :: zcoordsize = NZB + 2*NGUARD
  REAL, DIMENSION(zcoordsize) :: zcoords
  integer, dimension(1:MDIM) :: lo, hi, u_lo, u_hi

  ! get dimensions/limits and coordinates
  blkLimits = tileDesc%limits
  blkLimitsGC = tileDesc%blkLimitsGC
  lo(1:MDIM) = blkLimits(LOW,1:MDIM)
  hi(1:MDIM) = blkLimits(HIGH,1:MDIM)


  CALL Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xcoords)
  if (NDIM > 1) then
    CALL Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, ycoords)

    if (NDIM > 2) then
      CALL Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zcoords)
    end if
  end if

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
      do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)

        vx = 0.0d0
        vy = 0.0d0
        vz = 0.0d0
        ! gas kinetic energy
        ek = 0.5 * (vx*vx + vy*vy + vz*vz)

        rho = sim_rhoAmbient  ! fixed nH to be 100.
        temp = sim_tempAmbient
        p = (R/eos_singleSpeciesA)*rho*temp ! just dummy value, init with dens_temp
        e = p / (rt_gamma - 1.0)
        e = e / rho
        e = e + ek

        solnData(DENS_VAR,i,j,k) = rho
        solnData(TEMP_VAR,i,j,k) = temp
        solnData(PRES_VAR,i,j,k) = p
        solnData(EINT_VAR,i,j,k) = e - ek
        solnData(ENER_VAR,i,j,k) = e
        solnData(GAME_VAR,i,j,k) = rt_gamma
        solnData(GAMC_VAR,i,j,k) = rt_gamma
        solnData(VELX_VAR,i,j,k) = vx
        solnData(VELY_VAR,i,j,k) = vy
        solnData(VELZ_VAR,i,j,k) = vz
#ifdef YE_MSCALAR
        solnData(YE_MSCALAR,i,j,k) = sim_YeAmbient
#endif
      enddo
    enddo
  enddo
 
  return
end subroutine Simulation_initBlock
