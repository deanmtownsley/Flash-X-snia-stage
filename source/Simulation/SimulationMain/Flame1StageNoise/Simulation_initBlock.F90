!!****if* source/Simulation/SimulationMain/Flame1StageNoise/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            Grid_tile_t(IN)  :: tileDesc  )
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   solnData  - pointer to solution data
!!   tileDesc -  describes the tile or block to initialize
!!
!! AUTOGENROBODOC
!!
!!
!!***

!!  Initialization of flame in pressure equilibrium
!!  The flame front can be planar or spherical in multiple dimensions
!!  depending on the value of the pseudo_1d parameter.
!!  See description of parameters in Config file for more info.
!!
! Dean Townsley 2008

subroutine Simulation_initBlock(solnData,tileDesc)
  
  use Simulation_data
  use Grid_interface, ONLY : Grid_getCellCoords, &
                             Grid_getCellVolumes, &
                             Grid_subcellGeometry, &
                             Grid_getDeltas
  use Grid_tile, ONLY : Grid_tile_t
  use Flame_interface, ONLY : Flame_getProfile, Flame_rhJump
  use fl_effData, ONLY: fl_effDeltae

  implicit none
#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

  
  real, pointer :: solnData(:,:,:,:)
  type(Grid_tile_t), intent(in) :: tileDesc

  integer :: i, j, k

  integer, dimension(MDIM) :: cell
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  integer,dimension(LOW:HIGH,MDIM) :: tileLimits
  integer,dimension(LOW:HIGH,MDIM) :: grownTileLimits



  real, dimension(EOS_NUM) :: eosData
  real :: flam, velx
  real :: fsurf_x_position, fsurf_distance
  real :: ye, yi

!==============================================================================

  ! get the coordinate information for the current block
  tileLimits = tileDesc%limits
  grownTileLimits = tileDesc%grownLimits
  call tileDesc%deltas(deltas)  ! alternatively: call Grid_getDeltas(tileDesc%level, deltas)

  allocate(iCoords(grownTileLimits(LOW, IAXIS):grownTileLimits(HIGH, IAXIS)))
  allocate(jCoords(grownTileLimits(LOW, JAXIS):grownTileLimits(HIGH, JAXIS)))
  allocate(kCoords(grownTileLimits(LOW, KAXIS):grownTileLimits(HIGH, KAXIS)))

  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          iCoords)
  call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          jCoords)
  call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          kCoords)

  !-----------------------------------------------
  ! loop over all zones and init
  !-----------------------------------------------
  do k = grownTileLimits(LOW,KAXIS), grownTileLimits(HIGH,KAXIS)
     do j = grownTileLimits(LOW, JAXIS), grownTileLimits(HIGH, JAXIS)
        do i = grownTileLimits(LOW,IAXIS), grownTileLimits(HIGH, IAXIS)

           if (.not. sim_ignite) then
              ! no burned material, only unburned
              eosData(:) = sim_eosData_u(:)
              flam = 0.0
              velx=0.0
           else
              !-----------------------------------------------
              ! initialize, including a burned region
              !-----------------------------------------------

              ! find distance from flame surface (positive is in front of flame)
              if (sim_pseudo1d) then
                 ! planar flame surface with normal tilted up by angle theta from x direction
                 fsurf_x_position = sim_fracPerturb*(sim_xmax-sim_xmin) &
                                    - tan(PI*sim_theta/180.0)*(jCoords(j)-sim_yctrPerturb)
                 fsurf_distance = iCoords(i) - fsurf_x_position
              else
                 ! n-dimensional sphere centered at specified point with radius
                 !   frac_perturb* (x domain extent)
                 fsurf_distance = sqrt( (iCoords(i)-sim_xctrPerturb)**2 + &
                                        (jCoords(j)-sim_yctrPerturb)**2 + &
                                        (kCoords(k)-sim_zctrPerturb)**2 ) &
                                  - sim_fracPerturb*(sim_xmax-sim_xmin)
              endif

              ! determine local state in this zone
              if ( fsurf_distance > 1.5*sim_laminarWidth ) then
                 ! unburned material
                 eosData(:) = sim_eosData_u(:)
                 flam = 0.0
              else if ( fsurf_distance < -1.5*sim_laminarWidth ) then
                 ! fully burned
                 eosData(:) = sim_eosData_b(:)
                 flam = 1.0
              else
                 ! partially burned
                 call Flame_getProfile(fsurf_distance, flam)

                 ! calculate propertise for partially burned material
                 ! note, in fact ye_f and ye_a should be equal
                 yi = 1.0/sim_eosData_u(EOS_ABAR)*(1.0-flam) + (1.0/sim_eosData_b(EOS_ABAR))*flam
                 ye = sim_eosData_u(EOS_ZBAR)/sim_eosData_u(EOS_ABAR)*(1.0-flam) + (sim_eosData_b(EOS_ZBAR)/sim_eosData_b(EOS_ABAR))*flam
                 eosData(:) = sim_eosData_u(:)
                 eosData(EOS_ABAR) = 1.0/yi
                 eosData(EOS_ZBAR) = ye*eosData(EOS_ABAR)

                 ! put this in pressure equilibrium with unburned material
                 call Flame_rhJump(sim_eosData_u, eosData, flam*fl_effDeltae, 0.0, MODE_DENS_TEMP)

              endif

              ! init velocity field, nonzero only makes sense with a planar flame front
              if (sim_pseudo1d) then
                 velx = sim_flamespeed*sim_eosData_u(EOS_DENS)* &
                            (1.e0/sim_eosData_b(EOS_DENS) - 1.e0/eosData(EOS_DENS))
              else
                 velx = 0.0
              endif

           endif ! sim_ignite
           

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           
           solnData(DENS_VAR,i,j,k)= eosData(EOS_DENS)
           solnData(TEMP_VAR,i,j,k)= eosData(EOS_TEMP)
           solnData(FLAM_MSCALAR,i,j,k)= flam

           solnData(VELX_VAR,i,j,k) = velx
           solnData(VELY_VAR,i,j,k) = 0.0
           solnData(VELZ_VAR,i,j,k) = 0.0

           !  usually I would just call the EOS, but we happen to have all this data
           !  so we'll just put it in.
           solnData(ENER_VAR,i,j,k) = eosData(EOS_EINT)+0.5*velx**2
           solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
           solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
           solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
           solnData(GAME_VAR,i,j,k) = eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT))+1.0

        enddo
     enddo
  enddo
  
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

  return
  
end subroutine Simulation_initBlock






