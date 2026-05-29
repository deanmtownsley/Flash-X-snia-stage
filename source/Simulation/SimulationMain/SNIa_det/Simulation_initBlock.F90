!!
!! Dean M. Townsley 2012, .. , 2026
!!
!! Simulation grid initialization routine for SNIa_det (started from WDSurfHeDetFullStar)
!! See source/Simulation/Simulation_initBlock.F90 for API and notes.
!!
!!

subroutine Simulation_initBlock(solnData, tileDesc)
  
  use Simulation_data, ONLY: sim_fluffBoundaryRadius, sim_densFluffInner, sim_tempFluffInner, &
     sim_densFluffOuter, sim_tempFluffOuter, &
     sim_refFluffThresh, sim_refFluffMargin, sim_refFluffLevel, &
     sim_refNogenEnucThresh, sim_refNogenMargin, sim_refNogenLevel, &
     sim_refBurnedProductThresh, &
     sim_refEjectaPhaseStartTime, sim_refEjectaPhaseMaxRes, &
     sim_ign_numpnts, &
     sim_ignite, sim_ign_keep_pres, sim_ign_hemispherical, sim_ign_file, &
     sim_ign_x, sim_ign_y, sim_ign_z, sim_ign_r, sim_ign_temp_center, sim_ign_temp_edge, &
     sim_wd_npnts, sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_he4_tab, sim_wd_c12_tab, &
     sim_wd_o16_tab, sim_wd_ne20_tab, &
     sim_wd_dr_inv
  use sim_local_interface, ONLY : sim_interpolate1dWd, &
                         sim_find_ign_state_keeping_pres
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, &
    Grid_getCellCoords
  use Eos_interface, ONLY : Eos
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
#include "Multispecies.h"

  
  real,              pointer    :: solnData(:,:,:,:)
  type(Grid_tile_t), intent(in) :: tileDesc

  integer :: i, j, k

  integer,dimension(LOW:HIGH,MDIM) :: tileLimits
  integer,dimension(LOW:HIGH,MDIM) :: grownTileLimits
  integer, dimension(MDIM) :: cell
  integer :: isizeGC, jsizeGC, ksizeGC
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  real, dimension(MDIM) :: deltas

  real, dimension(EOS_NUM) :: zone_state
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: zone_abund
  logical :: zone_fluff
  real :: radius, ign_dist, ign_center_radius, test_dist
  real :: enuc, ign_temp
  integer :: ipi, testi


!==============================================================================

  ! get essential info about this block - index limits and cell coordinates
  tileLimits = tileDesc%limits
  grownTileLimits = tileDesc%grownLimits

  allocate( iCoords( grownTileLimits(LOW,IAXIS), grownTileLimits(HIGH,IAXIS) ) )
  allocate( jCoords( grownTileLimits(LOW,JAXIS), grownTileLimits(HIGH,JAXIS) ) )
  allocate( kCoords( grownTileLimits(LOW,KAXIS), grownTileLimits(HIGH,KAXIS) ) )
  call Grid_getCellCoords( IAXIS, CENTER, tileDesc%level, &
                           grownTileLimits(LOW,:), grownTileLimits(HIGH,:), iCoords)
  call Grid_getCellCoords( JAXIS, CENTER, tileDesc%level, &
                           grownTileLimits(LOW,:), grownTileLimits(HIGH,:), jCoords)
  call Grid_getCellCoords( KAXIS, CENTER, tileDesc%level, &
                           grownTileLimits(LOW,:), grownTileLimits(HIGH,:), kCoords)

  call tileDesc%deltas(deltas)

  !-----------------------------------------------
  ! loop over all zones and init
  !-----------------------------------------------
  do k = tileLimits(LOW,KAXIS), tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS), tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS), tileLimits(HIGH,IAXIS)

           !-----------------------------------------------
           !  determine state of material at this radius if unburned
           !  from external 1-d hyrdostatic model
           !-----------------------------------------------

           radius = iCoords(i)**2
           if (NDIM >= 2) radius = radius + jCoords(j)**2
           if (NDIM == 3) radius = radius + kCoords(k)**2
           radius = sqrt(radius)

           zone_abund(:) = 0.0
           call sim_interpolate1dWd(radius, deltas(IAXIS), zone_state(EOS_DENS), zone_state(EOS_TEMP), &
                                    zone_abund(HE4_SPEC), zone_abund(C12_SPEC), zone_abund(O16_SPEC), &
                                    zone_abund(NE20_SPEC), zone_fluff )

           enuc = 0.0

           if (sim_ignite) then

              !-----------------------------------------------
              ! initialize ignition point
              !-----------------------------------------------

              ! first find closest ignition point if multiple

              ipi = 1
              ign_dist = (iCoords(i) - sim_ign_x(ipi))**2
              if (NDIM >= 2) ign_dist = ign_dist + (jCoords(j) - sim_ign_y(ipi))**2
              if (NDIM == 3) ign_dist = ign_dist + (kCoords(k) - sim_ign_z(ipi))**2
              ign_dist = sqrt(ign_dist)
              ! check others if more than one
              do testi = 2, sim_ign_numpnts
                 test_dist = (iCoords(i) - sim_ign_x(testi))**2
                 if (NDIM >= 2) test_dist = test_dist + (jCoords(j) - sim_ign_y(testi))**2
                 if (NDIM == 3) test_dist = test_dist + (kCoords(k) - sim_ign_z(testi))**2
                 test_dist = sqrt(test_dist)
                 if ( test_dist < ign_dist ) then
                    ign_dist = test_dist
                    ipi = testi
                 endif
              enddo

              ! find radius of center of ignition point (for hemispherical)
              ign_center_radius = sim_ign_x(ipi)**2
              if (NDIM >= 2) ign_center_radius = ign_center_radius + sim_ign_y(ipi)**2
              if (NDIM == 3) ign_center_radius = ign_center_radius + sim_ign_z(ipi)**2
              ign_center_radius = sqrt(ign_center_radius)

              ! raise temperature if in ignition region
              if (     ( sim_ign_hemispherical &
                         .and. (radius >= ign_center_radius) .and.  (ign_dist <= sim_ign_r(ipi)) ) &
                  .or. ( (.not. sim_ign_hemispherical) .and. (ign_dist <= sim_ign_r(ipi))        ) ) then

                 ign_temp = max( zone_state(EOS_TEMP), sim_ign_temp_center(ipi) - &
                             (sim_ign_temp_center(ipi)-sim_ign_temp_edge(ipi))*ign_dist/sim_ign_r(ipi) )
                 if ( sim_ign_keep_pres ) then
                    ! would just call Eos but it doesn't have a T,P mode
                    call sim_find_ign_state_keeping_pres( ign_temp, zone_state, zone_abund)
                 else
                    zone_state(EOS_TEMP) = ign_temp
                 endif
                 ! to trigger refinement
                 enuc = 1.1*sim_refNogenEnucThresh

              endif

           endif ! sim_ignite

           call Eos(MODE_DENS_TEMP, 1, zone_state, zone_abund)
           

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, zone_state(EOS_DENS))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, zone_state(EOS_TEMP))

           call Grid_putPointData(blockId, CENTER, HE4_SPEC, EXTERIOR, cell, zone_abund(HE4_SPEC))
           call Grid_putPointData(blockId, CENTER, C12_SPEC, EXTERIOR, cell, zone_abund(C12_SPEC))
           call Grid_putPointData(blockId, CENTER, O16_SPEC, EXTERIOR, cell, zone_abund(O16_SPEC))
           call Grid_putPointData(blockId, CENTER, NE20_SPEC, EXTERIOR, cell, zone_abund(NE20_SPEC))

           call Grid_putPointData(blockId, CENTER, ENUC_VAR, EXTERIOR, cell, enuc)

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, cell, 0.0)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, cell, 0.0)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, cell, 0.0)

           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, zone_state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, zone_state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, zone_state(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, zone_state(EOS_GAMC))
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, &
                                       zone_state(EOS_PRES)/(zone_state(EOS_DENS)*zone_state(EOS_EINT))+1.0)

           if (zone_fluff) then
              call Grid_putPointData(blockId, CENTER, FLFF_MSCALAR, EXTERIOR, cell, 1.0)
           else
              call Grid_putPointData(blockId, CENTER, FLFF_MSCALAR, EXTERIOR, cell, 0.0)
           endif
        enddo
     enddo
  enddo
  
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

  return
  
end subroutine Simulation_initBlock






