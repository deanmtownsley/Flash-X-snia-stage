! Dean M. Townley 2009
!
! This is the static module data for the Simulation Unit for the SNIa_ddt setup
! Note that some of these are allocatable, and are allocated by Simulation_init()

module Simulation_data
#include "Flash.h"
#include "Eos.h"

  real,allocatable,dimension(:),save :: sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_he4_tab, &
                                        sim_wd_c12_tab, sim_wd_o16_tab, sim_wd_ne20_tab
  real, save :: sim_wd_dr_inv
  integer, save :: sim_wd_npnts

  real, save :: sim_fluffBoundaryRadius, sim_densFluffInner, sim_tempFluffInner, &
                                         sim_densFluffOuter, sim_tempFluffOuter

  logical, save :: sim_ignite, sim_ign_keep_pres, sim_ign_hemispherical, sim_ign_file
  real, allocatable, dimension(:), save :: sim_ign_x, sim_ign_y, sim_ign_z, sim_ign_r
  real, allocatable, dimension(:), save :: sim_ign_temp_center, sim_ign_temp_edge
  integer, save :: sim_ign_numpnts

  ! adaptive mesh marking parameters
  real, save :: sim_refFluffThresh, sim_refFluffMargin
  real, save :: sim_refNogenEnucThresh, sim_refNogenMargin
  integer, save :: sim_refFluffLevel, sim_refNogenLevel

  integer, save :: sim_refBurnedKeyProductIndex
  real, save :: sim_refBurnedProductThresh, sim_maxRadiusBurned

  real, save :: sim_refEjectaPhaseStartTime
  integer, save :: sim_refEjectaPhaseMaxRes


end module Simulation_data
