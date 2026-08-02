! Dean M. Townley 2009
!
! This is the static module data for the Simulation Unit for the SNIa_ddt setup
! Note that some of these are allocatable, and are allocated by Simulation_init()

module Simulation_data
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

!don: add ne20
  real,allocatable,dimension(:),save :: sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_c12_tab, sim_wd_ne20_tab, sim_wd_ne22_tab
  real, save :: sim_wd_dr_inv
  integer, save :: sim_wd_npnts

!don: add ne20
!cf: add sim_presFluff
  real, save :: sim_densFluff, sim_tempFluff, sim_presFluff, sim_xc12Fluff, sim_xne20Fluff, sim_xne22Fluff

  logical, save :: sim_ignite
  real, save :: sim_ignX, sim_ignY, sim_ignZ, sim_ignR

  logical, save :: sim_ignitionFile

!don: add sim_ignShell, a shell ignition flag
  logical, save :: sim_ignMPole, sim_ignMPoleSym, sim_ignShell
  real, save :: sim_ignMPoleA
  integer, save :: sim_ignMpoleMinL, sim_ignMpoleMaxL, sim_ignMpoleSeed
  integer, save :: sim_geom
  real, allocatable, dimension(:,:), save  :: mp_A, mp_delta

  logical, save :: sim_ignSin
  real, save    :: sim_ignSinN, sim_ignSinA

  real,    save           :: sim_laminarWidth

  real, save :: sim_refFluffDensThresh, sim_refFluffDensMargin
  real, save :: sim_refFluffThresh, sim_refFluffMargin
  real, save :: sim_refNogenEnucThresh, sim_refNogenFldtThresh, sim_refNogenMargin
  integer, save :: sim_refFluffDensLevel, sim_refFluffLevel, sim_refNogenLevel
  real, save :: sim_refCentRegionDist
  integer, save :: sim_refCentRegionLevel

  real, save :: sim_vrms_reduced, sim_vrms_center, sim_vrms_Tc, sim_vrms_T0, sim_vrms_alpha
  logical, save :: sim_read_turbfield
  integer, save :: sim_smooth_level
  character(len=4096), save :: sim_turbfield_filename
  real, save :: sim_turbfield_bbox(IAXIS:KAXIS,LOW:HIGH)

end module Simulation_data
