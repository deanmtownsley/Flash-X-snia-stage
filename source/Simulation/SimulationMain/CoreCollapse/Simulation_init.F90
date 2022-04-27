!!****if* source/Simulation/SimulationMain/CoreCollapse/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion
!!  problem.
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!  p_ambient       Initial ambient pressure
!!  rho_Ambient     Initial ambient density
!!  exp_energy      Explosion energy (distributed over 2^dimen central zones)
!!  t_init          Initial time since explosion
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***
subroutine Simulation_init()
  use Simulation_data
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Driver_interface, ONLY : Driver_abort, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_getGeometry
  use Tree, ONLY : lrefine_max
  use ut_modelInterface, ONLY : ut_modelRead
#ifdef EOS_HELMNSE
  use eos_hybridData, ONLY : eos_hybTransitionLow, eos_hybTransitionDens
  use Multispecies_interface, ONLY : Multispecies_getSumInv
  use Eos_interface, ONLY : eos 
#include "Eos.h"
#include "Multispecies.h"
#endif    

  implicit none

#include "constants.h"
#include "Simulation.h"

  integer :: i
  integer :: geom

  integer, parameter :: max_stored_vars = 200
  real :: var_temp(max_stored_vars)
#ifdef EOS_HELMNSE  
  real :: hybTrasitionLow, hybTrasitionDens  
  real, dimension(EOS_NUM) :: eosData
#endif  

  logical :: dr_restart

  character (len=256) :: current_line

  integer :: j, ipos, NUNK_VARS_stored
  !  integer :: var_key(max_stored_vars)
  integer :: var_key (NUNK_VARS)
  character (len=4) :: var_labels(max_stored_vars)

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call Grid_getGeometry(geom)

  call RuntimeParameters_get( 'model_file',   model_file)
  call RuntimeParameters_get( 'smlrho',   sim_smlrho)
  call RuntimeParameters_get( 'smallt',   sim_smallt)
  call RuntimeParameters_get( 'smallx',   sim_smallx)
  call RuntimeParameters_get( 'vel_wind', sim_windVel)
  call RuntimeParameters_get( 'mass_loss', sim_massLoss)
  call RuntimeParameters_get( 'vel_mult', sim_velMult)
  call RuntimeParameters_get( 'nsub', nsub)
  call RuntimeParameters_get( 'restart', sim_restart)
  call RuntimeParameters_get( 'burnUpdateEint', sim_burnUpdateEint)
  call RuntimeParameters_get( 'ener_exp', sim_expEner)

  call RuntimeParameters_get( 'point_mass', sim_pointMass)

  call RuntimeParameters_get( 'hole_radius', sim_holeRadius)

  call RuntimeParameters_get( 'rot_omega', sim_rotOmega)
  call RuntimeParameters_get( 'rot_a', sim_rotA)
  call RuntimeParameters_get( 'mag_e', sim_mage)

  call RuntimeParameters_get( 'mri_refine', sim_mriRefine)
  call RuntimeParameters_get( 'mri_time', sim_mriTime)

  call RuntimeParameters_get( 'mag_B0', sim_magB0)
  call RuntimeParameters_get( 'mag_R0', sim_magR0)

  call RuntimeParameters_get( 'hole_radius', sim_holeRad)

!!$  call RuntimeParameters_get("delep_Enu", sim_Enu)
!!$  call RuntimeParameters_get("delep_rhoOne", sim_rhoOne)
!!$  call RuntimeParameters_get("delep_rhoTwo", sim_rhoTwo)
!!$  call RuntimeParameters_get("delep_rhoThree", sim_rhoThree)
!!$  call RuntimeParameters_get("delep_yOne", sim_yOne)
!!$  call RuntimeParameters_get("delep_yTwo", sim_yTwo)
!!$  call RuntimeParameters_get("delep_yThree", sim_yThree)
!!$  call RuntimeParameters_get("delep_yc", sim_yc)

  call RuntimeParameters_get("use_PnotT", sim_usePnotT)

!!$  call RuntimeParameters_get('xmin',sim_xMin)
!!$  call RuntimeParameters_get('ymin',sim_yMin)
!!$  call RuntimeParameters_get('zmin',sim_zMin)
!!$  call RuntimeParameters_get('xmax',sim_xMax)
!!$  call RuntimeParameters_get('ymax',sim_yMax)
!!$  call RuntimeParameters_get('zmax',sim_zMax)

  !! Grab perturbation parameters
  call RuntimeParameters_get("use_perturb", sim_usePerturb)
  call RuntimeParameters_get("perturb_radMin0", sim_perturbRadMin0)
  call RuntimeParameters_get("perturb_radMax0", sim_perturbRadMax0)
  call RuntimeParameters_get("perturb_mag0", sim_perturbMag0)
  call RuntimeParameters_get("perturb_n0", sim_perturbN0)
  call RuntimeParameters_get("perturb_l0", sim_perturbL0)
  call RuntimeParameters_get("perturb_m0", sim_perturbM0)
  call RuntimeParameters_get("perturb_radMin1", sim_perturbRadMin1)
  call RuntimeParameters_get("perturb_radMax1", sim_perturbRadMax1)
  call RuntimeParameters_get("perturb_mag1", sim_perturbMag1)
  call RuntimeParameters_get("perturb_n1", sim_perturbN1)
  call RuntimeParameters_get("perturb_l1", sim_perturbL1)
  call RuntimeParameters_get("perturb_m1", sim_perturbM1)

  call RuntimeParameters_get("lrefine_center", sim_lrefCenter)
  if (sim_lrefCenter == 0) sim_lrefCenter = lrefine_max

  call RuntimeParameters_get("velRotFac", sim_velRotFac)

  call RuntimeParameters_get("bdryRad", sim_bdryRad)
  call RuntimeParameters_get("fullAngRefRad", sim_fullAngRefRad)

  call RuntimeParameters_get("tinitial", sim_tinitial)

  call RuntimeParameters_get("alwaysRefineShock", sim_alwaysRefineShock)
  sim_massLoss = sim_massLoss * 2.0e33 / 31556926.

  ! MRI Refinement flag
  sim_doMRIrefine = .FALSE.

  ! Random perturbation stuff
  call RuntimeParameters_get('use_randomperts', use_randomperts)
  call RuntimeParameters_get('rnd_seed', rnd_seed)
  call RuntimeParameters_get('rnd_scale', rnd_scale)

  do i = UNK_VARS_BEGIN, UNK_VARS_END
     call Simulation_mapIntToStr(i, unklabels(i), MAPBLOCK_UNK)
     call makeLowercase(unklabels(i))
  enddo

  ! Now call the model utility to get the model and coords
  call ut_modelRead(model_file, unklabels, NUNK_VARS, n1d_total, model_1d, xzn)

!!$  do i = 1, n1d_total
!!$     print *, xzn(i), sum(model_1d(i,SPECIES_BEGIN:SPECIES_END))
!!$  enddo

  !! Initialize sim_maxDens for use in Grid_markRefineDerefine
  sim_maxDens = maxval(model_1d(:,DENS_VAR))

#ifdef EOS_HELMNSE
  ! If we are using the Hybrid EOS with a network, we need to compute EINT from input TEMP
  ! using only Helmholtz
  call RuntimeParameters_get('eos_hybTransitionLow', hybTrasitionLow)
  call RuntimeParameters_get('eos_hybTransitionDens', hybTrasitionDens)
  ! cheat by setting the Hybrid EOS low transition density to huge
  eos_hybTransitionLow = HUGE(1.)
  eos_hybTransitionDens = HUGE(1.)
  ! Now loop over the zones in the 1D model and get the EINT for Helmholtz EOS
  ! Note that we require Ye to be in the initial model file for correctness in NSE,
  ! but SUMY/ABAR is NOT guaranteed to be there so we must compute it from the 
  ! input species
  do i=1,n1d_total
    call Multispecies_getSumInv(A, model_1d(i,SUMY_MSCALAR), model_1d(i,SPECIES_BEGIN:SPECIES_END))
    eosData(EOS_DENS) = model_1d(i,DENS_VAR)
    eosData(EOS_TEMP) = model_1d(i,TEMP_VAR)
    eosData(EOS_ABAR) = 1./model_1d(i,SUMY_MSCALAR)
    eosData(EOS_ZBAR) = model_1d(i,YE_MSCALAR)/model_1d(i,SUMY_MSCALAR)
    call Eos(MODE_DENS_TEMP, 1, eosData, massFrac=model_1d(i,SPECIES_BEGIN:SPECIES_END))
    model_1d(i,EINT_VAR) = eosData(EOS_EINT)
    model_1d(i,PRES_VAR) = eosData(EOS_PRES)
  end do
  ! Reset the transition densities to original values
  eos_hybTransitionLow  = hybTrasitionLow
  eos_hybTransitionDens = hybTrasitionDens  
#endif  

end subroutine Simulation_init
