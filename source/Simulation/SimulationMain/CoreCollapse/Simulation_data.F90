!!****if* source/Simulation/SimulationMain/CoreCollapse/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density

!!***

module Simulation_data

  implicit none
  save
#include "constants.h"
#include "Simulation.h"

  integer, save :: nvar_stored
  integer, parameter :: n1d_max = 10000 ! Max number of lines a file can have
  integer, save :: n1d_total ! Actual number of lines, calculated after input
  real, save :: sim_smlrho, sim_smallt,sim_smallx
  character(len=80),save :: model_file
  real, allocatable, save :: xzn(:)
  real,allocatable, save :: model_1d(:,:)
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save  :: sim_windVel, sim_massLoss, sim_velMult
  integer, save :: nsub
  character (len=4), save :: unklabels(UNK_VARS_BEGIN:UNK_VARS_END)

  integer, save :: sim_meshMe

  logical,save :: useCool
  real, save :: coolTemp

  logical, save :: sim_restart, sim_burnUpdateEint
  real, save :: sim_pointMass, sim_expEner, sim_holeRadius

  real, save :: sim_Enu
  real, save :: sim_rhoOne, sim_rhoTwo, sim_rhoThree
  real, save :: sim_yOne, sim_yTwo, sim_yc, sim_yThree

  logical, save :: sim_usePnotT

  real, save :: sim_shockRadTot, sim_shockRadNum
  logical, save :: sim_postBounce
  real, save :: sim_bounceTime
  real, save :: sim_massAccRate, sim_massAccNum

  real, save :: sim_rotOmega, sim_rotA, sim_mage
  real, save :: sim_magB0, sim_magR0

  real, save :: sim_holeRad

  real, save :: sim_maxDens

  logical, save :: sim_mriRefine
  real, save :: sim_mriTime

  logical, save :: sim_usePerturb
  real, save :: sim_perturbRadMin0, sim_perturbRadMax0, sim_perturbMag0
  integer, save :: sim_perturbN0, sim_perturbL0, sim_perturbM0
  real, save :: sim_perturbRadMin1, sim_perturbRadMax1, sim_perturbMag1
  integer, save :: sim_perturbN1, sim_perturbL1, sim_perturbM1

  integer, save :: sim_lrefCenter

  real, save :: sim_velRotFac

  logical, save :: sim_doMRIrefine

  real, save :: sim_avgShockR, sim_minShockR, sim_maxShockR

  logical, save :: sim_1Dprofile_2ms = .false.
  logical, save :: sim_1Dprofile_5ms = .false.
  logical, save :: sim_1Dprofile_10ms = .false.
  logical, save :: sim_1Dprofile_20ms = .false.
  logical, save :: sim_1Dprofile_30ms = .false.

  real, save :: sim_bdryRad, sim_fullAngRefRad

  real, save :: sim_tinitial

  logical, save :: sim_alwaysRefineShock

  logical :: use_randomperts
  integer :: rnd_seed
  real :: rnd_scale

end module Simulation_data
