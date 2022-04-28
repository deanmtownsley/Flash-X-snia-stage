!!****if* source/Particles/ParticlesMain/Particles_data
!!
!! NAME
!!    Particles_data
!!
!! SYNOPSIS
!!    Particles_data()
!!
!! DESCRIPTION
!!    Module to hold local variables and data types for Particles unit
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!    useParticles   BOOLEAN [TRUE]  Should particles be used in this simulation?
!!                                    in the setup.
!!    pt_maxPerProc  INTEGER [1000]   Maximum number of particles per processor. Allocates array space
!!                                   Particles are distributed per PROCESSOR rather than per BLOCK
!!    pt_dtFactor    REAL    [0.5]   Factor to make sure that time step is small enough that particles
!!                                   don't move farther than one block in each step
!!    pt_dtChangeTolerance REAL [0.4] For uncorrected Estimated Midpoint propagation (EstiMidpoint):
!!                                    Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities
!!    pt_small       REAL    [1.0E-10] Used for general comparisons of real values 
!!                                   For example, IF (abs(real1 - real2) .lt. pt_small) THEN
!!    pt_numParticlesWanted  INTEGER [100]  for withDensity -- requested number of particles per block
!!
!!***
!*******************************************************************************

#include "Flash.h"

module Particles_data

#ifdef FLASH_GRID_AMREX
 use pt_amrexParticleMod, ONLY : pt_amrexParticlecontainer_t
#endif
!===============================================================================

  implicit none

#include "constants.h"
#include "Particles.h"
#include "GridParticles.h"
!-------------------------------------------------------------------------------

!! The particles data structure is an real array of size (NPART_PROPS,MAX_PARTICLES_PER_PROCESSOR)
!! It is allocated at runtime, as MAX_PARTICLES_PER_PROCESSOR is a runtime parameter = pt_maxPerProc
   real, save, allocatable, dimension(:,:), target :: particles
!! Particle containers array to store the particles of each type.
#ifdef FLASH_GRID_AMREX
   type(pt_amrexParticlecontainer_t), save :: pt_containers(1:NPART_TYPES)
#endif

  integer pt_flashProp2AoSRProp(1:NPART_PROPS) 
#ifdef TYPE_PART_PROP
  integer pt_aoSRProp2FlashProp(1:NPART_PROPS-(2*MDIM+5)) 
#else
  integer pt_aoSRProp2FlashProp(1:NPART_PROPS-(2*MDIM+4)) 
#endif



  integer, save   :: pt_numLocal, pt_numLost
  logical, save   :: pt_keepLostParticles

! Run-time parameters, also described in Particles_init and above under PARAMETERS

  logical, save   :: useParticles = .true.   !if including Particles in simulation, default is true
  real, save      :: pt_dtFactor              ! a multiplying factor for time step limits
  real, save      :: pt_small                 ! a small number for velocity comparisons
  real, save      :: pt_dtChangeTolerance     ! percentage change threshold that
                                              ! controls when some schemes fall back to Euler

! Run-time parameters, local copies defined from some other Unit

  real, save      :: pt_xmin, pt_xmax, pt_ymin, pt_ymax, pt_zmin, pt_zmax ! from Grid/common
!                                  physical domain lower and upper bound in x/y/z
  integer, save    :: pt_geometry                            ! from Grid/common, an integer and string
  character(len=MAX_STRING_LENGTH), save :: pt_str_geometry  ! designation of the geometry configuration

! Run-time parameters, for grid initialization
  integer, save      :: pt_maxPerProc  ! for Lattice


! Guard Cell masks, initialized in Particles_init
  logical, save      :: pt_gcMaskForAdvance(NUNK_VARS) !This one is for Particles_advance
  logical, save      :: pt_gcMaskForWrite(NUNK_VARS+NDIM*NFACE_VARS) !This one for Particles_updateAttributes
  integer,save       :: pt_gcMaskSizeForAdvance=NUNK_VARS
  integer,save       :: pt_gcMaskSizeForWrite=NUNK_VARS+NDIM*NFACE_VARS
! Local variable to control timestepping
  logical, save      :: pt_restart

  integer,save         :: pt_globalMe, pt_globalNumProcs, pt_globalComm
  integer,save         :: pt_meshMe, pt_meshNumProcs, pt_meshComm
  logical, save        :: pt_posInitialized, pt_velInitialized, pt_resetTag
  integer,save         :: pt_logLevel
  integer,save         :: pt_numAtOnce
  integer,dimension(PT_MAX_ATTRIBUTES),save :: pt_mappedAttributes
  integer,dimension(PART_ATTR_DS_SIZE,PT_MAX_ATTRIBUTES), save :: pt_meshVar
  integer, save        :: pt_numMappedAttributes
  integer, save        :: pt_velNumAttrib
  integer, dimension(PART_ATTR_DS_SIZE,MDIM),save :: pt_velAttrib,pt_velPredAttrib
  integer, dimension(MDIM),save :: pt_posAttrib,pt_posPredAttrib
  integer, dimension(PART_TYPE_INFO_SIZE,NPART_TYPES),save :: pt_typeInfo
  


!! Paramters needed only with Lattice initialization
  integer, save      :: pt_numX, pt_numY, pt_numZ
  real, save         :: pt_initialXMin, pt_initialXMax, &
       pt_initialYMin, pt_initialYMax, pt_initialZMin, pt_initialZMax, &
       pt_initialRadius

!! Parameter needed only for With Density initialization

  integer, save      :: pt_pRand, pt_numParticlesWanted 
  real, save         :: pt_totalMass, pt_totalVolume, pt_averageDensity
  integer, dimension(GRPT_ALL),save      :: pt_indexList
  integer, save      :: pt_indexCount

  logical, save :: pt_reduceGcellFills

!! Parameter needed for generating the unique consecutive tags in parallel
  integer, save :: pt_startTagNumber=0

!! Parameters needed for IMC
  logical, save :: pt_is_eff_scattering, pt_is_grey
  integer, save    :: pt_maxnewnum, pt_max_rt_iterations
  real, save    :: pt_es_alpha
  real, save    :: pt_smlpush
  logical, save :: pt_half_rt_timesteps

  integer, parameter :: pt_STAY_ID =  1
  integer, parameter :: pt_SCAT_ID =  2
  integer, parameter :: pt_ESCAT_ID = 3
  integer, parameter :: pt_PION_ESCAT_ID = 4
  integer, parameter :: pt_ABS_ID =   5
  integer, parameter :: pt_CROSS_ID = 6

  integer, parameter :: pt_quadrant_I =   1
  integer, parameter :: pt_quadrant_II =  2
  integer, parameter :: pt_quadrant_III = 3
  integer, parameter :: pt_quadrant_IV =  4

  logical, save :: pt_FaceEmission
  logical, save :: pt_is_FacePlanck
  integer, save :: pt_FaceEmissionSide, pt_FaceEmissionAxis
  real, save    :: pt_FacePlanckTemp, pt_constFaceFlux
  logical, save :: pt_is_therm_face_vel
  logical, save :: pt_is_iso_face_vel
  logical, save :: pt_is_radial_face_vel
  integer, save :: pt_num_fmcps_tstep

  logical, save :: pt_PointEmission, pt_PointPulse
  real, save    :: pt_PointPulseErad
  real, save    :: pt_PointLuminosity, pt_PointSrcPosOffset
  integer, save :: pt_num_pmcps_tstep
  integer, save :: originprocID, originblkID

  ! use integer switch
  integer, save :: pt_initradfield_num

  logical, save :: pt_ThermalEmission
  logical, save :: pt_marshak_eos
  real, save    :: pt_grey_abs_opac, pt_grey_sca_opac, pt_grey_eps
  real, save    :: pt_dens_threshold
  real, save    :: pt_abs_threshold
  integer, save :: pt_num_tmcps_tstep

  real, save    :: pt_energy_min_eV, pt_energy_max_eV
  real, save    :: pt_energy_min_erg, pt_energy_max_erg

  logical, save :: pt_is_corrdl
  logical, save :: pt_is_deposit_urad, pt_is_deposit_energy,&
                   pt_is_deposit_momentum
  logical, save :: pt_is_thermally_coupled, pt_is_dynamically_coupled
  logical, save :: pt_is_scat_elastic, pt_is_escat_elastic
  logical, save :: pt_is_scat_iso, pt_is_escat_iso

  real, save    :: pt_temp_floor

  logical, save    :: pt_is_veldp

! for spectral transport
  real, dimension(THORNADO_NE) :: rt_Emid, rt_Etop, rt_Ebottom ! erg
  real, dimension(THORNADO_NE) :: rt_int4piE2dE ! erg^3
  real, dimension(THORNADO_NE) :: rt_int4piE3dE ! erg^4

  real, save    :: rt_xMin, rt_xMax, rt_yMin, rt_yMax, rt_zMin, rt_zMax
  real, save    :: rt_gamma

  real, save :: mH, kB, h_planck, clight, sigma, R, G_grav, a_rad
  real, save    :: ev2erg = 1.602E-12
  real, save    :: sigma_pi = 6.3d-18

end module Particles_data
