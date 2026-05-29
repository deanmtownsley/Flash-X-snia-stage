!!
!! Dean M. Townsley 2009
!!
!! Initialization of Simulation Unit for WDHeDetFullStar
!! See source/Simulation/Simulation_init.F90 for API spec and notes.
!! 

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, only : Logfile_stampMessage
  use Driver_data, only : dr_restart
  use Grid_interface, only : Grid_getGeometry
  use Logfile_interface, only : Logfile_stampMessage

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  character(len=256) :: initialWDFileName
  character(len=4096) :: ignitionFileName
  character(len=80) :: reportbuf
  character(len=256) :: burnedKeyProductStr
  integer :: istat, i
  real :: lastradius, radius
  integer :: ignfileunit=2

  integer :: try, l
  logical :: accept
  real    :: u, v, r
  real    :: deltaplus, deltaminus, costheta, P_lm1, P_l, hold

  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get( 'initialWDFile', initialWDFileName)

  call RuntimeParameters_get('fluff_boundary_radius', sim_fluffBoundaryRadius)
  call RuntimeParameters_get('dens_fluff_inner', sim_densFluffInner)
  call RuntimeParameters_get('temp_fluff_inner', sim_tempFluffInner)
  call RuntimeParameters_get('dens_fluff_outer', sim_densFluffOuter)
  call RuntimeParameters_get('temp_fluff_outer', sim_tempFluffOuter)

  call RuntimeParameters_get('refFluffThresh', sim_refFluffThresh)
  call RuntimeParameters_get('refFluffMargin', sim_refFluffMargin)
  call RuntimeParameters_get('refFluffLevel', sim_refFluffLevel)

  call RuntimeParameters_get('refNogenEnucThresh', sim_refNogenEnucThresh)
  call RuntimeParameters_get('refNogenMargin', sim_refNogenMargin)
  call RuntimeParameters_get('refNogenLevel', sim_refNogenLevel)

  call RuntimeParameters_get('refBurnedKeyProduct', burnedKeyProductStr)
  call Simulation_mapStrToInt( burnedKeyProductStr, sim_refBurnedKeyProductIndex, MAPBLOCK_UNK)
  call RuntimeParameters_get( 'refBurnedProductThresh', sim_refBurnedProductThresh)

  call RuntimeParameters_get('refEjectaPhaseStartTime', sim_refEjectaPhaseStartTime )
  call RuntimeParameters_get('refEjectaPhaseMaxRes', sim_refEjectaPhaseMaxRes )


  !----------------------------------------------------------
  ! Now get and set up information about ignition
  !----------------------------------------------------------
  call RuntimeParameters_get('ignite', sim_ignite)
  call RuntimeParameters_get('ign_keep_pres', sim_ign_keep_pres)
  call RuntimeParameters_get('ign_hemispherical', sim_ign_hemispherical)
  call RuntimeParameters_get('ign_file', sim_ign_file)
  if (sim_ign_file) then
     call Logfile_stampMessage('[Simulation_init] Reading ignition points from file')
     call RuntimeParameters_get('ign_file_name', ignitionFileName)
     open(unit=ignfileunit,file=ignitionFileName,status='OLD',iostat=istat)
     if (istat /= 0) call Driver_abortFlash('Unable to open ignition points file')
     ! one-line header ignored
     read(ignfileunit,*)
     read(ignfileunit,*) sim_ign_numpnts
     allocate(sim_ign_x(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_y(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_z(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_r(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_temp_center(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_temp_edge(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     do i = 1, sim_ign_numpnts
        read(ignfileunit,*) sim_ign_x(i), sim_ign_y(i), sim_ign_z(i), sim_ign_r(i), &
                                               sim_ign_temp_center(i), sim_ign_temp_edge(i)
     enddo
     close(ignfileunit)
  else
     ! just a single ignition point, read from parameter file
     sim_ign_numpnts = 1
     allocate(sim_ign_x(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_y(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_z(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_r(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_temp_center(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ign_temp_edge(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     call RuntimeParameters_get('ign_x', sim_ign_x(1))
     call RuntimeParameters_get('ign_y', sim_ign_y(1))
     call RuntimeParameters_get('ign_z', sim_ign_z(1))
     call RuntimeParameters_get('ign_r', sim_ign_r(1))
     call RuntimeParameters_get('ign_temp_center', sim_ign_temp_center(1))
     call RuntimeParameters_get('ign_temp_edge', sim_ign_temp_edge(1))
  endif

  !--------------------------------------------------------
  !  read in 1d initial wd profile
  !--------------------------------------------------------

  if (.not.dr_restart) then

    call Logfile_stampMessage('[Simulation_init] Reading initial 1-d WD profile')
    open(unit=2,file=initialWDFileName,status='OLD',iostat=istat)
    if (istat /= 0) call Driver_abortFlash('Unable to open initial WD profile')

    ! eat header
    read(2,*)
    read(2,*) sim_wd_npnts

    allocate(sim_wd_dens_tab(sim_wd_npnts),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
    allocate(sim_wd_temp_tab(sim_wd_npnts),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
    allocate(sim_wd_he4_tab(sim_wd_npnts),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
    allocate(sim_wd_c12_tab(sim_wd_npnts),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
    allocate(sim_wd_o16_tab(sim_wd_npnts),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
    allocate(sim_wd_ne20_tab(sim_wd_npnts),STAT=istat)
      if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
    lastradius = 0.0
    do i = 1, sim_wd_npnts
       lastradius = radius
       read(2,*) radius, sim_wd_dens_tab(i), sim_wd_temp_tab(i), sim_wd_he4_tab(i), &
                         sim_wd_c12_tab(i), sim_wd_o16_tab(i), sim_wd_ne20_tab(i)
    enddo
    close(2)
    sim_wd_dr_inv = 1.0/ (radius - lastradius)

!    do i = 1, sim_wd_npnts
!      print *, sim_wd_dens_tab(i), sim_wd_temp_tab(i), sim_wd_he4_tab(i), &
!                               sim_wd_c12_tab(i), sim_wd_o16_tab(i), sim_wd_ne20_tab(i)
!    enddo

  endif

end subroutine Simulation_init
