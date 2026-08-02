!!
!! Dean M. Townsley 2009
!!
!! Initialization of Simulation Unit for SNIa_ddt.
!! See source/Simulation/Simulation_init.F90 for API spec and notes.
!! 
!! This routine does the following:
!! Read in initial WD and store in Simulation Unit static module data area.
!! Retrieve refinement parameters from parameter file datastructures.
!! Retrieve parameters for initial condition and do further setup if necessary.
!!
!! If a randomized initial condition is chosen, the coefficients,
!! which are global, are calculated here.  They are simply calculated on
!! every processor because they are fairly cheap.
!! This initial condition is the randomized multipole used to generate
!! the sample in Townsley etal (2009, ApJ, in press, ArXiv:0906.4384).

subroutine Simulation_init()

  use Simulation_data
  use Flame_interface, ONLY : Flame_rhJumpReactive, Flame_getWidth, Flame_laminarSpeed
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use bn_paraInterface, ONLY : bn_paraFuelAshProperties
  use Logfile_interface, only : Logfile_stampMessage
  use Grid_interface, only : Grid_getGeometry
  use sim_local_interface, only : sim_LCGRandomIterate
  use Logfile_interface, only : Logfile_stampMessage

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  character(len=256) :: initialWDFileName
  character(len=80) :: reportbuf
  integer :: istat, i, j, jstart, jend, k
  real :: lastradius, radius

  integer :: try, l, m, sgn
  logical :: accept
  real    :: u, v, r, ir
  real    :: deltaplus, deltaminus, costheta, sintheta, phi
  real    :: P_lm1_m, P_l_m, hold
  real    :: factlmm, factlpm, fact
  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get( 'initialWDFile', initialWDFileName)

  call RuntimeParameters_get('dens_fluff', sim_densFluff)
  call RuntimeParameters_get('temp_fluff', sim_tempFluff)
!cf: add presFluff
  call RuntimeParameters_get('pres_fluff', sim_presFluff)
  call RuntimeParameters_get('temp_fluff', sim_tempFluff)
  call RuntimeParameters_get('xc12_fluff', sim_xc12Fluff)
!don: add ne20
  call RuntimeParameters_get('xne20_fluff', sim_xne20Fluff)
  call RuntimeParameters_get('xne22_fluff', sim_xne22Fluff)

  call RuntimeParameters_get('ignite', sim_ignite)
  call RuntimeParameters_get('x_match', sim_ignX)
  call RuntimeParameters_get('y_match', sim_ignY)
  call RuntimeParameters_get('z_match', sim_ignZ)
  call RuntimeParameters_get('r_match', sim_ignR)
  
  call RuntimeParameters_get('refFluffDensThresh', sim_refFluffDensThresh)
  call RuntimeParameters_get('refFluffDensMargin', sim_refFluffDensMargin)
  call RuntimeParameters_get('refFluffDensLevel', sim_refFluffDensLevel)

  call RuntimeParameters_get('refFluffThresh', sim_refFluffThresh)
  call RuntimeParameters_get('refFluffMargin', sim_refFluffMargin)
  call RuntimeParameters_get('refFluffLevel', sim_refFluffLevel)

  call RuntimeParameters_get('refNogenEnucThresh', sim_refNogenEnucThresh)
  call RuntimeParameters_get('refNogenFldtThresh', sim_refNogenFldtThresh)
  call RuntimeParameters_get('refNogenMargin', sim_refNogenMargin)
  call RuntimeParameters_get('refNogenLevel', sim_refNogenLevel)

  call RuntimeParameters_get('refCentRegionDist', sim_refCentRegionDist)
  call RuntimeParameters_get('refCentRegionLevel', sim_refCentRegionLevel)

!   write(*,*) 'alan: in Simulation_init'

  ! only need to get width of artificial flame once
  call Flame_getWidth(sim_laminarWidth)

  !--------------------------------------------------------
  !  read in 1d initial wd profile
  !--------------------------------------------------------

  call Logfile_stampMessage('[Simulation_init] Reading initial 1-d WD profile')
  open(unit=2,file=initialWDFileName,status='OLD',iostat=istat)
  if (istat /= 0) call Driver_abortFlash('Unable to open initial WD profile')

  ! eat header
  read(2,*)
  read(2,*) sim_wd_npnts
!alan
!  write(*,*)'no pts in model =',sim_wd_npnts

  allocate(sim_wd_dens_tab(sim_wd_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_temp_tab(sim_wd_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_c12_tab(sim_wd_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
!don: add ne20
  allocate(sim_wd_ne20_tab(sim_wd_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  allocate(sim_wd_ne22_tab(sim_wd_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for inital WD")
  radius = 0.0
  do i = 1, sim_wd_npnts
     lastradius = radius
!don: add ne20, note that the order of the columns has ne20 between c12 and ne22
     read(2,*) radius, sim_wd_dens_tab(i), sim_wd_temp_tab(i), sim_wd_c12_tab(i), sim_wd_ne20_tab(i), sim_wd_ne22_tab(i)
  enddo
  close(2)
  sim_wd_dr_inv = 1.0/ (radius - lastradius)
!alan
  !write(*,*)'after read, sim_wd_dr_inv =',sim_wd_dr_inv
  !pause "after read"
  !write(*,*)'after read, sim_wd_dens_tab =',sim_wd_dens_tab

  call RuntimeParameters_get('ignMPole', sim_ignMpole)
  call RuntimeParameters_get('ignMPoleA', sim_ignMpoleA)
  call RuntimeParameters_get('ignMPoleMinL', sim_ignMpoleMinL)
  call RuntimeParameters_get('ignMPoleMaxL', sim_ignMpoleMaxL)
  call RuntimeParameters_get('ignMPoleSym', sim_ignMPoleSym)
  call RuntimeParameters_get('ignMPoleSeed', sim_ignMPoleSeed)
  !don: read in shell ignition flag
  call RuntimeParameters_get('ignShell', sim_ignShell)
  write (*,*) 'Read parameter sim_ignShell=',sim_ignShell

  if (sim_ignMPole) then
     ! set up coefficients of multipole contributions
     ! check some things first to catch mistakes
     if (NDIM < 2) call Driver_abortFlash("[Simulation_init] multipole ignition not sensible in 1-dimension")
     call Grid_getGeometry(sim_geom)

     if (NDIM == 2) then
        if (sim_geom .ne. CYLINDRICAL) &
           call Driver_abortFlash("[Simulation_init] multipole ignition assumes cylindrical geometry in 2D")
        ! in 2D the ignition condition must be symmetric
        sim_ignMPoleSym = .true.
     endif

     if ( (NDIM == 3) .and. (sim_geom .ne. CARTESIAN) ) &
        call Driver_abortFlash("[Simulation_init] multipole ignition assumes cartesian geometry in 3D")

     ! the j-index is for m values ordered 0,-1,1,-2,2...
     jstart = 0
     if (sim_ignMPoleSym) then
        jend = 0 ! only m=0 for symmetric
     else
        jend = 2*sim_ignMpoleMaxL ! only support max_m = max_l
     endif

     allocate(mp_A(sim_ignMpoleMinL:sim_ignMpoleMaxL,jstart:jend), stat=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate mp_A in Simulation_init")
     allocate(mp_delta(sim_ignMpoleMinL:sim_ignMpoleMaxL,jstart:jend), &
           stat=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate mp_delta in Simulation_init")

     try = 0
     accept = .false.
     do while ( .not. accept )

        try = try + 1

        ! generate set of random coefficients
        do j = jstart, jend, 2
           m = j / 2
           do i = sim_ignMpoleMinL, sim_ignMpoleMaxL
              ! generate a normally distributed number, r
              ! via the Box-Muller method
              call sim_LCGRandomIterate(sim_ignMPoleSeed)
              u = real(sim_ignMPoleSeed)/2147483646.0
              ! avoid taking log of zero
              if (u==0) u = 1
              call sim_LCGRandomIterate(sim_ignMPoleSeed)
              v = real(sim_ignMPoleSeed)/2147483646.0
              r = sqrt(-2*log(u))*cos(2*PI*v)
              if (.not. sim_ignMPoleSym) then
                 ! need a complex coefficient, uniformly distributed
                 call sim_LCGRandomIterate(sim_ignMPoleSeed)
                 u = real(sim_ignMPoleSeed)/2147483646.0
                 if (u==0) u = 1
                 ! map into range -pi to pi
                 mp_delta(i,j) = 2*PI*u - PI
              else
                 mp_delta(i,j) = 0.0
              endif

              
              ! now set amplitude of this component, this will be multiplied by
              ! the legendre polynomial of order i to construct the initial 
              ! flame surface
              ! normalization is chosen to match spherical harmonics as
              ! given in Jackson.
              factlmm = 1.0   !!  (l - m)!
              if ( i - m > 0 ) then
                 do k = 1, i - m
                    factlmm = factlmm*k
                 enddo
              endif

              factlpm = 1.0   !!  (l + m)!
              if ( i + m > 0 ) then
                 do k = 1, i + m
                    factlpm = factlpm*k
                 enddo
              endif

              fact = sqrt((2*i+1)/4.0/PI * factlmm/factlpm)

              mp_A(i,j) = r*sim_ignMpoleA * fact

              ! Need to account for 3D, there are 2l+1 terms as opposed
              ! to just one.
              if (.not. sim_ignMPoleSym) mp_A(i,j) = mp_A(i,j) / sqrt(real(2*i+1))


              ! this will only happen for sim_ignMPoleSym = .false.
              if (m > 0) then
                 call sim_LCGRandomIterate(sim_ignMPoleSeed)
                 u = real(sim_ignMPoleSeed)/2147483646.0
                 if (u==0) u = 1
                 call sim_LCGRandomIterate(sim_ignMPoleSeed)
                 v = real(sim_ignMPoleSeed)/2147483646.0
                 r = sqrt(-2*log(u))*cos(2*PI*v)
                 ! need a complex coefficient
                 call sim_LCGRandomIterate(sim_ignMPoleSeed)
                 u = real(sim_ignMPoleSeed)/2147483646.0
                 if (u==0) u = 1
                 ! this is for negative m which is equivalent to the
                 ! complex conjugate (c.c.) of positive m up to a sign.
                 ! the coefficient itself is not c.c.'ed but it will
                 ! be multiplied by the c.c.'ed part of Ylm.
                 ! Account for sign change here.
                 mp_delta(i,j-1) = PI - 2*PI*u

                 if (mod(m,2) == 0) then 
                    sgn = 1
                 else 
                    sgn = -1
                 endif

                 mp_A(i,j-1) = sgn * r*sim_ignMpoleA * fact

                 ! Need to account for 3D, there are 2l+1 terms as opposed
                 ! to just one.
                 mp_A(i,j-1) = mp_A(i,j-1) / sqrt(real(2*i+1))
              endif

           ! another initialization which turned out to be non-stardard
           ! this was used for some early testing and is just kept for reference
           !call sim_LCGRandomIterate(sim_ignMPoleSeed)
           ! random phase angle, between 0 and pi since we will not use imaginary part
           !alpha = real(sim_ignMPoleSeed)/2147483646.0 * PI
           ! coefficients have amplitudes of equal complex amplitude, with random
           ! phase.  We then use real part.  Normalization is for spherical
           ! harmonics in Jackson
           !mp_A(i) = sqrt((2*i+1)/4.0/PI)*cos(alpha)*sim_ignMpoleA*2/(sim_ignMpoleMaxL-sim_ignMpoleMinL+1)

           !print *, i,alpha/PI, mp_A(i)
           enddo
        enddo

        ! now check that this is negative on + and - axis to help prevent pathologies
        if (sim_ignMPoleSym) then

           deltaplus = 0.0
           costheta = 1.0
           sintheta = 0.0

           ! Legendre polynomial
           m = 0 ! only possible value for symmetric
           do l = m, sim_ignMPoleMaxL
              if ( l == m ) then
                 P_l_m = 1.0
                 P_lm1_m = 0.0
              else
                 hold = P_l_m
                 P_l_m = (costheta*(2*l-1)*P_l_m - (l+m-1)*P_lm1_m)/real(l-m)
                 P_lm1_m = hold
              endif

              if ( l >= sim_ignMPoleMinL ) then
                    deltaplus = deltaplus + mp_A(l,m)*P_l_m
              endif
           enddo

           deltaminus = 0.0
           costheta = -1.0
           sintheta = 0.0
           ! Legendre polynomial
           do l = m, sim_ignMPoleMaxL
              if ( l == m ) then
                 P_l_m = 1.0
                 P_lm1_m = 0.0
              else
                 hold = P_l_m
                 P_l_m = (costheta*(2*l-1)*P_l_m - (l+m-1)*P_lm1_m)/real(l-m)
                 P_lm1_m = hold
              endif

              if ( l >= sim_ignMPoleMinL ) then
                 deltaminus = deltaminus + mp_A(l,m)*P_l_m
              endif
           enddo

           if (deltaplus < 0.0 .and. deltaminus < 0.0) accept = .true.

        else ! not symmetric so accept as is

           accept = .true.

        endif

     enddo


     ! need to report ending seed for statistically independent follow-on run
     write(reportbuf, *) '[Simulation_init] Multipole ignition coefficients initialized'
     call Logfile_stampMessage(reportbuf)
     write(reportbuf, *) '[Simulation_init] required iterations: ', try
     call Logfile_stampMessage(reportbuf)
     write(reportbuf, *) '[Simulation_init] endpoint ignMPoleSeed = ', sim_ignMPoleSeed
     call Logfile_stampMessage(reportbuf)
  endif

!   write(*,*)'alan: in end of Simulation_init'

  call RuntimeParameters_get('ignSin', sim_ignSin)
  call RuntimeParameters_get('ignSinN', sim_ignSinN)
  call RuntimeParameters_get('ignSinA', sim_ignSinA)

  call RuntimeParameters_get('sim_smooth_level', sim_smooth_level)
  call RuntimeParameters_get('sim_vrms_center', sim_vrms_center)
  call RuntimeParameters_get('sim_vrms_Tc', sim_vrms_Tc)
  call RuntimeParameters_get('sim_vrms_T0', sim_vrms_T0)

  ! calculate factors for scaling velocity field
  sim_vrms_reduced = sim_vrms_center / 10**(0.8)
  sim_vrms_alpha = log(10**(0.8)) / log( (sim_vrms_Tc-sim_vrms_T0)/sim_vrms_T0*2 )

  call RuntimeParameters_get('sim_turbfield_filename',sim_turbfield_filename)
  call RuntimeParameters_get('sim_read_turbfield', sim_read_turbfield)
  call RuntimeParameters_get('sim_turbfield_xmin',sim_turbfield_bbox(IAXIS,LOW) )
  call RuntimeParameters_get('sim_turbfield_xmax',sim_turbfield_bbox(IAXIS,HIGH) )
  call RuntimeParameters_get('sim_turbfield_ymin',sim_turbfield_bbox(JAXIS,LOW) )
  call RuntimeParameters_get('sim_turbfield_ymax',sim_turbfield_bbox(JAXIS,HIGH) )
  call RuntimeParameters_get('sim_turbfield_zmin',sim_turbfield_bbox(KAXIS,LOW) )
  call RuntimeParameters_get('sim_turbfield_zmax',sim_turbfield_bbox(KAXIS,HIGH) )

end subroutine Simulation_init
