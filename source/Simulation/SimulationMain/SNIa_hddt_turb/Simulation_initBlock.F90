!!
!! Dean M. Townsley 2009
!!
!! Simulation grid initialization routine for SNIa_ddt setup.
!! See source/Simulation/Simulation_initBlock.F90 for API and notes.
!!
!! This initializes the white dwarf and flame for Type Ia simulation
!! This initialization is generally intended for modelling DDTs
!! Many parameters used here (those in the Simulation Unit static
!! module data area, Simulation_data) were initialized in Simulation_init()
!!
!! See Townsley etal (2009, ApJ, in press, ArXiv:0906.4384) for further discussion of
!! the multipole randomized initial condition implemented here.
!!

subroutine Simulation_initBlock(blockID)
  
  use Simulation_data
  use sim_local_interface, ONLY : sim_interpolate1dWd
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, Grid_putPointData, &
    Grid_getCellCoords
  use Flame_interface, ONLY : Flame_getProfile, Flame_rhJump, Flame_rhJumpReactive
  use bn_paraInterface, ONLY : bn_paraFuelAshProperties
  use Eos_interface, ONLY : Eos

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  
  integer, intent(in) :: blockID

  integer :: i, j, k

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: cell
  integer :: isizeGC, jsizeGC, ksizeGC
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  real, dimension(MDIM) :: deltas

  real, dimension(EOS_NUM) :: state, unburned_state, nse_state
  real :: radius, ign_dist
  real :: P_l_m, P_lm1_m, hold, flame_radius, costheta, sintheta
  real :: theta, phi, fact
  integer :: l, m, n, p, nmax
  real :: fsurf_distance, rshell
!don: added ne20
  real :: flam, xc12initial, xne20initial, xne22initial, ye, dyi_qn, dqbar_qn
  real :: yi
  real :: qbar_nse
  real :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a
  real :: enuc
  logical :: zone_fluff

  real :: cgsMeVperAmu = 9.6485e17

!==============================================================================

  ! get essential info about this block - index limits and cell coordinates
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  isizeGC = blkLimitsGC(HIGH,IAXIS)
  allocate(iCoords(isizeGC))
  jsizeGC = blkLimitsGC(HIGH,JAXIS)
  allocate(jCoords(jsizeGC))
  ksizeGC = blkLimitsGC(HIGH,KAXIS)
  allocate(kCoords(ksizeGC))
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoords,isizeGC)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoords,jsizeGC)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCoords,ksizeGC)

  call Grid_getDeltas(blockID, deltas)

  !-----------------------------------------------
  ! loop over all zones and init
  !-----------------------------------------------
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           !-----------------------------------------------
           !  determine state of material at this radius if unburned
           !  from external 1-d hyrdostatic model
           !-----------------------------------------------

           radius = iCoords(i)**2
           if (NDIM >= 2) radius = radius + jCoords(j)**2
           if (NDIM == 3) radius = radius + kCoords(k)**2
           radius = sqrt(radius)

           call sim_interpolate1dWd(radius, deltas(IAXIS), unburned_state(EOS_DENS), unburned_state(EOS_TEMP), &
                &xc12initial, xne20initial, xne22initial, zone_fluff )
!alan
           !write(*,*)'after interpolation radius dens=',radius,unburned_state(EOS_DENS)
           call bn_paraFuelAshProperties(xc12initial, xne20initial, xne22initial, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

           !print *, xc12initial, xne20initial, xne22initial

           unburned_state(EOS_ABAR) = 1.0/yi_f
           unburned_state(EOS_ZBAR) = ye_f/yi_f

           !print *, unburned_state(EOS_ABAR), unburned_state(EOS_ZBAR)
           if (.not. sim_ignite) then
              ! no burned material, only unburned
              state(:) = unburned_state(:)
              call Eos(MODE_DENS_TEMP, 1, state)
              flam = 0.0
              ye = state(EOS_ZBAR)/state(EOS_ABAR)
              dyi_qn = 0.0
              dqbar_qn = 0.0
              enuc = 0.0
           else
              !-----------------------------------------------
              ! initialize, including a burned region
              !-----------------------------------------------

              ! find distance from flame surface (positive is in front of flame)

              ! default to a spherical region centered at specified coordinates
              ! distance from center of ignition region
              ign_dist = (iCoords(i) - sim_ignX)**2
              if (NDIM >= 2) ign_dist = ign_dist + (jCoords(j) - sim_ignY)**2
              if (NDIM == 3) ign_dist = ign_dist + (kCoords(k) - sim_ignZ)**2
              ign_dist = sqrt(ign_dist)

              flame_radius = sim_ignR
              rshell = 0.0

              if (sim_ignMpole .and. sim_ignSin)  &
                 call Driver_abortFlash("Simulation_initBlock: multipole and sinusoidal ignition are exclusive")

              if (sim_ignMpole) then
                 if (NDIM == 2) then
                    ! assume 2d is cylindrical
                    costheta = (jCoords(j)-sim_ignY)/ign_dist
                    phi = 0.0
                 else if (NDIM == 3) then
                    ! assume 3d is cartesian
                    costheta = (kCoords(k)-sim_ignZ)/ign_dist
                    phi = atan2( jCoords(j)-sim_ignY, iCoords(i)-sim_ignX )
                 endif
                 if (sim_ignMPoleSym) then
                    nmax = 0
                 else
                    nmax = 2*sim_ignMPoleMaxL
                 endif
                 sintheta = sqrt(1.0-costheta*costheta)
                 ! Legendre polynomial
                 do n = 0, nmax, 2
                    m = n / 2
                    do l = m, sim_ignMPoleMaxL
                       if ( l == m ) then
                          P_l_m = 1.0
                          if ( m > 0 ) then
                             fact = 1.0
                             do p = 1, m
                                P_l_m = -P_l_m*fact*sintheta
                                fact = fact + 2.0
                             enddo
                          endif
                          P_lm1_m = 0.0
                       else
                          hold = P_l_m
                          P_l_m = (costheta*(2*l-1)*P_l_m - (l+m-1)*P_lm1_m)/ &
                                  real(l-m)
                          P_lm1_m = hold
                       endif
                          
                       if ( l >= sim_ignMPoleMinL ) then
                          rshell = rshell +  &
                                mp_A(l,n)*P_l_m*cos(m*phi+mp_delta(l,n))
                          ! coefficients for -m are stored in n-1 index.
                          if (m /= 0) rshell = rshell +  &
                                mp_A(l,n-1)*P_l_m*cos(-m*phi+mp_delta(l,n-1))
                       endif
                    enddo
                 enddo
              endif

              if (sim_ignSin) then
                 ! sinusoidal (in theta) border between burned and unburned region
                 ! sim_ignitionSinN = number of periods between theta = 0 and pi
                    ! this is not a good place for this warning (reported for every block)
                    !... need to move
                    !if (myPE==MASTER_PE) then
                    !   write(6,*) "warning: sinusoidal ignition reduces to spherical in 1 dimension"
                    !endif
                 if (NDIM == 2) then
                    theta = acos( (jCoords(j)-sim_ignY)/ign_dist)
                    rshell = rshell - sim_ignSinA*cos(2*theta*sim_ignSinN)
                 else if (NDIM == 3) then
                    theta = acos( (kCoords(k)-sim_ignZ)/ign_dist)
                    rshell = rshell - sim_ignSinA*cos(2*theta*sim_ignSinN)
                 endif
              endif

              if (sim_ignShell) then
                fsurf_distance = abs(ign_dist - flame_radius) - rshell
              else
                flame_radius = flame_radius + rshell
                fsurf_distance = ign_dist - flame_radius
              endif

              ! determine local state in this zone
              ! assume deltas are equal
              if ( (fsurf_distance-0.5*deltas(IAXIS)) > 1.5*sim_laminarWidth ) then
                 ! whole cell unburned material
                 state(:) = unburned_state(:)
                 call Eos(MODE_DENS_TEMP, 1, state)
                 !print *, 'After calling Eos, state is:', state
                 flam = 0.0
                 ye = state(EOS_ZBAR)/state(EOS_ABAR)
                 dyi_qn = 0.0
                 dqbar_qn = 0.0
                 enuc = 0.0
              else if ( (fsurf_distance+0.5*deltas(IAXIS)) < -1.5*sim_laminarWidth ) then
                 ! fully burned to NSE
                 call Flame_rhJumpReactive(unburned_state, qbar_f, state, dqbar_qn, MODE_DENS_TEMP)
                 flam = 1.0
                 dyi_qn   = 1.0/state(EOS_ABAR)
                 ye    = dyi_qn*state(EOS_ZBAR)
                 enuc = 0.0
              else
                 ! partially burned
                 ! at least one cell will fall here (necessary to get initial refinement right)
                 call Flame_getProfile(fsurf_distance, flam)
                 ! calculate properties of NSE final state
                 call Flame_rhJumpReactive(unburned_state, qbar_f, nse_state, qbar_nse, MODE_DENS_TEMP)

                 ! calculate propertise for partially burned material
                 ! note, in fact ye_f and ye_a should be equal
                 yi = yi_f*(1.0-flam) + (1.0/nse_state(EOS_ABAR))*flam
                 ye = ye_f*(1.0-flam) + (nse_state(EOS_ZBAR)/nse_state(EOS_ABAR))*flam
                 state(:) = unburned_state(:)
                 state(EOS_ABAR) = 1.0/yi
                 state(EOS_ZBAR) = ye/yi

                 ! put this in pressure equilibrium with unburned material
                 call Flame_rhJump(unburned_state, state, flam*(qbar_nse-qbar_f)*cgsMeVperAmu, 0.0, MODE_DENS_TEMP)

                 dyi_qn   = flam * 1.0/nse_state(EOS_ABAR)
                 dqbar_qn = flam * qbar_nse

                 ! to trigger refinement
                 enuc = 1.1*sim_refNogenEnucThresh
              endif

           endif ! sim_ignite
           

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, state(EOS_DENS))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, state(EOS_TEMP))

           call Grid_putPointData(blockId, CENTER, FLAM_MSCALAR, EXTERIOR, cell, flam)

           call Grid_putPointData(blockId, CENTER, CI_MSCALAR, EXTERIOR, cell, xc12initial)
           call Grid_putPointData(blockId, CENTER, NE0I_MSCALAR, EXTERIOR, cell, xne20initial)
           call Grid_putPointData(blockId, CENTER, NE2I_MSCALAR, EXTERIOR, cell, xne22initial)

           call Grid_putPointData(blockId, CENTER, PHFA_MSCALAR, EXTERIOR, cell, flam)
           call Grid_putPointData(blockId, CENTER, PHAQ_MSCALAR, EXTERIOR, cell, flam)
           call Grid_putPointData(blockId, CENTER, PHQN_MSCALAR, EXTERIOR, cell, flam)

           call Grid_putPointData(blockId, CENTER, YE_MSCALAR, EXTERIOR, cell, ye)
           call Grid_putPointData(blockId, CENTER, DYQN_MSCALAR, EXTERIOR, cell, dyi_qn)
           call Grid_putPointData(blockId, CENTER, DQQN_MSCALAR, EXTERIOR, cell, dqbar_qn)

           call Grid_putPointData(blockId, CENTER, ENUC_VAR, EXTERIOR, cell, enuc)

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, cell, 0.0)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, cell, 0.0)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, cell, 0.0)

           !  usually I would just call the EOS, but we happen to have all this data
           !  so we'll just put it in.
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, state(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, state(EOS_GAMC))
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, &
                                       state(EOS_PRES)/(state(EOS_DENS)*state(EOS_EINT))+1.0)

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






