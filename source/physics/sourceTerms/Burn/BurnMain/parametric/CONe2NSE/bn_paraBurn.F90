! Dean Townsley 2007-8

subroutine bn_paraBurn(dens, temp, eint, pres, &
                       xc12init, xne22init, &
                       flame, flamedot, phi_fa, phi_aq, phi_qn, &
                       ye, dyi_qn, dqbar_qn, &
                       qdot, edotnu, dt, react_proximity, shock, &
                       ignite_detonation, phi_fa_det)
  
  use NSE_interface, ONLY:  NSE_finalAtDens, NSE_finalAtPres  
  use Eos_interface, ONLY: Eos_vector
  use Eos_data, ONLY: eos_smallt
  use bn_paraInterface, ONLY: bn_paraFuelAshProperties

  use Burn_data, ONLY : &
       c_l, N_A, m_p, m_n, m_e, bn_thermalReact, &
       yi_he4, yi_si28, q_si28, q_he4


  implicit none
  
#include "Simulation.h"
#include "constants.h"
#include "Eos.h"
  
  real, intent(in)    :: dens, temp, eint, pres
  real, intent(in)    :: xc12init, xne22init
  real, intent(in)    :: flame, flamedot
  real, intent(inout) :: phi_fa, phi_aq, phi_qn
  real, intent(inout) :: ye, dyi_qn, dqbar_qn
  real, intent(out)   :: qdot, edotnu
  real, intent(in)    :: dt, react_proximity, shock
  logical, intent(in) :: ignite_detonation
  real, intent(in)    :: phi_fa_det

  
  real, parameter :: inflame_threshold = 1.e-6
  real :: cgsMeVperGram = 9.6485e17

  real    :: dti
  
  ! for saving some initial states
  real :: phi_fa_i, phi_aq_i, phi_qn_i
  real :: eint_i, pres_i
  real :: qbar_i, ye_i, yi_i, dye_n

  real :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a

  ! for calculating time evolution of properties
  real :: phi_aq_dot, phi_qn_dot
  real :: r_c
  real :: yi_q_est
  real :: tau_nsqe_inv, tau_nse_inv, maxExp, testExp  
  real :: expq

  real :: qbar
  
  ! for predictions of nse final state
  real :: s
  real :: qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d
  real :: qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p
  real :: qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot

  ! for calculation of reaction rate
  real :: nasigmav, t9, t9a

  ! for estimating reaction temperatures
  real :: qbar_unburned, yi_unburned, qbar_burned, yi_burned, x_bnf, qbar_noflame
  real :: yi_noflame, ye_noflame, eint_noflame
  real,dimension(1, EOS_NUM) :: eosData_2D

  !-------------------------------------------------------------------------------------

  if (shock == 1) then
     qdot = 0.0
     edotnu = 0.0
     return
  endif

  !--------------------------------------------------------
  ! 0.  initialize some information
  !     including getting unburned properties and saving some initial data
  !--------------------------------------------------------
  !inverse timestep
  dti    = 1.e0/dt
  call bn_paraFuelAshProperties(xc12init,xne22init, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

  ! truncate advection errors which might have pushed out of limits
  phi_fa = max(0.0,min(1.0,phi_fa))
  phi_aq = max(0.0,min(phi_fa,phi_aq))
  phi_qn = max(0.0,min(phi_aq,phi_qn))

  dyi_qn = max(0.0,min(1.0,dyi_qn))
  
  !initial plasma properties
  qbar_i  = (1.0-phi_fa)*qbar_f + (phi_fa-phi_aq)*qbar_a + dqbar_qn
  ye_i  = ye
  yi_i  = (1.0-phi_fa)*yi_f + (phi_fa-phi_aq)*yi_a + dyi_qn

  dye_n = max( ye - (1.0-phi_qn)*ye_a, 0.0 )

  pres_i  = pres
  eint_i  = eint
  phi_fa_i  = phi_fa
  phi_aq_i  = phi_aq
  phi_qn_i  = phi_qn

  !--------------------------------------------------------------------------
  ! Calculate properties (qbar, yi, temp) of NSE final state
  !--------------------------------------------------------------------------
  ! We predict the properties of the final NSE state from the current 
  ! (partially burned) state, assuming the burn will occur at either
  ! constant density (isochoric) or constant pressure (isobaric).
  ! The properties have been pre-calculated and tabulated for speed.  Inside a
  ! flame we use an interpolation between the isobaric and isochoric predictions
  ! In addtion to qbar, yi, and temp of the final state, we also obtain
  ! edot (neutrino losses) and yedot (neutronization rate) for that state.
  ! These are used in the evolution of partially burned zones.
  !--------------------------------------------------------------------------
  if ( flame >= 0.9999 .or. flame < inflame_threshold ) then
     ! isochoric prediction
     call NSE_finalAtDens(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
       ye_i, dens, eint_i - qbar_i*cgsMeVperGram)
  
  else if (flame >= 0.99) then
     ! crossover region, linearly average
     call NSE_finalAtDens(qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d, &
            ye_i, dens, eint_i - qbar_i*cgsMeVperGram)
     call NSE_finalAtPres(qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p, &
          ye_i, pres_i, eint_i + pres_i/dens - qbar_i*cgsMeVperGram)
     s = (flame - 0.99)/(0.9999-0.99)
     qbar_finnse = s*qbar_finnse_d + (1.0-s)*qbar_finnse_p
     sumyi_finnse = 1.0/(s/sumyi_finnse_d + (1.0-s)/sumyi_finnse_p)
     temp_finnse = s*temp_finnse_d + (1.0-s)*temp_finnse_p
     edot = s*edot_d + (1.0-s)*edot_p
     yedot = s*yedot_d + (1.0-s)*yedot_p
       
  else
     ! isobaric prediction
     call NSE_finalAtPres(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
          ye_i, pres_i, eint_i + pres_i/dens - qbar_i*cgsMeVperGram)
  endif

  !-------------------------------------------------------------------
  !   calculate timescale calibrated in Calder etal 2007
  !-------------------------------------------------------------------

  testExp = 182.06e9/temp_finnse - 46.054e0
  ! protect from overflow for portability
  maxExp = log(HUGE(1.0))
  if (testExp .GE. maxExp) then
     tau_nsqe_inv = 0.0
  else
     tau_nsqe_inv = exp(-testExp)
  endif
  ! ajk     log(tau_nse)  = (dens**0.2e0)*exp(179.7e9/temp - 40.5e0)
  ! fang1   log(tau_nse)  = (dens**0.2e0)*exp(207.76e9/temp - 47.262e0)
  ! calder etal      log(tau_nse)  = exp(196.02e9/temp_finnse - 41.645e0)
  ! townsley etal '09 tau_nse = exp(-47.36 + 197.7/T_9)
  if ( flame >= 0.9999 .or. &
       ( bn_thermalReact .and. &
          ( flame < inflame_threshold .or. (phi_fa > flame + 0.1) ) ) ) then
     ! use local temperature to get NSE timescale outside flame
     ! (particularly in detonations)
!     testExp = 196.02e9/temp - 41.645e0
!si     testExp = 197.7e9/temp - 47.36e0
     testExp = 201.0e9/temp - 46.77e0
  else
     ! in flame with no apparent thermal burning, use estimated final temperature
     ! this is just a rate, so a discontinuous crossover shouldn't be a big deal
!     testExp = 196.02e9/temp_finnse - 41.645e0
!si     testExp = 197.7e9/temp_finnse - 47.36e0
     testExp = 201.0e9/temp_finnse - 46.77e0
  endif
  if (testExp .GE. maxExp) then
     tau_nse_inv = 0.0
  else
     tau_nse_inv = exp(-testExp)
  endif
 
  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
  ! if we are in the flame, only include thermal reactions if there is nearby
  ! thermal burning.  This prevents flame burning from causing thermal reactions
  ! without outside influence (such as an arriving detonation).
  if (bn_thermalReact .and. &
      .not. ( flame >= inflame_threshold .and. react_proximity > 1.0 ) ) then
     ! away from flame or in flame near other reactions, use thermal reactions

     ! first need a temperature
     if (flame < inflame_threshold .or. phi_fa > (1.e0-inflame_threshold) ) then
        ! away from flame or with nearly fully burned material ...
        ! burn based on nuclear reaction rate determined by t9
        t9 = temp/1.e9
     else
        ! in flame but not yet fully burned...
        ! estimate temperature of material if flame had not burned

        ! catch up corner case in which the flame just started burning
        !   (prevents a divide by zero risk below)
        phi_fa = max(phi_fa,inflame_threshold)
        ! x_bnf = fraction of material which was burned, but not by the flame
        x_bnf = max(0.e0, phi_fa-flame)
        qbar_unburned =  qbar_f
        yi_unburned =  yi_f
        qbar_burned = ( qbar_i - qbar_unburned*(1.0-phi_fa) ) / phi_fa
        yi_burned = ( yi_i - yi_unburned*(1.0-phi_fa) ) / phi_fa
        qbar_noflame = qbar_unburned*(1.e0-x_bnf) + qbar_burned*x_bnf
        yi_noflame = yi_unburned*(1.e0-x_bnf) + yi_burned*x_bnf
        !ye_noflame   = 0.5*(1-x_bnf) + ye_oldash_i*x_bnf
        ! just use current ye instead
        !ye_noflame = ye_i
        ye_noflame = ye_f
        eint_noflame = eint_i - cgsMeVperGram*(qbar_i-qbar_noflame)

        if (yi_noflame < 0.e0) then
           write (6,*) ' neg sumy', phi_fa, flame, dqbar_qn
        endif
        ! find a floor energy for local material
        eosData_2D(1, EOS_DENS) = dens
        eosData_2D(1, EOS_TEMP) = eos_smallt
        eosData_2D(1, EOS_ABAR) = 1.e0/yi_noflame
        eosData_2D(1, EOS_ZBAR) = ye_noflame*eosData_2D(1, EOS_ABAR)
        call Eos_vector(MODE_DENS_TEMP,1,eosData_2D)
        if ( eint_noflame < eosData_2D(1, EOS_EINT) ) then
           t9 = eos_smallt/1.e9
        else
           eosData_2D(1, EOS_EINT) = eint_noflame
           eosData_2D(1, EOS_TEMP) = temp ! guess value
           call Eos_vector(MODE_DENS_EI,1,eosData_2D)
           t9 = eosData_2D(1, EOS_TEMP)/1.e9
        endif
     endif

     ! thermally burn based on temperature just estimated
     ! From Caughlin & Fowler 1988, ADNDT 40 283
     ! note this is UNSCREENED, should be screened at some point
     t9a = t9/(1.0+0.0396*t9)
     nasigmav = 4.27e26*t9a**(5.0/6)*t9**(-1.5)*exp(-84.165*t9a**(-1.0/3.0)-2.12e-3*t9**3)
     ! integrate
     !     phi_fa_dot = dens*xc12init*(1.0-phi_fa)**2/12.0*nasigmav
     ! and add as evolution after flame
     r_c = dens*xc12init/12.0*nasigmav
     phi_fa = 1.0 - (1.0-phi_fa)/(1.0 + r_c*dt*(1.0-phi_fa))

     ! ignition override
     if (ignite_detonation) then
        ! fully burn this zone in this timestep (will hit limiter)
        phi_fa = phi_fa_det
     endif

  endif

  ! follow flame variable in any case
  ! (second operator split step if thermal burning is activated)
  phi_fa = max(0.0,min(1.e0,phi_fa+max(0.0e0,flamedot)*dt))
 

  !------------------------------------------------------
  ! 3 Update progress variables by direct integration, assuming
  !   phi_fa is constant (i.e. evolves much faster than others)

  expq = exp(-dt*tau_nsqe_inv)

  ! integration of     phi_aq_dot = (phi_fa-phi_aq)/tau_nsqe
  phi_aq = phi_fa - (phi_fa-phi_aq_i)*expq
  phi_aq = max(0.0,min(phi_fa,phi_aq))
  phi_aq_dot = (phi_aq-phi_aq_i)*dti

  ! integration of     phi_qn_dot = (phi_aq-phi_qn)^2/tau_nsqe
  ! assuming phi_aq is constant (i.e. operator split)
  ! need to handle case when phi_aq-phi_qn_i is zero
  if ((phi_aq-phi_qn_i) == 0.0) then
     phi_qn = phi_aq
  else
     phi_qn = phi_aq - 1.0/( 1.0/(phi_aq-phi_qn_i) + tau_nse_inv*dt )
  endif
  phi_qn = max(0.0,min(phi_aq,phi_qn))
  phi_qn_dot = (phi_qn-phi_qn_i)*dti

  ! now update partial properties
  ! operator split to mirror later treatment of dyi_qn
  dqbar_qn = dqbar_qn + phi_aq_dot*dt*236.537/28
  dqbar_qn = (phi_aq*236.537/28+phi_qn*(qbar_finnse-236.537/28))*(1.0-expq) + dqbar_qn*expq

  ! these estimating forms don't seem to work very well
!  yi_q_est = max( yi_si28, min( yi_he4, &
!             ( (q_si28 - qbar_finnse)*yi_he4 + (qbar_finnse-q_he4)*yi_si28) ) &
!              / (q_si28-q_he4) )
!  yi_q_est = max( yi_si28, min( yi_he4, &
!             ( (sumyi_finnse - 1.0/56.0)*yi_he4 + (yi_he4-sumyi_finnse)*yi_si28) ) &
!              / (yi_he4-1.0/56.0) )
  yi_q_est = yi_si28

  ! now update dyi_qn
  ! assuming phi_qn is constant (operator split again)
  ! two steps, first change due to change of phi_aq, then evolution
  !  better considering operator splitting of phi_qn and phi_aq evolution
  dyi_qn = dyi_qn + phi_aq_dot*dt*yi_q_est
  dyi_qn = (phi_aq*yi_q_est+phi_qn*(sumyi_finnse-yi_q_est))*(1.0-expq) + dyi_qn*expq

  ! implicit update not necessary, but use actual change in phi_qn from above
  !   ( this is like operator splitting it)
  dye_n = dye_n + ye_a*phi_qn_dot*dt + phi_qn*yedot*dt

  ! Construct final state
  qbar  = (1.0-phi_fa)*qbar_f + (phi_fa-phi_aq)*qbar_a + dqbar_qn
  ye    = (1.0-phi_qn)*ye_a + dye_n
  !sumyi = (1.0-phi_fa)*yi_f + (phi_fa-phi_aq)*yi_a + dyi_qn

  ! calculate energy release including neutrino losses
  ! and rest mass differences
  qdot = cgsMeVperGram*(qbar - qbar_i)*dti  &
       - phi_qn*( yedot*N_A*c_l*c_l*(m_p + m_e - m_n) + edot )

  edotnu = phi_qn*edot
  
end subroutine bn_paraBurn
