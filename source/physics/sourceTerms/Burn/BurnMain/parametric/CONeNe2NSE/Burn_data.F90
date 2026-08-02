!
! Dean Townsley 2008
!

module Burn_data

#include "constants.h"

  implicit none

  integer, save :: bn_meshMe, bn_meshComm, bn_meshNumProcs

  logical :: bn_useBurn
  logical :: bn_useShockBurn
  real    :: bn_enucDtFactor

  ! Logical file unit for saving detonation ignition points
  integer, save :: burn_lun = 31
  character(len=MAX_STRING_LENGTH), save :: bn_detIgnFileName 
  
  ! I don't strictly follow the flash convention by adding bn_ to everything
  ! because by having these in a module the namespace issue that resolves is
  ! already taken care of, and the convention is clumsy for things like these
  real, save         :: ye_he4, yi_he4, q_he4
  real, save         :: ye_c12, yi_c12, q_c12
  real, save         :: ye_o16, yi_o16, q_o16
!don: add ne20
  real, save         :: ye_ne20, yi_ne20, q_ne20
  real, save         :: ye_ne22, yi_ne22, q_ne22
  real, save         :: ye_mg24, yi_mg24, q_mg24
  real, save         :: ye_si28, yi_si28, q_si28

  real, save :: N_A, m_p, m_n, m_e, c_l
  
  logical, save      :: bn_thermalReact, bn_autoDDT
  real, save ::  bn_thermalReactInFlameThreshold

  ! table for detonation ignition points
  real, save :: pbIgnRho, pbIgnRhoFact, pbIgnPhfa, pbIgnDist, pbIgnRad, pbIgnSep
  integer, save :: pbIgnNum, pbIgnNumMax
  real, allocatable, dimension(:) :: pbIgnTime, pbIgnX, pbIgnY, pbIgnZ, pbIgnR

  ! variable for saving processor-local and global neutrino loss energy integrals
  real, save :: bn_neutLossThisProcStep, bn_neutLoss

end module
