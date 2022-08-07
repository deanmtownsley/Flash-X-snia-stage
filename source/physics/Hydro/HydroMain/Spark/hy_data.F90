Module hy_data
  implicit none
  save
  real :: hy_cfl
  logical :: hy_hydroComputeDtFirstCall
  logical :: hy_updateHydroFluxes
  integer :: hy_gcMaskSize
  logical :: hy_restart
  logical :: hy_useHydro, hy_telescoping
  integer :: hy_meshMe, hy_globalComm, hy_meshComm, hy_meshNumProcs
end Module hy_data
