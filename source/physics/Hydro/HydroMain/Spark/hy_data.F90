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

  !$omp target declare to &
  !$omp ( hy_cfl, &
  !$omp   hy_hydroComputeDtFirstCall, &
  !$omp   hy_updateHydroFluxes, &
  !$omp   hy_gcMaskSize, &
  !$omp   hy_restart, &
  !$omp   hy_useHydro, hy_telescoping, &
  !$omp   hy_meshMe, hy_globalComm, hy_meshComm, hy_meshNumProcs )

end Module hy_data
