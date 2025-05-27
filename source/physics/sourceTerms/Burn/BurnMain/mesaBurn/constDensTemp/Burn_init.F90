

subroutine Burn_init()

  use Burn_data, ONLY: bn_useBurn, bn_useShockBurn, bn_meshMe, bn_enucDtFactor, bn_map_fi_to_mi, &
                       bn_mesa_net_def_filename, bn_limiter_max_dlnT
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Driver_interface, ONLY: Driver_getMype, Driver_abort
  use Multispecies_interface, ONLY: Multispecies_getProperty
  use mesa_flash_nucreact, ONLY: init_mesa_internals, fill_nuclide_index_map

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Multispecies.h"

  integer :: i, ierr
  integer,dimension(NSPECIES) :: zlist, alist
  

  call Driver_getMype(MESH_COMM, bn_meshMe)
  call RuntimeParameters_get("useBurn", bn_useBurn )
  call RuntimeParameters_get('useShockBurn', bn_useShockBurn)
  call RuntimeParameters_get("enucDtFactor", bn_enucDtFactor)
  call RunTimeParameters_get("bn_mesa_net_def_filename", bn_mesa_net_def_filename)

  call RuntimeParameters_get("bn_limiter_max_dlnT", bn_limiter_max_dlnT)

  call init_mesa_internals(trim(bn_mesa_net_def_filename))

  ! get Z and A of species included in Flash order (needed to construct maps to mesa species)
  do i=1,NSPECIES
    call Multispecies_getProperty( SPECIES_BEGIN+i-1, Z, zlist(i) )
    call Multispecies_getProperty( SPECIES_BEGIN+i-1, A, alist(i) )
  end do
  call fill_nuclide_index_map( NSPECIES, zlist, alist, bn_map_fi_to_mi, ierr)
  if (ierr == 1) call Driver_abort("Unable to create map from flash species to MESA nuclides")

  return

end subroutine Burn_init
