

module Burn_data

  implicit none
#include "Flash.h"

  logical, save :: bn_useBurn
  logical, save :: bn_useShockBurn

  integer, save :: bn_meshMe
  real, save :: bn_enucDtFactor

  integer, dimension(NSPECIES), save :: bn_map_fi_to_mi ! array of mesa isotope index indexed by flash species index
  character (len=512), save :: bn_mesa_net_def_filename

  real, save :: bn_limiter_max_dlnT
end module Burn_data
