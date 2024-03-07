Module Orchestration_interfaceTypeDecl
   use,intrinsic :: iso_c_binding, ONLY : c_ptr, C_INT, C_DOUBLE
   use milhoja_types_mod,          ONLY : MILHOJA_INT
   use milhoja_runtime_mod,        ONLY : milhoja_runtime_taskFunction

   implicit none
#include "constants.h"

#define TYPENAME_PREFIXED(surname) TYPENAME_PREFIXED2(Orchestration,surname)
#include "Milhoja_tileCInfo.finc"

end Module Orchestration_interfaceTypeDecl
! Local Variables:
! f90-program-indent: 3
! indent-tabs-mode: nil
! End:
