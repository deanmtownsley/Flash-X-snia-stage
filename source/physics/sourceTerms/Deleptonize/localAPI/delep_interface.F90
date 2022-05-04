!!****ih* source/physics/sourceTerms/Deleptonize/localAPI/delep_interface
!!
!!   Deleptonize's private interfaces
!!***

Module delep_interface
#include "constants.h"

  interface
     subroutine delep_detectBounce (dotiling,dt,time)
       real,intent(IN) :: dt,time
       logical, intent(IN) :: dotiling
     end subroutine delep_detectBounce
  end interface

end Module delep_interface
