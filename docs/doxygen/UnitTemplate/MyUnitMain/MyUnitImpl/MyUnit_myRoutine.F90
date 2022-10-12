!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! @file

The file documentation block should be at the top of the file for legal reasons
and should not need modification.  C-preprocessor macros can go here between
documentation blocks if so desired.  The next documentation block should be
placed as close to the subroutine/function as possible.  The general goal is to
make that content as minimal as possible.  For example, we don't include @param
here since those are documented in the stubs.  The documentation in @details
should contain only implementation-specific details that are useful for
developers/maintainers.

@ingroup should refer to the group name of this file's concrete unit
implementation.  This is declared as the first item in the @defgroup statement
in the .dox in the same folder.

The argument to @stubref must always be the routine's name.

Aside from updating the name of the subroutine, no other changes should be made
to the @brief line.

#include "Simulation.h"

!> @ingroup MyUnitImpl
!! @stubref{MyUnit_myRoutine}
!!
!! @brief Concrete implementation of MyUnit_myRoutine.
!!
!! @details
!! <Add in detailed implementation-specific docs here if any.>
subroutine MyUnit_myRoutine(a, b, c)
...
end subroutine MyUnit_myRoutine
