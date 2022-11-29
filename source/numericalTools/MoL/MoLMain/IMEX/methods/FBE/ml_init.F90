!!****if* source/numericalTools/MoL/MoLMain/FBE/ml_init
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!  NAME
!!
!!      ml_init
!!
!!  SYNOPSIS
!!
!!      call ml_init()
!!
!!  DESCRIPTION
!!
!!      Initialize a method of lines unit implementation
!!
!!  ARGUMENTS
!!
!!
!!***
!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief ml_init implementation for FBE

!> @ingroup MoLFBE
!!
!! @brief Implements ml_init for FBE
!!
!! @stubref{ml_init}
subroutine ml_init()
   use MoL_data, only: ml_nscratch

   implicit none

   ! Just need MOL_RHS
   ml_nscratch = 0
end subroutine ml_init
