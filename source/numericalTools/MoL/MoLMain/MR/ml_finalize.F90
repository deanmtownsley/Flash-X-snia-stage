!!****if* source/numericalTools/MoL/MoLMain/MR/ml_finalize
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
!!      ml_finalize
!!
!!  SYNOPSIS
!!
!!      call ml_finalize()
!!
!!  DESCRIPTION
!!
!!      Finalize a method of lines unit implementation
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine ml_finalize()
    use mr_data

    implicit none

    if (allocated(FE)) deallocate(FE)
    if (allocated(FI)) deallocate(FI)
    if (allocated(FF)) deallocate(FF)

    if (allocated(gamK)) deallocate(gamK)
    if (allocated(wK))   deallocate(wK)

    if (allocated(gamBar)) deallocate(gamBar)
    if (allocated(wBar))   deallocate(wBar)

    if (allocated(cS)) deallocate(cS)

    if (allocated(AF)) deallocate(AF)
    if (allocated(bF)) deallocate(bF)
    if (allocated(cF)) deallocate(cF)
end subroutine ml_finalize
