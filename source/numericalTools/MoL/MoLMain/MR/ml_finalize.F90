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
   use mr_data, only: FE, FI, FF, mr_gamK, mr_wK, &
                      mr_gamBar, mr_wBar, mr_cS, mr_AF, mr_bF, mr_cF

   implicit none

   if (allocated(FE)) deallocate (FE)
   if (allocated(FI)) deallocate (FI)
   if (allocated(FF)) deallocate (FF)

   if (allocated(mr_gamK)) deallocate (mr_gamK)
   if (allocated(mr_wK)) deallocate (mr_wK)

   if (allocated(mr_gamBar)) deallocate (mr_gamBar)
   if (allocated(mr_wBar)) deallocate (mr_wBar)

   if (allocated(mr_cS)) deallocate (mr_cS)

   if (allocated(mr_AF)) deallocate (mr_AF)
   if (allocated(mr_bF)) deallocate (mr_bF)
   if (allocated(mr_cF)) deallocate (mr_cF)
end subroutine ml_finalize
