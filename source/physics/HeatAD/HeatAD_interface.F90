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
#include "constants.h"

module HeatAD_interface

   implicit none

   interface
      subroutine HeatAD_init(restart)
         implicit none
         logical, intent(in) :: restart
      end subroutine HeatAD_init
   end interface

   interface
      subroutine HeatAD_finalize()
         implicit none
      end subroutine HeatAD_finalize
   end interface

   interface
      subroutine HeatAD_solve(solndata, del, lo, hi, dt)
        real, pointer, dimension(:, :, :, :) :: solnData
        real, INTENT(IN) :: dt
        real,dimension(MDIM),intent(IN) ::  del
        integer, dimension(MDIM), intent(IN) :: lo, hi
      end subroutine HeatAD_solve
   end interface

   interface
      subroutine HeatAD_advection(solnData, facexData, faceyData, facezData, del, lo, hi)
        real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
        real,dimension(MDIM),intent(IN) :: del
        integer, dimension(MDIM+1) :: face
        integer, dimension(MDIM),intent(IN) :: lo, hi
      end subroutine HeatAD_advection
   end interface

   interface
      subroutine HeatAD_diffusion(solnData, del, lo, hi)
        real, pointer, dimension(:, :, :, :) :: solnData
        real,dimension(MDIM),intent(IN) ::  del
        integer, dimension(MDIM),intent(IN) :: lo, hi
      end subroutine HeatAD_diffusion
   end interface

   interface
      subroutine HeatAD_indicators()
         implicit none
      end subroutine HeatAD_indicators
   end interface

   interface
      subroutine HeatAD_reInitGridVars(solnData)
        implicit none
        real, dimension(:,:,:,:), pointer :: solnData
      end subroutine HeatAD_reInitGridVars
   end interface

   interface
      subroutine HeatAD_getGridVar(name, value)
         implicit none
         character(len=*), intent(in)  :: name
         integer, intent(out)          :: value
      end subroutine HeatAD_getGridVar
   end interface

   interface
      subroutine HeatAD_getScalarProp(name, value)
         implicit none
         character(len=*), intent(in)  :: name
         real, intent(out)             :: value
      end subroutine HeatAD_getScalarProp
   end interface

end module HeatAD_interface
