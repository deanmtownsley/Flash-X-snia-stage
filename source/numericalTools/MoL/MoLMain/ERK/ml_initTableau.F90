!!****if* source/numericalTools/MoL/MoLMain/ERK/ml_initTableau
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
!!      ml_initTableau
!!
!!  SYNOPSIS
!!
!!      call ml_initTableau()
!!
!!  DESCRIPTION
!!
!!      Initialize the tableau(s) used by the selected MoL integrator
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine ml_initTableau()
   use ml_data, only: ml_A, ml_b, ml_c, ml_stages, ml_order, ml_method
   use ml_interface, only: ml_error

   use erk_tableau, only: euler_init, rk2_heun_init, rk3_ssp_init, rk4_init

   implicit none

   if (allocated(ml_A)) deallocate (ml_A)
   if (allocated(ml_b)) deallocate (ml_b)
   if (allocated(ml_c)) deallocate (ml_c)

   select case (ml_method)

   case ("erk-euler")
      call euler_init(ml_A, ml_b, ml_c, ml_order, ml_stages)

   case ("erk-rk2-heun")
      call rk2_heun_init(ml_A, ml_b, ml_c, ml_order, ml_stages)

   case ("erk-rk3-ssp")
      call rk3_ssp_init(ml_A, ml_b, ml_c, ml_order, ml_stages)

   case ("erk-rk4")
      call rk4_init(ml_A, ml_b, ml_c, ml_order, ml_stages)

   case default
      call ml_error("Unkown ERK method")
   end select

end subroutine ml_initTableau
