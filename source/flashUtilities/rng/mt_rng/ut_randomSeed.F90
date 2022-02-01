!!****if* source/flashUtilities/rng/mt_rng/ut_randomSeed
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
!! NAME
!!
!!  ut_randomSeed
!!
!! SYNOPSIS
!!
!!  call ut_randomSeed(integer, optional, intent(OUT)  :: ut_size,
!!                     integer, optional, dimension(:), intent(IN)  :: ut_put,
!!                     integer, optional, dimension(:), intent(OUT)  :: ut_get)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ut_size : 
!!
!!   ut_put : 
!!
!!   ut_get : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

  subroutine ut_randomSeed(ut_size,ut_put, ut_get)

    implicit none

    integer, optional, dimension(:), intent(IN) :: ut_put
    integer, optional, dimension(:), intent(OUT) :: ut_get
    integer, optional, intent(OUT) :: ut_size
    logical :: noargs

    noargs = .not.(present(ut_size).or.present(ut_put).or.present(ut_get))
    if(noargs)then
       call ut_rand_init()
    else
       if(present(ut_size))then
          call ut_rand_init()
          call ut_get_size(ut_size)
       end if
       if(present(ut_get))call ut_get_seed(ut_get)
       if(present(ut_put))then
          print*,'calling put seed'
          call ut_put_seed(ut_put)
          print*,'done'
       end if
    end if
  end subroutine ut_randomSeed

