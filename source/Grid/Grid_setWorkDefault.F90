!!****f* source/Grid/Grid_setWorkDefault
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
!!  Grid_setWorkDefault
!!
!! SYNOPSIS
!!
!!  Grid_setWorkDefault(real(IN)   :: pwork,
!!                      real(IN)   :: lwork)
!!
!! DESCRIPTION
!!  Sets the default work value. These value is assigned to every
!!  block after regridding, unless gr_btExchangeWork = True.
!!  Used for uneven load distributions when using Paramesh.
!!
!! ARGUMENTS
!!  pwork - default work for non-leaf blocks
!!  lwork - default work for leaf blocks
!!
!!***

      subroutine Grid_setWorkDefault(pwork,lwork)

      implicit none
      real,intent(in)    :: pwork
      real,intent(in)    :: lwork

      !implementation in Paramesh

      return
      end subroutine Grid_setWorkDefault
