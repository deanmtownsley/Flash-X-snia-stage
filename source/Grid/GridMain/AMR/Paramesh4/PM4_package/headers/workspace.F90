!----------------------------------------------------------------------
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! workspace module


      Module workspace

#include "paramesh_preprocessor.fh"

      Private

! workspace arrays
      Public :: work, recvw, sendw, tempw, interp_mask_work, & 
     &          interp_mask_work_res
      Real, Allocatable, Save :: work(:,:,:,:,:)
      Real, Allocatable, Save :: recvw(:,:,:)
      Real, Allocatable, Save :: sendw(:,:,:)
      Real, Allocatable, Save :: tempw(:,:,:)
      Integer, Allocatable, Save :: interp_mask_work(:)
      Integer, Allocatable, Save :: interp_mask_work_res(:)

! common block storing the solution for cell-centered quantities.
      Public :: work1, recvw1, tempw1
      Real, Allocatable, Save :: work1(:,:,:,:)
      Real, Allocatable, Save :: recvw1(:,:,:,:)
      Real, Allocatable, Save :: tempw1(:,:,:)

! arrays used to store geometry information for the working block
      Public :: cell_vol_w
      Real, Allocatable :: cell_vol_w(:,:,:)

! Index arrays used to record destination data values for fine layer
! neighbor guardcells
      Integer,Public :: f2c_ind_work(2,3,27)

      End Module workspace


