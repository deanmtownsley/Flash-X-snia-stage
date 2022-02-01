! Helper function that indicates whether Grid metainformation coming
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
! from a checkpoint file, which has supposedly just just been used
! to initialize the following arrays owned by the Grid unit
! inplementation, is valid:
!
!    neigh
!    parent
!    child
!    surr_blks
logical function gr_pmIoTreeMetadataIsValid()
  use gr_specificData, ONLY : gr_gidIsValid
  use physicaldata,    ONLY : gsurrblks_set

  logical :: gsurrValid

  gsurrValid = (gsurrblks_set == +1)
  gr_pmIoTreeMetadataIsValid = gr_gidIsValid .AND. gsurrValid
end function gr_pmIoTreeMetadataIsValid
