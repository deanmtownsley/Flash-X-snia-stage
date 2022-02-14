!!****f* source/Grid/Grid_conserveField
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
!!  Grid_conserveField
!!
!! SYNOPSIS
!!
!!  call Grid_conserveField ()
!!
!! ARGUMENTS
!!  none
!!
!! DESCRIPTION
!!
!!  Corrects electric fields at refinement jump boundaries to make 
!!  sure electric fields at the common boundaries are consistent 
!!  at both refined and coarse meshes.
!!
!! NOTES
!!
!!  DEV: Currently only implemented for Paramesh4!
!!
!!  On the uniform grid this routine is not needed.
!!
!!***



subroutine Grid_conserveField ()


  implicit none

  return

end subroutine Grid_conserveField
