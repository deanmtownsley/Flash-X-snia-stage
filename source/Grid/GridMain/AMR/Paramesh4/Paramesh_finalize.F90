!!****if* source/Grid/GridMain/paramesh/Paramesh_finalize
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
!!  Paramesh_finalize
!! 
!! SYNOPSIS
!!  call Paramesh_finalize()
!!
!! DESCRIPTION
!!  Simple interface to finalize the mesh.
!!  Each mesh package will have a different copy of this function.
!!
!! USED BY
!!  Driver_finalizeAll via Grid_finalize
!!***
subroutine Paramesh_finalize()
  use paramesh_interfaces, ONLY : amr_close
  implicit none
  call amr_close()
end subroutine Paramesh_finalize

