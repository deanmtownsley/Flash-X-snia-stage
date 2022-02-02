!!****if* source/IO/IOMain/io_bcastScalars
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
!!  io_bcastScalars
!!
!! SYNOPSIS
!!
!!  io_bcastScalars()
!!
!! DESCRIPTION
!!
!! broadcasts scalars from global me to the other processors
!! calls nameValueLL_bcast to do the implementation
!! 
!! 
!! 
!!
!! ARGUMENTS
!!
!!        
!!
!!
!!
!!***

subroutine io_bcastScalars()

  use IO_data, only : io_scalar, io_globalMe
  
implicit none

  call nameValueLL_bcast(io_scalar, io_globalMe)

  return

end subroutine io_bcastScalars


  

