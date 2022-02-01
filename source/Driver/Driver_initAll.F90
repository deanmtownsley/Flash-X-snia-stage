!!****f* source/Driver/Driver_initAll
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
!!  Driver_initAll
!!
!! SYNOPSIS
!!
!!   Driver_initAll()
!!
!! DESCRIPTION
!!
!!  Performs Flash initializations, which includes:
!!
!!  Call all 'init' routines in units, in the right order.
!!  Order does matter, particularly when restarting from a 
!!  checkpoint file.
!!
!!  For the most part, Driver_initAll calls another units' init
!!  routines directly, like call IO_init or call Grid_init.  This
!!  routine also makes calls to other Driver initialization routines
!!  like Driver_initMaterialProperties or Driver_initSourceTerms.
!!  These routines then call the unit-specific initialization 
!!  routines.  This level of abstraction was added to simplify
!!  the initialization calls.
!!
!!
!!
!!***


!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"


subroutine Driver_initAll()

  
implicit none
  return
end subroutine Driver_initAll
