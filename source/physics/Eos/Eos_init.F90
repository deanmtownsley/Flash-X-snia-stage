!!****f* source/physics/Eos/Eos_init
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
!!  Eos_init
!!
!! 
!! SYNOPSIS
!!
!!  call Eos_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars needed
!!  by the EOS unit from the runtime parameters and physical
!!  constants facilities
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!  
!!   Particular implementations (Gamma,Helmholz,etc) of the unit
!!   define their own runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might overwrite these values with the 
!!   flash.par values for your specific run.  
!!
!!
!!***


subroutine Eos_init()
    
    implicit none
    ! stub for eos initialization.  A non-stub implementation of his routine
    ! will be supplied, if required,  by the main subunit of Eos, normally
    ! located under Eos/EosMain.

    return
end subroutine Eos_init


