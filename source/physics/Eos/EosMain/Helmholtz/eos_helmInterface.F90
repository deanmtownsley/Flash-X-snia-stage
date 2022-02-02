!!****ih* source/physics/Eos/EosMain/Helmholtz/eos_helmInterface
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
!!     eos_helmInterface
!!
!! SYNOPSIS
!!     use eos_helmInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! the Helmholtz Eos implementation.
!!
!!***

module eos_helmInterface
  implicit none
#include "Eos.h"

  interface
     subroutine eos_helm(eos_jlo,eos_jhi,mask)
       integer, intent(IN) :: eos_jlo, eos_jhi
       logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
     end subroutine eos_helm
  end interface
end module eos_helmInterface
