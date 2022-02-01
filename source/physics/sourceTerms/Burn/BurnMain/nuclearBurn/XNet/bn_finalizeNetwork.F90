!!****f* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/bn_finalizeNetwork
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
!!  bn_finalizeNetwork
!!
!!
!! SYNOPSIS
!! 
!!  call bn_finalizeNetwork()
!!
!!  
!! DESCRIPTION
!!
!!  Finalizes the bn_* submodule
!!
!!***


subroutine bn_finalizeNetwork()

  use bn_interface, ONLY : bn_xnetFinalize

  implicit none

  call bn_xnetFinalize

  return

end subroutine bn_finalizeNetwork
