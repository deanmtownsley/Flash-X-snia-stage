!!****if* source/Grid/localAPI/gr_mpolePotentials
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
!!  gr_mpolePotentials
!!
!! SYNOPSIS
!!
!!  gr_mpolePotentials  (integer, intent(in) :: ipotvar,
!!                       real,    intent(in) :: Poisson_factor )
!!
!! DESCRIPTION
!!
!!  Computes the potential field using the mass moments already
!!  calculated. On output tha variable indexed by ipotvar contains
!!  the potential. The calculations are entirely local to each
!!  processor, since each processor has a local copy of the moments.
!!
!!  This routine calls the appropriate subroutines according to
!!  the geometry specified.
!!
!! ARGUMENTS
!!
!!  ipotvar        : index to variable containing the potential
!!  Poisson_factor : the name says it all 
!!
!!***

subroutine gr_mpolePotentials (ipotvar,Poisson_factor)

  implicit none
    
  integer, intent (in) :: ipotvar
  real,    intent (in) :: Poisson_factor

  return
end subroutine gr_mpolePotentials
