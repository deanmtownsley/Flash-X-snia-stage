
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
Module gr_pmDivpres_mod

  public :: prol_fc_divpres
  logical, PARAMETER :: prol_fc_divpres = .false.
  ! No DIVPRES variables

Contains

  subroutine prol_fc_divpres_init(n,tf,i_divf_fc_vars)
    use Driver_interface, ONLY: Driver_abort
    implicit none
      
    integer, intent(in) :: n, tf, i_divf_fc_vars(tf,n)

    call Driver_abort(&
           "Attempting to initialize divpres module, not comfigured in!")

  end subroutine prol_fc_divpres_init


  subroutine gr_pmDivpresApply(recvfx,recvfy,recvfz, &
       ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
       lb)
    use paramesh_dimensions, ONLY: il_bnd1,jl_bnd1,kl_bnd1

    implicit none
    Real,intent(inout) :: recvfx(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    Real,intent(inout) :: recvfy(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    Real,intent(inout) :: recvfz(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    integer, intent(in) :: ia,ib,ja,jb,ka,kb
    integer, intent(in) :: idest,ioff,joff,koff
    integer, intent(in) :: lb


  end subroutine gr_pmDivpresApply

end Module gr_pmDivpres_mod
