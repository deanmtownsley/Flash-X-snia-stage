!!****if* source/Grid/GridMain/AMR/Paramesh4/Thornado/amr_1blk_cc_prol_dg
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
!!   amr_1blk_cc_prol_dg
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_dg (recv,ia,ib,ja,jb,ka,kb,idest,
!!                                 ioff,joff,koff,mype,ivar)
!!   Call amr_1blk_cc_prol_dg (real,
!!                                 integer, integer, integer, integer,
!!                                 integer, integer, integer, integer,
!!                                 integer, integer, integer, integer)
!!
!! ARGUMENTS
!!***Desriptions need to be updated***
!!  Real,    intent(inout) :: recv(:,:,:,:)
!!    Data array holding the data extracted from unk which will be prolonged
!!    and placed into the unk1 array.
!!
!!  Integer, Intent(in) :: ia,ib,ja,jb,ka,kb
!!    Integers which control the limits into unk1 where the prolonged data
!!    will be placed.
!!
!!  Integer, Intent(in) :: idest,ioff,joff,koff,mype
!!    idest controls which 'layer' into which the prolonged data will be
!!    placed in unk1.  ioff, joff and koff are offsets.  mype is is the
!!    local processor id.
!!
!!  Integer, intent(in) :: ivar,
!!    ivar is the varible number in unk which is prolonged.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   This routine takes data from the array recv, originally extracted
!!   from the solution array unk, and performs a prolongation operation
!!   on it, between the bounds ranges ia to ib, ja to jb, and ka to kb.
!!   The data in recv is from a parent block and the
!!   result of the prolongation operation is written directly into one
!!   layer of the working block array unk1(...,idest).
!!   The position of the child within the parent block is specified by
!!   the ioff, joff and koff arguments.
!!
!!   Created a custom prolongation routine that is called from
!!   amr_1blk_cc_prol_gen_unk_fun depending on the value of
!!   interp_mask_unk(ivar). The method is named amr_1blk_cc_prol_dg.
!!
!!   It is applied to all UNK variables whose corresponding element
!!   of interp_mask is set to 40.
!!
!!   The amr_1blk_cc_prol_dg routine accepts a "coarse" block (recv) as
!!   input and return a "fine" block (unk1).
!!   The routine will prepare an array from recv and then call a routine that 
!!   is in Thornado- the Refine_TwoMoment.
!!   It will subsequently transfer information as an output from the Thornado
!!   routine into the "fine" block array.
!!
!!
!! CALLS
!!
!!   No calls made to other PARAMESH routines.
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit the data in recv is interpolated
!!   and placed into the unk1 array.
!!
!! DESCRIPTION
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!!***

#include "paramesh_preprocessor.fh"

Subroutine amr_1blk_cc_prol_dg               &
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  &
         mype,ivar)

  !-----Use Statements
  Use paramesh_dimensions
  Use physicaldata
  Use tree
  Use prolong_arrays

  Use TwoMoment_MeshRefinementModule, Only : &
     RefineX_TwoMoment

  Implicit None

  !-----Input/Output Variables
  Real,    Intent(inout) :: recv(:,:,:,:)
  Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb
  Integer, Intent(in)    :: idest,ioff,joff,koff,mype
  Integer, Intent(in)    :: ivar

  !-----Local variables
  Integer :: ifl, ifu, jfl, jfu, kfl, kfu
  Integer :: icl, icu, jcl, jcu, kcl, kcu
  Integer :: offi, offj, offk
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: ic, jc, kc, ii, jj, kk, icc, jcc, kcc

  Integer :: iNodeX
  Integer :: crse_lo(3), crse_hi(3)
  Integer :: fine_lo(3), fine_hi(3)
  Integer :: iFineX, nFineX(3), nX(3)
  Integer :: u_lo(3), u_hi(3)

  Integer, Parameter :: largei = 100

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  Real    :: U_Crse(THORNADO_FLUID_NDOF)
  Real    :: U_Fine(THORNADO_FLUID_NDOF)

  !print*,'I am in prolongation dg!!!!'

  nFineX = 1
  nFineX(1:NDIM) = refine_factor

  ifl=ia
  ifu=ib
  jfl=ja
  jfu=jb
  kfl=ka
  kfu=kb

  offi = 0
  offj = 0
  offk = 0
  If (ioff > 0) offi = nxb/2
  If (joff > 0) offj = nyb*k2d/2
  If (koff > 0) offk = nzb*k3d/2

  kcl = ((kfl-nguard-1+largei)/2 + nguard - largei/2 )*k3d + 1 + offk
  kcu = ((kfu-nguard-1+largei)/2 + nguard - largei/2 )*k3d + 1 + offk
  jcl = ((jfl-nguard-1+largei)/2 + nguard - largei/2 )*k2d + 1 + offj
  jcu = ((jfu-nguard-1+largei)/2 + nguard - largei/2 )*k2d + 1 + offj
  icl =  (ifl-nguard-1+largei)/2 + nguard - largei/2       + 1 + offi
  icu =  (ifu-nguard-1+largei)/2 + nguard - largei/2       + 1 + offi

  fine_lo = [ ifl, jfl, kfl ] ! [ 5, 5, 5 ]
  fine_hi = [ ifu, jfu, kfu ] ! [ 20, 20, 20 ]

  crse_lo = [ icl, jcl, kcl ] ! [ 5, 5, 5 ] ; [ 13, 13, 13 ]
  crse_hi = [ icu, jcu, kcu ] ! [ 12, 12, 12 ] ; [ 20, 20 20 ]

  nX = 1
  nX(1:NDIM) = ( crse_hi(1:NDIM) - crse_lo(1:NDIM) + 1 ) / THORNADO_NNODESX ! [ 4, 4, 4 ]

  u_lo = 1
  u_hi = nX

  !-----Interpolation loop.

  ! loop over coarse (thornado) elements in parent block
  do kc = u_lo(3), u_hi(3)
     do jc = u_lo(2), u_hi(2)
        do ic = u_lo(1), u_hi(1) ! 1, 2, 3, 4

           ! get unk indices for parent block
           k1 = crse_lo(3) + THORNADO_NNODESX*(kc-1)*k3d
           j1 = crse_lo(2) + THORNADO_NNODESX*(jc-1)*k2d
           i1 = crse_lo(1) + THORNADO_NNODESX*(ic-1)    ! 5, 7, 9, 11
           ! 13, 15, 17, 19

           ! grab the data from coarse grid element quadrature points
           do iNodeX = 1, THORNADO_FLUID_NDOF
              kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
              jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
              ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1
              U_Crse(iNodeX) = recv(ivar,ii,jj,kk)
           end do

           ! unk offsets of first child element in child block
           k0 = 1 + refine_factor*THORNADO_NNODESX*(kc-1)*k3d
           j0 = 1 + refine_factor*THORNADO_NNODESX*(jc-1)*k2d
           i0 = 1 + refine_factor*THORNADO_NNODESX*(ic-1)     ! 1, 5, 9, 13

           ! loop over fine grid element for this parent element
           iFineX = 0
           do kcc = 1, nFineX(3)
              do jcc = 1, nFineX(2)
                 do icc = 1, nFineX(1)
                    iFineX = iFineX + 1

                    ! compute fine grid element quadrature points
                    U_Fine = RefineX_TwoMoment( iFineX, U_Crse )

                    ! calculate unk indices in child block
                    k = fine_lo(3) + (k0-1) + THORNADO_NNODESX*(kcc-1)
                    j = fine_lo(2) + (j0-1) + THORNADO_NNODESX*(jcc-1)
                    i = fine_lo(1) + (i0-1) + THORNADO_NNODESX*(icc-1) ! ic = 1 ; i0 = 1  -> 5, 7
                    ! ic = 2 ; i0 = 5  -> 9, 11
                    ! ic = 3 ; i0 = 9  -> 13, 15
                    ! ic = 4 ; i0 = 13 -> 17, 19

                    ! store the result in child block
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                       jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                       ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i
                       unk1(ivar,ii,jj,kk,idest) = U_Fine(iNodeX)
                    end do

                 end do
              end do
           end do

        end do
     end do
  end do

  Return
End Subroutine amr_1blk_cc_prol_dg
