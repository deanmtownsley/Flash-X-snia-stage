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
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!***

#include "paramesh_preprocessor.fh"

Subroutine amr_1blk_cc_prol_dg               &
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  &
         mype,ivar)

  !-----Use Statements
  Use paramesh_dimensions, ONLY: nxb, nyb, nzb
  Use physicaldata, ONLY: unk1


  use RadTrans_interface, ONLY: RadTrans_prolongDgData

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

  Integer, Parameter :: largei = 100

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  !print*,'I am in prolongation dg!!!!'

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
  If (joff > 0) offj = nyb*K2D/2
  If (koff > 0) offk = nzb*K3D/2

  kcl = ((kfl-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K3D + 1 + offk
  kcu = ((kfu-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K3D + 1 + offk
  jcl = ((jfl-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K2D + 1 + offj
  jcu = ((jfu-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K2D + 1 + offj
  icl =  (ifl-NGUARD-1+largei)/2 + NGUARD - largei/2       + 1 + offi
  icu =  (ifu-NGUARD-1+largei)/2 + NGUARD - largei/2       + 1 + offi


  call RadTrans_prolongDgData(recv(ivar,icl:icu,jcl:jcu,kcl:kcu), &
                              unk1(ivar,ifl:ifu,jfl:jfu,kfl:kfu,idest))

End Subroutine amr_1blk_cc_prol_dg
