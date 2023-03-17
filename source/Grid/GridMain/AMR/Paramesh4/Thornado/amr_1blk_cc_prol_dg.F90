!!****if* source/Grid/GridMain/AMR/Paramesh4/Thornado/amr_1blk_cc_prol_dg
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!! DESCRIPTION
!!
!!   This custom prolongation routine is called from
!!   amr_1blk_cc_prol_gen_unk_fun, when requested by the value of
!!   interp_mask_unk(ivar).
!!
!!   It is applied to all UNK variables whose corresponding element
!!   of interp_mask is set to 40.
!!
!!   The amr_1blk_cc_prol_dg routine accepts a "coarse" block (recv) as
!!   input and returns a "fine" block (unk1).
!!   The routine will prepare an array from recv and then call a routine that 
!!   is in Thornado- the Refine_TwoMoment.
!!   It will subsequently transfer information as an output from the Thornado
!!   routine into the "fine" block array.
!!
!! NOTES
!!
!!   In the case of curvilinear coordinates, the arrays cell_face_coord1, cell_face_coord2,
!!   cell_face_coord3 from the physicaldata module are assumed to contain the coordinates
!!   of cell faces for the full coarse block, i.e., for the data in the recv argument.
!!   The caller should therefore call amr_block_geometry (on the coarse block)
!!   before this routine is called.
!!
!! CALLS
!!
!!   RadTrans_prolongDgData
!!   No calls made to other PARAMESH routines.
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit the data in recv is interpolated
!!   and placed into the unk1 array.
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2022-09-22 Added computation of skip            - Klaus Weide
!!  2023-03-14 Pass cell_face_coordN when needed    - Klaus Weide
!!***

#include "paramesh_preprocessor.fh"
#include "constants.h"

Subroutine amr_1blk_cc_prol_dg               &
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  &
         mype,ivar)

  !-----Use Statements
  Use paramesh_dimensions, ONLY: nxb, nyb, nzb
  Use physicaldata, ONLY: unk1
  Use physicaldata, ONLY: cell_face_coord1, cell_face_coord2, cell_face_coord3
  Use physicaldata, ONLY: curvilinear

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
  Integer :: iclX, icuX, jclX, jcuX, kclX, kcuX ! extended range bounds
  Integer :: offi, offj, offk

  integer :: skip(MDIM)

  Integer, Parameter :: largei = 200

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

  ! The following variables are now unused; left for code comparision:
  kcl = ((kfl-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K3D + 1 + offk
  kcu = ((kfu-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K3D + 1 + offk
  jcl = ((jfl-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K2D + 1 + offj
  jcu = ((jfu-NGUARD-1+largei)/2 + NGUARD - largei/2 )*K2D + 1 + offj
  icl =  (ifl-NGUARD-1+largei)/2 + NGUARD - largei/2       + 1 + offi
  icu =  (ifu-NGUARD-1+largei)/2 + NGUARD - largei/2       + 1 + offi

!!$  !Example numbers are for (NXB,NYB,NZB)=(16,16,16), NGUARD=6,
!!$  !prolongation into a new child block (interior cells):
!!$  !  [ ia, ja, ka ] == [ 7, 7, 7 ]
!!$  !  [ ib, jb, kb ] == [ 22, 22, 22 ]
!!$  fine_lo = [ ifl, jfl, kfl ] ! [ 7, 7, 7 ]
!!$  fine_hi = [ ifu, jfu, kfu ] ! [ 22, 22, 22 ]
!!$
!!$  crse_lo = [ icl, jcl, kcl ] ! [ 7, 7, 7 ] ; [ 15, 15, 15 ]
!!$  crse_hi = [ icu, jcu, kcu ] ! [ 14, 14, 14 ] ; [ 22, 22 22 ]
!!$  !(Everything is fine here!)
!!$
!!$  !Example numbers are for (NXB,NYB,NZB)=(16,16,16), NGUARD=6,
!!$  !f/c interpolation at a left ; right (all-layers) gc region.
!!$  !  [ ia, ja, ka ] == [ 1, 7, 7 ]   ; [ 23, 7, 7 ]
!!$  !  [ ib, jb, kb ] == [ 6, 22, 22 ] ; [ 28, 22, 22 ]
!!$  fine_lo = [ ifl, jfl, kfl ] ! [ 1, 7, 7 ]   ; [ 23, 7, 7 ]
!!$  fine_hi = [ ifu, jfu, kfu ] ! [ 6, 22, 22 ] ; [ 28, 22, 22 ]
!!$
!!$  crse_lo = [ icl, jcl, kcl ] ! [ 4, 7, 7 ]   ; [ 23, 7, 15 ]
!!$  crse_hi = [ icu, jcu, kcu ] ! [ 6, 14, 14 ] ; [ 25, 14, 22 ]
!!$
!!$  !needs to be turned into:
!!$
!!$  crse_loX = [ iclX, jclX, kclX ] ! [ 3(!), 7, 7 ]   ; [ 23, 7, 15 ]
!!$  crse_hiX = [ icuX, jcuX, kcuX ] ! [ 6, 14, 14 ] ; [ 26(!), 14, 22 ]

  kclX = ( ((kfl-NGUARD-1+4*largei)/4 - largei)*2 + NGUARD  )*K3D + 1 + offk
  kcuX = ( ((kfu-NGUARD-1+4*largei)/4 - largei)*2 + NGUARD+1)*K3D + 1 + offk
  jclX = ( ((jfl-NGUARD-1+4*largei)/4 - largei)*2 + NGUARD  )*K2D + 1 + offj
  jcuX = ( ((jfu-NGUARD-1+4*largei)/4 - largei)*2 + NGUARD+1)*K2D + 1 + offj
  iclX =   ((ifl-NGUARD-1+4*largei)/4 - largei)*2 + NGUARD        + 1 + offi
  icuX =   ((ifu-NGUARD-1+4*largei)/4 - largei)*2 + NGUARD        + 2 + offi

  skip = (/ &
            modulo(ia-NGUARD-1, 4)    , &
            modulo(ja-NGUARD-1, 4)*K2D, &
            modulo(ka-NGUARD-1, 4)*K3D  &
         /)

  if (.NOT. curvilinear) then
     call RadTrans_prolongDgData(recv(ivar,iclX:icuX,jclX:jcuX,kclX:kcuX), &
                              unk1(ivar,ifl :ifu, jfl :jfu, kfl :kfu, idest), &
                              skip)
  else
     call RadTrans_prolongDgData(recv(ivar,iclX:icuX,jclX:jcuX,kclX:kcuX), &
                              unk1(ivar,ifl :ifu, jfl :jfu, kfl :kfu, idest), &
                              skip,                                        &
                              cell_face_coord1(iclX:icuX+1),               &
                              cell_face_coord2(jclX:jcuX+1),               &
                              cell_face_coord3(kclX:kcuX+1))
  end if

End Subroutine amr_1blk_cc_prol_dg
