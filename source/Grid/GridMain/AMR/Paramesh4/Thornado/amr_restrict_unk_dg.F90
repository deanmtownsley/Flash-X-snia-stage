!!****if* source/Grid/GridMain/AMR/Paramesh4/Thornado/amr_restrict_unk_dg
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
!!   amr_restrict_unk_dg
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_dg (datain, dataout, ioff, joff, koff)
!!   Call amr_restrict_unk_dg (real array, real array, integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned
!!   ioff, joff, koff : offsets of the restricted data that is returned within the coarse
!!                      block; each of these numbers should be either 0 or half
!!                      the number of interior cells in the block for the
!!                      relevant direction.
!!
!! INCLUDES
!! 
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! CALLS
!! 
!! RETURNS
!!
!!   Restricted data returned in array 'dataout'.  
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   DG quadrature-point data stored in 'unk'.
!!
!!   Data is passed in in the array 'datain' and returned in the array
!!   'dataout'.
!!   The restriction operates on a subset of variables in UNK selected
!!   by corresponding entries in the array interp_mask_unk_res.
!!
!! NOTES
!!
!!   In the case of curvilinear coordinates, the arrays cell_face_coord1, cell_face_coord2,
!!   cell_face_coord3 from the physicaldata module are assumed to contain the coordinates
!!   of cell faces for the full coarse block, i.e., for the data in the recv argument.
!!   The caller should therefore call amr_block_geometry (on the coarse block)
!!   before this routine is called. This is the case when the call ultimately comes
!!   from the routine mpi_amr_1blk_restrict and the 'curvilinear' flag is set.
!!
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2023-03-15 Pass cell_face_coordN when needed    - Klaus Weide
!!  2023-03-16 Accept ioff,joff,koff args           - Klaus Weide
!!  2023-04-05 Act on a set of variables, not all   - Austin Harris
!!  2023-04-07 Eliminate ivar argument              - Klaus Weide
!!***

#include "paramesh_preprocessor.fh"

Subroutine amr_restrict_unk_dg(datain,dataout,ioff,joff,koff)

  !-----Use Statements
  Use paramesh_dimensions, ONLY: nxb, nyb, nzb
  Use physicaldata, ONLY: cell_face_coord1, cell_face_coord2, cell_face_coord3
  Use physicaldata, ONLY: curvilinear
  Use physicaldata, ONLY: interp_mask_unk_res

  use RadTrans_interface, ONLY: RadTrans_restrictDgData

  Implicit None

  !-----Input/Output arguments.
  Real,    Intent(in)    :: datain(:,:,:,:)
  Real,    Intent(inout) :: dataout(:,:,:,:)
  Integer, Intent(in)    :: ioff,joff,koff

  !-----Local arrays and variables.
  Integer :: ifl,ifu
  Integer :: jfl,jfu
  Integer :: kfl,kfu
  Integer :: iclX, icuX, jclX, jcuX, kclX, kcuX ! coarse block bounds

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  !-----Being executable code.

  !print*,'I am in restriction dg!!!!'

  ifl = 1+NGUARD
  ifu = nxb+NGUARD
  jfl = 1+NGUARD*K2D
  jfu = nyb+NGUARD*K2D
  kfl = 1+NGUARD*K3D
  kfu = nzb+NGUARD*K3D

  if (.NOT. curvilinear) then
     call RadTrans_restrictDgData(datain (:,ifl:ifu,  jfl:jfu,  kfl:kfu), &
                                  dataout(:,ifl:ifu:2,jfl:jfu:2,kfl:kfu:2), &
                                  (interp_mask_unk_res == 40) )
  else
     ! Compute corresponding index bounds in the coarse block, into which the
     ! data we store in dataout will ultimately be copied.
     iclX = 1 +           ioff + NGUARD
     icuX = 1 +  nxb/2  + ioff + NGUARD
     jclX = 1 + (       + joff + NGUARD) * K2D
     jcuX = 1 + (nyb/2  + joff + NGUARD) * K2D
     kclX = 1 + (       + koff + NGUARD) * K3D
     kcuX = 1 + (nzb/2  + koff + NGUARD) * K3D
     call RadTrans_restrictDgData(datain (:,ifl:ifu,  jfl:jfu,  kfl:kfu),   &
                                  dataout(:,ifl:ifu:2,jfl:jfu:2,kfl:kfu:2), &
                                  (interp_mask_unk_res == 40), &
                                  cell_face_coord1(iclX:icuX+1),               &
                                  cell_face_coord2(jclX:jcuX+K2D),             &
                                  cell_face_coord3(kclX:kcuX+K3D))
  end if
End Subroutine amr_restrict_unk_dg
