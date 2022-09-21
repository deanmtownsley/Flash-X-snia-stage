!!****if* source/Grid/GridMain/AMR/Paramesh4/Thornado/amr_restrict_unk_dg
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
!!   amr_restrict_unk_dg
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_dg (datain, dataout, ivar)
!!   Call amr_restrict_unk_dg (real array, real array, integer, integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned
!!   Integer, Intent(in)    :: ivar  variable number in unk to restrict
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
!!   The last argument 'ivar' specifies which variable in 'unk' to apply
!!   the interpolation to.
!!
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!***

#include "paramesh_preprocessor.fh"

Subroutine amr_restrict_unk_dg(datain,dataout,ivar)

  !-----Use Statements
  Use paramesh_dimensions, ONLY: nxb, nyb, nzb

  use RadTrans_interface, ONLY: RadTrans_restrictDgData

  Implicit None

  !-----Input/Output arguments.
  Real,    Intent(in)    :: datain(:,:,:,:)
  Real,    Intent(inout) :: dataout(:,:,:,:)
  Integer, Intent(in)    :: ivar

  !-----Local arrays and variables.
  Integer :: ifl,ifu
  Integer :: jfl,jfu
  Integer :: kfl,kfu

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  !-----Being executable code.

  !print*,'I am in restriction dg!!!!'

  ifl = 1+NGUARD
  ifu = nxb+NGUARD
  jfl = 1+NGUARD*K2D
  jfu = nyb+NGUARD*K2D
  kfl = 1+NGUARD*K3D
  kfu = nzb+NGUARD*K3D

  call RadTrans_restrictDgData(datain (ivar,ifl:ifu,  jfl:jfu,  kfl:kfu), &
                               dataout(ivar,ifl:ifu:2,jfl:jfu:2,kfl:kfu:2))

End Subroutine amr_restrict_unk_dg
