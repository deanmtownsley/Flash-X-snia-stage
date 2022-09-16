!!****if* source/Grid/GridMain/AMR/Paramesh4/Thornado/amr_restrict_unk_dg
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
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
!!***

#include "paramesh_preprocessor.fh"

Subroutine amr_restrict_unk_dg(datain,dataout,ivar)

  !-----Use Statements
  Use paramesh_dimensions
  Use physicaldata

  Use TwoMoment_MeshRefinementModule, Only : &
    ProjectionMatrix_T, VolumeRatio
  Use ReferenceElementModuleX, Only : &
    WeightsX_q

  Implicit None

  !-----Input/Output arguments.
  Real,    Intent(in)    :: datain(:,:,:,:)
  Real,    Intent(inout) :: dataout(:,:,:,:)
  Integer, Intent(in)    :: ivar

  !-----Local arrays and variables.
  Integer :: ifl,ifu
  Integer :: jfl,jfu
  Integer :: kfl,kfu
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: ic, jc, kc, icc, jcc, kcc, ii, jj, kk

  Integer :: iNodeX
  Integer :: fine_lo(3), fine_hi(3)
  Integer :: iFineX, nFineX(3), nX(3)
  Integer :: u_lo(3), u_hi(3)

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  Real    :: U_Crse(THORNADO_FLUID_NDOF)
  Real    :: U_Fine(THORNADO_FLUID_NDOF)

  !-----Being executable code.

  !print*,'I am in restriction dg!!!!'

  nFineX = 1
  nFineX(1:NDIM) = refine_factor

  ifl = 1+nguard
  ifu = nxb+nguard
  jfl = 1+nguard*k2d
  jfu = nyb+nguard*k2d
  kfl = 1+nguard*k3d
  kfu = nzb+nguard*k3d

  fine_lo = [ ifl, jfl, kfl ] ! [ 5, 5, 5]
  fine_hi = [ ifu, jfu, kfu ] ! [ 20, 20, 20]

  nX = 1
  nX(1:NDIM) = ( fine_hi(1:NDIM) - fine_lo(1:NDIM) + 1 ) / THORNADO_NNODESX ! [ 8, 8, 8 ]

  u_lo = 1
  u_hi = nX / nFineX ! [ 4, 4, 4 ]

  ! loop over thornado elements in parent block
  do kc = u_lo(3), u_hi(3)
     do jc = u_lo(2), u_hi(2)
        do ic = u_lo(1), u_hi(1) ! 1, 2, 3, 4

           ! unk offsets of first child element in child block
           k0 = 1 + refine_factor*THORNADO_NNODESX*(kc-1)*k3d
           j0 = 1 + refine_factor*THORNADO_NNODESX*(jc-1)*k2d
           i0 = 1 + refine_factor*THORNADO_NNODESX*(ic-1)     ! 1, 5, 9, 13

           U_Crse = 0.0

           ! loop over child elements for this parent element
           iFineX = 0
           do kcc = 1, nFineX(3)
              do jcc = 1, nFineX(2)
                 do icc = 1, nFineX(1)
                    iFineX = iFineX + 1

                    ! unk offsets for child block
                    k = fine_lo(3) + (k0-1) + THORNADO_NNODESX*(kcc-1)
                    j = fine_lo(2) + (j0-1) + THORNADO_NNODESX*(jcc-1)
                    i = fine_lo(1) + (i0-1) + THORNADO_NNODESX*(icc-1) ! icc = 1 : 5, 9, 13, 17
                    ! icc = 2 : 7, 11, 15, 19

                    ! grab the data from fine grid element quadrature points
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                       jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                       ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i
                       U_Fine(iNodeX) = datain(ivar,ii,jj,kk)
                    end do

                    ! compute contribution to parent element from child element
                    U_Crse = U_Crse + VolumeRatio * MATMUL( ProjectionMatrix_T(:,:,iFineX), U_Fine ) / WeightsX_q

                 end do
              end do
           end do

           ! get unk indices for (temporary) parent block
           ! note the extra factor of 2 here is because Paramesh expects the parent block data
           ! to be returned in every other cell
           k1 = kfl + 2*THORNADO_NNODESX*(kc-1)*k3d
           j1 = jfl + 2*THORNADO_NNODESX*(jc-1)*k2d
           i1 = ifl + 2*THORNADO_NNODESX*(ic-1)     ! 1, 5, 9, 13

           ! store the result in parent block
           do iNodeX = 1, THORNADO_FLUID_NDOF
              kk = 2*mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
              jj = 2*mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
              ii = 2*mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1
              dataout(ivar,ii,jj,kk) = U_Crse(iNodeX)
           end do

        end do
     end do
  end do

  Return
End Subroutine amr_restrict_unk_dg
