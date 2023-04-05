!!****if* RadTrans/RadTransMain/TwoMoment/Thornado/RadTrans_restrictDgData_simple
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
!!  RadTrans_restrictDgData_simple
!!
!! SYNOPSIS
!!
!!  call RadTrans_restrictDgData(real(IN)    :: inData(:,:,:),
!!                               real(INOUT) :: outData(:,:,:))
!!
!!  call RadTrans_restrictDgData_simple(real(IN)    :: inData(:,:,:),
!!                                      real(INOUT) :: outData(:,:,:))
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   inData : real input array, may be a slice corresponding to a region of cells
!!            for one variable, taken from a larger array
!!
!!   outData : real output array, may be a slice corresponding to a region of cells
!!            for one variable from a larger array
!!
!! AUTOGENROBODOC
!!
!! NOTES
!!  The specific subroutine implemented here can be invoked by the generic name
!!  RadTrans_restrictDgData if the caller uses the generic interface definition
!!  in the RadTrans_interface module.
!!
!! SEE ALSO
!!  RadTrans_restrictDgData
!!
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2023-03-16 Named RadTrans_restrictDgData_simple - Klaus Weide
!!
!!***

#include "Simulation.h"

subroutine RadTrans_restrictDgData_simple(inData,outData,lmask)

  Use TwoMoment_MeshRefinementModule, Only : &
     CoarsenX_TwoMoment
  Use ArrayUtilitiesModule, Only : &
     CreatePackIndex

  implicit none
  real,intent(IN)    :: inData(:,:,:,:)
  real,intent(INOUT) :: outData(:,:,:,:)
  logical,intent(IN) :: lmask(:)

  !-----Local variables
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: ii, jj, kk, icc, jcc, kcc

  Integer :: iX1, iX2, iX3

  Integer :: iNodeX
  Integer :: iFineX, nFineX(3), nFine

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  real, allocatable :: U_Fine(:,:,:,:,:,:)
  real, allocatable :: U_Crse(:,:,:,:,:)

  integer :: nX(3)
  integer :: lo(3), hi(3)

  integer :: ivar_unk2dg(size(lmask))
  integer :: ivar_dg2unk(size(lmask))
  integer :: ivar, ivar_dg, nvar_dg

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    lmask ) &
  !$OMP MAP( alloc: ivar_unk2dg, ivar_dg2unk )
#elif defined( THORNADO_OACC   )
  !$ACC ENTER DATA &
  !$ACC COPYIN(     lmask ) &
  !$ACC CREATE(     ivar_unk2dg, ivar_dg2unk )
#endif

  CALL CreatePackIndex( lmask, nvar_dg, ivar_unk2dg, ivar_dg2unk )

  nFineX = 1
  nFineX(1:NDIM) = refine_factor !or use Thornado-native variables here

  nFine = product( nFineX )

  lo = [ lbound(inData,2), lbound(inData,3), lbound(inData,4) ]
  hi = [ ubound(inData,2), ubound(inData,3), ubound(inData,4) ]

  nX = 1
  nX(1:NDIM) = ( hi(1:NDIM) - lo(1:NDIM) + 1 ) / ( THORNADO_NNODESX * nFineX(1:NDIM) )

  allocate( U_Crse(THORNADO_FLUID_NDOF,nX(1),nX(2),nX(3),nvar_dg) )
  allocate( U_Fine(THORNADO_FLUID_NDOF,nFine,nX(1),nX(2),nX(3),nvar_dg) )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    inData, outData, nFineX, nX ) &
  !$OMP MAP( alloc: U_Crse, U_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC ENTER DATA &
  !$ACC COPYIN(     inData, outData, nFineX, nX ) &
  !$ACC CREATE(     U_Crse, U_Fine )
#endif

  ! loop over fine grid elements
#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(8) &
  !$OMP PRIVATE( i1, j1, k1, iFineX, ivar, &
  !$OMP          i0, j0, k0, i, j, k, ii, jj, kk )
#elif defined( THORNADO_OACC   )
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(8) &
  !$ACC PRIVATE( i1, j1, k1, iFineX, ivar, &
  !$ACC          i0, j0, k0, i, j, k, ii, jj, kk )
#elif defined( THORNADO_OMP    )
  !$OMP PARALLEL DO COLLAPSE(8) &
  !$OMP PRIVATE( i1, j1, k1, iFineX, ivar, &
  !$OMP          i0, j0, k0, i, j, k, ii, jj, kk )
#endif
  do ivar_dg = 1, nvar_dg

     do iX3 = 1, nX(3)
        do iX2 = 1, nX(2)
           do iX1 = 1, nX(1)

              do kcc = 1, nFineX(3)
                 do jcc = 1, nFineX(2)
                    do icc = 1, nFineX(1)

                       ! store the result in child block
                       do iNodeX = 1, THORNADO_FLUID_NDOF

                          ! which child element of the parent element we are considering
                          iFineX =     icc &
                                   + ( jcc - 1 ) * nFineX(1) & 
                                   + ( kcc - 1 ) * nFineX(1) * nFineX(2)

                          ! offsets for the parent (iX1,iX2,iX3) element in unk
                          i1 = ( iX1 - 1 ) * THORNADO_NNODESX + 1
                          j1 = ( iX2 - 1 ) * THORNADO_NNODESX + 1
                          k1 = ( iX3 - 1 ) * THORNADO_NNODESX + 1

                          ! offsets of first child element in output data
                          k0 = 1 + (nFineX(3) * (k1-1)) * K3D
                          j0 = 1 + (nFineX(2) * (j1-1)) * K2D
                          i0 = 1 +  nFineX(1) * (i1-1)

                          ! calculate unk indices in child block
                          k = k0 + THORNADO_NNODESX*(kcc-1)
                          j = j0 + THORNADO_NNODESX*(jcc-1)
                          i = i0 + THORNADO_NNODESX*(icc-1)

                          kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                          jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                          ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i

                          ivar = ivar_dg2unk(ivar_dg)
                          U_Fine(iNodeX,iFineX,iX1,iX2,iX3,ivar_dg) &
                             = inData(ivar,ii,jj,kk)

                       end do

                    end do
                 end do
              end do

           end do
        end do
     end do

  end do

  ! compute coarse grid element reconstruction
  CALL CoarsenX_TwoMoment( nX, nvar_dg, U_Fine, U_Crse )

  ! loop over coarse (thornado) elements in parent block
#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
  !$OMP PRIVATE( i1, j1, k1, kk, jj, ii, ivar )
#elif defined( THORNADO_OACC   )
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
  !$ACC PRIVATE( i1, j1, k1, kk, jj, ii, ivar )
#elif defined( THORNADO_OMP    )
  !$OMP PARALLEL DO COLLAPSE(5) &
  !$OMP PRIVATE( i1, j1, k1, kk, jj, ii, ivar )
#endif
  do ivar_dg = 1, nvar_dg

     do iX3 = 1, nX(3)
        do iX2 = 1, nX(2)
           do iX1 = 1, nX(1)

              do iNodeX = 1, THORNADO_FLUID_NDOF

                 i1 = ( iX1 - 1 ) * THORNADO_NNODESX + 1
                 j1 = ( iX2 - 1 ) * THORNADO_NNODESX + 1
                 k1 = ( iX3 - 1 ) * THORNADO_NNODESX + 1

                 kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
                 jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
                 ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1

                 ivar = ivar_dg2unk(ivar_dg)
                 outData(ivar,ii,jj,kk) &
                    = U_Crse(iNodeX,iX1,iX2,iX3,ivar_dg)

              end do

           end do
        end do
     end do

  end do

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from:    outData ) &
  !$OMP MAP( release: inData, nFineX, nX, U_Crse, U_Fine, lmask, ivar_unk2dg, ivar_dg2unk )
#elif defined( THORNADO_OACC   )
  !$ACC EXIT DATA &
  !$ACC COPYOUT(      outData ) &
  !$ACC DELETE(       inData, nFineX, nX, U_Crse, U_Fine, lmask, ivar_unk2dg, ivar_dg2unk )
#endif

  deallocate( U_Crse )
  deallocate( U_Fine )

end subroutine RadTrans_restrictDgData_simple
