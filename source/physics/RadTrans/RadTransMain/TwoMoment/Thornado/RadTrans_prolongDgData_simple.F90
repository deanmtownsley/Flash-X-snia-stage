!!****if* RadTrans/RadTransMain/TwoMoment/Thornado/RadTrans_prolongDgData_simple
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
!!  RadTrans_prolongDgData_simple
!!
!! SYNOPSIS
!!
!!  call RadTrans_prolongDgData(real(IN)   ,dimension(:,:,:) :: inData(:,:,:),
!!                              real(INOUT),dimension(:,:,:) :: outData(:,:,:),
!!                              integer(IN),dimension(MDIM)  :: skip(3))
!!  call RadTrans_prolongDgData_simple(real(IN)   ,dimension(:,:,:) :: inData(:,:,:),
!!                              real(INOUT),dimension(:,:,:) :: outData(:,:,:),
!!                              integer(IN),dimension(MDIM)  :: skip(3))
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
!!   skip : integer array, its values should be in the range
!!            0 ...  refine_factor*THORNADO_NNODESX - 1
!!          for the NDIM active spatial directions.
!!          For each spatial direction, it indicates by how much the first output
!!          element in that direction is offset wrt the first input element.
!!
!! AUTOGENROBODOC
!!
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2022-09-22 Added skip to the interface          - Klaus Weide
!!  2023-03-14 Named RadTrans_prolongDgData_simple  - Klaus Weide
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine RadTrans_prolongDgData_simple(inData,outData,skip,lmask)

  Use TwoMoment_MeshRefinementModule, Only : &
     RefineX_TwoMoment
  Use ArrayUtilitiesModule, Only : &
     CreatePackIndex

  implicit none
  real,intent(IN)    :: inData(:,:,:,:)
  real,intent(INOUT) :: outData(:,:,:,:)
  integer,intent(IN) :: skip(MDIM)
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
  nX(1:NDIM) = ( hi(1:NDIM) - lo(1:NDIM) + 1 ) / THORNADO_NNODESX

  allocate( U_Crse(THORNADO_FLUID_NDOF,nX(1),nX(2),nX(3),nvar_dg) )
  allocate( U_Fine(THORNADO_FLUID_NDOF,nFine,nX(1),nX(2),nX(3),nvar_dg) )

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    inData, outData, skip, nFineX, nX ) &
  !$OMP MAP( alloc: U_Crse, U_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC ENTER DATA &
  !$ACC COPYIN(     inData, outData, skip, nFineX, nX ) &
  !$ACC CREATE(     U_Crse, U_Fine )
#endif

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
                U_Crse(iNodeX,iX1,iX2,iX3,ivar_dg) &
                   = inData(ivar,ii,jj,kk)

             end do

          end do
       end do
    end do

  end do

  ! compute fine grid element projection
  CALL RefineX_TwoMoment( nX, nvar_dg, U_Crse, U_Fine )

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

                         ! offsets for the parent (iX1,iX2,iX3) element in unk
                         i1 = ( iX1 - 1 ) * THORNADO_NNODESX + 1
                         j1 = ( iX2 - 1 ) * THORNADO_NNODESX + 1
                         k1 = ( iX3 - 1 ) * THORNADO_NNODESX + 1

                         ! which child element of the parent element we are considering
                         iFineX =     icc &
                                  + ( jcc - 1 ) * nFineX(1) & 
                                  + ( kcc - 1 ) * nFineX(1) * nFineX(2)

                         ! offsets of first child element in output data
                         k0 = 1 + (nFineX(3) * (k1-1) - skip(3)) * K3D
                         j0 = 1 + (nFineX(2) * (j1-1) - skip(2)) * K2D
                         i0 = 1 +  nFineX(1) * (i1-1) - skip(1)

                         ! calculate unk indices in child block
                         k = k0 + THORNADO_NNODESX*(kcc-1)
                         j = j0 + THORNADO_NNODESX*(jcc-1)
                         i = i0 + THORNADO_NNODESX*(icc-1)

                         kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                         jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                         ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i

                         if (      kk.GE.lbound(outData,4) .AND. kk.LE.ubound(outData,4) &
                             .AND. jj.GE.lbound(outData,3) .AND. jj.LE.ubound(outData,3) &
                             .AND. ii.GE.lbound(outData,2) .AND. ii.LE.ubound(outData,2) ) then

                            ivar = ivar_dg2unk(ivar_dg)
                            outData(ivar,ii,jj,kk) &
                               = U_Fine(iNodeX,iFineX,iX1,iX2,iX3,ivar_dg)

                         end if
                      end do

                   end do
                end do
             end do

           end do
        end do
     end do

  end do

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from:    outData ) &
  !$OMP MAP( release: inData, skip, nFineX, nX, U_Crse, U_Fine, lmask, ivar_unk2dg, ivar_dg2unk )
#elif defined( THORNADO_OACC   )
  !$ACC EXIT DATA &
  !$ACC COPYOUT(      outData ) &
  !$ACC DELETE(       inData, skip, nFineX, nX, U_Crse, U_Fine, lmask, ivar_unk2dg, ivar_dg2unk )
#endif

  deallocate( U_Crse )
  deallocate( U_Fine )

end subroutine RadTrans_prolongDgData_simple
