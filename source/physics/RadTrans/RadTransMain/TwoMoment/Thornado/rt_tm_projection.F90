!!****f* source/physics/RadTrans/RadTransMain/TwoMoment/Thornado/rt_tm_projection
!!
!! NAME
!!
!!  rt_tm_projection
!!
!! SYNOPSIS
!!
!!  call rt_tm_projection ( real(IN),pointer :: solnData(:,:,:,:), 
!!                              integer(IN) :: nX(3),
!!                              integer(IN) :: lo(MDIM),hi(MDIM),
!!                              integer(IN) :: u_lo(MDIM),u_hi(MDIM) )
!!
!! DESCRIPTION
!!
!!  Interpolate fluid variables from DG node values from thornado 
!!  to cell-centered values for Flash.
!!  Copy radiation varibles from internal Thornado data structures.
!!
!! ARGUMENTS
!! solnData - pointer to cell-centered Flash data
!! nX - number of cells
!! lo - lower index limits (excluding guard cells)
!! hi - upper index limits (excluding guard cells)
!! u_lo - lower index limits (including guard cells)
!! u_hi - upper index limits (including guard cells)
!!
!!***

!!REORDER(4): solnData

subroutine rt_tm_projection(solnData,nX,lo,hi,u_lo,u_hi)

#include "Simulation.h"
#include "constants.h"

  use FluidFieldsModule, ONLY : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  use GeometryFieldsModule, ONLY : uGF, iGF_h_1, iGF_h_2, iGF_h_3
  use RadiationFieldsModule, ONLY : uCR
  use SubcellReconstructionModule, ONLY : ProjectionMatrix
  use UnitsModule, ONLY : Gram, Centimeter, Second, AtomicMassUnit

  implicit none

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(3), intent(in) :: nX
  integer, dimension(MDIM), intent(in) :: lo,hi,u_lo,u_hi

  integer :: i, j, k, ic, jc, kc, ii, jj, kk, n, ioff, ivar, i_sub, n_sub

  real :: D, S1, S2, S3, E, Nel, ekin
  integer :: iNode, iNodeX, iNodeE, iCF, iCR, iS, iE

  real, parameter :: conv_dens = Gram / Centimeter**3
  real, parameter :: conv_mom  = Gram / Centimeter**2 / Second
  real, parameter :: conv_enr  = Gram / Centimeter / Second**2
  real, parameter :: conv_ne   = Gram / Centimeter**3 / AtomicMassUnit
  real, parameter :: conv_J    = Gram/Second**2/Centimeter
  real, parameter :: conv_H    = Gram/Second**3

  ! Interpolate fluid vars from Thornado arrays into Flash arrays
#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(4)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG COLLAPSE(4) &
  !$ACC PRIVATE( i, j, k, ii, jj, kk, D, S1, S2, S3, E, Nel, ekin ) &
  !$ACC PRESENT( ProjectionMatrix, solnData, uCF, uGF, nX )
#elif defined(THORNADO_OMP)
  !$OMP PARALLEL DO COLLAPSE(4)
#endif
  do kc = 1, nX(3)
     do jc = 1, nX(2)
        do ic = 1, nX(1)
           do iNodeX = 1, THORNADO_FLUID_NDOF

              i = lo(IAXIS) + THORNADO_NNODESX*(ic-1)
              j = lo(JAXIS) + THORNADO_NNODESX*(jc-1)
              k = lo(KAXIS) + THORNADO_NNODESX*(kc-1)

              ii = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
              jj = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
              kk = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

              D   = 0.0
              S1  = 0.0
              S2  = 0.0
              S3  = 0.0
              E   = 0.0
              Nel = 0.0

#if defined(THORNADO_OMP_OL)
              !$OMP PARALLEL DO SIMD &
              !$OMP REDUCTION( +: D, S1, S2, S3, E, Nel )
#elif defined(THORNADO_OACC)
              !$ACC LOOP VECTOR &
              !$ACC REDUCTION( +: D, S1, S2, S3, E, Nel )
#endif
              do n = 1, THORNADO_FLUID_NDOF
                 D   = D   + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_D )
                 S1  = S1  + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_S1) &
                       / uGF(n,ic,jc,kc,iGF_h_1)
                 S2  = S2  + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_S2) &
                       / uGF(n,ic,jc,kc,iGF_h_2)
                 S3  = S3  + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_S3) &
                       / uGF(n,ic,jc,kc,iGF_h_3)
                 E   = E   + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_E )
                 Nel = Nel + ProjectionMatrix(iNodeX,n)*uCF(n,ic,jc,kc,iCF_Ne)
              end do

              solnData(DENS_VAR,  ii,jj,kk) = D   / conv_dens
              solnData(VELX_VAR,  ii,jj,kk) = S1  / conv_mom  / solnData(DENS_VAR,ii,jj,kk)
              solnData(VELY_VAR,  ii,jj,kk) = S2  / conv_mom  / solnData(DENS_VAR,ii,jj,kk)
              solnData(VELZ_VAR,  ii,jj,kk) = S3  / conv_mom  / solnData(DENS_VAR,ii,jj,kk)
              solnData(ENER_VAR,  ii,jj,kk) = E   / conv_enr  / solnData(DENS_VAR,ii,jj,kk)
              solnData(YE_MSCALAR,ii,jj,kk) = Nel / conv_ne   / solnData(DENS_VAR,ii,jj,kk)
#ifdef EINT_VAR
              ekin = 0.0
              do n = VELX_VAR, VELZ_VAR
                ekin = ekin + 0.5*solnData(n,ii,jj,kk)**2
              end do
              solnData(EINT_VAR,ii,jj,kk) = solnData(ENER_VAR,ii,jj,kk) - ekin
#endif
           end do
        end do
     end do
  end do

  ! Copy radiation vars from Thornado arrays into Flash arrays
#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(7)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG COLLAPSE(7) &
  !$ACC PRIVATE( i, j, k, ii, jj, kk, ioff, ivar, iNodeE, iNodeX ) &
  !$ACC PRESENT( solnData, uCR, nX, lo )
#elif defined(THORNADO_OMP)
  !$OMP PARALLEL DO COLLAPSE(7)
#endif
  do kc = 1, nX(3)
     do jc = 1, nX(2)
        do ic = 1, nX(1)
           do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS
              do iE = 1-THORNADO_SWE, THORNADO_NE+THORNADO_SWE
              do iNode = 1, THORNADO_RAD_NDOF

                 i = lo(IAXIS) + THORNADO_NNODESX*(ic-1)
                 j = lo(JAXIS) + THORNADO_NNODESX*(jc-1)
                 k = lo(KAXIS) + THORNADO_NNODESX*(kc-1)

                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 iNodeE = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                 iNodeX = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                 ii     = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
                 jj     = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                 kk     = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

                 ivar = ioff + iNodeE - 1

                 solnData(ivar,ii,jj,kk) = uCR(iNode,iE,ic,jc,kc,iCR,iS)

              end do
           end do ; end do ; end do
        end do
     end do
  end do


  return
end subroutine rt_tm_projection
