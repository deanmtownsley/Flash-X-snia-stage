!!****if* source/physics/Hydro/HydroMain/Spark/hy_reconstruct
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
!!  For this function the name and argument lists depend upon the variants used.
!!
!!
!!***

subroutine hy_reconstruct(i1,i2,i3,v,hy_rope,hy_uPlus, hy_uMinus, hy_flat)
#include "Simulation.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS
  
  implicit none
  
  integer, intent(IN) :: i1,i2,i3,v
  real,intent(IN) :: hy_rope(:,:,:,:), hy_flat(:,:,:)
  real,intent(OUT), dimension(:,:,:,:) :: hy_uMinus, hy_uPlus
  real :: invSumAlpha
  real, dimension(NRECON,3) :: W5p, W5m, betaWeno, Alpha5, omega
  integer :: g
  real, dimension(NRECON) :: mags, betaSum
  !! Set WENO5 coefficients once and for all
  !u_{1,i+1/2}= 2/6*u_{i-2} -7/6*u_{i-1} +11/6*u_{i}
  real, dimension(3), parameter :: coeff1p1(1:3) = (/ 2./6., -7./6., 11./6./)
  !u_{2,i+1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
  real, dimension(3), parameter :: coeff1p2(1:3) = (/-1./6.,  5./6.,  2./6./)
  !u_{3,i+1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
  real, dimension(3), parameter :: coeff1p3(1:3) = (/ 2./6.,  5./6., -1./6./)
  !(gamma1,gamma2,gamma3)
  real, dimension(3), parameter :: coeff2p(1:3)   = (/0.1, 0.6, 0.3/)
  !u_{1,i-1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
  real, dimension(3), parameter :: coeff1m1(1:3) = (/-1./6.,  5./6.,  2./6./)
  !u_{2,i-1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
  real, dimension(3), parameter :: coeff1m2(1:3) = (/ 2./6.,  5./6., -1./6./)
  !u_{3,i-1/2}=11/6*u_{i-2} -7/6*u_{i-1} + 2/6*u_{i}
  real, dimension(3), parameter :: coeff1m3(1:3) = (/ 11./6.,-7./6.,  2./6./)
  !(gamma1,gamma2,gamma3)
  real, dimension(3), parameter :: coeff2m(1:3)   = (/0.3, 0.6, 0.1/)
  
  real, parameter :: epsilon = 1.e-36
  real, parameter :: n13o12 = 13./12.
  
   
  ! Interpolation stencil for weno
  !! Calculate interface values at i+1/2
  W5p(v,1) = coeff1p1(1)*hy_rope(v,i1-2,i2,i3) + coeff1p1(2)*hy_rope(v,i1-1,i2,i3) + coeff1p1(3)*hy_rope(v,i1+0,i2,i3)
  W5p(v,2) = coeff1p2(1)*hy_rope(v,i1-1,i2,i3) + coeff1p2(2)*hy_rope(v,i1+0,i2,i3)   + coeff1p2(3)*hy_rope(v,i1+1,i2,i3)
  W5p(v,3) = coeff1p3(1)*hy_rope(v,i1+0,i2,i3)   + coeff1p3(2)*hy_rope(v,i1+1,i2,i3) + coeff1p3(3)*hy_rope(v,i1+2,i2,i3)
  
  !! Calculate interface values at i-1/2
  W5m(v,1) = coeff1m1(1)*hy_rope(v,i1-2,i2,i3) + coeff1m1(2)*hy_rope(v,i1-1,i2,i3) + coeff1m1(3)*hy_rope(v,i1+0,i2,i3)
  W5m(v,2) = coeff1m2(1)*hy_rope(v,i1-1,i2,i3) + coeff1m2(2)*hy_rope(v,i1+0,i2,i3)   + coeff1m2(3)*hy_rope(v,i1+1,i2,i3)
  W5m(v,3) = coeff1m3(1)*hy_rope(v,i1+0,i2,i3)   + coeff1m3(2)*hy_rope(v,i1+1,i2,i3) + coeff1m3(3)*hy_rope(v,i1+2,i2,i3)
  
  !! Calculate smoothness indicators at i+1/2
  betaWeno(v,1) = n13o12*(hy_rope(v,i1-2,i2,i3) - 2.*hy_rope(v,i1-1,i2,i3) +    hy_rope(v,i1+0,i2,i3)  )**2 &
       +            0.25*(hy_rope(v,i1-2,i2,i3) - 4.*hy_rope(v,i1-1,i2,i3) + 3.*hy_rope(v,i1,i2,i3)  )**2
  betaWeno(v,2) = n13o12*(hy_rope(v,i1-1,i2,i3) - 2.*hy_rope(v,i1,i2,i3)   +    hy_rope(v,i1+1,i2,i3))**2 &
       +            0.25*(hy_rope(v,i1-1,i2,i3)                      -    hy_rope(v,i1+1,i2,i3))**2
  betaWeno(v,3) = n13o12*(hy_rope(v,i1,i2,i3)   - 2.*hy_rope(v,i1+1,i2,i3) +    hy_rope(v,i1+2,i2,i3))**2 &
       +            0.25*(3.*hy_rope(v,i1,i2,i3)- 4.*hy_rope(v,i1+1,i2,i3) +    hy_rope(v,i1+2,i2,i3))**2
  
  !! Use problem-adaptive epsilong as in Tchekovskoy7, A3
  ! This does not seem to work with the WENO-Z indicators of Borges+08
  ! mags(v) = hy_rope(@M ind_m(2))**2 + hy_rope(@M ind_m(1))**2 + hy_rope(@M ind_p(0))**2 &
  !      + hy_rope(@M ind_p(1))**2 + hy_rope(@M ind_p(2))**2
  ! betaWeno(v,1) = betaWeno(v,1) + epsilon*mags(v) + TINY(1.0)
  ! betaWeno(v,2) = betaWeno(v,2) + epsilon*mags(v) + TINY(1.0)
  ! betaWeno(v,3) = betaWeno(v,3) + epsilon*mags(v) + TINY(1.0)
  
  !! This is WENO-Zv this is very similar to weno5 with wenoExp=1
  Alpha5(v,1) = coeff2p(1)*(1.+(abs(betaWeno(v,1)-betaWeno(v,3))/(betaWeno(v,1)+epsilon)))
  Alpha5(v,2) = coeff2p(2)*(1.+(abs(betaWeno(v,1)-betaWeno(v,3))/(betaWeno(v,2)+epsilon)))
  Alpha5(v,3) = coeff2p(3)*(1.+(abs(betaWeno(v,1)-betaWeno(v,3))/(betaWeno(v,3)+epsilon)))
  
  !! Normalize weights at i+1/2
  invSumAlpha = 1./(Alpha5(v,1)+Alpha5(v,2)+Alpha5(v,3))
  omega(v,1)  = Alpha5(v,1)*invSumAlpha
  omega(v,2)  = Alpha5(v,2)*invSumAlpha
  omega(v,3)  = Alpha5(v,3)*invSumAlpha
  
  !! Compute interface value at i+1/2
  hy_uPlus(v ,i1,i2,i3)  = omega(v,1)*W5p(v,1) + omega(v,2)*W5p(v,2) + omega(v,3)*W5p(v,3)
  !! Apply hy_flattening
  hy_uPlus(v ,i1,i2,i3) = hy_flat(i1,i2,i3)*hy_uPlus(v ,i1,i2,i3) + (1.-hy_flat(i1,i2,i3))*hy_rope(v,i1 ,i2,i3)
  
  !! Now move on to i-1/2
  !! This is WENO-Z
  Alpha5(v,1) = coeff2m(1)*(1.+(abs(betaWeno(v,1)-betaWeno(v,3))/(betaWeno(v,1)+epsilon)))
  Alpha5(v,2) = coeff2m(2)*(1.+(abs(betaWeno(v,1)-betaWeno(v,3))/(betaWeno(v,2)+epsilon)))
  Alpha5(v,3) = coeff2m(3)*(1.+(abs(betaWeno(v,1)-betaWeno(v,3))/(betaWeno(v,3)+epsilon)))
  
  !! Normalize weights at i-1/2
  invSumAlpha = 1./(Alpha5(v,1)+Alpha5(v,2)+Alpha5(v,3))
  omega(v,1)  = Alpha5(v,1)*invSumAlpha
  omega(v,2)  = Alpha5(v,2)*invSumAlpha
  omega(v,3)  = Alpha5(v,3)*invSumAlpha
  
  !! Compute interface value at i-1/2
  hy_uMinus(v ,i1,i2,i3) = omega(v,1)*W5m(v,1) + omega(v,2)*W5m(v,2) + omega(v,3)*W5m(v,3)
  !! Apply hy_flattening
  hy_uMinus(v ,i1,i2,i3) = hy_flat(i1,i2,i3)*hy_uMinus(v ,i1,i2,i3) + (1.-hy_flat(i1,i2,i3))*hy_rope(v, i1 ,i2,i3)
  !! Check for monotonicity
  if ( (hy_uPlus(v ,i1,i2,i3)-hy_rope(v, i1 ,i2,i3))*(hy_rope(v, i1 ,i2,i3)-hy_uMinus(v ,i1,i2,i3)) <= 0. ) then
     hy_uPlus(v ,i1,i2,i3)  = hy_rope(v, i1 ,i2,i3)
     hy_uMinus(v ,i1,i2,i3) = hy_rope(v, i1 ,i2,i3)
  end if
  
  
end  subroutine hy_reconstruct


