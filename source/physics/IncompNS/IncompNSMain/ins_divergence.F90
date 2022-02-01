!!
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
!!
!!***
#include "Simulation.h"
#include "constants.h"

SUBROUTINE ins_divergence(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
         dx,dy,dz,divv)

      ! This routine computes the divergence of the velocity field.

      implicit none

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv


      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
         ( uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
           uni(ix1:ix2,jy1:jy2,kz1:kz2) )/dx +  &       
         ( vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
           vni(ix1:ix2,jy1:jy2,kz1:kz2) )/dy


#if NDIM == MDIM
      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
           divv(ix1:ix2,jy1:jy2,kz1:kz2)     +  &
         ( wni(ix1:ix2,jy1:jy2,kz1+1:kz2+1) -   &
           wni(ix1:ix2,jy1:jy2,kz1:kz2) )/dz 
#endif

END SUBROUTINE ins_divergence
