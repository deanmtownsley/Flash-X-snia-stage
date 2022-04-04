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
  !***************************************************************
  ! This subroutine computes the discretization of the diffusion RHS of the 
  ! Helmholtz equation on a staggered uniform grid.
  !
  ! Input:  uni,vni     = velocity at timestep n
  !         ru1         = molecular viscosity !- kpd - Inverse Reynolds No
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         dx,dy       = grid spacing in x and y directions
  !
  ! Output: ru,rv    = u and v momentum for Helmholtz RHS
  !**************************************************************
SUBROUTINE ins_diffusion2d_vardens(uni,vni,ru1,ix1,ix2,jy1,jy2,dx,dy,ru,rv, &
                           visc,rhox,rhoy)
      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
      REAL, INTENT(IN):: ru1, dx, dy
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, visc, rhox, rhoy
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv

      INTEGER:: i, j
      REAL:: dx1, dy1, Mdens, th, cri, crc
      ! x-component variables
      REAL:: dudxp, dudxm, dudyp, dudym, dvdxp, dvdxm
      REAL:: txxp, txxm, tyyp, tyym
      REAL:: txyp, txym
      ! new y-component variables
      REAL:: dvdyp, dvdym
      INTEGER, parameter :: kz1 = 1

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy

      !++++++++++  U-COMPONENT  ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2+1 
            ! Diffusion Terms

             ! get derivatives at 1/2 locations

             ! velocity
             dudxp = (uni(i+1,j,kz1) - uni(i,j,kz1))*dx1
             dudxm = (uni(i,j,kz1)   - uni(i-1,j,kz1))*dx1
             dudyp = (uni(i,j+1,kz1) - uni(i,j,kz1))*dy1
             dudym = (uni(i,j,kz1)   - uni(i,j-1,kz1))*dy1

             !-Variable Viscosity Implementation (ru1 is 1/Re)
             txxp = ru1*visc(i  ,j,kz1)*(dudxp)
             txxm = ru1*visc(i-1,j,kz1)*(dudxm)
             tyyp = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i-1,j+1,kz1)+visc(i,j+1,kz1))*(dudyp)
             tyym = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i-1,j-1,kz1)+visc(i,j-1,kz1))*(dudym)

             Mdens = rhox(i,j,kz1) ! Inverse density

             ! calculate RHS for u-momentum
             ru(i,j,kz1) = ru(i,j,kz1)                               &
                           + Mdens*(txxp - txxm)*dx1                 & ! diffusion - normal terms 
                           + Mdens*(tyyp - tyym)*dy1                                  

          enddo
       enddo


    !++++++++++  V-COMPONENT  ++++++++++
       do j = jy1,jy2+1
          do i = ix1,ix2
             ! Diffusion Terms

             ! get derivatives at 1/2 locations
   
             ! velocity
             dvdxp = (vni(i+1,j,kz1) - vni(i,j,kz1))*dx1
             dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1))*dx1
             dvdyp = (vni(i,j+1,kz1) - vni(i,j,kz1))*dy1
             dvdym = (vni(i,j,kz1) - vni(i,j-1,kz1))*dy1

             !- Variable Viscosity Implementation (ru1 is 1/Re)
             txxp = ru1*0.25*(visc(i,j,kz1)+visc(i+1,j,kz1)+visc(i,j-1,kz1)+visc(i+1,j-1,kz1))*(dvdxp)
             txxm = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i,j-1,kz1)+visc(i-1,j-1,kz1))*(dvdxm)
             tyyp = ru1*visc(i,j,kz1)  *(dvdyp)
             tyym = ru1*visc(i,j-1,kz1)*(dvdym)

             Mdens = rhoy(i,j,kz1) ! Inverse density.

             ! calculate RHS for v-momentum
             rv(i,j,kz1) = rv(i,j,kz1)                                &
                           + Mdens* (txxp - txxm)*dx1                 &! diffusion - normal terms
                           + Mdens* (tyyp - tyym)*dy1 
          enddo
       enddo

END SUBROUTINE ins_diffusion2d_vardens

!-------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------!

SUBROUTINE ins_diffusion3d_vardens(uni,vni,wni,tv,ru1,      &
                        ix1,ix2,jy1,jy2,kz1,kz2, &
                        dx,dy,dz,ru,rv,rw,visc,  &
                        rhox,rhoy,rhoz)

  !*****************************************************************
  ! This subroutine computes the centered discretization of the RHS 
  ! of the momentum equation (viscous terms) on a 
  ! staggered uniform grid based on the Paramesh grid structure.
  !
  ! Input:  uni,vni,wni = velocity at timestep n
  !         tv          = eddy viscosity
  !         ru1         = molecular viscosity !- kpd - Inverse Reynolds No
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         kz1,kz2     = starting and ending z indices
  !         dx,dy,dz    = grid spacing in x, y, and z directions
  !
  ! Output: ru,rv,rw    = RHS of u, v, and w momentum equations
  !
  ! E. Balaras   July 1999
  ! P. Rabenold  August 2006
  !**************************************************************

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
      REAL, INTENT(IN):: ru1, dx, dy, dz
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, wni, tv, visc, rhox, rhoy
      REAL, DIMENSION(:,:,:), INTENT(IN):: rhoz
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv, rw

      INTEGER:: i, j, k
      REAL:: dx1, dy1, dz1, Mdens
      ! x-component variables
      REAL:: dudxp, dudxm, dudyp, dudym, dudzp, dudzm, dvdxp, dvdxm, &
             dwdxp, dwdxm
      REAL:: tvjp, tvjm, tvkp, tvkm
      REAL:: txxp, txxm, tyyp, tyym, tzzp, tzzm
      REAL:: txyp, txym, txzp, txzm
      ! additional y-component variables
      REAL:: dvdyp, dvdym, dvdzp, dvdzm, dwdyp, dwdym
      REAL:: tvip, tvim
      REAL:: tyzp, tyzm
      ! additional z-component variables
      REAL:: wzplus, wzminus
      REAL:: dwdzp, dwdzm
      REAL:: vvip,vvim,vvjp,vvjm,vvkp,vvkm

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy
      dz1 = 1.0/dz

      !++++++++++  U-COMPONENT (Variable Density)  ++++++++++
      do k = kz1,kz2
         do j = jy1,jy2
            do i = ix1,ix2+1
               ! Diffusion Terms

               ! get derivatives at 1/2 locations
               dudxp = (uni(i+1,j  ,k  ) - uni(i  ,j  ,k  ))*dx1
               dudxm = (uni(i  ,j  ,k  ) - uni(i-1,j  ,k  ))*dx1
               dudyp = (uni(i  ,j+1,k  ) - uni(i  ,j  ,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dudzp = (uni(i  ,j  ,k+1) - uni(i  ,j  ,k  ))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdxp = (vni(i  ,j+1,k  ) - vni(i-1,j+1,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dwdxp = (wni(i  ,j  ,k+1) - wni(i-1,j  ,k+1))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1

               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvjp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i-1,j+1,k  ))
               tvjm = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j  ,k+1) + tv(i-1,j  ,k+1))
               tvkm = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j,  k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = visc(i,j,k)
               vvim = visc(i-1,j,k)
               vvjp = 0.25*(visc(i-1,j  ,k  ) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j+1,k  ) + visc(i-1,j+1,k  ))
               vvjm = 0.25*(visc(i-1,j-1,k  ) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvkp = 0.25*(visc(i-1,j  ,k  ) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j  ,k+1) + visc(i-1,j  ,k+1))
               vvkm = 0.25*(visc(i-1,j  ,k-1) + visc(i  ,j,  k-1) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))


               ! flux of normal total stresses, mu*dU/dXj (ru1 is 1/Re)
               txxp = (ru1*vvip + 2.0*tv(i,j,k))   *(dudxp)    ! KPD 2*nut, but not nu 
               txxm = (ru1*vvim + 2.0*tv(i-1,j,k)) *(dudxm)    ! KPD 2*nut, but not nu 
               tyyp = (ru1*vvjp + tvjp)            *(dudyp)
               tyym = (ru1*vvjm + tvjm)            *(dudym)
               tzzp = (ru1*vvkp + tvkp)            *(dudzp)
               tzzm = (ru1*vvkm + tvkm)            *(dudzm)

               ! flux of cross SGS stresses, mu*dUi/dX
               txyp = tvjp*dvdxp                             ! KPD no mol visc, turb only 
               txym = tvjm*dvdxm                             ! KPD no mol visc, turb only 
               txzp = tvkp*dwdxp                             ! KPD no mol visc, turb only 
               txzm = tvkm*dwdxm                             ! KPD no mol visc, turb only 

               !- kpd - Mixture inverse density
               Mdens = rhox(i,j,k)

               ! calculate RHS for u-momentum
               ru(i,j,k) = ru(i,j,k)                                &
                           + Mdens*(txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens*(tyyp - tyym)*dy1                &
                           + Mdens*(tzzp - tzzm)*dz1                &
                           + Mdens*(txyp - txym)*dy1                &! TURBULENT cross terms
                           + Mdens*(txzp - txzm)*dz1
 
            enddo
         enddo
      enddo

      !++++++++++  V-COMPONENT  ++++++++++

      do k = kz1,kz2
         do j = jy1,jy2+1
            do i = ix1,ix2
               ! Diffusion Terms

               ! get derivatives at 1/2 locations
               dvdxp = (vni(i+1,j  ,k  ) - vni(i  ,j  ,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dvdyp = (vni(i  ,j+1,k  ) - vni(i  ,j  ,k  ))*dy1
               dvdym = (vni(i  ,j  ,k  ) - vni(i  ,j-1,k  ))*dy1
               dvdzp = (vni(i  ,j  ,k+1) - vni(i  ,j  ,k  ))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1
               dudyp = (uni(i+1,j  ,k  ) - uni(i+1,j-1,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dwdyp = (wni(i  ,j  ,k+1) - wni(i  ,j-1,k+1))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1

               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvip = 0.25*(tv(i  ,j-1,k  ) + tv(i+1,j-1,k  ) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i  ,j-1,k  ) + tv(i  ,j-1,k+1) + &
                            tv(i  ,j  ,k+1) + tv(i  ,j  ,k  ))
               tvkm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = 0.25*(visc(i  ,j-1,k  ) + visc(i+1,j-1,k  ) + &
                            visc(i+1,j  ,k  ) + visc(i  ,j  ,k  ))
               vvim = 0.25*(visc(i-1,j-1,k  ) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvjp = visc(i,j,k) 
               vvjm = visc(i,j-1,k)
               vvkp = 0.25*(visc(i  ,j-1,k  ) + visc(i  ,j-1,k+1) + &
                            visc(i  ,j  ,k+1) + visc(i  ,j  ,k  ))
               vvkm = 0.25*(visc(i  ,j-1,k-1) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i  ,j  ,k-1))

               ! flux of normal total stresses (ru1 is 1/Re)
               txxp = (ru1*vvip + tvip           )*(dvdxp)
               txxm = (ru1*vvim + tvim           )*(dvdxm)
               tyyp = (ru1*vvjp + 2.0*tv(i,j,k)  )*(dvdyp)
               tyym = (ru1*vvjm + 2.0*tv(i,j-1,k))*(dvdym)
               tzzp = (ru1*vvkp + tvkp           )*(dvdzp)
               tzzm = (ru1*vvkm + tvkm           )*(dvdzm)

               ! flux of cross SGS stresses
               txyp = tvip*dudyp
               txym = tvim*dudym
               tyzp = tvkp*dwdyp
               tyzm = tvkm*dwdym

               !- kpd - Mixture inverse density
               Mdens = rhoy(i,j,k)

               ! calculate RHS for v-momentum
               rv(i,j,k) = rv(i,j,k)                                 &
                           + Mdens* (txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens* (tyyp - tyym)*dy1                &
                           + Mdens* (tzzp - tzzm)*dz1                &
                           + Mdens* (txyp - txym)*dx1                &! diffusion - cross terms
                           + Mdens* (tyzp - tyzm)*dz1

            enddo
         enddo
      enddo

      !++++++++++  W-COMPONENT  ++++++++++
      
      do k = kz1,kz2+1
         do j = jy1,jy2
            do i = ix1,ix2
               ! Diffusion Terms

               ! get derivatives at 1/2 locations
               dwdxp = (wni(i+1,j  ,k  ) - wni(i  ,j  ,k  ))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1
               dwdyp = (wni(i  ,j+1,k  ) - wni(i  ,j  ,k  ))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1
               dwdzp = (wni(i  ,j  ,k+1) - wni(i  ,j  ,k  ))*dz1
               dwdzm = (wni(i  ,j  ,k  ) - wni(i  ,j  ,k-1))*dz1
               dudzp = (uni(i+1,j  ,k  ) - uni(i+1,j  ,k-1))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdzp = (vni(i  ,j+1,k  ) - vni(i  ,j+1,k-1))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1

              !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvip = 0.25*(tv(i  ,j  ,k-1) + tv(i+1,j  ,k-1) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j  ,k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvjp = 0.25*(tv(i  ,j  ,k-1) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i  ,j+1,k-1))
               tvjm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + & 
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = 0.25*(visc(i  ,j  ,k-1) + visc(i+1,j  ,k-1) + &
                            visc(i+1,j  ,k  ) + visc(i  ,j  ,k  ))
               vvim = 0.25*(visc(i-1,j  ,k-1) + visc(i  ,j  ,k-1) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvjp = 0.25*(visc(i  ,j  ,k-1) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j+1,k  ) + visc(i  ,j+1,k-1))
               vvjm = 0.25*(visc(i  ,j-1,k-1) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i  ,j  ,k-1))
               vvkp = visc(i,j,k)
               vvkm = visc(i,j,k-1)

               ! flux of normal total stresses (ru1 is 1/Re)
               txxp = (ru1*vvip + tvip           )*(dwdxp)
               txxm = (ru1*vvim + tvim           )*(dwdxm)
               tyyp = (ru1*vvjp + tvjp           )*(dwdyp)
               tyym = (ru1*vvjm + tvjm           )*(dwdym)
               tzzp = (ru1*vvkp + 2.0*tv(i,j,k)  )*(dwdzp)
               tzzm = (ru1*vvkm + 2.0*tv(i,j,k-1))*(dwdzm)

               ! flux of cross SGS stresses
               txzp = tvip*dudzp
               txzm = tvim*dudzm
               tyzp = tvjp*dvdzp
               tyzm = tvjm*dvdzm

               !- kpd - Mixture inverse density
               Mdens = rhoz(i,j,k)

               ! calculate RHS for w-momentum
               rw(i,j,k) = rw(i,j,k)                                 &
                           + Mdens* (txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens* (tyyp - tyym)*dy1                &
                           + Mdens* (tzzp - tzzm)*dz1                &
                           + Mdens* (txzp - txzm)*dx1                &! diffusion - cross terms
                           + Mdens* (tyzp - tyzm)*dy1
            enddo
         enddo
      enddo

END SUBROUTINE ins_diffusion3d_vardens
