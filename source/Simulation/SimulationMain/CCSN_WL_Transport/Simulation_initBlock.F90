!!****if* source/Simulation/SimulationMain/CCSN/Simulation_initBlock
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
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(Grid_tile_t(IN) :: tileDesc)
!!                       
!!
!! DESCRIPTION
!!
!! Initial conditions for Core Collapse SN problem
!!
!! ARGUMENTS
!!
!!  tileDesc -           describes the block (or tile) to initialize
!!
!! NOTES
!!  
!!  This problem is described in, e.g.,
!!  Couch, S.M. 2013, ApJ, 765, 29
!!  Couch, S.M. 2013, ApJ, 775, 35
!!  Couch, S.M. & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(solnData, tileDesc)

  use Driver_interface, ONLY : Driver_abort
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getDeltas, &
       Grid_getGeometry
  use Eos_interface, ONLY : Eos_wrapped, Eos_getAbarZbar, Eos
  use Multispecies_interface, ONLY : Multispecies_getSumFrac,  Multispecies_getSumInv

  use Simulation_data
  use chimera_model_module
  use model_interp_module

#if defined(THORNADO)
  ! thornado modules

  use KindModule, ONLY : DP, SqrtTiny
  use GeometryFieldsModule, ONLY: uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  use RadiationFieldsModule, ONLY : iCR_N, iCR_G1, iCR_G2, iCR_G3
  use UnitsModule, ONLY : Centimeter, Gram, Second, MeV, Kelvin
  use MeshModule, ONLY : NodeCoordinate, MeshE, MeshX
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch  
  use TwoMoment_UtilitiesModule_OrderV, ONLY : ComputeConserved_TwoMoment
#endif

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
#include "Multispecies.h"

  real,dimension(:,:,:,:),pointer :: solnData
  type(Grid_tile_t), intent(in)   :: tileDesc

  integer :: blockID
  real, allocatable, dimension(:) :: xCenter, xLeft, xRight
  real, allocatable, dimension(:) :: yCenter, yLeft, yRight
  real, allocatable, dimension(:) :: zCenter, zLeft, zRight
  real, dimension(MDIM) :: delta
  real :: dx, dy, dz

  real, dimension(EOS_NUM) :: eosData

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits
  integer, dimension(MDIM) :: rigid_axis
  integer :: iSize, jSize, kSize
  integer :: iSizeGC, jSizeGC, kSizeGC
  integer :: level

  integer :: meshGeom

  integer :: i, j, k, n, ivar

  real :: suminv
  real :: abar,zbar,sumY

  real :: dens_interp, velx_interp, vely_interp, velz_interp
  real :: temp_interp, pres_interp, eint_interp, spec_interp(NSPECIES)
  real :: ye_interp, var_interp
  real :: radCenter, thtCenter, phiCenter
  real :: radCenterVol, thtCenterVol, phiCenterVol

  ! hydro initialization variables
  integer, dimension(1:MDIM)   :: lo, hi
  integer, dimension(1:MDIM+1) :: u_lo, u_hi

  integer :: ii, jj, kk, ic, jc, kc
  integer :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3
  integer :: iS, iCR, iE, iNode, iNodeX, iNodeE, ioff
  integer :: nX(3), swX(3)
  real    :: xL(3), xR(3)
  real    :: xnode, ynode, znode, enode
  real    :: D, I1, I2, I3
  real    :: velr, velt, velp
  real    :: Gmdd11, Gmdd22, Gmdd33
  real    :: CN, G1, G2, G3

#if defined(THORNADO)
  real, parameter :: conv_x = Centimeter
  real, parameter :: UnitD  = Gram / Centimeter**3
  real, parameter :: UnitT  = Kelvin
  real, parameter :: UnitY  = 1.0
  real, parameter :: UnitV  = Centimeter/Second
  real, parameter :: conv_J = Gram/Second**2/Centimeter
  real, parameter :: conv_H = Gram/Second**3
  real, parameter :: conv_e = MeV
#endif


!!$#ifndef FLASH_GRID_AMREX
!!$  blockID = tileDesc%id    ! Maybe useful for debugging
!!$#endif
  blkLimits = tileDesc%limits
  level     = tileDesc%level

  !! allocate all needed space
  allocate(xCenter(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)))
  allocate(xLeft  (blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)))
  allocate(xRight (blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)))
  allocate(yCenter(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)))
  allocate(yLeft  (blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)))
  allocate(yRight (blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)))
  allocate(zCenter(blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
  allocate(zLeft  (blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
  allocate(zRight (blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

  xCenter(:) = 0.e0
  yCenter(:) = 0.e0
  zCenter(:) = 0.e0

  call Grid_getDeltas(level, delta)
  dx = delta(IAXIS)
  dy = delta(JAXIS)
  dz = delta(KAXIS)

! call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getGeometry(meshGeom)

  call Grid_getCellCoords(IAXIS,CENTER,     level, blkLimits(LOW,:), blkLimits(HIGH,:), xCenter)
  call Grid_getCellCoords(IAXIS,LEFT_EDGE,  level, blkLimits(LOW,:), blkLimits(HIGH,:), xLeft  )
  call Grid_getCellCoords(IAXIS,RIGHT_EDGE, level, blkLimits(LOW,:), blkLimits(HIGH,:), xRight )

  call Grid_getCellCoords(JAXIS,CENTER,     level, blkLimits(LOW,:), blkLimits(HIGH,:), yCenter)
  call Grid_getCellCoords(JAXIS,LEFT_EDGE,  level, blkLimits(LOW,:), blkLimits(HIGH,:), yLeft  )
  call Grid_getCellCoords(JAXIS,RIGHT_EDGE, level, blkLimits(LOW,:), blkLimits(HIGH,:), yRight )

  call Grid_getCellCoords(KAXIS,CENTER,     level, blkLimits(LOW,:), blkLimits(HIGH,:), zCenter)
  call Grid_getCellCoords(KAXIS,LEFT_EDGE,  level, blkLimits(LOW,:), blkLimits(HIGH,:), zLeft  )
  call Grid_getCellCoords(KAXIS,RIGHT_EDGE, level, blkLimits(LOW,:), blkLimits(HIGH,:), zRight )

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           ! yCenter and zCenter should be 0.0 when NDIM does not include them
           if ( meshGeom == SPHERICAL ) then

              radCenter = xCenter(i)
              radCenterVol = (1.0/3.0)*( xLeft(i)**3 &
                 + 0.5*(xRight(i)-xLeft(i))*(3.0*xLeft(i)*xRight(i) + (xRight(i)-xLeft(i))**2) )

              thtCenter = yCenter(j)
              thtCenterVol = 1.0 - 0.5 * ( cos(yLeft(j)) + cos(yRight(j)) ) 

              phiCenter = zCenter(k)
              phiCenterVol = phiCenter

           else if ( meshGeom == CYLINDRICAL ) then

              radCenter = sqrt(xCenter(i)**2 + yCenter(j)**2)
              radCenterVol = (1.0/3.0) * radCenter**3

              if ( yCenter(j) /= 0.0 ) then
                 thtCenter = atan( xCenter(i) / yCenter(j) )
              else
                 thtCenter = 0.5 * PI
              end if
              if ( thtCenter < 0.0 ) then
                 thtCenter = thtCenter + PI
              end if
              thtCenterVol = 1.0 - cos(thtCenter)

              phiCenter = zCenter(k)
              phiCenterVol = phiCenter

           else if ( meshGeom == CARTESIAN ) then

              radCenter = sqrt(xCenter(i)**2 + yCenter(j)**2 + zCenter(k)**2)
              radCenterVol = (1.0/3.0) * radCenter**3

              thtCenter = acos( zCenter(k) / radCenter )
              thtCenterVol = 1.0 - cos(thtCenter)

              phiCenter = mod( atan2( yCenter(j), xCenter(i) ), 2.0*PI )
              phiCenterVol = phiCenter

           else

              call Driver_abort("Geometry not supported")

           end if

           ! Use chimera data if within user-defined maximum radius
           if (radCenter <= sim_max_r .or. len_trim(progenitor_model_file) == 0) then

              ! Make sure we can interpolate to this point
              if (radCenter <= x_e_chim(nx_chim+1)) then

                 if ( NDIM == 1 ) then
                    ! Interpolate density in volume
                    call interp1d_linear(volx_c_chim, rho_c_chim(:,1,1), radCenterVol, dens_interp)
                    ! Interpolate other quantities in radius
                    call interp1d_linear(x_c_chim, et_c_chim(:,1,1), radCenter, eint_interp)
                    call interp1d_linear(x_c_chim, p_c_chim(:,1,1), radCenter, pres_interp)
                    call interp1d_linear(x_c_chim, t_c_chim(:,1,1), radCenter, temp_interp)
                    call interp1d_linear(x_c_chim, u_c_chim(:,1,1), radCenter, velx_interp)
                    call interp1d_linear(x_c_chim, v_c_chim(:,1,1), radCenter, vely_interp)
                    call interp1d_linear(x_c_chim, w_c_chim(:,1,1), radCenter, velz_interp)
                    do n = 1, NSPECIES
                       call interp1d_linear(x_c_chim, xn_c_chim(n,:,1,1), radCenter, spec_interp(n))
                    end do
#if defined (YE_MSCALAR)
                    call interp1d_linear(x_c_chim, ye_c_chim(:,1,1), radCenter, ye_interp)
#endif
                 else if ( NDIM == 2 ) then
                    ! Interpolate density in volume
                    call interp2d_linear(volx_c_chim, voly_c_chim, rho_c_chim(:,:,1), &
                       &                 radCenterVol, thtCenterVol, dens_interp)
                    ! Interpolate other quantities in radius
                    call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), et_c_chim(:,:,1), &
                       &                 radCenter, thtCenter, eint_interp)
                    call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), p_c_chim(:,:,1), &
                       &                 radCenter, thtCenter, pres_interp)
                    call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), t_c_chim(:,:,1), &
                       &                 radCenter, thtCenter, temp_interp)
                    call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), u_c_chim(:,:,1), &
                       &                 radCenter, thtCenter, velx_interp)
                    call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), v_c_chim(:,:,1), &
                       &                 radCenter, thtCenter, vely_interp)
                    call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), w_c_chim(:,:,1), &
                       &                 radCenter, thtCenter, velz_interp)
                    do n = 1, NSPECIES
                       call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), xn_c_chim(n,:,:,1), &
                          &                 radCenter, thtCenter, spec_interp(n))
                    end do
#if defined (YE_MSCALAR)
                    call interp2d_linear(x_c_chim, y_c_chim(1:ny_chim), ye_c_chim(:,:,1), &
                       &                 radCenter, thtCenter, ye_interp)
#endif
                 else if ( NDIM == 3 ) then
                    ! Interpolate density in volume
                    call interp3d_linear(volx_c_chim, voly_c_chim, volz_c_chim, rho_c_chim(:,:,:), &
                       &                 radCenterVol, thtCenterVol, phiCenterVol, dens_interp)
                    ! Interpolate other quantities in radius
                    call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), et_c_chim(:,:,:), &
                       &                 radCenter, thtCenter, phiCenter, eint_interp)
                    call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), p_c_chim(:,:,:), &
                       &                 radCenter, thtCenter, phiCenter, pres_interp)
                    call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), t_c_chim(:,:,:), &
                       &                 radCenter, thtCenter, phiCenter, temp_interp)
                    call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), u_c_chim(:,:,:), &
                       &                 radCenter, thtCenter, phiCenter, velx_interp)
                    call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), v_c_chim(:,:,:), &
                       &                 radCenter, thtCenter, phiCenter, vely_interp)
                    call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), w_c_chim(:,:,:), &
                       &                 radCenter, thtCenter, phiCenter, velz_interp)
                    do n = 1, NSPECIES
                       call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), xn_c_chim(n,:,:,:), &
                          &                 radCenter, thtCenter, phiCenter, spec_interp(n))
                    end do
#if defined (YE_MSCALAR)
                    call interp3d_linear(x_c_chim, y_c_chim(1:ny_chim), z_c_chim(1:nz_chim), ye_c_chim(:,:,:), &
                       &                 radCenter, thtCenter, phiCenter, ye_interp)
#endif
                 end if
                 solnData(DENS_VAR,i,j,k) = dens_interp
                 solnData(EINT_VAR,i,j,k) = eint_interp
                 solnData(PRES_VAR,i,j,k) = pres_interp
                 solnData(TEMP_VAR,i,j,k) = temp_interp
                 solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = spec_interp(:)
#if defined (YE_MSCALAR)
                 solnData(YE_MSCALAR,i,j,k) = ye_interp
#endif
                 if ( meshGeom == CARTESIAN ) then
                    solnData(VELX_VAR,i,j,k) = + velx_interp * sin( thtCenter ) * cos( phiCenter ) &
                       &                       + vely_interp * cos( thtCenter ) * cos( phiCenter ) &
                       &                       - velz_interp * sin( thtCenter ) * sin( phiCenter )
                    solnData(VELY_VAR,i,j,k) = + velx_interp * sin( thtCenter ) * sin( phiCenter ) &
                       &                       + vely_interp * cos( thtCenter ) * sin( phiCenter ) &
                       &                       + velz_interp * sin( thtCenter ) * cos( phiCenter )
                    solnData(VELZ_VAR,i,j,k) = + velx_interp * cos( thtCenter ) &
                       &                       - vely_interp * sin( thtCenter )
                 else if ( meshGeom == CYLINDRICAL ) then
                    solnData(VELX_VAR,i,j,k) = + velx_interp * sin( thtCenter ) &
                       &                       + vely_interp * cos( thtCenter )
                    solnData(VELY_VAR,i,j,k) = + velx_interp * cos( thtCenter ) &
                       &                       - vely_interp * sin( thtCenter )
                    solnData(VELZ_VAR,i,j,k) = velz_interp
                 else if ( meshGeom == SPHERICAL ) then
                    solnData(VELX_VAR,i,j,k) = velx_interp
                    solnData(VELY_VAR,i,j,k) = vely_interp
                    solnData(VELZ_VAR,i,j,k) = velz_interp
                 end if
              else

                 call Driver_abort("Beyond the Chimera data")

              end if

           else if (radCenter <= xzn(n1d_total)) then ! use progenitor data

              velx_interp = 0.0
              vely_interp = 0.0
              velz_interp = 0.0
              do ivar = 1, NUNK_VARS
                 if ( ivar == DENS_VAR ) then
                    call interp1d_linear(volxzn(1:n1d_total), model_1d(1:n1d_total,ivar), radCenterVol, var_interp)
                 else
                    call interp1d_linear(xzn(1:n1d_total), model_1d(1:n1d_total,ivar), radCenter, var_interp)
                 end if
                 ! hold on to velocities until later when we can convert them to the proper geometry
                 if ( ivar == VELX_VAR ) then
                    velx_interp = var_interp
                 else if ( ivar == VELY_VAR ) then
                    vely_interp = var_interp
                 else if ( ivar == VELZ_VAR ) then
                    velz_interp = var_interp
                 else
                    solnData(ivar,i,j,k) = var_interp
                 end if
              end do

              if ( meshGeom == CARTESIAN ) then
                 solnData(VELX_VAR,i,j,k) = + velx_interp * sin( thtCenter ) * cos( phiCenter ) &
                    &                       + vely_interp * cos( thtCenter ) * cos( phiCenter ) &
                    &                       - velz_interp * sin( thtCenter ) * sin( phiCenter )
                 solnData(VELY_VAR,i,j,k) = + velx_interp * sin( thtCenter ) * sin( phiCenter ) &
                    &                       + vely_interp * cos( thtCenter ) * sin( phiCenter ) &
                    &                       + velz_interp * sin( thtCenter ) * cos( phiCenter )
                 solnData(VELZ_VAR,i,j,k) = + velx_interp * cos( thtCenter ) &
                    &                       - vely_interp * sin( thtCenter )
              else if ( meshGeom == CYLINDRICAL ) then
                 solnData(VELX_VAR,i,j,k) = + velx_interp * sin( thtCenter ) &
                    &                       + vely_interp * cos( thtCenter )
                 solnData(VELY_VAR,i,j,k) = + velx_interp * cos( thtCenter ) &
                    &                       - vely_interp * sin( thtCenter )
                 solnData(VELZ_VAR,i,j,k) = velz_interp
              else if ( meshGeom == SPHERICAL ) then
                 solnData(VELX_VAR,i,j,k) = velx_interp
                 solnData(VELY_VAR,i,j,k) = vely_interp
                 solnData(VELZ_VAR,i,j,k) = velz_interp
              end if

           else

              call Driver_abort("Beyond the progenitor data")

           endif

#if NSPECIES > 0
#if defined (YE_MSCALAR)
           call renorm_mass_ye(solnData(YE_MSCALAR,i,j,k),solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
#else
           solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = &
              & max(sim_smallx,min(1.0,solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)))
           suminv = 1.0 / sum(solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
           solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = &
              & max(sim_smallx,min(1.0,suminv*solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)))
#endif

#endif
#if defined (SUMY_MSCALAR)
           call Eos_getAbarZbar(solnData(:,i,j,k),abar,zbar,sumY)
           solnData(SUMY_MSCALAR,i,j,k) = sumY
#endif

        end do
     end do
  end do

#if defined (THORNADO)
  !print*, 'Initializing radiation fields with thornado.'

  lo(1:MDIM) = tileDesc%limits(LOW,1:MDIM)
  hi(1:MDIM) = tileDesc%limits(HIGH,1:MDIM)

  nX = 1
  swX = 0
  xL = 0.0
  if ( meshGeom == CARTESIAN ) then
     xR = 1.0
  else if ( meshGeom == CYLINDRICAL ) then
     xR = [ 1.0, 1.0, 2.0*PI ]
  else if ( meshGeom == SPHERICAL ) then
     xR = [ 1.0, PI, 2.0*PI ]
  else
     call Driver_abort("Geometry not supported")
  end if

  nX(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
  swX(1:NDIM) = 2
  u_lo(2:4) = 1 - swX
  u_hi(2:4) = nX + swX
  u_lo(1)   = 1 - THORNADO_SWE
  u_hi(1)   = THORNADO_NE + THORNADO_SWE

  ! convert cm to m for Thornado
  xL(1) = xL(1) * conv_x
  xR(1) = xR(1) * conv_x
  if ( meshGeom /= SPHERICAL ) then
     xL(2) = xL(2) * conv_x
     xR(2) = xR(2) * conv_x
  end if
  if ( meshGeom == CARTESIAN ) then
     xL(3) = xL(3) * conv_x
     xR(3) = xR(3) * conv_x
  end if

  if ( meshGeom == SPHERICAL ) then

     call InitThornado_Patch &
          (nX, swX, xL, xR, THORNADO_NSPECIES, "spherical" )

  else if ( meshGeom == CYLINDRICAL ) then

     call InitThornado_Patch &
          (nX, swX, xL, xR, THORNADO_NSPECIES, "cylindrical" )

  else

     call InitThornado_Patch &
          (nX, swX, xL, xR, THORNADO_NSPECIES, "cartesian" )

  end if

  do iX3 = 1, nX(3)
     do iX2 = 1, nX(2)
        do iX1 = 1, nX(1)

           !print*, 'In cell: ', iX1, iX2, iX3

           i = lo(IAXIS) + THORNADO_NNODESX*(iX1-1)
           j = lo(JAXIS) + THORNADO_NNODESX*(iX2-1)
           k = lo(KAXIS) + THORNADO_NNODESX*(iX3-1)

           do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS ; do iE = u_lo(1), u_hi(1)

              ioff = THORNADO_BEGIN &
                 + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE) &
                           *THORNADO_NMOMENTS) &
                 + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                 + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

              do iNode = 1, THORNADO_RAD_NDOF

                 iNodeE  = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                 iNodeX  = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                 iNodeX1 = mod((iNodeX-1)                    ,THORNADO_NNODESX) + 1
                 ii      = iNodeX1 + i - 1

                 iNodeX2 = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + 1
                 jj      = iNodeX2 + j - 1

                 iNodeX3 = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + 1
                 kk      = iNodeX3 + k - 1
                 
                 ! calculate the indices
                 ivar = ioff + iNodeE - 1

                 !print*, 'ivar: ', ivar

                 ! calculate actual positions of the nodes used for the gaussian quadrature
                 xnode = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                 ynode = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
                 znode = NodeCoordinate( MeshX(3), iX3, iNodeX3 )
                 enode = NodeCoordinate( MeshE,    iE,  iNodeE  )

                 xnode = xnode / conv_x
                 if ( meshGeom /= SPHERICAL ) ynode = ynode / conv_x
                 if ( meshGeom == CARTESIAN ) znode = znode / conv_x

                 D = 1.0e-20
                 I1 = 0.0e0
                 I2 = 0.0e0
                 I3 = 0.0e0

                 if( meshGeom == SPHERICAL ) then

                    velr = solnData(VELX_VAR,ii,jj,kk) * UnitV
                    velt = 0.0
                    velp = 0.0

                    Gmdd11 = uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11)
                    Gmdd22 = uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22)
                    Gmdd33 = uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33)

                    call ComputeConserved_TwoMoment &
                       ( D, I1, I2, I3, &
                         CN, G1, G2, G3, &
                         velr, velt, velp, &
                         Gmdd11, Gmdd22, Gmdd33 )

                 end if

                 if( iCR == iCR_N ) solnData(ivar,ii,jj,kk) = CN

                 if( iCR == iCR_G1 ) solnData(ivar,ii,jj,kk) = G1

                 if( iCR == iCR_G2 ) solnData(ivar,ii,jj,kk) = G2

                 if( iCR == iCR_G3 ) solnData(ivar,ii,jj,kk) = G3

              end do

           end do ; end do ; end do ;

        end do
     end do
  end do

  call FreeThornado_Patch()
#endif

! call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yLeft)
  deallocate(yRight)
  deallocate(yCenter)
  deallocate(zLeft)
  deallocate(zRight)
  deallocate(zCenter)

  call Eos_wrapped(MODE_DENS_TEMP,tileDesc%limits,solnData)

  return
end subroutine Simulation_initBlock

#ifdef FLASH_MULTISPECIES
subroutine renorm_mass_ye(ye,xx)
! given an electron fraction (ye) and vector of mass fractions (xx),
! renormalize mass fractions such that they sum to 1 and reproduce the electron fraction
   use Multispecies_interface, ONLY : Multispecies_getProperty
   use Simulation_data, ONLY : sim_smallx
   implicit none

#include "Simulation.h"
#include "Multispecies.h"

   real, intent(in) :: ye
   real, intent(inout) :: xx(NSPECIES)
   real :: ny, zy, zny, zzy, aa, zz, nn, ainv, alpha, beta
   integer :: i

   ny = 0.0 ; zy = 0.0 ; zny = 0.0 ; zzy = 0.0
   do i = 1, NSPECIES
      xx(i) = max(sim_smallx,min(1.0,xx(i)))
      call Multispecies_getProperty(SPECIES_BEGIN + (i-1), A, aa)
      call Multispecies_getProperty(SPECIES_BEGIN + (i-1), Z, zz)
      call Multispecies_getProperty(SPECIES_BEGIN + (i-1), N, nn)
      ainv = 1.0 / aa
      ny = ny + xx(i) * nn * ainv
      zy = zy + xx(i) * zz * ainv
      zny = zny + xx(i) * zz * nn * ainv * ainv
      zzy = zzy + xx(i) * zz * zz * ainv * ainv
   end do
   beta = (ye*ny - zny) / (ny*zzy - zy*zny)
   alpha = (1.0 - beta*zy) / ny
   do i = 1, NSPECIES
      call Multispecies_getProperty(SPECIES_BEGIN + (i-1), Z, zz)
      call Multispecies_getProperty(SPECIES_BEGIN + (i-1), N, nn)
      call Multispecies_getProperty(SPECIES_BEGIN + (i-1), A, aa)
      xx(i) = xx(i) * (alpha*nn + beta*zz) / aa
      xx(i) = max(sim_smallx,min(1.0,xx(i)))
   end do

   return
end subroutine renorm_mass_ye
#endif
