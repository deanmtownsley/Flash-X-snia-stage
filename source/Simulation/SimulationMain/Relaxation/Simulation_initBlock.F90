!!****if* source/Simulation/SimulationMain/Relaxation/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: blockDesc  )
!!
!!
!! DESCRIPTION
!!
!!   Initialize solution data in one block for a deleptonization wave
!!
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!! PARAMETERS
!!
!!  
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(solnData, tileDesc)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abort
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_coordTransfm
  use Eos_interface, ONLY : Eos_multiDim

  use GeometryFieldsModule, ONLY : uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  use RadiationFieldsModule, ONLY : iCR_N, iCR_G1, iCR_G2, iCR_G3
  use UnitsModule, ONLY : Centimeter, Gram, Second, MeV, Kelvin, BoltzmannConstant
  use MeshModule, ONLY : NodeCoordinate, MeshE, MeshX
  use NeutrinoOpacitiesComputationModule, ONLY : ComputeEquilibriumDistributions_Point
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"

#if defined(THORNADO_ORDER_V)
  use TwoMoment_UtilitiesModule, ONLY : ComputeConserved_TwoMoment
#endif

  use rt_data, ONLY : rt_UpperBry1

  implicit none

  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in)     :: tileDesc

  real, allocatable, dimension(:) :: xCenter
  real, allocatable, dimension(:) :: yCenter
  real, allocatable, dimension(:) :: zCenter
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real :: rcc, tcc, pcc, velr, velt, velp

  integer, dimension(1:MDIM)   :: lo, hi
  integer, dimension(1:MDIM+1) :: u_lo, u_hi
  integer :: i, j, k, n, ii, jj, kk, ic, jc, kc
  integer :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3
  integer :: iS, iCR, iE, iNode, iNodeX, iNodeE, ioff, ivar
  integer :: nX(3), swX(3)
  real :: xL(3), xR(3)
  real :: xnode, ynode, znode, enode
  real :: E, kT, f_E
  real :: Nnu, Gnu1, Gnu2, Gnu3
  real :: Dnu, Inu1, Inu2, Inu3
  real :: dens_buffer, temp_buffer, ye_buffer, velr_buffer
  real :: dens_i, temp_i, velx_i, vely_i, velz_i, ye_i, mu_i

  real, parameter :: conv_x = Centimeter
  real, parameter :: UnitD  = Gram / Centimeter**3
  real, parameter :: UnitT  = Kelvin
  real, parameter :: UnitY  = 1.0
  real, parameter :: UnitV  = Centimeter/Second
  real, parameter :: conv_J = Gram/Second**2/Centimeter
  real, parameter :: conv_H = Gram/Second**3
  real, parameter :: conv_e = MeV

  ! get dimensions/limits and coordinates
  lo(1:MDIM) = tileDesc%limits(LOW,1:MDIM)
  hi(1:MDIM) = tileDesc%limits(HIGH,1:MDIM)

  !! allocate all needed space
  allocate(xCenter(lo(IAXIS):hi(IAXIS)))
  allocate(yCenter(lo(JAXIS):hi(JAXIS)))
  allocate(zCenter(lo(KAXIS):hi(KAXIS)))

  xCenter = 0.0
  yCenter = 0.0
  zCenter = 0.0

  call Grid_getCellCoords(IAXIS,CENTER,tileDesc%level,lo,hi,xCenter)
  if ( NDIM > 1 ) call Grid_getCellCoords(JAXIS,CENTER,tileDesc%level,lo,hi,yCenter)
  if ( NDIM > 2 ) call Grid_getCellCoords(KAXIS,CENTER,tileDesc%level,lo,hi,zCenter)

  nX = 1
  swX = 0
  xL = 0.0
  if ( sim_geometry == CARTESIAN ) then
     xR = 1.0
  else if ( sim_geometry == CYLINDRICAL ) then
     xR = [ 1.0, 1.0, 2.0*PI ]
  else if ( sim_geometry == SPHERICAL ) then
     xR = [ 1.0, PI, 2.0*PI ]
  else
     call Driver_abort("Geometry not supported")
  end if

  nX(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
  swX(1:NDIM) = NGUARD / THORNADO_NNODESX
  u_lo(2:4) = 1 - swX
  u_hi(2:4) = nX + swX
  u_lo(1)   = 1 - THORNADO_SWE
  u_hi(1)   = THORNADO_NE + THORNADO_SWE

  call tileDesc%boundBox(boundBox)
  xL(1:NDIM) = boundBox(LOW, 1:NDIM)
  xR(1:NDIM) = boundBox(HIGH,1:NDIM)

  ! convert cm to m for Thornado
  xL(1) = xL(1) * conv_x
  xR(1) = xR(1) * conv_x
  if ( sim_geometry /= SPHERICAL ) then
     xL(2) = xL(2) * conv_x
     xR(2) = xR(2) * conv_x
  end if
  if ( sim_geometry == CARTESIAN ) then
     xL(3) = xL(3) * conv_x
     xR(3) = xR(3) * conv_x
  end if

  call InitThornado_Patch &
       (nX, swX, xL, xR, THORNADO_NSPECIES, sim_str_geometry )

  do iX3 = 1, nX(3)
     do iX2 = 1, nX(2)
        do iX1 = 1, nX(1)

           i = lo(IAXIS) + THORNADO_NNODESX*(iX1-1)
           j = lo(JAXIS) + THORNADO_NNODESX*(iX2-1)
           k = lo(KAXIS) + THORNADO_NNODESX*(iX3-1)

           ! Initialize hydro data
           do iNodeX = 1, THORNADO_FLUID_NDOF

              ii = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
              jj = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
              kk = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

              if ( sim_use_model .or. sim_rad_option /= 4 ) then

                 ! compute the spherical polar cell-center coordinates
                 call Grid_coordTransfm &
                    ( xCenter(ii), yCenter(jj), zCenter(kk), &
                      rcc, tcc, pcc, geometryOut = SPHERICAL )

                 if ( sim_use_model ) then
                    call sim_interpProfile &
                       ( rcc, sim_dens_i, sim_temp_i, sim_ye_i, velr )
                 elseif ( sim_rad_option /= 4 ) then
                    call sim_analyticProfile &
                       ( rcc, sim_dens_i, sim_temp_i, sim_ye_i, velr )
                 end if

                 velt = 0.0
                 velp = 0.0
                 call sim_velFromSph &
                    ( rcc, tcc, pcc, velr, velt, velp, &
                      sim_geometry, sim_velx_i, sim_vely_i, sim_velz_i )

              end if

              solnData(DENS_VAR,ii,jj,kk)   = sim_dens_i
              solnData(TEMP_VAR,ii,jj,kk)   = sim_temp_i
              solnData(VELX_VAR,ii,jj,kk)   = sim_velx_i
              solnData(VELY_VAR,ii,jj,kk)   = sim_vely_i
              solnData(VELZ_VAR,ii,jj,kk)   = sim_velz_i
              solnData(YE_MSCALAR,ii,jj,kk) = sim_ye_i

#if NSPECIES>0
              do n = SPECIES_BEGIN,SPECIES_END
                 solnData(n,ii,jj,kk) = sim_xn_i(n)
              end do
#endif

           end do

           ! Initialize neutrino data
           do iS = 1, THORNADO_NSPECIES ; do iE = u_lo(1), u_hi(1)
              do iNode = 1, THORNADO_RAD_NDOF

                 iNodeE  = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                 iNodeX  = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                 iNodeX1 = mod((iNodeX-1)                    ,THORNADO_NNODESX) + 1
                 iNodeX2 = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + 1
                 iNodeX3 = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + 1

                 ii      = iNodeX1 + i - 1
                 jj      = iNodeX2 + j - 1
                 kk      = iNodeX3 + k - 1

                 ! calculate the indices
                 ivar = ioff + iNodeE - 1

                 ! calculate actual positions of the nodes used for the gaussian quadrature
                 xnode = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                 ynode = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
                 znode = NodeCoordinate( MeshX(3), iX3, iNodeX3 )
                 enode = NodeCoordinate( MeshE,    iE,  iNodeE  )

                 xnode = xnode / conv_x
                 if ( sim_geometry /= SPHERICAL ) ynode = ynode / conv_x
                 if ( sim_geometry == CARTESIAN ) znode = znode / conv_x

                 if ( sim_use_model .or. sim_rad_option /= 4 ) then

                    ! compute the spherical polar cell-center coordinates
                    call Grid_coordTransfm &
                       ( xnode, ynode, znode, &
                         rcc, tcc, pcc, geometryOut = SPHERICAL )

                    if ( sim_use_model ) then
                       call sim_interpProfile &
                          ( rcc, sim_dens_i, sim_temp_i, sim_ye_i, velr )
                    elseif( sim_rad_option /= 4 ) then ! analytical approx
                       call sim_analyticProfile &
                          ( rcc, sim_dens_i, sim_temp_i, sim_ye_i, velr )
                    end if

                    velt = 0.0
                    velp = 0.0
                    call sim_velFromSph &
                       ( rcc, tcc, pcc, velr, velt, velp, &
                         sim_geometry, sim_velx_i, sim_vely_i, sim_velz_i )

                 end if

                 dens_i = sim_dens_i * UnitD
                 temp_i = sim_temp_i * UnitT
                 velx_i = sim_velx_i * UnitV
                 vely_i = sim_vely_i * UnitV
                 velz_i = sim_velz_i * UnitV
                 ye_i   = sim_ye_i   * UnitY
                 mu_i   = sim_mu_i

                 Nnu  = 0.0
                 Gnu1 = 0.0
                 Gnu2 = 0.0
                 Gnu3 = 0.0

                 select case ( sim_rad_option )
                   case ( -1 ) ! Zero Distribution

                     Nnu = 1.0e-20
                     Gnu1 = 0.0e0

                   case ( 0 ) ! FD with sim_rintSwitch

                     if( rcc >= sim_rintSwitch ) then

                       call sim_interpProfile &
                            ( sim_rintSwitch, dens_buffer, temp_buffer, ye_buffer, velr_buffer)
                       dens_buffer = dens_buffer * UnitD
                       temp_buffer = temp_buffer * UnitT
                       ye_buffer   = ye_buffer   * UnitY
                       call ComputeEquilibriumDistributions_Point &
                            ( enode, dens_buffer, temp_buffer, ye_buffer, Nnu, iS )
                       Nnu = Nnu * ( sim_rintSwitch / rcc )**2

                     else
                       call ComputeEquilibriumDistributions_Point &
                            ( enode, dens_i, temp_i, ye_i, Nnu, iS )
                     end if

                     Gnu1 = 0.0e0

                   case ( 1 ) ! approximate neutrino spher

                     call sim_interp_rad( iE, iNodeE, rcc, iS, Nnu, Gnu1 )

                   case ( 2 ) ! Chimera Profile radiation

                     call sim_interpProfile_rad( enode/MeV, rcc, iS, Nnu, Gnu1 )

                   case ( 3 ) ! BoltzTran Profile radiation

                     call sim_interpProfile_rad( enode/MeV, rcc, iS, Nnu, Gnu1 )

                   case ( 4 ) ! Thornado Spectrum radiation

                     kT = BoltzmannConstant * temp_i
                     E  = enode

                     f_E = MAX( 0.99 * EXP( - ( E - 2.0 * kT )**2 &
                                               / ( 2.0*(1.0e1*MeV)**2 ) ), 1.0e-99 )
#if   defined(THORNADO_ORDER_1)
                     Nnu  = f_E
                     Gnu1 = 0.0
#elif defined(THORNADO_ORDER_V)
                     Dnu  = f_E * 0.50 * ( 1.0 - mu_i    )
                     Inu1 = f_E * 0.25 * ( 1.0 - mu_i**2 )
                     Inu2 = 0.0
                     Inu3 = 0.0
                     CALL ComputeConserved_TwoMoment &
                          ( Dnu, Inu1, Inu2, Inu3, &
                            Nnu, Gnu1, Gnu2, Gnu3, &
                            velx_i, vely_i, velz_i, &
                            uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                            uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                            uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )
#endif

                 end select

                 ! realizability check
                 ! --- non-negative zeroth moment
                 if( Nnu > rt_UpperBry1 ) then
                   !!$write(*,'(A,2ES12.3)') &
                   !!$  ' Triggered realizability check in Simulation_initBlock.F90 J=', Nnu
                   Nnu = max( 0.0e0, min(Nnu,rt_UpperBry1) )
                   if( Gnu1 < 0.0 ) then
                     Gnu1 = -(1.0d0-Nnu)*Nnu
                   else
                     Gnu1 = (1.0d0-Nnu)*Nnu
                   end if
                 end if
                 ! --- (1-J)*J >= abs(H)
                 if( abs(Gnu1) > (1.0d0-Nnu)*Nnu ) then
                   !!$write(*,'(A,2ES12.3)') &
                   !!$  ' Triggered realizability check in Simulation_initBlock.F90 gamma',&
                   !!$  Nnu, Gnu1
                   !!$write(*,'(A,2ES12.3)') '       @ (E,X)', enode/MeV, rcc
                   if( Gnu1 < 0.0 ) then
                     Gnu1 = -(1.0d0-Nnu)*Nnu
                   else
                     Gnu1 = (1.0d0-Nnu)*Nnu
                   end if
                   if( abs(Gnu1) > (1.0d0-Nnu)*Nnu )&
                      call Driver_abort("not preserve realizability in initializing radiation")
                 end if

                 ! J moment, iCR = 1
                 iCR  = iCR_N
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Nnu

                 ! H_x moment, iCR = 2
                 iCR  = iCR_G1
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Gnu1

                 ! H_y moment, iCR = 3
                 iCR  = iCR_G2
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Gnu2

                 ! H_z moment, iCR = 4
                 iCR  = iCR_G3
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Gnu3

              end do
           end do ; end do

        end do
     end do
  end do

  call FreeThornado_Patch()

  ! cleanup
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)

  call Eos_multiDim(MODE_DENS_TEMP,tileDesc%limits,solnData)

  return

contains

  subroutine sim_velFromSph( r, theta, phi, vr, vtheta, vphi, geometryOut, vx, vy, vz )

     implicit none

#include "constants.h"

     real, intent(in) :: r, theta, phi, vr, vtheta, vphi
     integer, intent(in) :: geometryOut
     real, intent(out) :: vx, vy, vz

     if ( geometryOut == CARTESIAN ) then
        vx = + vr     * sin( theta ) * cos( phi ) &
           & + vtheta * cos( theta ) * cos( phi ) &
           & - vphi   * sin( theta ) * sin( phi )
        vy = + vr     * sin( theta ) * sin( phi ) &
           & + vtheta * cos( theta ) * sin( phi ) &
           & + vphi   * sin( theta ) * cos( phi )
        vz = + vr     * cos( theta ) &
           & - vtheta * sin( theta )
     else if ( geometryOut == CYLINDRICAL ) then
        vx = + vr     * sin( theta ) &
           & + vtheta * cos( theta )
        vy = + vr     * cos( theta ) &
           & - vtheta * sin( theta )
        vz = vphi
     else if ( geometryOut == SPHERICAL ) then
        vx = vr
        vy = vtheta
        vz = vphi
     end if

     return
  end subroutine sim_velFromSph

  subroutine sim_interpProfile( rcc, dens, temp, ye, velr )

     use Simulation_data
     use Driver_interface, ONLY : Driver_abort
     use ut_interpolationInterface, ONLY : ut_hunt, ut_polint

     implicit none

     real, intent(in) :: rcc
     real, intent(out) :: dens, temp, ye, velr

     integer, parameter :: order = 1

     real :: err
     integer :: idx

     if ( rcc > xzn(n1d_max) ) then
        write(*,'(A,ES15.3)') ' location not covered by profile : r = ', rcc
        call Driver_abort('Profile mismatching.') 
     else

        call ut_hunt(xzn, n1d_max, rcc, idx)

        ! if we are below the minimum in radius, then just set the index to the first value
        idx = max(1,min(idx,n1d_max-1))

        ! AH: do we want to interpolate in volume-coordinates for density?
        ! AH: and interpolate rho*var for other values?
        call ut_polint(xzn(idx), model_1d(idx,DENS_VAR),   order+1, rcc, dens, err)
        call ut_polint(xzn(idx), model_1d(idx,TEMP_VAR),   order+1, rcc, temp, err)
        call ut_polint(xzn(idx), model_1d(idx,YE_MSCALAR), order+1, rcc, ye,   err)

        call ut_polint(xzn(idx), model_1d(idx,VELX_VAR),   order+1, rcc, velr, err)

     end if

     return
  end subroutine sim_interpProfile

  subroutine sim_interpProfile_rad( e_MeV, rcc, iS, Nnu, Gnu1 )

     use Simulation_data, ONLY : xznrad, ezn, model_1d_rad, n1d_nE
     use Driver_interface, ONLY : Driver_abort
     use ut_interpolationInterface, ONLY : ut_hunt

     implicit none

     integer, intent(in) :: iS
     real,    intent(in) :: e_MeV, rcc
     real,    intent(out) :: Nnu, Gnu1

     integer, parameter :: order = 1

     real    :: dx, de, p00, p10, p01, p11, gamma
     integer :: idx, ide, nX

     nX = size(xznrad)

     !! not able to handle region that out of the profile
     if ( rcc > xznrad(nX) ) then
        write(*,'(A,ES15.3)') ' location not covered by profile : r = ', rcc
        call Driver_abort('Profile mismatching.')
     end if
     !! not able to handle spectrum that exceeds the profile's coverage
     if ( e_MeV > ezn(size(ezn)) ) then
        write(*,'(A,ES15.3)') ' energy not covered by profile : e = ', e_MeV
        call Driver_abort('Profile mismatching.')
     end if
     !! if rcc < the inner most point in the profile
     !! use constant extrapolation
     if ( rcc < xznrad(1) ) then
        !!$print*, 'rcc < xznrad(1) : ', rcc, ' < ', xznrad(1)
        idx = 1
        call ut_hunt(ezn, n1d_nE, e_MeV, ide)

        ! use log interpolation
        de = log10(e_MeV/ezn(ide)) / log10(ezn(ide+1)/ezn(ide))
        dx = 0.0e0

        associate( Table1 => model_1d_rad(:,:,1,iS), &
                   Table2 => model_1d_rad(:,:,2,iS) )

        p00 = Table1( ide  , idx )
        P10 = Table1( ide+1, idx )
        p01 = Table1( ide  , idx )
        p11 = Table1( ide+1, idx )
 
        Nnu = BiLinear( p00, p10, p01, p11, de, dx )

        gamma = ( 1.0e0 - Nnu ) * Nnu

         p00 = Table2( ide  , idx   )
         P10 = Table2( ide+1, idx   )
         p01 = Table2( ide  , idx   )
         p11 = Table2( ide+1, idx   )

         Gnu1 = BiLinear( p00, p10, p01, p11, de, dx )
         if( abs(Gnu1) > gamma ) then
           if( Gnu1 < 0.0 ) then
             Gnu1 = - gamma
           else
             Gnu1 = gamma
           end if
         end if
        end associate
     end if

     if ( (rcc >= xznrad(1)) .and. ( rcc<= xznrad(nX)) .and. (e_MeV <= ezn(size(ezn))) ) then

        call ut_hunt(xznrad, nX, rcc, idx)
        ! if we are below the minimum in radius, then just set the index to the first value
        idx = max(1,min(idx,nX-1))

        if( e_MeV > ezn(size(ezn)) ) then
          ! if we are above the maximum in energy, then just set the index to the second last value
          ide = n1d_nE - 1
        else
          call ut_hunt(ezn, n1d_nE, e_MeV, ide)
          ! if we are below the minimum in energy, then just set the index to the first value
          ide = max(1,min(ide,n1d_nE-1))
        end if

        ! use log interpolation
        de = log10(e_MeV/ezn(ide)) / log10(ezn(ide+1)/ezn(ide))
        dx = log10(rcc  /xznrad(idx)) / log10(xznrad(idx+1)/xznrad(idx))

        associate( Table1 => model_1d_rad(:,:,1,iS), &
                   Table2 => model_1d_rad(:,:,2,iS) )

          p00 = Table1( ide  , idx   )
          P10 = Table1( ide+1, idx   )
          p01 = Table1( ide  , idx+1 )
          p11 = Table1( ide+1, idx+1 )

          Nnu = BiLinear( p00, p10, p01, p11, de, dx )
          if( Nnu < 0.0e0 .or. Nnu > 1.0e0 ) then
            !!$write(*,'(A,ES12.3,A,2I4,A,2ES12.3)') 'J (Nnu)=', Nnu, &
            !!$  ' hit the realizability bry [ide, idx] = [', &
            !!$  ide, idx,'] => [MeV,cm] =', e_MeV, rcc
            Nnu = min(p00, p10, p01, p11)
            write(*,'(A25,2ES15.3)') 'The four neighbors (e,x):', p00, p01
            write(*,'(A30,A15)') '', 'Nnu'
            write(*,'(A25,2ES15.3)') '', p10, p11
            Nnu = min(min(p00,p10),min(p01,p11))
            if( Nnu < 0.0 )then
              write(*,'(A,ES15.6)') 'Nnu = ', Nnu
              call Driver_abort('Fail to enforce J > 0')
            end if
          end if
          !!$if( Nnu > 1.0 ) then
          !!$  !write(*,'(A20,ES25.13,A,2I4)') 'Nnu exceed one !', Nnu-1.0d0, ' @ (indexE,indexX)', ide, idx
          !!$  !write(*,'(A25,2ES15.3)') 'The four neighbors (e,x):', p00, p01
          !!$  !write(*,'(A30,A15)') '', 'Nnu'
          !!$  !write(*,'(A25,2ES15.3)') '', p10, p11
          !!$  !write(*,'(A)') 'Enforce min(Nnu,1.0d0)'
          !!$  Nnu = min(Nnu,1.0e0)
          if( Nnu > 1.0e0 )then
            write(*,'(A,ES15.6)') 'Nnu = ', Nnu
            call Driver_abort('Fail to enforce J < 1.0')
          end if

          gamma = ( 1.0e0 - Nnu ) * Nnu

          p00 = Table2( ide  , idx   )
          P10 = Table2( ide+1, idx   )
          p01 = Table2( ide  , idx+1 )
          p11 = Table2( ide+1, idx+1 )

          Gnu1 = BiLinear( p00, p10, p01, p11, de, dx )
          if( abs(Gnu1) > gamma ) then
            if( Gnu1 < 0.0 ) then
              Gnu1 = - gamma
            else
              Gnu1 = gamma
            end if
          end if

        end associate

     end if

     return
  end subroutine sim_interpProfile_rad

  subroutine sim_interp_rad( iE, iNodeE, rcc, iS, Nnu, Gnu1 )

    use Simulation_data, ONLY: n1d_max, xzn, D_Nu_P, I1_Nu_P
    use ut_interpolationInterface, ONLY : ut_hunt

#include "Simulation.h"

    integer, intent(in)  :: iE, iNodeE, iS
    real,    intent(in)  :: rcc
    real,    intent(out) :: Nnu, Gnu1

    integer, parameter :: order = 1
    integer :: i_nodeE, idx
    real    :: err

    i_nodeE = (iE-1) * THORNADO_NNODESE + iNodeE

    if ( rcc > xzn(n1d_max) ) then
        write(*,'(A,ES15.3)') ' location not covered by profile : r = ', rcc, 'cm'
        call Driver_abort('Profile mismatching.')
    end if

    call ut_hunt(xzn, n1d_max, rcc, idx)
    idx = max( 1, idx )

    Nnu = Interpolate1D_Linear &
              ( rcc, xzn(idx), xzn(idx+1), &
                D_Nu_P (i_nodeE,idx,iS), D_Nu_P (i_nodeE,idx+1,iS) )
    Gnu1 = Interpolate1D_Linear &
              ( rcc, xzn(idx), xzn(idx+1), &
                I1_Nu_P (i_nodeE,idx,iS), I1_Nu_P (i_nodeE,idx+1,iS) )    

    return 
 
  end subroutine sim_interp_rad

  real function BiLinear &
    ( p00, p10, p01, p11, dX1, dX2 )

    real, intent(in) :: &
      p00, p10, p01, p11, dX1, dX2

    real :: ddX1, ddX2

    ddX1 = 1.0e0 - dX1
    ddX2 = 1.0e0 - dX2

    BiLinear =  ddX1 * ( ddX2 * p00 + dX2 * p01 ) &
               + dX1 * ( ddX2 * p10 + dX2 * p11 )

    return
  end function BiLinear

  subroutine sim_analyticProfile( rcc, dens, temp, ye, velr )

     use Simulation_data

     implicit none

     real, intent(in) :: rcc
     real, intent(out) :: dens, temp, ye, velr

     ! density
     dens = stepFunc( MaxD, MinD, H_D, R_D, rcc )

     ! temperature
     temp = stepFunc( MaxT, MinT, H_T, R_T, rcc )

     ! electron fraction
     ye   = stepFunc( MinY, MaxY, H_Y, R_Y, rcc )

     ! velocity
     velr = 0.0

     return
  end subroutine sim_analyticProfile

  real function stepFunc( y0, y1, h, r, rcc )

     implicit none

     real, intent(in) :: y0, y1, h, r, rcc
     real :: x

     x = (rcc - r) / h
     stepFunc = 0.5 * ( y0 * ( 1.0 - tanh(+x)) &
                      + y1 * ( 1.0 - tanh(-x)) )

     return
  end function stepFunc

  PURE REAL FUNCTION Interpolate1D_Linear( x, xL, xR, uL, uR )

    REAL, INTENT(in) :: x, xL, xR, uL, uR

    Interpolate1D_Linear &
      = ( ( xR - x ) * uL + ( x - xL ) * uR ) / ( xR - xL )

    RETURN
  END FUNCTION Interpolate1D_Linear

end subroutine Simulation_initBlock
