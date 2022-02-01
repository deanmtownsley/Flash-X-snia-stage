!!****if* source/Simulation/SimulationMain/DeleptonizationWave/Simulation_init
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
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!   Initialize general solution data for deleptonization wave radiation
!!   problem (steady-state)
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use ProgramHeaderModule, ONLY : nE, nDOF
  use RadiationFieldsModule, ONLY : nCR, nSpecies
  use Driver_interface, ONLY : Driver_getMype
  use Logfile_interface, ONLY : Logfile_stamp
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
     RuntimeParameters_mapStrToInt

  implicit none

#include "constants.h"
#include "Simulation.h"

  call Logfile_stamp( 'Entering simulation initialization' , '[Simulation_init]')
  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)

  call RuntimeParameters_get( 'restart', sim_restart)
  call RuntimeParameters_get("geometry",sim_str_geometry)
  call RuntimeParameters_mapStrToInt(sim_str_geometry, sim_geometry)

  call RuntimeParameters_get("sim_use_model", sim_use_model)
  call RuntimeParameters_get("sim_rad_option", sim_rad_option)
  call RuntimeParameters_get("sim_model_file", sim_model_file)
  call RuntimeParameters_get("sim_rintSwitch", sim_rintSwitch)

  sim_xn_i = 0.0e0
  if( NSPECIES > 0 ) sim_xn_i(SPECIES_BEGIN) = 1.0e0

  sim_velx_i = 0.0e0
  sim_vely_i = 0.0e0
  sim_velz_i = 0.0e0

  sim_nComp = nSpecies * nCR * nE * nDOF

  ! read in profile if not restart
  if ( .not. sim_restart ) then

    ! fluid field
    if ( sim_use_model ) call sim_readProfile
    ! radiation field
    select case ( sim_rad_option )
      case ( -1 )
        if (sim_meshMe == MASTER_PE) &
        print*, 'Using Zero distribution function to initial radiation field ...'
      case ( 0 )
        if (sim_meshMe == MASTER_PE) &
        print*, 'Using FD distribution to initial neutrino number density and zero out flux ...'
      case ( 1 )
        if (sim_meshMe == MASTER_PE) &
        print*, 'Computing approximate neutrino sphere radii to initial radiation field ...'
        call sim_computeSphereSolution
      case ( 2 )
        if (sim_meshMe == MASTER_PE) &
        print*, 'Looking into profile for initializing radiation field ...'
        call sim_readChimeraProfile_rad
      case ( 3 )
        if (sim_meshMe == MASTER_PE) &
        print*, 'Looking into profile for initializing radiation field ...'
        call sim_readBoltzTranProfile_rad
      case default
        if (sim_meshMe == MASTER_PE) &
        print*, 'Using default radtion field initialization: no neutrino ...'
        sim_rad_option = 0
    end select
  end if

  return

contains

  subroutine sim_readProfile

     use, intrinsic :: iso_fortran_env, only: iostat_end
     use Driver_interface, ONLY : Driver_abort
     use Simulation_interface, ONLY : Simulation_mapStrToInt
     use Simulation_data, ONLY: nvar_stored, sim_model_file

     implicit none

#include "constants.h"
#include "Simulation.h"

     character(128) :: current_line
     character(4) :: var_label
     integer :: fileUnit
     integer :: var_key(NUNK_VARS)
     integer :: ipos, i, j, k, ierr
     real, allocatable :: var_temp(:)

     ! open the file and read in the header
     open(newunit=fileUnit,file=trim(sim_model_file),status='old',iostat=ierr)
     if (ierr /= 0) call Driver_abort('Unable to open initial model file')
     if (sim_meshMe == MASTER_PE) &
     print *, 'Profile ', trim(sim_model_file), ' opened'
     read(fileUnit,'(A80)') current_line

     ! read in the number of variables line
     read(fileUnit,'(A80)') current_line
     ipos = index(current_line,'=') + 1
     read (current_line(ipos:),*) nvar_stored
     if (sim_meshMe == MASTER_PE) print *,"read nvar_stored", nvar_stored


     ! read in the number of points in profile
     read(fileUnit,'(A80)') current_line
     ipos = index(current_line,'=') + 1
     read (current_line(ipos:),*) n1d_max
     if (sim_meshMe == MASTER_PE) print *,"with ",n1d_max, " points"

     if (NUNK_VARS .NE. nvar_stored .AND. sim_meshMe == MASTER_PE) then
        print *, ' '
        print *, 'Warning: the number of variables stored in the'
        print *, 'input file is different than the number of'
        print *, 'variables in the current version of FLASH.'
        print *, ' '
        print *, 'The variables in the file that are also defined'
        print *, 'in FLASH will be read in.  Any missing variables'
        print *, 'will be initialized to zero'
        print *, ' '
     endif

     if (sim_meshMe == MASTER_PE) print *, "Vaiables in file:"
     var_key = NONEXISTENT
     do i = 1, nvar_stored
        read (fileUnit,*) var_label
        if (sim_meshMe == MASTER_PE) print *, var_label
        call Simulation_mapStrToInt(trim(adjustl(var_label)), j, MAPBLOCK_UNK)
        if ( j /= NONEXISTENT ) var_key(j) = i
     enddo

     allocate( xzn(n1d_max) )
     allocate( model_1d(n1d_max,NUNK_VARS) )
     xzn = 0.0
     model_1d = 0.0

     allocate( var_temp(nvar_stored) )
     do k = 1, n1d_max
        read(fileUnit,*,iostat=ierr) xzn(k), (var_temp(i),i=1,nvar_stored)
        If ( ierr == iostat_end ) exit

        ! put these in order, so model1d_var always contains the same variables
        ! in the same spots
        do j = 1, NUNK_VARS
           if (var_key(j) /= NONEXISTENT) model_1d(k,j) = var_temp(var_key(j))
        enddo
     enddo
     deallocate( var_temp )

     close(fileUnit, status = 'keep')

     if (sim_meshMe == MASTER_PE) &
     print *, 'Profile (', trim(sim_model_file),') read completed.'

  end subroutine sim_readProfile

  subroutine sim_readBoltzTranProfile_rad

    use Simulation_data, only : &
        sim_model_file, n1d_nE, ezn, xznrad, model_1d_rad

    character(40) :: radfilename
    integer :: fileUnit, ierr, nlength

    integer, parameter :: nE = 20
    integer, parameter :: nX = 306
    integer, parameter :: nM = 2
    integer, parameter :: nS = 2

    character(100) :: current_line
    integer :: iline, ii, iE, ipos, iS, iR, iX, ia
    integer :: totalline = 12529

    integer :: breakline
    real    :: buffer1, buffer2, buffer3

    nlength = len(trim(sim_model_file))
    radfilename = sim_model_file(1:nlength-1)//'r_Boltz'

    ! allocate and pass parameter
    n1d_nE = nE
    allocate( ezn(nE) )
    allocate( xznrad(nX) )
    allocate( model_1d_rad(nE,nX,nM,nS) )

    open(newunit=fileUnit,file=trim(radfilename),status='old',iostat=ierr)
    if (ierr /= 0) call Driver_abort('Unable to open initial model file')
    if (sim_meshMe == MASTER_PE) &
    print*, 'Profile ', trim(radfilename), ' opened'

    iline = 0
    iS = 0
    breakline = 0

    do ii = 1, totalline
      read(fileUnit,'(A80)') current_line
      ! -- seperate species
      if( current_line(1:13) == '   MOMENTS of' )then
         iS = iS + 1
         iE = 0
      end if
      ! -- read ezn values
      if( (current_line(1:4) == '  E=') .and. (iE < nE) ) then
        iE = iE + 1
        ! -- only load ezn once
        if ( iS == 1 ) then
          ipos = index(current_line,'=') + 1
          read (current_line(ipos:),*) ezn(iE)
        end if
      end if
      ! -- find break lines
      if( current_line(1:10) == '----------' ) then
        iline = iline + 1
        breakline = ii
      end if
      ! -- load data to model_1d_rad
      if( (ii >= breakline+2) .and. (ii <= breakline+1+nX) &
          .and. (breakline .ne. 0) ) then
        iX = ii-breakline-1
        read(current_line,*) &
        iR, xznrad(iX), buffer1, buffer2, &
        model_1d_rad(iE,iX,1,iS), model_1d_rad(iE,iX,2,iS), buffer3
        if( iR .ne. iX ) stop 'wrong profile reader'
        if( model_1d_rad(iE,iX,1,iS) < 0.0e0)then
           !!$write(*,'(A,4I4,A)') &
           !!$'  => triggeredd max(J,1.0e-40) and max(H1,0.0e0) @ (iE,iX,1,iS)=(',iE,iX,1,iS,')'
           model_1d_rad(iE,iX,1,iS) = max( model_1d_rad(iE,iX,1,iS), 1.0e-40 )
           model_1d_rad(iE,iX,2,iS) = max( model_1d_rad(iE,iX,2,iS), 0.0e0 )
           if( model_1d_rad(iE,iX,1,iS) < 0.0e0) stop 'failed enforce J > 1.0e-40'
        end if
      end if
      ! -- IMPORTANT --
      ! FRIST and LAST zone (307) data is not usable. DON'T USE.
      ! --------------
    end do

    close(fileUnit, status = 'keep')

    if (sim_meshMe == MASTER_PE) then
       print *, 'Profile ', trim(radfilename),' read completed :'
       write(*,'(A,2ES15.6)') ' ezn range [MeV]:', ezn(1), ezn(nE)
       write(*,'(A,2ES15.6)') ' xznrad range [CM] :', xznrad(1), xznrad(nX)
       write(*,'(A,2ES15.6)') ' model_1d_rad range:', minval(model_1d_rad), maxval(model_1d_rad)
       write(*,'(A,2ES15.6)') '   in which J range:', &
         minval(model_1d_rad(:,:,1,:)), maxval(model_1d_rad(:,:,1,:))
       write(*,'(A,2ES15.6)') '           H1 range:', &
         minval(model_1d_rad(:,:,2,:)), maxval(model_1d_rad(:,:,2,:))
       print*, ''
    end if

  end subroutine sim_readBoltzTranProfile_rad

  subroutine sim_readChimeraProfile_rad

     use Driver_interface, ONLY : Driver_abort
     use Simulation_data, ONLY : n1d_max, n1d_nE, ezn, model_1d_rad, &
                                 sim_model_file, xznrad, xzn

#include "Simulation.h"

     implicit none

     character(40) :: filename0, filename1
     character(128) :: current_line
     integer :: iS, pro_n, pro_i, nlength
     integer :: fileUnit, ierr, nline, ipos, iE, iX, il
     real :: x_c

     ! Chimera radiation field has the same radius as as fluid field
     allocate( xznrad(n1d_max) )
     xznrad = xzn

     ! open Chimera_G25+100ms_rPsi0_1 for nE and nline
     nlength = len(trim(sim_model_file))
     filename0 = sim_model_file(1:nlength-1)//'rPsi0_1'
     open(newunit=fileUnit,file=trim(filename0),status='old',iostat=ierr)
     if (ierr /= 0) &
       call Driver_abort('Unable to open initial model file: '//trim(filename0))
     read (fileUnit,'(A)') current_line
     read(fileUnit,'(A)') current_line
     ipos = index(current_line,'=') + 1
     read (current_line(ipos:),*) nline
     read(fileUnit,'(A80)') current_line
     ipos = index(current_line,'=') + 1
     read (current_line(ipos:),*) n1d_nE
     if (sim_meshMe == MASTER_PE) &
       print *, trim(filename0),' containts points = ', nline, 'nE = ',n1d_nE
     allocate( ezn(n1d_nE) )
     read(fileUnit,'(A80)') current_line
     read(fileUnit,*,iostat=ierr) (ezn(iE),iE=1,n1d_nE)
     close(fileUnit, status = 'keep')

     if( nline /= n1d_max ) call Driver_abort('Profiles mismatching.')
     allocate( model_1d_rad(n1d_nE,n1d_max,2,THORNADO_NSPECIES) )

     do iS = 1, THORNADO_NSPECIES
       write(filename0,'(A,I1)') sim_model_file(1:nlength-1)//'rPsi0_',iS
       write(filename1,'(A,I1)') sim_model_file(1:nlength-1)//'rPsi1_',iS

       open(newunit=fileUnit,file=trim(filename0),status='old',iostat=ierr)
       if (ierr /= 0) call Driver_abort('Unable to open initial model file:'//trim(filename0))
       if (sim_meshMe == MASTER_PE) print *, 'Profile ', trim(fileName0), ' opened'
       do il = 1, 7
         read(fileUnit,'(A)') current_line
       end do
       do iX = 1, nline
         read(fileUnit,*,iostat=ierr) pro_n, pro_i, x_c, (model_1d_rad(iE,iX,1,iS),iE=1,n1d_nE)
       end do
       close(fileUnit, status = 'keep')

     ! open the file and read in the header
       open(newunit=fileUnit,file=trim(filename1),status='old',iostat=ierr)
       if (ierr /= 0) call Driver_abort('Unable to open initial model file:'//trim(filename1))
       if (sim_meshMe == MASTER_PE) print *, 'Profile ', trim(fileName1), ' opened'
       do il = 1, 7
         read(fileUnit,'(A)') current_line
       end do
       do iX = 1, nline
         read(fileUnit,*,iostat=ierr) pro_n, pro_i, x_c, (model_1d_rad(iE,iX,2,iS),iE=1,n1d_nE)
       end do
       close(fileUnit, status = 'keep')

     end do

  end subroutine sim_readChimeraProfile_rad

  subroutine sim_computeSphereSolution

    use MeshModule, ONLY : NodeCoordinate, MeshE
    use UnitsModule, ONLY : MeV, Gram, Centimeter, Kelvin, Kilometer
    use Simulation_data, ONLY : n1d_max, xzn, D_Nu_P, I1_Nu_P
    use NeutrinoOpacitiesComputationModule, ONLY: &
        ComputeEquilibriumDistributions_Points, &
        ComputeNeutrinoOpacities_EC_Points
    use ut_interpolationInterface, ONLY : ut_hunt

#include "Simulation.h"

    integer :: iE, iNodeE, nE, nR, i, iS, iR, buff_int
    real :: enode, Tau
    real, allocatable :: R_P(:), D_P(:), T_P(:), Y_P(:)
    real, allocatable :: Chi(:,:,:), fEQ(:,:,:), R_Nu(:,:), E_Nu(:)

    real, parameter :: UnitD  = Gram / Centimeter**3
    real, parameter :: UnitT  = Kelvin
    real, parameter :: UnitY  = 1.0
    
    nR = n1d_max 
    nE = THORNADO_NE * THORNADO_NNODESE

    allocate( E_Nu(nE) )
    allocate( R_P(nR) )
    allocate( D_P(nR) )
    allocate( T_P(nR) )
    allocate( Y_P(nR) )
    allocate( R_Nu(nE,THORNADO_NSPECIES) )
    allocate( Chi(nE,nR,THORNADO_NSPECIES) ) 
    allocate( fEQ(nE,nR,THORNADO_NSPECIES) ) 
    allocate( D_Nu_P(nE,nR,THORNADO_NSPECIES) ) 
    allocate( I1_Nu_P(nE,nR,THORNADO_NSPECIES) ) 

    ! convert to thornado's unit
    R_P = xzn * Centimeter
    D_P = model_1d(:,DENS_VAR) * UnitD
    T_P = model_1d(:,TEMP_VAR) * UnitT
    Y_P = model_1d(:,YE_MSCALAR) * UnitY
    ! prevent too low temperature
    call ut_hunt(R_P,nR,5.0e3*Kilometer,buff_int)
    D_P = max( D_P, D_P(buff_int) )
    T_P = max( T_P, T_P(buff_int) )
    ! get thornado neutrino energies node value
    i = 1
    do iE = 1, THORNADO_NE
      do iNodeE = 1, THORNADO_NNODESE
        enode = NodeCoordinate( MeshE,iE,iNodeE)
        E_Nu(i) = enode
        i = i+1
      end do
    end do
    ! compute neutrino absorption opacities
    do iS = 1, THORNADO_NSPECIES
      call ComputeNeutrinoOpacities_EC_Points &
             ( 1, nE, 1, nR, E_Nu, D_P, T_P, Y_P, iS, Chi(:,:,iS) )
      call ComputeEquilibriumDistributions_Points &
             ( 1, nE, 1, nR, E_Nu, D_P, T_P, Y_P, fEQ(:,:,iS), iS )
    end do
    ! compute approximate neutrino sphere radii
    do iS = 1, THORNADO_NSPECIES; do iE = 1, nE
      Tau = 0.0e0
      do iR = nR-1, 1, -1
        if( Tau > 2.0e0 / 3.0e0 ) cycle
        Tau = Tau + 0.5e0 * ( R_P(iR+1) - R_P(iR) ) &
                          * ( Chi(iE,iR+1,iS) + Chi(iE,iR,iS) )
        R_Nu(iE,iS) = MAX( R_P(iR), 1.0e1 * Kilometer )
      end do 
    end do; end do
    ! use Chi and fEQ at local radius or neutrino sphere
    do iS = 1, THORNADO_NSPECIES; do iR = 1, nR; do iE = 1, nE
      call ut_hunt(R_P,nR,R_Nu(iE,iS),buff_int)    
      i = MIN( iR, buff_int )
      call ComputeSphereSolution &
               ( R_Nu(iE,iS), Chi(iE,i,iS), fEQ(iE,i,iS), &
                 R_P(iR), D_Nu_P(iE,iR,iS), I1_Nu_P(iE,iR,iS) )
    end do; end do; end do

    deallocate( R_P, D_P, T_P, Y_P )
    deallocate( E_Nu, R_Nu, Chi, fEQ )

    if (sim_meshMe == MASTER_PE) &
      print*, 'Compute Sphere Solution completed'
    return

  end subroutine sim_computeSphereSolution

  subroutine ComputeSphereSolution( R0, Chi, f0, R, D, I )

    REAL, INTENT(in)  :: R0, Chi, f0, R
    REAL, INTENT(out) :: D, I

    INTEGER, PARAMETER :: nMu = 2048
    INTEGER            :: iMu
    REAL               :: Mu(nMu), Distribution(nMu)

    DO iMu = 1, nMu

      Mu(iMu) = - 1.0e0 + 2.0e0 * DBLE(iMu-1)/DBLE(nMu-1)

      Distribution(iMu) = f_A( R, Mu(iMu), R0, f0, Chi )

    END DO

    D = 0.5e0 * TRAPEZ( nMu, Mu, Distribution )
    I = 0.5e0 * TRAPEZ( nMu, Mu, Distribution * Mu )

    D = max( D, 1.0e-100 )

  end subroutine ComputeSphereSolution

  real function TRAPEZ( n, x, y )

    INTEGER,  INTENT(in) :: n
    REAL    , INTENT(in) :: x(n), y(n)

    INTEGER :: i

    TRAPEZ = 0.0e0
    DO i = 1, n - 1
      TRAPEZ = TRAPEZ + 0.5e0 * ( x(i+1) - x(i) ) * ( y(i) + y(i+1) )
    END DO

    RETURN
  end function TRAPEZ

  REAL FUNCTION f_A( R, Mu, R0, f0, Chi )

    REAL, INTENT(in) :: R, Mu, R0, f0, Chi

    REAL :: s

    IF( R < R0 )THEN
      s = ( R * Mu + R0 * SQRT( 1.0e0 - ( R / R0 )**2 * ( 1.0e0 - Mu**2 ) ) )
    ELSE
      IF( Mu >= SQRT( 1.0e0 - ( R0 / R )**2 ) )THEN
        s = ( 2.0e0 * R0 * SQRT( 1.0e0 - ( R / R0 )**2 * ( 1.0e0 - Mu**2 ) ) )
      ELSE
        s = 0.0e0
      END IF
    END IF

    f_A = f0 * ( 1.0e0 - EXP( - Chi * s ) )

    RETURN
  END FUNCTION f_A

end subroutine Simulation_init
