!!****if* source/IO/IOMain/IO_writeIntegralQuantities
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
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flashx.dat. Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData
!!REORDER(4):scratchData

#include "Simulation.h"

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalComm
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator, &
                             Grid_getCellVolumes

  use IO_data, ONLY : io_globalMe, io_writeMscalarIntegrals
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,     ONLY : Grid_tile_t
  use Simulation_data, ONLY : sim_maxDens

#if defined(THORNADO)
  use rt_data, ONLY : rt_offGridFluxR
  use RadiationFieldsModule, ONLY : iNuE, iNuE_Bar, &
    iGR_N, iGR_J, iGR_H1, iGR_H2, iGR_H3, nGR
  use PhysicalConstantsModule, ONLY : SpeedOfLightCGS 
#endif

  implicit none

#include "Flashx_mpi.h"
#include "constants.h"
  
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: funit = 99
  integer :: error
  integer :: nGlobalSumUsed, iSum
  
  character (len=MAX_STRING_LENGTH), save :: fname

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

#if defined(THORNADO)
  integer, parameter ::  nGlobalThornado = 2*THORNADO_NMOMENTS ! Number of globally-summed thornado quantities
#else
  integer, parameter ::  nGlobalThornado = 0
#endif
#ifdef MAGP_VAR
  integer, parameter ::  nGlobalSumProp = 11 + nGlobalThornado  ! Number of globally-summed regular quantities
#else
  integer, parameter ::  nGlobalSumProp = 10 + nGlobalThornado  ! Number of globally-summed regular quantities
#endif
  integer, parameter ::  nGlobalSum = nGlobalSumProp + NMASS_SCALARS ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer, parameter :: iGR_N_NuE      = SCRATCH_GRID_VARS_BEGIN + (iNuE    -1)*nGR + iGR_N - 1 
  integer, parameter :: iGR_N_NuE_Bar  = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar-1)*nGR + iGR_N - 1

  integer, parameter :: iGR_J_NuE      = SCRATCH_GRID_VARS_BEGIN + (iNuE    -1)*nGR + iGR_J - 1
  integer, parameter :: iGR_J_NuE_Bar  = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar-1)*nGR + iGR_J - 1

  integer, parameter :: iGR_H1_NuE     = SCRATCH_GRID_VARS_BEGIN + (iNuE    -1)*nGR + iGR_H1 - 1
  integer, parameter :: iGR_H1_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar-1)*nGR + iGR_H1 - 1

  integer, parameter :: iGR_H2_NuE     = SCRATCH_GRID_VARS_BEGIN + (iNuE    -1)*nGR + iGR_H2 - 1
  integer, parameter :: iGR_H2_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar-1)*nGR + iGR_H2 - 1

  integer, parameter :: iGR_H3_NuE     = SCRATCH_GRID_VARS_BEGIN + (iNuE    -1)*nGR + iGR_H3 - 1
  integer, parameter :: iGR_H3_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar-1)*nGR + iGR_H3 - 1

  integer :: ivar
  integer :: i, j, k
  integer :: lo(1:MDIM)
  integer :: hi(1:MDIM)
  real    :: dvol
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, DIMENSION(:,:,:,:), POINTER :: scratchData

  real :: maxDensLocal, maxDensGlobal

  integer :: ioStat
  
  real, allocatable :: cellVolumes(:,:,:)

  nullify(solnData)
  nullify(scratchData)

  if (io_writeMscalarIntegrals) then
     nGlobalSumUsed = nGlobalSum
  else
     nGlobalSumUsed = nGlobalSumProp
  end if

  ! Sum quantities over all locally held leaf-node blocks.
  gsum(1:nGlobalSumUsed) = 0.
  lsum(1:nGlobalSumUsed) = 0.

  maxDensLocal = 0.

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while (itor%isValid())
     call itor%currentTile(tileDesc)

     lo = tileDesc%limits(LOW,  :)
     hi = tileDesc%limits(HIGH, :)
     allocate(cellVolumes(lo(IAXIS):hi(IAXIS), &
                          lo(JAXIS):hi(JAXIS), &
                          lo(KAXIS):hi(KAXIS)))
     call Grid_getCellVolumes(tileDesc%level, &
                              lbound(cellVolumes), ubound(cellVolumes), &
                              cellVolumes)

     call tileDesc%getDataPtr(solnData, CENTER)
     call tileDesc%getDataPtr(scratchData, SCRATCH)

     ! Sum contributions from the indicated blkLimits of cells.
     do       k = lo(KAXIS), hi(KAXIS)
        do    j = lo(JAXIS), hi(JAXIS)
           do i = lo(IAXIS), hi(IAXIS)

              dvol = cellVolumes(i, j, k)

              ! mass
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol
#endif


#ifdef DENS_VAR
#ifdef VELX_VAR
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * &
                   &                                solnData(VELX_VAR,i,j,k)*dvol

#endif
#ifdef VELY_VAR

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * &
                   &                                solnData(VELY_VAR,i,j,k)*dvol

#endif
#ifdef VELZ_VAR
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * &
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * &
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#ifdef MAGP_VAR
              ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * &
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

              lsum(5) = lsum(5) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
#endif


#ifdef VELX_VAR
#ifdef VELY_VAR
#ifdef VELZ_VAR
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * &
                   &                             (solnData(VELX_VAR,i,j,k)**2+ &
                   &                              solnData(VELY_VAR,i,j,k)**2+ &
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol

#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * &
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif

              !neutrino lepton number
              lsum(8) = lsum(8) + (scratchData(iGR_N_NuE,i,j,k) &
                                -  scratchData(iGR_N_NuE_Bar,i,j,k) )*dvol

              !NuE energy = \int (J_e + 2 v^i H_{ei}) dV
              lsum(9) = lsum(9) + (scratchData(iGR_J_NuE,i,j,k) + 2.0d0 &
                                * (solnData(VELX_VAR,i,j,k)*scratchData(iGR_H1_NuE,i,j,k) &
                                 + solnData(VELY_VAR,i,j,k)*scratchData(iGR_H2_NuE,i,j,k) &
                                 + solnData(VELZ_VAR,i,j,k)*scratchData(iGR_H3_NuE,i,j,k)) &
                                /SpeedOfLightCGS**2)*dvol

              !NuE_Bar energy = \int (J_ebar + 2 v^i H_{ebar i}) dV
              lsum(10) = lsum(10) + (scratchData(iGR_J_NuE_Bar,i,j,k) + 2.0d0 &
                                  * (solnData(VELX_VAR,i,j,k)*scratchData(iGR_H1_NuE_Bar,i,j,k) &
                                   + solnData(VELY_VAR,i,j,k)*scratchData(iGR_H2_NuE_Bar,i,j,k) &
                                   + solnData(VELZ_VAR,i,j,k)*scratchData(iGR_H3_NuE_Bar,i,j,k)) &
                                  /SpeedOfLightCGS**2)*dvol
#endif ! ifdef DENS_VAR

#ifdef MAGP_VAR
              ! magnetic energy
              lsum(11) = lsum(11) + solnData(MAGP_VAR,i,j,k)*dvol
#endif

#ifdef DENS_VAR
              if (io_writeMscalarIntegrals) then
                 iSum = nGlobalSumProp
!!$                 do ivar=MASS_SCALARS_BEGIN,MASS_SCALARS_END
                    lsum(iSum+1:iSum+NMASS_SCALARS) = &
                         lsum(iSum+1:iSum+NMASS_SCALARS) + &
                           solnData(DENS_VAR,i,j,k) * &
                           solnData(MASS_SCALARS_BEGIN: &
                                    MASS_SCALARS_END,i,j,k)*dvol
!!$                 end do
              end if

              ! Peak density
              maxDensLocal = max(maxDensLocal,solnData(DENS_VAR,i,j,k))
#endif
           enddo
        enddo
     enddo
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call tileDesc%releaseDataPtr(scratchData, SCRATCH)

     deallocate(cellVolumes)

     call itor%next()
  enddo
  call Grid_releaseTileIterator(itor)

#if defined(THORNADO)
  iSum = nGlobalSumProp - nGlobalThornado
  lsum(iSum+1:iSum+2*THORNADO_NMOMENTS) = rt_offGridFluxR
#endif

  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.

  call MPI_Reduce (lsum, gsum, nGlobalSumUsed, FLASH_REAL, MPI_SUM, &
       &                MASTER_PE, io_globalComm, error)

  lsum(1) = maxDensLocal
  call MPI_AllReduce ( lsum(1), lsum(2), 1, FLASH_REAL, MPI_MAX, &
     &  io_globalComm, error)
  maxDensGlobal = lsum(2)
  sim_maxDens = maxDensGlobal


  if (io_globalMe  == MASTER_PE) then

     ! create the file from scratch if it is a not a restart simulation,
     ! otherwise append to the end of the file

     !No matter what, we are opening the file. Check to see if already there
     ioStat = 0
     open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
     if(ioStat .NE. 0) then
        !print *, 'FILE FOUND'
        open(funit, file=trim(io_statsFileName), position='APPEND')
     endif

     if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then

#ifndef MAGP_VAR
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'x-momentum                ', &
             'y-momentum                ', &
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'e-type nu lepton number   ', &
             'NuE energy                ', &
             'NuE_Bar energy            ', &
             (msName(ivar),ivar=MASS_SCALARS_BEGIN,&
              min(MASS_SCALARS_END,&
                  MASS_SCALARS_BEGIN+nGlobalSumUsed-nGlobalSumProp-1))

#else

        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'x-momentum                ', &
             'y-momentum                ', &
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'e-type nu lepton number   ', &
             'NuE energy                ', &
             'NuE_Bar energy            ', &
             'MagEnergy                 ', &
             (msName(ivar),ivar=MASS_SCALARS_BEGIN,&
              min(MASS_SCALARS_END,&
                  MASS_SCALARS_BEGIN+nGlobalSumUsed-nGlobalSumProp-1))
#endif

10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11)
11      format('# simulation restarted')
     endif

     ! Write the global sums to the file.
     write (funit, 12) simtime, gsum(1:nGlobalSumUsed)

12   format (1x, 50(es25.18, :, 1x))

     close (funit)          ! Close the file.

  endif

#ifdef USEBARS
  call MPI_Barrier (io_globalComm, error)
#endif

  !=============================================================================

  return

  contains
    character(len=25) function msName(ivar)
      integer,intent(in) :: ivar
      character(len=25) :: str
      call Simulation_mapIntToStr(ivar,str,MAPBLOCK_UNK)
      msName = str
    end function msName
end subroutine IO_writeIntegralQuantities



