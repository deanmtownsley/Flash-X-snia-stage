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

subroutine IO_writeIntegralQuantities(isFirst, simTime)

   use IO_data, ONLY: io_restart, io_statsFileName, io_globalComm
   use Grid_interface, ONLY: Grid_getTileIterator, &
                             Grid_releaseTileIterator, &
                             Grid_getCellVolumes

   use IO_data, ONLY: io_globalMe, io_writeMscalarIntegrals
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t
   use Simulation_data, ONLY: sim_maxDens, sim_postBounce, sim_bounceTime

   ! use Eos_wlInterface, only: Eos_wlDetectBounce

#if defined(THORNADO)
   use rt_data, ONLY: rt_offGridFluxR
   use RadiationFieldsModule, ONLY: iNuE, iNuE_Bar, &
                                    iNuM, iNuM_Bar, iNuT, iNuT_Bar, &
                                    iGR_N, iGR_J, iGR_H1, iGR_H2, iGR_H3, nGR
   use PhysicalConstantsModule, ONLY: SpeedOfLightCGS
   use UnitsModule, ONLY: PlanckConstant, SpeedOfLight, &
                          Second, Erg
#endif

   implicit none

#include "Flashx_mpi.h"
#include "constants.h"

   real, intent(in) :: simTime

   integer, intent(in) :: isFirst

   integer :: funit = 99
   integer :: error
   integer :: nGlobalSumUsed, iSum

   character(len=MAX_STRING_LENGTH), save :: fname

   type(Grid_iterator_t) :: itor
   type(Grid_tile_t)     :: tileDesc

#if defined(THORNADO)
! Number of globally-summed thornado quantities:
! total electron type neutrino lepton number in domain,
! total NuE energy in domain,
! total NuE_Bar energy in domain,
! (total NuM energy in domain),
! (total NuM_Bar energy in domain),
! (total NuT energy in domain),
! (total NuT_Bar energy in domain),
! total neutrino lepton number that left grid,
! total neutrino energy that left grid,
! total neutrino (x,y,z)-momentum that left grid

   real, parameter :: hc3 = (PlanckConstant*SpeedOfLight)**3

   integer, parameter :: nOffGridThornado = 5

   integer, parameter :: iGR_N_NuE = SCRATCH_GRID_VARS_BEGIN + (iNuE - 1)*nGR + iGR_N - 1
   integer, parameter :: iGR_N_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar - 1)*nGR + iGR_N - 1

   integer, parameter :: iGR_J_NuE = SCRATCH_GRID_VARS_BEGIN + (iNuE - 1)*nGR + iGR_J - 1
   integer, parameter :: iGR_J_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar - 1)*nGR + iGR_J - 1

   integer, parameter :: iGR_H1_NuE = SCRATCH_GRID_VARS_BEGIN + (iNuE - 1)*nGR + iGR_H1 - 1
   integer, parameter :: iGR_H1_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar - 1)*nGR + iGR_H1 - 1

   integer, parameter :: iGR_H2_NuE = SCRATCH_GRID_VARS_BEGIN + (iNuE - 1)*nGR + iGR_H2 - 1
   integer, parameter :: iGR_H2_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar - 1)*nGR + iGR_H2 - 1

   integer, parameter :: iGR_H3_NuE = SCRATCH_GRID_VARS_BEGIN + (iNuE - 1)*nGR + iGR_H3 - 1
   integer, parameter :: iGR_H3_NuE_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuE_Bar - 1)*nGR + iGR_H3 - 1

#if THORNADO_NSPECIES == 6

   integer, parameter :: nVolumeThornado = 7
   integer, parameter :: nGlobalThornado = nVolumeThornado + nOffGridThornado

   integer, parameter :: iGR_J_NuM = SCRATCH_GRID_VARS_BEGIN + (iNuM - 1)*nGR + iGR_J - 1
   integer, parameter :: iGR_J_NuM_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuM_Bar - 1)*nGR + iGR_J - 1

   integer, parameter :: iGR_H1_NuM = SCRATCH_GRID_VARS_BEGIN + (iNuM - 1)*nGR + iGR_H1 - 1
   integer, parameter :: iGR_H1_NuM_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuM_Bar - 1)*nGR + iGR_H1 - 1

   integer, parameter :: iGR_H2_NuM = SCRATCH_GRID_VARS_BEGIN + (iNuM - 1)*nGR + iGR_H2 - 1
   integer, parameter :: iGR_H2_NuM_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuM_Bar - 1)*nGR + iGR_H2 - 1

   integer, parameter :: iGR_H3_NuM = SCRATCH_GRID_VARS_BEGIN + (iNuM - 1)*nGR + iGR_H3 - 1
   integer, parameter :: iGR_H3_NuM_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuM_Bar - 1)*nGR + iGR_H3 - 1

   integer, parameter :: iGR_J_NuT = SCRATCH_GRID_VARS_BEGIN + (iNuT - 1)*nGR + iGR_J - 1
   integer, parameter :: iGR_J_NuT_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuT_Bar - 1)*nGR + iGR_J - 1

   integer, parameter :: iGR_H1_NuT = SCRATCH_GRID_VARS_BEGIN + (iNuT - 1)*nGR + iGR_H1 - 1
   integer, parameter :: iGR_H1_NuT_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuT_Bar - 1)*nGR + iGR_H1 - 1

   integer, parameter :: iGR_H2_NuT = SCRATCH_GRID_VARS_BEGIN + (iNuT - 1)*nGR + iGR_H2 - 1
   integer, parameter :: iGR_H2_NuT_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuT_Bar - 1)*nGR + iGR_H2 - 1

   integer, parameter :: iGR_H3_NuT = SCRATCH_GRID_VARS_BEGIN + (iNuT - 1)*nGR + iGR_H3 - 1
   integer, parameter :: iGR_H3_NuT_Bar = SCRATCH_GRID_VARS_BEGIN + (iNuT_Bar - 1)*nGR + iGR_H3 - 1
#else
   integer, parameter :: nVolumeThornado = 3
   integer, parameter :: nGlobalThornado = nVolumeThornado + nOffGridThornado

#endif
#else
   integer, parameter ::  nGlobalThornado = 0
#endif
#ifdef MAGP_VAR
   integer, parameter ::  nGlobalSumProp = 8 + nGlobalThornado  ! Number of globally-summed regular quantities
#else
   integer, parameter ::  nGlobalSumProp = 7 + nGlobalThornado  ! Number of globally-summed regular quantities
#endif
   integer, parameter ::  nGlobalSum = nGlobalSumProp + NMASS_SCALARS ! Number of globally-summed quantities
   real :: gsum(nGlobalSum) !Global summed quantities
   real :: lsum(nGlobalSum) !Global summed quantities

   integer :: ivar
   integer :: i, j, k
   integer :: lo(1:MDIM)
   integer :: hi(1:MDIM)
   real    :: dvol
   real, DIMENSION(:, :, :, :), POINTER :: solnData
   real, DIMENSION(:, :, :, :), POINTER :: scratchData

   real :: maxDensLocal, maxDensGlobal

   integer :: ioStat

   real, allocatable :: cellVolumes(:, :, :)

   real :: postBounceTime

   nullify (solnData)
   nullify (scratchData)

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

      lo = tileDesc%limits(LOW, :)
      hi = tileDesc%limits(HIGH, :)
      allocate (cellVolumes(lo(IAXIS):hi(IAXIS), &
                            lo(JAXIS):hi(JAXIS), &
                            lo(KAXIS):hi(KAXIS)))
      call Grid_getCellVolumes(tileDesc%level, &
                               lbound(cellVolumes), ubound(cellVolumes), &
                               cellVolumes)

      call tileDesc%getDataPtr(solnData, CENTER)
      call tileDesc%getDataPtr(scratchData, SCRATCH)

      ! Sum contributions from the indicated blkLimits of cells.
      do k = lo(KAXIS), hi(KAXIS)
         do j = lo(JAXIS), hi(JAXIS)
            do i = lo(IAXIS), hi(IAXIS)

               dvol = cellVolumes(i, j, k)

               ! mass
#ifdef DENS_VAR
               lsum(1) = lsum(1) + solnData(DENS_VAR, i, j, k)*dvol
#endif

#ifdef DENS_VAR
#ifdef VELX_VAR
               ! momentum
               lsum(2) = lsum(2) + solnData(DENS_VAR, i, j, k)* &
                    &                                solnData(VELX_VAR, i, j, k)*dvol

#endif
#ifdef VELY_VAR

               lsum(3) = lsum(3) + solnData(DENS_VAR, i, j, k)* &
                    &                                solnData(VELY_VAR, i, j, k)*dvol

#endif
#ifdef VELZ_VAR
               lsum(4) = lsum(4) + solnData(DENS_VAR, i, j, k)* &
                    &                                solnData(VELZ_VAR, i, j, k)*dvol
#endif

               ! total energy
#ifdef ENER_VAR
               lsum(5) = lsum(5) + solnData(ENER_VAR, i, j, k)* &
                    &                                solnData(DENS_VAR, i, j, k)*dvol
#ifdef MAGP_VAR
               ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * &
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

               lsum(5) = lsum(5) + solnData(MAGP_VAR, i, j, k)*dvol
#endif
#endif

#ifdef VELX_VAR
#ifdef VELY_VAR
#ifdef VELZ_VAR
               ! kinetic energy
               lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR, i, j, k)* &
                    &                             (solnData(VELX_VAR, i, j, k)**2 + &
                    &                              solnData(VELY_VAR, i, j, k)**2 + &
                    &                              solnData(VELZ_VAR, i, j, k)**2)*dvol

#endif
#endif
#endif

#ifdef EINT_VAR
               ! internal energy
               lsum(7) = lsum(7) + solnData(DENS_VAR, i, j, k)* &
                    &                                solnData(EINT_VAR, i, j, k)*dvol
#endif

#if defined(THORNADO)
               !electron type neutrino lepton number
               lsum(8) = lsum(8) + (scratchData(iGR_N_NuE, i, j, k) &
                                    - scratchData(iGR_N_NuE_Bar, i, j, k))*dvol

               !NuE energy = \int (J_e + 2 v^i H_{ei}) dV
               lsum(9) = lsum(9) + (scratchData(iGR_J_NuE, i, j, k) + 2.0d0 &
                                    *(solnData(VELX_VAR, i, j, k)*scratchData(iGR_H1_NuE, i, j, k) &
                                      + solnData(VELY_VAR, i, j, k)*scratchData(iGR_H2_NuE, i, j, k) &
                                      + solnData(VELZ_VAR, i, j, k)*scratchData(iGR_H3_NuE, i, j, k)) &
                                    /SpeedOfLightCGS**2)*dvol

               !NuE_Bar energy = \int (J_ebar + 2 v^i H_{ebar i}) dV
               lsum(10) = lsum(10) + (scratchData(iGR_J_NuE_Bar, i, j, k) + 2.0d0 &
                                      *(solnData(VELX_VAR, i, j, k)*scratchData(iGR_H1_NuE_Bar, i, j, k) &
                                        + solnData(VELY_VAR, i, j, k)*scratchData(iGR_H2_NuE_Bar, i, j, k) &
                                        + solnData(VELZ_VAR, i, j, k)*scratchData(iGR_H3_NuE_Bar, i, j, k)) &
                                      /SpeedOfLightCGS**2)*dvol
#if THORNADO_NSPECIES == 6
               !NuM energy = \int (J_{\mu} + 2 v^i H_{\mu i}) dV
               lsum(11) = lsum(11) + (scratchData(iGR_J_NuM, i, j, k) + 2.0d0 &
                                      *(solnData(VELX_VAR, i, j, k)*scratchData(iGR_H1_NuM, i, j, k) &
                                        + solnData(VELY_VAR, i, j, k)*scratchData(iGR_H2_NuM, i, j, k) &
                                        + solnData(VELZ_VAR, i, j, k)*scratchData(iGR_H3_NuM, i, j, k)) &
                                      /SpeedOfLightCGS**2)*dvol

               !NuM_Bar energy = \int (J_{\mu bar} + 2 v^i H_{\mu bar i}) dV
               lsum(12) = lsum(12) + (scratchData(iGR_J_NuM_Bar, i, j, k) + 2.0d0 &
                                      *(solnData(VELX_VAR, i, j, k)*scratchData(iGR_H1_NuM_Bar, i, j, k) &
                                        + solnData(VELY_VAR, i, j, k)*scratchData(iGR_H2_NuM_Bar, i, j, k) &
                                        + solnData(VELZ_VAR, i, j, k)*scratchData(iGR_H3_NuM_Bar, i, j, k)) &
                                      /SpeedOfLightCGS**2)*dvol

               !NuT energy = \int (J_{\tau} + 2 v^i H_{\tau i}) dV
               lsum(13) = lsum(13) + (scratchData(iGR_J_NuT, i, j, k) + 2.0d0 &
                                      *(solnData(VELX_VAR, i, j, k)*scratchData(iGR_H1_NuT, i, j, k) &
                                        + solnData(VELY_VAR, i, j, k)*scratchData(iGR_H2_NuT, i, j, k) &
                                        + solnData(VELZ_VAR, i, j, k)*scratchData(iGR_H3_NuT, i, j, k)) &
                                      /SpeedOfLightCGS**2)*dvol

               !NuT_Bar energy = \int (J_{\tau bar} + 2 v^i H_{\tau bar i}) dV
               lsum(14) = lsum(14) + (scratchData(iGR_J_NuT_Bar, i, j, k) + 2.0d0 &
                                      *(solnData(VELX_VAR, i, j, k)*scratchData(iGR_H1_NuT_Bar, i, j, k) &
                                        + solnData(VELY_VAR, i, j, k)*scratchData(iGR_H2_NuT_Bar, i, j, k) &
                                        + solnData(VELZ_VAR, i, j, k)*scratchData(iGR_H3_NuT_Bar, i, j, k)) &
                                      /SpeedOfLightCGS**2)*dvol
#endif !
#endif
#endif ! ifdef DENS_VAR

#ifdef MAGP_VAR
               ! magnetic energy
               lsum(nGlobalSumProp - nGlobalThornado) = &
                  lsum(nGlobalSumProp - nGlobalThornado) + solnData(MAGP_VAR, i, j, k)*dvol
#endif

#ifdef DENS_VAR
               if (io_writeMscalarIntegrals) then
                  iSum = nGlobalSumProp
!!$                 do ivar=MASS_SCALARS_BEGIN,MASS_SCALARS_END
                  lsum(iSum + 1:iSum + NMASS_SCALARS) = &
                     lsum(iSum + 1:iSum + NMASS_SCALARS) + &
                     solnData(DENS_VAR, i, j, k)* &
                     solnData(MASS_SCALARS_BEGIN: &
                              MASS_SCALARS_END, i, j, k)*dvol
!!$                 end do
               end if

               ! Peak density
               maxDensLocal = max(maxDensLocal, solnData(DENS_VAR, i, j, k))
#endif
            end do
         end do
      end do
      call tileDesc%releaseDataPtr(solnData, CENTER)
      call tileDesc%releaseDataPtr(scratchData, SCRATCH)

      deallocate (cellVolumes)

      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)

#if defined(THORNADO)
   iSum = nGlobalSumProp - nOffGridThornado
   !lsum(iSum+1:iSum+2*THORNADO_NMOMENTS) = rt_offGridFluxR
   lsum(iSum + 1) = rt_offGridFluxR(1)*4.0d0*PI/hc3
   lsum(iSum + 2) = rt_offGridFluxR(5)*4.0d0*PI/hc3/Erg
   lsum(iSum + 3) = rt_offGridFluxR(6)*4.0d0*PI/hc3/Erg
   lsum(iSum + 4) = rt_offGridFluxR(7)*4.0d0*PI/hc3/Erg
   lsum(iSum + 5) = rt_offGridFluxR(8)*4.0d0*PI/hc3/Erg
#endif

   ! Now the MASTER_PE sums the local contributions from all of
   ! the processors and writes the total to a file.

   call MPI_Reduce(lsum, gsum, nGlobalSumUsed, FLASH_REAL, MPI_SUM, &
        &                MASTER_PE, io_globalComm, error)

   lsum(1) = maxDensLocal
   call MPI_AllReduce(lsum(1), lsum(2), 1, FLASH_REAL, MPI_MAX, &
      &  io_globalComm, error)
   maxDensGlobal = lsum(2)
   sim_maxDens = maxDensGlobal

   if (.not. sim_postBounce) &
      call sim_detectBounce(sim_postBounce, sim_bounceTime)

   postBounceTime = 0.0
   if (sim_postBounce) postBounceTime = simTime - sim_bounceTime

   if (io_globalMe == MASTER_PE) then

      ! create the file from scratch if it is a not a restart simulation,
      ! otherwise append to the end of the file

      !No matter what, we are opening the file. Check to see if already there
      ioStat = 0
      open (funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
      if (ioStat .NE. 0) then
         !print *, 'FILE FOUND'
         open (funit, file=trim(io_statsFileName), position='APPEND')
      end if

      if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then

         write (funit, 10) &
            '1_time                     ', &
            '2_postBounceTime           ', &
            '3_mass                     ', &
            '4_x-momentum               ', &
            '5_y-momentum               ', &
            '6_z-momentum               ', &
            '7_E_total                  ', &
            '8_E_kinetic                ', &
            '9_E_internal               ', &
#if defined(THORNADO)
            '10_e-type_nu_lepton_number  ', &
            '11_NuE_energy              ', &
            '12_NuE_Bar_energy          ', &
#if THORNADO_NSPECIES == 6
            '13_NuM_energy              ', &
            '14_NuM_Bar_energy          ', &
            '15_NuT_energy              ', &
            '16_NuT_Bar_energy          ', &
            '17_Nu_lepton_num_off_grid  ', &
            '18_Nu_energy_off_grid      ', &
            '19_Nu_x-momentum_off_grid  ', &
            '20_Nu_y-momentum_off_grid  ', &
            '21_Nu_z-momentum_off_grid  ', &
#ifdef MAGP_VAR
            '22_MagEnergy               ', &
#endif
#else
            '13_Nu_lepton_num_off_grid  ', &
            '14_Nu_energy_off_grid      ', &
            '15_Nu_x-momentum_off_grid  ', &
            '16_Nu_y-momentum_off_grid  ', &
            '17_Nu_z-momentum_off_grid  ', &
#ifdef MAGP_VAR
            '18_MagEnergy               ', &
#endif
#endif !THORNADO_NSPECIES == 6
#endif !defined(THORNADO)
#ifdef MAGP_VAR
            '10_MagEnergy                ', &
#endif
            (msName(ivar), ivar=MASS_SCALARS_BEGIN, &
             min(MASS_SCALARS_END, &
                 MASS_SCALARS_BEGIN + nGlobalSumUsed - nGlobalSumProp - 1))

10       format(2x, 50(a25, :, 1X))

      else if (isFirst .EQ. 1) then
         write (funit, 11)
11       format('# simulation restarted')
      end if

      ! Write the global sums to the file.
      write (funit, 12) simtime, postBounceTime, gsum(1:nGlobalSumUsed)

12    format(1x, 50(es25.16e4, :, 1x))

      close (funit)          ! Close the file.

   end if

#ifdef USEBARS
   call MPI_Barrier(io_globalComm, error)
#endif

   !=============================================================================

   return

contains
   character(len=25) function msName(ivar)
      integer, intent(in) :: ivar
      character(len=25) :: str
      call Simulation_mapIntToStr(ivar, str, MAPBLOCK_UNK)
      msName = str
   end function msName
end subroutine IO_writeIntegralQuantities

