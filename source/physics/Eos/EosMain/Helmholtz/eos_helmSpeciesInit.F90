!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!> @ingroup physics_Eos
!!
!! @brief  Initialize the Helmholtz EOS.  The table data is read in on processor 0
!!  and broadcast to all other processors during initialization.
!!  This routine first checks for a binary copy of the table (helm_table.bdat), 
!!  and then for the ASCII version (helm_table.dat).  If the binary table 
!!  file was not found, it is created by this routine for subsequent use.
!!

#ifdef DEBUG_ALL
#define DEBUG_EOS
#endif


#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

subroutine eos_helmSpeciesInit()

  use Eos_data, ONLY : eos_type, eos_meshMe, &
       eos_eintSwitch, eos_smallt
  use eos_helmData 
  use Driver_interface, ONLY : Driver_abort
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

#include "Flashx_mpi_implicitNone.fh"



  ! vector_eos.fh computes the vector length from nxb, nyb, nzb, so 
  ! this information must be provided

  integer:: unitEos =2
  integer :: i, j
  real :: tstp, dstp
  integer :: istat, ierr


  !get the runtime parameters

  call RuntimeParameters_get('smallt', eos_smallt)
  call RuntimeParameters_get('eos_tolerance', eos_tol)
  call RuntimeParameters_get('eos_maxNewton', eos_maxNewton)
  call RuntimeParameters_get('eos_coulombMult', eos_coulombMult)
  call RuntimeParameters_get('eos_useMultiSpecies',eos_useMultiSpecies)
  call RuntimeParameters_get("eos_forceConstantInput",eos_forceConstantInput)

#ifdef DEBUG_EOS
  print *, 'in eos_helmSpeciesInit'
#endif
  call RuntimeParameters_get('eos_coulombAbort', eos_coulombAbort)
#ifdef DEBUG_EOS
  print *, 'done with RuntimeParameters_get (eos_coulombAbort)'
#endif

#ifndef EINT_VAR
  if (eos_eintSwitch > 0.0) then
     call Driver_abort("[Eos_init] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif

  if (eos_meshMe==MASTER_PE) then
     open(unit=unitEos,file='helm_table.bdat',status='old',iostat=istat)
     close(unit=unitEos)
     print*,'about to open file'
     istat=1
     if (istat.ne.0) then
        write(*,*) '[Eos_init] Cannot open helm_table.bdat!'
        write(*,*) '[Eos_init] Trying old helm_table.dat!'
        
        open(unit=unitEos,file='helm_table.dat',status='old',iostat=istat)

        if (istat .ne. 0) then
           write(*,*) '[Eos_init] ERROR: opening helm_table.dat!'
           call Driver_abort("[Eos_init] ERROR: opening helm_table.dat")
        endif

        !..read the helmholtz free energy table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_table(i,j,eos_f),eos_table(i,j,eos_fd),eos_table(i,j,eos_ft),&
                   eos_table(i,j,eos_fdd),eos_table(i,j,eos_ftt), & 
                   eos_table(i,j,eos_fdt),eos_table(i,j, eos_fddt),eos_table(i,j,eos_fdtt),eos_table(i,j,eos_fddtt)
           enddo
        enddo

        !..read the pressure derivative with density table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_table(i,j,eos_dpdf),eos_table(i,j,eos_dpdfd),&
                   eos_table(i,j,eos_dpdft),eos_table(i,j,eos_dpdfdt)
           enddo
        enddo

        !..read the electron chemical potential table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_table(i,j,eos_ef),eos_table(i,j,eos_efd),&
                   eos_table(i,j,eos_eft),eos_table(i,j,eos_efdt)
           enddo
        enddo

        !..read the number density table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_table(i,j,eos_xf),eos_table(i,j,eos_xfd),&
                   eos_table(i,j,eos_xft),eos_table(i,j,eos_xfdt)
           enddo
        enddo

        !..close up the data file and write a message
        close(unitEos)

        !..dump binary version of table for later use
        istat = EOSIMAX*EOSJMAX*EOST
        call eos_writeHfet(istat, eos_table)

        !..read binary version of table
     else
        istat = EOSIMAX*EOSJMAX*EOST
        call eos_readHfet(istat,eos_table)
     endif
  endif

  !..broadcast to rest of processors
  istat = EOSIMAX*EOSJMAX*EOST
  call MPI_BCAST(eos_table,      istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  eos_tlo   = 3.0e0
  eos_thi   = 13.0e0
  tstp      = (eos_thi - eos_tlo)/float(EOSJMAX-1)
  eos_tstpi = 1.0e0/tstp
  eos_dlo   = -12.0e0
  eos_dhi   = 15.0e0
  dstp      = (eos_dhi - eos_dlo)/float(EOSIMAX-1)
  eos_dstpi = 1.0e0/dstp
  do j=1,EOSJMAX
     eos_temps(j,eos_t) = 10.0e0**(eos_tlo + (j-1)*tstp)
     do i=1,EOSIMAX
        eos_rhos(i,eos_d) = 10.0e0**(eos_dlo + (i-1)*dstp)
     enddo
  enddo


  !..store the temperature and density differences and their inverses 
  do j=1,EOSJMAX-1
     eos_temps(j,eos_dt)   = eos_temps(j+1,eos_t) - eos_temps(j,eos_t)
     eos_temps(j,eos_dtSqr)  = eos_temps(j,eos_dt)*eos_temps(j,eos_dt)
     eos_temps(j,eos_dtInv)  = 1.0e0/eos_temps(j,eos_dt)
     eos_temps(j,eos_dtSqrInv) = 1.0e0/eos_temps(j,eos_dtSqr)
  enddo
  do i=1,EOSIMAX-1
     eos_rhos(i,eos_dd)   = eos_rhos(i+1,eos_d) - eos_rhos(i,eos_d)
     eos_rhos(i,eos_ddSqr)  = eos_rhos(i,eos_dd)*eos_rhos(i,eos_dd)
     eos_rhos(i,eos_ddInv)  = 1.0e0/eos_rhos(i,eos_dd)
     eos_rhos(i,eos_ddSqrInv) = 1.0e0/eos_rhos(i,eos_ddSqr)
  enddo
  eos_type=EOS_HLM

  return
end subroutine eos_helmSpeciesInit
