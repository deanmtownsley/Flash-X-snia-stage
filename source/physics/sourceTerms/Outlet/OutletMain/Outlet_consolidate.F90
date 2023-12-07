!!***if* source/physics/sourceTerms/Outlet/OutletMain/Outlet_consolidate
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  Outlet_consolidate
!!
!!
!! SYNOPSIS
!!  Outlet_consolidate()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine Outlet_consolidate()

   use Outlet_data, ONLY: out_QOut, out_QAux, out_volOut, out_volAux, &
                          out_QOutLiq, out_QOutGas, out_QAuxLiq, out_QAuxGas, &
                          out_volOutLiq, out_volOutGas, out_volAuxLiq, out_volAuxGas, out_meshMe

   implicit none
   include "Flashx_mpi.h"

   integer :: ierr

   out_QOut = 0.
   out_QOutLiq = 0.
   out_QOutGas = 0.

   out_volOut = 0.
   out_volOutLiq = 0.
   out_volOutGas = 0.

#ifdef OUTLET_PHASED

   call MPI_Allreduce(out_QAuxLiq, out_QOutLiq, (HIGH-LOW+1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(out_QAuxGas, out_QOutGas, (HIGH-LOW+1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(out_volAuxLiq, out_volOutLiq, (HIGH-LOW+1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(out_volAuxGas, out_volOutGas, (HIGH-LOW+1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   out_QOutLiq = out_QOutLiq/(out_volOutLiq+1e-13)
   out_QOutGas = out_QOutGas/(out_volOutGas+1e-13)

   if (out_meshMe .eq. MASTER_PE) then
      write (*, *) 'Outlet Liq Velocity LOW,', out_QOutLiq(LOW, :)
      write (*, *) 'Outlet Liq Velocity HIGH,', out_QOutLiq(HIGH, :)
      write (*, *) '--------------------------------------------------------'
      write (*, *) 'Outlet Gas Velocity LOW,', out_QOutGas(LOW, :)
      write (*, *) 'Outlet Gas Velocity HIGH,', out_QOutGas(HIGH, :)
   end if

#else

   call MPI_Allreduce(out_QAux, out_QOut, (HIGH-LOW+1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(out_volAux, out_volOut, (HIGH-LOW+1)*MDIM, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   out_QOut = out_QOut/(out_volOut+1e-13)

   if (out_meshMe .eq. MASTER_PE) then
      write (*, *) 'Outlet Velocity LOW,', out_QOut(LOW, :)
      write (*, *) 'Outlet Velocity HIGH,', out_QOut(HIGH, :)
   end if

#endif

   out_QAux = 0.
   out_QAuxLiq = 0.
   out_QAuxGas = 0.

   out_volAux = 0.
   out_volAuxLiq = 0.
   out_volAuxGas = 0.

end subroutine Outlet_consolidate
