!!****h* source/Driver/Driver_interface
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
!!  Driver_interface
!!
!! SYNOPSIS
!!
!!  use  Driver_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Driver module
!! that defines its public interfaces.
!!***
Module Driver_interface

   implicit none

# include "Simulation.h"
# include "constants.h"

   interface
      subroutine Driver_abort(errorMessage)
         implicit none
         character(len=*), intent(in) :: errorMessage
      end subroutine Driver_abort
   end interface

   interface
      subroutine Driver_dbgBreak()
      end subroutine Driver_dbgBreak
   end interface

   interface
      subroutine Driver_computeDt(nbegin, nstep, &
                                  simTime, dtOld, dtNew)
         implicit none
         integer, intent(IN) :: nbegin, nstep
         real, intent(IN) :: simTime    !! current simulation time
         real, intent(IN) :: dtOld      !! last time step we used
         real, intent(OUT):: dtNew      !! the new timestep we get. to be returned.
      end subroutine Driver_computeDt
   end interface

   interface
      subroutine Driver_driftBlock(src_file, src_line, blk, ptr, gds)
         implicit none
         character(len=*), intent(in) :: src_file
         integer, intent(in) :: src_line
         integer, intent(in) :: blk
         real, intent(in) :: ptr(:, :, :, :)
         integer, intent(in) :: gds
      end subroutine Driver_driftBlock
   end interface

   interface
      subroutine Driver_driftUnk(src_file, src_line, flags)
         implicit none
         character(len=*), intent(in) :: src_file
         integer, intent(in) :: src_line
         integer, intent(in), optional :: flags
      end subroutine Driver_driftUnk
   end interface

   interface
      subroutine Driver_evolveAll()
         implicit none
      end subroutine Driver_evolveAll
   end interface

   interface
      subroutine Driver_finalizeAll()
         implicit none
      end subroutine Driver_finalizeAll
   end interface

   interface
      subroutine Driver_finalizeSourceTerms(restart)
         implicit none
         logical, intent(in) :: restart
      end subroutine Driver_finalizeSourceTerms
   end interface

   interface
      subroutine Driver_getDt(dt)
         implicit none
         real, intent(out) :: dt
      end subroutine Driver_getDt
   end interface

   interface
      subroutine Driver_getElapsedWCTime(elapsedWCTime)
         implicit none
         real, intent(out) :: elapsedWCTime
      end subroutine Driver_getElapsedWCTime
   end interface

   interface
      subroutine Driver_getNStep(nstep)
         implicit none
         integer, intent(out) :: nstep
      end subroutine Driver_getNStep
   end interface

   interface
      subroutine Driver_getSimTime(simulationTime, simGeneration)
         implicit none
         real, intent(out) :: simulationTime
         integer, intent(out), OPTIONAL :: simGeneration
      end subroutine Driver_getSimTime
   end interface

   interface
      subroutine Driver_init()
         implicit none
      end subroutine Driver_init
   end interface

   interface
      subroutine Driver_notifyGridChange()
         implicit none
      end subroutine Driver_notifyGridChange
   end interface

   interface
      subroutine Driver_initAll()
         implicit none
      end subroutine Driver_initAll
   end interface

   interface
      subroutine Driver_initMaterialProperties()
         implicit none
      end subroutine Driver_initMaterialProperties
   end interface

   interface
      subroutine Driver_initNumericalTools()
      end subroutine Driver_initNumericalTools
   end interface

   interface
      subroutine Driver_initParallel()
      end subroutine Driver_initParallel
   end interface

   interface
      subroutine Driver_initSourceTerms(restart)
         implicit none
         logical, intent(in) :: restart
      end subroutine Driver_initSourceTerms
   end interface

   interface
      subroutine Driver_sendOutputData()
         implicit none
      end subroutine Driver_sendOutputData
   end interface

   interface
      subroutine Driver_sourceTerms(dt, pass)
         implicit none
         real, intent(IN)    :: dt
         integer, OPTIONAL, intent(IN):: pass
      end subroutine Driver_sourceTerms
   end interface

   interface
      subroutine Driver_verifyInitDt()
         implicit none
      end subroutine Driver_verifyInitDt
   end interface

   interface
      subroutine Driver_getTimeStamp(dateStr)
         implicit none
         character(len=40), intent(OUT)     :: dateStr
      end subroutine Driver_getTimeStamp
   end interface

   interface
      subroutine Driver_putTimeStamp(dateStr)
         implicit none
         character(len=40), intent(IN)     :: dateStr
      end subroutine Driver_putTimeStamp
   end interface

   interface
      subroutine Driver_checkMPIErrorCode(errorCode)
         implicit none
         integer, intent(IN) :: errorCode
      end subroutine Driver_checkMPIErrorCode
   end interface

   interface
      subroutine Driver_superTimeStep(dt, nuSTS, nstepSTS, nstepTotalSTS, dt_subSTS)
         implicit none
         real, intent(IN)    :: dt, nuSTS
         integer, intent(IN) :: nstepSTS, nstepTotalSTS
         real, intent(OUT)   :: dt_subSTS
      end subroutine Driver_superTimeStep
   end interface

   interface
      subroutine Driver_getMype(communicatorType, mype, axis)
         integer, INTENT(IN) :: communicatorType
         integer, INTENT(OUT) :: mype
         integer, optional, intent(IN) :: axis
      end subroutine Driver_getMype
   end interface

   interface
      subroutine Driver_getNumProcs(communicatorType, numProcs, axis)
         integer, INTENT(IN) :: communicatorType
         integer, INTENT(OUT) :: numProcs
         integer, optional, intent(IN) :: axis
      end subroutine Driver_getNumProcs
   end interface

   interface
      subroutine Driver_getComm(communicatorType, communicator, axis)
         integer, INTENT(IN) :: communicatorType
         integer, INTENT(OUT) :: communicator
         integer, optional, intent(IN) :: axis
      end subroutine Driver_getComm
   end interface

   interface
      subroutine Driver_setupParallelEnv()

         implicit none
      end subroutine Driver_setupParallelEnv
   end interface

   interface
      subroutine Driver_mpiThreadSupport(mpiThreadSupport)
         implicit none
         logical, intent(OUT) :: mpiThreadSupport
      end subroutine Driver_mpiThreadSupport
   end interface

   interface
      subroutine Driver_logMemoryUsage(callsite)
         implicit none
         character(len=*), intent(IN) :: callsite
      end subroutine Driver_logMemoryUsage
   end interface

   interface
      subroutine Driver_diagnostics(blockCount, blockList, dt)
         implicit none
         real, intent(IN) :: dt
         integer, intent(IN) :: blockCount
         integer, dimension(blockCount), intent(IN):: blockList
      end subroutine Driver_diagnostics
   end interface

   interface Driver_envGetScalar
      subroutine Driver_envGetScalarInt(name, value)
         implicit none
         character(len=*), intent(in)          :: name
         integer, intent(out)                  :: value
      end subroutine Driver_envGetScalarInt
   end interface

end Module Driver_interface
