!!****h* source/physics/IncompNS/IncompNS_interface
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!
!!  IncompNS_interface
!!
!! SYNOPSIS
!!
!!  use IncompNS_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Incompressible Navier-Stokes (INS)
!! module that defines its public interfaces.
!!
!!***

Module IncompNS_interface

  implicit none

#include "Simulation.h"
#include "constants.h"

  interface ! IncompNS_computeDt

    subroutine IncompNS_computeDt(ins_mindt,ins_minloc)
      implicit none
      real, intent(INOUT) :: ins_mindt
      integer, intent(INOUT) :: ins_minloc(5)
    end subroutine IncompNS_computeDt

  end interface


  interface !IncompNS
    subroutine IncompNS (timeEndAdv,  dt, dtOld, sweepOrder)
      implicit none 
      integer, INTENT(IN) :: sweepOrder
      real,    INTENT(IN) :: timeEndAdv, dt, dtOld  
    end subroutine IncompNS
  end interface


  interface !IncompNS_init
    subroutine IncompNS_init (restart)
      implicit none
      logical, intent(IN) :: restart
    end subroutine IncompNS_init
  end interface


  interface !IncompNS_finalize
    subroutine IncompNS_finalize ()
      implicit none
    end subroutine IncompNS_finalize
  end interface


  interface !IncompNS_sendOutputData
    subroutine IncompNS_sendOutputData ()
      implicit none
    end subroutine IncompNS_sendOutputData
  end interface

  interface
    subroutine IncompNS_stats()
      implicit none
    end subroutine
  end interface

  interface
     subroutine IncompNS_statsIOExport(expt_flag)
       implicit none
       logical, intent(in) :: expt_flag
     end subroutine IncompNS_statsIOExport
  end interface

  interface
     subroutine IncompNS_velomgToCenter()
      end subroutine IncompNS_velomgToCenter
  end interface

  interface IncompNS_getScalarProp
     subroutine IncompNS_getScalarPropReal(name, value)
     implicit none
     character(len=*), intent(in)  :: name
     real, intent(out)             :: value
     end subroutine IncompNS_getScalarPropReal

     subroutine IncompNS_getScalarPropLogical(name, value)
     implicit none
     character(len=*), intent(in)  :: name
     logical, intent(out)          :: value
     end subroutine IncompNS_getScalarPropLogical
  end interface
  
  interface
     subroutine IncompNS_getVectorProp(name, value)
     implicit none
     character(len=*), intent(in)         :: name
     real, dimension(MDIM), intent(out)   :: value
     end subroutine IncompNS_getVectorProp
  end interface

  interface
     subroutine IncompNS_setVectorProp(name, value)
     implicit none
     character(len=*), intent(in)         :: name
     real, dimension(MDIM), intent(in)    :: value
     end subroutine IncompNS_setVectorProp
  end interface

  interface
     subroutine IncompNS_getGridVar(name, value)
     implicit none
     character(len=*), intent(in)  :: name
     integer, intent(out)          :: value
     end subroutine IncompNS_getGridVar
  end interface

  interface
     subroutine IncompNS_predictor(dt)
     implicit none
     real,    INTENT(IN) :: dt
     end subroutine IncompNS_predictor
  end interface

  interface
     subroutine IncompNS_solvePoisson(dt)
     implicit none
     real,    INTENT(IN) :: dt
     end subroutine IncompNS_solvePoisson
  end interface

  interface
     subroutine IncompNS_corrector(dt)
     implicit none
     real,    INTENT(IN) :: dt
     end subroutine IncompNS_corrector
  end interface

  interface
     subroutine IncompNS_divergence()
     implicit none
     end subroutine IncompNS_divergence
  end interface

  interface
     subroutine IncompNS_advection()
     implicit none
     end subroutine IncompNS_advection
  end interface

  interface
     subroutine IncompNS_diffusion()
     implicit none
     end subroutine IncompNS_diffusion
  end interface

  interface
     subroutine IncompNS_indicators()
     implicit none
     end subroutine IncompNS_indicators
  end interface

  interface
     subroutine IncompNS_reInitGridVars()
     implicit none
     end subroutine IncompNS_reInitGridVars
  end interface

end Module IncompNS_interface
