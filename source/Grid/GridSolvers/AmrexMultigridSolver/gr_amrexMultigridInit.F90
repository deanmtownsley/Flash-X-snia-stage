subroutine gr_amrexMultigridInit()
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

    use Timers_interface, ONLY : Timers_start, Timers_stop
    use amrex_amr_module, ONLY : amrex_get_finest_level
    use amrex_fort_module, ONLY : amrex_real
    use gr_amrexMultigridData
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Grid_data, ONLY: gr_meshMe

#include "constants.h"

    call Timers_start("gr_multigridAmrexInit")

    call RuntimeParameters_get("gr_amrexMG_composite_solve",gr_amrexMG_composite_solve)
    call RuntimeParameters_get("gr_amrexMG_max_iter",gr_amrexMG_max_iter)
    call RuntimeParameters_get("gr_amrexMG_linop_maxorder",gr_amrexMG_linop_maxorder)
    call RuntimeParameters_get("gr_amrexMG_Tol",gr_amrexMG_Tol)

    if(gr_meshMe .eq. MASTER_PE) then
      write(*,*) 'gr_amrexMG_composite_solve   =',gr_amrexMG_composite_solve
      write(*,*) 'gr_amrexMG_max_iter          =',gr_amrexMG_max_iter
      write(*,*) 'gr_amrexMG_linop_maxorder    =',gr_amrexMG_linop_maxorder
      write(*,*) 'gr_amrexMG_Tol               =',gr_amrexMG_Tol
    end if 

    gr_amrexMG_maxLevel = amrex_get_finest_level()

!   Allocate space for multifab array storing phi (solution) and rhs
    allocate(gr_amrexMG_solution(0:gr_amrexMG_maxLevel))
    allocate(gr_amrexMG_rhs(0:gr_amrexMG_maxLevel))
    allocate(gr_amrexMG_ba(0:gr_amrexMG_maxLevel))
    allocate(gr_amrexMG_dm(0:gr_amrexMG_maxLevel))

    call Timers_stop("gr_multigridAmrexInit")

end subroutine gr_amrexMultigridInit
