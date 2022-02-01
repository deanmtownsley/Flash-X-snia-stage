!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/Grid_solvePoisson
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
!!  NAME 
!!
!! Grid_solvePoisson
!!
!!  SYNOPSIS
!!
!!  call Grid_solvePoisson()
!!
!!
!!  DESCRIPTION 
!! This routine solves the Poisson equation from the 
!! Amrex Linear Solvers using the variables from Unk multifab
!! for rhs and unknown phi 
!!
!!
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (must be a cell centered variable)
!!  iSrc     - the index of the source variable (must be a cell centered variable)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!             Only the first 2*NDIM elements are significant. They are interpreted
!!             in the order (X left, X right, Y left, Y right, Z left, Z right).
!!             Valid values are:
!!               GRID_PDE_BND_PERIODIC (1)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous or constant Dirichlet)
!!               GRID_PDE_BND_NEUMANN (3) (homogeneous or constant Neumann)
!!
!!  bcValues - the values to boundary conditions, currently not used (treated as 0)
!!  poisfact      - scaling factor to be used in calculation, currently not used (treated as 1)
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Currently, solver only works for periodic, Neumann and Dirichlet BCs
!!  Other BCs to be implemented later
!!  Relative and absolute tolerences for multigrid sovle - 1.e-10, 0.0
!!
!!***

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)
    use Timers_interface, ONLY : Timers_start, Timers_stop
    use Driver_interface, ONLY : Driver_abort
    use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
         GRID_PDE_BND_NEUMANN,   &
         GRID_PDE_BND_DIRICHLET
    use amrex_multigrid_module, ONLY : amrex_multigrid_build, amrex_multigrid_destroy
    use amrex_poisson_module, ONLY : amrex_poisson_build, amrex_poisson_destroy
    use amrex_lo_bctypes_module, ONLY : amrex_lo_periodic, amrex_lo_dirichlet, amrex_lo_neumann
    use amrex_amr_module, ONLY : amrex_geom
    use amrex_fort_module,     ONLY : amrex_real
    use gr_physicalMultifabs,  ONLY : unk
    !!
    use amrex_multifab_module, ONLY : amrex_multifab, amrex_multifab_destroy, amrex_multifab_build_alias
    use amrex_boxarray_module, ONLY : amrex_boxarray, amrex_boxarray_destroy
    use amrex_distromap_module, ONLY : amrex_distromap, amrex_distromap_destroy
    use Grid_data, only: gr_globalMe
    use gr_amrexMultigridData

    implicit none  
    integer, intent(in)    :: iSoln, iSrc
    integer, intent(in)    :: bcTypes(6)
    real, intent(in)       :: bcValues(2,6)
    real, intent(inout)    :: poisfact
    integer                :: amrexPoissonBcTypes(6)

    integer                :: i,ilev    
    real(amrex_real)       :: err

#include "Simulation.h"
#include "constants.h"   
  
    call Timers_start("Grid_solvePoisson")
    if(poisfact .ne. 1) then
      if(gr_globalMe .eq. MASTER_PE) print*,"[WARNING] [Grid_solvePoisson] Variable poisfact &
                                                                    was set not 1. It is being ignored and set to 1!!!"
      poisfact =1.
    endif

    !---------------------------------------------------------------------------------
    !-------1. SETUP AMREX MULTIFABS--------------------------------------------------
    !---------------------------------------------------------------------------------
    do ilev = 0, gr_amrexMG_maxLevel
        call amrex_multifab_build_alias(gr_amrexMG_solution(ilev), unk(ilev), iSoln, 1)
        call amrex_multifab_build_alias(gr_amrexMG_rhs(ilev)     , unk(ilev), iSrc , 1)
        call gr_amrexMG_solution(ilev)%setVal(0.0_amrex_real)
    end do

    gr_amrexMG_ba=gr_amrexMG_rhs%ba
    gr_amrexMG_dm=gr_amrexMG_rhs%dm

    !---------------------------------------------------------------------------------
    !-------2. SETUP MULTIGRID AND SOLVE-------------- -------------------------------
    !---------------------------------------------------------------------------------
    if(gr_amrexMG_composite_solve) then !Composite solve
       ! Build poisson object with the geometry amrex_geom, boxarray unk%ba  and distromap unk%dm
       call amrex_poisson_build(gr_amrexMG_poisson, amrex_geom(0:gr_amrexMG_maxLevel), gr_amrexMG_rhs%ba, gr_amrexMG_rhs%dm, &
                              metric_term=.false., agglomeration=gr_amrexMG_agglomeration, consolidation=gr_amrexMG_consolidation)
       
       call gr_amrexMG_poisson % set_maxorder(gr_amrexMG_linop_maxorder)

       !Select BCs to send to AMReX poisson solver
       do i=1,6
         select case (bcTypes(i))
         case (GRID_PDE_BND_PERIODIC)
            amrexPoissonBcTypes(i)=amrex_lo_periodic
         case (GRID_PDE_BND_NEUMANN)
            amrexPoissonBcTypes(i)=amrex_lo_neumann
         case (GRID_PDE_BND_DIRICHLET)
            amrexPoissonBcTypes(i)=amrex_lo_dirichlet
         case default
            call Driver_abort('BC not implemented for AMReX poisson solver!')
         end select
       end do
       call gr_amrexMG_poisson % set_domain_bc([amrexPoissonBcTypes(1),amrexPoissonBcTypes(3),amrexPoissonBcTypes(5)], &
                                             [amrexPoissonBcTypes(2),amrexPoissonBcTypes(4),amrexPoissonBcTypes(6)])

       do ilev = 0, gr_amrexMG_maxLevel
         !solution multifab's ghost cells at physical boundaries have been set to bc values.
         call gr_amrexMG_poisson % set_level_bc(ilev, gr_amrexMG_solution(ilev))
       end do

       call amrex_multigrid_build(gr_amrexMG_multigrid, gr_amrexMG_poisson)
       call gr_amrexMG_multigrid % set_verbose(gr_amrexMG_verbose)
       call gr_amrexMG_multigrid % set_bottom_verbose(gr_amrexMG_bottom_verbose)
       call gr_amrexMG_multigrid % set_max_iter(gr_amrexMG_max_iter)
       call gr_amrexMG_multigrid % set_max_fmg_iter(gr_amrexMG_max_fmg_iter)

       if(gr_globalMe .eq. MASTER_PE) print*, "Calling multigrid solve, maxlev", gr_amrexMG_maxLevel
       err = gr_amrexMG_multigrid % solve(gr_amrexMG_solution, gr_amrexMG_rhs, gr_amrexMG_Tol, 0.0_amrex_real)
       if(gr_globalMe .eq. MASTER_PE) print*, err
       call amrex_multigrid_destroy(gr_amrexMG_multigrid)
       call amrex_poisson_destroy(gr_amrexMG_poisson)

    else ! level-by-level solver
       do ilev = 0, gr_amrexMG_maxLevel
         call amrex_poisson_build(gr_amrexMG_poisson, [amrex_geom(ilev)], [gr_amrexMG_ba(ilev)], [gr_amrexMG_dm(ilev)], &
             metric_term=.false., agglomeration=gr_amrexMG_agglomeration, consolidation=gr_amrexMG_consolidation)
       
         do i=1,6
          select case (bcTypes(i))
           case (GRID_PDE_BND_PERIODIC)
           amrexPoissonBcTypes(i)=amrex_lo_periodic
          case (GRID_PDE_BND_NEUMANN)
           amrexPoissonBcTypes(i)=amrex_lo_neumann
          case (GRID_PDE_BND_DIRICHLET)
           amrexPoissonBcTypes(i)=amrex_lo_dirichlet
          case default
           call Driver_abort('Only periodic BC implemented for AMReX poissonsolver!')
          end select
         end do
         call gr_amrexMG_poisson %set_domain_bc([amrexPoissonBcTypes(1),amrexPoissonBcTypes(3),amrexPoissonBcTypes(5)],&
                                                [amrexPoissonBcTypes(2),amrexPoissonBcTypes(4),amrexPoissonBcTypes(6)])

         if (ilev > 0) then
           ! use coarse level data to set up bc at corase/fine boundary
           call gr_amrexMG_poisson % set_coarse_fine_bc(gr_amrexMG_solution(ilev-1), gr_amrexMG_ref_ratio)
         end if
         ! Note that to the linear solver, the level is ZERO.  In
         ! this problem, when lev > 0, solution(lev) is going to
         ! be ignored because fine level grids are completed
         ! surrounded by coarse level.  If fine level grids do touch
         ! phyical domain, the multifab must have bc values at
         ! physical boundaries stored in ghost cells.
         call gr_amrexMG_poisson % set_level_bc(0, gr_amrexMG_solution(ilev))

         call amrex_multigrid_build(gr_amrexMG_multigrid, gr_amrexMG_poisson)
         call gr_amrexMG_multigrid % set_verbose(gr_amrexMG_verbose)
         call gr_amrexMG_multigrid % set_bottom_verbose(gr_amrexMG_bottom_verbose)
         call gr_amrexMG_multigrid % set_max_iter(gr_amrexMG_max_iter)
         call gr_amrexMG_multigrid % set_max_fmg_iter(gr_amrexMG_max_fmg_iter)

         err = gr_amrexMG_multigrid % solve([gr_amrexMG_solution(ilev)], [gr_amrexMG_rhs(ilev)],gr_amrexMG_Tol, 0.0_amrex_real)
         if(gr_globalMe .eq. MASTER_PE) print*, err
         call amrex_poisson_destroy(gr_amrexMG_poisson)
         call amrex_multigrid_destroy(gr_amrexMG_multigrid)
       end do
    endif

    !---------------------------------------------------------------------------------
    !-------3. CLEANUP MULTIFABS---------------------------------------------
    !---------------------------------------------------------------------------------
    do ilev = 0, gr_amrexMG_maxLevel
     call amrex_multifab_destroy(gr_amrexMG_solution(ilev))
     call amrex_multifab_destroy(gr_amrexMG_rhs(ilev))
     call amrex_boxarray_destroy(gr_amrexMG_ba(ilev))
     call amrex_distromap_destroy(gr_amrexMG_dm(ilev))
    end do
     
    call Timers_stop("Grid_solvePoisson")

end subroutine Grid_solvePoisson
