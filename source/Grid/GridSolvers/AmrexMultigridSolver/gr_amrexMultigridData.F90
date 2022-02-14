!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexMultigridData
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
!!  gr_amrexMultigridData
!!
!! SYNOPSIS
!!  use gr_amrexMultigridData
!!
!! DESCRIPTION
!!
!!***

#include "Simulation.h"

module gr_amrexMultigridData
  
  use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t
  use amrex_multifab_module, ONLY : amrex_multifab
  use amrex_boxarray_module,     ONLY : amrex_boxarray
  use amrex_distromap_module,     ONLY : amrex_distromap
  use amrex_geometry_module, ONLY : amrex_geometry
  use amrex_fort_module,     ONLY : amrex_real
  use amrex_multigrid_module, ONLY : amrex_multigrid
  use amrex_poisson_module, ONLY: amrex_poisson
  
  implicit none

  ! parameters

  integer, save :: gr_amrexMG_maxLevel
  integer, save :: gr_amrexMG_ref_ratio = 2
  integer, save :: gr_amrexMG_n_cell
  integer, save :: gr_amrexMG_max_grid_size

  logical, save :: gr_amrexMG_composite_solve = .FALSE.

  ! prob_type 1 here is Poisson with homogeneous Dirichlet boundary.
  ! prob_type 2 here is ABecLaplacian with homogeneous Neumann boundary.
  integer, save :: gr_amrexMG_prob_type = 1

  integer, save :: gr_amrexMG_verbose = 2
  integer, save :: gr_amrexMG_bottom_verbose = 0
  integer, save :: gr_amrexMG_max_iter = 100
  integer, save :: gr_amrexMG_max_fmg_iter = 0
  integer, save :: gr_amrexMG_linop_maxorder = 2
  logical, save :: gr_amrexMG_agglomeration = .true.
  logical, save :: gr_amrexMG_consolidation = .true.

  ! data
  type(amrex_geometry), allocatable, save :: gr_amrexMG_geom(:)
  type(amrex_boxarray), allocatable, save :: gr_amrexMG_ba(:)
  type(amrex_distromap), allocatable, save :: gr_amrexMG_dm(:)

  type(amrex_multifab), allocatable, save :: gr_amrexMG_solution(:)
  type(amrex_multifab), allocatable, save :: gr_amrexMG_rhs(:)
  type(amrex_multifab), allocatable, save :: gr_amrexMG_exact_solution(:)
  type(amrex_multifab), allocatable, save :: gr_amrexMG_acoef(:)
  type(amrex_multifab), allocatable, save :: gr_amrexMG_bcoef(:)
  real(amrex_real), save :: gr_amrexMG_ascalar, gr_amrexMG_bscalar

  type(amrex_poisson) :: gr_amrexMG_poisson
  type(amrex_multigrid) :: gr_amrexMG_multigrid

  real(amrex_real), save :: gr_amrexMG_Tol
  
end module gr_amrexMultigridData
