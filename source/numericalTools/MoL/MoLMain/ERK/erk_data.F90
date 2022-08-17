module erk_data

   implicit none

   character(len=:), allocatable, save :: erk_method

   integer, save :: erk_order
   integer, save :: erk_stages

   ! RK tableau
   real, dimension(:, :), allocatable, target, save :: erk_A
   real, dimension(:), allocatable, target, save :: erk_b
   real, dimension(:), allocatable, target, save :: erk_c

   ! Indices for intermediate RHS states
   integer, allocatable, save :: erk_K(:)

end module erk_data
