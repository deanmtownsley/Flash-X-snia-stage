! Edited by Nethra Rajavel to work with newest version of mesa and Flash 4.6

module mesa_flash_nucreact
!   use chem_def
!   use chem_lib
!   use eos_def
!   use eos_lib
!   use net_def
!   use net_lib
   use rates_def, only: extended_screening
!   use rates_lib, only: rates_init
!   use const_def, only: Qconv, secyer, kerg, avo
!   use const_lib
!   use utils_lib
!   use mtx_def
!   use crlibm_lib

   implicit none

   ! indexes for abundances in array passed to mesa net
   ! a few for convenience
   !integer :: He4_mi,  C12_mi, Si28_mi, Ni56_mi
   integer, dimension(:), pointer :: net_iso, chem_id, which_rates

   integer :: net_handle, eos_handle, num_nuclides, num_reactions

   ! various burner parameters, some initialized in the init function
   double precision, dimension(:), pointer :: rate_factors
   ! TODO just take out decsol switching if we're only going to use one option
   integer :: decsol_switch = 5000, decsol_choice
   character (len=32) :: small_mtx_decsol, large_mtx_decsol
   double precision :: weak_rate_factor = 1
   integer :: screening_mode = extended_screening, max_steps = 1000
   double precision :: theta_e_for_graboske_et_al = 1
   ! MESA test problem uses eps=1e-6, flashaprox13 uses 1e-5
   double precision :: stptry = 1d-16, eps = 1d-5, odescal = 1d-6
   logical :: okay_to_reuse_rate_screened = .true., use_pivoting = .true., trace = .false., burn_dbg = .false.

   double precision, dimension(:), pointer :: mass_excesses
   integer, dimension(:), pointer :: aion

contains

subroutine init_mesa_internals( net_def_filename )

   use chem_def, only: chem_isos
   use chem_lib
   use eos_lib
   use rates_def
   use rates_lib, only: rates_init
   use const_lib
   use math_lib !previously crlibm_lib
   use net_def, only: Net_General_Info
   use net_lib
   use mtx_lib, only: decsol_option

   implicit none

   character (len=*), intent(in) :: net_def_filename
   integer:: ierr, i
   ! note there is also a net_info living in net_def
   type (Net_General_Info), pointer  :: mesaflash_net_info

   ! blank argument makes it use the environmental variable MESA_DIR
   call const_init('',ierr)
   if (ierr /= 0) then
      write(*,*) 'const_init failed'
      stop 1
   end if
   call math_init !previously crlibm_init
   call chem_init('isotopes.data', ierr)
   if (ierr /= 0) then
      write(*,*) 'chem_init failed'
      stop 1
   end if


   ! had to init eos even though it is not used (yet),
   ! otherwise some calls in the rates lib don't link
   ! since this is a linking error there is probably a missing 'use'
   ! statement somewhere, but I'm not going to go hunting
   ! TODO hunt it down
   !call eos_init('mesa', '', '', '', .true., ierr)
   call eos_init('', .true., ierr)
   if (ierr /= 0) then
      write(*,*) 'eos_init failed'
      stop 1
   end if


   !call rates_init('reactions.list', '', 'rate_tables', .false., '','','', ierr)
   call rates_init('reactions.list', '', 'rate_tables', .false., .false., '','','', ierr)   
   if (ierr /= 0) then
      write(*,*) 'rates_init failed'
      stop 1
   end if


   ! -------------------------
   ! now initialize network we want and save various info about it

   call net_init(ierr)
   if (ierr /= 0) then
      write(*,*) 'net_init failed'
      stop 1
   end if
   net_handle = alloc_net_handle(ierr)
   eos_handle = alloc_eos_handle(ierr) !added by N.R
   if (ierr /= 0) then
      write(*,*) 'alloc_net_handle failed'
      stop 2
   end if
   call net_start_def(net_handle, ierr)
   if (ierr /= 0) then
      write(*,*) 'net_start_def failed'
      stop 2
   end if
   call read_net_file(net_def_filename, net_handle, ierr)
   if (ierr /= 0) then
      write(*,*) 'read_net_file failed'
      stop 2
   end if
   call net_finish_def(net_handle, ierr)
   if (ierr /= 0) then
      write(*,*) 'net_finish_def failed'
      stop 2
   end if

   call net_ptr(net_handle, mesaflash_net_info, ierr)
   if (ierr /= 0) then
      write(*,*) 'getting net pointer failed'
      stop 2
   end if

   num_nuclides = mesaflash_net_info% num_isos
   num_reactions = mesaflash_net_info% num_reactions

   allocate( rate_factors(num_reactions) )
   rate_factors(:) = 1

   large_mtx_decsol = 'klu'
   small_mtx_decsol = 'lapack'
   if ( num_nuclides >= decsol_switch) then
      decsol_choice = decsol_option(large_mtx_decsol, ierr)
      if (ierr /= 0 ) then
         write(*,*) 'error: unkown large_mtx_decsol: ' // trim(large_mtx_decsol)
         stop 2
      end if
   else
      decsol_choice = decsol_option(small_mtx_decsol, ierr)
      if (ierr /= 0 ) then
         write(*,*) 'error: unkown small_mtx_decsol: ' // trim(small_mtx_decsol)
         stop 2
      end if
   end if

   ! TODO may need to set more rates here see set_which_rates in burn test

   !allocate( which_rates(rates_reaction_id_max) )
   !which_rates(:) = rates_JR_if_available
   !call net_set_which_rates(net_handle, which_rates, ierr)
   !if (ierr /= 0) then
   !   write(*,*) 'net_set_which_rates failed'
   !   stop 2
   !end if
   ! TODO not sure if I need to do this...
   call net_setup_tables(net_handle, 'rate_tables', ierr)
   !                                             ^ cache_suffix
   if (ierr /= 0) then
      write(*,*) 'net_setup_tables failed'
      stop 2
   end if

   allocate(net_iso(num_chem_isos), chem_id(num_chem_isos))

   call get_chem_id_table(net_handle, num_nuclides, chem_id, ierr)
   if (ierr /= 0) then
      write(*,*) 'get_chem_id_table failed'
      stop 2
   end if
   call get_net_iso_table(net_handle, net_iso,ierr)
   if (ierr /= 0) then
      write(*,*) 'get_net_iso_table failed'
      stop 2
   end if

   ! initialize mass excesses for later computation of energy release
   allocate(aion(num_nuclides), mass_excesses(num_nuclides))
   do i=1,num_nuclides
      aion(i) = chem_isos% Z_plus_N(chem_id(i))
      mass_excesses(i) = get_mass_excess(chem_isos, chem_id(i))
   end do

   ! get some indices for convenience
   !He4_mi = net_iso( chem_get_iso_id('he4') )
   !C12_mi = net_iso( chem_get_iso_id('c12') )
   !Si28_mi = net_iso( chem_get_iso_id('si28') )
   !Ni56_mi = net_iso( chem_get_iso_id('ni56') )


end subroutine init_mesa_internals


! construct an index map from flash species index to mesa "isotope" index
subroutine fill_nuclide_index_map( num_flash_species, zlist, alist, map_fi_to_mi, ierr)

   use chem_lib, only: lookup_ZN
   integer, intent(in) :: num_flash_species
   integer, dimension(num_flash_species), intent(in) :: zlist, alist
   integer, dimension(num_flash_species), intent(out) :: map_fi_to_mi
   integer, intent(out) :: ierr

   integer :: fi ! flash index
   integer :: mi ! mesa index
   integer :: id

   if (num_flash_species /= num_nuclides) then
      write(*,*) 'num_flash_species /= num_nuclides'
      write(*,*) 'num_flash_species =', num_flash_species
      write(*,*) 'num_nuclides = ', num_nuclides
      ierr = 1
      return
   end if
   
   do fi=1,num_flash_species
      id = lookup_ZN( zlist(fi), alist(fi)-zlist(fi) )
      if (id <= 0 ) then
         ! failed to find this species in mesa
         write(*,*) 'failed to find mesa chem id for species with flash index ', fi, 'Z=', zlist(fi), 'A=', alist(fi)
         ierr = 1
         return
      end if
      mi = net_iso( lookup_ZN( zlist(fi), alist(fi)-zlist(fi) ) )
      if (mi <= 0 ) then
         ! failed to find this species in mesa
         write(*,*) 'failed to find mesa index for species with flash index ', fi
         ierr = 1
         return
      end if
      map_fi_to_mi(fi) = mi
   end do
   ierr = 0

end subroutine fill_nuclide_index_map
      
   
      




subroutine burn_const_dens_temp( dt, de_limit, dens, temp, abund )

   use rates_def, only: std_reaction_Qs, std_reaction_neuQs
   use net_lib, only: net_1_zone_burn!, net_1_zone_burn_work_size, net_work_size
   use math_lib, only: log10 
   use const_def, only: mev_amu
   use chem_def, only: num_categories
   ! time interval to integrate, density, and temperature
   double precision, intent(inout) :: dt
   double precision, intent(in) :: de_limit, dens, temp
   double precision, dimension(:), intent(inout) :: abund ! abundances in mesa order

   double precision :: xin(num_nuclides)
   double precision, dimension(:), pointer :: xout
   integer, parameter :: num_times = 1
   double precision, dimension(:), pointer :: times, log10Ts_f1, log10Rhos_f1, etas_f1
   double precision, dimension(:,:), pointer :: log10Ts_f, log10Rhos_f, etas_f
   double precision, dimension(:), pointer :: dxdt_source_term

   integer :: burn_work_array_len, net_work_array_len
   double precision, dimension(:), pointer :: burn_work_array, net_work_array
   integer :: nfcn, njac, nstep, naccpt, nrej, ierr
   double precision :: avg_eps_nuc, eps_neu_total
   double precision, target :: eps_nuc_categories(num_categories)
   double precision :: init_mass_excess, prev_substep_t, overstep_t, frac, norm
   double precision :: prev_substep_Y(num_nuclides), overstep_Y(num_nuclides)
   double precision :: over_mass_excess, prev_mass_excess
   integer :: i

   ! repack abundances for mesa net
   xin(:)=abund(:)

   ! setup temperature, density
   allocate( times(num_times), log10Ts_f1(4*num_times), log10Rhos_f1(4*num_times), etas_f1(4*num_times) )
   times(1) = dt

   log10Ts_f(1:4,1:num_times) => log10Ts_f1(1:4*num_times)
   log10Rhos_f(1:4,1:num_times) => log10Rhos_f1(1:4*num_times)
   etas_f(1:4,1:num_times) => etas_f1(1:4*num_times)

   log10Ts_f(1,1) = log10(temp)
   log10Ts_f(2:4,1) = 0
   log10Rhos_f(1,1) = log10(dens)
   log10Rhos_f(2:4,1) = 0
   etas_f(1:4,1) = 0

   dxdt_source_term => null()

   !burn_work_array_len = net_1_zone_burn_work_size(net_handle,ierr)
   !if (ierr /=0) then
   !   write(*,*) 'Unable to get burn work size'
   !   stop 1
   !endif
   !net_work_array_len = net_work_size(net_handle,ierr)
   !if (ierr /=0) then
   !   write(*,*) 'Unable to get net work size'
   !   stop 1
   !endif
   !allocate( net_work_array(net_work_array_len), burn_work_array(burn_work_array_len) )
   allocate( xout(num_nuclides) )
   do i=1,num_nuclides
      prev_substep_Y(i) = xin(i)/aion(i)
   end do
   init_mass_excess = dot_product(mass_excesses, prev_substep_Y)
   prev_substep_t = 0.0
   stptry = dt
   call net_1_zone_burn( &
     net_handle, eos_handle, num_nuclides, num_reactions, 0d0, dt, xin, &
     num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
     rate_factors, weak_rate_factor, &
     std_reaction_Qs, std_reaction_neuQs, &
     screening_mode,&
     stptry, max_steps, eps, odescal, &
     use_pivoting, trace, burn_dbg, mesaflash_burn_finish_substep, &
     xout, eps_nuc_categories, avg_eps_nuc, eps_neu_total, &
     nfcn, njac, nstep, naccpt, nrej, ierr)

   if (ierr == 15) then
      ! stopped due to overstepping de_limit
      ! interpolate to correct final state
      prev_mass_excess = dot_product(mass_excesses, prev_substep_Y)
      over_mass_excess = dot_product(mass_excesses, overstep_Y)
      frac = (de_limit/mev_amu - (init_mass_excess-prev_mass_excess)) / (prev_mass_excess-over_mass_excess)
      if (frac < 0.0 .or. frac > 1.0) then
         write(*,*) 'error in energy limiter, de_limit not between steps'
      end if
      norm = 0.0
      do i=1,num_nuclides
         xout(i) = ((1.0-frac)*prev_substep_Y(i) + frac*overstep_Y(i))*aion(i)
         norm = norm + xout(i)
      end do
      do i=1,num_nuclides
         xout(i) = xout(i)/norm
      end do
      dt = (1.0-frac)*prev_substep_t + frac*overstep_t
   elseif (ierr /= 0) then
      write(*,*) 'error in net_1_zone_burn'
      stop 1
   end if

   !deallocate(burn_work_array, net_work_array)
   deallocate(times, log10Ts_f1, log10Rhos_f1, etas_f1)

   abund(:) = xout(:)

   deallocate( xout )
   
   contains

   ! callback function
   subroutine mesaflash_burn_finish_substep(step, time, Y, ierr)
      integer, intent(in) :: step
      double precision, intent(in) :: time, Y(:)
      integer, intent(out) :: ierr
      double precision :: tot_mass_excess, de
      integer :: i

      tot_mass_excess = dot_product(mass_excesses, Y)

      de = init_mass_excess - tot_mass_excess

      if (de <= de_limit/mev_amu) then
         ! just save result and return
         prev_substep_Y(:) = Y(:)
         prev_substep_t = time
         ierr = 0
      else
         ! have exceeded limit, store data and return error to force exit
         ! have to store the data because it isn't actually returned normally
         ! when exiting via this error condition
         ierr = 15
         overstep_Y(:) = Y(:)
         overstep_t = time
      end if
   end subroutine mesaflash_burn_finish_substep

end subroutine burn_const_dens_temp

end module mesa_flash_nucreact
