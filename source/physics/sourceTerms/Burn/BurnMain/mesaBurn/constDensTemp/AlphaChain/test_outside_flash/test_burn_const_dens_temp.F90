! program to develop flash interface for mesa burn modules

! for now we will assume hard-coding of alpha-chain (aprox13 nuclides)

program test_burn_const_dens_temp
   use mesa_flash_nucreact, only: init_mesa_internals, burn_const_dens_temp, net_iso
   use chem_lib, only: chem_get_iso_id

   implicit none

   double precision :: dt, dens, temp, abund(13)
   integer :: C12_mi

   call init_mesa_internals('alpha_chain_plus_3he_co.net')


   ! mesa burn test problem:
   ! a simple integration for 1e4 seconds at dens=6e9, temp=9.6e9
   ! starting from pure Si28
  
   ! but flash's approx13 doesn't like those conditions

   ! my test problem for flash
   ! integrate at dens=1e8 and temp=6e9 for 1e-5 seconds
   ! starting from pure C12

   dens = 1.0d8
   temp = 6.0d9
   dt = 1e-5
   ! set abundance to pure Si
   ! use same ord
   C12_mi = net_iso( chem_get_iso_id('c12') )
   abund(:)= 0.0d0
   abund(C12_mi) = 1.0d0

   write(*,*) 'initial abundances', abund
   call burn_const_dens_temp( dt, dens, temp, abund )

   ! print final abundances
   write(*,*) 'final abundances', abund

end program test_burn_const_dens_temp


