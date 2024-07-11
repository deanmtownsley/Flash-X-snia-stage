! program to generate Config file list of species for flash for a mesa network definition

! Dean M. Townsley, 2015

! need to change to allow specification of nuclear network definition file

program generate_flash_species_config
   use net_def, only: get_net_ptr, Net_General_Info
   use chem_def, only: chem_isos
   use mesa_flash_nucreact, only: init_mesa_internals, net_handle

   type (Net_General_Info), pointer :: g
   integer :: ierr, i, id
   character (len=4096) :: net_def_filename

   if ( command_argument_count() /= 1 ) then
      write (*,*) 'Please specify one argument that is name of mesa net def file.'
      stop
   end if

   call get_command_argument( 1, value=net_def_filename, status=ierr)
   if (ierr/=0) then
      write (*,*) 'Net def file name too long.'
      stop
   endif

   call init_mesa_internals(trim(net_def_filename))

   ierr = 0
   call get_net_ptr(net_handle, g, ierr)
   if (ierr /=0) then
      write (*,*) 'error getting mesa net pointer'
      stop
   end if

   write (*,*) '# species list for Flash Config file from mesa net definition file ..'
   do i=1, g% num_isos
      id = g% chem_id(i)
      if ( id > 0 ) then
         write (*,*) 'SPECIES ', trim(chem_isos% name(id))
      else
         write (*,*) 'Error: nuclide with nonpositive id: ', id,'. What is it?'
         stop
      end if
   end do

end program


