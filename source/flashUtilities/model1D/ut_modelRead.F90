!!****if* source/flashUtilities/Model1D/ut_modelRead
!!
!! NAME
!!
!!  ut_modelRead
!!
!!
!! SYNOPSIS
!!
!!  ut_modelRead(  
!!                    string(IN)  :: model_file
!!                    real(IN)    :: requestVars
!!                    integer(IN) :: requestSize
!!                    integer(OUT)       :: modelZones
!!                    real, pointer(OUT) :: modelVars
!!                    real(OUT)          :: modelCoords
!!               )
!!
!! DESCRIPTION
!!   
!!    Reads in an initial 1D model in FLASH format and
!!    returns it, along with coordinates, to calling
!!    routine.
!!
!! ARGUMENTS
!!
!!    model_file  : name of file to be read. Must be in "FLASH" 1D model format
!!    requestVars : vector of strings with LOWERCASE variable names that are desired
!!    requestSize : number of variables (size of requestVars)
!!    modelZones  : array size for each variable in modelVars. Will be determined automatically.
!!    modelVars   : output pointer to the variables. All variable NOT found in file will be zero
!!    modelCoords : array of coordinates for modelVars
!!
!! NOTES
!!
!!    The order in which you send in the requestVars is the 
!!    order that they will be returned.
!!
!!***

subroutine ut_modelRead(model_file, requestVars, requestSize, modelZones, modelVars, modelCoords)
   use Driver_interface, ONLY : Driver_getMype
   
   implicit none 
   
#include "constants.h"
   
   character (len=80), INTENT(IN) :: model_file 
   integer, INTENT(IN) :: requestSize
   character (len=4), INTENT(IN) :: requestVars(requestSize)
   integer, INTENT(OUT) :: modelZones
   real, allocatable, INTENT(OUT) :: modelVars(:,:)
   real, allocatable, INTENT(OUT) :: modelCoords(:)
   
   ! Local variable declarations
   character (len=256) :: current_line
   integer :: i,j, ipos, nvar_stored
   character (len=4), ALLOCATABLE :: var_labels(:)
   integer :: var_key (requestSize)
   real, allocatable :: var_temp(:)
   integer :: meshMe, stat
   
   call Driver_getMype(MESH_COMM, meshMe)
   
   ! open the file and read in the header
   call read_header()
   if (meshMe == MASTER_PE) print *, 'file opened ', model_file
   if (meshMe == MASTER_PE) print *,"read nvar_stored", nvar_stored
   if (requestSize .NE. nvar_stored .AND. meshMe == MASTER_PE) then
      print *, ' '
      print *, 'Warning: the number of variables stored in the'
      print *, 'input file is different than the number of'
      print *, 'variables in the current version of FLASH.'
      print *, ' '
      print *, 'The variables in the file that are also defined'
      print *, 'in FLASH will be read in.  Any missing variables'
      print *, 'will be initialized to zero'
      print *, ' '
   endif
   
   do j = 1, requestSize
      var_key(j) = NONEXISTENT
      
      do i = 1, nvar_stored      
         if (requestVars(j) == var_labels(i)) then
            var_key(j) = i
         endif      
      enddo
      
      if (var_key(j) == NONEXISTENT) then
         if(meshMe == MASTER_PE) then
            print *, 'Warning, variable: ', requestVars(j), &
            ' not found in the input file.'
            print *, 'initializing ', requestVars(j), ' to 0'
            print *, ' '
         endif
      endif
      
   enddo
   
   ! Read file once through to get number of lines
   
   modelZones = 0 
   do  
      read(2,*,iostat=stat)
      if (stat /= 0) exit 
      modelZones = modelZones + 1
   end do
   
   ! Close the model and deallocate var_labels 
   close(2) 
   DEALLOCATE(var_labels) 
   
   ! Allocate stuff
   allocate(modelVars(modelZones, requestSize))
   allocate(modelCoords(modelZones))
   allocate(var_temp(nvar_stored))
   
   ! Read model for real
   call read_header() 
   do i = 1, modelZones
      read(2,*,iostat=stat) modelCoords(i), (var_temp(j),j=1,nvar_stored)
      if (stat /= 0) exit
      ! put these in order, so model1d_var always contains the same variables
      ! in the same spots
      do j = 1, requestSize
         if (var_key(j) /= NONEXISTENT) then
            modelVars(i,j) = var_temp(var_key(j))
         else
            modelVars(i,j) = 0.0
         endif         
      enddo      
   enddo
   
   close(2)
   deallocate(var_labels)
   deallocate(var_temp)

   if (meshMe .EQ. MASTER_PE) then
     print *, 'file read completed'
     print *, modelZones, 'points read in'
  endif

contains
   
   subroutine read_header()
      
      open(unit=2,file=model_file,status='old')
      read (2,'(a80)') current_line   
      
      ! read in the number of variables per line
      read (2,'(a80)') current_line
      !write(*,*) current_line
      ipos = index(current_line,'=') + 1
      read (current_line(ipos:),*) nvar_stored      
      
      ! Allocate var_labels
      allocate(var_labels(nvar_stored))
      
      do i = 1, nvar_stored
         read (2,'(a4)') var_labels(i)
         call makeLowercase(var_labels(i))
      enddo
   
   end subroutine read_header 

end subroutine ut_modelRead