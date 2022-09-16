!!****ih* source/flashUtilities/Model1D/ut_model_interface
!!
!! NAME
!!
!!  ut_model_interface
!!
!! SYNOPSIS
!!
!!  use ut_model_interface
!!
!! DESCRIPTION
!!
!!  Interface module for 1D model reading utilities.
!!
!!***
module ut_modelInterface
  implicit none      
  
  interface
  subroutine ut_modelRead(model_file, requestVars, requestSize, modelZones, modelVars, modelCoords)
    implicit none     
    character (len=80), INTENT(IN) :: model_file 
    integer, INTENT(IN) :: requestSize
    character (len=4), INTENT(IN) :: requestVars(requestSize)
    integer, INTENT(OUT) :: modelZones
    real, allocatable, INTENT(OUT) :: modelVars(:,:)
    real, allocatable, INTENT(OUT) :: modelCoords(:)
  end subroutine ut_modelRead
  end interface
  
end module ut_modelInterface
