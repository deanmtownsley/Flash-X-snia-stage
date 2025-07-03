!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant-diff/Conductivity_init
!!
!! NAME
!!
!!  Conductivity_init
!!
!! SYNOPSIS
!!
!!  call Conductivity_init()
!!                          
!!                         
!!
!! DESCRIPTION
!!
!! Initializes the data in Conductivity data to a constant value "diff_constant"
!! from RuntimeParameters
!!
!! ARGUMENTS
!!
!!  none
!!
!!
!!***

subroutine Conductivity_init()

  use Conductivity_data, ONLY : cond_diffConstant
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none


  ! Everybody should know this
  call cond_commonInit()

  
  call RuntimeParameters_get("diff_constant", cond_diffConstant)

end subroutine Conductivity_init
