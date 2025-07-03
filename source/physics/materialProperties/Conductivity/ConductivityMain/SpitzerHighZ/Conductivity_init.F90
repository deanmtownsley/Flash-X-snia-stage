!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/Conductivity_init
!!
!! NAME
!!
!!  Conductivity_init
!!
!! SYNOPSIS
!!
!!  call Conductivity_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  none
!!
!!
!!***
subroutine Conductivity_init

  implicit none

  call cond_commonInit()

end subroutine Conductivity_init

