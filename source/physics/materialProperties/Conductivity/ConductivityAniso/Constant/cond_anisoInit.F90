!!****if* source/physics/materialProperties/Conductivity/ConductivityAniso/Constant/cond_anisoInit
!!
!! NAME
!!
!!  cond_anisoInit
!!
!! SYNOPSIS
!!
!!  call cond_anisoInit()
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
subroutine cond_anisoInit
  use cond_anisoData, ONLY: cond_constantParallel, cond_constantPerpendicular, cond_constantCross
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  call RuntimeParameters_get("cond_constantParallel", cond_constantParallel)
  call RuntimeParameters_get("cond_constantPerpendicular", cond_constantPerpendicular)
  call RuntimeParameters_get("cond_constantCross", cond_constantCross)

end subroutine cond_anisoInit

