!!****f* source/physics/RadTrans/RadTrans_getDbgContext
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
!!  RadTrans_getDbgContext
!!
!! SYNOPSIS
!!
!!  call RadTrans_getDbgContext(type(RadTrans_dbgContext_t),intent(OUT)  :: context)
!!
!! DESCRIPTION
!!
!! Stub
!!
!! ARGUMENTS
!!
!!   context : the context
!!!!
!!
!!***

subroutine RadTrans_getDbgContext(context)
  use RadTrans_interface, ONLY: RadTrans_dbgContext_t
  implicit none
  type(RadTrans_dbgContext_t),intent(OUT) :: context
  
  context%step = 0
  context%group = 0

end subroutine RadTrans_getDbgContext

subroutine RadTrans_getDbgContextPtr(context)
  use RadTrans_interface, ONLY: RadTrans_dbgContext_t
  implicit none
  type(RadTrans_dbgContext_t),pointer :: context
  
  nullify(context)

end subroutine RadTrans_getDbgContextPtr
