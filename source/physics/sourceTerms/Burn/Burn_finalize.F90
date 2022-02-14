!!****f* source/physics/sourceTerms/Burn/Burn_finalize
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
!!  Burn_finalize
!!
!!
!! SYNOPSIS
!! 
!!  call Burn_finalize()
!!
!!  
!! DESCRIPTION
!!
!!  Finalizes the Burn module.
!!
!! NOTES
!!  
!!  There is no implementation that does anything.
!!
!!  The NSE arrays used with parametricBurn are deallocated by
!!  NSE_finalize, which should be called directly frm Driver_finalizeAll.
!!
!!***


subroutine Burn_finalize()


  implicit none

  return

end subroutine Burn_finalize
