!!****if* source/Grid/localAPI/gr_mpoleDumpMoments
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
!!  gr_mpoleDumpMoments
!!
!! SYNOPSIS
!!  gr_mpoleDumpMoments()
!!
!! DESCRIPTION
!!
!!  Utility routine to output the regular Moment_R array and the irregular
!!  Moment_I array to a text file. The information is written out to a file
!!  named basenm_momentDump.txt, where basenm is the runtime parameter for
!!  output file names. The file is appended at each time for each iteration.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!  
!! NOTES
!!
!!  At the end of each iteration time, the line --- finished present iteration ---
!!  is inserted to make post-processing easier. The original distinctive
!!  phrase "Chakka Khan Chakka Khan I feel for you this is the end of this timestep"
!!  has NOT been retained. Sorry, but I am open for discussion about including it
!!  again.
!!
!!***

subroutine gr_mpoleDumpMoments ()

  implicit none

  return
end subroutine gr_mpoleDumpMoments
