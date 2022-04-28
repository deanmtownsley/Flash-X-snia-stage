!!****if* source/IO/IOMain/IO_checkForPlot
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
!!    IO_checkForPlot
!!
!! SYNOPSIS
!!
!!    IO_checkForPlot(logical(out) : wrotePlot)
!!
!! DESCRIPTION
!!
!!   This routine sets the argument to true when called on a cycle
!!   following a plot.
!!
!! ARGUMENTS
!!
!!  wrotePlot : set to true if plot written, otherwise false
!!
!!***
subroutine IO_checkForPlot(wrotePlot)
  use IO_data, ONLY: io_wrotePlot
  implicit none
  
  logical, intent(out) :: wrotePlot

  wrotePlot = io_wrotePlot

end subroutine IO_checkForPlot
