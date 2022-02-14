!!****if* source/monitors/Debugger/DebuggerMain/Debugger_data
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
!!  Debugger_data
!!
!! SYNOPSIS
!!
!!  use Debugger_data
!!
!! DESCRIPTION
!!
!!  Holds the data needed by the Debugger unit
!!
!!***

module Debugger_data
  implicit none
  logical, save :: dbg_doHeapCheck
  integer, save :: dbg_globalMe
end module Debugger_data
