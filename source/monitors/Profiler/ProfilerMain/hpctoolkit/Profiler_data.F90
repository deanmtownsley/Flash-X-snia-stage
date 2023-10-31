!!****if* source/monitors/Profiler/ProfilerMain/Profiler_data
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
!!  Profiler_data
!!
!! SYNOPSIS
!!
!!  use Profiler_data
!!
!! DESCRIPTION
!!
!!  Holds the data needed by the Profiler unit
!!
!!***

module Profiler_data
   implicit none
   logical, save :: prf_profilerIsOn
   logical, save :: prf_evolutionOnly
   character(len=50), save :: prf_groupName
   integer, save :: prf_meshMe
end module Profiler_data
