subroutine sim_readProfile
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
  use Simulation_data, ONLY : sim_nProfile, &
       sim_profFileName,                    &
       sim_profileInitial,                  &
       sim_rProf, sim_vProf, sim_rhoProf, sim_pProf
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none

  integer :: i, j
  real    :: errIgnored

  call Timers_start("readProfile")
  ! Read in the solution file to use for interpolation
  open(unit=10, file=sim_profFileName, status="OLD", action="READ")
  do i=1, 17
     read(10,*)
  end do
  do i=1,sim_nProfile
     read(10,*) (sim_profileInitial(j,i+1),j=1,4)
  end do
  close(10)


  ! We have left one space for an initiial extra point. Now fill it it.
  sim_profileInitial(1,1) = 0.0 ! add a first point at r=0.0
  sim_profileInitial(2,1) = 0.0 ! velocity is zero at 0.0
  call ut_polint(sim_profileInitial(1,2:3),sim_profileInitial(3,2:3), 2, 0.0, &
                 sim_profileInitial(3,1), errIgnored) ! density
  call ut_polint(sim_profileInitial(1,2:3),sim_profileInitial(4,2:3), 2, 0.0, &
                 sim_profileInitial(4,1), errIgnored) ! pressure
  call Timers_stop("readProfile")

end subroutine sim_readProfile
