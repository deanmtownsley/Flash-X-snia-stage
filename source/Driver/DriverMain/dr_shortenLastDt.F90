!!****if* source/Driver/DriverMain/dr_shortenLastDt
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
!!  dr_shortenLastDt
!!
!! SYNOPSIS
!!
!!  call dr_shortenLastDt(real(INOUT):: dt,
!!                        real (IN)  :: simTime,
!!                        real (IN)  :: tMax,
!!                      logical(OUT) :: shorteningDt,
!!                      integer(IN)  :: dtPerStep)
!!
!! DESCRIPTION
!!
!!  At the last timestep, taking a full timestep will many times exceed
!!  the maximum time of the simulation. This routine adjusts the dt value
!!  for the last time step so that it aligns exactly with tmax
!!
!!
!! ARGUMENTS
!!
!!   dt           - delta t
!!   simTime      - evolution time
!!   tMax         - the maximum evolution time
!!   shorteningDt - logical argument indicating whether to shorten dt
!!   dtPerStep    - How many advances (each by dt) Does the Driver do
!!                  in one iteration? (e.g., 2 for the default Split
!!                  driver and 1 for the default Unsplit driver)
!!
!!***

subroutine dr_shortenLastDt(dt, simTime, tMax, shorteningDt, dtPerStep)
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Driver_data, ONLY: dr_globalMe, dr_shortenLastStepBeforeTMax

  implicit none

  real,intent(inout) :: dt
  real,intent(in) :: simTime, tMax
  logical,intent(out) :: shorteningDt
  integer,intent(in) :: dtPerStep  !2 for split driver, 1 otherwise

  real :: origDt, remainingTime

  shorteningDt = .FALSE.
  if (dr_shortenLastStepBeforeTMax) then
     if (simTime + dtPerStep*dt > tMax) then
        remainingTime = tMax - simTime
        origDt = dt
        dt = remainingTime / dtPerStep
99      format ('Shortening dt from',(1PG26.19),' to',(1PG26.19),' for reaching tmax=',(1PG26.19))
        if (dr_globalMe==0) print 99,origDt,dt,tMax
        call Logfile_stampMessage( 'Shortening the last timestep to reach tmax exactly:')
        shorteningDt = .TRUE.
     end if
  end if

end subroutine dr_shortenLastDt



