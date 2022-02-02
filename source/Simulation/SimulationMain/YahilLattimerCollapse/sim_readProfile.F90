!!****if* source/Simulation/SimulationMain/YahilLattimerCollapse/sim_readProfile
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
!!  DESCRIPTION
!!
!!   Read in a Yahil profile which contains the position, density, and 
!!   velocity and return the user specified profile
!!
!!  AUTHOR AND DATE
!! 
!!   Ran Chu (rchu@vols.utk.edu)   5/13/2019
!!
!!***

subroutine sim_readProfile

  use Simulation_data, ONLY : sim_nProfile, &
       sim_profFileName,                    &
       sim_profileInitial,                  &
       sim_Prof_nVariable,                  &
       sim_Prof_iX, sim_Prof_iD, sim_Prof_iV, sim_Prof_iM, sim_Prof_iE, &
       sim_pCentral, sim_rhoCentral, sim_kappa, &
       sim_gammaInitial, sim_collapsetime, &
       sim_rProf, sim_vProf, sim_rhoProf, sim_eProf, sim_pProf, &
       sim_mProf, model_1d
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abort

  implicit none

  integer :: i, j, iostat
  real    :: errIgnored

  REAL    :: C_X, C_D, C_V, C_M
  REAL    :: GravitationalConstant = 6.673d-8

  call Timers_start("readProfile")

  ! open the file and read in the header
  open(unit=10, file=sim_profFileName, status="OLD", action="READ")
  do i=1, 1
     read(10,*)
  end do
  
!!$RC:  if (sim_meshMe == MASTER_PE) then
     print *, &
       'file opened, reading profile: ', sim_profFileName
!!$RC:  end if

  do i=1,sim_nProfile
     read(10,*,IOSTAT=iostat) (sim_profileInitial(j,i),j=1,4)

     if( iostat < 0 ) then
       write(*,'(A2,A)') '', 'Warning: File contains less than 2049 entries'
       exit
     else if ( iostat > 0 ) then
       call Driver_abort('Error: error reading file')
     end if
  end do
  close(10)

  sim_kappa = sim_pCentral / sim_rhoCentral**( sim_gammaInitial ) 

  C_X = SQRT( sim_kappa ) &
        * GravitationalConstant**( 0.5d0 * ( 1.0d0 - sim_gammaInitial ) ) &
        * sim_collapsetime**( 2.0d0 - sim_gammaInitial )

  C_D = 1.d0 / ( sim_collapsetime * sim_collapsetime * GravitationalConstant )

  C_V = SQRT( sim_kappa ) &
        * GravitationalConstant**( 0.5d0 * ( 1.0d0 - sim_gammaInitial ) ) &
        * sim_collapsetime**( 1.0d0 - sim_gammaInitial )

  C_M = sim_kappa**( 3/2 ) &
        * GravitationalConstant**( 0.5d0 * ( 1.0d0 - 3.0d0 * sim_gammaInitial ) ) &
        * sim_collapsetime**( 4.0d0 - 3.0d0 * sim_gammaInitial )

  sim_rProf   = sim_profileInitial(sim_Prof_iX,:) * C_X

  sim_rhoProf = sim_profileInitial(sim_Prof_iD,:) * C_D

  sim_vProf   = sim_profileInitial(sim_Prof_iV,:) * C_V

  sim_mProf   = sim_profileInitial(sim_Prof_iM,:) * C_M

  sim_eProf   = sim_kappa * sim_rhoProf ** sim_gammaInitial &
                / ( sim_gammaInitial - 1.d0 )

  sim_profileInitial(sim_Prof_iX,:) = sim_rProf
  sim_profileInitial(sim_Prof_iD,:) = sim_rhoProf
  sim_profileInitial(sim_Prof_iV,:) = sim_vProf
  sim_profileInitial(sim_Prof_iM,:) = sim_mProf
  sim_profileInitial(sim_Prof_iE,:) = sim_eProf

  sim_pProf = sim_kappa * sim_rhoProf**sim_gammaInitial 

  call Timers_stop("readProfile")

end subroutine sim_readProfile
