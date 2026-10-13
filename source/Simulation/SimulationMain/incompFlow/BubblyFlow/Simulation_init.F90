!!****if* source/Simulation/SimulationMain/incompFlow/BubblyFlow/Simulation_init
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
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!   none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine Simulation_init()

   use Driver_interface, ONLY: Driver_getMype
   use Simulation_data, ONLY: sim_xMin, sim_yMin, &
                              sim_xMax, sim_yMax, &
                              sim_zMin, sim_zMax, &
                              sim_meshMe, sim_numBubbles, &
                              sim_xBubble, sim_yBubble, sim_zBubble

   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use HDF5

   implicit none

   integer(HID_T)                 :: file
   integer                        :: h5err
   integer(HID_T)                 :: dset
   integer(HSIZE_T), dimension(3) :: dims
 
   call Driver_getMype(MESH_COMM, sim_meshMe)

   call RuntimeParameters_get('xmin', sim_xMin)
   call RuntimeParameters_get('ymin', sim_yMin)
   call RuntimeParameters_get('zmin', sim_zMin)
   call RuntimeParameters_get('xmax', sim_xMax)
   call RuntimeParameters_get('ymax', sim_yMax)
   call RuntimeParameters_get('zmax', sim_zMax)
   call RuntimeParameters_get('sim_numBubbles', sim_numBubbles)

   allocate(sim_xBubble(sim_numBubbles))
   allocate(sim_yBubble(sim_numBubbles))
   allocate(sim_zBubble(sim_numBubbles))

   call h5open_f(h5err)
   call h5fopen_f("bubble_dist_hdf5", H5F_ACC_RDONLY_F, file, h5err) !H5F_ACC_RDONLY_F
   if (h5err < 0) call Driver_abort('Unable to open bubble distribution file')

   dims = (/sim_numBubbles, 1, 1/)

   call h5dopen_f(file, "xBubble", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read xBubble')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, sim_xBubble, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "yBubble", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read yBubble')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, sim_yBubble, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "zBubble", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read zBubble')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, sim_zBubble, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5fclose_f(file, h5err)
   call h5close_f(h5err)

end subroutine Simulation_init
