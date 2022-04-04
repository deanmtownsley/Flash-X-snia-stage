!!***if* source/Simulation/SimulationForcing/incompFlow/Heater/sim_heaterRead
!!
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
!!
!!***
subroutine sim_heaterRead(heaterID, heaterFile)

   use sim_heaterData, ONLY: sim_nucSeedRadius, sim_heaterType, sim_heaterInfo
   use Driver_interface, ONLY: Driver_abort
   use HDF5

   implicit none
   integer, intent(in)          :: heaterID
   character(len=*), intent(in) :: heaterFile

   !-------------------------------------------------------------------------------
   integer(HID_T)                 :: file
   integer                        :: h5err
   integer(HID_T)                 :: dset
   integer(HSIZE_T), dimension(3) :: dims
   type(sim_heaterType), pointer  :: heater

   !--------------------------------------------------------------------------------
   heater => sim_heaterInfo(heaterID)

   call h5open_f(h5err)
   call h5fopen_f(trim(heaterFile), H5F_ACC_RDONLY_F, file, h5err) !H5F_ACC_RDONLY_F
   if (h5err < 0) call Driver_abort('Unable to open heater file')

   dims = (/1, 1, 1/)

   call h5dopen_f(file, "heater/xMin", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/xMin')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%xMin, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/xMax", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/xMax')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%xMax, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/yMin", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/yMin')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%yMin, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/yMax", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/yMax')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%yMax, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/zMin", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/zMin')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%zMin, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/zMax", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/zMax')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%zMax, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/wallTemp", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/wallTemp')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%wallTemp, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/advAngle", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/advAngle')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%advAngle, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/rcdAngle", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/rcdAngle')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%rcdAngle, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/velContact", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/velContact')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%velContact, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "heater/nucWaitTime", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read heater/nucWaitTime')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%nucWaitTime, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, 'site/num', dset, h5err) ! dset handle to database required
   if (h5err < 0) call Driver_abort('Unable to read site/num')
   call h5dread_f(dset, H5T_NATIVE_INTEGER, heater%numSites, dims, h5err)
   call h5dclose_f(dset, h5err)

   allocate (heater%xSite(heater%numSites))
   allocate (heater%ySite(heater%numSites))
   allocate (heater%zSite(heater%numSites))
   allocate (heater%siteRadii(heater%numSites))
   allocate (heater%siteIsAttachedCurr(heater%numSites))
   allocate (heater%siteIsAttachedPrev(heater%numSites))
   allocate (heater%siteTimeStamp(heater%numSites))

   dims = (/heater%numSites, 1, 1/)

   call h5dopen_f(file, "site/x", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read site/x')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%xSite, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "site/y", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read site/y')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%ySite, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "site/z", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read site/z')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%zSite, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "init/radii", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read init/radii')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, heater%siteRadii, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5fclose_f(file, h5err)
   call h5close_f(h5err)

   heater%siteIsAttachedCurr = .false.
   heater%siteIsAttachedPrev = .false.
   heater%siteTimeStamp = 0.
   heater%seedRadius = sim_nucSeedRadius
   heater%seedHeight = heater%seedRadius*cos(heater%rcdAngle*acos(-1.0)/180)

   return
end subroutine sim_heaterRead
