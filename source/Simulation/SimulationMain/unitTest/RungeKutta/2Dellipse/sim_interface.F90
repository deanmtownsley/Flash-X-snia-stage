!!****ih* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/sim_interface
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
!!  sim_interface
!!
!! SYNOPSIS
!!
!!  use sim_interface
!!
!! DESCRIPTION
!!
!!  The interface for the current simulation unit.
!!
!!***

Module sim_interface

  interface
    subroutine sim_calculateInitialData ()
    end subroutine sim_calculateInitialData
  end interface

  interface
     function sim_distancePoint2Ellipse2Dxy (a, b, Rot, Cx, Cy, Px, Py)
       real, intent (in)  :: a, b
       real, intent (in)  :: Rot
       real, intent (in)  :: Cx, Cy
       real, intent (in)  :: Px, Py
       real               :: sim_distancePoint2Ellipse2Dxy (1:2)
     end function sim_distancePoint2Ellipse2Dxy
  end interface

  interface
    function sim_ODEfunction (t,y)
      real, intent (in) :: t
      real, intent (in) :: y (:)
      real              :: sim_ODEfunction (1:size (y))
    end function sim_ODEfunction
  end interface

end Module sim_interface
