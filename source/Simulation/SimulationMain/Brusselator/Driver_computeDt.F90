subroutine Driver_computeDt(nbegin, nstep, simTime, dtOld, dtNew)
   use Simulation_data, only: sim_k

   implicit none

   integer, intent(in)  :: nbegin, nstep
   real, intent(in)  :: simTime    !! current simulation time
   real, intent(in)  :: dtOld      !! last time step we used
   real, intent(out) :: dtNew      !! the new timestep we get. to be returned.

   dtNew = 0.1d0/2d0**sim_k
end subroutine Driver_computeDt
