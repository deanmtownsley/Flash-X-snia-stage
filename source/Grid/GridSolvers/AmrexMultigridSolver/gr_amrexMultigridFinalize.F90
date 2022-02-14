subroutine gr_amrexMultigridFinalize()
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

    use Timers_interface, ONLY : Timers_start, Timers_stop
    use gr_amrexMultigridData

    call Timers_start("gr_multigridAmrexFinalize")
    deallocate(gr_amrexMG_solution)
    deallocate(gr_amrexMG_rhs)
    deallocate(gr_amrexMG_ba)
    deallocate(gr_amrexMG_dm)
    call Timers_stop("gr_multigridAmrexFinalize")

end subroutine gr_amrexMultigridFinalize
