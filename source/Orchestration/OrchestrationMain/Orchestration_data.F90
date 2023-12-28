module Orchestration_data
   implicit none


!!!!!----- MPI & OpenMP Information    
   integer, save :: or_globalComm
   integer, save :: or_globalMe

   integer, save :: or_meshComm
   integer, save :: or_meshMe
   integer, save :: or_meshNumProcs

end module Orchestration_data
! Local Variables:
! Mode: F90
! f90-program-indent: 3
! indent-tabs-mode: nil
! End:
