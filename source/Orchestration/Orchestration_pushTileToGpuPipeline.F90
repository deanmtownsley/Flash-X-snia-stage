!! @copyright Copyright 2024 UChicago Argonne, LLC and contributors
subroutine Orchestration_pushTileToGpuPipeline(prototype_Cptr, nThreads, &
                                            tileCInfo)
   use iso_c_binding, ONLY : C_PTR

   use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t

   implicit none

   type(C_PTR),                            intent(IN) :: prototype_Cptr
   integer,                                intent(IN) :: nThreads
   type(Orchestration_tileCInfo_t),intent(IN),target :: tileCInfo
end subroutine Orchestration_pushTileToGpuPipeline
! Local Variables:
! f90-program-indent: 3
! f90-do-indent: 3
! f90-type-indent: 3
! indent-tabs-mode: nil
! End:
