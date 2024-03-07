!! @copyright Copyright 2024 UChicago Argonne, LLC and contributors
subroutine Orchestration_setupPipelineForCpuTasks(MH_taskFunction, nThreads)
    use Orchestration_interfaceTypeDecl, ONLY : milhoja_runtime_taskFunction

    implicit none

    procedure(milhoja_runtime_taskFunction)            :: MH_taskFunction
    integer,                                intent(IN) :: nThreads

end subroutine Orchestration_setupPipelineForCpuTasks

