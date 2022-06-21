module ml_imex

    implicit none

contains

subroutine ml_imex_ssp_222(tableau)
    use ml_rkInterface, only: imex_tableau_t

    implicit none

    type(imex_tableau_t), intent(out) :: tableau

    real, parameter :: g = 1d0 - 1d0/sqrt(2d0)

    tableau%p  = 2
    tableau%sE = 2
    tableau%sI = 2

    allocate(tableau%AE(2,2))
    allocate(tableau%bE(2))
    allocate(tableau%cE(2))

    allocate(tableau%AI(2,2))
    allocate(tableau%bI(2))
    allocate(tableau%cI(2))

    tableau%AE      = 0d0
    tableau%AE(2,1) = 1d0

    tableau%bE(1) = 1d0/2d0
    tableau%bE(2) = 1d0/2d0

    tableau%cE(1) = 0d0
    tableau%cE(2) = 1d0

    tableau%AI      = 0d0
    tableau%AI(1,1) = g
    tableau%AI(2,1) = 1d0 - 2d0*g
    tableau%AI(2,2) = g

    tableau%bI(1) = 1d0/2d0
    tableau%bI(2) = 1d0/2d0

    tableau%cI(1) = g
    tableau%cI(2) = 1d0 - g
end subroutine ml_imex_ssp_222

subroutine ml_imex_ssp_322(tableau)
    use ml_rkInterface, only: imex_tableau_t

    implicit none

    type(imex_tableau_t), intent(out) :: tableau

    tableau%p  = 2
    tableau%sE = 3
    tableau%sI = 3

    allocate(tableau%AE(3,3))
    allocate(tableau%bE(3))
    allocate(tableau%cE(3))

    allocate(tableau%AI(3,3))
    allocate(tableau%bI(3))
    allocate(tableau%cI(3))

    tableau%AE      = 0d0
    tableau%AE(3,2) = 1d0

    tableau%bE(1) = 0d0
    tableau%bE(2) = 1d0/2d0
    tableau%bE(3) = 1d0/2d0

    tableau%cE(1) = 0d0
    tableau%cE(2) = 1d0
    tableau%cE(3) = 1d0

    tableau%AI      = 0d0
    tableau%AI(1,1) = 1d0/2d0
    tableau%AI(2,1) = -1d0/2d0
    tableau%AI(2,2) = 1d0/2d0
    tableau%AI(3,2) = 1d0/2d0
    tableau%AI(3,3) = 1d0/2d0

    tableau%bI(1) = 0d0
    tableau%bI(2) = 1d0/2d0
    tableau%bI(3) = 1d0/2d0

    tableau%cI(1) = 1d0
    tableau%cI(2) = 0d0
    tableau%cI(3) = 1d0
end subroutine ml_imex_ssp_322

end module ml_imex
