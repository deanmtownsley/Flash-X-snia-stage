module ml_erk

    implicit none

contains

subroutine ml_erk_euler(tableau)
    use ml_rkInterface, only: erk_tableau_t

    implicit none

    type(erk_tableau_t), intent(out) :: tableau

    tableau%p = 1
    tableau%s = 1

    allocate(tableau%A(1,1))
    allocate(tableau%b(1))
    allocate(tableau%c(1))

    tableau%A(1,1) = 0d0

    tableau%b(1) = 1d0

    tableau%c(1) = 0d0
end subroutine ml_erk_euler

subroutine ml_erk_rk2(tableau, alp)
    use ml_rkInterface, only: erk_tableau_t

    implicit none

    type(erk_tableau_t), intent(out) :: tableau
    real, intent(in) :: alp

    tableau%p = 2
    tableau%s = 2

    allocate(tableau%A(2,2))
    allocate(tableau%b(2))
    allocate(tableau%c(2))

    tableau%A      = 0d0
    tableau%A(2,1) = alp

    tableau%b(1) = 1d0 - 1d0/(2d0*alp)
    tableau%b(2) = 1d0/(2d0*alp)

    tableau%c(1) = 0d0
    tableau%c(2) = alp
end subroutine ml_erk_rk2

subroutine ml_erk_rk3_ssp(tableau)
    use ml_rkInterface, only: erk_tableau_t

    implicit none

    type(erk_tableau_t), intent(out) :: tableau

    tableau%p = 3
    tableau%s = 3

    allocate(tableau%A(3,3))
    allocate(tableau%b(3))
    allocate(tableau%c(3))

    tableau%A      = 0d0
    tableau%A(2,1) = 1d0
    tableau%A(3,1) = 1d0/4d0
    tableau%A(3,2) = 1d0/4d0

    tableau%b(1) = 1d0/6d0
    tableau%b(2) = 1d0/6d0
    tableau%b(3) = 2d0/3d0

    tableau%c(1) = 0d0
    tableau%c(2) = 1d0
    tableau%c(3) = 1d0/2d0
end subroutine ml_erk_rk3_ssp

subroutine ml_erk_rk4(tableau)
    use ml_rkInterface, only: erk_tableau_t

    implicit none

    type(erk_tableau_t), intent(out) :: tableau

    tableau%p = 4
    tableau%s = 4

    allocate(tableau%A(4,4))
    allocate(tableau%b(4))
    allocate(tableau%c(4))

    tableau%A      = 0d0
    tableau%A(2,1) = 1d0/2d0
    tableau%A(3,2) = 1d0/2d0
    tableau%A(4,3) = 1d0

    tableau%b(1) = 1d0/6d0
    tableau%b(2) = 1d0/3d0
    tableau%b(3) = 1d0/3d0
    tableau%b(4) = 1d0/6d0

    tableau%c(1) = 0d0
    tableau%c(2) = 1d0/2d0
    tableau%c(3) = 1d0/2d0
    tableau%c(4) = 1d0
end subroutine ml_erk_rk4

end module ml_erk
