module imex_tableau

    implicit none

    ! ================= !
    !  IMEX-SSP(2,2,2)  !
    ! ================= !
    integer,                      save :: ssp_222_p
    integer,                      save :: ssp_222_s
    real, dimension(2,2), target, save :: ssp_222_AE, ssp_222_AI
    real, dimension(2),   target, save :: ssp_222_bE, ssp_222_bI
    real, dimension(2),   target, save :: ssp_222_cE, ssp_222_cI

contains

subroutine ssp_222_init()
    implicit none

    real, parameter :: g = 1d0 - 1d0/sqrt(2d0)

    ssp_222_p = 2
    ssp_222_s = 2

    ! Explicit
    ssp_222_AE      = 0d0
    ssp_222_AE(2,1) = 1d0

    ssp_222_bE(1) = 1d0/2d0
    ssp_222_bE(2) = 1d0/2d0

    ssp_222_cE(1) = 0d0
    ssp_222_cE(2) = 1d0

    ! Implicit
    ssp_222_AI      = 0d0
    ssp_222_AI(1,1) = g
    ssp_222_AI(2,1) = 1d0 - 2d0*g
    ssp_222_AI(2,2) = g

    ssp_222_bI(1) = 1d0/2d0
    ssp_222_bI(2) = 1d0/2d0

    ssp_222_cI(1) = g
    ssp_222_cI(2) = 1d0 - g
end subroutine ssp_222_init

subroutine imex_tableau_init()
    implicit none

    call ssp_222
end subroutine imex_tableau_init

subroutine imex_tableau_get(method, p, s, AE, bE, cE, AI, bI, cI)
    use ml_interface, only: ml_error

    implicit none

    character(len=*), intent(in)  :: method
    integer,          intent(out) :: p, s
    real, pointer                 :: AE(:,:), bE(:), cE(:)
    real, pointer                 :: AI(:,:), bI(:), cI(:)

    select case(trim(method))

    case ("imex-ssp-222")
        p =  ssp_222_p
        s =  ssp_222_s

        AE => ssp_222_AE
        bE => ssp_222_bE
        cE => ssp_222_cE

        AI => ssp_222_AI
        bI => ssp_222_bI
        cI => ssp_222_cI

    case default
        call ml_error("Unkown IMEX method")
    end select
end subroutine imex_tableau_get

end module imex_tableau
