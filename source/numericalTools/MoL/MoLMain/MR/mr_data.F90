module mr_data

    implicit none

    character(len=:), allocatable, save :: mr_slowMethod
    character(len=:), allocatable, save :: mr_fastMethod

    integer, save :: mr_nstages_slow, mr_nstages_fast
    integer, save :: mr_nsubcycle

    integer, save :: mr_nscratch_slow, mr_nscratch_fast

    ! RHS indexinng
    integer, dimension(:), allocatable, save :: FE, FI, FF

    ! Initial state for fast evolution
    integer, save :: FAST_INITIAL

    ! Butcher tableau for slow evolution
    integer, save :: kmax
    real, dimension(:,:,:), allocatable, save :: gamK,   wK
    real, dimension(:,:),   allocatable, save :: gamBar, wBar
    real, dimension(:), allocatable, save :: cS

    ! Butcher tableau for fast evolution
    real, dimension(:,:), allocatable, save :: AF
    real, dimension(:),   allocatable, save :: bF, cF

end module mr_data
