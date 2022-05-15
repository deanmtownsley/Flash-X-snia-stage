subroutine get_scratch_indices(dirLims, i_new, j_new, k_new, i, j, k)
    implicit none
    integer ,intent(in) :: i, j ,k
    integer ,intent(in), dimension(LOW:HIGH,MDIM)  ::  dirLims
    integer ,intent(OUT) :: i_new, j_new, k_new
    !$omp declare target
    i_new = i - dirLims(LOW,1)
    j_new = j - dirLims(LOW,2)
    k_new = k - dirLims(LOW,3)
  
  end subroutine get_scratch_indices