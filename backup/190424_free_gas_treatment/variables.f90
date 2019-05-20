module variables 

implicit none

! steady-state parameters 
    integer :: n_totcyc, n_inact, n_act
    integer :: curr_cyc, curr_act
    integer :: n_history
    
    integer :: E_mode !> 0 for MG // 1 for CE
    
    
    integer :: tally_switch !> o for off // 1 for on
    real(8) :: keff, k_col, k_tl
    real(8), allocatable:: kprt(:)
    
        
    !==============================================================================
    integer ::    & 
    & ncore,      & !Total number of cores
    & icore,      & !My core id
    & score,      & !Server rank
    & ierr,       & !Error
    & core        !Global communicator

    !==============================================================================
    ! real(8) :: macro_t,  macro_a, macro_f, macro_nuf
    ! !$OMP THREADPRIVATE(macro_t)
    ! !$OMP THREADPRIVATE(macro_a)
    ! !$OMP THREADPRIVATE(macro_f)
    ! !$OMP THREADPRIVATE(macro_nuf)
    
    ! Error message
    character(len=80) :: err_msg !> error message when exception handler is called

end module 
