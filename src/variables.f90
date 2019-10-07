module variables 

implicit none

! steady-state parameters 
    integer :: n_totcyc, n_inact, n_act
    integer :: curr_cyc, curr_act
    integer :: n_history
    
    integer :: E_mode !> 0 for MG // 1 for CE
    
    
    integer :: tally_switch !> o for off // 1 for on
    real(8) :: keff, k_col, k_tl
	real(8) :: k_vrc, fiss_vrc, loss_vrc, fiss_last, keff_vrc
    real(8), allocatable:: kprt(:)
    real(8) :: Nominal_Power
    real(8) :: cyc_power, avg_power
	
! depletion parameters 
	logical :: do_burn
	
	
	
	
	
    !==============================================================================
    integer ::    & 
    & ncore,      & !Total number of cores
    & icore,      & !My core id
    & score,      & !Server rank
    & ierr,       & !Error
    & core        !Global communicator

    ! Error message
    character(len=80) :: err_msg !> error message when exception handler is called

end module 
