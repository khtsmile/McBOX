subroutine premc
    use constants
    use input_reader
    use ace_xs, only : setugrid
    use bank_header, only: source_bank
    use simulation, only: bank_initialize 
    use ENTROPY
    
    implicit none
    
    !===========================================================================
    !Read input geometry and cross section data
    call init_var
    call read_ctrl
    if (E_mode == 0) then 
        if(icore==score) print '(A30)', '    Multi-group Energy Mode...' 
        call read_MG_XS
    elseif (E_mode == 1) then 
        if(icore==score) print '(A29)', '    Continuous Energy Mode...' 
        call read_inventory
        call read_CE_mat
    endif
    call read_geom 
    if(tally_switch > 0) call read_tally
    
    !===========================================================================
    !Set lethargy grid for hash-based energy look-up
    call setugrid

    ! ==========================================================================
    call ENTRP_INIT
    
    !===========================================================================
    !Source bank initialize
    allocate(source_bank(n_history)) 
    call bank_initialize(source_bank)

end subroutine
