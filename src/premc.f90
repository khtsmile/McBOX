subroutine premc
use constants
use variables, only: E_mode, do_burn
use ace_xs, only : setugrid
use bank_header, only: source_bank
use simulation, only: bank_initialize 
use depletion_module
use input_reader

implicit none

integer :: i 
!==============================================================================
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

call read_depletion
if (do_burn) then 
    if(icore==score)  print '(A28)', '    Reading Depletion Lib...' 
	call getdepletionlibrary
endif 

!==============================================================================
!Set material library for burnup and equilibrium Xe135 search 
call setmat

!==============================================================================
!Set lethargy grid for hash-based energy look-up
call setugrid



end subroutine
