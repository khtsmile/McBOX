subroutine premc
use constants
use input_reader
use ace_xs, only : setugrid

implicit none
!integer, allocatable :: seed(:)
integer :: n, i, i_iso, iso
integer :: seed_src


!==============================================================================
!List of output files
!if(icore==score) then
!  open(prt_keff,  file="out1_keff",action="write",status="replace") !position='append')
!  open(prt_ntpy,  file="out2_ntpy",action="write",status="replace") !position='append')
!  open(prt_pop,   file="out3_pop",action="write",status="replace") !position='append')
!  open(prt_tsize, file="out4_tsize",action="write",status="replace") !position='append')
!  open(monitor,   file="out5_monitor",action="write",status="replace") !position='append')
!  open(prt_pwr,   file="out6_pwr",action="write",status="replace") !position='append')
!  open(prt_pk,    file="out7_pk",action="write",status="replace") !position='append')
!  open(prt_pk2,   file="out8_pk2",action="write",status="replace") !position='append')
!  open(prt_pk3,   file="out9_pk3",action="write",status="replace") !position='append')
!  open(prt_time,  file="time_analysis",action="write",status="replace") !position='append')
!  open(prt_rSF,   file="conv_analysis",action="write",status="replace") !position='append')
!  open(prt_pCMFD, file="pCMFD_keff",action="write",status="replace") !position='append')
!  open(prt_ratio, file="ratio_estimator",action="write",status="replace") !position='append')
!
!  open(debug1,    file="debug1",action="write",status="replace") !position='append')
!  open(debug2,    file="debug2",action="write",status="replace") !position='append')
!
!endif
!
!
!!==============================================================================
!!Set random seed
!call random_seed(size=n)
!call system_clock(seed_src)
!allocate(seed(1:n))
!seed = abs(mod(seed_src*181*(icore-83)*359,104729))
!call random_seed(put=seed)
!if(icore==score) print *, "Setting random seed"

!==============================================================================
!Read input geometry and cross section data
call init_var
call read_ctrl
if (E_mode == 0) then 
    call read_MG_XS
elseif (E_mode == 1) then 
    print '(A29)', '    Continuous Energy Mode...' 
    call read_inventory
    call read_CE_mat
endif
call read_geom 
if(tally_switch > 0) call read_tally

!do i = 1, size(materials)
!    print *, materials(i)%mat_name
!    do i_iso = 1, materials(i)%n_iso
!        iso = materials(i)%ace_idx(i_iso)
!        print *, i_iso, ace(iso)%library 
!    enddo 
!    print *, '' 
!enddo 
!stop 
!==============================================================================
!Get events for TMC 
!call getevent


!==============================================================================
!Allocate allocatable variable and initilize them
!call setdat


!==============================================================================
!Allocating dynamic variables for TMC
!call settmc


!==============================================================================
!Set lethargy grid for hash-based energy look-up
call setugrid


!==============================================================================
!Set multi-group cross section tally
!call setmgtal


!==============================================================================
!Estimate heavy metal mass
!call heavymetal



end subroutine
