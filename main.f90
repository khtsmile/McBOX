program main
use constants
use variables
use ENTROPY,    only : mprupon, entrp0
use FMFD_HEADER, only: p_fmfd, k_fmfd, e_fmfd, fmfdon, cmfdon
use simulation 
use DEPLETION_MODULE
use omp_lib
use mpi

implicit none

integer :: i, ierg, iso 
integer :: provide 
real(8) :: time1, time2, time3, time4
real(8) :: k_sum 
logical :: isopened
real(8) :: tt1, tt2, tt3
integer :: jj, kk

!> Preparation for parallelization ===============================================
!call omp_set_num_threads(1)
call MPI_Init_thread(MPI_THREAD_SINGLE, provide, ierr)
core = MPI_COMM_WORLD
call MPI_COMM_RANK(core,icore,ierr)
call MPI_COMM_SIZE(core,ncore,ierr)

!> PreMC : Read input / Initialize / Set Random Seed etc. ========================
call premc
call TIME_MEASURE

!> Stead-state Simlulation Start =================================================
call START_MSG

BURNUP : do

    if ( do_burn ) then
        call BURNUP_MSG
        call INIT_BURNUP
    end if

    ! steady-state calculation
    time3 = omp_get_wtime()
    Do curr_cyc = 1, n_totcyc
        curr_act = curr_cyc - n_inact
        !> history wise transport simulation
        time1 = omp_get_wtime()
        call simulate_history(curr_cyc)
        time2 = omp_get_wtime()
        t_tot(curr_cyc) = time2-time1
        call RUN_MSG
        if ( curr_cyc == n_totcyc ) exit
    Enddo
    time4 = omp_get_wtime()
    call END_MSG

    !> Gather Burnup Tallies
    call MPI_reduce_burnup()

    !> Make & Solve depletion matrix
    call depletion

    !> Check burnup loop exit condition
    if ( do_burn ) then
        if ( istep_burnup > nstep_burnup ) exit BURNUP
    else
        exit BURNUP
    end if

    call MPI_BARRIER(core,ierr)

end do BURNUP

deallocate(source_bank)
inquire(unit=prt_flux, opened=isopened)
if ( isopened ) close(prt_flux)
inquire(unit=prt_powr, opened=isopened)
if ( isopened ) close(prt_powr)
close(prt_keff)

call MPI_FINALIZE(ierr)

contains

function AVG(val)
    real(8):: avg
    real(8), intent(in):: val(:)

    avg = sum(val)/size(val)

end function

function STD(val)
    real(8):: std
    real(8), intent(in):: val(:)
    integer:: length
    real(8):: avg

    length = size(val)
    avg = sum(val)/length
    std = sqrt(dot_product((val-avg),(val-avg))/(length*length))*1E5
    if ( isnan(std) ) std = 0

end function

subroutine TIME_MEASURE
    use SIMULATION_HEADER, only: t_MC, t_det, t_tot
    implicit none
    allocate(t_MC(n_totcyc))
    allocate(t_det(n_totcyc))
    allocate(t_tot(n_totcyc))
    t_MC  = 0
    t_det = 0

end subroutine

! =============================================================================
! START_MSG
! =============================================================================
subroutine START_MSG

if ( icore==score ) then  

    write(*,*)
    write(*,*)
    write(*,10), '  > Num of Threads per Node   ', omp_get_max_threads()
    write(*,10), '  > Num of MPI Nodes          ', ncore
    write(*,10), '  > Num of Histories per Cycle', ngen
    write(*,11), '  > Skip Cycles:',n_inact , &
                 '  /  Active Cycles:', n_totcyc-n_inact
    write(*,*)
    if (tally_switch > 0) then 
        write(*,*), ' > Tally is On :: See tally.inp'
    else 
        write(*,*), ' > Tally is OFF' 
    endif 

    if ( mprupon ) then
        write(*,*), ' > m-PRUP is on', rampup, ccrt, scrt
    else
        write(*,*), ' > m-PRUP is OFF'
    end if

    if ( fmfdon ) then
        if ( cmfdon ) then
        write(*,*), ' > FMFD with CMFD is On'
        else
        write(*,*), ' > FMFD is On'
        end if
    else
        write(*,*), ' > FMFD is OFF'
    end if

    write(*,*)
    write(*,*), '   Transport Simulation Starts...' 

endif

10 format(A30,I9)
11 format(A16,I5,A19,I5)
12 format(A16)
13 format(A25,I4)
14 format(A22,I3)
15 format(A27,I4)
16 format(A22,I3)

end subroutine

! =============================================================================
! BURNUP_MSG
! =============================================================================
subroutine BURNUP_MSG
    write(*,10), '   =========================================='
    write(*,11), '      Burnup step', istep_burnup
    write(*,12), burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
    write(*,10), '   =========================================='

    10 format(A45)
    11 format(A17,I4)
    12 format(F14.2,A16)

end subroutine


! =============================================================================
! RUN_MSG
! =============================================================================
subroutine RUN_MSG
use ENTROPY, only: up_sign
implicit none
    
if ( icore == score ) then
    if ( curr_cyc <= n_inact ) then
    if ( up_sign ) then
    write(*,10), curr_cyc, time2-time1, "sec", entrp0, " | ", "keff", keff, &
                 "//", ngen
    up_sign = .false.
    else
    write(*,11), curr_cyc, time2-time1, "sec", entrp0, " | ", "keff", keff
    end if
    else
        kprt(curr_act) = keff
        write(*,12), curr_cyc, time2-time1, "sec", entrp0, " | ", &
            "keff", keff, &
            "avg", sum(kprt(1:curr_act))/dble(curr_act), &
            "SD", std(kprt(1:curr_act))
    end if
end if

10 format(i8,f9.2,1x,a,f10.5,1x,a,1x,a,f9.5,2x,a,i10)
11 format(i8,f9.2,1x,a,f10.5,1x,a,1x,a,f9.5)
12 format(i8,f9.2,1x,a,f10.5,1x,a,1x,2(a,f9.5,3x),a,f9.3)

end subroutine


! =============================================================================
! END_MSG
! =============================================================================
subroutine END_MSG
    use tally
    implicit none

if ( icore == score ) then
    write(*,*)
    write(*,*), '   Simulation of Burnup Step Terminated...'
    write(*,10), "    - Elapsed time    : ", &
        time4 - time3, 'sec', (time4-time3)/60, 'min'
    write(*,11), "    - Step Final keff : ", &
        sum(kprt(1:n_act))/dble(n_act), "+/-", STD(kprt(1:n_act))
    write(*,*)

    10 format(A,F10.3,A4,F8.2,A4)
    11 format(A,F10.6,A4,F8.3)


    ! multiplication factor
    if ( fmfdon ) then
    do ii = 1, n_inact
        write(*,12), k_fmfd(ii)
    end do
    do ii = n_inact+1, n_totcyc
        write(*,13), k_fmfd(ii), AVG(k_fmfd(n_inact+1:ii)), STD(k_fmfd(n_inact+1:ii))
    end do
    write(*,*)
    end if

    ! power distribution
    if ( fmfdon ) then
    do ii = 1, n_totcyc
    write(*,16), entrp3(ii)
    end do
    write(*,*)
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        p_fmfd(0,ii,jj,kk) = AVG(p_fmfd(1:,ii,jj,kk))
        e_fmfd(ii,jj,kk) = STD(p_fmfd(1:,ii,jj,kk))
    end do
    end do
    end do
    write(*,14), p_fmfd(0,:,:,:)
    write(*,*)
    write(*,14), e_fmfd(:,:,:)
    write(*,*)
    end if

    ! computing time
    t_MC = t_tot - t_det
    do ii = 1, n_totcyc
        write(*,15), t_MC(ii), t_det(ii), t_tot(ii)
    end do

    12 format(F10.6)
    13 format(2F10.6,F10.2)
    14 format(102ES15.7)
    15 format(3ES15.7)
    16 format(F10.5)

end if



end subroutine

end program 
