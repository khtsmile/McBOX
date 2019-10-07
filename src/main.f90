program main
use constants
use variables
use CMFD,         only : CMFD_lat, CMFD_type, n_acc
use simulation 
use depletion_module

use omp_lib
use mpi


implicit none

integer :: i, ierg, iso 
integer :: provide 
real(8) :: time1, time2, time3, time4
real(8) :: k_sum 
logical :: isopened
real(8), allocatable :: kprt_vrc(:)

!> Preparation for parallelization ===============================================
call MPI_Init_thread(MPI_THREAD_SINGLE, provide, ierr)
core = MPI_COMM_WORLD
call MPI_COMM_RANK(core,icore,ierr)
call MPI_COMM_SIZE(core,ncore,ierr)
score = 0 !server core=0

!> PreMC : Read input / Initialize / Set Random Seed etc. ========================
call premc

!> Node Information and Calculation condition ====================================
if(icore==score) then  
    print *, ''
    print *, ''
    print '(A30, I9)', '  > Num of Threads per Node   ', omp_get_max_threads()
    print '(A30, I9)', '  > Num of MPI Nodes          ', ncore
    print '(A30, I9)', '  > Num of Histories per Cycle', n_history
    print '(A16,I5,A19,I5)', '  > Skip Cycles:',n_inact ,'  /  Active Cycles:', n_totcyc-n_inact
    print *, '' 
    if (CMFD_lat <= 0) then 
        print '(A16)', '  > CMFD is OFF ' 
    elseif (CMFD_lat > 0 .and. CMFD_type == 1 ) then 
        print '(A25,I4)', '  > CMFD is ON :: Lattice', CMFD_lat
        print '(A22,I3)', '  > CMFD Accumulations', n_acc
    elseif (CMFD_lat > 0 .and. CMFD_type == 2 ) then 
        print '(A27,I4)', '  > p-CMFD is ON :: Lattice', CMFD_lat
        print '(A22,I3)', '  > CMFD Accumulations', n_acc
    endif 
    print *, ''
    if (tally_switch > 0) then 
        print *, ' > Tally is On :: See tally.inp'
    else 
        print *, ' > Tally is OFF' 
    endif 
    print *, '' 
    print *, '   Transport Simulation Starts...' 
endif

allocate(kprt_vrc(size(kprt)))
!> Initiate Burnup Step Loop =========================================================
BURNUP : Do 
	!==============================================================================
	!Source bank initialize
	allocate(source_bank(n_history)) 
	call bank_initialize(source_bank)
	
	if (icore==score .and. do_burn == .true.) then 
		print '(a45)', '   =========================================='
		print '(a17,i4)', '      Burnup step', istep_burnup
		print '(f14.2,a16)', burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
		print '(a45)', '   =========================================='
	endif 
	
	k_sum = 0; 
	!> Initialize Burnup Tallies
	call Init_burnup() 
	
	!> Stead-state Simlulation Start =================================================
	time3 = omp_get_wtime()
	Do curr_cyc=1, n_totcyc
		curr_act = curr_cyc - n_inact
		!> history wise transport simulation
		time1 = omp_get_wtime()
		call simulate_history()
		time2 = omp_get_wtime()

		kprt(curr_cyc) = keff
		kprt_vrc(curr_cyc) = keff_vrc
		
		if ( icore == score ) then
		if ( curr_cyc <= n_inact+1 ) then
			write(*,1), curr_cyc, time2-time1, "sec", keff
		else
			write(*,2), curr_cyc, time2-time1, "sec", keff, "|", &
				sum(kprt(n_inact+1:curr_cyc))/dble(curr_cyc-n_inact), &
				"+/-", std(kprt(n_inact+1:curr_cyc))
		end if
		end if
	Enddo

	time4 = omp_get_wtime()

	if ( icore == score ) then 
		print *, '   Simulation of Burnup Step Terminated...' 

		write(*,3), '   elapsed time', time4 - time3, 'sec'
		print *
		write(*,4), "    Step Final keff : ", &
		sum(kprt(n_inact+1:n_totcyc))/dble(n_totcyc-n_inact), &
		"+/-", STD(kprt(n_inact+1:n_totcyc))
		
		write(prt_keff,'(a,4f12.6)'), "    Final keff : ", &
		sum(kprt(n_inact+1:n_totcyc))/dble(n_totcyc-n_inact), STD(kprt(n_inact+1:n_totcyc)), &
		sum(kprt_vrc(n_inact+1:n_totcyc))/dble(n_totcyc-n_inact), STD(kprt_vrc(n_inact+1:n_totcyc))
		
	endif

	!> Gather Burnup Tallies
	call MPI_reduce_burnup()

	!> Make & Solve depletion matrix
	call depletion()
	
	
	!Check burnup loop exit condition
	if(do_burn) then
		if(istep_burnup > nstep_burnup) exit BURNUP
	else
		exit BURNUP
	end if
	
	
	call MPI_BARRIER(core,ierr)
	
	deallocate(source_bank)
	
	
enddo BURNUP


1 format(i5,f9.3,1x,a,1x,f10.6)
2 format(i5,f9.3,1x,a,1x,f10.6,2x,a,1x,f10.6,1x,a,f12.6)
3 format(A16,F10.4,A5)
4 format(a,F10.6,1x,a,F10.3)




inquire(unit=prt_flux, opened=isopened)
if ( isopened ) close(prt_flux)
inquire(unit=prt_powr, opened=isopened)
if ( isopened ) close(prt_powr)

close(prt_keff)

call MPI_FINALIZE(ierr)

contains

function STD(val)
    real(8):: std
    real(8), intent(in):: val(:)
    integer:: length
    real(8):: avg

    length = size(val)
    avg = sum(val)/length
    std = sqrt(dot_product((val-avg),(val-avg))/(length*length))*1E5

end function



end program 
