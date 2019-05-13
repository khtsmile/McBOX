program main
use constants
use variables
use CMFD, 		only : CMFD_lat, CMFD_type, n_acc
use simulation 
use omp_lib
use mpi

implicit none

integer :: i, ierg, iso 
integer :: provide 
real(8) :: time1, time2, time3, time4
real(8) :: k_sum 
logical :: isopened

!> Preparation for parallelization ===============================================
call MPI_Init_thread(MPI_THREAD_SINGLE, provide, ierr)
core = MPI_COMM_WORLD
call MPI_COMM_RANK(core,icore,ierr)
call MPI_COMM_SIZE(core,ncore,ierr)
score = 0 !server core=0

!> PreMC : Read input / Initialize / Set Random Seed etc. ========================
call premc

if(icore==score) then 
	print *, ''
	print *, ''
	print '(A30, I8)', '  > Num of Threads per Node   ', omp_get_max_threads()
	print '(A30, I8)', '  > Num of MPI Nodes          ', ncore
	print '(A30, I8)', '  > Num of Histories per Cycle', n_history
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

!> Stead-state Simlulation Start =================================================
time3 = omp_get_wtime()
Do curr_cyc=1, n_totcyc
	curr_act = curr_cyc - n_inact
	!> history wise transport simulation
	time1 = omp_get_wtime()
	call simulate_history()
	time2 = omp_get_wtime()
	
	if(icore==score) print '(A3,I4,A2,I4,A15, 1F10.7, A21,1F7.3,A3)', &
			'   ',curr_cyc,' /',n_totcyc,' cycle k eff : ', keff, '      elased time : ', time2-time1, 'sec'
	
Enddo
if(icore==score) print *, '   All Simulation Ends...' 

time4 = omp_get_wtime()
if(icore==score) print '(A16, F10.4, A5)', '   elapsed time', time4 - time3, 'sec'


deallocate(source_bank)
inquire(unit=prt_flux, opened=isopened)
if ( isopened ) close(prt_flux)
inquire(unit=prt_powr, opened=isopened)
if ( isopened ) close(prt_powr)

close(prt_keff)

call MPI_FINALIZE(ierr)

end program 
