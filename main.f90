program main
use constants
use variables
use simulation 
use omp_lib
use mpi

use ace_xs
!use bank_header,             only : source_bank, fission_bank

implicit none

integer :: i, ierg, iso 
real(8) :: time1, time2, micro(5)
real(8) :: k_sum 
integer:: ii

!==============================================================================
!Preparation for parallelization
!call omp_set_num_threads(1)
call MPI_INIT(ierr)
core = MPI_COMM_WORLD
call MPI_COMM_RANK(core,icore,ierr)
call MPI_COMM_SIZE(core,ncore,ierr)
score = 0 !server core=0

!open(99, file="testfile.out",action="write", status="replace")

!**Set random seed
!**Read input file
!**Allocate dynamic variable
!**Get ace-format library
!**Initializing variable
!**Multi-group cross section tally (for all cell types)
call premc
derg = 100
cnt = 0
!call system_clock(count_rate=cnt_rate)
!call system_clock(count_max=cnt_max)

!open(prt_keff, file="keff.out",action="write", status="old")

allocate(source_bank(n_history)) 
call bank_initialize(source_bank)

!do i = 1, size(source_bank) 
!    write(wt_coord,*) source_bank(i)%xyz(:) 
!enddo 

!ierg = 24913
! ierg = 1
! iso = 2
! print *, ace(iso)%library, ace(iso)%E(ierg), 'MeV'
! print *, ''
! print *, ace(iso)%sigt(ierg)
! print *, ace(iso)%sigd(ierg)
! print *, ace(iso)%sigel(ierg)
! if (allocated(ace(iso)%sigf)) then 
!     print *, ace(iso)%sigf(ierg)
! else 
!     print *, 'xs_fis_total is not set '
! endif
! print *, getnu (iso,ace(iso)%E(ierg))
! print *, 'sig inel ', ace(iso)%sigt(ierg)- ace(iso)%sigd(ierg)-ace(iso)%sigel(ierg)-ace(iso)%sigf(ierg)
!                 
! do i = 1, ace(iso)%NXS(4) 
!     if (ace(iso)%TY(i)==19) then 
!         print *, ace(iso)%MT(i),ace(iso)%TY(i), ace(iso)%sig_MT(i)%cx(ierg), '***'
!     else 
!         print *, ace(iso)%MT(i),ace(iso)%TY(i), ace(iso)%sig_MT(i)%cx(ierg)
!     endif 
! enddo 
! 
! stop 

open(wt_coord, file="coordinates",action="write", status="replace")

print *, ''
print *, ''
print '(A18, I8)',  '    Num of threads', omp_get_max_threads()
print '(A20, I10)', '    Num of particles', n_history
print *, '   Transport Simulation Starts...' 

Do curr_cyc=1, n_totcyc

    curr_act = curr_cyc - n_inact
    !> history wise transport simulation
    time1 = omp_get_wtime()
    call simulate_history()
    time2 = omp_get_wtime()
!    print '(A3,I3,A2,I3,A15, 1F10.7, A21,1F6.3,A3)', &
!            '   ',curr_cyc,' /',n_totcyc,' cycle k eff : ', keff, '      elased time : ', time2-time1, 'sec'

    kprt(curr_cyc) = keff
    if ( curr_cyc <= n_inact+1 ) then
        write(*,1), curr_cyc, time2-time1, "sec", keff
    else
        write(*,2), curr_cyc, time2-time1, "sec", keff, "|", &
            sum(kprt(n_inact+1:curr_cyc))/dble(curr_cyc-n_inact), &
            "+/-", std(kprt(n_inact+1:curr_cyc))
    end if

!    do ii = 1, 500
!        print*, derg*(ii-1), cnt(ii)
!    end do
!    stop

    write(prt_keff,*) keff

Enddo
print *, '   All Simulation Ends...' 

1 format(i5,f9.3,1x,a,1x,f10.6)
2 format(i5,f9.3,1x,a,1x,f10.6,2x,a,1x,f10.6,1x,a,f12.6)

!do i = 1, size(source_bank) 
!    write(wt_coord,*) source_bank(i)%xyz(1:2)
!enddo 


print *, ''
!rewind(prt_keff)
!k_sum = 0 
!do i = 1, n_totcyc
!    read(prt_keff,*) keff
!    if (i<=n_inact) cycle 
!    k_sum = k_sum + keff
!enddo 
!print '(A22, F10.6)', '   Final keff mean :: ', k_sum / dble(n_totcyc-n_inact)
write(*,3), "    Final keff : ", &
    sum(kprt(n_inact+1:n_totcyc))/dble(n_totcyc-n_inact), &
    "+/-", STD(kprt(n_inact+1:n_totcyc))
3 format(a,F10.6,1x,a,F10.3)

call MPI_FINALIZE(ierr)

contains

function std(val)
    real(8):: std
    real(8), intent(in):: val(:)
    integer:: length
    real(8):: avg

    length = size(val)
    avg = sum(val)/length
    std = sqrt(dot_product((val-avg),(val-avg))/(length*length))*1E5

end function



end program 
