module simulation
    use omp_lib
    use mpi
    use constants
    use tracking,           only : transport, transport_DT
    use variables
    use particle_header
    use randoms
    use simulation_header
    use bank_header,        only : source_bank, fission_bank, temp_bank, thread_bank, bank_idx
    use geometry,           only : find_cell, find_cell_xyz, cell_contains
    use geometry_header,    only : cells, surfaces, universes, sgrid, lattices, lattice_coord
    use surface_header,     only : surface
    use XS_header
    use material_header
    use ENTROPY
    use MPRUP,              only : GENSIZE
    use ace_header,         only : ace
    use tally,              only : TallyCoord, TallyFlux, TallyPower, tally1, &
                                   tally2
    use FMFD,               only : FMFD_initialize_thread, FMFD_initialize, &
                                   FMFD_solve, n_skip, n_acc, fsd, FMFD_ID, &
                                   process_FMFD, NORM_FMFD, fmfdon, &
                                   nfm, k_fmfd
                                    
    
    implicit none 
    
    contains 
    
subroutine simulate_history(cyc)
    use FMFD, only: fm_thread
    implicit none
    integer, intent(in):: cyc
    integer :: i, j, k, isize, i_surf
    integer :: a,b,c, n
    type(particle) :: p
    logical :: found 
    integer :: tid, my_id
    integer :: ista, iend
    real(8) :: Jtemp
    real(8), allocatable :: shape(:)
    integer :: i_xyz(3), id(3)
    real(8) :: rcv_buf
    integer :: realex, intex, restype, ndata, idata
    integer, dimension(0:4) :: blocklength, displacement, oldtype 
    integer, allocatable :: ircnt(:), idisp(:) 
    real(8) :: time1, time2

    if (allocated(fission_bank)) call move_alloc(fission_bank, source_bank)
    if ( icore == score ) then
        call SHENTROPY(source_bank)
        if ( mprupon .or. genup ) call GENSIZE(cyc)
    end if
    allocate(fission_bank(0))
    isize = size(source_bank)
    
    !> Distribute source_bank to slave nodes 
    call para_range(1, isize, ncore, icore, ista, iend)        
    k_col = 0; k_tl = 0;
    if ( fmfdon ) call FMFD_initialize()

    !$omp parallel private(p) shared(source_bank, fission_bank, temp_bank)
      thread_bank(:)%wgt = 0; bank_idx = 0
      if ( fmfdon ) call FMFD_initialize_thread()
      !$omp do reduction(+:k_col, k_tl)
        !do i=1, isize 
        do i= ista, iend 
            !print*, i
            call p%initialize()
            call p%set(source_bank(i))
            
            do while ((p%alive == .true.).and.(p%wgt > wgt_min))
                call transport(p,cyc)
                !call transport_DT(p)
            enddo 
            
            !if buffer is almost full -> add to the fission bank
            if (bank_idx > int(size(thread_bank)*0.01*(80-OMP_GET_THREAD_NUM()))) then 
              !$omp critical
                isize = size(fission_bank)
                if(allocated(temp_bank)) deallocate(temp_bank)
                allocate(temp_bank(1:isize+bank_idx)) 
                if (isize>0) temp_bank(1:isize) = fission_bank(:)
                deallocate(fission_bank)
                temp_bank(isize+1:isize+bank_idx) = thread_bank(1:bank_idx)
                call move_alloc(temp_bank, fission_bank)
              !$omp end critical
                bank_idx = 0
            endif

!            if ( mod(i,10000) == 0 ) then
!            !if ( mod(i,1) == 0 .and. i >= 8922 ) then
!            !if ( i > 920 ) then
!            write(*,1), fm_thread(34,1:17,1)%J1(2)
!            1 format(es15.7)
!            print*, "*cycle*", i
!            print*
!            pause
!            end if

        enddo
      !$omp end do

      !$omp critical
        !> gather thread fission bank
        isize = size(fission_bank)
        if(allocated(temp_bank)) deallocate(temp_bank)
        allocate(temp_bank(1:isize+bank_idx)) 
        if (isize>0) temp_bank(1:isize) = fission_bank(:)
        deallocate(fission_bank)
        temp_bank(isize+1:isize+bank_idx) = thread_bank(1:bank_idx)
        call move_alloc(temp_bank, fission_bank)
        
        !> normalize thread FMFD parameters (can be done outside critical)
        if ( fmfdon ) call NORM_FMFD()
        
      !$omp end critical
    !$omp end parallel

    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !> Process tallied FMFD parameters ==========================================
    if ( fmfdon ) call PROCESS_FMFD()
    
    !> Gather keff from the slave nodes =========================================
    call MPI_REDUCE(k_col,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    k_col = rcv_buf
    
    call MPI_REDUCE(k_tl,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    k_tl = rcv_buf
    
    
    !> Calculate k_eff ==========================================================
    k_col = k_col / real(ngen,8)
    k_tl  = k_tl  / real(ngen,8) 
    !keff  = (k_tl + k_col) / 2.0d0 ; 
    keff = k_col
    !print*, "MC", keff
    
    if (icore == score) write(prt_keff,*) keff, k_col, k_tl
    
    call MPI_BCAST(keff, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr) 
    
    !> Gather fission_bank from the slave nodes =================================        
    ndata = size(fission_bank)
    allocate(ircnt(1:ncore))
    allocate(idisp(1:ncore))
    
    do i = 1, ncore
        idata = ndata
        call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
        ircnt(i) = idata
    enddo 
    idisp(1) = 0
    do i = 2, ncore 
        idisp(i) = idisp(i-1) + ircnt(i-1)
    enddo
    
    blocklength(0) = 1
    blocklength(1) = 3
    blocklength(2) = 3
    blocklength(3) = 1
    blocklength(4) = 1
    
    call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr) 
    displacement(0) = 0; displacement(1) = realex; displacement(2) = 4*realex
    displacement(3) = 7*realex; displacement(4) = 8*realex
    
    oldtype(0:3) = MPI_double_precision
    oldtype(4)   = MPI_INTEGER8
    allocate(temp_bank(1:sum(ircnt)))
            
    call MPI_TYPE_STRUCT (5, blocklength, displacement, oldtype,restype, ierr)
    call MPI_TYPE_COMMIT (restype, ierr)
    call MPI_ALLGATHERV(fission_bank,ndata,restype,temp_bank,ircnt,idisp,restype,MPI_COMM_WORLD,ierr)
    call MPI_TYPE_FREE(restype, ierr)
    
    deallocate(ircnt); deallocate(idisp); deallocate(fission_bank) 
    call move_alloc(temp_bank, fission_bank)
    
    !> Normalize source weight  =================================================
    isize = size(fission_bank)
    fission_bank(:)%wgt = real(ngen,8)/real(isize,8)
    
    !> Normalize tally (flux & power) ===========================================
    !if ( tally_switch > 0 .and. cyc > n_inact ) then 
    if ( tally_switch > 0 ) then 
        !> Calculate tallied flux 
        isize = size(TallyFlux) 
        do i = 1, isize
            TallyFlux(i)  = TallyFlux(i)/(real(ngen,8)*TallyCoord(i)%vol)
            TallyPower(i) = TallyPower(i)/(real(ngen,8)*TallyCoord(i)%vol)
        enddo

        call MPI_REDUCE(TallyFlux,tally1,isize,MPI_DOUBLE_PRECISION,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(TallyPower,tally2,isize,MPI_DOUBLE_PRECISION,MPI_SUM,score,MPI_COMM_WORLD,ierr)

!        if (icore == score) then 
!            tally1(:) = tally1(:) * dble(isize) / sum(tally1)
!            tally2(:) = tally2(:) * Nominal_Power / sum(tally2)
            !> Print tally1 if active cycle
!            write(prt_flux, 100) tally1(:)
!            write(prt_powr, 100) tally2(:)
!            100 format(<isize>ES15.7)
!        endif 
    endif
    
    !> Solve FMFD and apply FSD shape feedback ==================================
    if ( fmfdon .and. cyc > n_skip) then 
        if (icore == score) then
            call CPU_TIME(time1)
            k_fmfd(cyc) = keff
            call FMFD_SOLVE(k_fmfd(cyc),fsd)
            call CPU_TIME(time2)
            t_det(cyc) = time2-time1
        end if
        call MPI_BCAST(fsd, nfm(1)*nfm(2)*nfm(3), &
            MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr) 
        ! find fission_bank lattice index and apply shape
        do i = 1, isize 
            id = FMFD_ID(fission_bank(i)%xyz)
            fission_bank(i)%wgt = fission_bank(i)%wgt * fsd(id(1),id(2),id(3))
        enddo
      
    endif 
    
    !> initialize the global tally parameters ===================================
    k_col = 0.0d0; k_tl = 0.0d0; 
    
    
end subroutine 


subroutine bank_initialize(this)
    class(bank)    :: this(:)
    type(surface), pointer :: surfptr
    integer        :: i, j, cell_idx, mat_idx, univ_idx, iso
    integer        :: n_hist, i_hist, i_cell
    real(8)        :: min(3), max(3), xyz(3), e1, e2
    logical        :: found
    
    
    if(icore==score) print *, '   Initializing Source...'
    
    if (.not. allocated(sgrid)) then 
        print *, 'not allocated'
        do i = 1, size(cells(i_cell)%pos_surf_idx)
            surfptr => surfaces(cells(i_cell)%pos_surf_idx(i))
            if (surfptr%surf_type == pz) then 
                max(3) = surfptr%parmtrs(1)
            elseif (surfptr%surf_type == sqcz) then 
                min(1) = -surfptr%parmtrs(3)
                min(2) = -surfptr%parmtrs(3)
                max(1) =  surfptr%parmtrs(3)
                max(2) =  surfptr%parmtrs(3)
            elseif (surfptr%surf_type == cylz) then
                min(1) = -surfptr%parmtrs(3)
                min(2) = -surfptr%parmtrs(3)
                max(1) =  surfptr%parmtrs(3)
                max(2) =  surfptr%parmtrs(3)
            endif 
        enddo 
    else 
        min(1:3) = sgrid(1:3)
        max(1:3) = sgrid(4:6)
    endif 
            
    
    do i = 1, size(this)
        this(i) % wgt = 1
        found         = .false.
        this(i) % uvw = rand_vec()
        
        ! multigroup MC
        if (E_mode == 0) then
            this(i) % G = 1
            search_MG: do while (found == .false.) 
                do j = 1, 3
                    this(i) % xyz(j) = rang()*(max(j)-min(j)) + min(j)
                enddo
                univ_idx = 0; xyz = this(i)%xyz
                call find_cell_xyz(xyz, univ_idx, cell_idx)
                if (cells(cell_idx)%mat_idx == 0) then 
                    found = .false. 
                    cycle search_MG
                endif
                mat_idx = cells(cell_idx)%mat_idx
                !print '(I3, I3, A20)', i, cell_idx, XS_MG(mat_idx)%mat_id
                do j = 1, size(XS_MG(mat_idx)%sig_fis)  
                    if(XS_MG(mat_idx)%sig_fis(j) > 0.0001) then 
                        found = .true.  
                        exit search_MG
                    else 
                        found = .false.
                        cycle search_MG
                    endif 
                enddo 
                
            enddo search_MG
        

        ! continuous energy MC
        elseif (E_mode == 1) then
                    
            e1 = rang()*32.d0 + 1.d0
            e2 = rang()
            this(i)%E = fes(e1) + e2*(fes(e1+1)-fes(e1))
                            
            search_CE: do while ( found == .false.) 
                do j = 1, 3
                    this(i) % xyz(j) = rang()*(max(j)-min(j)) + min(j)
                enddo
                univ_idx = 0; xyz = this(i)%xyz
                call find_cell_xyz(xyz, univ_idx, cell_idx)
                ! within a material region
                if (cells(cell_idx)%mat_idx == 0) then 
                    found = .false. 
                    cycle search_CE
                endif
                mat_idx = cells(cell_idx)%mat_idx
                ! within a fissionable material region
                do j = 1, materials(mat_idx)%n_iso
                    iso = materials(mat_idx)%ace_idx(j)
                    if(ace(iso)%jxs(2) /= 0) then 
                        found = .true.  
                    endif 
                enddo 
                
            enddo search_CE 
            
            
        endif
        
        
        !this(i) % xyz   = (/0.001, -10, 1/)
        !this(i) % uvw    = (/-1, 0, 0/)!rand_vec()
        !this(i) % G     = 1
        
    enddo 
    
    if(icore==score) print *, '   Source Initializing Complete...'
    

end subroutine

    
subroutine para_range(n1, n2, nprocs, irank, ista, iend)
    integer :: iwork1, iwork2 
    integer, intent(in) :: n1, n2, nprocs, irank 
    integer, intent(inout) :: ista, iend
    
    iwork1 = (n2 - n1 + 1) / nprocs
    iwork2 = MOD(n2 - n1 + 1, nprocs)
    ista = irank * iwork1 + n1 + MIN(irank, iwork2)
    iend = ista + iwork1 - 1
    if (iwork2 > irank) iend = iend + 1

end subroutine 

end module
