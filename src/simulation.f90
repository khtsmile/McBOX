module simulation
    use omp_lib
    use mpi
    use constants
    use tracking,           only : transport, transport_DT, transport_VRC
    use variables
    use particle_header
    use randoms
    use bank_header,        only : bank, source_bank, fission_bank, temp_bank, thread_bank, bank_idx
    use geometry,           only : find_cell, find_cell_xyz, cell_contains
    use geometry_header,    only : cells, surfaces, universes, sgrid, lattices, lattice_coord
    use surface_header,     only : surface
    use XS_header
    use material_header
    use depletion_module,    only : inline_xenon

    use ace_header,         only : ace
    use tally,              only : TallyCoord, TallyFlux, TallyPower
    use CMFD,               only : CMFD_initialize_thread, CMFD_initialize, &
                                   idx_lat, CMFD_solve, CMFD_lat, n_skip, n_acc, &
                                   process_cmfd_par, normalize_CMFD_par 
                                    
    
    implicit none 
    
    contains 
    
    subroutine simulate_history()
    integer :: i, j, k, isize, i_surf
    integer :: a,b,c, n
    type(particle) :: p
    logical :: found 
    integer :: tid, my_id
    integer :: ista, iend
    real(8) :: Jtemp
    real(8), allocatable :: shape(:)
    integer :: i_xyz(3)
    real(8) :: rcv_buf
    integer :: realex, intex, restype, ndata, idata
    integer, dimension(0:4) :: blocklength, displacement, oldtype 
    integer, allocatable :: ircnt(:), idisp(:) 
    
    if (allocated(fission_bank)) & 
    call move_alloc(fission_bank, source_bank)
    allocate(fission_bank(0))
    isize = size(source_bank)
    !> Distribute source_bank to slave nodes 
    call para_range(1, isize, ncore, icore, ista, iend)        
        
        !> Initialize tally parameters
        k_col = 0; k_tl = 0; k_vrc = 0; fiss_vrc = 0; loss_vrc = 0;  
        cyc_power = 0; 

        !> Initialize CMFD tally 
        if (CMFD_lat > 0) call CMFD_initialize()
        !$omp parallel private(p) shared(source_bank, fission_bank, temp_bank)
          thread_bank(:)%wgt = 0; bank_idx = 0
          if (CMFD_lat > 0) call CMFD_initialize_thread()
          !$omp do reduction(+:k_col, k_tl, k_vrc)
            !do i=1, isize 
            do i= ista, iend 
                call p%initialize()
                call p%set(source_bank(i))
                
                do while ((p%alive == .true.).and.(p%wgt > wgt_min))
                    call transport(p)
                    !call transport_DT(p)
                    !call transport_VRC(p)
                enddo 
                
                !if buffer is almost full -> add to the fission bank
                if (bank_idx > int(size(thread_bank)*0.01*(90-OMP_GET_THREAD_NUM()))) then 
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
            
            !> normalize thread CMFD parameters (can be done outside critical)
            call normalize_CMFD_par()
            
          !$omp end critical
        !$omp end parallel
        
        !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !> Process tallied CMFD parameters ==========================================
        call process_CMFD_par()
        
        !> Gather tallies from the slave nodes =========================================
        call MPI_REDUCE(k_col, rcv_buf, 1, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
        k_col = rcv_buf
        
        call MPI_REDUCE(k_tl, rcv_buf, 1, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
        k_tl = rcv_buf
        
        call MPI_REDUCE(fiss_vrc, rcv_buf, 1, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
        fiss_vrc = rcv_buf 
        call MPI_REDUCE(loss_vrc, rcv_buf, 1, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
        loss_vrc = rcv_buf 
        !call MPI_REDUCE(k_vrc, rcv_buf, 1, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
        !k_vrc = rcv_buf !+ k_col
        
        call MPI_REDUCE(cyc_power, rcv_buf, 1, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
        cyc_power = rcv_buf
        if(curr_cyc > n_inact) avg_power = avg_power + cyc_power
        
        !> Calculate k_eff ==========================================================
        k_col = k_col / real(n_history,8)
        k_tl  = k_tl  / real(n_history,8) 
        k_vrc = fiss_vrc / real(n_history,8) 
        !k_vrc = fiss_vrc / loss_vrc
        keff_vrc = k_vrc
        !k_vrc = fiss_vrc / fiss_last
        !k_vrc = fiss_vrc / real(n_history,8)
        !fiss_last = fiss_vrc 
        !if (icore == score) print *, k_col, k_vrc, fiss_vrc / real(n_history,8) 
        !keff  = (k_tl + k_col) / 2.0d0 ; 
        keff = k_col
        !keff = k_tl
        !if (icore == score) write(prt_keff,'(5F12.7)') k_col, k_tl, k_vrc, fiss_vrc / real(n_history,8) 
        
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
        
        blocklength(0) = 1; blocklength(1) = 3; blocklength(2) = 3; blocklength(3) = 1; blocklength(4) = 1
        
        call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr) 
        displacement(0) = 0; displacement(1) = realex; displacement(2) = 4*realex
        displacement(3) = 7*realex; displacement(4) = 8*realex
        
        oldtype(0:3) = MPI_double_precision
        oldtype(4)   = MPI_INTEGER8
        allocate(temp_bank(1:sum(ircnt)))
                
        call MPI_TYPE_STRUCT (5, blocklength, displacement, oldtype,restype, ierr)
        call MPI_TYPE_COMMIT (restype, ierr)
        call MPI_ALLGATHERV(fission_bank, ndata, restype, temp_bank,ircnt,idisp,restype, MPI_COMM_WORLD, ierr)
        call MPI_TYPE_FREE(restype, ierr)
        
        deallocate(ircnt); deallocate(idisp); deallocate(fission_bank) 
        call move_alloc(temp_bank, fission_bank)
        
        !> Normalize source weight  =================================================
        isize = size(fission_bank)
        fission_bank(:)%wgt = real(n_history,8)/real(isize,8)
        
        !> Solve CMFD and apply FSD shape feedback ==================================
        if (CMFD_lat > 0 .and. curr_cyc > n_skip) then 
            a = lattices(idx_lat)%n_xyz(1)
            b = lattices(idx_lat)%n_xyz(2)
            c = lattices(idx_lat)%n_xyz(3)
            allocate(shape(a*b*c)) 
            
            if (icore == score) then 
                call CMFD_solve(lattices(idx_lat)%n_xyz(1),lattices(idx_lat)%n_xyz(2), &
                                lattices(idx_lat)%n_xyz(3), shape) 
            endif 
            call MPI_BCAST(shape, a*b*c, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr) 
            ! find fission_bank lattice index and apply shape
            do i = 1, isize 
                i_xyz = lattice_coord (lattices(idx_lat), fission_bank(i)%xyz)
                j = i_xyz(1) + a*(i_xyz(2)-1) + a*b*(i_xyz(3)-1) 
                fission_bank(i)%wgt = fission_bank(i)%wgt * shape(j)
            enddo
            deallocate(shape)
          
        endif 
        
        !> Normalize tally (flux & power) ===========================================
        if (tally_switch>0) then 
            !> Calculate tallied flux 
            isize = size(TallyFlux) 
            do i = 1, isize
                TallyFlux(i)  = TallyFlux(i)/(real(n_history,8)*TallyCoord(i)%vol)
                TallyPower(i) = TallyPower(i)/(real(n_history,8)*TallyCoord(i)%vol)
            enddo
            call MPI_REDUCE(TallyFlux, TallyFlux, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(TallyPower, TallyPower, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
            
            if (icore == score) then 
                TallyFlux(:) = TallyFlux(:) * dble(isize) / sum(TallyFlux)
                TallyPower(:) = TallyPower(:) * Nominal_Power / sum(TallyPower)
                !> Print TallyFlux if active cycle
                if ( curr_act > 0 ) then 
                    write(prt_flux, 100) TallyFlux(:)
                    write(prt_powr, 100) TallyPower(:)
                endif 
            100    format(<isize>ES15.7)
            endif 
            !TallyFlux(:) =0
            !TallyPower(:)=0
        endif
        
    
        !> Inline Equilibrium Xe-135 
        !call inline_xenon()
        
        
        !> initialize the global tally parameters ===================================
        !k_col = 0.0d0; k_tl = 0.0d0; 
        
        
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
            this(i) % wgt     = 1
            found             = .false.
            this(i) % uvw    = rand_vec()
            
            
            if (E_mode == 0) then !> multigroup MC
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
            
            elseif (E_mode == 1) then !> continuous energy MC
                        
                e1 = rang()*32.d0 + 1.d0
                e2 = rang()
                this(i)%E = fes(e1) + e2*(fes(e1+1)-fes(e1))
                                
                search_CE: do while (found == .false.) 
                    do j = 1, 3
                        this(i) % xyz(j) = rang()*(max(j)-min(j)) + min(j)
                    enddo
                    univ_idx = 0; xyz = this(i)%xyz
                    call find_cell_xyz(xyz, univ_idx, cell_idx)
                    if (cells(cell_idx)%mat_idx == 0) then 
                        found = .false. 
                        cycle search_CE
                    endif
                    mat_idx = cells(cell_idx)%mat_idx
                    
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
