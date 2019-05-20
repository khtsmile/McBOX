module simulation
    use omp_lib

    use constants
    use tracking,             only : transport 
    use variables
    use particle_header
    use randoms
    use simulation_header
    use bank_header,         only : source_bank, fission_bank, temp_bank, thread_bank, bank_idx
    use geometry,             only : find_cell, find_cell_xyz, cell_contains
    use geometry_header,     only : cells, surfaces, universes, sgrid, lattices
    use surface_header, only : surface
    use XS_header
    use material_header
    use ace_header, only : ace
    use tally,             only: TallyCoord, TallyFlux
    
    
    implicit none 
    
    contains 
    
    subroutine simulate_history()
        integer :: i, j, isize, jsize, nthread
        type(particle) :: p
        logical :: found 
        integer :: tid, my_id
        !real(8) :: time1, time2, time3
        
        if (allocated(fission_bank)) & 
            call move_alloc(fission_bank, source_bank)
        allocate(fission_bank(0))
        
        isize = size(source_bank)
        
        !time1 = omp_get_wtime()
        
        !found = .false.
        !call p%initialize()
        !p%n_coord = 1
        !p%coord(1)%xyz = xyz
        !!call find_cell(p, found, j)
        !print *, cell_contains(cells(18), p)
        !stop
        
        k_col = 0 
        !$omp parallel private(p) shared(source_bank, fission_bank, temp_bank)
          thread_bank(:)%wgt = 0; bank_idx = 0
          
          !$omp do reduction(+:k_col)
            do i=1, isize 
                call p%initialize()
                call p%set(source_bank(i))
                
                !print *, ''
                !print *, '=================================================='
                !print '(A25,I4, F10.4)', 'Initial neutron energy   ', i, p%E 
                !print '(A25,4F10.4)', 'Initial neutron position ', &
                !p%coord(1)%xyz(:), sqrt(p%coord(1)%xyz(1)**2 + &
                !p%coord(1)%xyz(2)**2 + p%coord(1)%xyz(3)**2 )
                !write(wt_coord, *) p % coord(1) % xyz

                do while ((p%alive == .true.).and.(p%wgt > wgt_min))
                    call transport(p)
                enddo 
                
                !if buffer is almost full -> add to the fission bank
                if (bank_idx > int(size(thread_bank)*0.01*(98-omp_get_thread_num()))) then 
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
          
          !$omp critical
            isize = size(fission_bank)
            if(allocated(temp_bank)) deallocate(temp_bank)
            allocate(temp_bank(1:isize+bank_idx)) 
            if (isize>0) temp_bank(1:isize) = fission_bank(:)
            deallocate(fission_bank)
            temp_bank(isize+1:isize+bank_idx) = thread_bank(1:bank_idx)
            call move_alloc(temp_bank, fission_bank)
          !$omp end critical
          
        !$omp end parallel
        
        !time2 = omp_get_wtime()
        
        !print *, 'elapsed time', time2 - time1
        
        isize = size(fission_bank) 
        do i = 1, isize
            fission_bank(i)%wgt = real(n_history,8)/real(isize,8)
        enddo
        
        !print *, 'bank size  ', isize
        k_col = k_col / real(n_history,8)
        keff  = k_col
        k_col = 0.0d0
        
        !do i = 1, size(TallyFlux)
        !    TallyFlux(i) = TallyFlux(i)/(real(n_history,8)*TallyCoord(i)%vol)
        !enddo
        
        !print *, TallyFlux(:)
        !TallyFlux(:)=0
        
        !open(wt_coord, file="coordinates",action="write", status="old")
        !do i = 1, size(fission_bank) 
        !    write(wt_coord,*) fission_bank(i)%xyz(:) 
        !enddo 
        !close(wt_coord)
        
    end subroutine 

    
    subroutine bank_initialize(this)
        class(bank)    :: this(:)
        type(surface), pointer :: surfptr
        integer        :: i, j, cell_idx, mat_idx, univ_idx, iso
        integer        :: n_hist, i_hist, i_cell
        real(8)        :: min(3), max(3), xyz(3), e1, e2
        logical        :: found
        
        
        print *, '   Initializing Source...'
        
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
                            !write(prt_keff,*) xyz(:)
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
        
        print *, '   Source Initializing Complete...'
        
    
    end subroutine
    
        
    
end module
