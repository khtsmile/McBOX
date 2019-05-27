module input_reader
    use constants
    use geometry_header
    use surface_header
    use variables
    use XS_header 
    use material_header
    use tally,                only: TallyCoord, TallyFlux, TallyPower
    use CMFD,                 only: getXYZ, CMFD_lat, n_skip, n_acc, CMFD_type
    
    use ace_header
    use ace_module 
    
    implicit none
    
    contains  

    
    
    subroutine init_var
        allocate(universes(0:0))
        universes(0)%univ_type     = 0
        universes(0)%univ_id     = 0
        universes(0)%xyz(:)        = 0
        universes(0)%ncell         = 0
        
        keff = 1
    end subroutine
    
    subroutine read_geom
        
        implicit none
        
        integer :: i, j, k, ix, iy,iz, idx, n, level
        integer :: i_cell, i_univ, i_lat, i_surf 
        integer :: ntemp, itemp
        integer :: ierr
        real(8) :: dtemp, xyz(3)
        character(100) :: line
        character(50) :: option, temp, test, mat_id
        character(1)  :: pnum
        character(20) :: title
        character(30) :: filename
        
        logical :: found 
        
        !Read geom.inp
        open(rd_geom, file="./inputfile/geom.inp",action="read", status="old")
    
        ierr = 0;
        do while (ierr.eq.0)
            read (rd_geom, FMT='(A)', iostat=ierr) line 
            if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
            
            
            !> option identifier 
            j = 0
            do while (j.le.len(line))
                j = j+1 
                if (line(j:j).eq.' ') exit
                option = line(1:j)        
            enddo 
            !<================================================================>!
            !> READ AND ADD INFO TO THE LINKED LIST
            select case (option)
            case ("title") 
                title = line(j+1:)
                filename = trim(title)//'_keff.out'
                
                inquire(file=filename, exist=found)
                if (found) then
                  if (icore == score) open(prt_keff, file=filename, status="old")
                else
                  if (icore == score) open(prt_keff, file=filename, status="new")
                end if
                
            case ("surf")
                isize = 0
                if (allocated(surfaces)) isize = size(surfaces) 
                isize = isize+1
                allocate(surfaces_temp(isize))
                surfaces_temp(1:isize-1) = surfaces(:) 
                
                call read_surf(surfaces_temp(isize), line)
                
                if(allocated(surfaces)) deallocate(surfaces)
                call move_alloc(surfaces_temp, surfaces)            
                
            case ("cell")
                isize = 0
                if (allocated(cells)) isize = size(cells) 
                isize = isize+1
                allocate(cells_temp(isize))
                if (isize > 1) cells_temp(1:isize-1) = cells(:) 
                call read_cell (cells_temp(isize), line) 
                if(allocated(cells)) deallocate(cells)
                call move_alloc(cells_temp, cells)            
                            
            case ("pin")
                isize = 0
                if (allocated(universes)) isize = size(universes)-1
                isize = isize+1
                allocate(universes_temp(0:isize))
                if (isize > 1) universes_temp(0:isize-1) = universes(:) 
                call read_pin (universes_temp(isize), line) 
                universes_temp(isize)%xyz(:) = 0
                
                univptr => universes_temp(isize)
                
                allocate(univptr%r(univptr%ncell-1))
                allocate(univptr%cell(univptr%ncell))
                if(allocated(universes)) deallocate(universes)
                call move_alloc(universes_temp, universes)

                !> generage cell from pin
                isize = 0
                if (allocated(cells)) isize = size(cells) 
                isize = isize + univptr%ncell
                allocate(cells_temp(isize))
                if (isize > 1) cells_temp(1:isize-univptr%ncell) = cells(:) 
                if(allocated(cells)) deallocate(cells)
                call move_alloc(cells_temp, cells)
                
                do i = 1, univptr%ncell-1
                    j = size(cells)-univptr%ncell+i
                    univptr%cell(i) = j
                    read(rd_geom,*) mat_id,  univptr%r(i)
                    
                    if (E_mode == 0) cells(j)%mat_idx = find_mat_idx(XS_MG,mat_id)
                    if (E_mode == 1) cells(j)%mat_idx = find_CE_mat_idx (materials, mat_id)
                enddo
                j = size(cells)
                univptr%cell(i) = j
                read(rd_geom,*) mat_id
                
                if (E_mode == 0) cells(j)%mat_idx = find_mat_idx(XS_MG,mat_id)
                if (E_mode == 1) cells(j)%mat_idx = find_CE_mat_idx (materials, mat_id)
                call gen_cells_from_pin (univptr, cells(j-univptr%ncell+1:j)) 
                
                
                !> generage surface from pin
                if (univptr%ncell > 1) then 
                    isize = 0
                    if (allocated(surfaces)) isize = size(surfaces) 
                    isize = isize + univptr%ncell -1
                    allocate(surfaces_temp(isize))
                    if (isize > 1) surfaces_temp(1:isize-univptr%ncell+1) = surfaces(:) 
                    if(allocated(surfaces)) deallocate(surfaces)
                    call move_alloc(surfaces_temp, surfaces)
                    j = size(surfaces)
                    call gen_surfs_from_pin (univptr, surfaces(j-univptr%ncell+2:j)) 
                endif 
            
            case ("lat")
                isize = 0
                if (allocated(lattices)) isize = size(lattices) 
                isize = isize+1
                allocate(lattices_temp(isize))
                if (isize > 1) lattices_temp(1:isize-1) = lattices(:) 
                
                lat_ptr => lattices_temp(isize)
                call read_lat(lat_ptr, line, rd_geom) 
                allocate(lat_ptr%lat(lat_ptr%n_xyz(1),lat_ptr%n_xyz(2),lat_ptr%n_xyz(3))) 
                
                iy = 0; 
                
                do while (iy.lt.(lat_ptr%n_xyz(2)*lat_ptr%n_xyz(3)))
                    read (rd_geom, FMT='(A)', iostat=ierr) line 
                    if ((len_trim(line)==0).or.(scan(line,"%"))/=0) then 
                        cycle  
                    else 
                        iy = iy+1 
                        temp = adjustl(line); temp = trim(temp); 
                        ix = 0; idx = 0 
                        do while (ix.lt.lat_ptr%n_xyz(1)) 
                            idx = idx+1
                            if (temp(idx:idx).eq.' ') then 
                                cycle 
                            else
                                call process_line(temp, idx, option) 
                            endif 
                            ix = ix+1 
                            read (option, *) lat_ptr%lat(ix,mod(iy,lat_ptr%n_xyz(2)),int(iy/lat_ptr%n_xyz(2))+1)
                        enddo 
                    endif 
                enddo 
                if(allocated(lattices)) deallocate(lattices)
                call move_alloc(lattices_temp, lattices)    
                
                
                
            case ('bc') 
                call read_bc (surfaces, line)
                
            case ('sgrid') 
                allocate(sgrid(6))
                call read_sgrid(line)
                                
            case default 
                print *, 'NO SUCH OPTION ::', option 
                stop
            end select
            
            
        enddo
        
        
        ! ===================================================================================== !
        !> add pure universes from cells(:)
        do i = 1, size(cells)
            !> if univ_id = 0 then add to base universe / else add to the tail
            if (cells(i)%univ_id == 0 ) then 
                universes(0)%ncell = universes(0)%ncell+1
            else
                found = .false.
                do j = 1, size(universes(1:))
                    if (universes(j)%univ_id == cells(i)%univ_id) then 
                        found = .true.
                        idx = j
                        exit 
                    endif
                enddo 
                if (found == .false.) then 
                    isize = size(universes(1:))
                    allocate(universes_temp(0:isize+1))
                    do j = 1, isize 
                        universes_temp(j) = universes(j)
                    enddo
                    obj_univ%univ_type     = 0
                    obj_univ%univ_id     = cells(i)%univ_id
                    obj_univ%xyz(:)        = 0 
                    obj_univ%ncell         = 1
                    universes_temp(isize+1) = obj_univ
                    call move_alloc(universes_temp, universes)
                elseif (.not. allocated(universes(idx)%r)) then 
                    universes(idx)%ncell = universes(idx)%ncell+1
                endif
            endif
        enddo         
        
        
        !> 3. Update surface info to subpin cells 
        do i = 1, size(cells)
            read(cells(i)%cell_id,*) temp
            if (temp(1:1) == 'p') then                  !> sub_pin cell 
                read(temp(2:2),'(I)') itemp
                read(temp(4:4),'(I)') j
                
                do i_univ = 0, size(universes) 
                    if (itemp == universes(i_univ)%univ_id) exit
                enddo 
                
                !print *, i_univ, universes(i_univ)%univ_id , universes(1)%ncell
                if (universes(i_univ)%ncell > 1) then
                  if (j == 1)then 
                      allocate(cells(i)%neg_surf_idx(1))
                      allocate(cells(i)%pos_surf_idx(0))
                      cells(i)%neg_surf_idx(1) = find_surf_idx(surfaces,temp(1:2)//'s1')
                  elseif (j == universes(i_univ)%ncell) then
                      allocate(cells(i)%neg_surf_idx(0))
                      allocate(cells(i)%pos_surf_idx(1))
                      write(line,*) j-1
                      cells(i)%pos_surf_idx(1) = find_surf_idx(surfaces,temp(1:2)//'s'//adjustl(line))
                  else 
                      allocate(cells(i)%pos_surf_idx(1))
                      write(line,*) j-1
                      !print *, 'check1'
                      cells(i)%pos_surf_idx(1) = find_surf_idx(surfaces,temp(1:2)//'s'//adjustl(line))
                      allocate(cells(i)%neg_surf_idx(1))
                      write(line,*) j
                      cells(i)%neg_surf_idx(1) = find_surf_idx(surfaces,temp(1:2)//'s'//adjustl(line))
                      !print *, 'check2'
                  endif
                endif
                cells(i)%operand_flag = 1
                
            else                                         !> ordinary cell
                ix=0; iy=0;
                cells(i)%nsurf = size(cells(i)%list_of_surface_IDs)
                do j = 1, cells(i)%nsurf
                    if (cells(i)%list_of_surface_IDs(j) < 0) ix = ix+1
                    if (cells(i)%list_of_surface_IDs(j) > 0) iy = iy+1
                enddo 
                allocate(cells(i)%neg_surf_idx(ix))
                allocate(cells(i)%pos_surf_idx(iy))
                ix=1; iy=1;
                do j = 1, cells(i)%nsurf
                    write(line, *) abs(cells(i)%list_of_surface_IDs(j))
                    if (cells(i)%list_of_surface_IDs(j) < 0) then 
                        cells(i)%neg_surf_idx(ix) = find_surf_idx(surfaces,adjustl(line))
                        ix = ix+1
                    
                    elseif (cells(i)%list_of_surface_IDs(j) > 0) then
                        cells(i)%pos_surf_idx(iy) = find_surf_idx(surfaces,adjustl(line))
                        iy = iy+1
                    endif
                enddo 
                
                !> cell translation
                if (cells(i)%nsurf == 1) then 
                    write(line,*) abs(cells(i)%list_of_surface_IDs(1))
                    idx = find_surf_idx(surfaces,adjustl(line))
                    !if (surfaces(idx)%surf_type == sqcx) 
                    !if (surfaces(idx)%surf_type == sqcy)
                    if (surfaces(idx)%surf_type == sqcz) then 
                        allocate(cells(i)%translation(3))
                        cells(i)%translation(1) = surfaces(idx)%parmtrs(1)
                        cells(i)%translation(2) = surfaces(idx)%parmtrs(2)
                        cells(i)%translation(3) = 0
                    endif
                    !if (surfaces(idx)%surf_type == cylx)
                    !if (surfaces(idx)%surf_type == cyly)
                    if (surfaces(idx)%surf_type == cylz) then 
                        allocate(cells(i)%translation(3))
                        cells(i)%translation(1) = surfaces(idx)%parmtrs(1)
                        cells(i)%translation(2) = surfaces(idx)%parmtrs(2)
                        cells(i)%translation(3) = 0
                    endif
                    if (surfaces(idx)%surf_type == sph) then 
                        allocate(cells(i)%translation(3))
                        cells(i)%translation(1) = surfaces(idx)%parmtrs(1)
                        cells(i)%translation(2) = surfaces(idx)%parmtrs(2)
                        cells(i)%translation(3) = surfaces(idx)%parmtrs(3)
                    endif
                    
                endif
            endif 
        enddo 
        
        !> 4. Add cells to pin universe
        !> Add the cells to the universe cell list   
        do i = 0, size(universes(1:))
            idx = 1
            if (.not.allocated(universes(i)%cell)) allocate(universes(i)%cell(universes(i)%ncell))
            do k = 1, size(cells)
                if (cells(k)%univ_id == universes(i)%univ_id) then 
                    universes(i)%cell(idx) = k
                    idx = idx+1
                endif 
            enddo 
        enddo
        
        !> 5. Change lattice universe name to universe index
        do i = 1, size(lattices) 
            do ix = 1, lattices(i)%n_xyz(1)
                do iy = 1, lattices(i)%n_xyz(2)
                    do iz = 1, lattices(i)%n_xyz(3)
                        itemp = find_univ_idx(universes,lattices(i)%lat(ix,iy,iz) )
                        lattices(i)%lat(ix,iy,iz) = itemp 
                    enddo 
                enddo 
            enddo 
        enddo 
        
        
        do i = 1, size(cells) 
            associate(this => cells(i))
                if (this%fill < 0) then 
                    this%filltype = FILL_MATERIAL
                elseif (in_the_list_univ(universes, this%fill)) then 
                    this%filltype = FILL_UNIVERSE
                elseif (in_the_list_lat(lattices, this%fill)) then
                    this%filltype = FILL_LATTICE
                else 
                    print *, 'ERROR : WRONG SHIT FILLING THIS CELL', cells(i)
                    stop
                endif
            end associate
        enddo 
        
        
        
        !> CHECK THE INPUT READ RESULT
        !call check_input_result(universes,lattices, cells,surfaces)
                
        !> READ DONE
        close(rd_geom)
        if(icore==score) print '(A25)', '    GEOM  READ COMPLETE...' 
        
    end subroutine

    

    


    subroutine read_pin (Pinobj, line)
        type(universe) :: Pinobj
        character(*) :: line
        character(50):: temp, option
        integer :: i, j, k
        
        j = 0; k = 0 
        do while (j.lt.len(line)) 
            j = j+1; k = k+1 
            if (line(j:j).ne.' ') then 
                call process_line (line, j, temp ) 
                
                if (k.eq.1) read(temp, *) option 
                if (k.eq.2) read(temp, *) Pinobj%univ_id
                if (k.eq.3) read(temp, *) Pinobj%ncell
            endif 
        enddo
        
    end subroutine 
    
    
    subroutine read_cell (Cellobj, line) 
        class(Cell) :: Cellobj
        character(*) :: line
        character(50):: temp, option, flag, mat_id
        integer :: i, j, k
    
        j = 0; k = 0 
        do while (j.lt.len(line)) 
            j = j+1; k = k+1 
            if (line(j:j).ne.' ') then 
                call process_line (line, j, temp ) 
                if (k.eq.1) read(temp, *) option 
                if (k.eq.2) read(temp, *) Cellobj%cell_id
                if (k.eq.3) read(temp, *) Cellobj%univ_id
                if (k.eq.4) then 
                    if (temp(1:4).eq.'fill') then 
                        Cellobj%mat_idx = -1
                        read(temp(5:len(trim(temp))), *) Cellobj%fill
                    else 
                        read(temp, *) mat_id
                        if (mat_id == 'outside') then 
                            Cellobj%mat_idx = 0
                        else 
                            !Cellobj%mat_idx = find_mat_idx(XS_MG,mat_id)
                            if (E_mode == 0) Cellobj%mat_idx = find_mat_idx(XS_MG,mat_id)
                            if (E_mode == 1) Cellobj%mat_idx = find_CE_mat_idx (materials, mat_id)
                        endif
                        Cellobj%fill = -1
                    endif 
                    temp = line(j+1: len(line)); temp = trim(temp)
                    read(temp(1:1), *) flag
                    if (adjustl(flag) == '&') then 
                        Cellobj%operand_flag = 1
                    elseif (adjustl(flag) == '|') then 
                        Cellobj%operand_flag = -1
                    else 
                        Cellobj%operand_flag = 0
                    endif 
                    call surfcount_for_cell(trim(temp(2:)),j, Cellobj)
                    
                endif 
            endif 
        enddo
    
    end subroutine 

    subroutine read_sgrid (line)
        character(*) :: line
        character(30):: temp, option, surf_id
        integer :: i, j, k, idx
    
        j = 0; k = 0 
        do while (j.lt.len(line)) 
            j = j+1; k = k+1 
            if (line(j:j).ne.' ') then 
                call process_line (line, j, temp ) 
                
                if (k.eq.1) read(temp, *) option
                if (k.eq.2) read(temp, *) sgrid(1) 
                if (k.eq.3) read(temp, *) sgrid(2) 
                if (k.eq.4) read(temp, *) sgrid(3) 
                if (k.eq.5) read(temp, *) sgrid(4) 
                if (k.eq.6) read(temp, *) sgrid(5) 
                if (k.eq.7) read(temp, *) sgrid(6) 
                
            endif 
        enddo
        
    end subroutine 
    
    subroutine surfcount_for_cell (line, idx, Cellobj)
        class(Cell) :: Cellobj
        character(*) :: line
        character(50):: temp
        integer      :: i, idx,idx_temp, nsurf
        
        line = adjustl(line)     
    
        idx = index(line, ' ') 
        
        nsurf = 0 
        do i = 1, len(line) 
            if ((line(i:i).eq.' ').and.(line(i-1:i-1).ne.' ')) nsurf = nsurf+1         
        enddo 
        !print *, nsurf
        
        allocate(Cellobj%list_of_surface_IDs(nsurf))
        
        idx_temp = 1; idx = 1; temp = line
        do i = 1, nsurf
            temp = adjustl(temp(idx:len(temp))) 
            idx = index(temp, ' ')
            read(temp(1:idx-1),*) Cellobj%list_of_surface_IDs(i)
        enddo 
    
    end subroutine 
    
    
    subroutine read_lat (Latobj, line, rd_geom)
        type(lattice) :: Latobj
        character(*) :: line
        character(50):: temp, option
        integer, intent(in) :: rd_geom
        integer :: i, j, k, ix, iy, ierr, idx

            j = 0; k = 0 
            do while (j.lt.len(line)) 
                j = j+1; k = k+1 
                if (line(j:j).ne.' ') then 
                    call process_line (line, j, temp ) 
                    
                    if (k.eq.1) read(temp, *) option
                    if (k.eq.2) read(temp, *) Latobj%lat_id    
                    if (k.eq.3) read(temp, *) Latobj%lat_type    
                    if (k.eq.4) read(temp, *) Latobj%xyz(1) 
                    if (k.eq.5) read(temp, *) Latobj%xyz(2) 
                    if (k.eq.6) read(temp, *) Latobj%xyz(3) 
                    
                    if (k.eq.7) read(temp, *) Latobj%n_xyz(1)
                    if (k.eq.8) read(temp, *) Latobj%n_xyz(2)
                    if (k.eq.9) read(temp, *) Latobj%n_xyz(3)
                    
                    if (k.eq.10) read(temp, *) Latobj%pitch(1)
                    if (k.eq.11) read(temp, *) Latobj%pitch(2)
                    if (k.eq.12) read(temp, *) Latobj%pitch(3)
                endif 
            enddo


    end subroutine
    
    subroutine gen_cells_from_pin (Pinobj, cellobj) 
        type(universe) :: Pinobj
        type(cell):: cellobj(:)
        integer :: i 
        character(1) :: num, pin_id
                
        do i = 1, size(cellobj)
            write(num, '(1I1)') i
            write(pin_id, '(I1)') pinobj%univ_id
            
            cellobj(i)%cell_id     = 'p'//pin_id//'c'//num             
            cellobj(i)%idx         = 0                                
            cellobj(i)%univ_id     = pinobj%univ_id
            cellobj(i)%nsurf     = 2
            
            cellobj(i)%fill = -1
            !cellobj(i)%list_of_surface_indices(:) =             !> TO BE EDITED
            
        enddo 
        
    end subroutine    
    
    subroutine gen_surfs_from_pin(Pinobj, surfobj) 
        type(universe) :: Pinobj
        type(surface):: surfobj(:)
        integer :: i 
        character(1) :: num, pin_id
        
        
        do i = 1, size(surfobj) 
            write(num, '(1I1)') i
            write(pin_id, '(I1)') pinobj%univ_id
            
            surfobj(i)%surf_id         = 'p'//pin_id//'s'//num                     
            surfobj(i)%surf_type     = cylz
            surfobj(i)%bc             = 0
            surfobj(i)%parmtrs(1)     = 0
            surfobj(i)%parmtrs(2)     = 0
            surfobj(i)%parmtrs(3)     = pinobj%r(i)
        enddo 
    end subroutine

    subroutine read_MG_XS 
        implicit none 
        integer :: i, j, i_mat, n_mat,i_group, j_group, ierr
        character(50) :: mat_id, line, option 
        
        open(rd_xs, file="./inputfile/MG_XS.inp",action="read", status="old")
        ierr = 0; n_mat = 0 
        !allocate(XS_MG(1)); allocate(XS_MG_temp(1));
        do while (ierr.eq.0)
            read (rd_xs, FMT='(A)', iostat=ierr) line 
            if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
            
            j = 0; 
            do while (j.le.len(line))
                j = j+1 
                if (line(j:j).eq.' ') exit
                option = line(1:j)        
            enddo 
            
            
            select case (option)
            case('group') 
                read(line(j+1:),*) n_group
                
            case('mat')  
                        
                n_mat = n_mat+1
                read(line(j+1:),*) mat_id
                
                allocate(XS_MG_temp(n_mat))
                if (n_mat > 1) XS_MG_temp(1:n_mat-1) = XS_MG(:) 
                
                
                XS_MG_ptr => XS_MG_temp(n_mat)
                allocate(XS_MG_ptr%sig_tr(n_group)) 
                allocate(XS_MG_ptr%sig_abs(n_group)) 
                allocate(XS_MG_ptr%sig_cap(n_group)) 
                allocate(XS_MG_ptr%sig_fis(n_group)) 
                allocate(XS_MG_ptr%nu(n_group)) 
                allocate(XS_MG_ptr%chi(n_group)) 
                allocate(XS_MG_ptr%sig_scat(n_group,n_group)) 
                
                XS_MG_temp(n_mat)%mat_id = mat_id
                
                
                read(rd_xs, *) (XS_MG_temp(n_mat)%sig_tr(i_group), i_group = 1, n_group)
                read(rd_xs, *) (XS_MG_temp(n_mat)%sig_abs(i_group), i_group = 1, n_group)
                read(rd_xs, *) (XS_MG_temp(n_mat)%sig_cap(i_group), i_group = 1, n_group)
                read(rd_xs, *) (XS_MG_temp(n_mat)%sig_fis(i_group), i_group = 1, n_group)
                read(rd_xs, *) (XS_MG_temp(n_mat)%nu(i_group), i_group = 1, n_group)
                read(rd_xs, *) (XS_MG_temp(n_mat)%chi(i_group), i_group = 1, n_group)
                do j_group = 1, n_group 
                    read(rd_xs, *) (XS_MG_temp(n_mat)%sig_scat(j_group,i_group), i_group = 1, n_group)
                enddo 
                if(allocated(XS_MG)) deallocate(XS_MG)
                call move_alloc(XS_MG_temp, XS_MG)
                
                
            case default
                print *, 'WRONG OPTION :: SHOULD BE mat OPTION'
                STOP 
                
            end select
            
        enddo
        
        !do i_mat = 1, n_mat 
        !    print *, i_mat, XS_MG(i_mat)%mat_id
        !    print '(7F10.7)', XS_MG(i_mat)%sig_scat(1,:)
        !enddo 
        
        close(rd_xs)

        !> READ DONE
        if(icore==score) print '(A25)', '    XS    READ COMPLETE...' 
    end subroutine    
    
    
    subroutine read_ctrl        
        integer :: i, j, k, idx, n, level
        integer :: ierr
        character(100) :: line
        character(50)  :: option, temp, test
        
        
        !Read ctrl.inp
        open(rd_ctrl, file="./inputfile/ctrl.inp",action="read", status="old")
    
        ierr = 0; CMFD_lat = -1
        do while (ierr.eq.0)
            read (rd_ctrl, FMT='(A)', iostat=ierr) line 
            if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
            !> option identifier 
            j = 0
            do while (j.le.len(line))
                j = j+1 
                if (line(j:j).eq.' ') exit
                option = line(1:j)        
            enddo 
            select case (option)
            case ("pop")
                read(line(j+1:), *) n_history, n_totcyc, n_act
                n_inact = n_totcyc - n_act
                
            case ("energy") 
                read(line(j+1:), *) E_mode
                
            case ("nugrid")
                read(line(j+1:), *) nugrid
            
            case ("tally") 
                read(line(j+1:), *) tally_switch
                if (tally_switch > 0 .and. icore == score) then 
                    open(prt_spec, file="flux.out",action="write", status="replace")
                    open(prt_spec, file="power.out",action="write", status="replace")
                endif 
                
                
            case ("CMFD") 
                read(line(j+1:), *) CMFD_lat, n_skip, n_acc
                CMFD_type = 1
            case ("pCMFD") 
                read(line(j+1:), *) CMFD_lat, n_skip, n_acc
                CMFD_type = 2 
            case ("power") 
                read(line(j+1:), *) Nominal_Power
            end select
        enddo
        close(rd_ctrl)
        if(icore==score) print '(A25)', '    CTRL  READ COMPLETE...' 
         
    end subroutine
    
    subroutine read_CE_mat
        implicit none 
        integer :: i, j,i_iso, n_mat, ierr, idx
        character(50) :: iso_id, line, option, mat_name
        
        open(rd_mat, file="./inputfile/CE_mat.inp",action="read", status="old")
        ierr = 0; n_mat = 0 
        do while (ierr.eq.0)
            read (rd_mat, FMT='(A)', iostat=ierr) line 
            if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
            
            j = 0; 
            do while (j.le.len(line))
                j = j+1 
                if (line(j:j).eq.' ') exit
                option = line(1:j)        
            enddo 
            
            
            select case (option)
            case ("mat")
                n_mat = n_mat+1
                allocate(materials_temp(n_mat))
                if (n_mat > 1) materials_temp(1:n_mat-1) = materials(:) 
                CE_mat_ptr => materials_temp(n_mat)
                
                line = line(j+1:);    j = index(line(:),' ')
                read(line(1:j-1),*) CE_mat_ptr%mat_name
                
                line = adjustl(line(j:)); j = index(line(:),' ')
                read(line(1:j-1),*) CE_mat_ptr%density_gpcc
                
                line = adjustl(line(j:)); j = index(line(:),' ')
                read(line(1:j-1),*) CE_mat_ptr%n_iso
                
                
                allocate(CE_mat_ptr%ace_idx(CE_mat_ptr%n_iso))
                allocate(CE_mat_ptr%numden(CE_mat_ptr%n_iso)) 
                
                do i_iso = 1, CE_mat_ptr%n_iso 
                    read(rd_mat, *) line, CE_mat_ptr%numden(i_iso)
                    line = adjustl(line); j = index(line(:),' ')
                    read(line(1:j-1),*) iso_id
                    CE_mat_ptr%ace_idx(i_iso) = find_ACE_iso_idx (ace, iso_id)
                enddo 
                
                if(allocated(materials)) deallocate(materials)
                call move_alloc(materials_temp, materials)            
                
            end select
        enddo
        close(rd_mat)
        if(icore==score) print '(A27)', '    CE MAT READ COMPLETE...' 
        
    end subroutine
    
    subroutine read_inventory 
        implicit none 
        integer :: i, j, iso, ierr
        character(50) :: mat_id, line, option 
        
        ! ======================================== !
        !         ACE read start (all nodes)
        ! ======================================== !
        open(rd_inven, file="./inputfile/inventory.inp",action="read", status="old")
        read(rd_inven,'(A)') library_path(:)
        read(rd_inven,*) num_iso
        allocate(ace(1:num_iso))

        do iso =1, num_iso
            read(rd_inven,*) i, ace(iso)%library
        end do
        
        call set_ace
        
        close(rd_inven)
        if(icore==score) print '(A29)', '    ACE Lib. READ COMPLETE...' 
        
    end subroutine
    
    subroutine read_tally
        integer :: i, j, k, idx, n, level, i_univ, i_lat, i_pin, n_pin
        integer :: a,b,c, i_xyz(3), i_save
        integer :: ierr
        character(100) :: line
        character(50)  :: option, temp, test
        
        
        !Read tally.inp
        open(rd_tally, file="./inputfile/tally.inp",action="read", status="old")
    
        ierr = 0;
        do while (ierr.eq.0)
            read (rd_tally, FMT='(A)', iostat=ierr) line 
            if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
            !> option identifier             
            j = index(adjustl(line),' ') -1 
            option = line(1:j)
            
            if (trim(option) /= "tally") then 
                print *, 'ERROR Tally Read :: tally.inp     ', option
                stop
            end if

            line = adjustl(line(j+1:))
            j = index(line,' ')-1
            option = trim(line(1:j))

            select case (option)
            case ("cell")
                read(line(j+1:), *) n
                allocate(TallyCoord(1:n))
                allocate(TallyFlux(1:n))
                allocate(TallyPower(1:n))
                TallyCoord(:)%flag = 1
                do i = 1, n
10                    read (rd_tally, FMT='(A)', iostat=ierr) line 
                    if ((len_trim(line)==0).or.(scan(line,"%"))/=0) goto 10  
                    level = 0; idx = index(line,'>')
                    do while (idx /= 0)
                        level = level + 1 
                        j = index(line,' '); read(line(1:j),*) temp ; line = line(j:) 
                        line = adjustl(line);j = index(line,' '); read(line(1:j),*) i_univ;  line = line(j:)
                        line = adjustl(line);j = index(line,' '); read(line(1:j),*) i_lat  ; line = line(j:)
                        line = adjustl(line);j = index(line,' '); read(line(1:j),*) TallyCoord(i)%coord(level)%lattice_x; line = line(j:)
                        line = adjustl(line);j = index(line,' '); read(line(1:j),*) TallyCoord(i)%coord(level)%lattice_y; line = line(j:)
                        line = adjustl(line);j = index(line,' '); read(line(1:j),*) TallyCoord(i)%coord(level)%lattice_z; line = line(j:)
                        
                        TallyCoord(i)%coord(level)%cell     = find_cell_idx(cells, adjustl(temp))
                        if (i_univ == 0) then 
                            TallyCoord(i)%coord(level)%universe = 0 
                        else 
                            TallyCoord(i)%coord(level)%universe = find_univ_idx(universes, i_univ)
                        endif
                        
                        if (i_lat == 0) then 
                            TallyCoord(i)%coord(level)%lattice = 0 
                        else 
                            TallyCoord(i)%coord(level)%lattice  =  find_lat_idx(lattices, i_lat)
                        endif
                        
                        idx = index(line,'>')
                        line = line(idx+1:); line = adjustl(line)
                    enddo 
                    idx = index(line, 'vol')
                    line = adjustl(line(4:)); read(line(1:),*) TallyCoord(i)%vol
                    TallyCoord(i)%n_coord = level
                enddo 
                TallyFlux(:)  = 0
                TallyPower(:) = 0 

            case("pin") 
                read(line(j+1:), *) n
                allocate(TallyCoord(1:n)); TallyCoord(n)%n_coord = 0 
                allocate(TallyFlux(1:n))
                allocate(TallyPower(1:n))
                TallyCoord(:)%flag = 0
                i = 1; 
                do while ( i <= n) 
                11  read (rd_tally, FMT='(A)', iostat=ierr) line
                    if ((len_trim(line)==0).or.(scan(line,"%"))/=0) goto 11

                    level = 0; idx = index(line,'>')
                    do while (idx /= 0)
                        level = level + 1 
                        !print *, 'level', level, line 
                        j = index(line,' '); read(line(1:j),*) temp ; line = line(j:)                         
                        line = adjustl(line);j = index(line,' '); read(line(1:j),*) i_univ;  line = line(j:)
                        line = adjustl(line);j = index(line,' '); read(line(1:j),*) i_lat  ; line = line(j:)
                        
                        if (i_lat /= 0) i_lat = find_lat_idx(lattices, i_lat)
                        line = adjustl(line);
                        if (line(1:3) == 'all') then 
                            a = lattices(i_lat)%n_xyz(1)
                            b = lattices(i_lat)%n_xyz(2)
                            c = lattices(i_lat)%n_xyz(3)
                            n_pin = a*b*c
                            i_save = i
                            do i_pin = 1, n_pin
                                i_xyz = getXYZ(i_pin, a,b,c)
                                
                                TallyCoord(i)%coord(1:level) = TallyCoord(i_save)%coord(1:level)
                                
                                TallyCoord(i)%coord(level)%lattice_x = i_xyz(1)
                                TallyCoord(i)%coord(level)%lattice_y = i_xyz(2)
                                TallyCoord(i)%coord(level)%lattice_z = i_xyz(3)
                                
                                if (i_lat == 0) then 
                                    TallyCoord(i)%coord(level)%lattice = 0 
                                else 
                                    TallyCoord(i)%coord(level)%lattice  =  i_lat
                                endif
                                TallyCoord(i)%n_coord = level
                                !> Reset the cell & univ designation (for pin only)
                                TallyCoord(i)%coord(level)%cell = 0
                                TallyCoord(i)%coord(level)%universe = 0 
                                
                                i = i + 1
                                
                            enddo 
                            idx = 0 
                        else
                            line = adjustl(line);j = index(line,' '); read(line(1:j),*) TallyCoord(i)%coord(level)%lattice_x; line = line(j:)
                            line = adjustl(line);j = index(line,' '); read(line(1:j),*) TallyCoord(i)%coord(level)%lattice_y; line = line(j:)
                            line = adjustl(line);j = index(line,' '); read(line(1:j),*) TallyCoord(i)%coord(level)%lattice_z; line = line(j:)
                            
                            TallyCoord(i)%coord(level)%cell     = find_cell_idx(cells, adjustl(temp))
                            
                            if (i_univ == 0) then 
                                TallyCoord(i)%coord(level)%universe = 0 
                            else 
                                TallyCoord(i)%coord(level)%universe = find_univ_idx(universes, i_univ)
                            endif
                            
                            if (i_lat == 0) then 
                                TallyCoord(i)%coord(level)%lattice = 0 
                            else 
                                TallyCoord(i)%coord(level)%lattice  =  i_lat
                            endif
                            
                            idx = index(line,'>')
                            line = line(idx+1:); line = adjustl(line)
                            
                            if (idx == 0) then 
                                i = i + 1
                                TallyCoord(i)%n_coord = level
                                !> Reset the cell designation (for pin only)
                                TallyCoord(i)%coord(level)%cell = 0
                            endif
                        endif 
                    enddo 
                    !> pin volume is not needed (all same)
                    TallyCoord(:)%vol = 1
                enddo                     
                TallyFlux(:)  = 0
                TallyPower(:) = 0
                
            end select
        enddo
        
        close(rd_tally)
        if(icore==score) print '(A25)', '    TALLY READ COMPLETE...' 
        
    end subroutine
    subroutine check_input_result(universes,lattices, cells,surfaces)
        type(surface) :: surfaces(:)
        type(lattice) :: lattices(:)
        type(universe):: universes(0:)
        type(cell)      :: cells(:)
        integer :: i, j, iy, iz
        
        print *, ' ========== SURFACE CHECK =========='
        print *, ' INDEX     ID      surf_type'
        do i = 1, size(surfaces) 
            print '(I5, A10, I10)', i, trim(surfaces(i)%surf_id), surfaces(i)%surf_type
        enddo 
        print *, ''
        print *, ''
        
        
        print *, ' ========== CELL CHECK =========='
        do i = 1, size(cells) 
            print *, 'cell number :', i 
            print *, 'cell id:', cells(i)%cell_id!, cells(i)%univ_id, cells(i)%mat_id
            print *, 'surf list: ', cells(i)%list_of_surface_IDs(:)
            print *, 'neg : ', surfaces(cells(i)%neg_surf_idx(:))%surf_id
            print *, 'pos : ', surfaces(cells(i)%pos_surf_idx(:))%surf_id
            print *, 'operand : ',cells(i)%operand_flag
            print *, 'fill type:', cells(i)%fill_type()
            print *, 'material idx', cells(i)%mat_idx
            print *, ''
            print *, ''
        enddo 
        
        
        print *, ' ========== UNIVERSE CHECK =========='
        do i = 0,size(universes(1:))
            print *, ''
            print *, 'univ_id = ', universes(i)%univ_id, '   # of cell', universes(i)%ncell 
            do j = 1, universes(i)%ncell
                print '(a4,I2,2A7)', 'cell',j,'    ', cells(universes(i)%cell(j))%cell_id
            enddo 
        enddo 
        print *, ''
        print *, ''
        
        
        print *, ' ========== LATTICE CHECK =========='
        do i = 1, size(lattices) 
            print *, ''
            print *, 'lattice id = ', lattices(i)%lat_id
            do iz = 1, lattices(i)%n_xyz(3)
                do iy = 1, lattices(i)%n_xyz(2)
                    print '(17I2)', (universes(lattices(i)%lat(j,iy,iz))%univ_id, j = 1, lattices(i)%n_xyz(1))
                enddo 
                
                print *, '' 
            enddo 
        enddo
        
        
        
            
    end subroutine 
    
end module
