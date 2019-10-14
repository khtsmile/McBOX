module input_reader
    use constants
    use geometry_header
    use surface_header
    use variables
    use XS_header 
    use material_header
    use tally,    only: TallyCoord, TallyFlux, TallyPower, CoordStruct
    use ENTROPY,  only: mprupon, rampup, crt1, crt2, elength
    use FMFD,     only: n_skip, n_acc, fm0, fm1, fm2, nfm, dfm, fcr, fcz, &
                        fmfdon, cmfdon
    use CMFD,     only: ncm
    use geometry, only: getXYZ
    use read_functions
    use depletion_module
    use ace_header
    use ace_module
    
    implicit none
    
    contains 
    
! =============================================================================
! INIT_VAR
! =============================================================================
subroutine init_var
    allocate(universes(0:0))
    universes(0)%univ_type = 0
    universes(0)%univ_id   = 0
    universes(0)%xyz(:)    = 0
    universes(0)%ncell     = 0
    
    keff = 1D0
    avg_power = 0D0

end subroutine
    
! =============================================================================
! READ_GEOM reads the input about the geometry
! =============================================================================
subroutine read_geom
    
    implicit none
    
    integer :: i, j, k, ix, iy,iz, idx, n, level
    integer :: i_cell, i_univ, i_lat, i_surf 
    integer :: ntemp, itemp
    integer :: ierr
    real(8) :: dtemp, xyz(3)
    character(100) :: line
    character(100) :: option, temp, test, mat_id
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
            allocate(surfaces_temp(1:isize))
            if ( isize > 1 ) surfaces_temp(1:isize-1) = surfaces(:) 
            
            call read_surf(surfaces_temp(isize), line)
            
            if(allocated(surfaces)) deallocate(surfaces)
            call move_alloc(surfaces_temp, surfaces)
            
        case ("cell")
            isize = 0
            if (allocated(cells)) isize = size(cells) 
            isize = isize+1
            allocate(cells_temp(1:isize))
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
            
            allocate(univptr%r(1:univptr%ncell-1))
            allocate(univptr%cell(1:univptr%ncell))
            if(allocated(universes)) deallocate(universes)
            call move_alloc(universes_temp, universes)

            !> generage cell from pin
            if (allocated(cells)) then 
                isize = size(cells) 
                isize = isize + univptr%ncell
                allocate(cells_temp(1:isize))
                if (isize > 1) cells_temp(1:isize-univptr%ncell) = cells(:) 
                deallocate(cells)
            else 
                isize = 0
                isize = isize + univptr%ncell
                allocate(cells_temp(1:isize))
            endif
            
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
                allocate(surfaces_temp(1:isize))
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
            allocate(lattices_temp(1:isize))
            if (isize > 1) lattices_temp(1:isize-1) = lattices(:) 
            
            lat_ptr => lattices_temp(isize)
            call read_lat(lat_ptr, line, rd_geom) 
            allocate(lat_ptr%lat(1:lat_ptr%n_xyz(1),1:lat_ptr%n_xyz(2),1:lat_ptr%n_xyz(3))) 
            
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
                        if (mod(iy,lat_ptr%n_xyz(2)) == 0) then 
                            read (option, *) lat_ptr%lat(ix,lat_ptr%n_xyz(2), &
                                                CEILING(dble(iy)/dble(lat_ptr%n_xyz(2))))
                        else 
                            read (option, *) lat_ptr%lat(ix,mod(iy,lat_ptr%n_xyz(2)), &
                                                CEILING(dble(iy)/dble(lat_ptr%n_xyz(2))))
                        endif 
                    enddo 
                endif 
            enddo 
            if(allocated(lattices)) deallocate(lattices)
            call move_alloc(lattices_temp, lattices)    
            
            
            
        case ('bc') 
            call read_bc (surfaces, line)
            
        case ('sgrid') 
            allocate(sgrid(1:6))
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
            if ( .not. found ) then 
                isize = size(universes(1:))
                allocate(universes_temp(0:isize+1))
                do j = 1, isize 
                    universes_temp(j) = universes(j)
                enddo
                obj_univ%univ_type = 0
                obj_univ%univ_id   = cells(i)%univ_id
                obj_univ%xyz(:)    = 0 
                obj_univ%ncell     = 1
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

! =============================================================================
! READ_PIN processes the input for the pin geometry
! =============================================================================
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
    
! =============================================================================
! READ_CELL processes the input for the cell geometry
! =============================================================================
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

! =============================================================================
! READ_SGRID
! =============================================================================
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
    
! =============================================================================
! SURFCOUNT_FOR_CELL
! =============================================================================
subroutine surfcount_for_cell (line, idx, Cellobj)
    class(Cell) :: Cellobj
    character(*) :: line
    character(50):: temp
    integer      :: i, idx,idx_temp, nsurf
    
    line = adjustl(line)     

    idx = index(line, ' ') 
    
    nsurf = 0
    do i = 1, len(line) 
        if ((line(i:i).eq.' ').and.(i == 1)) then
            nsurf = nsurf+1
        elseif ((line(i:i).eq.' ').and.(i > 1)) then 
            if (line(i-1:i-1).ne.' ') nsurf = nsurf+1
        endif
    enddo 
            
    allocate(Cellobj%list_of_surface_IDs(nsurf))
    
    idx_temp = 1; idx = 1; temp = line
    do i = 1, nsurf
        temp = adjustl(temp(idx:len(temp))) 
        idx = index(temp, ' ')
        read(temp(1:idx-1),*) Cellobj%list_of_surface_IDs(i)
    enddo 

end subroutine 
    
    
! =============================================================================
! READ_LAT
! =============================================================================
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
    
    
! =============================================================================
! READ_CTRL reads the input for the calculation control
! =============================================================================
subroutine READ_CTRL        
    implicit none
    logical :: file_exists
    integer :: Open_Error, File_Error
    character(4) :: Card
    character:: Card_Type

    file_exists = .false.
    inquire(file="./inputfile/ctrl.inp",exist=file_exists)
    if(file_exists==.false.) then
      print *, "FATAL ERROR :: NO CTRL.INP FILE "
      stop
    end if 
    
    open(unit=rd_ctrl,file="./inputfile/ctrl.inp",status='old', action='read',iostat=Open_Error)
    Read_File : do
        read(rd_ctrl,*,iostat=File_Error) Card
        if (File_Error/=0) exit Read_File
        if (Card=="CARD" .or. Compare_String(Card,"card")) then
            backspace(rd_ctrl)
            read(rd_ctrl,*,iostat=File_Error) Card,Card_Type
            call Small_to_Capital(Card_Type)
            if (icore==score) print *, "ctrl.inp :: CARD ", Card_Type," is being read..."
            call Read_Card(rd_ctrl,Card_Type)
        end if
    end do Read_File
    close(rd_ctrl)
        
!        !Read ctrl.inp
!        open(rd_ctrl, file="./inputfile/ctrl.inp",action="read", status="old")
!    
!        ierr = 0; 
!        do while (ierr.eq.0)
!            read (rd_ctrl, FMT='(A)', iostat=ierr) line 
!            if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
!            !> option identifier 
!            j = 0
!            do while (j.le.len(line))
!                j = j+1 
!                if (line(j:j).eq.' ') exit
!                option = line(1:j)        
!            enddo 
!            select case (option)
!            case ("pop")
!                read(line(j+1:), *) n_history, n_inact, n_act
!                n_totcyc = n_inact + n_act
!                ngen = n_history
!                allocate(kprt(n_act))
!                
!            case ("energy") 
!                read(line(j+1:), *) E_mode
!                
!            case ("nugrid")
!                read(line(j+1:), *) nugrid
!            
!            case ("tally") 
!                read(line(j+1:), *) tally_switch
!                if ( tally_switch > 0 .and. icore == score ) then
!                open(prt_spec,file="flux.out",action="write",status="replace")
!                open(prt_spec,file="power.out",action="write",status="replace")
!                end if
!
!            case ("DBRC","dbrc")
!                n_iso0K = 1
!                call READ_DBRC(trim(line(j+1:)))
!            case ("entropy")
!                read(line(j+1:), *) en0(:), en1(:), nen(:)
!!            case ("CMFD") 
!!                read(line(j+1:), *) CMFD_lat, n_skip, n_acc
!!                CMFD_type = 1
!!            case ("pCMFD") 
!!                read(line(j+1:), *) CMFD_lat, n_skip, n_acc
!!                CMFD_type = 2 
!            case ("FMFD","fmfd")
!                call FMFD_INITIAL
!                call FMFD_READ(rd_ctrl)
!            case ("power") 
!                read(line(j+1:), *) Nominal_Power
!            case ("PRUP","prup")
!                call PRUP_INITIAL
!                call READ_PRUP(adjustl(line(5:)))
!            end select
!        enddo
!        close(rd_ctrl)

    if ( icore == score ) print '(A25)', '    CTRL  READ COMPLETE...' 
    
end subroutine READ_CTRL
    
! =============================================================================
! Read_Card reads the type of input card
! =============================================================================
    subroutine Read_Card(File_Number,Card_Type)
        use ENTROPY, only: en0, en1, nen
        implicit none
        integer :: i, j 
        integer,intent(in)::File_Number
        character(*),intent(inout)::Card_Type
        integer::File_Error
        character(30):: Char_Temp
        character(80):: line, lib1,lib2
        character(1)::Equal
        integer :: n
        logical :: switch
        
        File_Error=0
        n = 0 
        select case(Card_Type)
        case('A') 
            Read_Card_A : do
                if(File_Error/=0) call Card_Error(Card_Type,Char_Temp)
                read(File_Number,*,iostat=File_Error) Char_Temp
                Call Small_to_Capital(Char_Temp)
                Card_A_Inp : select case(Char_Temp)
                case("ENERGY_MODE")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, E_mode
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("NOMINAL_POWER")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, Nominal_Power
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("NUGRID")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, nugrid
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("HISTORY_PER_CYCLE")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, ngen
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("NUMBER_INACTIVE")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, n_inact
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("NUMBER_ACTIVE")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, n_act
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                    n_totcyc = n_act + n_inact
                    allocate(kprt(n_totcyc))

                case("DBRC")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, switch
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("DBRC_E_MIN")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, DBRC_E_MIN
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("DBRC_E_MAX")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, DBRC_E_MAX
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("N_ISO0K")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, n_iso0K
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                    allocate(ace0K(n_iso0K))
                case("DBRC_LIB") 
                    backspace(File_Number)
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, ace0K(1)%library
                    do i = 2, n_iso0K
                        read(File_Number,*,iostat=File_Error) ace0K(i)%library
                    enddo 

                case("ENTROPY")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, en0, en1, nen
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                case("PRUP")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, mprupon
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                    if ( mprupon ) call PRUP_INITIAL
                case("IGEN")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, ngen
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                case("DGEN")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, rampup
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                case("CRT1")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, crt1
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                case("CRT2")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, crt2
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                case("EACC")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, elength
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)

                case("FMFD")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, fmfdon
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                    if ( fmfdon ) call FMFD_INITIAL
                case("FMFD_GRID")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, fm0, fm1, nfm
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                    fm2 = fm1 - fm0
                    dfm = fm2 / dble(nfm)
                    call FMFD_ERR0
                case("FMFD_ACC")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, n_acc
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                    call FMFD_ERR0
                case("FMFD_SKIP")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, n_skip
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                    call FMFD_ERR0
                case("ONE_CMFD")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, fcr, fcz
                    if ( Equal /= "=" ) call Card_Error (Card_Type,Char_Temp)
                    ncm(1:2) = nfm(1:2) / fcr
                    ncm(3)   = nfm(3)   / fcz
                    call FMFD_ERR0
                    call FMFD_ERR1
                    cmfdon = .true.

                case("NUMBER_CMFD_SKIP")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, n_skip
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("NUMBER_CMFD_ACC")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, n_acc
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("TALLY")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, tally_switch
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                    if( tally_switch > 0 .and. icore == score ) then
                        open(prt_flux,file="flux.out",action="write",status="replace")
                        open(prt_powr,file="power.out",action="write",status="replace")
                    end if
                end select Card_A_Inp
                if (Char_Temp=="ENDA") Exit Read_Card_A
            end do Read_Card_A
        
        case('D')
            Read_Card_D : do
                read(File_Number,*,iostat=File_Error) Char_Temp
                if (Char_Temp(1:3)=="MAT" .or. Compare_String(Char_Temp(1:3),"mat")) then
                    ! add a new material slot
                    n = n+1
                    allocate(materials_temp(n))
                    if (n > 1) materials_temp(1:n-1) = materials(:) 
                    CE_mat_ptr => materials_temp(n)
                
                    if(File_Error/=0) call Card_Error(Card_Type,Char_Temp)
                    Read_Mat : do 
                        read(File_Number,*,iostat=File_Error) Char_Temp
                        Call Small_to_Capital(Char_Temp)
                        Card_D_Inp : select case(Char_Temp)
                        case("MAT_NAME")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%mat_name
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                        case("DENSITY_GPCC")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%density_gpcc
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                        case("VOL")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%vol
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                        case("FISSIONABLE")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%fissionable
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                        case("DEPLETABLE")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%depletable
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                        case("SAB")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%sab
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                            if (CE_mat_ptr%sab) then 
                                backspace(File_Number)
                                read(File_Number,*,iostat=File_Error) & 
                                    Char_Temp, Equal, line, lib1,lib2
                                call READ_SAB_MAT(j,lib1,lib2)
                            endif
                            
                        case("N_ISO")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%n_iso
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)    
                            allocate(CE_mat_ptr%ace_idx(1:CE_mat_ptr%n_iso))
                            allocate(CE_mat_ptr%numden(1:CE_mat_ptr%n_iso)) 
                            allocate(CE_mat_ptr%temp(1:CE_mat_ptr%n_iso)) 
                            
                        case("ISOTOPES")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, line, CE_mat_ptr%numden(1)
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                            CE_mat_ptr%ace_idx(1) = find_ACE_iso_idx (ace, line)
                            do i = 2, CE_mat_ptr%n_iso
                                read(File_Number,*,iostat=File_Error) line, CE_mat_ptr%numden(i)
                                CE_mat_ptr%ace_idx(i) = find_ACE_iso_idx (ace, line)
                            enddo 
                            
                        case("DOPPLER")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CE_mat_ptr%db
                            if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                    
                        case("TEMPERATURE")
                            backspace(File_Number)
                            read(File_Number,*,iostat=File_Error) Char_Temp, Equal, line, CE_mat_ptr%temp(1)
                            do i = 2, CE_mat_ptr%n_iso
                                read(File_Number,*,iostat=File_Error) line, CE_mat_ptr%temp(i)
                            enddo 
                            CE_mat_ptr%temp = CE_mat_ptr%temp * K_B

                        end select Card_D_Inp
                    
                        if (Char_Temp(1:7)=="END_MAT") Exit Read_Mat
                        
                    enddo Read_Mat
                    if(allocated(materials)) deallocate(materials)
                    call move_alloc(materials_temp, materials)
                    
                end if
                
                if (Char_Temp=="ENDD") Exit Read_Card_D
            end do Read_Card_D            
            n_materials = n
            
        case('E') ! Depletion input 
            Read_Card_E : do
                if(File_Error/=0) call Card_Error(Card_Type,Char_Temp)
                read(File_Number,*,iostat=File_Error) Char_Temp
                Call Small_to_Capital(Char_Temp)
                Card_E_Inp : select case(Char_Temp)
                ! 01_01. DO_BURN Title
                case("DO_BURN")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, do_burn
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                !case("REAL_POWER")    
                !    backspace(File_Number)
                !    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, RealPower
                !    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                ! 01_02. Read Data Format
                case("MATRIX_EXPONENTIAL_SOLVER")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, Matrix_Exponential_Solver
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                ! 01_03. Read Energy Group
                case("CRAM_ORDER")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, CRAM_ORDER
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("NSTEP_BURNUP")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, NSTEP_BURNUP
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                    allocate(burn_step(0:NSTEP_BURNUP))
                case("BURNUP_TIME")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal, burn_step(1)
                    do i = 2, NSTEP_BURNUP
                        read(File_Number,*,iostat=File_Error) burn_step(i)
                    enddo 
                    burn_step(0) = 0.0d0
                    burn_step = burn_step * 86400.d0 !Unit in [sec]                
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                case("LIBRARY_PATH")
                    backspace(File_Number)
                    read(File_Number,*,iostat=File_Error) Char_Temp, Equal
                    if(Equal/="=") call Card_Error(Card_Type,Char_Temp)
                    read(File_Number,'(A)',iostat=File_Error) dep_lib(:)
                    dep_lib = adjustl(dep_lib)
                end select Card_E_Inp
                if (Char_Temp=="ENDE") Exit Read_Card_E
            end do Read_Card_E
            
        case default 
            if (icore==score) print *, 'No such card type defined ::', card_type
            stop
        end select
        
        RealPower = Nominal_Power
        
    end subroutine Read_Card

    
    subroutine read_depletion 
        
        logical :: file_exists
        integer :: Open_Error, File_Error
        character(4)::Card
        character::Card_Type    

        file_exists = .false.
        inquire(file="./inputfile/depletion.inp",exist=file_exists)
        if(file_exists==.false.) then
          do_burn = .false.
          return
        end if 
        
        open(unit=rd_dep,file="./inputfile/depletion.inp",status='old', action='read',iostat=Open_Error)
        Read_File : do
            read(rd_dep,*,iostat=File_Error) Card
            if (File_Error/=0) exit Read_File
            if (Card=="CARD" .or. Compare_String(Card,"card")) then
                backspace(rd_dep)
                read(rd_dep,*,iostat=File_Error) Card,Card_Type
                call Small_to_Capital(Card_Type)
                if (icore==score) print *, "depletion.inp :: CARD ", Card_Type," is being read..."
                call Read_Card(rd_dep,Card_Type)
            end if
        end do Read_File
        close(rd_dep)
        
    end subroutine read_depletion

! =========================================================================
! FMFD_INITIAL
! =========================================================================
subroutine FMFD_INITIAL
    use FMFD,   only: fmfdon, n_acc
    implicit none

    fmfdon = .true.
    n_acc  = 2
    n_skip = 1

end subroutine

! =========================================================================
! FMFD_ERRS
! =========================================================================
subroutine FMFD_ERR0
    if ( .not. fmfdon ) then
    print*, "FMFD method is not ready"
    stop
    end if
end subroutine

subroutine FMFD_ERR1
    if ( mod(nfm(1),fcr) /= 0 .or. mod(nfm(3),fcz) /= 0 ) then
        print*, "CMFD mesh grid /= FMFD mesh grid"
        stop
    end if
end subroutine

!! =========================================================================
!! FMFD_READ
!! =========================================================================
!subroutine FMFD_READ(rd)
!    use FMFD_HEADER, only: fm0, fm1, fm2, nfm, dfm, ncm, cmfdon, fcr, fcz
!    implicit none
!    integer, intent(in):: rd
!    integer :: j
!    integer :: ierr
!    character(100) :: line
!    character(50)  :: option
!
!    ierr = 0
!    do while (ierr.eq.0)
!        read (rd, FMT='(A)', iostat=ierr) line 
!        if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
!        !> option identifier 
!        j = 0
!        do while (j.le.len(line))
!            j = j+1 
!            if (line(j:j).eq.' ') exit
!            option = line(1:j)        
!        enddo 
!        select case (option)
!        case ("grid")   ! essential
!            read(line(j+1:), *) fm0(:), fm1(:), nfm(:)
!            fm2(:) = fm1(:) - fm0(:)
!            dfm(:) = fm2(:) / nfm(:)
!        case ("acc")    ! optional 
!            read(line(j+1:), *) n_acc
!        case ("CMFD","cmfd")
!            cmfdon = .true.
!            read(line(j+1:), *) fcr, fcz
!            ncm(1:2) = nfm(1:2) / fcr
!            ncm(3) = nfm(3) / fcz
!
!            if ( mod(nfm(1),fcr) /= 0 .or. mod(nfm(3),fcz) /= 0 ) then
!                print*, "CMFD mesh grid /= FMFD mesh grid"
!                stop
!            end if
!        case default
!            backspace(rd)
!            return
!        end select
!    enddo
!
!end subroutine

! =========================================================================
! PRUP_INTIAL
! =========================================================================
subroutine PRUP_INITIAL
    use ENTROPY, only: rampup, ccrt, scrt, elength, mprupon, shannon, &
                       dshannon, entrp1, entrp2, genup
    implicit none
    genup    = .true.
    rampup   = ngen
    ccrt     = 1
    scrt     = 1
    elength  = 2
    n_inact  = 1000
    n_totcyc = n_inact + n_act
    allocate(shannon(2*elength))
    shannon  = 0
    dshannon = 0
    entrp1   = 0
    entrp2   = 0

end subroutine

! =========================================================================
! SET_PRUP
! =========================================================================
subroutine SET_PRUP
    use ENTROPY, only: rampup, ccrt, scrt, crt1, crt2, crt3, elength, &
                       shannon
    implicit none

    ! type of criteria
    if ( crt1 > 0 ) then
        ccrt = 2
    end if
    if ( crt2 > 0 ) then
        scrt = 2
    end if

    ! initial criteria
    if ( ccrt == 1 ) crt1 = 9D-2/sqrt(dble(ngen))
    if ( scrt == 1 ) crt2 = 3D-2/sqrt(dble(ngen))
    crt3 = crt2

end subroutine

! =========================================================================
! READ_DBRC
! =========================================================================
subroutine READ_DBRC(line)
    implicit none
    character(len=*):: line
    character(len=80):: left
    integer:: j, i

    ! minimum E
    j = index(line(:),' ')
    read(line(1:j-1),*) DBRC_E_min

    ! maximum E
    line = adjustl(line(j+1:)); j = index(line(:),' ')
    read(line(1:j-1),*) DBRC_E_max

    ! first library
    line = adjustl(line(j+1:)); j = index(line(:),' ')
    left(:) = line(:)

    ! other libraries
    do
        line = adjustl(line(j+1:))
        if ( len_trim(line) == 0 ) then
            exit
        else
            j = index(line(:),' ')
            n_iso0K = n_iso0K + 1
        end if
    end do

    ! library reading
    allocate(ace0K(n_iso0K))
    read(left,*) (ace0K(i)%library, i = 1, n_iso0K)

end subroutine

subroutine read_CE_mat
     
    logical :: file_exists
    integer :: Open_Error, File_Error
    character(4)::Card
    character::Card_Type    
    
    file_exists = .false.
    inquire(file="./inputfile/CE_mat.inp",exist=file_exists)
    if(file_exists==.false.) then
        if (icore==score) print *, "INPUT ERROR :: NO CE_mat.inp file"
        stop 
    end if 
    open(unit=rd_mat,file="./inputfile/CE_mat.inp",status='old', action='read',iostat=Open_Error)
    Read_File : do
        read(rd_mat,*,iostat=File_Error) Card
        if (File_Error/=0) exit Read_File
        if (Card=="CARD" .or. Compare_String(Card,"card")) then
            backspace(rd_mat)
            read(rd_mat,*,iostat=File_Error) Card,Card_Type
            call Small_to_Capital(Card_Type)
            if (icore==score) print *, "CE_mat.inp :: CARD ", Card_Type," is being read..."
            call Read_Card(rd_mat,Card_Type)
        end if
    end do Read_File
    close(rd_mat)
    
    if ( icore == score ) print '(A27)', '    CE MAT READ COMPLETE...' 
    
end subroutine

! =============================================================================    
! READ_SAB_MAT reads which isotope is considered by thermal scattering; S(a,b)
! =============================================================================
subroutine READ_SAB_MAT(j,lib1,lib2)
    integer, intent(inout):: j
    !character(*), intent(inout):: line
    character(80):: lib1    ! which isotope
    character(80):: lib2    ! which library
    integer:: i, k

    ! find a isotope and connect to the corresponding library
    do i = 1, num_iso
        if ( trim(ace(i)%library) == trim(lib1) ) then
            do k = 1, sab_iso
                if ( trim(sab(k)%library) == trim(lib2) ) then
                    ace(i)%sab_iso = k
                end if
            end do
        end if
    end do

end subroutine READ_SAB_MAT

    
    subroutine read_inventory 
        implicit none 
        integer :: i, j, iso, ierr
        character(50) :: mat_id, option 
        
        ! ======================================== !
        !                 ACE read start
        ! ======================================== !
        open(rd_inven, file="./inputfile/inventory.inp",action="read", status="old")
        read(rd_inven,'(A)') library_path(:)
        read(rd_inven,*) num_iso
        allocate(ace(1:num_iso))

        do iso =1, num_iso
            read(rd_inven,*) i, ace(iso)%library
        end do
        close(rd_inven)

        call SET_SAB
        call SET_DBRC
        call set_ace
        if(icore==score) print '(A29)', '    ACE Lib. READ COMPLETE...' 
                
    
    end subroutine


! =============================================================================
! SET_SAB
! =============================================================================
subroutine SET_SAB
    integer:: iso
    integer:: dummy
    integer:: i

    sab_iso = 0
    do iso = 1, num_iso
    if ( IS_SAB(trim(ace(iso)%library)) ) sab_iso = sab_iso + 1
    end do

    if ( allocated(ace) ) deallocate(ace)

    open(rd_inven, file="./inputfile/inventory.inp",action="read", status="old")
    read(rd_inven,'(A)') library_path(:)
    read(rd_inven,*) dummy
    num_iso = num_iso - sab_iso
    allocate(ace(1:num_iso))
    allocate(sab(1:sab_iso))

    do iso = 1, num_iso
        read(rd_inven,*) i, ace(iso)%library
    end do
    do iso = 1, sab_iso
        read(rd_inven,*), i, sab(iso)%library
    end do
    close(rd_inven)

end subroutine

! =============================================================================
! SET_DBRC
! =============================================================================
subroutine SET_DBRC
    integer:: ii, jj

    do ii = 1, n_iso0K
    do jj = 1, num_iso
        if ( trim(ace(jj)%library(1:5))  &
            == trim(ace0K(ii)%library(1:5)) ) then
            ace(jj)%resonant = ii
            exit
        end if
    end do
    end do

end subroutine



function IS_SAB(name_of_lib) result(lib_type)
    character(len=*), intent(in):: name_of_lib  ! name of library
    logical:: lib_type  ! type of library
    integer:: length    ! length of character

    length = len_trim(name_of_lib)
    lib_type = .false.
    if ( name_of_lib(length-5:length-5) == 't' ) lib_type = .true. 

!    select case(name_of_lib(length-5:length-5))
!    case('c'); lib_type = 1 ! continuous energy
!    case('t'); lib_type = 4 ! thermal scattering; S(a,b)
!    end select

end function






    
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
                TallyFlux(:) =0; TallyPower(:)=0
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



    ! ===========================================================================================================
    ! ===========================================================================================================
    ! ===========================================================================================================
    ! ===========================================================================================================
    ! ===========================================================================================================
!    subroutine read_CE_mat
!        implicit none 
!        integer :: i, j,i_iso, n_mat, ierr, idx
!        character(80) :: iso_id, line, option, mat_name
!        
!        open(rd_mat, file="./inputfile/CE_mat.inp",action="read", status="old")
!        ierr = 0; n_mat = 0 
!        do while (ierr.eq.0)
!            read (rd_mat, FMT='(A)', iostat=ierr) line 
!            if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
!            
!            j = 0; 
!            do while (j.le.len(line))
!                j = j+1 
!                if (line(j:j).eq.' ') exit
!                option = line(1:j)        
!            enddo 
!            
!            
!            select case (option)
!            case ("mat")
!                n_mat = n_mat+1
!                allocate(materials_temp(n_mat))
!                if (n_mat > 1) materials_temp(1:n_mat-1) = materials(:) 
!                CE_mat_ptr => materials_temp(n_mat)
!                
!                line = line(j+1:);    j = index(line(:),' ')
!                read(line(1:j-1),*) CE_mat_ptr%mat_name
!                
!                line = adjustl(line(j:)); j = index(line(:),' ')
!                read(line(1:j-1),*) CE_mat_ptr%density_gpcc
!                
!                line = adjustl(line(j:)); j = index(line(:),' ')
!                read(line(1:j-1),*) CE_mat_ptr%n_iso
!
!                line = adjustl(line(j:)); j = index(line(:),' ')
!                if ( line(1:j-1) == "sab" ) call READ_SAB_MAT(j,line)
!                
!                allocate(CE_mat_ptr%ace_idx(CE_mat_ptr%n_iso))
!                allocate(CE_mat_ptr%numden(CE_mat_ptr%n_iso)) 
!                
!                do i_iso = 1, CE_mat_ptr%n_iso 
!                    read(rd_mat, *) line, CE_mat_ptr%numden(i_iso)
!                    line = adjustl(line); j = index(line(:),' ')
!                    read(line(1:j-1),*) iso_id
!                    CE_mat_ptr%ace_idx(i_iso) = find_ACE_iso_idx (ace, iso_id)
!                enddo 
!                
!                if(allocated(materials)) deallocate(materials)
!                call move_alloc(materials_temp, materials)            
!                
!            end select
!        enddo
!        close(rd_mat)
!        if ( icore == score ) print '(A27)', '    CE MAT READ COMPLETE...' 
!        
!    end subroutine
!
!! =============================================================================    
!! READ_SAB_MAT reads which isotope is considered by thermal scattering; S(a,b)
!! =============================================================================
!subroutine READ_SAB_MAT(j,line)
!    integer, intent(inout):: j
!    character(*), intent(inout):: line
!    character(80):: lib1    ! which isotope
!    character(80):: lib2    ! which library
!    integer:: i, k
!
!
!    CE_mat_ptr%sab = .true.
!
!    line = adjustl(line(j:)); j = index(line(:),' ')
!    lib1 = line(1:j-1)
!    line = adjustl(line(j:)); j = index(line(:),' ')
!    lib2 = line(1:j-1)
!
!    ! find a isotope and connect to the corresponding library
!    do i = 1, num_iso
!    if ( trim(ace(i)%library) == trim(lib1) ) then
!    do k = 1, sab_iso
!    if ( trim(sab(k)%library) == trim(lib2) ) then
!        ace(i)%sab_iso = k
!    end if
!    end do
!    end if
!    end do
!
!end subroutine READ_SAB_MAT
!
!    
!    subroutine read_inventory 
!        implicit none 
!        integer :: i, j, iso, ierr
!        character(50) :: mat_id, option 
!        
!        ! ======================================== !
!        !                 ACE read start
!        ! ======================================== !
!        open(rd_inven, file="./inputfile/inventory.inp",action="read", status="old")
!        read(rd_inven,'(A)') library_path(:)
!        read(rd_inven,*) num_iso
!        allocate(ace(1:num_iso))
!
!        do iso =1, num_iso
!            read(rd_inven,*) i, ace(iso)%library
!        end do
!        close(rd_inven)
!
!        call SET_SAB
!        call SET_DBRC
!        call set_ace
!        if(icore==score) print '(A29)', '    ACE Lib. READ COMPLETE...' 
!                
!    
!    end subroutine
!
!
!! =============================================================================
!! SET_SAB
!! =============================================================================
!subroutine SET_SAB
!    integer:: iso
!    integer:: dummy
!    integer:: i
!
!    sab_iso = 0
!    do iso = 1, num_iso
!    if ( IS_SAB(trim(ace(iso)%library)) ) sab_iso = sab_iso + 1
!    end do
!
!    if ( allocated(ace) ) deallocate(ace)
!
!    open(rd_inven, file="./inputfile/inventory.inp",action="read", status="old")
!    read(rd_inven,'(A)') library_path(:)
!    read(rd_inven,*) dummy
!    num_iso = num_iso - sab_iso
!    allocate(ace(1:num_iso))
!    allocate(sab(1:sab_iso))
!
!    do iso = 1, num_iso
!        read(rd_inven,*) i, ace(iso)%library
!    end do
!    do iso = 1, sab_iso
!        read(rd_inven,*), i, sab(iso)%library
!    end do
!    close(rd_inven)
!
!end subroutine
!
!! =============================================================================
!! SET_DBRC
!! =============================================================================
!subroutine SET_DBRC
!    integer:: ii, jj
!
!    do ii = 1, n_iso0K
!    do jj = 1, num_iso
!        if ( trim(ace(jj)%library(1:5))  &
!            == trim(ace0K(ii)%library(1:5)) ) then
!            ace(jj)%resonant = ii
!            exit
!        end if
!    end do
!    end do
!
!end subroutine
!
!
!
!function IS_SAB(name_of_lib) result(lib_type)
!    character(len=*), intent(in):: name_of_lib  ! name of library
!    logical:: lib_type  ! type of library
!    integer:: length    ! length of character
!
!    length = len_trim(name_of_lib)
!    lib_type = .false.
!    if ( name_of_lib(length-5:length-5) == 't' ) lib_type = .true. 
!
!!    select case(name_of_lib(length-5:length-5))
!!    case('c'); lib_type = 1 ! continuous energy
!!    case('t'); lib_type = 4 ! thermal scattering; S(a,b)
!!    end select
!
!end function
!    
!
!subroutine read_tally
!    use TALLY, only: tally1, tally2
!    implicit none
!    integer :: i, j, k, idx, n, level, i_univ, i_lat, i_pin, n_pin
!    integer :: a,b,c, i_xyz(3), i_save
!    integer :: ierr
!    character(100) :: line
!    character(50)  :: option, temp, test
!    type(CoordStruct), pointer:: TC
!    
!    
!    !Read tally.inp
!    open(rd_tally, file="./inputfile/tally.inp",action="read", status="old")
!
!    ierr = 0;
!    do while (ierr.eq.0)
!        read (rd_tally, FMT='(A)', iostat=ierr) line 
!        if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle  
!        !> option identifier             
!        j = index(adjustl(line),' ') -1 
!        option = line(1:j)
!        
!        if (trim(option) /= "tally") then 
!            print *, 'ERROR Tally Read :: tally.inp     ', option
!            stop
!        end if
!
!        line = adjustl(line(j+1:))
!        j = index(line,' ')-1
!        option = trim(line(1:j))
!
!        select case (option)
!        case ("cell")
!            read(line(j+1:), *) n
!            allocate(TallyCoord(1:n))
!            allocate(TallyFlux(1:n))
!            allocate(TallyPower(1:n))
!            allocate(tally1(1:n))
!            allocate(tally2(1:n))
!            tally1 = 0; tally2 = 0
!            TallyCoord(:)%flag = 1
!            do i = 1, n
!                if ( associated(TC) ) nullify(TC)
!                TC => TallyCoord(i)
!                read (rd_tally, FMT='(A)', iostat=ierr) line 
!                if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle
!                level = 0; idx = index(line,'>')
!                do while (idx /= 0)
!                    level = level + 1 
!                    j = index(line,' '); read(line(1:j),*) temp ; line = line(j:) 
!                    line = adjustl(line);j = index(line,' ')
!                    read(line(1:j),*) i_univ;  line = line(j:)
!                    line = adjustl(line);j = index(line,' ')
!                    read(line(1:j),*) i_lat  ; line = line(j:)
!                    line = adjustl(line);j = index(line,' ')
!                    read(line(1:j),*) TC%coord(level)%lattice_x; line = line(j:)
!                    line = adjustl(line);j = index(line,' ')
!                    read(line(1:j),*) TC%coord(level)%lattice_y; line = line(j:)
!                    line = adjustl(line);j = index(line,' ')
!                    read(line(1:j),*) TC%coord(level)%lattice_z; line = line(j:)
!                    
!                    TC%coord(level)%cell     = find_cell_idx(cells, adjustl(temp))
!                    if (i_univ == 0) then 
!                        TC%coord(level)%universe = 0 
!                    else 
!                        TC%coord(level)%universe = find_univ_idx(universes, i_univ)
!                    endif
!                    
!                    if (i_lat == 0) then 
!                        TC%coord(level)%lattice = 0 
!                    else 
!                        TC%coord(level)%lattice = find_lat_idx(lattices, i_lat)
!                    endif
!                    
!                    idx = index(line,'>')
!                    line = line(idx+1:); line = adjustl(line)
!                enddo 
!                idx = index(line, 'vol')
!                line = adjustl(line(4:)); read(line(1:),*) TC%vol
!                TC%n_coord = level
!            enddo 
!            TallyFlux(:)  = 0
!            TallyPower(:) = 0 
!
!        case("pin") 
!            read(line(j+1:), *) n
!            allocate(TallyCoord(1:n)); TallyCoord(n)%n_coord = 0 
!            allocate(TallyFlux(1:n))
!            allocate(TallyPower(1:n))
!            allocate(tally1(1:n))
!            allocate(tally2(1:n))
!            tally1 = 0; tally2 = 0
!            TallyCoord(:)%flag = 0
!            i = 1; 
!            do while ( i <= n )
!                !if ( associated(TC) ) nullify(TC)
!                TC => TallyCoord(i)
!                read (rd_tally, FMT='(A)', iostat=ierr) line
!                if ((len_trim(line)==0).or.(scan(line,"%"))/=0) cycle
!
!                level = 0; idx = index(line,'>')
!                do while (idx /= 0)
!                    level = level + 1 
!                    !print *, 'level', level, line 
!                    j = index(line,' '); read(line(1:j),*) temp ; line = line(j:)                         
!                    line = adjustl(line);j = index(line,' '); read(line(1:j),*) i_univ;  line = line(j:)
!                    line = adjustl(line);j = index(line,' '); read(line(1:j),*) i_lat  ; line = line(j:)
!                    
!                    if (i_lat /= 0) i_lat = find_lat_idx(lattices, i_lat)
!                    line = adjustl(line);
!                    if (line(1:3) == 'all') then 
!                        a = lattices(i_lat)%n_xyz(1)
!                        b = lattices(i_lat)%n_xyz(2)
!                        c = lattices(i_lat)%n_xyz(3)
!                        n_pin = a*b*c
!                        i_save = i
!                        do i_pin = 1, n_pin
!                            i_xyz = getXYZ(i_pin, a,b,c)
!                            
!                            TallyCoord(i)%coord(1:level) = TallyCoord(i_save)%coord(1:level)
!                            
!                            TallyCoord(i)%coord(level)%lattice_x = i_xyz(1)
!                            TallyCoord(i)%coord(level)%lattice_y = i_xyz(2)
!                            TallyCoord(i)%coord(level)%lattice_z = i_xyz(3)
!                            
!                            if (i_lat == 0) then 
!                                TallyCoord(i)%coord(level)%lattice = 0 
!                            else 
!                                TallyCoord(i)%coord(level)%lattice  =  i_lat
!                            endif
!                            TallyCoord(i)%n_coord = level
!                            !> Reset the cell & univ designation (for pin only)
!                            TallyCoord(i)%coord(level)%cell = 0
!                            TallyCoord(i)%coord(level)%universe = 0 
!                            
!                            i = i + 1
!                            
!                        enddo 
!                        idx = 0 
!                    else
!                        line = adjustl(line);j = index(line,' ')
!                        read(line(1:j),*) TC%coord(level)%lattice_x; line = line(j:)
!                        line = adjustl(line);j = index(line,' ')
!                        read(line(1:j),*) TC%coord(level)%lattice_y; line = line(j:)
!                        line = adjustl(line);j = index(line,' ')
!                        read(line(1:j),*) TC%coord(level)%lattice_z; line = line(j:)
!                        
!                        TC%coord(level)%cell     = find_cell_idx(cells, adjustl(temp))
!                        
!                        if (i_univ == 0) then 
!                            TC%coord(level)%universe = 0 
!                        else 
!                            TC%coord(level)%universe = find_univ_idx(universes, i_univ)
!                        endif
!                        
!                        if (i_lat == 0) then 
!                            TC%coord(level)%lattice = 0 
!                        else 
!                            TC%coord(level)%lattice  =  i_lat
!                        endif
!                        
!                        idx = index(line,'>')
!                        line = line(idx+1:); line = adjustl(line)
!                        
!                        if (idx == 0) then 
!                            i = i + 1
!                            TC%n_coord = level
!                            !> Reset the cell designation (for pin only)
!                            TC%coord(level)%cell = 0
!                        endif
!                    endif 
!                enddo 
!                !> pin volume is not needed (all same)
!                TallyCoord(:)%vol = 1
!            enddo                     
!            TallyFlux(:)  = 0
!            TallyPower(:) = 0
!            
!        end select
!    enddo
!    
!    close(rd_tally)
!    if(icore==score) print '(A25)', '    TALLY READ COMPLETE...' 
!    
!end subroutine
!
!subroutine check_input_result(universes,lattices, cells,surfaces)
!    type(surface) :: surfaces(:)
!    type(lattice) :: lattices(:)
!    type(universe):: universes(0:)
!    type(cell)      :: cells(:)
!    integer :: i, j, iy, iz
!    
!    print *, ' ========== SURFACE CHECK =========='
!    print *, ' INDEX     ID      surf_type'
!    do i = 1, size(surfaces) 
!        print '(I5, A10, I10)', i, trim(surfaces(i)%surf_id), surfaces(i)%surf_type
!    enddo 
!    print *, ''
!    print *, ''
!    
!    
!    print *, ' ========== CELL CHECK =========='
!    do i = 1, size(cells) 
!        print *, 'cell number :', i 
!        print *, 'cell id:', cells(i)%cell_id!, cells(i)%univ_id, cells(i)%mat_id
!        print *, 'surf list: ', cells(i)%list_of_surface_IDs(:)
!        print *, 'neg : ', surfaces(cells(i)%neg_surf_idx(:))%surf_id
!        print *, 'pos : ', surfaces(cells(i)%pos_surf_idx(:))%surf_id
!        print *, 'operand : ',cells(i)%operand_flag
!        print *, 'fill type:', cells(i)%fill_type()
!        print *, 'material idx', cells(i)%mat_idx
!        print *, ''
!        print *, ''
!    enddo 
!    
!    
!    print *, ' ========== UNIVERSE CHECK =========='
!    do i = 0,size(universes(1:))
!        print *, ''
!        print *, 'univ_id = ', universes(i)%univ_id, '   # of cell', universes(i)%ncell 
!        do j = 1, universes(i)%ncell
!            print '(a4,I2,2A7)', 'cell',j,'    ', cells(universes(i)%cell(j))%cell_id
!        enddo 
!    enddo 
!    print *, ''
!    print *, ''
!    
!    
!    print *, ' ========== LATTICE CHECK =========='
!    do i = 1, size(lattices) 
!        print *, ''
!        print *, 'lattice id = ', lattices(i)%lat_id
!        do iz = 1, lattices(i)%n_xyz(3)
!            do iy = 1, lattices(i)%n_xyz(2)
!                print '(17I2)', (universes(lattices(i)%lat(j,iy,iz))%univ_id, j = 1, lattices(i)%n_xyz(1))
!            enddo 
!            
!            print *, '' 
!        enddo 
!    enddo
!    
!    
!    
!        
!end subroutine 
    
end module
