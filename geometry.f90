module geometry
    use constants,             only: TINY_BIT, MAX_COORD
    use surface_header 
    use geometry_header
    use particle_header
    use omp_lib
    
    implicit none
    
    
    
    contains 
    
    function cell_contains(c, p) result(in_cell)
        type(Cell), intent(in) :: c
        type(Particle), intent(in) :: p
        logical :: in_cell
        integer :: i,j, n 
        
        j = p % n_coord
        if (c%operand_flag >= 0) then   !> and 
            in_cell = .true.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),p%coord(j)%xyz) == .false.) in_cell = .false.
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), p%coord(j)%xyz) == .true.) in_cell = .false.
            enddo
        else     !> or
            in_cell = .false.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),p%coord(j)%xyz) == .true.) in_cell = .true.
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), p%coord(j)%xyz) == .false.) in_cell = .true.
            enddo
        endif
        
    end function 
    
    function cell_xyz(c, xyz) result(in_cell)
        type(Cell), intent(in) :: c
        real(8),    intent(in) :: xyz(3)
        logical :: in_cell
        integer :: i,j, n 
        
        if (c%operand_flag >= 0) then   !> and 
            in_cell = .true.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),xyz) == .false.) in_cell = .false.
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), xyz) == .true.) in_cell = .false.
            enddo
        else     !> or
            in_cell = .false.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),xyz) == .true.) in_cell = .true.
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), xyz) == .false.) in_cell = .true.
            enddo
        endif
        
    end function 
    
    
    !===============================================================================
    ! FIND_CELL determines what cell a source particle is in within a particular
    ! universe. If the base universe is passed, the particle should be found as long
    ! as it's within the geometry
    !===============================================================================
    !recursive subroutine find_cell(p, found, search_cells, cell_idx)
    recursive subroutine find_cell(p, found, cell_idx)
        type(Particle), intent(inout) :: p
        logical,        intent(inout) :: found
        !integer,        optional      :: search_cells(:)
        integer,         optional      :: cell_idx
        integer :: i                    ! index over cells
        integer :: j, k, idx            ! coordinate level index
        integer :: offset               ! instance # of a distributed cell
        integer :: distribcell_index
        integer :: i_xyz(3)             ! indices in lattice
        integer :: n                    ! number of cells to search
        integer :: i_cell               ! index in cells array
        integer :: i_univ_cell            ! index of cell in universe cell list
        integer :: i_universe           ! index in universes array

        
        do j = p % n_coord + 1, MAX_COORD
          call p % coord(j) % reset()
        enddo
        j = p % n_coord

        i_universe = p % coord(j) % universe
        
        
        if ( j == 1 ) then
            !p%univ => universes(0)
            p%coord(1)% universe = 0 
            idx = 0
        else
            idx = p % coord(j) % universe
            !p%univ => universes(idx)
            !p%coord(p%n_coord)% universe = idx
            !> coordinate translation for the new universe
            !> call last_univ_coord_tranlation(universes,p,j)
            !p%coord(j)%xyz = p%coord(j-1)%xyz에다가 연산.
            
        endif
        
        !n = p%univ%ncell; 
        n = universes(idx)%ncell; found = .false.
        !print *, 'num of cell', n
        CELL_LOOP: do i = 1, n
            ! select cells based on whether we are searching a universe or a provided
            ! list of cells (this would be for lists of neighbor cells)
            i_cell = universes(idx)%cell(i)
            ! Move on to the next cell if the particle is not inside this cell
            !print *, i, cells(i_cell)%cell_id
            !print *, p%coord(1)%xyz
            if(cell_contains(cells(i_cell), p)) then
                ! Set cell on this level
                p % coord(j) % cell = i_cell
                found = .true.
                !print *, 'in univ ',universes(idx)%univ_id , 'found in cell:', cells(i_cell)%cell_id
                !print *, i_cell
                !print *, p%coord(j)%xyz
                exit
            endif
            
        end do CELL_LOOP
        
        
                
        if ( found ) then
            associate(c => cells(i_cell))
                CELL_TYPE: if (c % filltype == FILL_MATERIAL) then
                    ! ======================================================================
                    ! AT LOWEST UNIVERSE, TERMINATE SEARCH
                    
                    p % coord(j) % cell = i_cell      ! index in cells(:) array
                    p % material = c%mat_idx
                    
                    if ( p%material == 0 ) then 
                        p%wgt = 0 
                        p%alive = .false. 
                    endif
                    
                    if ( present(cell_idx) ) cell_idx = i_cell
                    
                    
                elseif (c % filltype == FILL_UNIVERSE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
                    !print *, 'FILL_UNIVERSE : ', universes(find_univ_idx(universes, c % fill))%univ_id
                    !print *, 'cell translation', c % translation
                    ! Store lower level coordinates
                    
                    p % coord(j + 1) % xyz = p % coord(j) % xyz
                    p % coord(j + 1) % uvw = p % coord(j) % uvw
                    
                    ! Move particle to next level and set universe
                    j = j + 1
                    p % n_coord = j
                                        
                    !p % coord(j) % universe = universes%find_idx(cells(i_cell) % fill)
                    p % coord(j) % universe = find_univ_idx(universes, c % fill)
                    !print *, 'universe', c % fill
                    
                    ! Apply translation
                    if (allocated(c % translation)) then
                        p % coord(j) % xyz = p % coord(j) % xyz - c % translation
                    end if
                    
                    call find_cell(p, found, cell_idx)
                    !j = p % n_coord
                
                
                elseif (c % filltype == FILL_LATTICE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
                    !print *, 'FILL_LATTICE : ', lattices(find_lat_idx(lattices, c % fill))%lat_id 
                    
                    associate (latptr => lattices(find_lat_idx(lattices, c % fill)))
                        ! Determine lattice indices
                        !i_xyz = lat % get_indices(p % coord(j) % xyz + TINY_BIT * p % coord(j) % uvw)
                        i_xyz = lattice_coord(latptr,p % coord(j) % xyz)  !latptr%lat_pos(p % coord(j) % xyz)
                        !print *, p % coord(j) % xyz
                        !print *, 'lattice coordinate :', i_xyz(1:2)
                        ! Store lower level coordinates
                        p % coord(j + 1) % xyz = get_local_xyz(latptr, p % coord(j) % xyz, i_xyz)
                        p % coord(j + 1) % uvw = p % coord(j) % uvw
                        
                        ! set particle lattice indices
                        p % coord(j + 1) % lattice   = find_lat_idx(lattices, c % fill)
                        !print *, 'lattice   ', c%cell_id, c % fill, p % coord(j + 1) % lattice
                        p % coord(j + 1) % lattice_x = i_xyz(1)
                        p % coord(j + 1) % lattice_y = i_xyz(2)
                        p % coord(j + 1) % lattice_z = i_xyz(3)
                        
                        !print *, 'find cell', lattices%find_idx(c % fill), i_xyz(:)
                        
                        ! Set the next lowest coordinate level.
                        p % coord(j + 1) % universe = latptr % lat(i_xyz(1),i_xyz(2),i_xyz(3))
                        !print *, 'lattices universe', universes(p % coord(j + 1) % universe)%univ_id
                        
                    end associate
                    
                    ! Move particle to next level and search for the lower cells.
                    j = j + 1
                    p % n_coord = j
                    
                    call find_cell(p, found, cell_idx)
                    !j = p % n_coord
                    
                    
                end if CELL_TYPE
            end associate
        end if        
        
        
    end subroutine 
    
    !===============================================================================
    ! TRANSFORM_COORD transforms the xyz coordinate from univ1 to univ2 
    !===============================================================================
    subroutine transform_coord ()
        
        
    end subroutine
    
    
    
    
    !===============================================================================
    ! DISTANCE_TO_BOUNDARY calculates the distance to the nearest boundary for a
    ! particle 'p' traveling in a certain direction. For a cell in a subuniverse
    ! that has a parent cell, also include the surfaces of the edge of the universe.
    !===============================================================================
    subroutine distance_to_boundary(p, dist, surface_crossed)!, level)!, lattice_translation, next_level)
        type(Particle), intent(inout) :: p
        real(8),        intent(out)   :: dist
        !type(LocalCoord), intent(out) :: level(MAX_COORD) 
        integer,        intent(out)   :: surface_crossed
        !integer,        intent(out)   :: lattice_translation(3)
        !integer,        intent(out)   :: next_level
        
        integer :: idx
        integer :: i, j, j_idx
        integer :: i_xyz(3)           ! lattice indices
        integer :: level_surf_cross   ! surface crossed on current level
        integer :: level_lat_trans(3) ! lattice translation on current level
        real(8) :: xyz_t(3)           ! local particle coordinates
        real(8) :: d_lat              ! distance to lattice boundary
        real(8) :: d_surf             ! distance to surface
        real(8) :: xyz_cross(3)       ! coordinates at projected surface crossing
        real(8) :: surf_uvw(3)        ! surface normal direction

        integer :: idx_surf
        integer :: i_xyz_out(3)
        real(8) :: dist_temp
        integer :: univ_idx

        !do j = 1, MAX_COORD
        !  call level(j) % reset()
        !enddo
        
        ! inialize distance to infinity (huge)
        dist = INFINITY
        d_lat = INFINITY
        d_surf = INFINITY
        dist_temp = INFINITY
        i_xyz_out(:) = 0 
        
        !lattice_translation(:) = [0, 0, 0]

        !next_level = 0
        surface_crossed = -1 
        
        !> surface 탐색 순서
        ! 1. universe 내부 cell boundary
        ! 2. 만약 위 거리가 toolong 보다 길면 -> 상위 universe 탐색 (lattice or universe)
        p%coord(:)%dist = INFINITY            ! reset level distance
        univ_idx = 0 
        j = 0 
        LEVEL_LOOP: do 
            j = j+1

            ! get pointer to cell on this level
            idx = p % coord(j) % cell
            
            ! =======================================================================
            ! FIND MINIMUM DISTANCE TO SURFACE IN THIS CELL
            call cell_distance(cells(idx), p%coord(j)%xyz, p%coord(j)%uvw, surfaces, p%coord(j)%dist, idx_surf)
            
            ! =======================================================================
            ! FIND MINIMUM DISTANCE TO LATTICE SURFACES        
            !if ((p%coord(j) % dist > TOOLONG).and.(p % coord(j) % lattice /= NONE)) then 
            if (p % coord(j) % lattice /= NONE) then
                i_xyz(1) = p % coord(j) % lattice_x
                i_xyz(2) = p % coord(j) % lattice_y
                i_xyz(3) = p % coord(j) % lattice_z
                
                call lat_distance(lattices(p % coord(j) % lattice), surfaces, p % coord(j) % xyz, &
                                    p % coord(j) % uvw, i_xyz, dist_temp, idx_surf)
                                    
                if (dist_temp < p % coord(j) % dist) p % coord(j) % dist = dist_temp
                
            endif 
            
            
            if (j == 1) then 
                j_idx = 1
                surface_crossed = idx_surf
            elseif ((p%coord(j)%dist < p%coord(j_idx)%dist+TINY_BIT).and.&
                    (p%coord(j)%dist > p%coord(j_idx)%dist-TINY_BIT)) then  !> similar 
                if (surfaces(j)%bc == 2)  then 
                    surface_crossed = idx_surf
                    j_idx = j
                endif 
            elseif ((p%coord(j)%dist < p%coord(j_idx)%dist)) then 
                surface_crossed = idx_surf
                j_idx = j
            endif 
            
            if (j >= p%n_coord ) exit
        enddo LEVEL_LOOP
        
        dist = p%coord(j_idx)%dist
        
        
        !dist = p%coord(1)%dist; level(1) = p%coord(1)
        !do i = 1, p%n_coord 
        !    if (dist > p%coord(i)%dist) then 
        !        dist = p%coord(i)%dist
        !        level(1:i) =p%coord(1:i) 
        !    enddo
        !enddo
        
        !dist = minval(p%coord(:)%dist)
        
        
    end subroutine
    
    
    
    
    !===============================================================================
    ! NEIGHBOR_LISTS builds a list of neighboring cells to each surface to speed up
    ! searches when a cell boundary is crossed.
    !===============================================================================
    subroutine neighbor_lists()
        
    end subroutine 
    
    
    
    recursive subroutine find_cell_xyz(xyz, idx_univ, cell_idx)
        real(8),         intent(inout) :: xyz(3)
        logical                          :: found
        integer,         optional      :: cell_idx
        integer :: i                    ! index over cells
        integer :: j, k, idx            ! coordinate level index
        integer :: offset               ! instance # of a distributed cell
        integer :: distribcell_index
        integer :: i_xyz(3)             ! indices in lattice
        integer :: n                    ! number of cells to search
        integer :: i_cell               ! index in cells array
        integer :: i_univ_cell            ! index of cell in universe cell list
        integer :: i_universe           ! index in universes array
        integer,         intent(inout) :: idx_univ
        
        
        n = universes(idx_univ)%ncell; found = .false.
        !print *, 'num of cell', n
        CELL_LOOP: do i = 1, n

            ! select cells based on whether we are searching a universe or a provided
            ! list of cells (this would be for lists of neighbor cells)
            i_cell = universes(idx_univ)%cell(i)
            ! Move on to the next cell if the particle is not inside this cell
            if(cell_xyz(cells(i_cell), xyz)) then
                found = .true.
                exit
            endif
            !print *, 'cell not found' 
            !stop
        end do CELL_LOOP
        
                
        if ( found ) then
            associate(c => cells(i_cell))
                CELL_TYPE: if (c % fill_type() == FILL_MATERIAL) then
                    ! ======================================================================
                    ! AT LOWEST UNIVERSE, TERMINATE SEARCH
                    !print *, 'FILL_MATERIAL : ', c%cell_id,  c%mat_idx
                    
                    if (present(cell_idx)) then 
                        cell_idx = i_cell
                    endif
                    
                elseif (c % fill_type() == FILL_UNIVERSE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
                    !print *, 'FILL_UNIVERSE : ', universes(universes%find_idx(c % fill))%univ_id
                    !print *, 'cell translation', c % translation
                    ! Store lower level coordinates
                    !p % coord(j + 1) % xyz = p % coord(j) % xyz
                    !p % coord(j + 1) % uvw = p % coord(j) % uvw

                    ! Move particle to next level and set universe
                    idx_univ = find_univ_idx(universes, c % fill)!universes%find_idx(c % fill)

                    ! Apply translation
                    if (allocated(c % translation)) then
                        xyz = xyz - c % translation
                    end if

                    call find_cell_xyz(xyz, idx_univ, cell_idx)
            
            
                elseif (c % fill_type() == FILL_LATTICE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
                    !print *, 'FILL_LATTICE : ', lattices(lattices%find_idx(c % fill))%lat_id 
                    
                    associate (latptr => lattices(find_lat_idx(lattices,c % fill)))
                        ! Determine lattice indices
                        !i_xyz = lat % get_indices(p % coord(j) % xyz + TINY_BIT * p % coord(j) % uvw)
                        i_xyz = lattice_coord(latptr, xyz)!latptr%lat_pos(xyz)
                        !print *, 'lattice coordinate :', i_xyz(1:2)
                        ! Store lower level coordinates
                        xyz = get_local_xyz(latptr , xyz, i_xyz)
                        ! Set the next lowest coordinate level.
                        !print *, latptr % lat(i_xyz(1),i_xyz(2))
                        !print *, idx_univ
                        idx_univ = latptr % lat(i_xyz(1),i_xyz(2),i_xyz(3))
                        
                        !print *, idx_univ
                    end associate
                    
                    ! Move particle to next level and search for the lower cells.                    
                    call find_cell_xyz(xyz, idx_univ, cell_idx)
                    
                    
                end if CELL_TYPE
            end associate
            
        else 
            print *, 'cell not found'
            stop
        end if        
        
        
    end subroutine     
    
    
    
end module
