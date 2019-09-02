module tracking
    use omp_lib
    use variables
    use constants
    use surface_header
    use geometry_header,    only: cell, lattices
    use geometry
    use particle_header,    only: particle
    use randoms,            only: rang
    use physics,            only: collision_MG, collision_MG_DT
    use XS_header 
    use tally,              only: TallyCoord, TallyFlux, FindTallyBin, TallyPower
    use ace_xs,             only: getMacroXS
    use material_header,    only: materials
    use ace_reactions,      only: collision_CE
    use FMFD,               only: FMFD_distance, FMFD_TRK, FMFD_COL, &
                                FMFD_SURF, fmfdon
    
    implicit none
    logical:: inside1
    real(8):: xyz0(3)
    real(8):: d0, d1


contains
 
!===============================================================================
! TRANSPORT - the main logic for moving a particle through geometry.
!===============================================================================

subroutine transport(p,cyc)

    type(Particle), intent(inout) :: p
    integer, intent(in):: cyc
    
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: d_FMFD                 ! distance to FMFD grid
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: macro_xs(5)
    real(8) :: xyz(3)
    integer :: i_cell, i_bin(4), i_lat, i_surf
    integer :: i_xyz(3), idx_xyz
    logical :: inside_FMFD
    integer :: income_FMFD
    real(8) :: ddiff
    logical :: fm_crossed
    
    found_cell = .false.
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)

    
    !> Surface distance(boundary)
    call distance_to_boundary(p, d_boundary, surface_crossed)
    
    !> Sample a distance to collision
    if (E_mode == 0) then 
        macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) &
                    + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
        macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
        macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
        d_collision = -log(rang())/macro_xs(1)
    elseif (E_mode == 1) then 
        macro_xs = getMacroXS(materials(p%material), p%E)
        d_collision = -log(rang())/macro_xs(1)
    endif 
    
    ! ===================================================
    !> CMFD distance 
    d_FMFD = INFINITY
    if ( fmfdon ) call FMFD_DISTANCE (p,i_xyz,d_FMFD,inside_FMFD,income_FMFD,i_surf)
    
    !> minimum distance
    ddiff = abs(d_boundary-d_FMFD)/d_boundary
    if ( ddiff < TINY_BIT ) then
        d_FMFD = d_boundary
        fm_crossed = .true.
    else if ( d_boundary < 5E-5 .and. ddiff < 1E-8 ) then
        d_FMFD = d_boundary
        fm_crossed = .true.
    else
        fm_crossed = .false.
    end if
    distance = min(d_boundary, d_collision, d_FMFD)
!    print*, "dst", d_boundary, d_collision, d_FMFD
!    print*, "xyz0", p%coord(1)%xyz
!    print*, "xyz1", p%coord(1)%xyz+distance*p%coord(1)%uvw
!    print*, "ijk", i_xyz, inside_FMFD
!    if ( .not. fm_crossed .and. inside_FMFD ) then
!        print*, d_boundary, d_FMFD, abs(d_boundary-d_FMFD)/d_boundary
!        pause
!    end if
!    pause
!    print*, "uvw", sign(1.,p%coord(1)%uvw(1)), &
!                   sign(1.,p%coord(1)%uvw(2)), &
!                   sign(1.,p%coord(1)%uvw(3))

!    if ( inside1 == .false. .and. inside_FMFD == .true. ) then
!    !if ( d0 /= d1 .and. d1 < 1D100 ) then
!        print*, inside1, inside_FMFD
!        print*, "xyz", p%coord(1)%xyz
!        print*, "uvw", sign(1.,p%coord(1)%uvw(1)), sign(1.,p%coord(1)%uvw(2)), sign(1.,p%coord(1)%uvw(3))
!        print*, "ddd", d_boundary, d_FMFD
!        print*, "d1", p%coord(1)%xyz+d_boundary*p%coord(1)%uvw
!        print*, "d2", p%coord(1)%xyz+d_FMFD*p%coord(1)%uvw
!        print*, "=?", d0 == d1
!        !pause
!    !end if
!    end if
!    print'(a,l,3f15.7,3x,a,2es20.12)', "*", inside1, xyz0, "/", d_boundary, d_FMFD
!    d0 = d_boundary
!    d1 = d_FMFD
!    inside1 = inside_FMFD
!    xyz0 = p%coord(1)%xyz

    !> Track-length estimator
    !$omp atomic 
    k_tl = k_tl + distance*p%wgt*macro_xs(4) 
    
    !> Flux Tally ===========================================================================
    !if ( tally_switch > 0 .and. cyc > n_inact ) then 
    if ( tally_switch > 0 ) then 
        i_bin = FindTallyBin(p)
        if ( i_bin(1) > 0 ) then 
            !$omp atomic
            TallyFlux(i_bin(1)) = TallyFlux(i_bin(1)) + distance*p%wgt
    !> ==== Power Tally =====================================================================
            !$omp atomic
            TallyPower(i_bin(1)) = TallyPower(i_bin(1)) + distance*p%wgt*macro_xs(5)
        endif
    endif 
    
    
    !> CMFD Tally (track length) 
    if ( fmfdon .and. inside_FMFD ) call &
    FMFD_TRK(p%wgt,distance,macro_xs,i_xyz,i_bin(2:4),p%coord(1)%xyz)

    !> Advance particle
    do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
    enddo
    
    if ( distance == d_collision ) then ! collision 
        if (E_mode == 0) then 
            call collision_MG(p)
        else !(E_mode == 1) 
            call collision_CE(p)
        endif
        !if ( fmfdon .and. inside_FMFD ) call FMFD_COL(p%wgt,macro_xs,i_xyz)

    elseif  ( distance == d_FMFD ) then 
        call FMFD_SURF(inside_FMFD, income_FMFD,i_surf, i_xyz, p%wgt, &
                       surfaces(surface_crossed)%bc)
        if ( fm_crossed ) then
            call cross_surface(p, surface_crossed)
        else
            p%coord(1)%xyz = p%coord(1)%xyz + TINY_BIT * p%coord(1)%uvw
        end if
    else
        call cross_surface(p, surface_crossed)
    endif
    
end subroutine transport

!===============================================================================
! CROSS_SURFACE handles all surface crossings, whether the particle leaks out of
! the geometry, is reflected, or crosses into a new lattice or cell
!===============================================================================

subroutine cross_surface(p, surface_crossed)
    type(Particle), intent(inout) :: p
    
    integer :: surface_crossed
    real(8) :: xyz(3)     ! Saved global coordinate
    integer :: i_surface  ! index in surfaces
    logical :: rotational ! if rotational periodic BC applied
    logical :: found      ! particle found in universe?
    class(Surface), pointer :: surf
    class(Surface), pointer :: surf2 ! periodic partner surface
    integer ::i, i_cell, i_cell_prev
    
    
    !p % n_coord = 1
    !call find_cell(p, found)
    if (surfaces(surface_crossed)%bc == 1) then     !> Vacuum BC
        p % wgt = 0 
        p % alive = .false.
    elseif (surfaces(surface_crossed)%bc == 2) then !> Reflective BC 
        p % n_coord = 1
        p % coord(1) % xyz(:) = p % coord(1) % xyz(:) - TINY_BIT * p % coord(1) % uvw(:)
        !call find_cell(p, found, i_cell)
        
        call reflective_bc(p%coord(1)%uvw, p%coord(1)%xyz, surface_crossed)
        !p%last_material = p%material
        !p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
        !print *, 'reflected'
        
    else
        p % n_coord = 1
        p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
        !print *, 'pass',surfaces(surface_crossed)%surf_id 
    endif
    !print *, 'next', p % coord(1) % xyz
    call find_cell(p, found)
    
    !print *, cells(p%coord(1)%cell)%cell_id
    !if (p%coord(p%n_coord)%cell == 0 ) print *, 'here'
    
end subroutine cross_surface



subroutine reflective_bc (uvw, xyz, surface_crossed)
    real(8), intent(inout) :: uvw(3)
    real(8), intent(in)       :: xyz(3)
    real(8) :: xyz_(3), r
    integer :: surface_crossed
    integer :: surf_type
    integer :: flag
    
    surf_type = surfaces(surface_crossed)%surf_type
    
    select case(surf_type) 
    
    case (1) !> px
        uvw(1) = -uvw(1) 
    case (2) !> py
        uvw(2) = -uvw(2) 
    case (3) !> pz
        uvw(3) = -uvw(3) 
        
    case (6) !> sqcz  (TO BE EDITTED)
        xyz_(3)   = xyz(3) 
        xyz_(1:2) = xyz(1:2) - surfaces(surface_crossed)%parmtrs(1:2) 
        r = surfaces(surface_crossed)%parmtrs(3) 
        
        if ((xyz_(2) > -r-TINY_BIT).and.(xyz_(2) < -r + TINY_BIT)) then 
            uvw(2) = -uvw(2)
        elseif ((xyz_(1) > -r-TINY_BIT).and.(xyz_(1) < -r + TINY_BIT)) then 
            uvw(1) = -uvw(1)
        elseif ((xyz_(2) > r-TINY_BIT).and.(xyz_(2) < r + TINY_BIT)) then 
            uvw(2) = -uvw(2) 
        elseif ((xyz_(1) > r-TINY_BIT).and.(xyz_(1) < r + TINY_BIT)) then 
            uvw(1) = -uvw(1) 
        else 
            print *, 'particle is not on the surface'
            stop 
        end if
        
    case (9) !> cylz
        xyz_(3)   = xyz(3) 
        xyz_(1:2) = xyz(1:2) - surfaces(surface_crossed)%parmtrs(1:2) 
        
        uvw(1) = uvw(1) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2))*xyz_(1)&
                    /(surfaces(surface_crossed)%parmtrs(3))**2
        uvw(2) = uvw(2) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2))*xyz_(2)&
                    /(surfaces(surface_crossed)%parmtrs(3))**2
        uvw(3) = uvw(3) 
        
    case (10) !> sph
        xyz_(1:3) = xyz(1:3) - surfaces(surface_crossed)%parmtrs(1:3) 
        
        uvw(1) = uvw(1) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2) + xyz_(3)*uvw(3))*xyz_(1) &
                    /(surfaces(surface_crossed)%parmtrs(4))**2
        uvw(2) = uvw(2) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2) + xyz_(3)*uvw(3))*xyz_(2) &
                    /(surfaces(surface_crossed)%parmtrs(4))**2
        uvw(3) = uvw(3) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2) + xyz_(3)*uvw(3))*xyz_(3) &
                    /(surfaces(surface_crossed)%parmtrs(4))**2
        
    end select 
    
    
end subroutine



!===============================================================================
! transport_DT handles the Woodcock Delta-tracking algorithm 
!===============================================================================
subroutine transport_DT(p) 

    type(Particle), intent(inout) :: p
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: d_FMFD                  ! distance to CMFD grid
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    
    real(8) :: macro_major              ! the global majorant cross-section
    real(8), allocatable :: macro_tot(:)
    real(8) :: macro_xs(5)
    real(8) :: xyz(3)
    integer :: i_cell, idx_surf
    integer :: i_xyz(3), idx_xyz, idx
    integer :: bc
    integer :: n_mat 
    
    
    found_cell = .false.
    i_xyz(:) = -1; idx_xyz = -1
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
    
    !> Determine macro_major 
    if (E_mode == 0) then 
        n_mat = size( XS_MG ) 
        allocate(macro_tot(n_mat))
        do i = 1, n_mat 
            macro_tot(i) = (sum(XS_MG(i)%sig_scat(p%g,:)) + XS_MG(i)%sig_abs(p%g))
        enddo 
    else 
        n_mat = size( materials ) 
        allocate(macro_tot(n_mat))
        do i = 1, n_mat 
            macro_xs = getMacroXS(materials(i), p%E)
            macro_tot(i) = macro_xs(1) 
        enddo 
    endif 
    macro_major = maxval(macro_tot, n_mat)
    deallocate(macro_tot)
    
    !> Sample a distance to collision
    d_collision = -log(rang())/macro_major
    
    !> Sample distances from special boundaries in univ 0
    idx = p % coord(1) % cell
    call cell_distance(cells(idx), p%coord(1)%xyz, p%coord(1)%uvw, surfaces, d_boundary, idx_surf)
    
    distance = min(d_boundary, d_collision)
    !print *, d_boundary, d_collision
    
    !> Determine Virtual / Real collsion OR Cross-surface
    if (d_collision < d_boundary) then ! collision 
        ! Advance particle
        p % n_coord = 1
        p % coord(1) % xyz = p % coord(1) % xyz + distance * p % coord(1) % uvw
        call find_cell(p, found_cell, i_cell)
        
        if (E_mode == 0) then 
            macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
        elseif (E_mode == 1) then 
            macro_xs = getMacroXS(materials(p%material), p%E)
        endif 
        
        ! Reject? 
        if (rang() < macro_xs(1)/macro_major ) then !> real collision
            if (E_mode == 0) then 
                !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
                call collision_MG_DT(p, macro_major)
            else 
                !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
                call collision_CE(p)
            endif
        endif 
    else
        call cross_surface(p, idx_surf)
    endif
    
    p%n_coord = 1 
    
end subroutine



end module tracking
