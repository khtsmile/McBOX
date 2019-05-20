module tracking
    use omp_lib
    use variables
    use constants
    use surface_header
    use geometry_header,    only: cell
    use geometry
    use particle_header,    only: particle
    use randoms,             only: rang
    use physics,            only: collision_MG
    use XS_header 
    use tally,                 only: TallyCoord, TallyFlux, FindTallyBin
    use ace_xs,             only: getMacroXS
    use material_header,    only: materials
    use ace_reactions,         only: collision_CE
    
    implicit none

contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

    subroutine transport(p)
    
        type(Particle), intent(inout) :: p
        
        integer :: i 
        integer :: j                      ! coordinate level
        integer :: next_level             ! next coordinate level to check
        integer :: surface_crossed        ! surface which particle is on
        integer :: lattice_translation(3) ! in-lattice translation vector
        integer :: n_event                ! number of collisions/crossings
        real(8) :: d_boundary             ! distance to nearest boundary
        real(8) :: d_collision            ! sampled distance to collision
        real(8) :: distance               ! distance particle travels
        logical :: found_cell             ! found cell which particle is in?
        real(8) :: macro_xs(4)
        real(8) :: xyz(3)
        integer :: i_cell, i_bin
        
        
        found_cell = .false.
        
        call find_cell(p, found_cell, i_cell)
        
        call distance_to_boundary(p, d_boundary, surface_crossed)    
        
        ! Sample a distance to collision
        if (E_mode == 0) then 
            d_collision = -log(rang())/(sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
        elseif (E_mode == 1) then 
            macro_xs = getMacroXS(materials(p%material), p%E)
            d_collision = -log(rang())/macro_xs(1)
        endif 
        
        distance = min(d_boundary, d_collision)
                
        !> Flux Tally ===========================================================================
        !if (tally_switch>0) then 
        !    i_bin = FindTallyBin(p)
        !    if (i_bin > 0) then 
        !        !$omp atomic
        !        TallyFlux(i_bin) = TallyFlux(i_bin) + distance*p%wgt/(real(n_history,8)*TallyCoord(i_bin)%vol)
        !    endif
        !endif 
        !> ======================================================================================
        
        
        do j = 1, p % n_coord
            p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
        enddo
        !write(wt_coord, *) p % coord(1) % xyz
        !print *, p%E
        if (d_boundary > d_collision) then ! collision 
            if (E_mode == 0) call collision_MG(p)
            if (E_mode == 1) call collision_CE(p)
            
        else ! surface crossing
            !write(wt_coord, *) p % coord(1) % xyz
            !write(wt_coord, *) p%E, p%wgt
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
            !print *, 'escaped'
        elseif (surfaces(surface_crossed)%bc == 2) then !> Reflective BC 
        
            !i_cell_prev = p%coord(p%n_coord)%cell
            !print *,'xyz', p%coord(1)%xyz
            p % n_coord = 1
            p % coord(1) % xyz(:) = p % coord(1) % xyz(:) - TINY_BIT * p % coord(1) % uvw(:)
            !call find_cell(p, found, i_cell)
            !if(i_cell == 0 ) then 
            !    print *, 'leak'
            !    print *, p % coord(1) % xyz(:)
            !    stop 
            !endif
            !print *,'uvw', p%coord(1)%uvw
            
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
    
end module tracking
