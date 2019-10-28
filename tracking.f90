module tracking
    use omp_lib
    use variables
    use constants
    use surface_header
    use geometry_header,    only: cell, lattices, cells, cell_distance
    use geometry,           only: cross_surface, distance_to_boundary, find_cell
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
    use DEPLETION_MODULE,   only: tally_burnup
    use VRC,                only: trace_psudoray
    use TH_HEADER,          only: th_on
    use TEMPERATURE,        only: TH_INSIDE, TH_COL
    
    implicit none

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
    integer :: i_xyz(3), idx_xyz, j_xyz(3)
    logical :: inside_FMFD
    integer :: income_FMFD
    logical :: inside_th
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

    ! =========================================================================
    !> TH distance
!    d_TH = INFINITY
!    if ( th_on ) call TH_DISTANCE(p,j_xyz,d_TH)
    
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

    !> Track-length estimator
    !$omp atomic 
    k_tl = k_tl + distance*p%wgt*macro_xs(4) 

    ! ==================== EX k-eff tally ====================
    ! !$omp atomic 
    ! fiss_vrc = fiss_vrc + p%wgt*macro_xs(4)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
    ! !$omp atomic
    ! loss_vrc = loss_vrc + p%wgt*macro_xs(2)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
    ! =========================================================
    
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

    !> Cycle-power Tally ====================================================================
    ! *** no need ***
    if(curr_cyc > n_inact) then 
        !$omp atomic
        cyc_power = cyc_power + distance*p%wgt*macro_xs(5)
    endif 

    !> Burn-up Tally ========================================================================
    call tally_burnup (p%material, distance, p%wgt, p%E)
    
    !> FMFD Tally (track length) 
    if ( fmfdon .and. inside_FMFD ) call FMFD_TRK(p%wgt,distance,macro_xs,i_xyz)

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
        if ( th_on ) then
            call TH_INSIDE(p%coord(1)%xyz(:),j_xyz(:),inside_th)
            if ( inside_th ) call TH_COL(p%wgt,macro_xs(1),macro_xs(4),j_xyz(:))
        end if

    elseif  ( distance == d_FMFD ) then 
        call FMFD_SURF(inside_FMFD, income_FMFD,i_surf, i_xyz, &
                        p%coord(1)%uvw, p%wgt, surfaces(surface_crossed)%bc)
        if ( fm_crossed ) then
            call cross_surface(p, surface_crossed)
        else
            p%coord(1)%xyz = p%coord(1)%xyz + TINY_BIT * p%coord(1)%uvw
        end if
    else
        p%n_cross = p%n_cross + 1 
        call cross_surface(p, surface_crossed)
        if(p%alive == .false.) then 
            !$omp atomic
            loss_vrc = loss_vrc + p%wgt! * exp(-macro_xs(1)*d_boundary)
        endif 
    endif
    
end subroutine transport

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


!===============================================================================
! transport_VRC handles the Volumetric-Ray-Casting Method 
! in combination with the Delta-tracking 
!===============================================================================
subroutine transport_VRC(p)

    type(Particle), intent(inout) :: p
    type(Particle) :: p_psudo
    
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: d_CMFD                  ! distance to CMFD grid
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: macro_xs(5)
    real(8) :: xyz(3)
    integer :: i_cell, i_bin, i_lat, i_surf
    integer :: i_xyz(3), idx_xyz
    logical :: inside_CMFD
    integer :: bc
    real(8) :: val
    
    
    found_cell = .false.
    i_xyz(:) = -1; idx_xyz = -1
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
    
    call distance_to_boundary(p, d_boundary, surface_crossed)
    
    ! Sample a distance to collision
    if (E_mode == 0) then 
        d_collision = -log(rang())/(sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
        macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
        macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
        
    elseif (E_mode == 1) then 
        macro_xs = getMacroXS(materials(p%material), p%E)
        d_collision = -log(rang())/macro_xs(1)
    endif 
            
    distance = min(d_boundary, d_collision)
    !> Track-length estimator
    !$omp atomic 
    k_tl = k_tl + distance*p%wgt*macro_xs(4)
    
    
    if (p % vrc_traced == .false.) then 
        p_psudo = p
        call trace_psudoray(p_psudo)
        call p_psudo%clear()
        p % vrc_traced = .true.
    endif
    
    
    !> Advance particle
    do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
    enddo
            
    if (distance == d_collision) then ! collision     
        if (E_mode == 0) then 
            !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
            call collision_MG(p)
        else !(E_mode == 1) 
            !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
            call collision_CE(p)
        endif
        p_psudo = p
        call trace_psudoray(p_psudo)
        !call p_psudo%clear()
        
    else
        call cross_surface(p, surface_crossed)
        if(p%alive == .false.) then 
            !$omp atomic
            loss_vrc = loss_vrc + p%wgt! * exp(-macro_xs(1)*d_boundary)
        endif 
    endif
    
end subroutine transport_VRC


!subroutine transport_VRC(p) 
!
!    type(Particle), intent(inout) :: p
!    integer :: i 
!    integer :: j                      ! coordinate level
!    integer :: next_level             ! next coordinate level to check
!    integer :: surface_crossed        ! surface which particle is on
!    integer :: lattice_translation(3) ! in-lattice translation vector
!    real(8) :: d_boundary             ! distance to nearest boundary
!    real(8) :: d_collision            ! distance to collision
!    real(8) :: d_CMFD                  ! distance to CMFD grid
!    real(8) :: distance               ! distance particle travels
!    logical :: found_cell             ! found cell which particle is in?
!    
!    real(8) :: macro_major              ! the global majorant cross-section
!    real(8), allocatable :: macro_tot(:)
!    real(8) :: macro_xs(5)
!    real(8) :: xyz(3)
!    integer :: i_cell, idx_surf
!    integer :: i_xyz(3), idx_xyz, idx
!    integer :: bc
!    integer :: n_mat 
!    type(Particle) :: p_psudo
!    
!    found_cell = .false.
!    i_xyz(:) = -1; idx_xyz = -1
!    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
!    
!    !> Determine macro_major 
!    if (E_mode == 0) then 
!        n_mat = size( XS_MG ) 
!        allocate(macro_tot(n_mat))
!        do i = 1, n_mat 
!            macro_tot(i) = (sum(XS_MG(i)%sig_scat(p%g,:)) + XS_MG(i)%sig_abs(p%g))
!        enddo 
!    else 
!        n_mat = size( materials ) 
!        allocate(macro_tot(n_mat))
!        do i = 1, n_mat 
!            macro_xs = getMacroXS(materials(i), p%E)
!            macro_tot(i) = macro_xs(1) 
!        enddo 
!    endif 
!    macro_major = maxval(macro_tot, n_mat)
!    deallocate(macro_tot)
!    
!    !> Sample a distance to collision
!    d_collision = -log(rang())/macro_major
!    
!    !> Sample distances from special boundaries in univ 0
!    idx = p % coord(1) % cell
!    call cell_distance(cells(idx), p%coord(1)%xyz, p%coord(1)%uvw, surfaces, d_boundary, idx_surf)
!    
!    distance = min(d_boundary, d_collision)
!	
!	if (p % vrc_traced == .false.) then 
!		p_psudo = p
!		call trace_psudoray(p_psudo)
!		call p_psudo%clear()
!		p % vrc_traced = .true.
!	endif
!	
!	print *, p%coord(1)%xyz(1)
!    !> Determine Virtual / Real collsion OR Cross-surface
!    if (d_collision < d_boundary) then ! collision 
!        ! Advance particle
!        p % n_coord = 1
!        p % coord(1) % xyz = p % coord(1) % xyz + distance * p % coord(1) % uvw
!        call find_cell(p, found_cell, i_cell)
!        
!        if (E_mode == 0) then 
!            macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
!        elseif (E_mode == 1) then 
!            macro_xs = getMacroXS(materials(p%material), p%E)
!        endif 
!        
!        ! Reject? 
!        if (rang() < macro_xs(1)/macro_major ) then !> real collision
!			print *, 'real collision'
!
!            if (E_mode == 0) then 
!                !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
!                call collision_MG_DT(p, macro_major)
!            else 
!                !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
!                call collision_CE(p)
!            endif
!			p_psudo = p
!			call trace_psudoray(p_psudo)
!			call p_psudo%clear()
!        endif 
!    else 
!        p % n_coord = 1
!        p % coord(1) % xyz = p % coord(1) % xyz + distance * p % coord(1) % uvw
!        call cross_surface(p, idx_surf)
!		print *, 'real surf cross', p%coord(1)%xyz(1), idx_surf
!    endif
!    
!    p%n_coord = 1 
!    

end module tracking
