module FMFD
    use omp_lib 
    use mpi 
    use geometry_header,    only : lattices, find_lat_idx, lattice_coord
    use variables
    use constants,          only : INFINITY, TiNY_BIT, prt_keff
    use particle_header,    only : particle
    
    implicit none

    type :: FMFD_parameters
        real(8) :: phi
        real(8) :: sig_t 
        real(8) :: sig_a 
        real(8) :: nusig_f 
        real(8) :: J(6) 
        real(8) :: J_pn(6)
        real(8) :: J_pp(6)
    end type
    logical :: fmfdon = .false.
    integer :: CMFD_lat = -1
    integer :: n_skip, n_acc
    integer :: CMFD_type
    real(8) :: CMFD_area(6), CMFD_volume
    integer :: idx_lat
    real(8) :: pitch(3)
    
    type :: FMFD_accumulation
        type(FMFD_parameters), allocatable :: par(:)
    endtype 
    
    type(FMFD_parameters), allocatable :: fm(:), fmavg(:), fmthread(:)
    !$OMP THREADPRIVATE(fmthread)
    type(FMFD_accumulation),target, allocatable :: fmacc(:)

    ! FMFD grid
    real(8):: fm0(3)    ! (x0,y0,z0)
    real(8):: fm1(3)    ! (x1,y1,z1)
    real(8):: fm2(3)    ! (x1-x0,y1-y0,z1-z0)
    integer:: nfm(3)    ! (nx,ny,nz)
    real(8):: dfm(3)    ! (dx,dy,dz)
    
    contains 

    subroutine CMFD_distance (p,i_xyz,idx_xyz,d_CMFD,inside_CMFD,i_surf)
        type(particle), intent(in) :: p
        integer, intent(inout) :: i_xyz(3), idx_xyz
        real(8), intent(inout) :: d_CMFD
        logical, intent(inout) :: inside_CMFD 
        integer, intent(inout) :: i_surf
        real(8) :: xyz(3), uvw(3)
        real(8) :: d_temp(6)
        integer :: i,j, i_coord
        integer :: idx_temp, idx_surf
        real(8) :: J_temp
        integer :: a, b, c
        
        d_CMFD = INFINITY
        if (CMFD_lat < 0) return
        
        ! Find lattice index from the particle
        do i = 1, p%n_coord 
            if (p%coord(i)%lattice == idx_lat) then 
                i_xyz(1) = p%coord(i)%lattice_x
                i_xyz(2) = p%coord(i)%lattice_y
                i_xyz(3) = p%coord(i)%lattice_z
                
                xyz(:) = p%coord(i)%xyz(:)
                uvw(:) = p%coord(i)%uvw(:)
                i_coord = i
                exit
            endif
        enddo
        
        !> Check if the particle is outside the CMFD grid
        inside_CMFD = .true. 
        do i = 1, 3
            if (abs(xyz(i)) > pitch(i)/2.0d0) then 
                inside_CMFD = .false.
                exit 
            endif 
        enddo 
        d_CMFD = INFINITY
        if (inside_CMFD) then
            d_temp(1) = (-pitch(1)/2.0 - xyz(1))/uvw(1)
            d_temp(6) = ( pitch(1)/2.0 - xyz(1))/uvw(1)
            d_temp(2) = (-pitch(2)/2.0 - xyz(2))/uvw(2)
            d_temp(5) = ( pitch(2)/2.0 - xyz(2))/uvw(2)
            d_temp(3) = (-pitch(3)/2.0 - xyz(3))/uvw(3)
            d_temp(4) = ( pitch(3)/2.0 - xyz(3))/uvw(3)
            
             
            do i = 1, 6 
                if (d_temp(i) < 0) cycle
                if (d_CMFD > d_temp(i)) then 
                    d_CMFD = d_temp(i)
                    i_surf = i
                endif 
            enddo 
        
        else !> the particle is outside the CMFD grid
            ! 여기서 xyz는 윗단계 universe의 xyz
            call distance_cuboid(lattices(idx_lat)%n_xyz, &
                        p%coord(i_coord-1)%xyz, p%coord(i_coord-1)%uvw, i_xyz, d_CMFD,i_surf)
            
        endif
        a = lattices(idx_lat)%n_xyz(1)
        b = lattices(idx_lat)%n_xyz(2)
        c = lattices(idx_lat)%n_xyz(3)
        
        idx_xyz = a*b*(i_xyz(3)-1) + a*(i_xyz(2)-1) + i_xyz(1)

    end subroutine
    
! =========================================================================
! FMFD_DISTANCE by Inny
! =========================================================================
subroutine FMFD_distance (p, i_xyz, idx_xyz, d_CMFD, inside_CMFD,i_surf)
    type(particle), intent(in) :: p
    integer, intent(inout) :: i_xyz(3), idx_xyz
    real(8), intent(inout) :: d_CMFD
    logical, intent(inout) :: inside_CMFD 
    integer, intent(inout) :: i_surf
    real(8) :: xyz(3), uvw(3)
    real(8) :: d_temp(6)
    integer :: i,j, i_coord
    integer :: idx_temp, idx_surf
    real(8) :: J_temp
    integer :: a, b, c
    
    if (CMFD_lat < 0) return

    ! Find lattice index in FMFD grid
    i_xyz = FMFD_ID(p%coord(1)%xyz(:))
    
    ! check if the particle is inside the FMFD grid
    inside_CMFD = .true.
    do i = 1, 3
        if ( i_xyz(i) < 0 .or. i_xyz(i) > nfm(i) ) then
            inside_CMFD = .false.
            exit
        end if
    end do

    ! find the closest distance to the gird surface
    d_CMFD = INFINITY
    if ( inside_CMFD ) then
        uvw(:) = p%coord(1)%uvw(:)
        d_temp(1) = ((dfm(1)*(i_xyz(1)-1))-xyz(1))/uvw(1)   ! x0
        d_temp(2) = ((dfm(1)*(i_xyz(1)  ))-xyz(1))/uvw(1)   ! x1
        d_temp(3) = ((dfm(2)*(i_xyz(2)-1))-xyz(2))/uvw(2)   ! y0
        d_temp(4) = ((dfm(2)*(i_xyz(2)  ))-xyz(2))/uvw(2)   ! y1
        d_temp(5) = ((dfm(3)*(i_xyz(3)-1))-xyz(3))/uvw(3)   ! y0
        d_temp(6) = ((dfm(3)*(i_xyz(3)  ))-xyz(3))/uvw(3)   ! y1

        do i = 1, 6
        if ( d_temp(i) > 0 .and. d_CMFD > d_temp(i) ) then
            d_CMFD = d_temp(i)
            i_surf = i
        end if
        end do
    end if

    idx_xyz = nfm(1)*nfm(2)*(i_xyz(3)-1) + nfm(1)*(i_xyz(2)-1) + i_xyz(1)

!    ! Find lattice index from the particle
!    do i = 1, p%n_coord 
!        if (p%coord(i)%lattice == idx_lat) then 
!            i_xyz(1) = p%coord(i)%lattice_x
!            i_xyz(2) = p%coord(i)%lattice_y
!            i_xyz(3) = p%coord(i)%lattice_z
!            
!            xyz(:) = p%coord(i)%xyz(:)
!            uvw(:) = p%coord(i)%uvw(:)
!            i_coord = i
!            exit
!        endif
!    enddo

!    !> Check if the particle is outside the CMFD grid
!    inside_CMFD = .true. 
!    do i = 1, 3
!        if (abs(xyz(i)) > pitch(i)/2.0d0) then 
!            inside_CMFD = .false.
!            exit 
!        endif 
!    enddo 
!    d_CMFD = INFINITY
!    if (inside_CMFD) then
!        d_temp(1) = (-pitch(1)/2.0 - xyz(1))/uvw(1)
!        d_temp(6) = ( pitch(1)/2.0 - xyz(1))/uvw(1)
!        d_temp(2) = (-pitch(2)/2.0 - xyz(2))/uvw(2)
!        d_temp(5) = ( pitch(2)/2.0 - xyz(2))/uvw(2)
!        d_temp(3) = (-pitch(3)/2.0 - xyz(3))/uvw(3)
!        d_temp(4) = ( pitch(3)/2.0 - xyz(3))/uvw(3)
!        
!         
!        do i = 1, 6 
!            if (d_temp(i) < 0) cycle
!            if (d_CMFD > d_temp(i)) then 
!                d_CMFD = d_temp(i)
!                i_surf = i
!            endif 
!        enddo 
!    
!    else !> the particle is outside the CMFD grid
!        ! 여기서 xyz는 윗단계 universe의 xyz
!        call distance_cuboid(lattices(idx_lat)%n_xyz, &
!                    p%coord(i_coord-1)%xyz, p%coord(i_coord-1)%uvw, i_xyz, d_CMFD,i_surf)
!        
!    endif
!    a = lattices(idx_lat)%n_xyz(1)
!    b = lattices(idx_lat)%n_xyz(2)
!    c = lattices(idx_lat)%n_xyz(3)
!    
!    idx_xyz = a*b*(i_xyz(3)-1) + a*(i_xyz(2)-1) + i_xyz(1)

end subroutine

! =============================================================================
! FMFD_ID finds the x, y, z indice in the FMFD mesh grid
! =============================================================================
function FMFD_ID(fmxyz) result(fmid)
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: fmid(3)               ! indice
       
    fmid(:) = int((fmxyz(:)-fm0(:))/dfm(:))+1

end function
    
    
    
    subroutine distance_cuboid(n_xyz, xyz, uvw, i_xyz, d_CMFD,i_surf)
        real(8), intent(in) :: xyz(3), uvw(3)
        integer, intent(in) :: n_xyz(3)
        integer, intent(inout) :: i_xyz(3)
        real(8), intent(inout) :: d_CMFD
        integer, intent(inout) :: i_surf
        real(8) :: d(6), xyz_(3), temp, r(3), xyz_next(3) 
        real(8) :: dist
        integer :: i,j
        
        i_surf = 0 
        xyz_(:) = xyz(:) - lattices(idx_lat)%xyz(:)
        
        
        do i = 1, 3 
            r(i) = n_xyz(i)*pitch(i) / 2.0d0
        enddo
                
        d(1)   = (-r(1)-xyz_(1))/uvw(1)
        temp = xyz_(2)+d(1)*uvw(2)
        if ((temp < -r(2)).or.(temp > r(2))) d(1) = INFINITY
        temp = xyz_(3)+d(1)*uvw(3)
        if ((temp < -r(3)).or.(temp > r(3))) d(1) = INFINITY
        
        d(6)   = ( r(1)-xyz_(1))/uvw(1)
        temp = xyz_(2)+d(6)*uvw(2)
        if ((temp < -r(2)).or.(temp > r(2))) d(6) = INFINITY
        temp = xyz_(3)+d(6)*uvw(3)
        if ((temp < -r(3)).or.(temp > r(3))) d(6) = INFINITY
        
        d(2)   = (-r(2)-xyz_(2))/uvw(2)
        temp = xyz_(1)+d(2)*uvw(1)
        if ((temp < -r(1)).or.(temp > r(1))) d(2) = INFINITY
        temp = xyz_(3)+d(2)*uvw(3)
        if ((temp < -r(3)).or.(temp > r(3))) d(2) = INFINITY
        
        d(5)   = ( r(2)-xyz_(2))/uvw(2)
        temp = xyz_(1)+d(5)*uvw(1)
        if ((temp < -r(1)).or.(temp > r(1))) d(5) = INFINITY
        temp = xyz_(3)+d(5)*uvw(3)
        if ((temp < -r(3)).or.(temp > r(3))) d(5) = INFINITY
        
        d(3)   = (-r(3)-xyz_(3))/uvw(3)
        temp = xyz_(1)+d(3)*uvw(1)
        if ((temp < -r(1)).or.(temp > r(1))) d(3) = INFINITY
        temp = xyz_(2)+d(3)*uvw(2)
        if ((temp < -r(2)).or.(temp > r(2))) d(3) = INFINITY
        
        d(4)   = ( r(3)-xyz_(3))/uvw(3)
        temp = xyz_(1)+d(4)*uvw(1)
        if ((temp < -r(1)).or.(temp > r(1))) d(4) = INFINITY
        temp = xyz_(2)+d(4)*uvw(2)
        if ((temp < -r(2)).or.(temp > r(2))) d(4) = INFINITY
        
        
        d_CMFD = INFINITY
        do i = 1, 6 
            if (d(i) < 0) cycle
            if (d_CMFD > d(i)) then 
                d_CMFD = d(i)
                i_surf = i
            endif 
        enddo 
        
        if (i_surf == 0) return 
        
        xyz_next(:) = xyz_(:) + (d_CMFD + 10*tiny_bit) * uvw(:)
        
        i_xyz = lattice_coord (lattices(idx_lat), xyz_next)
        
        
        
    end subroutine
    
    ! ========================================================= !
    !  CMFD_Tally Cell-wise    :: flux, Sig_t, Sig_a, nuSig_f
    !                Surface-wise :: (net)current
    ! ========================================================= !
    subroutine CMFD_tally(wgt,distance,macro_xs,idx_xyz)
        real(8), intent(in) :: wgt
        real(8), intent(in) :: distance
        real(8), intent(in) :: macro_xs(5)
        integer, intent(in) :: idx_xyz
        real(8) :: flux
        
        flux = wgt * distance
        
        fmthread(idx_xyz) % phi = & 
        fmthread(idx_xyz) % phi + flux
        fmthread(idx_xyz) % sig_t = &
        fmthread(idx_xyz) % sig_t + flux*macro_xs(1)
        fmthread(idx_xyz) % sig_a = &
        fmthread(idx_xyz) % sig_a + flux*macro_xs(2)
        fmthread(idx_xyz) % nusig_f = &
        fmthread(idx_xyz) % nusig_f + flux*macro_xs(4)
        
    end subroutine
    
    subroutine CMFD_tally_col(wgt, macro_xs, idx_xyz, inside_CMFD)
        real(8), intent(in) :: wgt
        real(8), intent(in) :: macro_xs(5)
        integer, intent(in) :: idx_xyz
        logical, intent(in) :: inside_CMFD
        real(8) :: flux
        
        if (.not. inside_CMFD) return 
        
        flux = wgt / macro_xs(1)
        
        fmthread(idx_xyz) % phi = & 
        fmthread(idx_xyz) % phi + flux
        fmthread(idx_xyz) % sig_t = &
        fmthread(idx_xyz) % sig_t + flux*macro_xs(1)
        fmthread(idx_xyz) % sig_a = &
        fmthread(idx_xyz) % sig_a + flux*macro_xs(2)
        fmthread(idx_xyz) % nusig_f = &
        fmthread(idx_xyz) % nusig_f + flux*macro_xs(4)
        
    end subroutine
    
    
    
    subroutine CMFD_curr_tally (inside_CMFD, i_surf, idx_xyz, wgt, uvw, bc) 
        logical, intent(in) :: inside_CMFD
        integer, intent(in) :: i_surf, idx_xyz
        real(8), intent(in) :: wgt, uvw(3)
        integer, intent(in) :: bc
        
        if (inside_CMFD == .true. ) then 
            if (i_surf <= 3) then 
                fmthread(idx_xyz) % J_pn(i_surf) = &
                fmthread(idx_xyz) % J_pn(i_surf) + wgt / CMFD_area(i_surf)
                if (bc == 2) then 
                    fmthread(idx_xyz) % J_pp(i_surf) = &
                    fmthread(idx_xyz) % J_pp(i_surf) + wgt / CMFD_area(i_surf)
                endif
            else 
                fmthread(idx_xyz) % J_pp(i_surf) = &
                fmthread(idx_xyz) % J_pp(i_surf) + wgt / CMFD_area(i_surf)
                if (bc == 2) then 
                    fmthread(idx_xyz) % J_pn(i_surf) = &
                    fmthread(idx_xyz) % J_pn(i_surf) + wgt / CMFD_area(i_surf)
                endif
            endif
        else 
            if (i_surf > 3) then 
                fmthread(idx_xyz) % J_pn(i_surf) = &
                fmthread(idx_xyz) % J_pn(i_surf) + wgt / CMFD_area(i_surf)
            else 
                fmthread(idx_xyz) % J_pp(i_surf) = &
                fmthread(idx_xyz) % J_pp(i_surf) + wgt / CMFD_area(i_surf)
            endif
        endif 
                
    end subroutine
    
    ! ========================================================= !
    !  FMFD_initialize :: initialize CMFD tally bins
    ! ========================================================= !
    subroutine FMFD_initialize()
        integer :: i 
        real(8) :: pitch_temp(6)
        type(FMFD_accumulation), pointer :: acc
        
        allocate(fm(nfm(1)*nfm(2)*nfm(3)))
        allocate(fmavg(nfm(1)*nfm(2)*nfm(3)))
                
        allocate(fmacc(n_acc)) 
        do i = 1, n_acc 
            acc => fmacc(i)
            allocate(acc%par(nfm(1)*nfm(2)*nfm(3)))
        enddo
        
        do i = 1, 3 
            if (lattices(idx_lat)%n_xyz(i) == 1) then 
                pitch(i) = INFINITY
                pitch_temp(i) = 1 
            else 
                pitch(i) = lattices(idx_lat)%pitch(i)
                pitch_temp(i) = lattices(idx_lat)%pitch(i)
            endif
        enddo 
        
        CMFD_area(1) = pitch_temp(2)*pitch_temp(3)
        CMFD_area(6) = pitch_temp(2)*pitch_temp(3)
        CMFD_area(2) = pitch_temp(1)*pitch_temp(3)
        CMFD_area(5) = pitch_temp(1)*pitch_temp(3)
        CMFD_area(3:4) = pitch_temp(1)*pitch_temp(2)
        CMFD_volume = pitch_temp(1)*pitch_temp(2)*pitch_temp(3)

        fm(:) % phi     = 0 
        fm(:) % sig_t   = 0 
        fm(:) % sig_a   = 0 
        fm(:) % nusig_f = 0 
        do i = 1, 6
            fm(:) % J(i)    = 0 
            fm(:) % J_pn(i) = 0 
            fm(:) % J_pp(i) = 0 
        enddo 
            
    end subroutine
    
    subroutine FMFD_initialize_thread() 
        integer :: i 
        
        allocate(fmthread(nfm(1)*nfm(2)*nfm(3)))
        
        fmthread(:) % phi     = 0 
        fmthread(:) % sig_t   = 0 
        fmthread(:) % sig_a   = 0 
        fmthread(:) % nusig_f = 0 
        do i = 1, 6
        fmthread(:) % J(i)    = 0 
        fmthread(:) % J_pn(i) = 0 
        fmthread(:) % J_pp(i) = 0 
        enddo 
        
    end subroutine
    
    ! ========================================================= !
    !  CMFD_solve 
    ! ========================================================= !
    subroutine FMFD_solve(a,b,c,shape)
        integer, intent(in) :: a,b,c 
        real(8), intent(inout) :: shape(:) 
        
        integer :: i, j, k, idx, size
        
        real(8) :: sig_a(a*b*c), nusig_f(a*b*c), D(a*b*c)
        real(8) :: F(a*b*c)
        real(8) :: M(a*b*c,7), D_tilda(a*b*c,6), D_hat(a*b*c,6), Js(a*b*c,6), J_pp(a*b*c,6),J_pn(a*b*c,6)
        real(8) :: keff_CMFD, phi(a*b*c), phi0(a*b*c), phi_CMFD(a*b*c), phi_noacc(a*b*c)
                
        do i = 1, a*b*c
            sig_a(i)   = fmavg(i)%sig_a
            nusig_f(i) = fmavg(i)%nusig_f 
            D(i)       = 1. / (3. * fmavg(i)%sig_t) 
            Js(i,:)    = fmavg(i)%J(:) 
            J_pn(i,:)  = fmavg(i)%J_pn(:) 
            J_pp(i,:)  = fmavg(i)%J_pp(:) 
            phi(i)     = fmavg(i)%phi
            
            phi_noacc(i) = fm(i)%phi
        enddo 
        
        
        size = a*b*c
        !> calculate D_tilda & D_hat 
        call D_hat_calculation(D_hat, Js, phi,D, a,b,c,size, pitch(1),pitch(2),pitch(3))
         
         
        if (CMFD_type == 1) then 
            !> CMFD 
            M = getM(D, D_hat, sig_a, a,b,c,pitch(1),pitch(2),pitch(3))
        else!if (CMFD_type == 2) then 
            !> p-CMFD
            M = getMp (phi, D, J_pp, J_pn, sig_a, a,b,c,pitch(1),pitch(2),pitch(3))
        endif 
        F = nusig_f
        
        keff_CMFD = keff
        phi_CMFD  = phi 
        call powerIter (keff_CMFD, M, F, phi_CMFD, nusig_f, size, a,b)
        
        
        !print *, 'CMFD keff  ',keff_CMFD
        !write(prt_keff,*) keff_CMFD, k_col, k_tl
        
        !> CMFD feedback (modulation)
        phi_noacc(:)= phi_noacc(:)*size/sum(phi_noacc)
        phi_CMFD(:) = phi_CMFD(:)*size/sum(phi_CMFD)
        shape(:)     = phi_CMFD(:)/phi_noacc(:)
        
    end subroutine
    
    
    
    ! ========================================================= !
    !  Miscel. functions for CMFD_solve
    ! ========================================================= !
    
    ! function which produces M matrix from xs & D_hat
    function getM (D, D_hat, sig_a, a,b,c,dx,dy,dz) result (M)
        integer :: i, j, a,b,c, index, index0
        real(8) :: dx, dy, dz
        real(8), dimension(:), intent(in):: sig_a(*), D(*)
        real(8), dimension(:,:) :: M(a*b*c,7), D_tilda(a*b*c,6), D_hat(a*b*c,6)
        integer, dimension(3) :: xyz
        
        M(:,:) = 0 
        index0 = a*b*c
        
        ! D_tilda set 
        do i = 1, index0
            xyz = getXYZ(i, a,b,c)

            if (xyz(1).ne.1) then 
                D_tilda(i,1) = 2.0*D(i)*D(i-1)/(D(i)+D(i-1))/dx
            else 
                D_tilda(i,1) = 0
            endif
            if (xyz(1).ne.a) then 
                D_tilda(i,6) = 2.0*D(i)*D(i+1)/(D(i)+D(i+1))/dx
            else 
                D_tilda(i,6) = 0
            endif 
            if (xyz(2).ne.1) then 
                D_tilda(i,2) = 2.0*D(i)*D(i-a)/(D(i)+D(i-a))/dy
            else 
                D_tilda(i,2) = 0
            endif
            if (xyz(2).ne.b) then 
                D_tilda(i,5) = 2.0*D(i)*D(i+a)/(D(i)+D(i+a))/dy
            else 
                D_tilda(i,5) = 0
            endif 
            if (xyz(3).ne.1) then 
                D_tilda(i,3) = 2*D(i)*D(i-b*a)/(D(i)+D(i-b*a))/dz
            else 
                D_tilda(i,3) = 0
            endif 
            if (xyz(3).ne.c) then 
                D_tilda(i,4) = 2.0*D(i)*D(i+b*a)/(D(i)+D(i+b*a))/dz
            else 
                D_tilda(i,4) = 0
            endif
            
        enddo 
        
        ! M matrix set 
        do i = 1, index0
            M(i,1) = -(D_tilda(i,1)+D_hat(i,1))/dx
            M(i,6) = (-D_tilda(i,6)+D_hat(i,6))/dx
            M(i,2) = -(D_tilda(i,2)+D_hat(i,2))/dy
            M(i,5) = (-D_tilda(i,5)+D_hat(i,5))/dy
            M(i,3) = -(D_tilda(i,3)+D_hat(i,3))/dz
            M(i,4) = (-D_tilda(i,4)+D_hat(i,4))/dz
            
            M(i,7) = (-(-D_tilda(i,1)+D_hat(i,1))+(D_tilda(i,6)+D_hat(i,6)))/dx &
                    +(-(-D_tilda(i,2)+D_hat(i,2))+(D_tilda(i,5)+D_hat(i,5)))/dy &
                    +(-(-D_tilda(i,3)+D_hat(i,3))+(D_tilda(i,4)+D_hat(i,4)))/dz + sig_a(i)
                        
            xyz = getXYZ(i, a,b,c)
            
            if (xyz(1) == 1) M(i,1) = 0
            if (xyz(1) == a) M(i,6) = 0
            if (xyz(2) == 1) M(i,2) = 0
            if (xyz(2) == b) M(i,5) = 0
            if (xyz(3) == 1) M(i,3) = 0
            if (xyz(3) == c) M(i,4) = 0
        enddo
        
    end function getM
    
    
    subroutine powerIter (k_eff, M, F, phi,nusigf, size, a,b)
        integer :: i, j, size, a, b, iter_max, iter
        real(8) :: k_eff, k_eff0
        real(8) :: suma, sumb, err
        
        real(8), dimension(size) :: F, phi, phi0, Q, nusigf
        real(8), dimension(size,7) :: M
        
        k_eff = 1 
        err = 1; iter_max = 1000; iter = 0;
        do while ((err.gt.1.0d-10).and.(iter.lt.iter_max))
            iter = iter + 1
            !print *, iter, k_eff, err
            k_eff0 = k_eff; phi0 = phi;
            Q(:) = F(:)*phi(:)/k_eff
            phi(:) = 0
            !call conjgrad_hepta(M,Q,phi,size,a,b)
            phi = BiCGStab_hepta(M,Q,a,b)
            suma = 0; sumb = 0; 
            do i = 1, size
                suma = suma + F(i)*phi(i)*F(i)*phi(i)
                sumb = sumb + F(i)*phi(i)*F(i)*phi0(i)
            enddo
            k_eff = k_eff*suma/sumb
            
            err  = error_calculation(phi, phi0, size)
        enddo ! while loop
        
        
    end subroutine 
    
    function error_calculation(phi1, phi10, dim) result (error_max)

        integer :: dim, i
        real(8) :: phi1(dim), phi10(dim)
        real(8) :: error, error_Max

        
        error=abs(phi1(1)-phi10(1))/phi1(1)
        error_Max=error
        do i=2, dim
            error=abs(phi1(i)-phi10(i))/phi1(i)
            if(error.GT.error_Max) then
                error_Max=error
            end if
        end do
    
    end function error_calculation
    
    function getXYZ(index, a, b, c) result(xyz)
        integer, intent(in) :: index, a,b,c
        integer, dimension(3) :: xyz

        if ((index.le.0).or.(index.gt.(a*b*c))) then
            print *, 'function getXYZ() :: INDEX OUT OF RANGE'
            stop
        endif
        
        xyz(3) = ceiling(real(index, 8)/real(a*b,8))
        xyz(2) = ceiling(real(index - (xyz(3)-1)*a*b, 8)/a)
        xyz(1) = index - (xyz(3)-1)*a*b - (xyz(2)-1)*a
    
    end function getXYZ

    
    subroutine D_hat_calculation(Dhat, Js, phi,D, a,b,c,size, dx,dy,dz)
        integer :: i, j, jj, size, a,b,c
        integer, dimension(3) :: xyz

        real(8) :: D_tilda, dx, dy, dz
        real(8), intent(inout):: Dhat(size,6)
        real(8), intent(in):: phi(size), D(size), Js(size,6)
        
        xyz(:) = 0 
        do j = 1, size
            xyz = getXYZ(j, a,b,c)
            do jj = 1, 6
                if (jj.eq.1) then
                    if ( xyz(1) == 1 ) then 
                        Dhat(j,jj) = Js(j,jj)/phi(j) 
                    else
                       D_tilda    = 2*D(j)*D(j-1)/(D(j)+D(j-1))/dx
                       Dhat(j,jj) = ((Js(j,jj) + D_tilda*(phi(j)-phi(j-1)))/(phi(j)+phi(j-1)))
                    endif 

                elseif (jj.eq.6) then
                    if (xyz(1) == a) then
                        Dhat(j,jj) = Js(j,jj)/phi(j)
                    else
                       D_tilda    = 2*D(j)*D(j+1)/(D(j)+D(j+1))/dx
                       Dhat(j,jj) = ((Js(j,jj) + D_tilda*(phi(j+1)-phi(j)))/(phi(j+1)+phi(j)))
                    endif

                elseif (jj == 2) then
                    if (xyz(2) == 1) then
                       Dhat(j,jj) = Js(j,jj)/phi(j) 
                    else
                       D_tilda    = 2*D(j)*D(j-a)/(D(j)+D(j-a))/dy
                       Dhat(j,jj) = ((Js(j,jj) + D_tilda*(phi(j)-phi(j-a)))/(phi(j)+phi(j-a)))
                    endif

                elseif (jj == 5) then
                    if (xyz(2) == b) then
                        Dhat(j,jj) = Js(j,jj)/phi(j) 
                    else
                       D_tilda    = 2*D(j)*D(j+a)/(D(j)+D(j+a))/dy
                       Dhat(j,jj) = ((Js(j,jj) + D_tilda*(phi(j+a)-phi(j)))/(phi(j+a)+phi(j)))
                    endif

                elseif (jj == 3) then
                    if (xyz(3) == 1) then
                        Dhat(j,jj) = Js(j,jj)/phi(j) 
                    else 
                       D_tilda    = 2*D(j)*D(j-b*a)/(D(j)+D(j-b*a))/dz
                       Dhat(j,jj) = ((Js(j,jj) + D_tilda*(phi(j)-phi(j-a*b)))/(phi(j)+phi(j-a*b)))
                    endif 

                elseif (jj == 4) then
                    if (xyz(3) == c) then
                        Dhat(j,jj) = Js(j,jj)/phi(j) 
                    else
                       D_tilda    = 2*D(j)*D(j+b*a)/(D(j)+D(j+b*a))/dz
                       Dhat(j,jj) = ((Js(j,jj) + D_tilda*(phi(j+a*b)-phi(j)))/(phi(j+a*b)+phi(j)))
                    endif 

                endif 

            enddo 
        enddo 

    end subroutine     
    
    
    
    ! function which produces p-CMFD M matrix
    function getMp (phi, D, J_pp, J_pn, sig_a, a,b,c,dx,dy,dz) result (M)
        integer, intent(in) :: a,b,c
        real(8), intent(in) :: phi(a*b*c), sig_a(a*b*c), D(a*b*c), J_pp(a*b*c,6), J_pn(a*b*c,6), dx, dy, dz
        real(8) :: M(a*b*c,7)
        real(8) :: D_tilda(a*b*c,6), D_hat_p(a*b*c,6), D_hat_n(a*b*c,6)
        integer :: xyz(3)
        integer :: i, j, index, index0
        
        M(:,:) = 0
        index0 = a*b*c
        
        ! D_tilda set 
        do i = 1, index0
            xyz = getXYZ(i,a,b,c)
            M(i,7) = 0 
            if (xyz(1).ne.1) then 
                D_tilda(i,1) = 2.0*D(i)*D(i-1)/(D(i)+D(i-1))/dx
                D_hat_p(i,1) = -(1./2./phi(i-1))*(2*J_pp(i,1) + D_tilda(i,1)*(phi(i)-phi(i-1)))
                D_hat_n(i,1) =  (1./2./phi(i))*(2*J_pn(i,1)   - D_tilda(i,1)*(phi(i)-phi(i-1)))
                M(i,7) = M(i,7) + (D_tilda(i,1)+D_hat_n(i,1))/dx
            else 
                D_tilda(i,1) = 0
                D_hat_p(i,1) = -J_pp(i,1)/phi(i)
                D_hat_n(i,1) =  J_pn(i,1)/phi(i)
                M(i,7) = M(i,7) + (D_hat_p(i,1) + D_hat_n(i,1))/dx
            endif
            if (xyz(1).ne.a) then 
                D_tilda(i,6) = 2.0*D(i)*D(i+1)/(D(i)+D(i+1))/dx
                D_hat_p(i,6) = -(1./2./phi(i))*(2*J_pp(i,6) + D_tilda(i,6)*(phi(i+1)-phi(i)))
                D_hat_n(i,6) =  (1./2./phi(i+1))*(2*J_pn(i,6) - D_tilda(i,6)*(phi(i+1)-phi(i)))
                M(i,7) = M(i,7) + (D_tilda(i,6)-D_hat_p(i,6))/dx
            else 
                D_tilda(i,6) = 0
                D_hat_p(i,6) = -J_pp(i,6)/phi(i)
                D_hat_n(i,6) =  J_pn(i,6)/phi(i)
                M(i,7) = M(i,7) - (D_hat_p(i,6) + D_hat_n(i,6))/dx
            endif 
            
            if (xyz(2).ne.1) then 
                D_tilda(i,2) = 2.0*D(i)*D(i-a)/(D(i)+D(i-a))/dy
                D_hat_p(i,2) = -(1./2./phi(i-a))*(2*J_pp(i,2) + D_tilda(i,2)*(phi(i)-phi(i-a)))
                D_hat_n(i,2) =  (1./2./phi(i))  *(2*J_pn(i,2) - D_tilda(i,2)*(phi(i)-phi(i-a)))
                M(i,7) = M(i,7) + (D_tilda(i,2)+D_hat_n(i,2))/dy 
            else 
                D_tilda(i,2) = 0
                D_hat_p(i,2) = -J_pp(i,2)/phi(i)
                D_hat_n(i,2) =  J_pn(i,2)/phi(i)
                M(i,7) = M(i,7) + (D_hat_p(i,2) + D_hat_n(i,2))/dy
            endif
            if (xyz(2).ne.b) then 
                D_tilda(i,5) = 2.0*D(i)*D(i+a)/(D(i)+D(i+a))/dy
                D_hat_p(i,5) = -(1./2./phi(i))*(2*J_pp(i,5)   + D_tilda(i,5)*(phi(i+a)-phi(i)))
                D_hat_n(i,5) =  (1./2./phi(i+a))*(2*J_pn(i,5) - D_tilda(i,5)*(phi(i+a)-phi(i)))
                M(i,7) = M(i,7) + (D_tilda(i,5)-D_hat_p(i,5))/dy
            else 
                D_tilda(i,5) = 0
                D_hat_p(i,5) = -J_pp(i,5)/phi(i)
                D_hat_n(i,5) =  J_pn(i,5)/phi(i)
                M(i,7) = M(i,7) - (D_hat_p(i,5) + D_hat_n(i,5))/dy
            endif 
            
            if (xyz(3).ne.1) then 
                D_tilda(i,3) = 2*D(i)*D(i-b*a)/(D(i)+D(i-b*a))/dz
                D_hat_p(i,3) = -(1./2./phi(i-a*b))*(2*J_pp(i,3) + D_tilda(i,3)*(phi(i)-phi(i-a*b)))
                D_hat_n(i,3) =  (1./2./phi(i))*(2*J_pn(i,3)     - D_tilda(i,3)*(phi(i)-phi(i-a*b)))
                M(i,7) = M(i,7) + (D_tilda(i,3)+D_hat_n(i,3))/dz
            else 
                D_tilda(i,3) = 0
                D_hat_p(i,3) = -J_pp(i,3)/phi(i)
                D_hat_n(i,3) =  J_pn(i,3)/phi(i)
                M(i,7) = M(i,7) + (D_hat_p(i,3) + D_hat_n(i,3))/dz
            endif 
            if (xyz(3).ne.c) then 
                D_tilda(i,4) = 2.0*D(i)*D(i+b*a)/(D(i)+D(i+b*a))/dz
                D_hat_p(i,4) = -(1./2./phi(i))*(2*J_pp(i,4)     + D_tilda(i,4)*(phi(i+a*b)-phi(i)))
                D_hat_n(i,4) =  (1./2./phi(i+a*b))*(2*J_pn(i,4) - D_tilda(i,4)*(phi(i+a*b)-phi(i)))
                M(i,7) = M(i,7) + (D_tilda(i,4)-D_hat_p(i,4))/dz
            else 
                D_tilda(i,4) = 0
                D_hat_p(i,4) = -J_pp(i,4)/phi(i)
                D_hat_n(i,4) =  J_pn(i,4)/phi(i)
                M(i,7) = M(i,7) - (D_hat_p(i,4) + D_hat_n(i,4))/dz
            endif
        
            ! M matrix set
            M(i,1) = -(D_tilda(i,1)-D_hat_p(i,1))/dx
            M(i,6) = -(D_tilda(i,6)+D_hat_n(i,6))/dx
            M(i,2) = -(D_tilda(i,2)-D_hat_p(i,2))/dy
            M(i,5) = -(D_tilda(i,5)+D_hat_n(i,5))/dy
            M(i,3) = -(D_tilda(i,3)-D_hat_p(i,3))/dz
            M(i,4) = -(D_tilda(i,4)+D_hat_n(i,4))/dz
            
            M(i,7) = M(i,7) + sig_a(i)
            
            xyz = getXYZ(i, a,b,c)
            
            if (xyz(1) == 1) M(i,1) = 0
            if (xyz(1) == a) M(i,6) = 0
            if (xyz(2) == 1) M(i,2) = 0
            if (xyz(2) == b) M(i,5) = 0
            if (xyz(3) == 1) M(i,3) = 0
            if (xyz(3) == c) M(i,4) = 0
            
        enddo
        
    end function getMp
    
    
    ! ============================================== !
    !         Hepta diagonal matrix solvers            !
    ! ============================================== !
    function BiCGStab_hepta(M,Q,a,b) result(x)
        real(8), intent(in)                   :: M (:,:)
        real(8), intent(in)                   :: Q ( : )
        real(8), dimension(1:size(Q, dim=1))  :: x
        real(8), dimension(1:size(Q, dim=1))  :: r, rs, v, p, s, t
        real(8), parameter                    :: e = 1d-10
        real(8)                               :: rho      , rho_prev
        real(8)                               :: alpha    , omega   , beta
        real(8)                               :: norm_r   , norm_b       
        real(8)                               :: summesion, temp
        integer, intent(in)                      :: a, b
        integer                               :: it=0, i, index0

        x   = 0.0
        r   = Q
        rs  = r
        rho = 1.0; alpha = 1.0; omega = 1.0
        v   = 0.0; p  = 0.0

        norm_r = sqrt(dot_product(r,r))
        norm_b = sqrt(dot_product(Q,Q))
        
        index0 = size(Q, dim=1)

        do while(norm_r .GT. e*norm_b)
            rho_prev = rho                                      
            rho      = dot_product(rs,r)                        
            beta     = (rho/rho_prev) * (alpha/omega)           
        
            p        = r + beta * (p - omega*v)                 
        
            v(:) = 0
            do i = 1, index0  !> v = matmul(M,p)
                if (M(i,1).ne.0) v(i) = v(i) + M(i,1)*p(i-1)
                if (M(i,6).ne.0) v(i) = v(i) + M(i,6)*p(i+1)
                if (M(i,2).ne.0) v(i) = v(i) + M(i,2)*p(i-a)
                if (M(i,5).ne.0) v(i) = v(i) + M(i,5)*p(i+a)
                if (M(i,3).ne.0) v(i) = v(i) + M(i,3)*p(i-a*b)
                if (M(i,4).ne.0) v(i) = v(i) + M(i,4)*p(i+a*b)
                                 v(i) = v(i) + M(i,7)*p(i)
            enddo 
            
            alpha    = rho/dot_product(rs,v)
            s        = r - alpha*v
            t(:) = 0
            do i = 1, index0  !> t = matmul(M,s)
                if (M(i,1).ne.0) t(i) = t(i) + M(i,1)*s(i-1)
                if (M(i,6).ne.0) t(i) = t(i) + M(i,6)*s(i+1)
                if (M(i,2).ne.0) t(i) = t(i) + M(i,2)*s(i-a)
                if (M(i,5).ne.0) t(i) = t(i) + M(i,5)*s(i+a)
                if (M(i,3).ne.0) t(i) = t(i) + M(i,3)*s(i-a*b)
                if (M(i,4).ne.0) t(i) = t(i) + M(i,4)*s(i+a*b)
                                 t(i) = t(i) + M(i,7)*s(i)
            enddo 
            
            
            omega    = dot_product(t,s)/dot_product(t,t)
            x        = x + alpha*p + omega*s
            r        = s - omega*t
            norm_r   = sqrt(dot_product(r,r))
            norm_b   = sqrt(dot_product(Q,Q))
        
            it = it + 1
        end do   
        
        return
    end function BiCGStab_hepta     

    subroutine conjgrad_hepta(M,Q,x,index,a, b)  
        integer :: index, i, iter, j, a, b
        real(8), dimension(:,:) :: M(index, 7)
        real(8), dimension(index)   :: Q, x, r, p, Ap
        real(8) :: rsold, alpha, rsnew
        
        r = Q
        p = r
        rsold = dot_product(r,r)
    
        rsnew = 1
        iter = 0 
        do while (sqrt(rsnew).gt.1.0d-10)
            iter = iter+1 
            ! Ap = A*p
            Ap(:) = 0
            do i = 1, index
                if (M(i,1).ne.0) Ap(i) = Ap(i) + M(i,1)*p(i-1)
                if (M(i,6).ne.0) Ap(i) = Ap(i) + M(i,6)*p(i+1)
                if (M(i,2).ne.0) Ap(i) = Ap(i) + M(i,2)*p(i-a)
                if (M(i,5).ne.0) Ap(i) = Ap(i) + M(i,5)*p(i+a)
                if (M(i,3).ne.0) Ap(i) = Ap(i) + M(i,3)*p(i-a*b)
                if (M(i,4).ne.0) Ap(i) = Ap(i) + M(i,4)*p(i+a*b)
                                 Ap(i) = Ap(i) + M(i,7)*p(i)
            enddo 
            
            alpha = rsold/dot_product(p,Ap)
            
            do j = 1, index
                x(j) = x(j) + alpha*p(j)
                r(j) = r(j) - alpha*Ap(j)
            enddo
            
            rsnew = dot_product(r,r)
            
            do j = 1, index
                p(j) = r(j) + (rsnew/rsold)*p(j)
            enddo
            
            rsold = rsnew
        enddo
        
    end subroutine conjgrad_hepta
    
    

    subroutine normalize_FMFD()
        integer:: i
    
        if (CMFD_lat <= 0) return
        
        fmthread(:) % phi     = fmthread(:) % phi     / dble(ngen) / CMFD_volume
        fmthread(:) % sig_t   = fmthread(:) % sig_t   / dble(ngen) / CMFD_volume
        fmthread(:) % sig_a   = fmthread(:) % sig_a   / dble(ngen) / CMFD_volume
        fmthread(:) % nusig_f = fmthread(:) % nusig_f / dble(ngen) / CMFD_volume
        do i = 1, 6
        fmthread(:) % J_pp(i) = fmthread(:) % J_pp(i) / dble(ngen)
        fmthread(:) % J_pn(i) = fmthread(:) % J_pn(i) / dble(ngen)
        end do
        
        !> gather thread CMFD parameters
        fm(:) % phi     = fm(:) % phi     + fmthread(:)%phi
        fm(:) % sig_t   = fm(:) % sig_t   + fmthread(:)%sig_t 
        fm(:) % sig_a   = fm(:) % sig_a   + fmthread(:)%sig_a 
        fm(:) % nusig_f = fm(:) % nusig_f + fmthread(:)%nusig_f 
        do i = 1, 6
        fm(:) % J_pp(i) = fm(:) % J_pp(i) + fmthread(:)%J_pp(i)
        fm(:) % J_pn(i) = fm(:) % J_pn(i) + fmthread(:)%J_pn(i)
        end do

    end subroutine
    
    
    subroutine process_FMFD() 

        integer :: a, b, c, n
        integer :: i, xyz(3), idx, i_surf, i_acc, i_curr
        !> MPI derived type reduce parameters 
        integer, dimension(7) :: CMFD_blocklength, CMFD_displacement, CMFD_datatype 
        integer :: intex, realex, restype, CMFD_op  ! MPI int & real extern / new type
        type(FMFD_parameters), allocatable :: CMFD_MPI_slut(:)
        
        if (CMFD_lat <= 0) return 
        
        data CMFD_blocklength /1, 1, 1, 1, 6, 6, 6 /
        
        a = lattices(idx_lat)%n_xyz(1)
        b = lattices(idx_lat)%n_xyz(2)
        c = lattices(idx_lat)%n_xyz(3)
        allocate(CMFD_MPI_slut(a*b*c))
        n = a*b*c 
        
        !> Gather CMFD parameters from the slave nodes 
        call MPI_TYPE_EXTENT(MPI_REAL8, realex, ierr) 
        CMFD_displacement(1) = 0
        CMFD_displacement(2) = realex
        CMFD_displacement(3) = 2*realex
        CMFD_displacement(4) = 3*realex
        CMFD_displacement(5) = 4*realex
        CMFD_displacement(6) = (4+6)*realex
        CMFD_displacement(7) = (4+2*6)*realex
        
        CMFD_datatype(:) = MPI_REAL8
        
        call MPI_TYPE_STRUCT (7, CMFD_blocklength, CMFD_displacement, CMFD_datatype,restype, ierr)
        call MPI_TYPE_COMMIT (restype, ierr)
        call MPI_Op_create(CMFD_sum, .true. , CMFD_op, ierr)
        call MPI_REDUCE(fm, CMFD_MPI_slut, n, restype, CMFD_op, score, MPI_COMM_WORLD, ierr)
        fm = CMFD_MPI_slut
         
        call MPI_TYPE_FREE(restype, ierr)
        deallocate(CMFD_MPI_slut)
        call MPI_Op_Free(CMFD_Op, ierr)
        
        if (icore /= score) return 
        
        do idx = 1, n 
            xyz = getXYZ(idx, a,b,c) 
            if (xyz(1) /= 1) fm(idx)%J_pp(1) = fm(idx-1)%J_pp(6)
            if (xyz(1) /= a) fm(idx)%J_pn(6) = fm(idx+1)%J_pn(1)
            if (xyz(2) /= 1) fm(idx)%J_pp(2) = fm(idx-a)%J_pp(5)
            if (xyz(2) /= b) fm(idx)%J_pn(5) = fm(idx+a)%J_pn(2)
            if (xyz(3) /= 1) fm(idx)%J_pp(3) = fm(idx-a*b)%J_pp(4)
            if (xyz(3) /= c) fm(idx)%J_pn(4) = fm(idx+a*b)%J_pn(3)
            
            do i_surf = 1, 6
                fm(idx) % J(i_surf) = &
                    (fm(idx) % J_pp(i_surf) - fm(idx) % J_pn(i_surf))
            enddo 
        enddo 
        
        fm(:) % sig_t   = fm(:) % sig_t   / fm(:) % phi
        fm(:) % sig_a   = fm(:) % sig_a   / fm(:) % phi
        fm(:) % nusig_f = fm(:) % nusig_f / fm(:) % phi
        

        !> save the next accumulation
        do i = n_acc, 2, -1
            fmacc(i)%par(:) = fmacc(i-1)%par(:)
        enddo 
        fmacc(1)%par(:) = fm(:)
        
        
        !> average with the accumulated parameters
        fmavg(:)%phi      = 0 
        fmavg(:)%sig_t     = 0 
        fmavg(:)%sig_a     = 0 
        fmavg(:)%nusig_f    = 0 
        do i_curr = 1, 6 
          fmavg(:)%J(i_curr) = 0 
          fmavg(:)%J_pn(i_curr) = 0 
          fmavg(:)%J_pp(i_curr) = 0 
        enddo 
        
        do i = 1, a*b*c
            do i_acc = 1, n_acc
                fmavg(i)%phi     = fmavg(i)%phi     + fmacc(i_acc)%par(i)%phi
                fmavg(i)%sig_t   = fmavg(i)%sig_t   + fmacc(i_acc)%par(i)%sig_t
                fmavg(i)%sig_a   = fmavg(i)%sig_a   + fmacc(i_acc)%par(i)%sig_a
                fmavg(i)%nusig_f = fmavg(i)%nusig_f + fmacc(i_acc)%par(i)%nusig_f
                do i_curr = 1, 6 
                fmavg(i)%J(i_curr)    = &
                    fmavg(i)%J(i_curr)    + fmacc(i_acc)%par(i)%J(i_curr)
                fmavg(i)%J_pn(i_curr) = &
                    fmavg(i)%J_pn(i_curr) + fmacc(i_acc)%par(i)%J_pn(i_curr)
                fmavg(i)%J_pp(i_curr) = &
                    fmavg(i)%J_pp(i_curr) + fmacc(i_acc)%par(i)%J_pp(i_curr)
                enddo 
            enddo 
        enddo
        
        fmavg(:)%phi     = fmavg(:)%phi     / dble(n_acc)
        fmavg(:)%sig_t   = fmavg(:)%sig_t   / dble(n_acc)
        fmavg(:)%sig_a   = fmavg(:)%sig_a   / dble(n_acc)
        fmavg(:)%nusig_f = fmavg(:)%nusig_f / dble(n_acc)
        do i = 1, 6
            fmavg(:)%J(i)    = fmavg(:)%J(i)    / dble(n_acc)
            fmavg(:)%J_pn(i) = fmavg(:)%J_pn(i) / dble(n_acc)
            fmavg(:)%J_pp(i) = fmavg(:)%J_pp(i) / dble(n_acc)
        enddo
        
    end subroutine 
    
    subroutine CMFD_sum (inpar, inoutpar, len, partype) 
        integer :: len, jj, partype
        type(FMFD_parameters), intent(in) :: inpar(len)
        type(FMFD_parameters), intent(inout) :: inoutpar(len) 
        
        inoutpar(:)%phi     = inoutpar(:)%phi     + inpar(:)%phi        
        inoutpar(:)%sig_t   = inoutpar(:)%sig_t   + inpar(:)%sig_t     
        inoutpar(:)%sig_a   = inoutpar(:)%sig_a   + inpar(:)%sig_a     
        inoutpar(:)%nusig_f = inoutpar(:)%nusig_f + inpar(:)%nusig_f
        do jj = 1, 6 
            inoutpar(:)%J(jj)    = inoutpar(:)%J(jj)    + inpar(:)%J(jj)
            inoutpar(:)%J_pn(jj) = inoutpar(:)%J_pn(jj) + inpar(:)%J_pn(jj)
            inoutpar(:)%J_pp(jj) = inoutpar(:)%J_pp(jj) + inpar(:)%J_pp(jj)
        enddo 

    end subroutine  

    
    
end module     
