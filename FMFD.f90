module FMFD
    use omp_lib 
    use mpi 
    use geometry_header,    only : lattices, find_lat_idx, lattice_coord
    use variables
    use constants,          only : INFINITY, TiNY_BIT, prt_keff
    use particle_header,    only : particle
    
    implicit none
    ! x0 / x1 / y0 / y1 / z0 / z1

    type :: FMFD_parameters
        real(8) :: phi
        real(8) :: sig_t 
        real(8) :: sig_a 
        real(8) :: nusig_f 
        real(8) :: Jn(6) 
        real(8) :: J0(6)
        real(8) :: J1(6)
    end type
    logical :: fmfdon = .false.
    integer :: n_skip, n_acc
    integer :: FMFD_type        ! 1 FMFD / 2 p-FMFD / 3 1-node CMFD
    real(8) :: a_fm(6), v_fm
    
    type :: FMFD_accumulation
        type(FMFD_parameters), allocatable :: fm(:,:,:)
    endtype 
    
    type(FMFD_parameters), allocatable :: fm(:,:,:)
    type(FMFD_parameters), allocatable :: fm_avg(:,:,:)
    type(FMFD_parameters), allocatable :: fm_thread(:,:,:)
    !$OMP THREADPRIVATE(fm_thread)
    type(FMFD_accumulation), allocatable :: acc(:)

    ! FMFD grid
    real(8):: fm0(3)    ! (x0,y0,z0)
    real(8):: fm1(3)    ! (x1,y1,z1)
    real(8):: fm2(3)    ! (x1-x0,y1-y0,z1-z0)
    integer:: nfm(3)    ! (nx,ny,nz)
    real(8):: dfm(3)    ! (dx,dy,dz)

    ! fission source distribution
    real(8), allocatable:: fsd_mc(:,:,:)
    real(8), allocatable:: fsd(:,:,:)
    
    contains 


! =============================================================================
! FMFD_initialize_1
! =============================================================================
subroutine FMFD_allocation()
    integer:: i
    type(FMFD_accumulation), allocatable:: acc(:)

    ! parameters allocation
    allocate(fm(nfm(1),nfm(2),nfm(3)))
    allocate(fm_avg(nfm(1),nfm(2),nfm(3)))
    allocate(fsd_mc(nfm(1),nfm(2),nfm(3)))
    allocate(fsd(nfm(1),nfm(2),nfm(3)))
    allocate(acc(n_acc)) 
    do i = 1, n_acc 
        allocate(acc(i)%fm(nfm(1),nfm(2),nfm(3)))
    enddo
    
    ! area and volume of mesh cell
    a_fm(:) = dfm(1)*dfm(3)
    a_fm(5) = dfm(1)*dfm(2)
    a_fm(6) = a_fm(5)
    v_fm    = dfm(1)*dfm(2)*dfm(3)

end subroutine


! =============================================================================
! FMFD_INITIALIZE_2 initializes parameters used in the FMFD calculaiton
! =============================================================================
subroutine FMFD_initialize()
    integer :: i 

    ! initialization
    fm(:,:,:) % phi     = 0 
    fm(:,:,:) % sig_t   = 0 
    fm(:,:,:) % sig_a   = 0 
    fm(:,:,:) % nusig_f = 0 
    do i = 1, 6
        fm(:,:,:) % Jn(i) = 0 
        fm(:,:,:) % J0(i) = 0 
        fm(:,:,:) % J1(i) = 0 
    enddo 
        
end subroutine

! =============================================================================
! FMFD_INITIALIZE_THREAD initializes thread-wise parameters
! =============================================================================
subroutine FMFD_initialize_thread() 
    integer :: i 
    
    allocate(fm_thread(nfm(1),nfm(2),nfm(3)))
    
    fm_thread(:,:,:) % phi     = 0 
    fm_thread(:,:,:) % sig_t   = 0 
    fm_thread(:,:,:) % sig_a   = 0 
    fm_thread(:,:,:) % nusig_f = 0 
    do i = 1, 6
    fm_thread(:,:,:) % Jn(i)   = 0 
    fm_thread(:,:,:) % J0(i)   = 0 
    fm_thread(:,:,:) % J1(i)   = 0 
    enddo 
    
end subroutine

!    subroutine CMFD_distance (p,i_xyz,idx_xyz,d_CMFD,inside_CMFD,i_surf)
!        type(particle), intent(in) :: p
!        integer, intent(inout) :: i_xyz(3), idx_xyz
!        real(8), intent(inout) :: d_CMFD
!        logical, intent(inout) :: inside_CMFD 
!        integer, intent(inout) :: i_surf
!        real(8) :: xyz(3), uvw(3)
!        real(8) :: d_temp(6)
!        integer :: i,j, i_coord
!        integer :: idx_temp, idx_surf
!        real(8) :: J_temp
!        integer :: a, b, c
!        
!        d_CMFD = INFINITY
!        if (CMFD_lat < 0) return
!        
!        ! Find lattice index from the particle
!        do i = 1, p%n_coord 
!            if (p%coord(i)%lattice == idx_lat) then 
!                i_xyz(1) = p%coord(i)%lattice_x
!                i_xyz(2) = p%coord(i)%lattice_y
!                i_xyz(3) = p%coord(i)%lattice_z
!                
!                xyz(:) = p%coord(i)%xyz(:)
!                uvw(:) = p%coord(i)%uvw(:)
!                i_coord = i
!                exit
!            endif
!        enddo
!        
!        !> Check if the particle is outside the CMFD grid
!        inside_CMFD = .true. 
!        do i = 1, 3
!            if (abs(xyz(i)) > pitch(i)/2.0d0) then 
!                inside_CMFD = .false.
!                exit 
!            endif 
!        enddo 
!        d_CMFD = INFINITY
!        if (inside_CMFD) then
!            d_temp(1) = (-pitch(1)/2.0 - xyz(1))/uvw(1)
!            d_temp(6) = ( pitch(1)/2.0 - xyz(1))/uvw(1)
!            d_temp(2) = (-pitch(2)/2.0 - xyz(2))/uvw(2)
!            d_temp(5) = ( pitch(2)/2.0 - xyz(2))/uvw(2)
!            d_temp(3) = (-pitch(3)/2.0 - xyz(3))/uvw(3)
!            d_temp(4) = ( pitch(3)/2.0 - xyz(3))/uvw(3)
!            
!             
!            do i = 1, 6 
!                if (d_temp(i) < 0) cycle
!                if (d_CMFD > d_temp(i)) then 
!                    d_CMFD = d_temp(i)
!                    i_surf = i
!                endif 
!            enddo 
!        
!        else !> the particle is outside the CMFD grid
!            ! 여기서 xyz는 윗단계 universe의 xyz
!            call distance_cuboid(lattices(idx_lat)%n_xyz, &
!                        p%coord(i_coord-1)%xyz, p%coord(i_coord-1)%uvw, i_xyz, d_CMFD,i_surf)
!            
!        endif
!        a = lattices(idx_lat)%n_xyz(1)
!        b = lattices(idx_lat)%n_xyz(2)
!        c = lattices(idx_lat)%n_xyz(3)
!        
!        idx_xyz = a*b*(i_xyz(3)-1) + a*(i_xyz(2)-1) + i_xyz(1)
!
!    end subroutine
    
! =========================================================================
! FMFD_DISTANCE
! =========================================================================
subroutine FMFD_DISTANCE (p,i_xyz,d_FMFD,inside_FMFD,income_FMFD,i_surf)
    type(particle), intent(in) :: p
    integer, intent(inout) :: i_xyz(3)
    real(8), intent(inout) :: d_FMFD
    logical, intent(inout) :: inside_FMFD 
    integer, intent(inout) :: income_FMFD
    integer, intent(inout) :: i_surf
    real(8) :: xyz(3), uvw(3)
    real(8) :: d_temp(6)
    integer :: i,j, i_coord
    integer :: idx_temp, idx_surf
    real(8) :: J_temp
    integer :: a, b, c
    
    ! Find lattice index in FMFD grid
    i_xyz = FMFD_ID(p%coord(1)%xyz(:))
    
    ! check if the particle is inside the FMFD grid
    inside_FMFD = .true.
    do i = 1, 3
        if ( i_xyz(i) < 1 .or. i_xyz(i) > nfm(i) ) then
            inside_FMFD = .false.
            exit
        end if
    end do
    
    ! check if the particle is coming to the FMFD grid
    if ( .not. inside_FMFD ) call INCOMING(i_xyz(:),income_FMFD)

    ! find the closest distance to the grid surface
    d_FMFD = INFINITY
    uvw(:) = p%coord(1)%uvw(:)
    if ( inside_FMFD .or. income_FMFD /= 0 ) then
        d_temp(1) = ((dfm(1)*(i_xyz(1)-1))-xyz(1))/uvw(1)   ! x0
        d_temp(2) = ((dfm(1)*(i_xyz(1)  ))-xyz(1))/uvw(1)   ! x1
        d_temp(3) = ((dfm(2)*(i_xyz(2)-1))-xyz(2))/uvw(2)   ! y0
        d_temp(4) = ((dfm(2)*(i_xyz(2)  ))-xyz(2))/uvw(2)   ! y1
        d_temp(5) = ((dfm(3)*(i_xyz(3)-1))-xyz(3))/uvw(3)   ! z0
        d_temp(6) = ((dfm(3)*(i_xyz(3)  ))-xyz(3))/uvw(3)   ! z1

        ! inside
        do i = 1, 6
        if ( d_temp(i) > 0 .and. d_FMFD > d_temp(i) ) then
            d_FMFD = d_temp(i)
            i_surf = i
        end if
        end do

        ! income
        if ( i_surf /= income_FMFD ) income_FMFD = 0

    end if



end subroutine

! =============================================================================
! FMFD_ID finds the x, y, z indice in the FMFD mesh grid
! =============================================================================
function FMFD_ID(fmxyz) result(fmid)
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: fmid(3)               ! indice
       
    fmid(:) = int((fmxyz(:)-fm0(:))/dfm(:))+1

end function

! =============================================================================
! ISINCOMING determines if the particle is coming to the FMFD grid
! =============================================================================
subroutine INCOMING(xyz,income)
    integer, intent(in):: xyz(:)
    integer, intent(out):: income

    income = 0
    ! x0
    if ( xyz(1) == 0 ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
            income = 1
        end if
        end if
    end if
    ! x1
    if ( xyz(1) == nfm(1)+1 ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
            income = 2
        end if
        end if
    end if
    ! y0
    if ( xyz(2) == 0 ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
            income = 3
        end if
        end if
    end if
    ! y1
    if ( xyz(2) == nfm(2)+1 ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
            income = 4
        end if
        end if
    end if
    ! z0
    if ( xyz(3) == 0 ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
            income = 5
        end if
        end if
    end if
    ! z1
    if ( xyz(3) == nfm(3)+1 ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
            income = 6
        end if
        end if
    end if

end subroutine

    
!! =============================================================================    
!! DISTANCE_CUBOID ?
!! =============================================================================    
!subroutine distance_cuboid(n_xyz, xyz, uvw, i_xyz, d_CMFD,i_surf)
!    real(8), intent(in) :: xyz(3), uvw(3)
!    integer, intent(in) :: n_xyz(3)
!    integer, intent(inout) :: i_xyz(3)
!    real(8), intent(inout) :: d_CMFD
!    integer, intent(inout) :: i_surf
!    real(8) :: d(6), xyz_(3), temp, r(3), xyz_next(3) 
!    real(8) :: dist
!    integer :: i,j
!    
!    i_surf = 0 
!    xyz_(:) = xyz(:) - lattices(idx_lat)%xyz(:)
!    
!    
!    do i = 1, 3 
!        r(i) = n_xyz(i)*pitch(i) / 2.0d0
!    enddo
!            
!    d(1)   = (-r(1)-xyz_(1))/uvw(1)
!    temp = xyz_(2)+d(1)*uvw(2)
!    if ((temp < -r(2)).or.(temp > r(2))) d(1) = INFINITY
!    temp = xyz_(3)+d(1)*uvw(3)
!    if ((temp < -r(3)).or.(temp > r(3))) d(1) = INFINITY
!    
!    d(6)   = ( r(1)-xyz_(1))/uvw(1)
!    temp = xyz_(2)+d(6)*uvw(2)
!    if ((temp < -r(2)).or.(temp > r(2))) d(6) = INFINITY
!    temp = xyz_(3)+d(6)*uvw(3)
!    if ((temp < -r(3)).or.(temp > r(3))) d(6) = INFINITY
!    
!    d(2)   = (-r(2)-xyz_(2))/uvw(2)
!    temp = xyz_(1)+d(2)*uvw(1)
!    if ((temp < -r(1)).or.(temp > r(1))) d(2) = INFINITY
!    temp = xyz_(3)+d(2)*uvw(3)
!    if ((temp < -r(3)).or.(temp > r(3))) d(2) = INFINITY
!    
!    d(5)   = ( r(2)-xyz_(2))/uvw(2)
!    temp = xyz_(1)+d(5)*uvw(1)
!    if ((temp < -r(1)).or.(temp > r(1))) d(5) = INFINITY
!    temp = xyz_(3)+d(5)*uvw(3)
!    if ((temp < -r(3)).or.(temp > r(3))) d(5) = INFINITY
!    
!    d(3)   = (-r(3)-xyz_(3))/uvw(3)
!    temp = xyz_(1)+d(3)*uvw(1)
!    if ((temp < -r(1)).or.(temp > r(1))) d(3) = INFINITY
!    temp = xyz_(2)+d(3)*uvw(2)
!    if ((temp < -r(2)).or.(temp > r(2))) d(3) = INFINITY
!    
!    d(4)   = ( r(3)-xyz_(3))/uvw(3)
!    temp = xyz_(1)+d(4)*uvw(1)
!    if ((temp < -r(1)).or.(temp > r(1))) d(4) = INFINITY
!    temp = xyz_(2)+d(4)*uvw(2)
!    if ((temp < -r(2)).or.(temp > r(2))) d(4) = INFINITY
!    
!    
!    d_CMFD = INFINITY
!    do i = 1, 6 
!        if (d(i) < 0) cycle
!        if (d_CMFD > d(i)) then 
!            d_CMFD = d(i)
!            i_surf = i
!        endif 
!    enddo 
!    
!    if (i_surf == 0) return 
!    
!    xyz_next(:) = xyz_(:) + (d_CMFD + 10*tiny_bit) * uvw(:)
!    
!    i_xyz = lattice_coord (lattices(idx_lat), xyz_next)
!    
!end subroutine

! =============================================================================
! FMFD_TRK calculates FMFD parameters such as flux, group contstans by 
! track-length estiamtor
! =============================================================================
subroutine FMFD_TRK(wgt,distance,macro_xs,id)
    real(8), intent(in) :: wgt
    real(8), intent(in) :: distance
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(3)
    real(8) :: flux
    
    flux = wgt * distance
    
    fm_thread(id(1),id(2),id(3)) % phi = & 
    fm_thread(id(1),id(2),id(3)) % phi + flux
    fm_thread(id(1),id(2),id(3)) % sig_t = &
    fm_thread(id(1),id(2),id(3)) % sig_t + flux*macro_xs(1)
    fm_thread(id(1),id(2),id(3)) % sig_a = &
    fm_thread(id(1),id(2),id(3)) % sig_a + flux*macro_xs(2)
    fm_thread(id(1),id(2),id(3)) % nusig_f = &
    fm_thread(id(1),id(2),id(3)) % nusig_f + flux*macro_xs(4)
    
end subroutine

! =============================================================================
! FMFD_COL calculates FMFD parameters such as flux, group contstans by 
! collision estimator
! =============================================================================
subroutine FMFD_COL(wgt, macro_xs,id)
    real(8), intent(in) :: wgt
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(3)
    real(8) :: flux
    
    flux = wgt / macro_xs(1)
    
    fm_thread(id(1),id(2),id(3)) % phi = & 
    fm_thread(id(1),id(2),id(3)) % phi + flux
    fm_thread(id(1),id(2),id(3)) % sig_t = &
    fm_thread(id(1),id(2),id(3)) % sig_t + flux*macro_xs(1)
    fm_thread(id(1),id(2),id(3)) % sig_a = &
    fm_thread(id(1),id(2),id(3)) % sig_a + flux*macro_xs(2)
    fm_thread(id(1),id(2),id(3)) % nusig_f = &
    fm_thread(id(1),id(2),id(3)) % nusig_f + flux*macro_xs(4)
    
end subroutine


! =============================================================================
! FMFD_SURF calculates FMFD surface parameters like net and particle current
! =============================================================================
subroutine FMFD_SURF (inside,income, is, id, wgt, bc) 
    logical, intent(in) :: inside
    integer, intent(in) :: income
    integer, intent(in) :: is, id(3)
    real(8), intent(in) :: wgt
    integer, intent(in) :: bc
    

    ! boundary surfaces
!    if ( inside ) then
!        ! x0
!        if ( id(1) == 1 .and. is == 1 .and. fmbc(1) == 2 ) then
!            fm_thread(id(1),id(2),id(3))%J1(1) = &
!            fm_thread(id(1),id(2),id(3))%J1(1) + wgt
!        end if
!        ! x1
!        if ( id(1) == nfm(1) .and. is == 2 .and. fmbc(2) == 2 ) then
!            fm_thread(id(1),id(2),id(3))%J0(2) = &
!            fm_thread(id(1),id(2),id(3))%J0(2) + wgt
!        end if
!        ! y0
!        if ( id(2) == 1 .and. is == 3 .and. fmbc(3) == 2 ) then
!            fm_thread(id(1),id(2),id(3))%J1(3) = &
!            fm_thread(id(1),id(2),id(3))%J1(3) + wgt
!        end if
!        ! y1
!        if ( id(2) == nfm(2) .and. is == 4 .and. fmbc(4) == 2 ) then
!            fm_thread(id(1),id(2),id(3))%J0(4) = &
!            fm_thread(id(1),id(2),id(3))%J0(4) + wgt
!        end if
!        ! z0
!        if ( id(3) == 1 .and. is == 5 .and. fmbc(5) == 2 ) then
!            fm_thread(id(1),id(2),id(3))%J1(5) = &
!            fm_thread(id(1),id(2),id(3))%J1(5) + wgt
!        end if
!        ! y1
!        if ( id(3) == nfm(3) .and. is == 6 .and. fmbc(6) == 2 ) then
!            fm_thread(id(1),id(2),id(3))%J0(6) = &
!            fm_thread(id(1),id(2),id(3))%J0(6) + wgt
!        end if
!    end if

    ! inner nodes
    if ( inside ) then 
        select case(is)
        case(1,3,5)
            fm_thread(id(1),id(2),id(3))%J0(is) = &
            fm_thread(id(1),id(2),id(3))%J0(is) + wgt
        case(2,4,6)
            fm_thread(id(1),id(2),id(3))%J1(is) = &
            fm_thread(id(1),id(2),id(3))%J1(is) + wgt
        end select
        return
    end if

    ! boundary nodes
    select case(income)
    case(1)
        fm_thread(id(1)+1,id(2),id(3))%J1(1) = &
        fm_thread(id(1)+1,id(2),id(3))%J1(1) + wgt
    case(2)
        fm_thread(id(1)-1,id(2),id(3))%J0(2) = &
        fm_thread(id(1)-1,id(2),id(3))%J0(2) + wgt
    case(3)
        fm_thread(id(1),id(2)+1,id(3))%J1(3) = &
        fm_thread(id(1),id(2)+1,id(3))%J1(3) + wgt
    case(4)
        fm_thread(id(1),id(2)-1,id(3))%J0(4) = &
        fm_thread(id(1),id(2)-1,id(3))%J0(4) + wgt
    case(5)
        fm_thread(id(1),id(2),id(3)+1)%J1(5) = &
        fm_thread(id(1),id(2),id(3)+1)%J1(5) + wgt
    case(6)
        fm_thread(id(1),id(2),id(3)-1)%J0(6) = &
        fm_thread(id(1),id(2),id(3)-1)%J0(6) + wgt
    end select
            
end subroutine

! =============================================================================
! NORM_FMFD normalizes cycle-wise FMFD parameters
! =============================================================================
subroutine NORM_FMFD()
    integer:: i
    real(8):: aa, bb    ! parameters

!    ! volume quantity normalization
!    aa = dble(ngen)*v_fm
!    fm_thread(:,:,:) % phi     = fm_thread(:,:,:) % phi     / aa
!    fm_thread(:,:,:) % sig_t   = fm_thread(:,:,:) % sig_t   / aa
!    fm_thread(:,:,:) % sig_a   = fm_thread(:,:,:) % sig_a   / aa
!    fm_thread(:,:,:) % nusig_f = fm_thread(:,:,:) % nusig_f / aa
!
!    ! surface quantity normalization
!    bb = dble(ngen)*a_fm(1)
!    do i = 1, 2
!    fm_thread(:,:,:) % J0(i) = fm_thread(:,:,:) % J0(i) / bb
!    fm_thread(:,:,:) % J1(i) = fm_thread(:,:,:) % J1(i) / bb
!    end do
!    bb = dble(ngen)*a_fm(2)
!    do i = 3, 4
!    fm_thread(:,:,:) % J0(i) = fm_thread(:,:,:) % J0(i) / bb
!    fm_thread(:,:,:) % J1(i) = fm_thread(:,:,:) % J1(i) / bb
!    end do
!    bb = dble(ngen)*a_fm(3)
!    do i = 5, 6
!    fm_thread(:,:,:) % J0(i) = fm_thread(:,:,:) % J0(i) / bb
!    fm_thread(:,:,:) % J1(i) = fm_thread(:,:,:) % J1(i) / bb
!    end do
    
    !> gather thread FMFD parameters
    fm(:,:,:) % phi     = fm(:,:,:) % phi     + fm_thread(:,:,:)%phi
    fm(:,:,:) % sig_t   = fm(:,:,:) % sig_t   + fm_thread(:,:,:)%sig_t 
    fm(:,:,:) % sig_a   = fm(:,:,:) % sig_a   + fm_thread(:,:,:)%sig_a 
    fm(:,:,:) % nusig_f = fm(:,:,:) % nusig_f + fm_thread(:,:,:)%nusig_f 
    do i = 1, 6
    fm(:,:,:) % J0(i) = fm(:,:,:) % J0(i) + fm_thread(:,:,:)%J0(i)
    fm(:,:,:) % J1(i) = fm(:,:,:) % J1(i) + fm_thread(:,:,:)%J1(i)
    end do

    ! volume quantity normalization
    aa = dble(ngen)*v_fm
    fm(:,:,:) % phi     = fm(:,:,:) % phi     / aa
    fm(:,:,:) % sig_t   = fm(:,:,:) % sig_t   / aa
    fm(:,:,:) % sig_a   = fm(:,:,:) % sig_a   / aa
    fm(:,:,:) % nusig_f = fm(:,:,:) % nusig_f / aa

    ! surface quantity normalization
    do i = 1, 6
    bb = dble(ngen)*a_fm(i)
    fm(:,:,:) % J0(i) = fm(:,:,:) % J0(i) / bb
    fm(:,:,:) % J1(i) = fm(:,:,:) % J1(i) / bb
    end do

end subroutine


! =============================================================================
! PROCESS_FMFD deals with MPI process and average quantities
! =============================================================================
subroutine PROCESS_FMFD() 
    !> MPI derived type reduce parameters 
    integer, dimension(7) :: FMFD_blocklength, FMFD_displacement, FMFD_datatype 
    integer :: intex, realex, restype, FMFD_op  ! MPI int & real extern / new type
    type(FMFD_parameters), allocatable :: FMFD_MPI_slut(:,:,:)
    integer :: i, j, k, l
    
    data FMFD_blocklength /1, 1, 1, 1, 6, 6, 6 /
    
    !> Gather FMFD parameters from the slave nodes 
    call MPI_TYPE_EXTENT(MPI_REAL8, realex, ierr) 
    FMFD_displacement(1) = 0
    FMFD_displacement(2) = realex
    FMFD_displacement(3) = 2*realex
    FMFD_displacement(4) = 3*realex
    FMFD_displacement(5) = 4*realex
    FMFD_displacement(6) = (4+6)*realex
    FMFD_displacement(7) = (4+2*6)*realex
    
    FMFD_datatype(:) = MPI_REAL8
    
    call MPI_TYPE_STRUCT (7, FMFD_blocklength, FMFD_displacement, FMFD_datatype,restype, ierr)
    call MPI_TYPE_COMMIT (restype, ierr)
    call MPI_Op_create(FMFD_SUM, .true. , FMFD_op, ierr)
    call MPI_REDUCE(fm, FMFD_MPI_slut, nfm(1)*nfm(2)*nfm(3), restype, FMFD_op, score, MPI_COMM_WORLD, ierr)
    fm = FMFD_MPI_slut
     
    call MPI_TYPE_FREE(restype, ierr)
    deallocate(FMFD_MPI_slut)
    call MPI_Op_Free(FMFD_Op, ierr)
    
    if ( icore /= score ) return 

!    ! data swapping
!    do i = 1, nfm(1)
!    do j = 1, nfm(2)
!    do k = 1, nfm(3)
!        if ( i /= 1 )      fm(i,j,k)%J1(1) = fm(i-1,j,k)%J1(2)
!        if ( i /= nfm(1) ) fm(i,j,k)%J0(2) = fm(i+1,j,k)%J0(1)
!        if ( j /= 1 )      fm(i,j,k)%J1(3) = fm(i-1,j,k)%J0(4)
!        if ( j /= nfm(2) ) fm(i,j,k)%J0(4) = fm(i+1,j,k)%J1(3)
!        if ( k /= 1 )      fm(i,j,k)%J1(5) = fm(i-1,j,k)%J0(6)
!        if ( k /= nfm(3) ) fm(i,j,k)%J0(6) = fm(i+1,j,k)%J1(5)
!        fm(i,j,k)%Jn(:)  = fm(i,j,k)%J1(:) - fm(i,j,k)%J0(:)
!    end do
!    end do
!    end do

    ! net current
    do i = 1, 6
    fm(:,:,:) % Jn(i) = fm(:,:,:) % J1(i) - fm(:,:,:) % J0(i)
    end do
    
    ! group constant
    fm(:,:,:) % sig_t   = fm(:,:,:) % sig_t   / fm(:,:,:) % phi
    fm(:,:,:) % sig_a   = fm(:,:,:) % sig_a   / fm(:,:,:) % phi
    fm(:,:,:) % nusig_f = fm(:,:,:) % nusig_f / fm(:,:,:) % phi

    ! save the next accumulation
    do i = n_acc, 2, -1
        acc(i)%fm(:,:,:) = acc(i-1)%fm(:,:,:)
    enddo 
    acc(1)%fm(:,:,:) = fm(:,:,:)
    
    !> average with the accumulated parameters
    fm_avg(:,:,:)%phi     = 0 
    fm_avg(:,:,:)%sig_t   = 0 
    fm_avg(:,:,:)%sig_a   = 0 
    fm_avg(:,:,:)%nusig_f = 0 
    do i = 1, 6 
      fm_avg(:,:,:)%Jn(i) = 0 
      fm_avg(:,:,:)%J0(i) = 0 
      fm_avg(:,:,:)%J1(i) = 0 
    enddo 
    
    ! accumulation
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
    do l = 1, n_acc
       fm_avg(i,j,k)%phi     = fm_avg(i,j,k)%phi     + acc(l)%fm(i,j,k)%phi
       fm_avg(i,j,k)%sig_t   = fm_avg(i,j,k)%sig_t   + acc(l)%fm(i,j,k)%sig_t
       fm_avg(i,j,k)%sig_a   = fm_avg(i,j,k)%sig_a   + acc(l)%fm(i,j,k)%sig_a
       fm_avg(i,j,k)%nusig_f = fm_avg(i,j,k)%nusig_f + acc(l)%fm(i,j,k)%nusig_f

       fm_avg(i,j,k)%Jn(:) = fm_avg(i,j,k)%Jn(:) + acc(l)%fm(i,j,k)%Jn(:)
       fm_avg(i,j,k)%J0(:) = fm_avg(i,j,k)%J0(:) + acc(l)%fm(i,j,k)%J0(:)
       fm_avg(i,j,k)%J1(:) = fm_avg(i,j,k)%J1(:) + acc(l)%fm(i,j,k)%J1(:)
    end do
    end do
    end do
    end do
    
    ! average
    fm_avg(:,:,:)%phi     = fm_avg(:,:,:)%phi     / dble(n_acc)
    fm_avg(:,:,:)%sig_t   = fm_avg(:,:,:)%sig_t   / dble(n_acc)
    fm_avg(:,:,:)%sig_a   = fm_avg(:,:,:)%sig_a   / dble(n_acc)
    fm_avg(:,:,:)%nusig_f = fm_avg(:,:,:)%nusig_f / dble(n_acc)
    do i = 1, 6
    fm_avg(:,:,:)%Jn(i) = fm_avg(:,:,:)%Jn(i) / dble(n_acc)
    fm_avg(:,:,:)%J0(i) = fm_avg(:,:,:)%J0(i) / dble(n_acc)
    fm_avg(:,:,:)%J1(i) = fm_avg(:,:,:)%J1(i) / dble(n_acc)
    enddo
    
end subroutine 


! =============================================================================
! FMFD_SOLVE solves FMFD eigenvalue problem
! =============================================================================
subroutine FMFD_SOLVE(keff,fsd)
    real(8), intent(in):: keff              ! multiplication factor
    real(8), intent(inout) :: fsd(:,:,:)    ! fission source distribution
    real(8) :: M(nfm(1),nfm(2),nfm(3),7)    ! FMFD matrix
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: &
        sig_a, &      ! absorption
        nusig_f, &    ! nu X fission
        D, &          ! diffusion coefficient
        F, &          ! matrix component for source
        phi0, &       ! neutron flux (previous)
        phi1          ! neutron flux (current)
    real(8), dimension(nfm(1),nfm(2),nfm(3),6):: &
        D_tilda, &    ! D tilda
        D_hat, &      ! correction factor
        Jn, &         ! net current
        J0, &         ! partial current -
        J1            ! partial current +
    real(8) :: k_fmfd
    integer :: i, j, k
            
    ! copy parameters
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        sig_a(i,j,k)   = fm_avg(i,j,k)%sig_a
        nusig_f(i,j,k) = fm_avg(i,j,k)%nusig_f 
        D(i,j,k)       = 1D0 / (3D0 * fm_avg(i,j,k)%sig_t) 
        Jn(i,j,k,:)    = fm_avg(i,j,k)%Jn(:) 
        J0(i,j,k,:)    = fm_avg(i,j,k)%J0(:) 
        J1(i,j,k,:)    = fm_avg(i,j,k)%J1(:) 
        phi1(i,j,k)    = fm_avg(i,j,k)%phi
    enddo 
    enddo
    enddo
    
    !> calculate D_tilda & D_hat 
    call D_TILDA_CALCULATION(D,D_tilda)
    call D_HAT_CALCULATION(D_tilda,Jn,phi1,D_hat)
     
    ! matrix composition
    if ( FMFD_type == 1 ) then 
        !> CMFD 
        M = getM(D_tilda, D_hat, sig_a)
    else!if (CMFD_type == 2) then 
        !> p-CMFD
        !M = getMp (phi, D, J_pp, J_pn, sig_a, a,b,c,pitch(1),pitch(2),pitch(3))
    endif 
    F = nusig_f
    
    k_fmfd = keff
    call POWER (k_fmfd, M, F, phi0, phi1, nusig_f)
    
    !print *, 'CMFD keff  ',keff_CMFD
    !write(prt_keff,*) keff_CMFD, k_col, k_tl
    
    !> CMFD feedback (modulation)
    fsd_MC(:,:,:) = fsd_MC(:,:,:) / sum(fsd_MC)
    phi1(:,:,:)   = phi1(:,:,:) / sum(phi1)
    fsd(:,:,:)    = phi1(:,:,:) / fsd_MC(:,:,:)
    
end subroutine

! =============================================================================
! D_TILDA_CALCULATION
! =============================================================================
subroutine D_TILDA_CALCULATION(D,Dt)
    implicit none
    real(8), intent(in) :: D(:,:,:)
    real(8), intent(out):: Dt(:,:,:,:)
    integer:: i, j, k
    
    ! initialization
    Dt = 0

    ! inner region
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 ) &      ! x0
        Dt(i,j,k,1) = 2D0*D(i,j,k)*D(i-1,j,k)/(D(i,j,k)+D(i-1,j,k))/dfm(1)
        if ( i /= nfm(1) ) & ! x1
        Dt(i,j,k,2) = 2D0*D(i+1,j,k)*D(i,j,k)/(D(i+1,j,k)+D(i,j,k))/dfm(1)
        if ( j /= 1 ) &      ! y0
        Dt(i,j,k,3) = 2D0*D(i,j,k)*D(i,j-1,k)/(D(i,j,k)+D(i,j-1,k))/dfm(2)
        if ( j /= nfm(2) ) & ! y1
        Dt(i,j,k,4) = 2D0*D(i,j+1,k)*D(i,j,k)/(D(i,j+1,k)+D(i,j,k))/dfm(2)
        if ( k /= 1 ) &      ! y0
        Dt(i,j,k,5) = 2D0*D(i,j,k)*D(i,j,k-1)/(D(i,j,k)+D(i,j,k-1))/dfm(3)
        if ( k /= nfm(3) ) & ! y1
        Dt(i,j,k,6) = 2D0*D(i,k,k+1)*D(i,j,k)/(D(i,j,k+1)+D(i,j,k))/dfm(3)
    end do
    end do
    end do

end subroutine


! =============================================================================
! D_HAT_CALCULATION calculates correction factors
! =============================================================================
subroutine D_HAT_CALCULATION(Dt,Jn,phi,Dh)
    real(8), intent(in)   :: Dt(:,:,:,:), Jn(:,:,:,:), phi(:,:,:)
    real(8), intent(inout):: Dh(:,:,:,:)
    integer :: i, j, k

    ! inner region
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 ) &      ! x0
        Dh(i,j,k,1) = (Jn(i,j,k,1)+Dt(i,j,k,1)*(phi(i,j,k)-phi(i-1,j,k)))/ &
                    (phi(i,j,k)+phi(i-1,j,k))
        if ( i /= nfm(1) ) & ! x1
        Dh(i,j,k,2) = (Jn(i,j,k,2)+Dt(i,j,k,2)*(phi(i+1,j,k)-phi(i,j,k)))/ &
                    (phi(i+1,j,k)+phi(i,j,k))
        if ( j /= 1 ) &      ! y0
        Dh(i,j,k,3) = (Jn(i,j,k,3)+Dt(i,j,k,3)*(phi(i,j,k)-phi(i,j-1,k)))/ &
                    (phi(i,j,k)+phi(i,j-1,k))
        if ( j /= nfm(2) ) & ! y1
        Dh(i,j,k,4) = (Jn(i,j,k,4)+Dt(i,j,k,4)*(phi(i,j+1,k)-phi(i,j,k)))/ &
                    (phi(i,j+1,k)+phi(i,j,k))
        if ( k /= 1 ) &      ! y0
        Dh(i,j,k,5) = (Jn(i,j,k,5)+Dt(i,j,k,5)*(phi(i,j,k)-phi(i,j,k-1)))/ &
                    (phi(i,j,k)+phi(i,j,k-1))
        if ( k /= nfm(3) ) & ! y1
        Dh(i,j,k,6) = (Jn(i,j,k,6)+Dt(i,j,k,6)*(phi(i,j,k+1)-phi(i,j,k)))/ &
                    (phi(i,j,k+1)+phi(i,j,k))
    end do
    end do
    end do

    ! boundary
    i = 1;      Dh(i,:,:,1) = Jn(i,:,:,1)/phi(i,:,:)
    i = nfm(1); Dh(i,:,:,2) = Jn(i,:,:,2)/phi(i,:,:)
    j = 1;      Dh(:,j,:,3) = Jn(:,j,:,3)/phi(:,j,:)
    j = nfm(2); Dh(:,j,:,4) = Jn(:,j,:,4)/phi(:,j,:)
    k = 1;      Dh(:,:,k,5) = Jn(:,:,k,5)/phi(:,:,k)
    k = nfm(3); Dh(:,:,k,6) = Jn(:,:,k,6)/phi(:,:,k)

end subroutine


! ========================================================= !
!  Miscel. functions for CMFD_solve
! ========================================================= !

! function which produces M matrix from xs & D_hat
function getM (Dt, Dh, sig_a) result(M)
    real(8), intent(in) :: Dt(:,:,:,:)
    real(8), intent(in) :: Dh(:,:,:,:)
    real(8), intent(in) :: sig_a(:,:,:)
    real(8):: M(nfm(1),nfm(2),nfm(3),7)
    integer :: i, j, k
    
    ! initialization
    M = 0 
    
    ! M matrix set 
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 )      M(i,j,k,1) = -(Dt(i,j,k,1)+Dh(i,j,k,1))/dfm(1)
        if ( i /= nfm(1) ) M(i,j,k,2) = (-Dt(i,j,k,2)+Dh(i,j,k,2))/dfm(2)
        if ( j /= 1 )      M(i,j,k,3) = -(Dt(i,j,k,3)+Dh(i,j,k,3))/dfm(2)
        if ( j /= nfm(2) ) M(i,j,k,4) = (-Dt(i,j,k,4)+Dh(i,j,k,4))/dfm(2)
        if ( k /= 1 )      M(i,j,k,5) = -(Dt(i,j,k,5)+Dh(i,j,k,5))/dfm(3)
        if ( k /= nfm(3) ) M(i,j,k,6) = (-Dt(i,j,k,6)+Dh(i,j,k,6))/dfm(3)
        
        M(i,j,k,7) = &
            +(Dt(i,j,k,1)-Dh(i,j,k,1)+Dt(i,j,k,2)+Dh(i,j,k,2))/dfm(1) &
            +(Dt(i,j,k,3)-Dh(i,j,k,3)+Dt(i,j,k,4)+Dh(i,j,k,4))/dfm(2) &
            +(Dt(i,j,k,5)-Dh(i,j,k,5)+Dt(i,j,k,6)+Dh(i,j,k,6))/dfm(3) &
            +sig_a(i,j,k)
    enddo
    enddo
    enddo
    
end function getM


subroutine POWER (k_eff1, M, F, nusigf, phi0, phi1)
    real(8), intent(inout):: k_eff1
    real(8), intent(in) :: M(:,:,:,:), F(:,:,:)
    real(8), intent(in) :: nusigf(:,:,:)
    real(8), intent(inout):: phi0(:,:,:), phi1(:,:,:)
    integer :: iter_max = 1000
    integer :: iter
    real(8) :: k_eff0
    real(8) :: err
    
    err = 1; iter_max = 1000; iter = 0
    do while ( ( err > 1.0d-10 ) .and. (iter < iter_max) )
        iter = iter + 1
        k_eff0 = k_eff1
        phi0 = phi1
        phi1 = BiCGStab_hepta(M(:,:,:,:),F(:,:,:)/phi0(:,:,:)/k_eff0)
        k_eff1 = k_eff1*sum(F*phi1*F*phi1)/sum(F*phi1*F*phi0)
        err = maxval((phi1(:,:,:)-phi0(:,:,:))/phi1(:,:,:))
    enddo
    
end subroutine 

! ============================================== !
!         Hepta diagonal matrix solvers            !
! ============================================== !
function BiCGStab_hepta(M,Q) result(x)
    real(8), intent(in) :: M (:,:,:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: x, r, rs, v, p, s, t
    real(8), parameter :: e = 1d-10
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    real(8) :: summesion, temp
    integer :: it = 0
    integer :: i, j, k

    x     = 0.0
    r     = Q
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sqrt(sum(r*r))
    norm_b = sqrt(sum(Q*Q))
    
    do while ( norm_r .GT. e*norm_b )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
    
        p        = r + beta * (p - omega*v)
    
!        do i = 1, index0  !> v = matmul(M,p)
!            if (M(i,1).ne.0) v(i) = v(i) + M(i,1)*p(i-1)
!            if (M(i,6).ne.0) v(i) = v(i) + M(i,6)*p(i+1)
!            if (M(i,2).ne.0) v(i) = v(i) + M(i,2)*p(i-a)
!            if (M(i,5).ne.0) v(i) = v(i) + M(i,5)*p(i+a)
!            if (M(i,3).ne.0) v(i) = v(i) + M(i,3)*p(i-a*b)
!            if (M(i,4).ne.0) v(i) = v(i) + M(i,4)*p(i+a*b)
!                             v(i) = v(i) + M(i,7)*p(i)
!        enddo 

        v(:,:,:) = 0
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,1)*p(i-1,j,k)
            if ( i /= nfm(1) ) v(i,j,k) = v(i,j,k) + M(i,j,k,2)*p(i+1,j,k)
            if ( j /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,3)*p(i,j-1,k)
            if ( j /= nfm(2) ) v(i,j,k) = v(i,j,k) + M(i,j,k,4)*p(i,j+1,k)
            if ( k /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,5)*p(i,j,k-1)
            if ( k /= nfm(3) ) v(i,j,k) = v(i,j,k) + M(i,j,k,6)*p(i,j,k+1)
                               v(i,j,k) = v(i,j,k) + M(i,j,k,7)*p(i,j,k)
        end do
        end do
        end do
        
!        t(:) = 0
!        do i = 1, index0  !> t = matmul(M,s)
!            if (M(i,1).ne.0) t(i) = t(i) + M(i,1)*s(i-1)
!            if (M(i,6).ne.0) t(i) = t(i) + M(i,6)*s(i+1)
!            if (M(i,2).ne.0) t(i) = t(i) + M(i,2)*s(i-a)
!            if (M(i,5).ne.0) t(i) = t(i) + M(i,5)*s(i+a)
!            if (M(i,3).ne.0) t(i) = t(i) + M(i,3)*s(i-a*b)
!            if (M(i,4).ne.0) t(i) = t(i) + M(i,4)*s(i+a*b)
!                             t(i) = t(i) + M(i,7)*s(i)
!        enddo 
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,1)*s(i-1,j,k)
            if ( i /= nfm(1) ) t(i,j,k) = t(i,j,k) + M(i,j,k,2)*s(i+1,j,k)
            if ( j /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,3)*s(i,j-1,k)
            if ( j /= nfm(2) ) t(i,j,k) = t(i,j,k) + M(i,j,k,4)*s(i,j+1,k)
            if ( k /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,5)*s(i,j,k-1)
            if ( k /= nfm(3) ) t(i,j,k) = t(i,j,k) + M(i,j,k,6)*s(i,j,k+1)
                               t(i,j,k) = t(i,j,k) + M(i,j,k,7)*s(i,j,k)
        end do
        end do
        end do
        
        omega    = sum(t*s)/sum(t*t)
        x        = x + alpha*p + omega*s
        r        = s - omega*t
        norm_r   = sqrt(sum(r*r))
        norm_b   = sqrt(sum(Q*Q))
    
        it = it + 1
    end do   
    
end function BiCGStab_hepta     


! function which produces p-CMFD M matrix
!function getMp (phi, D, J_pp, J_pn, sig_a, a,b,c,dx,dy,dz) result (M)
!    integer, intent(in) :: a,b,c
!    real(8), intent(in) :: phi(a*b*c), sig_a(a*b*c), D(a*b*c), J_pp(a*b*c,6), J_pn(a*b*c,6), dx, dy, dz
!    real(8) :: M(a*b*c,7)
!    real(8) :: D_tilda(a*b*c,6), D_hat_p(a*b*c,6), D_hat_n(a*b*c,6)
!    integer :: xyz(3)
!    integer :: i, j, index, index0
!    
!    M(:,:) = 0
!    index0 = a*b*c
!    
!    ! D_tilda set 
!    do i = 1, index0
!        xyz = getXYZ(i,a,b,c)
!        M(i,7) = 0 
!        if (xyz(1).ne.1) then 
!            D_tilda(i,1) = 2.0*D(i)*D(i-1)/(D(i)+D(i-1))/dx
!            D_hat_p(i,1) = -(1./2./phi(i-1))*(2*J_pp(i,1) + D_tilda(i,1)*(phi(i)-phi(i-1)))
!            D_hat_n(i,1) =  (1./2./phi(i))*(2*J_pn(i,1)   - D_tilda(i,1)*(phi(i)-phi(i-1)))
!            M(i,7) = M(i,7) + (D_tilda(i,1)+D_hat_n(i,1))/dx
!        else 
!            D_tilda(i,1) = 0
!            D_hat_p(i,1) = -J_pp(i,1)/phi(i)
!            D_hat_n(i,1) =  J_pn(i,1)/phi(i)
!            M(i,7) = M(i,7) + (D_hat_p(i,1) + D_hat_n(i,1))/dx
!        endif
!        if (xyz(1).ne.a) then 
!            D_tilda(i,6) = 2.0*D(i)*D(i+1)/(D(i)+D(i+1))/dx
!            D_hat_p(i,6) = -(1./2./phi(i))*(2*J_pp(i,6) + D_tilda(i,6)*(phi(i+1)-phi(i)))
!            D_hat_n(i,6) =  (1./2./phi(i+1))*(2*J_pn(i,6) - D_tilda(i,6)*(phi(i+1)-phi(i)))
!            M(i,7) = M(i,7) + (D_tilda(i,6)-D_hat_p(i,6))/dx
!        else 
!            D_tilda(i,6) = 0
!            D_hat_p(i,6) = -J_pp(i,6)/phi(i)
!            D_hat_n(i,6) =  J_pn(i,6)/phi(i)
!            M(i,7) = M(i,7) - (D_hat_p(i,6) + D_hat_n(i,6))/dx
!        endif 
!        
!        if (xyz(2).ne.1) then 
!            D_tilda(i,2) = 2.0*D(i)*D(i-a)/(D(i)+D(i-a))/dy
!            D_hat_p(i,2) = -(1./2./phi(i-a))*(2*J_pp(i,2) + D_tilda(i,2)*(phi(i)-phi(i-a)))
!            D_hat_n(i,2) =  (1./2./phi(i))  *(2*J_pn(i,2) - D_tilda(i,2)*(phi(i)-phi(i-a)))
!            M(i,7) = M(i,7) + (D_tilda(i,2)+D_hat_n(i,2))/dy 
!        else 
!            D_tilda(i,2) = 0
!            D_hat_p(i,2) = -J_pp(i,2)/phi(i)
!            D_hat_n(i,2) =  J_pn(i,2)/phi(i)
!            M(i,7) = M(i,7) + (D_hat_p(i,2) + D_hat_n(i,2))/dy
!        endif
!        if (xyz(2).ne.b) then 
!            D_tilda(i,5) = 2.0*D(i)*D(i+a)/(D(i)+D(i+a))/dy
!            D_hat_p(i,5) = -(1./2./phi(i))*(2*J_pp(i,5)   + D_tilda(i,5)*(phi(i+a)-phi(i)))
!            D_hat_n(i,5) =  (1./2./phi(i+a))*(2*J_pn(i,5) - D_tilda(i,5)*(phi(i+a)-phi(i)))
!            M(i,7) = M(i,7) + (D_tilda(i,5)-D_hat_p(i,5))/dy
!        else 
!            D_tilda(i,5) = 0
!            D_hat_p(i,5) = -J_pp(i,5)/phi(i)
!            D_hat_n(i,5) =  J_pn(i,5)/phi(i)
!            M(i,7) = M(i,7) - (D_hat_p(i,5) + D_hat_n(i,5))/dy
!        endif 
!        
!        if (xyz(3).ne.1) then 
!            D_tilda(i,3) = 2*D(i)*D(i-b*a)/(D(i)+D(i-b*a))/dz
!            D_hat_p(i,3) = -(1./2./phi(i-a*b))*(2*J_pp(i,3) + D_tilda(i,3)*(phi(i)-phi(i-a*b)))
!            D_hat_n(i,3) =  (1./2./phi(i))*(2*J_pn(i,3)     - D_tilda(i,3)*(phi(i)-phi(i-a*b)))
!            M(i,7) = M(i,7) + (D_tilda(i,3)+D_hat_n(i,3))/dz
!        else 
!            D_tilda(i,3) = 0
!            D_hat_p(i,3) = -J_pp(i,3)/phi(i)
!            D_hat_n(i,3) =  J_pn(i,3)/phi(i)
!            M(i,7) = M(i,7) + (D_hat_p(i,3) + D_hat_n(i,3))/dz
!        endif 
!        if (xyz(3).ne.c) then 
!            D_tilda(i,4) = 2.0*D(i)*D(i+b*a)/(D(i)+D(i+b*a))/dz
!            D_hat_p(i,4) = -(1./2./phi(i))*(2*J_pp(i,4)     + D_tilda(i,4)*(phi(i+a*b)-phi(i)))
!            D_hat_n(i,4) =  (1./2./phi(i+a*b))*(2*J_pn(i,4) - D_tilda(i,4)*(phi(i+a*b)-phi(i)))
!            M(i,7) = M(i,7) + (D_tilda(i,4)-D_hat_p(i,4))/dz
!        else 
!            D_tilda(i,4) = 0
!            D_hat_p(i,4) = -J_pp(i,4)/phi(i)
!            D_hat_n(i,4) =  J_pn(i,4)/phi(i)
!            M(i,7) = M(i,7) - (D_hat_p(i,4) + D_hat_n(i,4))/dz
!        endif
!    
!        ! M matrix set
!        M(i,1) = -(D_tilda(i,1)-D_hat_p(i,1))/dx
!        M(i,6) = -(D_tilda(i,6)+D_hat_n(i,6))/dx
!        M(i,2) = -(D_tilda(i,2)-D_hat_p(i,2))/dy
!        M(i,5) = -(D_tilda(i,5)+D_hat_n(i,5))/dy
!        M(i,3) = -(D_tilda(i,3)-D_hat_p(i,3))/dz
!        M(i,4) = -(D_tilda(i,4)+D_hat_n(i,4))/dz
!        
!        M(i,7) = M(i,7) + sig_a(i)
!        
!        xyz = getXYZ(i, a,b,c)
!        
!        if (xyz(1) == 1) M(i,1) = 0
!        if (xyz(1) == a) M(i,6) = 0
!        if (xyz(2) == 1) M(i,2) = 0
!        if (xyz(2) == b) M(i,5) = 0
!        if (xyz(3) == 1) M(i,3) = 0
!        if (xyz(3) == c) M(i,4) = 0
!        
!    enddo
!    
!end function getMp



!subroutine conjgrad_hepta(M,Q,x,index,a, b)  
!    integer :: index, i, iter, j, a, b
!    real(8), dimension(:,:) :: M(index, 7)
!    real(8), dimension(index)   :: Q, x, r, p, Ap
!    real(8) :: rsold, alpha, rsnew
!    
!    r = Q
!    p = r
!    rsold = dot_product(r,r)
!
!    rsnew = 1
!    iter = 0 
!    do while (sqrt(rsnew).gt.1.0d-10)
!        iter = iter+1 
!        ! Ap = A*p
!        Ap(:) = 0
!        do i = 1, index
!            if (M(i,1).ne.0) Ap(i) = Ap(i) + M(i,1)*p(i-1)
!            if (M(i,6).ne.0) Ap(i) = Ap(i) + M(i,6)*p(i+1)
!            if (M(i,2).ne.0) Ap(i) = Ap(i) + M(i,2)*p(i-a)
!            if (M(i,5).ne.0) Ap(i) = Ap(i) + M(i,5)*p(i+a)
!            if (M(i,3).ne.0) Ap(i) = Ap(i) + M(i,3)*p(i-a*b)
!            if (M(i,4).ne.0) Ap(i) = Ap(i) + M(i,4)*p(i+a*b)
!                             Ap(i) = Ap(i) + M(i,7)*p(i)
!        enddo 
!        
!        alpha = rsold/dot_product(p,Ap)
!        
!        do j = 1, index
!            x(j) = x(j) + alpha*p(j)
!            r(j) = r(j) - alpha*Ap(j)
!        enddo
!        
!        rsnew = dot_product(r,r)
!        
!        do j = 1, index
!            p(j) = r(j) + (rsnew/rsold)*p(j)
!        enddo
!        
!        rsold = rsnew
!    enddo
!    
!end subroutine conjgrad_hepta



subroutine CMFD_SUM (inpar, inoutpar, partype) 
    type(FMFD_parameters), intent(in) :: inpar(nfm(1),nfm(2),nfm(3))
    type(FMFD_parameters), intent(inout) :: inoutpar(nfm(1),nfm(2),nfm(3))
    integer :: ii, partype
    
    inoutpar(:,:,:)%phi     = inoutpar(:,:,:)%phi     + inpar(:,:,:)%phi
    inoutpar(:,:,:)%sig_t   = inoutpar(:,:,:)%sig_t   + inpar(:,:,:)%sig_t
    inoutpar(:,:,:)%sig_a   = inoutpar(:,:,:)%sig_a   + inpar(:,:,:)%sig_a
    inoutpar(:,:,:)%nusig_f = inoutpar(:,:,:)%nusig_f + inpar(:,:,:)%nusig_f
    do ii = 1, 6 
    inoutpar(:,:,:)%Jn(ii) = inoutpar(:,:,:)%Jn(ii) + inpar(:,:,:)%Jn(ii)
    inoutpar(:,:,:)%J0(ii) = inoutpar(:,:,:)%J0(ii) + inpar(:,:,:)%J0(ii)
    inoutpar(:,:,:)%J1(ii) = inoutpar(:,:,:)%J1(ii) + inpar(:,:,:)%J1(ii)
    enddo 

end subroutine  

    
    
end module     
