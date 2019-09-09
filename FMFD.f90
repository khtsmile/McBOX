module FMFD
    use omp_lib 
    use mpi 
    use geometry_header,    only : lattices, find_lat_idx, lattice_coord
    use variables
    use constants,          only : INFINITY, TiNY_BIT, prt_keff
    use particle_header,    only : particle
    use FMFD_HEADER
    
    implicit none
    
    contains 


! =============================================================================
! FMFD_initialize_1
! =============================================================================
subroutine FMFD_allocation()
    integer:: i, j

    ! parameters allocation
    allocate(fm(nfm(1),nfm(2),nfm(3)))
    allocate(fm_avg(nfm(1),nfm(2),nfm(3)))
    allocate(fsd_MC(nfm(1),nfm(2),nfm(3)))
    allocate(fsd(nfm(1),nfm(2),nfm(3)))
    allocate(acc(n_acc)) 
    do i = 1, n_acc 
        allocate(acc(i)%fm(nfm(1),nfm(2),nfm(3)))
        acc(i)%fm(:,:,:)%phi     = 0
        acc(i)%fm(:,:,:)%sig_t   = 0
        acc(i)%fm(:,:,:)%sig_a   = 0
        acc(i)%fm(:,:,:)%nusig_f = 0
        do j = 1, 6
        acc(i)%fm(:,:,:)%Jn(j)   = 0
        acc(i)%fm(:,:,:)%J0(j)   = 0
        acc(i)%fm(:,:,:)%J1(j)   = 0
        acc(i)%fm(:,:,:)%sphi(j) = 0
        end do
    enddo
    
    ! area and volume of mesh cell
    a_fm(:) = dfm(1)*dfm(3)
    a_fm(5) = dfm(1)*dfm(2)
    a_fm(6) = a_fm(5)
    v_fm    = dfm(1)*dfm(2)*dfm(3)


    ! CMFD parameters
    if ( cmfdon ) then
        allocate(cm_t   (ncm(1),ncm(2),ncm(3)), &
                 cmD    (ncm(1),ncm(2),ncm(3)), &
                 cm_a   (ncm(1),ncm(2),ncm(3)), &
                 cm_nf  (ncm(1),ncm(2),ncm(3)), &
                 cm_s   (ncm(1),ncm(2),ncm(3)), &
                 cm_phi0(ncm(1),ncm(2),ncm(3)), &
                 cm_phi1(ncm(1),ncm(2),ncm(3)), &
                 deltf0 (nfm(1),nfm(2),nfm(3),6), &
                 deltf1 (nfm(1),nfm(2),nfm(3),6), &
                 deltc0 (ncm(1),ncm(2),ncm(3),6), &
                 deltc1 (ncm(1),ncm(2),ncm(3),6), &
                 jsrc   (nfm(1),nfm(2),nfm(3),6), &
                 fsrc   (nfm(1),nfm(2),nfm(3),6), &
                 cmJ0   (ncm(1),ncm(2),ncm(3),6), &
                 cmJ1   (ncm(1),ncm(2),ncm(3),6), &
                 cmJn   (ncm(1),ncm(2),ncm(3),6), &
                 cmF    (ncm(1),ncm(2),ncm(3),6), &
                 cmDt   (ncm(1),ncm(2),ncm(3),6), &
                 cmDh   (ncm(1),ncm(2),ncm(3),6), &
                 Mcm    (ncm(1),ncm(2),ncm(3),7), &
                 Mfm    (nfm(1),nfm(2),nfm(3),7))
        cm_t    = 0
        cmD     = 0
        cm_a    = 0
        cm_nf   = 0
        cm_s    = 0
        cm_phi0 = 0
        cm_phi1 = 0
        deltf0  = 0
        deltf1  = 0
        deltc0  = 0
        deltc1  = 0
        jsrc    = 0
        fsrc    = 0
        cmJ0    = 0
        cmJ1    = 0
        cmJn    = 0
        cmF     = 0
        cmDt    = 0
        cmDh    = 0
        Mcm     = 0
        Mfm     = 0
    end if

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
        fm(:,:,:) % Jn(i)   = 0 
        fm(:,:,:) % J0(i)   = 0 
        fm(:,:,:) % J1(i)   = 0 
        fm(:,:,:) % sphi(i) = 0 
    enddo 

    fsd_MC = 0
        
end subroutine

! =============================================================================
! FMFD_INITIALIZE_THREAD initializes thread-wise parameters
! =============================================================================
subroutine FMFD_initialize_thread() 
    integer :: i 
    
    if ( .not. allocated(fm_thread) ) allocate(fm_thread(nfm(1),nfm(2),nfm(3)))
    
    fm_thread(:,:,:) % phi     = 0 
    fm_thread(:,:,:) % sig_t   = 0 
    fm_thread(:,:,:) % sig_a   = 0 
    fm_thread(:,:,:) % nusig_f = 0 
    do i = 1, 6
    fm_thread(:,:,:) % Jn(i)   = 0 
    fm_thread(:,:,:) % J0(i)   = 0 
    fm_thread(:,:,:) % J1(i)   = 0 
    fm_thread(:,:,:) % sphi(i) = 0
    enddo 
    
end subroutine
    

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
    real(8) :: xyz(3), uvw(3), xyz1(3)
    real(8) :: d_temp(6)
    integer :: i,j, i_coord
    integer :: idx_temp, idx_surf
    real(8) :: J_temp
    integer :: a, b, c

    ! distance initialization
    d_FMFD = INFINITY
    
    ! Find lattice index in FMFD grid
    xyz(:) = p%coord(1)%xyz(:)
    uvw(:) = p%coord(1)%uvw(:)
    i_xyz  = FMFD_ID(p%coord(1)%xyz(:))
!    print*, "ijk", i_xyz
    
    ! the particle is inside the FMFD grid
    inside_FMFD = INSIDE(xyz)
    income_FMFD = 0

    if ( inside_FMFD ) then
        d_temp(1) = ((dfm(1)*(i_xyz(1)-1)+fm0(1))-xyz(1))/uvw(1)   ! x0
        d_temp(2) = ((dfm(1)*(i_xyz(1)  )+fm0(1))-xyz(1))/uvw(1)   ! x1
        d_temp(3) = ((dfm(2)*(i_xyz(2)-1)+fm0(2))-xyz(2))/uvw(2)   ! y0
        d_temp(4) = ((dfm(2)*(i_xyz(2)  )+fm0(2))-xyz(2))/uvw(2)   ! y1
        d_temp(5) = ((dfm(3)*(i_xyz(3)-1)+fm0(3))-xyz(3))/uvw(3)   ! z0
        d_temp(6) = ((dfm(3)*(i_xyz(3)  )+fm0(3))-xyz(3))/uvw(3)   ! z1
        
        do i = 1, 6
        if ( d_temp(i) > 0 .and. d_FMFD > d_temp(i) ) then
            d_FMFD = d_temp(i)
            i_surf = i
        end if
        end do

    ! the particle is outside the FMFD grid
    else
        d_temp(1) = (fm0(1)-xyz(1))/uvw(1)   ! x0
        d_temp(2) = (fm1(1)-xyz(1))/uvw(1)   ! x1
        d_temp(3) = (fm0(2)-xyz(2))/uvw(2)   ! x0
        d_temp(4) = (fm1(2)-xyz(2))/uvw(2)   ! y1
        d_temp(5) = (fm0(3)-xyz(3))/uvw(3)   ! z0
        d_temp(6) = (fm1(3)-xyz(3))/uvw(3)   ! z1

        do i = 1, 6
        if ( d_temp(i) > 0 .and. d_FMFD > d_temp(i) ) then
            xyz1(:) = xyz(:) + d_temp(i)*uvw(:)
            select case(i)
            case(1,2)
                if ( xyz1(2) < fm0(2) .or. fm1(2) < xyz1(2) ) cycle
                if ( xyz1(3) < fm0(3) .or. fm1(3) < xyz1(3) ) cycle
            case(3,4)
                if ( xyz1(3) < fm0(3) .or. fm1(3) < xyz1(3) ) cycle
                if ( xyz1(1) < fm0(1) .or. fm1(1) < xyz1(1) ) cycle
            case(5,6)
                if ( xyz1(1) < fm0(1) .or. fm1(1) < xyz1(1) ) cycle
                if ( xyz1(2) < fm0(2) .or. fm1(2) < xyz1(2) ) cycle
            end select
            i_xyz = FMFD_ID(xyz1)
            d_FMFD = d_temp(i)
            income_FMFD = i
        end if
        end do
    end if

end subroutine

! =============================================================================
! FMFD_ID finds the x, y, z indice in the FMFD mesh grid
! =============================================================================
function FMFD_ID(fmxyz) result(fmid)
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: fmid(3)               ! indice
       
    fmid(:) = floor((fmxyz(:)-fm0(:))/dfm(:))+1

end function

! =============================================================================
! FMFD_ID finds the x, y, z indice in the FMFD mesh grid
! =============================================================================
function INSIDE(fmxyz) result(inside_FMFD)
    real(8), intent(in):: fmxyz(:)  ! coordinate
    integer:: inside_FMFD
    integer:: ii
       
    inside_FMFD = .true.
    do ii = 1, 3
        if ( fmxyz(ii) < fm0(ii) .or. fmxyz(ii) >= fm1(ii) ) then
            inside_FMFD = .false.
            exit
        end if
    end do

end function

! =============================================================================
! ISINCOMING determines if the particle is coming to the FMFD grid
! =============================================================================
subroutine INCOMING(xyz,income)
    integer, intent(in):: xyz(:)
    integer, intent(out):: income

    ! x0
    if ( xyz(1) == 0 ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
            income = 2
        end if
        end if
    end if
    ! x1
    if ( xyz(1) == nfm(1)+1 ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
            income = 1
        end if
        end if
    end if
    ! y0
    if ( xyz(2) == 0 ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
            income = 4
        end if
        end if
    end if
    ! y1
    if ( xyz(2) == nfm(2)+1 ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
            income = 3
        end if
        end if
    end if
    ! z0
    if ( xyz(3) == 0 ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
            income = 6
        end if
        end if
    end if
    ! z1
    if ( xyz(3) == nfm(3)+1 ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
            income = 5
        end if
        end if
    end if

end subroutine

! =============================================================================
! FMFD_TRK calculates FMFD parameters such as flux, group contstans by 
! track-length estiamtor
! =============================================================================
subroutine FMFD_TRK(wgt,distance,macro_xs,id,idmc,mcxyz)
    implicit none
    type(Particle):: p
    real(8), intent(in) :: wgt
    real(8), intent(in) :: distance
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(1:3)
    integer, intent(in) :: idmc(1:3)
    real(8), intent(in) :: mcxyz(1:3)
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

!    if ( id(1) == 34 .and. id(2) == 34 ) then
!        print*, id
!        print*, fm_thread(id(1),id(2),id(3))%phi
!        pause
!    end if
    
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
subroutine FMFD_SURF (inside,income, is, id, uvw, wgt, bc)
    logical, intent(in) :: inside
    integer, intent(in) :: income
    integer, intent(in) :: is, id(3)
    real(8), intent(in) :: uvw(3)
    real(8), intent(in) :: wgt
    integer, intent(in) :: bc
    integer:: dir
    
    ! inner nodes
    if ( inside ) then 
        ! surface partial current
        select case(is)
        case(1,3,5)
            fm_thread(id(1),id(2),id(3))%J0(is) = &
            fm_thread(id(1),id(2),id(3))%J0(is) + wgt
            if ( cmfdon ) &
            fm_thread(id(1),id(2),id(3))%sphi(is) = &
            fm_thread(id(1),id(2),id(3))%sphi(is) + wgt/abs(uvw((is+1)/2))
        case(2,4,6)
            fm_thread(id(1),id(2),id(3))%J1(is) = &
            fm_thread(id(1),id(2),id(3))%J1(is) + wgt
            if ( cmfdon ) &
            fm_thread(id(1),id(2),id(3))%sphi(is) = &
            fm_thread(id(1),id(2),id(3))%sphi(is) + wgt/abs(uvw(is/2))
        end select

        ! boundary condition
        if ( bc == 2 ) then
        select case(is)
        case(1,3,5)
            fm_thread(id(1),id(2),id(3))%J1(is) = &
            fm_thread(id(1),id(2),id(3))%J1(is) + wgt
            if ( cmfdon ) &
            fm_thread(id(1),id(2),id(3))%sphi(is) = &
            fm_thread(id(1),id(2),id(3))%sphi(is) + wgt/abs(uvw((is+1)/2))
        case(2,4,6)
            fm_thread(id(1),id(2),id(3))%J0(is) = &
            fm_thread(id(1),id(2),id(3))%J0(is) + wgt
            if ( cmfdon ) &
            fm_thread(id(1),id(2),id(3))%sphi(is) = &
            fm_thread(id(1),id(2),id(3))%sphi(is) + wgt/abs(uvw((is+1)/2))
        end select
        end if
        return
    end if

    ! boundary nodes
    select case(income)
    case(1)
        fm_thread(1,id(2),id(3))%J1(1) = &
        fm_thread(1,id(2),id(3))%J1(1) + wgt
        if ( cmfdon ) &
        fm_thread(1,id(2),id(3))%sphi(1) = &
        fm_thread(1,id(2),id(3))%sphi(1) + wgt/abs(uvw(1))
    case(2)
        fm_thread(nfm(1),id(2),id(3))%J0(2) = &
        fm_thread(nfm(1),id(2),id(3))%J0(2) + wgt
        if ( cmfdon ) &
        fm_thread(nfm(1),id(2),id(3))%sphi(2) = &
        fm_thread(nfm(1),id(2),id(3))%sphi(2) + wgt/abs(uvw(1))
    case(3)
        fm_thread(id(1),1,id(3))%J1(3) = &
        fm_thread(id(1),1,id(3))%J1(3) + wgt
        if ( cmfdon ) &
        fm_thread(id(1),1,id(3))%sphi(3) = &
        fm_thread(id(1),1,id(3))%sphi(3) + wgt/abs(uvw(2))
    case(4)
        fm_thread(id(1),nfm(2),id(3))%J0(4) = &
        fm_thread(id(1),nfm(2),id(3))%J0(4) + wgt
        if ( cmfdon ) &
        fm_thread(id(1),nfm(2),id(3))%sphi(4) = &
        fm_thread(id(1),nfm(2),id(3))%sphi(4) + wgt/abs(uvw(2))
    case(5)
        fm_thread(id(1),id(2),1)%J1(5) = &
        fm_thread(id(1),id(2),1)%J1(5) + wgt
        if ( cmfdon ) &
        fm_thread(id(1),id(2),1)%sphi(5) = &
        fm_thread(id(1),id(2),1)%sphi(5) + wgt/abs(uvw(3))
    case(6)
        fm_thread(id(1),id(2),nfm(3))%J0(6) = &
        fm_thread(id(1),id(2),nfm(3))%J0(6) + wgt
        if ( cmfdon ) &
        fm_thread(id(1),id(2),nfm(3))%sphi(6) = &
        fm_thread(id(1),id(2),nfm(3))%sphi(6) + wgt/abs(uvw(3))
    end select
            
end subroutine

! =============================================================================
! NORM_FMFD normalizes cycle-wise FMFD parameters
! =============================================================================
subroutine NORM_FMFD()
    implicit none

    !> gather thread FMFD parameters
    fm(:,:,:) % phi     = fm(:,:,:) % phi     + fm_thread(:,:,:)%phi
    fm(:,:,:) % sig_t   = fm(:,:,:) % sig_t   + fm_thread(:,:,:)%sig_t 
    fm(:,:,:) % sig_a   = fm(:,:,:) % sig_a   + fm_thread(:,:,:)%sig_a 
    fm(:,:,:) % nusig_f = fm(:,:,:) % nusig_f + fm_thread(:,:,:)%nusig_f 
    do ii = 1, 6
    fm(:,:,:) % J0(ii)   = fm(:,:,:) % J0(ii)   + fm_thread(:,:,:)%J0(ii)
    fm(:,:,:) % J1(ii)   = fm(:,:,:) % J1(ii)   + fm_thread(:,:,:)%J1(ii)
    fm(:,:,:) % sphi(ii) = fm(:,:,:) % sphi(ii) + fm_thread(:,:,:)%sphi(ii)
    end do

end subroutine


! =============================================================================
! PROCESS_FMFD deals with MPI process and average quantities
! =============================================================================
subroutine PROCESS_FMFD() 
    !> MPI derived type reduce parameters 
    integer, dimension(8) :: FMFD_blocklength, FMFD_displacement, FMFD_datatype 
    integer :: intex, realex, restype, FMFD_op  ! MPI int & real extern / new type
    type(FMFD_parameters), allocatable :: FMFD_MPI_slut(:,:,:)
    integer :: i, j, k, l
    real(8) :: aa, bb

    data FMFD_blocklength /1, 1, 1, 1, 6, 6, 6, 6 /

    allocate(FMFD_MPI_slut(nfm(1),nfm(2),nfm(3)))
    
    !> Gather FMFD parameters from the slave nodes 
    call MPI_TYPE_EXTENT(MPI_REAL8, realex, ierr) 
    FMFD_displacement(1) = 0
    FMFD_displacement(2) = realex
    FMFD_displacement(3) = 2*realex
    FMFD_displacement(4) = 3*realex
    FMFD_displacement(5) = 4*realex
    FMFD_displacement(6) = (4+6)*realex
    FMFD_displacement(7) = (4+2*6)*realex
    FMFD_displacement(8) = (4+3*6)*realex
    
    FMFD_datatype(:) = MPI_REAL8
    
    call MPI_TYPE_STRUCT (8, FMFD_blocklength, FMFD_displacement, FMFD_datatype,restype, ierr)
    call MPI_TYPE_COMMIT (restype, ierr)
    call MPI_TYPE_COMMIT (restype, ierr)
    call MPI_Op_create(FMFD_SUM, .true. , FMFD_op, ierr)
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
    call MPI_REDUCE(fm(i,j,k), FMFD_MPI_slut(i,j,k), 1, restype, FMFD_op, score, MPI_COMM_WORLD, ierr)
    end do
    end do
    end do
    fm = FMFD_MPI_slut
     
    call MPI_TYPE_FREE(restype, ierr)
    deallocate(FMFD_MPI_slut)
    call MPI_Op_Free(FMFD_Op, ierr)
    
    if ( icore /= score ) return 

    ! current swapping
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 )       fm(i,j,k)%J1(1) = fm(i-1,j,k)%J1(2)
        if ( i /= nfm(1) )  fm(i,j,k)%J0(2) = fm(i+1,j,k)%J0(1)
        if ( j /= 1 )       fm(i,j,k)%J1(3) = fm(i,j-1,k)%J1(4)
        if ( j /= nfm(2) )  fm(i,j,k)%J0(4) = fm(i,j+1,k)%J0(3)
        if ( k /= 1 )       fm(i,j,k)%J1(5) = fm(i,j,k-1)%J1(6)
        if ( k /= nfm(3) )  fm(i,j,k)%J0(6) = fm(i,j,k+1)%J0(5)
    end do
    end do
    end do

    ! surface flux swapping
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 ) then
        fm(i,j,k)%sphi(1) = fm(i,j,k)%sphi(1) + fm(i-1,j,k)%sphi(2)
        fm(i-1,j,k)%sphi(2) = fm(i,j,k)%sphi(1)
        end if
        if ( j /= 1 ) then
        fm(i,j,k)%sphi(3) = fm(i,j,k)%sphi(3) + fm(i,j-1,k)%sphi(4)
        fm(i,j-1,k)%sphi(4) = fm(i,j,k)%sphi(3)
        end if
        if ( k /= 1 ) then
        fm(i,j,k)%sphi(5) = fm(i,j,k)%sphi(5) + fm(i,j,k-1)%sphi(6)
        fm(i,j,k-1)%sphi(6) = fm(i,j,k)%sphi(5)
        end if
    end do
    end do
    end do

    ! group constant
    fm(:,:,:) % sig_t   = fm(:,:,:) % sig_t   / fm(:,:,:) % phi
    fm(:,:,:) % sig_a   = fm(:,:,:) % sig_a   / fm(:,:,:) % phi
    fm(:,:,:) % nusig_f = fm(:,:,:) % nusig_f / fm(:,:,:) % phi
    fm(:,:,:) % phi     = fm(:,:,:) % phi     / (dble(ngen)*v_fm)

    ! surface quantity normalization
    do i = 1, 6
    bb = dble(ngen)*a_fm(i)
    fm(:,:,:) % J0(i)   = fm(:,:,:) % J0(i) / bb
    fm(:,:,:) % J1(i)   = fm(:,:,:) % J1(i) / bb
    fm(:,:,:) % sphi(i) = fm(:,:,:) % sphi(i) / bb
    end do

    ! net current
    do i = 1, 6
    fm(:,:,:) % Jn(i) = fm(:,:,:) % J1(i) - fm(:,:,:) % J0(i)
    end do

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
      fm_avg(:,:,:)%sphi(i) = 0
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

       fm_avg(i,j,k)%Jn(:)   = fm_avg(i,j,k)%Jn(:)   + acc(l)%fm(i,j,k)%Jn(:)
       fm_avg(i,j,k)%J0(:)   = fm_avg(i,j,k)%J0(:)   + acc(l)%fm(i,j,k)%J0(:)
       fm_avg(i,j,k)%J1(:)   = fm_avg(i,j,k)%J1(:)   + acc(l)%fm(i,j,k)%J1(:)
       fm_avg(i,j,k)%sphi(:) = fm_avg(i,j,k)%sphi(:) + acc(l)%fm(i,j,k)%sphi(:)
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
    fm_avg(:,:,:)%Jn(i)   = fm_avg(:,:,:)%Jn(i)   / dble(n_acc)
    fm_avg(:,:,:)%J0(i)   = fm_avg(:,:,:)%J0(i)   / dble(n_acc)
    fm_avg(:,:,:)%J1(i)   = fm_avg(:,:,:)%J1(i)   / dble(n_acc)
    fm_avg(:,:,:)%sphi(i) = fm_avg(:,:,:)%sphi(i) / dble(n_acc)
    enddo

end subroutine 


! =============================================================================
! FMFD_SOLVE solves FMFD eigenvalue problem
! =============================================================================
subroutine FMFD_SOLVE(keff,fsd)
    use CMFD, only: ONE_NODE_CMFD
    implicit none
    real(8), intent(in):: keff              ! multiplication factor
    real(8), intent(inout) :: fsd(:,:,:)    ! fission source distribution
    real(8) :: M(nfm(1),nfm(2),nfm(3),7)    ! FMFD matrix
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: &
        sig_t, &      ! total
        sig_a, &      ! absorption
        nusig_f, &    ! nu X fission
        D, &          ! diffusion coefficient
        phi0, &       ! neutron flux (previous)
        phi1, &       ! neutron flux (current)
        fsd_FM        ! fission source distribution for FMFD
    real(8), dimension(nfm(1),nfm(2),nfm(3),6):: &
        D_tilda, &    ! D tilda
        D_hat, &      ! correction factor
        Jn, &         ! net current
        J0, &         ! partial current -
        J1, &         ! partial current +
        sphi          ! surface flux
    integer :: i, j, k

    ! copy parameters
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        sig_t(i,j,k)   = fm_avg(i,j,k)%sig_t
        sig_a(i,j,k)   = fm_avg(i,j,k)%sig_a
        nusig_f(i,j,k) = fm_avg(i,j,k)%nusig_f 
        D(i,j,k)       = 1D0 / (3D0 * fm_avg(i,j,k)%sig_t) 
        Jn(i,j,k,:)    = fm_avg(i,j,k)%Jn(:) 
        J0(i,j,k,:)    = fm_avg(i,j,k)%J0(:) 
        J1(i,j,k,:)    = fm_avg(i,j,k)%J1(:) 
        phi1(i,j,k)    = fm_avg(i,j,k)%phi
        sphi(i,j,k,:)  = fm_avg(i,j,k)%sphi(:)
    enddo 
    enddo
    enddo


    if ( cmfdon ) then
        do i = 1, 6
            sphi(:,:,:,1) = 2D0*J1(:,:,:,1)+2D0*J0(:,:,:,1)
            sphi(:,:,:,2) = 2D0*J1(:,:,:,2)+2D0*J0(:,:,:,2)
            sphi(:,:,:,3) = 2D0*J1(:,:,:,3)+2D0*J0(:,:,:,3)
            sphi(:,:,:,4) = 2D0*J1(:,:,:,4)+2D0*J0(:,:,:,4)
            sphi(:,:,:,5) = 2D0*J1(:,:,:,5)+2D0*J0(:,:,:,5)
            sphi(:,:,:,6) = 2D0*J1(:,:,:,6)+2D0*J0(:,:,:,6)
        end do

        call ONE_NODE_CMFD(keff,sig_t,sig_a,nusig_f,D,phi1,J0,J1,Jn,sphi)

    else
        !> calculate D_tilda & D_hat 
        call D_TILDA_CALCULATION(D,D_tilda)
        call D_HAT_CALCULATION(D_tilda,Jn,phi1,D_hat)

    
        ! matrix composition
        if ( FMFD_type == 1 ) then 
            !> CMFD 
            M = getM(D_tilda, D_hat, sig_a)
        else!if (CMFD_type == 2) then 
            !> p-CMFD
            !M = getMp (phi,D,J_pp,J_pn,sig_a,a,b,c,pitch(1),pitch(2),pitch(3))
        endif 
         
        k_fmfd = keff
        call POWER (k_fmfd, M, phi0, phi1, nusig_f)
    end if

    !print *, 'CMFD keff  ',keff_CMFD
    !write(prt_keff,*) keff_CMFD, k_col, k_tl

!    do jj = 34, 1, -1
!    write(8,1), (phi1(ii,jj,1)/sum(phi1), ii = 1, 34)
!    end do
!    write(8,1), Jn(:,:,1,1)/sum(phi1(:,:,1))
!    write(8,*)
!    write(8,1), Jn(:,:,1,2)/sum(phi1(:,:,1))
!    write(8,*)
!    write(8,1), Jn(:,:,1,3)/sum(phi1(:,:,1))
!    write(8,*)
!    write(8,1), Jn(:,:,1,4)/sum(phi1(:,:,1))
!    write(8,*)
!    write(8,1), Jn(:,:,1,5)/phi1(:,:,1)
!    write(8,*)
!    write(8,1), Jn(:,:,1,6)/phi1(:,:,1)
!    write(8,*)
!    1 format(34es15.7)
    
    !> CMFD feedback (modulation)
    fsd_MC(:,:,:) = fsd_MC(:,:,:) / sum(fsd_MC)
    fsd_FM(:,:,:) = nusig_f(:,:,:)*phi1(:,:,:)
    fsd_FM(:,:,:) = fsd_FM(:,:,:) / sum(fsd_FM)
    where ( fsd_MC(:,:,:) /= 0 ) &
    fsd(:,:,:)    = fsd_FM(:,:,:) / fsd_MC(:,:,:)

    print*, "k_fmfd", k_fmfd
!    write(8,*)
!    do j = 1, nfm(2)
!        write(8,1), (phi1(i,j,1), i = 1, nfm(1))
!    end do
!    write(8,*)
!    do j = 1, nfm(2)
!        write(8,1), (fsd(i,j,1), i = 1, nfm(1))
!    end do
!    write(8,*)
!    1 format(34es15.7)
    stop
    
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
        Dt(i,j,k,6) = 2D0*D(i,j,k+1)*D(i,j,k)/(D(i,j,k+1)+D(i,j,k))/dfm(3)
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
    do k = 1, nfm(3)
    do j = 1, nfm(2)
    do i = 1, nfm(1)
        if ( i /= 1 )      M(i,j,k,3) = -(Dt(i,j,k,1)+Dh(i,j,k,1))/dfm(1) ! x0
        if ( i /= nfm(1) ) M(i,j,k,5) = -(Dt(i,j,k,2)-Dh(i,j,k,2))/dfm(1) ! x1
        if ( j /= 1 )      M(i,j,k,2) = -(Dt(i,j,k,3)+Dh(i,j,k,3))/dfm(2) ! y0
        if ( j /= nfm(2) ) M(i,j,k,6) = -(Dt(i,j,k,4)-Dh(i,j,k,4))/dfm(2) ! y1
        if ( k /= 1 )      M(i,j,k,1) = -(Dt(i,j,k,5)+Dh(i,j,k,5))/dfm(3) ! z0
        if ( k /= nfm(3) ) M(i,j,k,7) = -(Dt(i,j,k,6)-Dh(i,j,k,6))/dfm(3) ! z1
        
        M(i,j,k,4) = &
            +(Dt(i,j,k,1)-Dh(i,j,k,1)+Dt(i,j,k,2)+Dh(i,j,k,2))/dfm(1) &
            +(Dt(i,j,k,3)-Dh(i,j,k,3)+Dt(i,j,k,4)+Dh(i,j,k,4))/dfm(2) &
            +(Dt(i,j,k,5)-Dh(i,j,k,5)+Dt(i,j,k,6)+Dh(i,j,k,6))/dfm(3) &
            +sig_a(i,j,k)

    enddo
    enddo
    enddo

end function getM


subroutine POWER (k_eff, M, phi0, phi1, nusig_f)
    use SOLVERS, only: BICGStab_hepta
    implicit none
    real(8), intent(inout):: k_eff
    real(8), intent(in) :: M(:,:,:,:), nusig_f(:,:,:)
    real(8), intent(inout):: phi0(:,:,:), phi1(:,:,:)
    real(8), parameter:: ONE = 1D0
    real(8) :: F(1:nfm(1),1:nfm(2),1:nfm(3))
    integer :: iter_max = 1E5
    integer :: iter = 0
    real(8) :: err
    
    err = ONE
    do while ( ( err > 1.0d-10 ) .and. (iter < iter_max) )
        iter = iter + 1
        phi0 = phi1
        F = nusig_f(:,:,:)*phi1(:,:,:)/k_eff
        !call SOR(M(:,:,:,:),F(:,:,:),phi1(:,:,:))
        phi1 = BiCGStab_hepta(M(:,:,:,:),F(:,:,:))
        k_eff = k_eff*sum(nusig_f*phi1*nusig_f*phi1) &
               / sum(nusig_f*phi1*nusig_f*phi0)
        err = maxval(abs(phi1(:,:,:)-phi0(:,:,:))/phi1(:,:,:))
    enddo
    
end subroutine 


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



subroutine FMFD_SUM (inpar, inoutpar, partype) 
    type(FMFD_parameters), intent(in) :: inpar(nfm(1),nfm(2),nfm(3))
    type(FMFD_parameters), intent(inout) :: inoutpar(nfm(1),nfm(2),nfm(3))
    integer :: ii, partype
    
    inoutpar(:,:,:)%phi     = inoutpar(:,:,:)%phi     + inpar(:,:,:)%phi
    inoutpar(:,:,:)%sig_t   = inoutpar(:,:,:)%sig_t   + inpar(:,:,:)%sig_t
    inoutpar(:,:,:)%sig_a   = inoutpar(:,:,:)%sig_a   + inpar(:,:,:)%sig_a
    inoutpar(:,:,:)%nusig_f = inoutpar(:,:,:)%nusig_f + inpar(:,:,:)%nusig_f
    do ii = 1, 6 
    inoutpar(:,:,:)%Jn(ii)   = inoutpar(:,:,:)%Jn(ii)   + inpar(:,:,:)%Jn(ii)
    inoutpar(:,:,:)%J0(ii)   = inoutpar(:,:,:)%J0(ii)   + inpar(:,:,:)%J0(ii)
    inoutpar(:,:,:)%J1(ii)   = inoutpar(:,:,:)%J1(ii)   + inpar(:,:,:)%J1(ii)
    inoutpar(:,:,:)%sphi(ii) = inoutpar(:,:,:)%sphi(ii) + inpar(:,:,:)%sphi(ii)
    enddo 

end subroutine  

    
    
end module     
