module ENTROPY
    use variables
    use GEOMETRY_HEADER
    use CONSTANTS
    use PARTICLE_HEADER, only: Particle
    implicit none
    integer:: ii
    real(8):: entrp     ! entropy at this cycle
    real(8):: entrp1    ! MPRUP entropy at previous generation
    real(8):: entrp2    ! MPRUP entropy at current generation
    real(8):: dentrp    ! entropy difference
    integer:: n_ben                     ! No of Cell entropy
    type entropy_bins
        integer:: n_sub
        integer:: id                    ! lattice ID
        real(8), allocatable:: entrp(:,:,:)
    end type
    type(entropy_bins), allocatable:: bens(:)
    !real(8), allocatable:: entrp(:,:,:)

    ! =========================================================================
    ! Modified paricle ramp-up method (MPRUP)
    logical:: mprupon = .false.         ! MPRUP on ?
    logical:: genup   = .false.         ! generation size up ?
    integer:: rampup                    ! rampup generation size
    real(8), allocatable:: shannon(:)   ! cycle-wise shannon
    real(8):: dshannon                  ! averaged difference
    integer:: elength                   ! accumulation length
    integer:: ccrt, scrt                ! type of criteria, judged by
                                        ! 1 - relative error
                                        ! 2 - no. of cycles / histories
    real(8):: crt1, crt2, crt3          ! convergence / stopping / finishing
    real(8):: en0(3)                    ! (x0,y0,z0)
    real(8):: en1(3)                    ! (x1,y1,z1)
    real(8):: en2(3)                    ! (x2,y2,z2)
    integer:: nen(3)                    ! (nx,ny,nz)
    real(8):: den(3)                    ! (dx,dy,dz)
    logical:: entrp_grid = .false.

contains

! =============================================================================
! ENTRP_INIT initializes parameters and options for Shannon entropy calcluation
! =============================================================================
subroutine ENTRP_INIT
    implicit none
    type(Cell), pointer:: c
    type(lattice), pointer:: l
    integer, allocatable:: temp_bens(:) ! temporary bins for entropy

    ! -------------------------------------------------------------------------
    ! user defined entropy mesh grid
    if ( entrp_grid ) then
        en2(:) = en1(:) - en0(:)
        den(:) = en2(:) / nen(:)

        n_ben = 1
        allocate(bens(1))
        allocate(bens(1)%entrp(nen(1),nen(2),nen(3)))
        bens(1)%entrp(:,:,:) = 0D0

!        allocate(entrp(nen(1),nen(2),nen(3)))
!        entrp = 0D0

    ! -------------------------------------------------------------------------
    ! assembly size entropy mesh grid
    else
        ! search assembly lattice
        allocate(temp_bens(100))
        n_ben = 0
        do ii = 1, size(cells)
            c => cells(ii)
            if ( c%filltype == FILL_LATTICE .and. c%univ_id == BASE ) then
                n_ben = n_ben + 1
                temp_bens(n_ben) = find_lat_idx(lattices,c%fill)
            end if
        end do
        nullify(c)
    
        ! allocate entorpy bins
        allocate(bens(n_ben))
        bens(1:n_ben)%id = temp_bens(1:n_ben)
        do ii = 1, n_ben
            l => lattices(bens(ii)%id)
            allocate(bens(ii)%entrp(l%n_xyz(1),l%n_xyz(2),l%n_xyz(3)))
            bens(ii)%entrp(:,:,:) = 0D0
        end do
        nullify(l)
    end if

end subroutine


! =============================================================================
! SHENTROPY calculates Shannon entropy at every cycle
! =============================================================================
subroutine SHENTROPY(sb)
    use BANK_HEADER,        only: Bank
    use GEOMETRY,           only: FIND_CELL
    implicit none
    type(Bank), intent(in):: sb(:)
    type(Particle):: p
    integer:: exyz(0:3)
    logical:: found
    integer:: ii, jj

    ! initialization
    do ii = 1, n_ben
        bens(ii)%entrp(:,:,:) = 0
    end do

    ! local entropy
    do ii = 1, size(sb)
        call p%initialize()
        call p%set(sb(ii))
        call FIND_CELL(p,found)
        exyz(:) = FIND_ENTRP_LAT(p)
        bens(exyz(0))%entrp(exyz(1),exyz(2),exyz(3)) = &
        bens(exyz(0))%entrp(exyz(1),exyz(2),exyz(3)) + p%wgt
    end do

    ! global entropy
    entrp = 0
    do ii = 1, n_ben
        entrp = entrp + sum(bens(ii)%entrp(:,:,:))
    end do
    do ii = 1, n_ben
        bens(ii)%entrp(:,:,:) = bens(ii)%entrp(:,:,:) / entrp
    end do
    do ii = 1, n_ben
        where ( bens(ii)%entrp /= 0 ) &
            bens(ii)%entrp = -bens(ii)%entrp*log(bens(ii)%entrp)/log(2D0)
    end do
    entrp = 0
    do ii = 1, n_ben
        entrp = entrp + sum(bens(ii)%entrp)
    end do

    ! MPRUP method
    if ( mprupon ) then
        shannon = eoshift(shannon,-1,entrp,1)
        dshannon = abs(sum(shannon(1:elength)-shannon(elength+1:2*elength))) &
                 / sum(shannon(1:elength))
    end if

end subroutine

! =============================================================================
! FIND_ENTRP_LAT
! =============================================================================
function FIND_ENTRP_LAT(p) result(exyz)
    implicit none
    type(Particle), intent(in):: p
    integer:: exyz(0:3)   ! 1x 2y 3z 4i (i=cell index)
    integer:: i_ben
    integer:: ncord
    integer:: icord
    real(8):: xyz(3)


    if ( entrp_grid ) then
        xyz(:) = p%coord(1)%xyz(:)
        exyz(0) = 1
        exyz(1:3) = floor((xyz(1:3)-en0(1:3))/den(1:3))+1

    else
        ! search lattice index
        do i_ben = 1, n_ben
        do icord = 1, p%n_coord
            if ( p%coord(icord)%lattice == bens(i_ben)%id ) then
                exyz(0) = i_ben
                exyz(1) = p%coord(icord)%lattice_x
                exyz(2) = p%coord(icord)%lattice_y
                exyz(3) = p%coord(icord)%lattice_z
                return
            end if
        end do
        end do
   end if

end function

end module
