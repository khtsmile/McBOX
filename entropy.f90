module ENTROPY
    use GEOMETRY_HEADER
    use CONSTANTS
    use PARTICLE_HEADER, only: Particle
    implicit none
    integer:: ii
    real(8):: entrp
    integer:: n_ben                     ! No of Cell entropy
    type entropy_bins
        integer:: n_sub
        integer:: id                    ! lattice ID
        real(8), allocatable:: entrp(:,:,:)
    end type
    type(entropy_bins), allocatable:: bens(:)

contains

! =============================================================================
! ENTRP_INIT initializes parameters and options for Shannon entropy calcluation
! =============================================================================
subroutine ENTRP_INIT
    implicit none
    type(Cell), pointer:: c
    type(lattice), pointer:: l
    integer, allocatable:: temp_bens(:) ! temporary bins for entropy

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

end subroutine


! =============================================================================
! SHENTROPY calculates Shannon entropy at every cycle
! =============================================================================
subroutine SHENTROPY(sb)
    use BANK_HEADER,    only: Bank
    use GEOMETRY,       only: FIND_CELL
    implicit none
    type(Bank), intent(in):: sb(:)
    type(Particle):: p
    integer:: exyz(0:3)
    logical:: found

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

end function

end module
