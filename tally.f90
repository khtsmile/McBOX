module tally
    use geometry_header 
    use constants
    use particle_header, only : Particle
    
    implicit none
    
    
    type :: LocalCoord
        integer :: cell      = NONE
        integer :: universe  = NONE
        integer :: lattice   = NONE
        integer :: lattice_x = NONE
        integer :: lattice_y = NONE
        integer :: lattice_z = NONE
      contains 
        procedure :: reset => reset_coord
    end type LocalCoord

    type :: CoordStruct
        integer          :: n_coord
        type(LocalCoord) :: coord(MAX_COORD) 
        real(8)          :: vol
        integer          :: flag ! 0 if pin 
                                 ! 1 if cell
    end type
    
    !> Variables
    type(CoordStruct), allocatable, target :: TallyCoord(:)
    real(8), allocatable :: TallyFlux(:), TallyPower(:)
    real(8), allocatable :: tally1(:), tally2(:)
    
    
    contains
    
    elemental subroutine reset_coord(this)
        class(LocalCoord), intent(inout) :: this
        
        this % cell     = NONE
        this % universe = NONE
        this % lattice     = NONE
        this % lattice_x = NONE
        this % lattice_y = NONE
        this % lattice_z = NONE
        
    end subroutine reset_coord


    function FindTallyBin(p) result (idx)
        type(particle), intent(in) :: p 
        integer :: idx(4)
        integer :: i_bin, i_coord, i, j
        integer, dimension(6) :: A, B, C 
        
        idx = -1
        Bin: do i_bin = 1, size(TallyCoord)
            !if (p%n_coord /= TallyCoord(i_bin)%n_coord) cycle Bin
            do i_coord = 1, TallyCoord(i_bin)%n_coord-1! p%n_coord-1
                
                A(1) = p%coord(i_coord)%cell
                A(2) = p%coord(i_coord)%universe 
                A(3) = p%coord(i_coord)%lattice
                A(4) = p%coord(i_coord)%lattice_x
                A(5) = p%coord(i_coord)%lattice_y
                A(6) = p%coord(i_coord)%lattice_z
                
                B(1) = TallyCoord(i_bin)%coord(i_coord)%cell
                B(2) = TallyCoord(i_bin)%coord(i_coord)%universe
                B(3) = TallyCoord(i_bin)%coord(i_coord)%lattice
                B(4) = TallyCoord(i_bin)%coord(i_coord)%lattice_x
                B(5) = TallyCoord(i_bin)%coord(i_coord)%lattice_y
                B(6) = TallyCoord(i_bin)%coord(i_coord)%lattice_z
                
                C(:) = abs(A(:)-B(:))
                if (sum(C) /= 0) cycle Bin
                
            enddo 
            i_coord = TallyCoord(i_bin)%n_coord !p%n_coord
            A(1) = p%coord(i_coord)%cell     * TallyCoord(i_bin)%flag
            A(2) = p%coord(i_coord)%universe * TallyCoord(i_bin)%flag
            A(3) = p%coord(i_coord)%lattice
            A(4) = p%coord(i_coord)%lattice_x
            A(5) = p%coord(i_coord)%lattice_y
            A(6) = p%coord(i_coord)%lattice_z
            
            B(1) = TallyCoord(i_bin)%coord(i_coord)%cell
            B(2) = TallyCoord(i_bin)%coord(i_coord)%universe
            B(3) = TallyCoord(i_bin)%coord(i_coord)%lattice
            B(4) = TallyCoord(i_bin)%coord(i_coord)%lattice_x
            B(5) = TallyCoord(i_bin)%coord(i_coord)%lattice_y
            B(6) = TallyCoord(i_bin)%coord(i_coord)%lattice_z
            
            C(:) = abs(A(:)-B(:))
            if (sum(C) /= 0) cycle Bin
            idx(1) = i_bin
            idx(2:4) = B(4:6)
            exit Bin 

        enddo Bin

    end function 

end module 
