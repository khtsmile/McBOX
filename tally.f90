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
	end type
	
	type(CoordStruct), allocatable :: TallyCoord(:)
	real(8), 		   allocatable :: TallyFlux(:)
	
	
	contains
	
	elemental subroutine reset_coord(this)
		class(LocalCoord), intent(inout) :: this
		
		this % cell 	= NONE
		this % universe = NONE
		this % lattice 	= NONE
		this % lattice_x = NONE
		this % lattice_y = NONE
		this % lattice_z = NONE
		
	end subroutine reset_coord


	function FindTallyBin(p) result (idx)
		type(particle), intent(in) :: p 
		integer :: idx
		integer :: i_bin, i_coord, i, j
		
		idx = -1
		Bin: do i_bin = 1, size(TallyCoord)
			if (p%n_coord /= TallyCoord(i_bin)%n_coord) cycle Bin
			
			do i_coord = 1, p%n_coord
				if (p%coord(i_coord)%cell /= TallyCoord(i_bin)%coord(i_coord)%cell) cycle Bin
				if (p%coord(i_coord)%universe /= TallyCoord(i_bin)%coord(i_coord)%universe) cycle Bin
				if (p%coord(i_coord)%lattice /= TallyCoord(i_bin)%coord(i_coord)%lattice) cycle Bin
				if (p%coord(i_coord)%lattice_x /= TallyCoord(i_bin)%coord(i_coord)%lattice_x) cycle Bin
				if (p%coord(i_coord)%lattice_y /= TallyCoord(i_bin)%coord(i_coord)%lattice_y) cycle Bin
				if (p%coord(i_coord)%lattice_z /= TallyCoord(i_bin)%coord(i_coord)%lattice_z) cycle Bin
			enddo 
			
			idx = i_bin; exit
		enddo Bin
		
		
	end function 


end module 