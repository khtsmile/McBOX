module material_header 
	use omp_lib
	
	implicit none 

	type Material_CE
		character(len=20)    :: mat_name        ! User-defined name
		integer              :: n_iso = 0  		! number of isotopes (nuclides)
		integer, allocatable :: ace_idx(:)      ! index in nuclides array
		real(8), allocatable :: numden(:) 		! nuclide atom density in atom/b-cm
		real(8)              :: density_gpcc    ! total density in g/cm^3
		
		! Does this material contain fissionable nuclides? Is it depletable?
		logical :: fissionable = .false.
		logical :: depletable = .false.
        ! Isotopes for S(a,b)?
        logical:: sab = .false.
		
	end type Material_CE
	
	integer :: n_materials ! # of materials
	
	type(Material_CE), allocatable, target :: materials(:), materials_temp(:)
	type(Material_CE), pointer :: CE_mat_ptr

	
	contains
	
	function find_CE_mat_idx (this, mat_id) result (idx) 
		type(Material_CE) :: this(:) 
		character(*) :: mat_id 
		integer :: i, idx, n 
		
		n = size(this)
		do i = 1, n
			if (trim(this(i)%mat_name) == mat_id) then 
				idx = i 
				return 
			endif
		enddo 
		print *, "no such CE mat id : ", mat_id 
		stop 
		
	end function 

	
	
	
	
	
end module 
