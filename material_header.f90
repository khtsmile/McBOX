module material_header 
    
    implicit none 

    type Material_CE
        character(len=20)    :: mat_name       ! User-defined name
        integer              :: n_iso = 0      ! number of isotopes (nuclides)
        integer, allocatable :: ace_idx(:)     ! index in nuclides array
        real(8), allocatable :: numden(:)      ! nuclide atom density (#/b-cm)
        real(8), allocatable :: full_numden(:) ! AD for all inventory (#/b-cm)
        real(8), allocatable :: temp(:)        ! temperature (MeV)
        real(8)              :: density_gpcc ! total density in g/cm^3
        real(8)              :: vol          ! material volume (cm^3)
        
        ! Does this material contain fissionable nuclides? Is it depletable?
        logical :: fissionable = .false.
        logical :: depletable = .false.
        real(8) :: flux = 0.0d0
        real(8) :: pwr  = 0.0d0
        real(8), allocatable :: ogxs(:,:)

        !ogxs(:,1) = One-group Volume-integrated (n,g)  // mt 102 (ENDF)
        !ogxs(:,2) = One-group Volume-integrated (n,2n) // mt 16 (ENDF)
        !if idx_mt_iso is fissionable
        !  ogxs(:,3) = One-group Volume-integrated (n,3n) // mt 17 (ENDF)
        !  ogxs(:,4) = One-group Volume-integrated (n,f)  // mt 18 (ENDF)
        !end if
        !if idx_mt_iso is not fissionable
        !  ogxs(:,3) = One-group Volume-integrated (n,alpha) // mt 107 (ENDF)
        !  ogxs(:,4) = One-group Volume-integrated (n,p)     // mt 103 (ENDF)
        !end if
        
        ! Isotopes for S(a,b)?
        logical :: sab = .false.

        ! Doppler broadening
        logical :: db = .false.
        
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
