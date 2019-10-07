module bank_header 

    use constants
    use randoms 
    use geometry_header, only: cells, universes, surfaces
    
    
    implicit none
    
    type :: Bank
        real(8) :: wgt           ! weight of bank site
        real(8) :: xyz(3)        ! location of bank particle
        real(8) :: uvw(3)        ! diretional cosines
        real(8) :: E             ! energy for CE
        integer :: G             ! energy group if in MG mode.
        
      contains 
        !procedure :: initialize => bank_initialize
    end type Bank
    
    ! Source and fission bank
    type(Bank), allocatable, target :: source_bank(:)
    type(Bank), allocatable, target :: fission_bank(:)
    type(Bank), allocatable         :: temp_bank(:)
    type(Bank)                      :: thread_bank(1000)
    !$OMP THREADPRIVATE(thread_bank)
    integer                         :: bank_idx
    !$OMP THREADPRIVATE(bank_idx)
    
    contains 
    
end module
