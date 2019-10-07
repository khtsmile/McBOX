module VRC 
	use constants, 			only: wgt_min 
    use particle_header,    only: particle
	use material_header,	only: materials 
	use XS_header, 			only: XS_MG
	use variables, 			only: E_mode, k_vrc, fiss_vrc, loss_vrc
	use ace_xs, 			only: getMacroXS
	use geometry, 			only: distance_to_boundary, find_cell, cross_surface
	use randoms, 			only: rang
    use omp_lib
	
	implicit none 
	
	type RayType
		real(8) :: wgt
		real(8) :: xyz(3)
		real(8) :: uvw(3) 
		real(8) :: E 
		integer :: G 
		real(8), allocatable :: Sigt(:)
		real(8), allocatable :: nuSigf(:) 
		
		contains
        procedure :: reset => reset_ray
	endtype 
	type(RayType) :: RayBuffer(1000) 
	
	contains 
	
	!==============================================================
	! TRACE_PSUDORAY 
	!==============================================================
	subroutine trace_psudoray(p)
        type(Particle), intent(inout) :: p
		real(8) :: val
        integer :: j                      ! coordinate level
        integer :: surface_crossed        ! surface which particle is on
        real(8) :: d_boundary             ! distance to nearest boundary
        logical :: found_cell             ! found cell which particle is in?
        real(8) :: macro_xs(5)
        integer :: i_cell
		real(8) :: wgt0, wgt_s, wgt
		integer :: iter 
		
		wgt0 = p%wgt
		val = 0 
		!print *, 'psudo-ray start'
		!found_cell = .false.
		call find_cell(p, found_cell, i_cell)
        do while (p%alive == .true.)
			! calculate reduced weight 
			!if (.not. found_cell) call find_cell(p, found_cell, i_cell)
			!print *, p%coord(1)%xyz(1), XS_MG(p%material)%mat_id
			if (E_mode == 0) then 
				macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
				macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
				macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
			elseif (E_mode == 1) then 
				macro_xs = getMacroXS(materials(p%material), p%E)
			endif 
			
			! Sample a distance to boundary
			call distance_to_boundary(p, d_boundary, surface_crossed)
			p%wgt = p%wgt*exp(-val)
			val = d_boundary * macro_xs(1)
			!> Volumetric-ray-casting estimator
			!$omp atomic
			fiss_vrc = fiss_vrc + p%wgt*macro_xs(4)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
			!$omp atomic
			loss_vrc = loss_vrc + p%wgt*macro_xs(2)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
			!> Advance particle
			do j = 1, p % n_coord
				p % coord(j) % xyz = p % coord(j) % xyz + d_boundary * p % coord(j) % uvw
			enddo
			
			call cross_surface(p, surface_crossed)

			!if (p%wgt < 0.000000001) then 
			!	p%alive = .false.
			!endif 
			! Russian Roulette for psudo-ray
			if (p%wgt < 0.0001) THEN !call Russian_Roulette(p)
				wgt_s = 2*0.0001
				if ((p%wgt/wgt_s).ge.rang()) then
					p%wgt = wgt_s
				else
					p%alive = .false.
				endif 
			endif 
			
		enddo 
		
	end subroutine	
	
	
	
	!subroutine trace_psudoray(p, val)
    !    type(Particle), intent(inout) :: p
	!	real(8), intent(inout) :: val
	!	
    !    integer :: j                      ! coordinate level
    !    integer :: surface_crossed        ! surface which particle is on
    !    real(8) :: d_boundary             ! distance to nearest boundary
    !    logical :: found_cell             ! found cell which particle is in?
    !    real(8) :: macro_xs(5)
    !    integer :: i_cell
	!	real(8) :: wgt0, wgt_s
	!	integer :: iter 
	!	
	!	wgt0 = p%wgt
    !    do while (p%alive == .true.)
	!		! calculate reduced weight 
	!		p%wgt = wgt0 * exp(-val)
	!		found_cell = .false.
	!		call find_cell(p, found_cell, i_cell)
	!		!print *, XS_MG(p%material)%mat_id, p%wgt, exp(-val)
	!		
	!		if (E_mode == 0) then 
	!			macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
	!			macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
	!			macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
	!		elseif (E_mode == 1) then 
	!			macro_xs = getMacroXS(materials(p%material), p%E)
	!		endif 
    !
	!		! Sample a distance to boundary
	!		call distance_to_boundary(p, d_boundary, surface_crossed)
	!		val = val + d_boundary * macro_xs(1)
	!		
	!		!> Volumetric-ray-casting estimator
	!		!$omp atomic
	!		fiss_vrc = fiss_vrc + p%wgt*macro_xs(4)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
	!		!$omp atomic
	!		loss_vrc = loss_vrc + p%wgt*macro_xs(2)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
	!		
	!		!> Advance particle
	!		do j = 1, p % n_coord
	!			p % coord(j) % xyz = p % coord(j) % xyz + d_boundary * p % coord(j) % uvw
	!		enddo
	!		!print *, p%coord(1)%xyz(1:2) 
	!		call cross_surface(p, surface_crossed)
	!		!if (p%alive == .false.) then ! particle leak
	!		!	!$omp atomic
	!		!	loss_vrc = loss_vrc + wgt0 * exp(-val)!p%wgt!
	!		!endif
	!		if (p%wgt < 0.0000001) then 
	!			p%alive = .false.
	!			!print *, 'killed ************'
	!		endif 
	!		! Russian Roulette for psudo-ray
	!		!if (p%wgt < wgt_min) THEN !call Russian_Roulette(p)
	!		!	wgt_s = 2*wgt_min
	!		!	if ((p%wgt/wgt_s).ge.rang()) then
	!		!		p%wgt = wgt_s
	!		!	else
	!		!		p%alive = .false.
	!		!	endif 
	!		!endif 
	!		
	!	enddo 
	!	
	!end subroutine
	
	
	!==============================================================
	! CREATE_RAY 
	!==============================================================
	subroutine create_Ray (p, ray) 
        type(Particle), intent(in) :: p
		type(RayType), intent(inout) :: ray
		integer :: i, n_mat
		real(8) :: macro_xs(5)
		
		ray%wgt    = p%wgt
		ray%xyz(:) = p%coord(1)%xyz(:) 
		ray%uvw(:) = p%coord(1)%uvw(:) 
		
		if (E_mode == 0) then 
			ray%G = p%G
			n_mat = size(XS_MG)
			allocate(ray%sigt(n_mat)) 
			allocate(ray%nuSigf(n_mat)) 
			
			do i = 1, n_mat
				ray%sigt(i)  = (sum(XS_MG(i)%sig_scat(p%g,:)) + XS_MG(i)%sig_abs(p%g))
				ray%nuSigf(i) = XS_MG(i)%sig_fis(p%g)*XS_MG(i)%nu(p%g)
			enddo 
			
		else 
			ray%E = p%E
			n_mat = size(materials)
			allocate(ray%sigt(n_mat)) 
			allocate(ray%nuSigf(n_mat))
			
			do i = 1, n_mat
				macro_xs = getMacroXS (materials(i), ray%E)
				ray%sigt(i)  = macro_xs(1)
				ray%nuSigf(i) = macro_xs(4)
			enddo 
			
		endif
		
	end subroutine
	
	
	!==============================================================
	! RESET_RAY clears data of a ray
	!==============================================================
    elemental subroutine reset_ray(this)
        class(RayType), intent(inout) :: this
		
        this % xyz = 0 
        this % uvw = 0 
        this % E   = 0 
        this % G   = 0 
		if (allocated(this%sigt)) deallocate(this%sigt)
		if (allocated(this%nuSigf)) deallocate(this%nuSigf)
		
    end subroutine reset_ray
	
end module 