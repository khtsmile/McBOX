module physics
	use omp_lib
	use variables,			only: keff, k_col, k_tl
	use constants
	use particle_header 
	use XS_header 
	use bank_header, 		only: fission_bank, thread_bank, bank_idx
	use randoms
	
	implicit none 
	
	contains
	
	subroutine collision_MG(p)
		type(Particle), intent(inout) :: p
		!type(MG_XS),		  pointer :: XSptr
		real(8) :: sig_tot, temp, rnum, wgt_s
		integer :: i_group, idx_group, n_group, n, bsize
		p % n_collision = p % n_collision + 1
		p % n_coord = 1
		
		
		
		!XSptr => XS_MG(p%material)
		!print *, p%material
		
		sig_tot = sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g)
		
		!> Collision estimator 
		!$omp atomic
		k_col = k_col + p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)

		
		!print *, k_col, p%wgt, (XSptr%nu(p%g)*XSptr%sig_fis(p%g)/sig_tot)
		
		!> Fission bank add
		n = int(p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/keff) + rang())
		!print *, p%wgt, p%wgt*(XSptr%nu(p%g)*XSptr%sig_fis(p%g)/sig_tot)*(1/keff) + rnum,rnum
        
		if (n > 0) then
			bank_idx = bank_idx + 1
			!print *, omp_get_thread_num(), bank_idx
			thread_bank(bank_idx)%xyz = p%coord(1)%xyz
			thread_bank(bank_idx)%uvw = rand_vec()
			rnum = rang()
			do i_group = 1, size(XS_MG(p%material)%chi(:))
				if (rnum < sum(XS_MG(p%material)%chi(1:i_group))/sum(XS_MG(p%material)%chi(:))) then 
					thread_bank(bank_idx)%G = i_group
					exit 
				endif 
			enddo


			!bsize = size(thread_bank)
			!if(.not. allocated(temp_bank)) allocate(temp_bank(1:bsize+1)) 
			!if (bsize>0) temp_bank(1:bsize) = thread_bank(:)
            !temp_bank(bsize+1)%xyz = p%coord(1)%xyz
			!temp_bank(bsize+1)%uvw = rand_vec()
			!rnum = rang()
			!do i_group = 1, size(XS_MG(p%material)%chi(:))
			!	if (rnum < sum(XS_MG(p%material)%chi(1:i_group))/sum(XS_MG(p%material)%chi(:))) then 
			!		temp_bank(bsize+1)%G = i_group
			!		exit 
			!	endif
			!enddo
			!deallocate(thread_bank)
			!print *, omp_get_thread_num(), bsize, size(temp_bank), size(thread_bank)
			!print *, thread_bank(:)%xyz(1)
			!call move_alloc(temp_bank, thread_bank)
		endif
		
		rnum = rang()
		do i_group = 1, size(XS_MG(p%material)%sig_scat(p%g,:))
			if (rnum < sum(XS_MG(p%material)%sig_scat(p%g,1:i_group))/sum(XS_MG(p%material)%sig_scat(p%g,:))) then 
				idx_group = i_group
				exit
			endif 
		enddo 
		
		p % wgt = p % wgt * sum(XS_MG(p%material)%sig_scat(p%g,:))/sig_tot
		p % g   = idx_group
		p % coord(1)% uvw(:) = rand_vec()
		
		if (p%wgt < wgt_min) THEN !call Russian_Roulette(p)
			wgt_s = 2*wgt_min
			if ((p%wgt/wgt_s).ge.rang()) then
				p%wgt = wgt_s
			else
				p%wgt = 0 
				p%alive = .false.
			endif 
		endif 
		
		
	end subroutine 
	
	!subroutine Russian_Roulette(p)
	!	type(Particle), intent(inout) :: p
	!	real(8) :: wgt_s
	!	
	!	wgt_s = 2*wgt_min
	!	if ((p%wgt/wgt_s).ge.rang()) then
	!		p%wgt = wgt_s
	!	else
	!		p%wgt = 0 
	!		p%alive = .false.
	!	endif 
	!	
	!end subroutine
	

end module 