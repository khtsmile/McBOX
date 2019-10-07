module physics
    use omp_lib
    use variables,            only: keff, k_col, k_tl
    use constants
    use particle_header 
    use XS_header 
    use bank_header,         only: fission_bank, thread_bank, bank_idx
    use randoms
    
    implicit none 
    
    contains
    
    subroutine collision_MG(p)
        type(Particle), intent(inout) :: p
        real(8) :: sig_tot, temp, rnum, wgt_s, uvw_temp(3)
        integer :: i, i_group, idx_group, n_group, n, bsize
        p % n_collision = p % n_collision + 1
        p % n_coord = 1
        
        
        sig_tot = sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g)
        
        !> Collision estimator 
        !$omp atomic
        k_col = k_col + p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)

        !> Fission bank add
        n = int(p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/keff) + rang())
        
        if (n > 0) then
            bank_idx = bank_idx + 1
            thread_bank(bank_idx)%xyz = p%coord(1)%xyz
            thread_bank(bank_idx)%uvw = rand_vec()
            rnum = rang()
            do i_group = 1, size(XS_MG(p%material)%chi(:))
                if (rnum < sum(XS_MG(p%material)%chi(1:i_group))/sum(XS_MG(p%material)%chi(:))) then 
                    thread_bank(bank_idx)%G = i_group
                    exit
                endif 
            enddo
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
        p % last_uvw(:) = p % coord(1)% uvw(:)
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
    
    subroutine collision_MG_DT(p, macro_major)
        type(Particle), intent(inout) :: p
        real(8), intent(in) :: macro_major
        real(8) :: sig_tot, temp, rnum, wgt_s, uvw_temp(3)
        integer :: i, i_group, idx_group, n_group, n, bsize, i_source
        p % n_collision = p % n_collision + 1
        p % n_coord = 1
        
        sig_tot = sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g)
        
        !> Collision estimator 
        !$omp atomic
        k_col = k_col + p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)

        
        rnum = rang()
        if (rnum < (XS_MG(p%material)%sig_abs(p%g) - XS_MG(p%material)%sig_fis(p%g)) / sig_tot) then 
            p%wgt   = 0
            p%alive = .false.
            return
            
        elseif (rnum < (XS_MG(p%material)%sig_abs(p%g)) / sig_tot) then 
            !> Fission bank add
            !n = int(p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/keff) + rang())
            n = int(p%wgt*XS_MG(p%material)%nu(p%g)*(1./keff) + rang())
            do i_source = 1, n
                bank_idx = bank_idx + 1
                thread_bank(bank_idx)%xyz = p%coord(1)%xyz
                thread_bank(bank_idx)%uvw = rand_vec()
                rnum = rang()
                do i_group = 1, size(XS_MG(p%material)%chi(:))
                    if (rnum < sum(XS_MG(p%material)%chi(1:i_group))/sum(XS_MG(p%material)%chi(:))) then 
                        thread_bank(bank_idx)%G = i_group
                        exit
                    endif 
                enddo
            enddo
            p%wgt   = 0
            p%alive = .false.
            return
             
        
        else
            rnum = rang()
            do i_group = 1, size(XS_MG(p%material)%sig_scat(p%g,:))
                if (rnum < sum(XS_MG(p%material)%sig_scat(p%g,1:i_group))/sum(XS_MG(p%material)%sig_scat(p%g,:))) then 
                    idx_group = i_group
                    exit
                endif 
            enddo 
            !p % wgt  = p % wgt * sum(XS_MG(p%material)%sig_scat(p%g,:))/sig_tot
            p % g    = idx_group
            p % coord(1) % uvw(:) = rand_vec()
        endif

    end subroutine

end module 
