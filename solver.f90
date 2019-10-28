module SOLVERS
    use FMFD_HEADER, only: nfm, ncm, fcr, fcz
    implicit none

    contains

! =============================================================================
!         Hepta diagonal matrix solvers
! =============================================================================
function BiCGStab_hepta(M,Q) result(x)
    real(8), intent(in) :: M (:,:,:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: x, r, rs, v, p, s, t
    real(8), parameter :: e = 1D-8
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    real(8) :: summesion, temp
    integer :: it = 0
    integer :: i, j, k

    x     = 0.0
    r     = Q
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sqrt(sum(r*r))
    norm_b = sqrt(sum(Q*Q))
    
    do while ( norm_r .GT. e*norm_b )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)

        v(:,:,:) = 0
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,3)*p(i-1,j,k) ! x0
            if ( i /= nfm(1) ) v(i,j,k) = v(i,j,k) + M(i,j,k,5)*p(i+1,j,k) ! x1
            if ( j /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,2)*p(i,j-1,k) ! y0
            if ( j /= nfm(2) ) v(i,j,k) = v(i,j,k) + M(i,j,k,6)*p(i,j+1,k) ! y1
            if ( k /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,1)*p(i,j,k-1) ! z0
            if ( k /= nfm(3) ) v(i,j,k) = v(i,j,k) + M(i,j,k,7)*p(i,j,k+1) ! z1
                               v(i,j,k) = v(i,j,k) + M(i,j,k,4)*p(i,j,k)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,3)*s(i-1,j,k)
            if ( i /= nfm(1) ) t(i,j,k) = t(i,j,k) + M(i,j,k,5)*s(i+1,j,k)
            if ( j /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,2)*s(i,j-1,k)
            if ( j /= nfm(2) ) t(i,j,k) = t(i,j,k) + M(i,j,k,6)*s(i,j+1,k)
            if ( k /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,1)*s(i,j,k-1)
            if ( k /= nfm(3) ) t(i,j,k) = t(i,j,k) + M(i,j,k,7)*s(i,j,k+1)
                               t(i,j,k) = t(i,j,k) + M(i,j,k,4)*s(i,j,k)
        end do
        end do
        end do
        
        omega  = sum(t*s)/sum(t*t)
        x      = x + alpha*p + omega*s
        r      = s - omega*t
        norm_r = sqrt(sum(r*r))
        norm_b = sqrt(sum(Q*Q))
    
        it = it + 1
    end do   
    
end function BiCGStab_hepta     


! =============================================================================
! BICG_G
! =============================================================================
function BiCG_G(M,Q) result(x)
    real(8), intent(in) :: M (:,:,:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8), dimension(ncm(1),ncm(2),ncm(3)):: x, r, rs, v, p, s, t
    real(8), parameter :: e = 1D-8
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    real(8) :: summesion, temp
    integer :: it = 0
    integer :: i, j, k

    x     = 0.0
    r     = Q
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sqrt(sum(r*r))
    norm_b = sqrt(sum(Q*Q))
    
    do while ( norm_r .GT. e*norm_b )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)

        v(:,:,:) = 0
        do i = 1, ncm(1)
        do j = 1, ncm(2)
        do k = 1, ncm(3)
            if ( i /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,3)*p(i-1,j,k) ! x0
            if ( i /= ncm(1) ) v(i,j,k) = v(i,j,k) + M(i,j,k,5)*p(i+1,j,k) ! x1
            if ( j /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,2)*p(i,j-1,k) ! y0
            if ( j /= ncm(2) ) v(i,j,k) = v(i,j,k) + M(i,j,k,6)*p(i,j+1,k) ! y1
            if ( k /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,1)*p(i,j,k-1) ! z0
            if ( k /= ncm(3) ) v(i,j,k) = v(i,j,k) + M(i,j,k,7)*p(i,j,k+1) ! z1
                               v(i,j,k) = v(i,j,k) + M(i,j,k,4)*p(i,j,k)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do i = 1, ncm(1)
        do j = 1, ncm(2)
        do k = 1, ncm(3)
            if ( i /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,3)*s(i-1,j,k)
            if ( i /= ncm(1) ) t(i,j,k) = t(i,j,k) + M(i,j,k,5)*s(i+1,j,k)
            if ( j /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,2)*s(i,j-1,k)
            if ( j /= ncm(2) ) t(i,j,k) = t(i,j,k) + M(i,j,k,6)*s(i,j+1,k)
            if ( k /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,1)*s(i,j,k-1)
            if ( k /= ncm(3) ) t(i,j,k) = t(i,j,k) + M(i,j,k,7)*s(i,j,k+1)
                               t(i,j,k) = t(i,j,k) + M(i,j,k,4)*s(i,j,k)
        end do
        end do
        end do
        
        omega  = sum(t*s)/sum(t*t)
        x      = x + alpha*p + omega*s
        r      = s - omega*t
        norm_r = sqrt(sum(r*r))
        norm_b = sqrt(sum(Q*Q))
    
        it = it + 1
    end do   
    
end function BiCG_G

! =============================================================================
! BICG_L
! =============================================================================
function BiCG_L(M,Q) result(x)
    real(8), intent(in) :: M (1:nfm(1),1:nfm(2),1:nfm(3),1:7)
    real(8), intent(in) :: Q (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8) :: x (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8), dimension(1:fcr,1:fcr,1:fcz):: xx, r, rs, v, p, s, t
    real(8), parameter :: e = 1D-8
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    integer :: ii, jj, kk, mm, nn, oo
    integer :: id(3), id0(3), id1(3)

    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; id1(1) = id0(1) + fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; id1(2) = id0(2) + fcr
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz; id1(3) = id0(3) + fcz

    xx    = 0.0
    r(1:fcr,1:fcr,1:fcz) = Q(id0(1)+1:id1(1),id0(2)+1:id1(2),id0(3)+1:id1(3))
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sqrt(sum(r*r))
    norm_b = sqrt(sum(r*r))*e

    do while ( norm_r .GT. norm_b )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)

        v(:,:,:) = 0
        do mm = 1, fcr; id(1) = id0(1) + mm
        do nn = 1, fcr; id(2) = id0(2) + nn
        do oo = 1, fcz; id(3) = id0(3) + oo
            if ( mm /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm-1,nn,oo) & ! x0
                                * M(id(1),id(2),id(3),3)
            if ( mm /= fcr ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm+1,nn,oo) & ! x1
                                * M(id(1),id(2),id(3),5)
            if ( nn /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn-1,oo) & ! y0
                                * M(id(1),id(2),id(3),2)
            if ( nn /= fcr ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn+1,oo) & ! y1
                                * M(id(1),id(2),id(3),6)
            if ( oo /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo-1) & ! z0
                                * M(id(1),id(2),id(3),1)
            if ( oo /= fcz ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo+1) & ! z1
                                * M(id(1),id(2),id(3),7)
                             v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo) &
                                * M(id(1),id(2),id(3),4)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do mm = 1, fcr; id(1) = id0(1) + mm
        do nn = 1, fcr; id(2) = id0(2) + nn
        do oo = 1, fcz; id(3) = id0(3) + oo
            if ( mm /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + v(mm-1,nn,oo) & ! x0
                                * M(id(1),id(2),id(3),3)
            if ( mm /= fcr ) t(mm,nn,oo) = t(mm,nn,oo) + v(mm+1,nn,oo) & ! x1
                                * M(id(1),id(2),id(3),5)
            if ( nn /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + v(mm,nn-1,oo) & ! y0
                                * M(id(1),id(2),id(3),2)
            if ( nn /= fcr ) t(mm,nn,oo) = t(mm,nn,oo) + v(mm,nn+1,oo) & ! y1
                                * M(id(1),id(2),id(3),6)
            if ( oo /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + v(mm,nn,oo-1) & ! z0
                                * M(id(1),id(2),id(3),1)
            if ( oo /= fcz ) t(mm,nn,oo) = t(mm,nn,oo) + v(mm,nn,oo+1) & ! z1
                                * M(id(1),id(2),id(3),7)
                             t(mm,nn,oo) = t(mm,nn,oo) + v(mm,nn,oo) &
                                * M(id(1),id(2),id(3),4)
        end do
        end do
        end do
        
        omega  = sum(t*s)/sum(t*t)
        xx     = xx + alpha*p + omega*s
        r      = s - omega*t
        norm_r = sqrt(sum(r*r))

    end do   
    x(id0(1)+1:id1(1),id0(2)+1:id1(2),id0(3)+1:id1(3)) = xx(1:fcr,1:fcr,1:fcz)
    end do
    end do
    end do
    
end function BiCG_L

! =============================================================================
! SOR
! =============================================================================
subroutine SOR(m0,ss,ff)
    implicit none
    real(8), intent(in   ):: m0(:,:,:,:)
    real(8), intent(in   ):: ss(:,:,:)
    real(8), intent(inout):: ff(:,:,:)
    integer:: ii, jj, kk, ee, mm
    real(8):: temp
    integer:: id(3)
    real(8):: relax = 1.4D0
    integer:: n_inner = 5

    do mm=1, n_inner
    do kk=1, nfm(3)
    do jj=1, nfm(2)
    do ii=1, nfm(1)
       temp = ss(ii,jj,kk)
       if ( ii /= 1 )      temp = temp - m0(ii,jj,kk,3)*ff(ii-1,jj,kk)
       if ( ii /= nfm(1) ) temp = temp - m0(ii,jj,kk,5)*ff(ii+1,jj,kk)
       if ( jj /= 1 )      temp = temp - m0(ii,jj,kk,2)*ff(ii,jj-1,kk)
       if ( jj /= nfm(2) ) temp = temp - m0(ii,jj,kk,6)*ff(ii,jj+1,kk)
       if ( kk /= 1 )      temp = temp - m0(ii,jj,kk,1)*ff(ii,jj,kk-1)
       if ( kk /= nfm(3) ) temp = temp - m0(ii,jj,kk,7)*ff(ii,jj,kk+1)
       ff(ii,jj,kk) = (1D0-relax)*ff(ii,jj,kk)+relax*temp/m0(ii,jj,kk,4)
    end do
    end do
    end do
    end do

end subroutine


! =============================================================================
! CG is a matrix solver by the SOR
! =============================================================================
subroutine SORG(m0,ss,ff)
    implicit none
    real(8), intent(in   ):: m0(1:ncm(1),1:ncm(2),1:ncm(3),1:7)
    real(8), intent(in   ):: ss(1:ncm(1),1:ncm(2),1:ncm(3))
    real(8), intent(inout):: ff(1:ncm(1),1:ncm(2),1:ncm(3))
    integer:: ii, jj, kk, ee, mm
    real(8):: temp
    integer:: id(3)
    real(8):: relax = 1.4D0
    integer:: n_inner = 5

    do mm=1, n_inner
    do kk=1, ncm(3)
    do jj=1, ncm(2)
    do ii=1, ncm(1)
       temp = ss(ii,jj,kk)
       if ( ii /= 1 )      temp = temp - m0(ii,jj,kk,3)*ff(ii-1,jj,kk)
       if ( ii /= ncm(1) ) temp = temp - m0(ii,jj,kk,5)*ff(ii+1,jj,kk)
       if ( jj /= 1 )      temp = temp - m0(ii,jj,kk,2)*ff(ii,jj-1,kk)
       if ( jj /= ncm(2) ) temp = temp - m0(ii,jj,kk,6)*ff(ii,jj+1,kk)
       if ( kk /= 1 )      temp = temp - m0(ii,jj,kk,1)*ff(ii,jj,kk-1)
       if ( kk /= ncm(3) ) temp = temp - m0(ii,jj,kk,7)*ff(ii,jj,kk+1)
       ff(ii,jj,kk) = (1D0-relax)*ff(ii,jj,kk)+relax*temp/m0(ii,jj,kk,4)
    end do
    end do
    end do
    end do

end subroutine

subroutine SORL(m0,ss,ff)
    implicit none
    real(8), intent(in   ):: m0(1:nfm(1),1:nfm(2),1:nfm(3),1:7)
    real(8), intent(in   ):: ss(1:nfm(1),1:nfm(2),1:nfm(3))
    real(8), intent(inout):: ff(1:nfm(1),1:nfm(2),1:nfm(3))
    integer:: ii, jj, kk, mm, nn, oo, ll
    real(8):: temp
    integer:: id0(3), id(3)
    real(8):: relax = 1.4D0
    integer:: n_inner = 5

    do ll=1, n_inner
    do kk=1, ncm(3); id0(3) = (kk-1)*fcz
    do jj=1, ncm(2); id0(2) = (jj-1)*fcr
    do ii=1, ncm(1); id0(1) = (ii-1)*fcr
      do mm = 1, fcr; id(1) = id0(1)+mm
      do nn = 1, fcr; id(2) = id0(2)+nn
      do oo = 1, fcz; id(3) = id0(3)+oo
        temp = ss(id(1),id(2),id(3))
        if ( mm /= 1 )   temp = temp - m0(id(1),id(2),id(3),3)*ff(id(1)-1,id(2),id(3))
        if ( mm /= fcr ) temp = temp - m0(id(1),id(2),id(3),5)*ff(id(1)+1,id(2),id(3))
        if ( nn /= 1 )   temp = temp - m0(id(1),id(2),id(3),2)*ff(id(1),id(2)-1,id(3))
        if ( nn /= fcr ) temp = temp - m0(id(1),id(2),id(3),6)*ff(id(1),id(2)+1,id(3))
        if ( oo /= 1 )   temp = temp - m0(id(1),id(2),id(3),1)*ff(id(1),id(2),id(3)-1)
        if ( oo /= fcz ) temp = temp - m0(id(1),id(2),id(3),7)*ff(id(1),id(2),id(3)+1)
        ff(id(1),id(2),id(3)) = &
            (1D0-relax)*ff(id(1),id(2),id(3))+relax*temp/m0(id(1),id(2),id(3),4)
      end do
      end do
      end do
    end do
    end do
    end do
    end do

end subroutine

end module
!    !$omp parallel default(shared) private(ii,jj,kk,ll,id0,id,temp)
!    !$omp do
!    !$omp end do
!    !$omp end parallel
