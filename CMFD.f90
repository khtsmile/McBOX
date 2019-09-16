module CMFD
    use FMFD_HEADER
    implicit none


    contains


! =============================================================================
! ONE_NODE_CMFD
! =============================================================================
subroutine ONE_NODE_CMFD(keff,fm_t,fm_a,fm_nf,fmD,fm_phi1,fmJ0,fmJ1,fmJn,fmF)
    use SOLVERS, only: BICG_G, SORL, BiCG_L, SORG
    implicit none
    real(8), intent(inout):: keff
    real(8), intent(in), dimension(:,:,:):: fm_t, fm_a, fm_nf, fmD
    real(8), intent(inout):: fm_phi1(:,:,:)
    real(8), intent(inout), dimension(:,:,:,:):: fmJ0, fmJ1, fmJn, fmF
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: &
        fm_phi0, &  ! neutron flux
        fm_s        ! neutron source
    real(8), dimension(nfm(1),nfm(2),nfm(3),6):: &
        fmDt, &     ! D tilda
        fmDh        ! D hat
    real(8) :: error, k_pre
    integer :: global, local

    call L2G(fm_phi1,fm_t,fm_a,fm_nf,fmJn,fmF)
    call L_DTILDA(fmD,fmDt)
    call L_DHAT(fmDt,fm_phi1,fmJn,fmDh)
    call D_BC(fmD,fmDt)
    call L_BC(fmDt,fmDh)
    call L_MATRIX(fmDt,fmDh,fm_a)

    do
    ! ------------------------------- GLOBAL
    error = 1D0
    call G_DHAT(cmJn,cmDt,cm_phi1,cmF,cmDh)
    call G_MATRIX(cmDt,cmDh)
    k_pre = keff

    do global = 1, 5
    cm_phi0 = cm_phi1
    cm_s = cm_nf*cm_phi0/keff
    cm_phi1 = BiCG_G(Mcm,cm_s)
    !call SORG(Mcm,cm_s,cm_phi1)
    keff = keff*sum(cm_nf*cm_phi1*cm_nf*cm_phi1) &
           / sum(cm_nf*cm_phi0*cm_nf*cm_phi1)
    end do
    error = abs(keff-k_pre)/keff
    !print*, keff, error
    if ( error < 1D-8 .or. isnan(keff) ) exit
    ! ------------------------------- LOCAL
    call G_INJ(cmDt,cmDh,cm_phi1)

    do local = 1, 2
    call G2L(fm_phi0,fm_phi1,fmJ0,fmJ1)
    call L_SOURCE(fm_phi0,fm_phi1,keff,fm_nf,fm_s,fmJ0,fmJ1)
!    fm_phi1(:,:,:) = BICG_L(Mfm(:,:,:,:),fm_s(:,:,:))
    call SORL(Mfm,fm_s,fm_phi1)
    call L_OUTJ(fm_phi0,fm_phi1,fmF,fmJ0,fmJ1,fmJn)
    call L_REFJ(fmF,fmJ0,fmJ1,fmJn)
    call G_XS(fm_t,fm_a,fm_nf,fm_phi1)
    end do
    end do



end subroutine

! =============================================================================
! L2G homogenizes the reactor parameters from local to global
! =============================================================================
subroutine L2G(phi,fm_tot,fm_abso,fm_nufiss,fmJn,fmF)
    implicit none
    real(8), intent(in), dimension(:,:,:):: phi, fm_tot, fm_abso, fm_nufiss
    real(8), intent(in), dimension(:,:,:,:):: fmJn, fmF
    real(8):: ssum(2)

    ! -------------------------------------------------------------------------
    ! homogenization
    do ii = 1, ncm(1); id(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id(2) = (jj-1)*fcr
    do kk = 1, ncm(3); id(3) = (kk-1)*fcz
        cm_phi1(ii,jj,kk) = sum(phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))
        cm_t(ii,jj,kk) = sum(fm_tot(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz)*phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))/cm_phi1(ii,jj,kk)
        cm_a(ii,jj,kk) = sum(fm_abso(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz)*phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))/cm_phi1(ii,jj,kk)
        cm_nf(ii,jj,kk) = sum(fm_nufiss(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz)*phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))/cm_phi1(ii,jj,kk)
    end do
    end do
    end do
    cmD = 1D0 / (3D0 * cm_t)
    where ( cm_t == 0 ) cmD = 0
    cm_phi1 = cm_phi1 / (fcr*fcr*fcz)

    ! interface diffusion coefficient
    do ii = 1, ncm(1)
    do jj = 1, ncm(2)
    do kk = 1, ncm(3)
        cmDt(ii,jj,kk,1) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(1))
        cmDt(ii,jj,kk,2) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(1))
        cmDt(ii,jj,kk,3) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(2))
        cmDt(ii,jj,kk,4) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(2))
        cmDt(ii,jj,kk,5) = 2D0*cmD(ii,jj,kk)/(fcz*dfm(3))
        cmDt(ii,jj,kk,6) = 2D0*cmD(ii,jj,kk)/(fcz*dfm(3))
    end do
    end do
    end do

    ! -------------------------------------------------------------------------
    ! surface average
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        ! x-direction
        ssum = 0;       id(1) = id0(1)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),1)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),1)
        end do
        end do
        cmJn(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,1)  = ssum(2) / (fcr*fcz)
        if ( ii /= 1 ) then
        cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)
        cmF(ii-1,jj,kk,2)  = cmF(ii,jj,kk,1)
        end if
        ! y-direction
        ssum = 0;       id(2) = id0(2)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),3)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),3)
        end do
        end do
        cmJn(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,3)  = ssum(2) / (fcr*fcz)
        if ( jj /= 1 ) then
        cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)
        cmF(ii,jj-1,kk,4)  = cmF(ii,jj,kk,3)
        end if
        ! z-direction
        ssum = 0;       id(3) = id0(3)+1
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),5)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),5)
        end do
        end do
        cmJn(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,5)  = ssum(2) / (fcr*fcr)
        if ( kk /= 1 ) then
        cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)
        cmF(ii,jj,kk-1,6)  = cmF(ii,jj,kk,5)
        end if
    end do
    end do
    end do
    ! Closure
    !   x-direction
    ii = ncm(1);       id(1) = nfm(1)
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),2)
        ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),2)
    end do
    end do
    cmJn(ii,jj,kk,2) = ssum(1) / (fcr*fcz)
    cmF(ii,jj,kk,2)  = ssum(2) / (fcr*fcz)
    end do
    end do
    !   y-direction
    jj = ncm(2);       id(2) = nfm(2)
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),4)
        ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),4)
    end do
    end do
    cmJn(ii,jj,kk,4) = ssum(1) / (fcr*fcz)
    cmF(ii,jj,kk,4)  = ssum(2) / (fcr*fcz)
    end do
    end do
    !   z-direction
    kk = ncm(3);       id(3) = nfm(3)
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),6)
        ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),6)
    end do
    end do
    cmJn(ii,jj,kk,6) = ssum(1) / (fcr*fcr)
    cmF(ii,jj,kk,6)  = ssum(2) / (fcr*fcr)
    end do
    end do

end subroutine


! =============================================================================
! D_TILDA
! =============================================================================
subroutine L_DTILDA(D,Dt)
    implicit none
    real(8), intent(in) :: D(:,:,:)
    real(8), intent(out):: Dt(:,:,:,:)
    
    ! inner region
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        if ( ii /= 1 ) then         ! x0
            Dt(ii,jj,kk,1) = 2D0*D(ii,jj,kk)*D(ii-1,jj,kk) &
                /(D(ii,jj,kk)+D(ii-1,jj,kk))/dfm(1)
            deltf0(ii,jj,kk,1) = Dt(ii,jj,kk,1)
        end if
        if ( ii /= nfm(1) ) then    ! x1
            Dt(ii,jj,kk,2) = 2D0*D(ii+1,jj,kk)*D(ii,jj,kk) &
                /(D(ii+1,jj,kk)+D(ii,jj,kk))/dfm(1)
            deltf0(ii,jj,kk,2) = Dt(ii,jj,kk,2)
        end if
        if ( jj /= 1 ) then         ! y0
            Dt(ii,jj,kk,3) = 2D0*D(ii,jj,kk)*D(ii,jj-1,kk) &
                /(D(ii,jj,kk)+D(ii,jj-1,kk))/dfm(2)
            deltf0(ii,jj,kk,3) = Dt(ii,jj,kk,3)
        end if
        if ( jj /= nfm(2) ) then    ! y1
            Dt(ii,jj,kk,4) = 2D0*D(ii,jj+1,kk)*D(ii,jj,kk) &
                /(D(ii,jj+1,kk)+D(ii,jj,kk))/dfm(2)
            deltf0(ii,jj,kk,4) = Dt(ii,jj,kk,4)
        end if
        if ( kk /= 1 ) then         ! z0
            Dt(ii,jj,kk,5) = 2D0*D(ii,jj,kk)*D(ii,jj,kk-1) &
                /(D(ii,jj,kk)+D(ii,jj,kk-1))/dfm(3)
            deltf0(ii,jj,kk,5) = Dt(ii,jj,kk,5)
        end if
        if ( kk /= nfm(3) ) then    ! z1
            Dt(ii,jj,kk,6) = 2D0*D(ii,jj,kk+1)*D(ii,jj,kk) &
                /(D(ii,jj,kk+1)+D(ii,jj,kk))/dfm(3)
            deltf0(ii,jj,kk,6) = Dt(ii,jj,kk,6)
        end if
    end do
    end do
    end do

end subroutine

! =============================================================================
! L_DHAT
! =============================================================================
subroutine L_DHAT(Dt,phi,Jn,Dh)
    implicit none
    real(8):: Dt(:,:,:,:), phi(:,:,:), Jn(:,:,:,:), Dh(:,:,:,:)

    do kk = 1, nfm(3)
    do jj = 1, nfm(2)
    do ii = 1, nfm(1)
        if ( ii /= 1 )      Dh(ii,jj,kk,1) = (Jn(ii,jj,kk,1)+Dt(ii,jj,kk,1) &
            *(phi(ii,jj,kk)-phi(ii-1,jj,kk)))/(phi(ii,jj,kk)+phi(ii-1,jj,kk))
        if ( ii /= nfm(1) ) Dh(ii,jj,kk,2) = (Jn(ii,jj,kk,2)+Dt(ii,jj,kk,2) &
            *(phi(ii+1,jj,kk)-phi(ii,jj,kk)))/(phi(ii+1,jj,kk)+phi(ii,jj,kk))
        if ( jj /= 1 )      Dh(ii,jj,kk,3) = (Jn(ii,jj,kk,3)+Dt(ii,jj,kk,3) &
            *(phi(ii,jj,kk)-phi(ii,jj-1,kk)))/(phi(ii,jj,kk)+phi(ii,jj-1,kk))
        if ( jj /= nfm(2) ) Dh(ii,jj,kk,4) = (Jn(ii,jj,kk,4)+Dt(ii,jj,kk,4) &
            *(phi(ii,jj+1,kk)-phi(ii,jj,kk)))/(phi(ii,jj+1,kk)+phi(ii,jj,kk))
        if ( kk /= 1 )      Dh(ii,jj,kk,5) = (Jn(ii,jj,kk,5)+Dt(ii,jj,kk,5) &
            *(phi(ii,jj,kk)-phi(ii,jj,kk-1)))/(phi(ii,jj,kk)+phi(ii,jj,kk-1))
        if ( kk /= nfm(3) ) Dh(ii,jj,kk,6) = (Jn(ii,jj,kk,6)+Dt(ii,jj,kk,6) &
            *(phi(ii,jj,kk+1)-phi(ii,jj,kk)))/(phi(ii,jj,kk+1)+phi(ii,jj,kk))
    end do
    end do
    end do

    ! Boundary condition
    ii = 1;      Dh(ii,:,:,1) = Jn(ii,:,:,1)/phi(ii,:,:)
    ii = nfm(1); Dh(ii,:,:,2) = Jn(ii,:,:,2)/phi(ii,:,:)
    jj = 1;      Dh(:,jj,:,3) = Jn(:,jj,:,3)/phi(:,jj,:)
    jj = nfm(2); Dh(:,jj,:,4) = Jn(:,jj,:,4)/phi(:,jj,:)
    kk = 1;      Dh(:,:,kk,5) = Jn(:,:,kk,5)/phi(:,:,kk)
    kk = nfm(3); Dh(:,:,kk,6) = Jn(:,:,kk,6)/phi(:,:,kk)

end subroutine

! =============================================================================
! D_BC
! =============================================================================
subroutine D_BC(D,Dt)
    implicit none
    real(8):: D(:,:,:), Dt(:,:,:,:)

    ! diffusion coefficient at boundary
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        id(1) = id0(1)+1    ! x0
        Dt(id(1),:,:,1) = 2D0*D(id(1),:,:)/dfm(1); deltf0(id(1),:,:,1) = 0D0
        id(1) = id0(1)+fcr  ! x1
        Dt(id(1),:,:,2) = 2D0*D(id(1),:,:)/dfm(1); deltf0(id(1),:,:,2) = 0D0
        id(2) = id0(2)+1    ! y0
        Dt(:,id(2),:,3) = 2D0*D(:,id(2),:)/dfm(2); deltf0(:,id(2),:,3) = 0D0
        id(2) = id0(2)+fcr  ! y1
        Dt(:,id(2),:,4) = 2D0*D(:,id(2),:)/dfm(2); deltf0(:,id(2),:,4) = 0D0
        id(3) = id0(3)+1    ! z0
        Dt(:,:,id(3),5) = 2D0*D(:,:,id(3))/dfm(3); deltf0(:,:,id(3),5) = 0D0
        id(3) = id0(3)+fcz  ! z1
        Dt(:,:,id(3),6) = 2D0*D(:,:,id(3))/dfm(3); deltf0(:,:,id(3),6) = 0D0
    end do
    end do
    end do

end subroutine


! =============================================================================
! L_BC
! =============================================================================
subroutine L_BC(Dt,Dh)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:)

    deltf1 = Dh

    ! interface boundary
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( ii /= 1 ) then;        id(1) = id0(1)+1    ! x0
        deltf1(id(1),:,:,1) = -(Dt(id(1),:,:,1)-Dh(id(1),:,:,1)) &
                            / (1D0+2D0*Dt(id(1),:,:,1))
        end if
        if ( ii /= ncm(1) ) then;   id(1) = id0(1)+fcr  ! x1
        deltf1(id(1),:,:,2) = (Dt(id(1),:,:,2)+Dh(id(1),:,:,2)) &
                            / (1D0+2D0*Dt(id(1),:,:,2))
        end if
        if ( jj /= 1 ) then;        id(2) = id0(2)+1    ! y0
        deltf1(:,id(2),:,3) = -(Dt(:,id(2),:,3)-Dh(:,id(2),:,3)) &
                            / (1D0+2D0*Dt(:,id(2),:,3))
        end if
        if ( jj /= ncm(2) ) then;   id(2) = id0(2)+fcr  ! y1
        deltf1(:,id(2),:,4) = (Dt(:,id(2),:,4)+Dh(:,id(2),:,4)) &
                            / (1D0+2D0*Dt(:,id(2),:,4))
        end if
        if ( kk /= 1 ) then;        id(3) = id0(3)+1    ! z0
        deltf1(:,:,id(3),5) = -(Dt(:,:,id(3),5)-Dh(:,:,id(3),5)) &
                            / (1D0+2D0*Dt(:,:,id(3),5))
        end if
        if ( kk /= ncm(3) ) then;   id(3) = id0(3)+fcz  ! z1
        deltf1(:,:,id(3),6) = (Dt(:,:,id(3),6)+Dh(:,:,id(3),6)) &
                            / (1D0+2D0*Dt(:,:,id(3),6))
        end if
    end do
    end do
    end do

end subroutine

! =============================================================================
! L_MATRIX
! =============================================================================
subroutine L_MATRIX(Dt,Dh,abso)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:), abso(:,:,:)
    real(8):: deno(nfm(1),nfm(2),nfm(3))  ! denominator
    real(8):: deno1

    ! Matrix formulation

    ! -------------------------------------------------------------------------
    !   migration term
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)

        if ( kk /= 1   ) Mfm(ii,jj,kk,1) = &
                -(deltf0(ii,jj,kk,5)+deltf1(ii,jj,kk,5))/dfm(3)
        if ( jj /= 1   ) Mfm(ii,jj,kk,2) = &
                -(deltf0(ii,jj,kk,3)+deltf1(ii,jj,kk,3))/dfm(2)
        if ( ii /= 1   ) Mfm(ii,jj,kk,3) = &
                -(deltf0(ii,jj,kk,1)+deltf1(ii,jj,kk,1))/dfm(1)
        if ( ii /= fcr ) Mfm(ii,jj,kk,5) = &
                -(deltf0(ii,jj,kk,2)-deltf1(ii,jj,kk,2))/dfm(1)
        if ( jj /= fcr ) Mfm(ii,jj,kk,6) = &
                -(deltf0(ii,jj,kk,4)-deltf1(ii,jj,kk,4))/dfm(2)
        if ( kk /= fcz ) Mfm(ii,jj,kk,7) = &
                -(deltf0(ii,jj,kk,6)-deltf1(ii,jj,kk,6))/dfm(3)
        
        Mfm(ii,jj,kk,4) = &
            +(deltf0(ii,jj,kk,1)-deltf1(ii,jj,kk,1))/dfm(1) &
            +(deltf0(ii,jj,kk,2)+deltf1(ii,jj,kk,2))/dfm(1) &
            +(deltf0(ii,jj,kk,3)-deltf1(ii,jj,kk,3))/dfm(2) &
            +(deltf0(ii,jj,kk,4)+deltf1(ii,jj,kk,4))/dfm(2) &
            +(deltf0(ii,jj,kk,5)-deltf1(ii,jj,kk,5))/dfm(3) &
            +(deltf0(ii,jj,kk,6)+deltf1(ii,jj,kk,6))/dfm(3) &
            +abso(ii,jj,kk)

    end do
    end do
    end do

    ! -------------------------------------------------------------------------
    !   source term
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        if ( ii /= 1 ) then;        id(1) = id0(1)+1    ! x0
            deno(id(1),:,:) = 1D0+2D0*Dt(id(1),:,:,1)
            jsrc(id(1),:,:,1) = 4D0*Dt(id(1),:,:,1)/deno(id(1),:,:)
            fsrc(id(1),:,:,1) = Dh(id(1),:,:,1)/deno(id(1),:,:)
        end if
        if ( ii /= ncm(1) ) then;   id(1) = id0(1)+fcr  ! x1
            deno(id(1),:,:) = 1D0+2D0*Dt(id(1),:,:,2)
            jsrc(id(1),:,:,2) = 4D0*Dt(id(1),:,:,2)/deno(id(1),:,:)
            fsrc(id(1),:,:,2) = Dh(id(1),:,:,2)/deno(id(1),:,:)
        end if
        if ( jj /= 1 ) then;        id(2) = id0(2)+1    ! y0
            deno(:,id(2),:) = 1D0+2D0*Dt(:,id(2),:,3)
            jsrc(:,id(2),:,3) = 4D0*Dt(:,id(2),:,3)/deno(:,id(2),:)
            fsrc(:,id(2),:,3) = Dh(:,id(2),:,3)/deno(:,id(2),:)
        end if
        if ( jj /= ncm(2) ) then;   id(2) = id0(2)+fcr  ! y1
            deno(:,id(2),:) = 1D0+2D0*Dt(:,id(2),:,4)
            jsrc(:,id(2),:,4) = 4D0*Dt(:,id(2),:,4)/deno(:,id(2),:)
            fsrc(:,id(2),:,4) = Dh(:,id(2),:,4)/deno(:,id(2),:)
        end if
        if ( kk /= 1 ) then;        id(3) = id0(3)+1
            deno(:,:,id(3)) = 1D0+2D0*Dt(:,:,id(3),5)
            jsrc(:,:,id(3),5) = 4D0*Dt(:,:,id(3),5)/deno(:,:,id(3))
            fsrc(:,:,id(3),5) = Dh(:,:,id(3),5)/deno(:,:,id(3))
        end if
        if ( kk /= ncm(3) ) then;   id(3) = id0(3)+fcz
            deno(:,:,id(3)) = 1D0+2D0*Dt(:,:,id(3),6)
            jsrc(:,:,id(3),6) = 4D0*Dt(:,:,id(3),6)/deno(:,:,id(3))
            fsrc(:,:,id(3),6) = Dh(:,:,id(3),6)/deno(:,:,id(3))
        end if
    end do
    end do
    end do

end subroutine

! =============================================================================
! G_DHAT
! =============================================================================
subroutine G_DHAT(Jn,Dt,vphi,sphi,Dh)
    implicit none
    real(8), intent(in) :: Jn(:,:,:,:), Dt(:,:,:,:), vphi(:,:,:), sphi(:,:,:,:)
    real(8), intent(out):: Dh(:,:,:,:)

    ! x0 +
    Dh(:,:,:,1) = (Jn(:,:,:,1)+Dt(:,:,:,1) &
        *(vphi(:,:,:)-sphi(:,:,:,1)))/(vphi(:,:,:)+sphi(:,:,:,1))
    ! x1 -
    Dh(:,:,:,2) = (Jn(:,:,:,2)+Dt(:,:,:,2) &
        *(sphi(:,:,:,2)-vphi(:,:,:)))/(sphi(:,:,:,2)+vphi(:,:,:))
    ! y0 +
    Dh(:,:,:,3) = (Jn(:,:,:,3)+Dt(:,:,:,3) &
        *(vphi(:,:,:)-sphi(:,:,:,3)))/(vphi(:,:,:)+sphi(:,:,:,3))
    ! y1 -
    Dh(:,:,:,4) = (Jn(:,:,:,4)+Dt(:,:,:,4) &
        *(sphi(:,:,:,4)-vphi(:,:,:)))/(sphi(:,:,:,4)+vphi(:,:,:))
    ! z0 +
    Dh(:,:,:,5) = (Jn(:,:,:,5)+Dt(:,:,:,5) &
        *(vphi(:,:,:)-sphi(:,:,:,5)))/(vphi(:,:,:)+sphi(:,:,:,5))
    ! z1 -
    Dh(:,:,:,6) = (Jn(:,:,:,6)+Dt(:,:,:,6) &
        *(sphi(:,:,:,6)-vphi(:,:,:)))/(sphi(:,:,:,6)+vphi(:,:,:))

    ! boundary condition
    ii = 1;      Dh(ii,:,:,1) = Jn(ii,:,:,1) / vphi(ii,:,:)
    ii = ncm(1); Dh(ii,:,:,2) = Jn(ii,:,:,2) / vphi(ii,:,:)
    jj = 1;      Dh(:,jj,:,3) = Jn(:,jj,:,3) / vphi(:,jj,:)
    jj = ncm(2); Dh(:,jj,:,4) = Jn(:,jj,:,4) / vphi(:,jj,:)
    kk = 1;      Dh(:,:,kk,5) = Jn(:,:,kk,5) / vphi(:,:,kk)
    kk = ncm(3); Dh(:,:,kk,6) = Jn(:,:,kk,6) / vphi(:,:,kk)

end subroutine


! =============================================================================
! G_MATRIX
! =============================================================================
subroutine G_MATRIX(Dt,Dh)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:)
    real(8):: deno    ! denominator of the parameter

    ! diffusion coefficient
    do ii = 1, ncm(1)
    do jj = 1, ncm(2)
    do kk = 1, ncm(3)
        if ( ii /= 1 ) then         ! x0
        deno = Dt(ii-1,jj,kk,2)+Dt(ii,jj,kk,1)+Dh(ii,jj,kk,1)-Dh(ii-1,jj,kk,2)
        deltc0(ii,jj,kk,1) = (Dt(ii-1,jj,kk,2)*Dt(ii,jj,kk,1) &
            +Dh(ii,jj,kk,1)*Dh(ii-1,jj,kk,2))/deno
        deltc1(ii,jj,kk,1) = (Dt(ii-1,jj,kk,2)*Dh(ii,jj,kk,1) &
            +Dt(ii,jj,kk,1)*Dh(ii-1,jj,kk,2))/deno
        end if
        if ( ii /= ncm(1) ) then    ! x1
        deno = Dt(ii,jj,kk,2)+Dt(ii+1,jj,kk,1)+Dh(ii+1,jj,kk,1)-Dh(ii,jj,kk,2)
        deltc0(ii,jj,kk,2) = (Dt(ii,jj,kk,2)*Dt(ii+1,jj,kk,1) &
            +Dh(ii+1,jj,kk,1)*Dh(ii,jj,kk,2))/deno
        deltc1(ii,jj,kk,2) = (Dt(ii,jj,kk,2)*Dh(ii+1,jj,kk,1) &
            +Dt(ii+1,jj,kk,1)*Dh(ii,jj,kk,2))/deno
        end if
        if ( jj /= 1 ) then         ! y0
        deno = Dt(ii,jj-1,kk,4)+Dt(ii,jj,kk,3)+Dh(ii,jj,kk,3)-Dh(ii,jj-1,kk,4)
        deltc0(ii,jj,kk,3) = (Dt(ii,jj-1,kk,4)*Dt(ii,jj,kk,3) &
            +Dh(ii,jj,kk,3)*Dh(ii,jj-1,kk,4))/deno
        deltc1(ii,jj,kk,3) = (Dt(ii,jj-1,kk,4)*Dh(ii,jj,kk,3) &
            +Dt(ii,jj,kk,3)*Dh(ii,jj-1,kk,4))/deno
        end if
        if ( jj /= ncm(2) ) then    ! y1
        deno = Dt(ii,jj,kk,4)+Dt(ii,jj+1,kk,3)+Dh(ii,jj+1,kk,3)-Dh(ii,jj,kk,4)
        deltc0(ii,jj,kk,4) = (Dt(ii,jj,kk,4)*Dt(ii,jj+1,kk,3) &
            +Dh(ii,jj+1,kk,3)*Dh(ii,jj,kk,4))/deno
        deltc1(ii,jj,kk,4) = (Dt(ii,jj,kk,4)*Dh(ii,jj+1,kk,3) &
            +Dt(ii,jj+1,kk,3)*Dh(ii,jj,kk,4))/deno
        end if
        if ( kk /= 1 ) then         ! z0
        deno = Dt(ii,jj,kk-1,6)+Dt(ii,jj,kk,5)+Dh(ii,jj,kk,5)-Dh(ii,jj,kk-1,6)
        deltc0(ii,jj,kk,5) = (Dt(ii,jj,kk-1,6)*Dt(ii,jj,kk,5) &
            +Dh(ii,jj,kk,5)*Dh(ii,jj,kk-1,6))/deno
        deltc1(ii,jj,kk,5) = (Dt(ii,jj,kk-1,6)*Dh(ii,jj,kk,5) &
            +Dt(ii,jj,kk,5)*Dh(ii,jj,kk-1,6))/deno
        end if
        if ( kk /= ncm(3) ) then    ! z1
        deno = Dt(ii,jj,kk,6)+Dt(ii,jj,kk+1,5)+Dh(ii,jj,kk+1,5)-Dh(ii,jj,kk,6)
        deltc0(ii,jj,kk,6) = (Dt(ii,jj,kk,6)*Dt(ii,jj,kk+1,5) &
            +Dh(ii,jj,kk+1,5)*Dh(ii,jj,kk,6))/deno
        deltc1(ii,jj,kk,6) = (Dt(ii,jj,kk,6)*Dh(ii,jj,kk+1,5) &
            +Dt(ii,jj,kk+1,5)*Dh(ii,jj,kk,6))/deno
        end if
    end do
    end do
    end do
    ! boundary condition (J/phi)
    ii = 1;      deltc1(ii,:,:,1) = Dh(ii,:,:,1)
    ii = ncm(1); deltc1(ii,:,:,2) = Dh(ii,:,:,2)
    jj = 1;      deltc1(:,jj,:,3) = Dh(:,jj,:,3)
    jj = ncm(2); deltc1(:,jj,:,4) = Dh(:,jj,:,4)
    kk = 1;      deltc1(:,:,kk,5) = Dh(:,:,kk,5)
    kk = ncm(3); deltc1(:,:,kk,6) = Dh(:,:,kk,6)

    ! cell components
    do kk = 1, ncm(3)
    do jj = 1, ncm(2)
    do ii = 1, ncm(1)
        ! conventional FDM
        if ( kk /= 1 )      Mcm(ii,jj,kk,1) = &
            -(deltc0(ii,jj,kk,5)+deltc1(ii,jj,kk,5))/(dfm(3)*fcz)
        if ( jj /= 1 )      Mcm(ii,jj,kk,2) = &
            -(deltc0(ii,jj,kk,3)+deltc1(ii,jj,kk,3))/(dfm(2)*fcr)
        if ( ii /= 1 )      Mcm(ii,jj,kk,3) = &
            -(deltc0(ii,jj,kk,1)+deltc1(ii,jj,kk,1))/(dfm(1)*fcr)
        if ( ii /= ncm(1) ) Mcm(ii,jj,kk,5) = &
            -(deltc0(ii,jj,kk,2)-deltc1(ii,jj,kk,2))/(dfm(1)*fcr)
        if ( jj /= ncm(2) ) Mcm(ii,jj,kk,6) = &
            -(deltc0(ii,jj,kk,4)-deltc1(ii,jj,kk,4))/(dfm(2)*fcr)
        if ( kk /= ncm(3) ) Mcm(ii,jj,kk,7) = &
            -(deltc0(ii,jj,kk,6)-deltc1(ii,jj,kk,6))/(dfm(3)*fcz)
        
        Mcm(ii,jj,kk,4)= &
            +(deltc0(ii,jj,kk,1)-deltc1(ii,jj,kk,1))/(dfm(1)*fcr) &
            +(deltc0(ii,jj,kk,2)+deltc1(ii,jj,kk,2))/(dfm(1)*fcr) &
            +(deltc0(ii,jj,kk,3)-deltc1(ii,jj,kk,3))/(dfm(2)*fcr) &
            +(deltc0(ii,jj,kk,4)+deltc1(ii,jj,kk,4))/(dfm(2)*fcr) &
            +(deltc0(ii,jj,kk,5)-deltc1(ii,jj,kk,5))/(dfm(3)*fcz) &
            +(deltc0(ii,jj,kk,6)+deltc1(ii,jj,kk,6))/(dfm(3)*fcz) &
            +cm_a(ii,jj,kk)

    end do
    end do
    end do

end subroutine


! =============================================================================
! G_INJ
! =============================================================================
subroutine G_INJ(Dt,Dh,phi)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:), phi(:,:,:)
    real(8):: deno

    ! incoming partial current for FMFD boundary condition
    do kk = 1, ncm(3)
    do jj = 1, ncm(2)
    do ii = 1, ncm(1)
        ! x-direction ---------------------------------------------------------
        if ( ii /= 1 ) then
        deno = Dt(ii-1,jj,kk,2)+Dt(ii,jj,kk,1)+Dh(ii,jj,kk,1)-Dh(ii-1,jj,kk,2)
        cmF(ii,jj,kk,1) = ((Dt(ii,jj,kk,1)-Dh(ii,jj,kk,1))*phi(ii,jj,kk) &
            +(Dt(ii-1,jj,kk,2)+Dh(ii-1,jj,kk,2))*phi(ii-1,jj,kk))/deno
        cmF(ii-1,jj,kk,2) = cmF(ii,jj,kk,1)

        cmJn(ii,jj,kk,1) = -Dt(ii,jj,kk,1)*(phi(ii,jj,kk)-cmF(ii,jj,kk,1)) &
                           +Dh(ii,jj,kk,1)*(phi(ii,jj,kk)+cmF(ii,jj,kk,1))
        cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)

        cmJ0(ii-1,jj,kk,2) = 25D-2*cmF(ii-1,jj,kk,2)-5D-1*cmJn(ii-1,jj,kk,2)
        cmJ1(ii,jj,kk,1)   = 25D-2*cmF(ii,jj,kk,1)+5D-1*cmJn(ii,jj,kk,1)
        end if
        ! y-direction ---------------------------------------------------------
        if ( jj /= 1 ) then
        deno = Dt(ii,jj-1,kk,4)+Dt(ii,jj,kk,3)+Dh(ii,jj,kk,3)-Dh(ii,jj-1,kk,4)
        cmF(ii,jj,kk,3) = ((Dt(ii,jj,kk,3)-Dh(ii,jj,kk,3))*phi(ii,jj,kk) &
            +(Dt(ii,jj-1,kk,4)+Dh(ii,jj-1,kk,4))*phi(ii,jj-1,kk))/deno
        cmF(ii,jj-1,kk,4) = cmF(ii,jj,kk,3)

        cmJn(ii,jj,kk,3) = -Dt(ii,jj,kk,3)*(phi(ii,jj,kk)-cmF(ii,jj,kk,3)) &
                           +Dh(ii,jj,kk,3)*(phi(ii,jj,kk)+cmF(ii,jj,kk,3))
        cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)

        cmJ0(ii,jj-1,kk,4) = 25D-2*cmF(ii,jj-1,kk,4)-5D-1*cmJn(ii,jj-1,kk,4)
        cmJ1(ii,jj,kk,3)   = 25D-2*cmF(ii,jj,kk,3)+5D-1*cmJn(ii,jj,kk,3)
        end if
        ! z-direction ---------------------------------------------------------
        if ( kk /= 1 ) then
        deno = Dt(ii,jj,kk-1,6)+Dt(ii,jj,kk,5)+Dh(ii,jj,kk,5)-Dh(ii,jj,kk-1,6)
        cmF(ii,jj,kk,5) = ((Dt(ii,jj,kk,5)-Dh(ii,jj,kk,5))*phi(ii,jj,kk) &
            +(Dt(ii,jj,kk-1,6)+Dh(ii,jj,kk-1,6))*phi(ii,jj,kk-1))/deno
        cmF(ii,jj,kk-1,6) = cmF(ii,jj,kk,5)

        cmJn(ii,jj,kk,5) = -Dt(ii,jj,kk,5)*(phi(ii,jj,kk)-cmF(ii,jj,kk,5)) &
                           +Dh(ii,jj,kk,5)*(phi(ii,jj,kk)+cmF(ii,jj,kk,5))
        cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)

        cmJ0(ii,jj,kk-1,6) = 25D-2*cmF(ii,jj,kk-1,6)-5D-1*cmJn(ii,jj,kk-1,6)
        cmJ1(ii,jj,kk,5)   = 25D-2*cmF(ii,jj,kk,5)+5D-1*cmJn(ii,jj,kk,5)
        end if
    end do
    end do
    end do

end subroutine

! =============================================================================
! G2L carries out the flux and current modulation (from GLOBAL to LOCAL)
! =============================================================================
subroutine G2L(phi0,phi1,fmJ0,fmJ1)
    implicit none
    real(8), intent(inout):: phi1(:,:,:), fmJ0(:,:,:,:), fmJ1(:,:,:,:)
    real(8), intent(out)::   phi0(:,:,:)
    real(8):: ssum

    phi0(:,:,:) = phi1(:,:,:)

    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr

    ! flux modulation
    phi1(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz) = &
    phi1(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz) &
    /sum(phi1(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz)) &
    *(fcr*fcr*fcz)*cm_phi1(ii,jj,kk)

    ! partial current modulation
    ! x0
    ssum = 0;       id(1) = id0(1)+1
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum = ssum + fmJ1(id(1),id(2),id(3),1)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum /= 0 ) fmJ1(id(1),id(2),id(3),1) = &
            fmJ1(id(1),id(2),id(3),1)/ssum*(fcr*fcz)*cmJ1(ii,jj,kk,1)
    end do
    end do
    ! x1
    ssum = 0;       id(1) = id0(1)+fcr
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum = ssum + fmJ0(id(1),id(2),id(3),2)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum /= 0 ) fmJ0(id(1),id(2),id(3),2) = &
            fmJ0(id(1),id(2),id(3),2)/ssum*(fcr*fcz)*cmJ0(ii,jj,kk,2)
    end do
    end do
    ! y0
    ssum = 0;       id(2) = id0(2)+1
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum = ssum + fmJ1(id(1),id(2),id(3),3)
    end do 
    end do 
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        if ( ssum /= 0 ) fmJ1(id(1),id(2),id(3),3) = &
            fmJ1(id(1),id(2),id(3),3)/ssum*(fcr*fcz)*cmJ1(ii,jj,kk,3)
    end do
    end do 
    ! y1
    ssum = 0;       id(2) = id0(2)+fcr
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum = ssum + fmJ0(id(1),id(2),id(3),4)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        if ( ssum /= 0 ) fmJ0(id(1),id(2),id(3),4) = &
            fmJ0(id(1),id(2),id(3),4)/ssum*(fcr*fcz)*cmJ0(ii,jj,kk,4)
    end do
    end do
    ! z0
    ssum = 0;       id(3) = id0(3)+1
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum = ssum + fmJ1(id(1),id(2),id(3),5)
    end do
    end do
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum /= 0 ) fmJ1(id(1),id(2),id(3),5) = &
            fmJ1(id(1),id(2),id(3),5)/ssum*(fcr*fcr)*cmJ1(ii,jj,kk,5)
    end do
    end do
    ! z1
    ssum = 0;       id(3) = id0(3)+fcz
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum = ssum + fmJ0(id(1),id(2),id(3),6)
    end do
    end do
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum /= 0 ) fmJ0(id(1),id(2),id(3),6) = &
            fmJ0(id(1),id(2),id(3),6)/ssum*(fcr*fcr)*cmJ0(ii,jj,kk,6)
    end do
    end do

    end do
    end do
    end do

end subroutine


! =============================================================================
! L_SOURCE
! =============================================================================
subroutine L_SOURCE(phi0,phi1,keff,fm_nf,fm_s,fmJ0,fmJ1)
    implicit none
    real(8), intent(inout):: fm_s(:,:,:)
    real(8), intent(in):: phi0(:,:,:), phi1(:,:,:), fm_nf(:,:,:), keff
    real(8), intent(in):: fmJ0(:,:,:,:), fmJ1(:,:,:,:)
    real(8):: fsource

    ! neutron source (fission + scattering)
    fm_s(:,:,:) = fm_nf(:,:,:)*phi1(:,:,:)/keff

    ! interface BC
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( ii /= 1 ) then; id(1) = id0(1)+1           ! x0
            fm_s(id(1),:,:) = fm_s(id(1),:,:) &
                +(jsrc(id(1),:,:,1)*fmJ1(id(1),:,:,1) &
                +fsrc(id(1),:,:,1)*phi0(id(1)-1,:,:))/dfm(1)
        end if
        if ( ii /= ncm(1) ) then; id(1) = id0(1)+fcr    ! x1
            fm_s(id(1),:,:) = fm_s(id(1),:,:) & 
                +(jsrc(id(1),:,:,2)*fmJ0(id(1),:,:,2) &
                -fsrc(id(1),:,:,2)*phi0(id(1)+1,:,:))/dfm(1)
        end if
    end do
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
        if ( jj /= 1 ) then; id(2) = id0(2)+1           ! y0
            fm_s(:,id(2),:) = fm_s(:,id(2),:) &
                +(jsrc(:,id(2),:,3)*fmJ1(:,id(2),:,3) &
                +fsrc(:,id(2),:,3)*phi0(:,id(2)-1,:))/dfm(2)
        end if
        if ( jj /= ncm(2) ) then; id(2) = id0(2)+fcr    ! y1
            fm_s(:,id(2),:) = fm_s(:,id(2),:) &
                +(jsrc(:,id(2),:,4)*fmJ0(:,id(2),:,4) &
                -fsrc(:,id(2),:,4)*phi0(:,id(2)+1,:))/dfm(2)
        end if
    end do
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        if ( kk /= 1 ) then; id(3) = id0(3)+1           ! z0
            fm_s(:,:,id(3)) = fm_s(:,:,id(3)) &
                +(jsrc(:,:,id(3),5)*fmJ1(:,:,id(3),5) &
                +fsrc(:,:,id(3),5)*phi0(:,:,id(3)-1))/dfm(3)
        end if
        if ( kk /= ncm(3) ) then; id(3) = id0(3)+fcz    ! z1
            fm_s(:,:,id(3)) = fm_s(:,:,id(3)) &
                +(jsrc(:,:,id(3),6)*fmJ0(:,:,id(3),6) &
                -fsrc(:,:,id(3),6)*phi0(:,:,id(3)+1))/dfm(3)
        end if
    end do

end subroutine

! =============================================================================
! L_OUTJ
! =============================================================================
subroutine L_OUTJ(phi0,phi1,fmF,fmJ0,fmJ1,fmJn)
    implicit none
    real(8), intent(in):: phi0(:,:,:), phi1(:,:,:)
    real(8), intent(inout):: fmF(:,:,:,:), fmJ0(:,:,:,:)
    real(8), intent(inout):: fmJ1(:,:,:,:), fmJn(:,:,:,:)
    real(8):: netJ(nfm(1),nfm(2),nfm(3))

    ! outgoing partial current
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( ii /= 1 ) then; id(1) = id0(1)+1           ! x0
            fmJn(id(1),:,:,1) = +jsrc(id(1),:,:,1)*fmJ1(id(1),:,:,1) &
                +deltf1(id(1),:,:,1)*phi1(id(1),:,:) &
                +fsrc(id(1),:,:,1)*phi0(id(1)-1,:,:)
            fmF(id(1),:,:,1) = 4D0*fmJ1(id(1),:,:,1)-2D0*fmJn(id(1),:,:,1)
            fmJ0(id(1),:,:,1) = 25D-2*fmF(id(1),:,:,1)-5D-1*fmJn(id(1),:,:,1)
        end if
        if ( ii /= ncm(1) ) then; id(1) = id0(1)+fcr    ! x1
            fmJn(id(1),:,:,2) = -jsrc(id(1),:,:,2)*fmJ0(id(1),:,:,2) &
                +deltf1(id(1),:,:,2)*phi1(id(1),:,:) &
                +fsrc(id(1),:,:,2)*phi0(id(1)+1,:,:)
            fmF(id(1),:,:,2) = 4D0*fmJ0(id(1),:,:,2)+2D0*fmJn(id(1),:,:,2)
            fmJ1(id(1),:,:,2) = 25D-2*fmF(id(1),:,:,2)+5D-1*fmJn(id(1),:,:,2)
        end if
    end do
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
        if ( jj /= 1 ) then; id(2) = id0(2)+1           ! y0
            fmJn(:,id(2),:,3) = +jsrc(:,id(2),:,3)*fmJ1(:,id(2),:,3) &
                +deltf1(:,id(2),:,3)*phi1(:,id(2),:) &
                +fsrc(:,id(2),:,3)*phi0(:,id(2)-1,:)
            fmF(:,id(2),:,3) = 4D0*fmJ1(:,id(2),:,3)-2D0*fmJn(:,id(2),:,3)
            fmJ0(:,id(2),:,3) = 25D-2*fmF(:,id(2),:,3)-5D-1*fmJn(:,id(2),:,3)
        end if
        if ( jj /= ncm(2) ) then; id(2) = id0(2)+fcr    ! y1
            fmJn(:,id(2),:,4) = -jsrc(:,id(2),:,4)*fmJ0(:,id(2),:,4) &
                +deltf1(:,id(2),:,4)*phi1(:,id(2),:) &
                +fsrc(:,id(2),:,4)*phi0(:,id(2)+1,:)
            fmF(:,id(2),:,4) = 4D0*fmJ0(:,id(2),:,4)+2D0*fmJn(:,id(2),:,4)
            fmJ1(:,id(2),:,4) = 25D-2*fmF(:,id(2),:,4)+5D-1*fmJn(:,id(2),:,4)
        end if
    end do
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        if ( kk /= 1 ) then; id(3) = id0(3)+1           ! z0
            fmJn(:,:,id(3),5) = +jsrc(:,:,id(3),5)*fmJ1(:,:,id(3),5) &
                +deltf1(:,:,id(3),5)*phi1(:,:,id(3)) &
                +fsrc(:,:,id(3),5)*phi0(:,:,id(3)-1)
            fmF(:,:,id(3),5) = 4D0*fmJ1(:,:,id(3),5)-2D0*fmJn(:,:,id(3),5)
            fmJ0(:,:,id(3),5) = 25D-2*fmF(:,:,id(3),5)-5D-1*fmJn(:,:,id(3),5)
        end if
        if ( kk /= ncm(3) ) then; id(3) = id0(3)+fcz    ! z1
            fmJn(:,:,id(3),6) = -jsrc(:,:,id(3),6)*fmJ0(:,:,id(3),6) &
                +deltf1(:,:,id(3),6)*phi1(:,:,id(3)) &
                +fsrc(:,:,id(3),6)*phi0(:,:,id(3)+1)
            fmF(:,:,id(3),6) = 4D0*fmJ0(:,:,id(3),6)+2D0*fmJn(:,:,id(3),6)
            fmJ1(:,:,id(3),6) = 25D-2*fmF(:,:,id(3),6)+5D-1*fmJn(:,:,id(3),6)
        end if
    end do

    ! boundary surface
    ii = 1;      fmJn(ii,:,:,1) = deltf1(ii,:,:,1)*phi1(ii,:,:)
    ii = nfm(1); fmJn(ii,:,:,2) = deltf1(ii,:,:,2)*phi1(ii,:,:)
    jj = 1;      fmJn(:,jj,:,3) = deltf1(:,jj,:,3)*phi1(:,jj,:)
    jj = nfm(2); fmJn(:,jj,:,4) = deltf1(:,jj,:,4)*phi1(:,jj,:)
    kk = 1;      fmJn(:,:,kk,5) = deltf1(:,:,kk,5)*phi1(:,:,kk)
    kk = nfm(3); fmJn(:,:,kk,6) = deltf1(:,:,kk,6)*phi1(:,:,kk)

    ! data swapping
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        ! x0
        if ( ii /= 1 ) then
        id(1) = id0(1)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ0(id(1)-1,id(2),id(3),2) = fmJ0(id(1),id(2),id(3),1)
            !fmJn(id(1),id(2),id(3),1) = fmJ1(id(1)-1,id(2),id(3),2)-fmJ0(id(1),id(2),id(3),1)
        end do
        end do
        end if
        ! x1
        if ( ii /= ncm(1) ) then
        id(1) = id0(1)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ1(id(1)+1,id(2),id(3),1) = fmJ1(id(1),id(2),id(3),2)
            !fmJn(id(1),id(2),id(3),2) = fmJ1(id(1),id(2),id(3),2)-fmJ0(id(1)+1,id(2),id(3),1)
        end do
        end do
        end if
        ! y0
        if ( jj /= 1 ) then
        id(2) = id0(2)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJ0(id(1),id(2)-1,id(3),4) = fmJ0(id(1),id(2),id(3),3)
            !fmJn(id(1),id(2),id(3),3) = fmJ1(id(1),id(2)-1,id(3),4)-fmJ0(id(1),id(2),id(3),3)
        end do
        end do
        end if
        ! y1
        if ( jj /= ncm(2) ) then
        id(2) = id0(2)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJ1(id(1),id(2)+1,id(3),3) = fmJ1(id(1),id(2),id(3),4)
            !fmJn(id(1),id(2),id(3),4) = fmJ1(id(1),id(2),id(3),4)-fmJ0(id(1),id(2)+1,id(3),3)
        end do
        end do
        end if
        ! z0
        if ( kk /= 1 ) then
        id(3) = id0(3)+1
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ0(id(1),id(2),id(3)-1,6) = fmJ0(id(1),id(2),id(3),5)
            !fmJn(id(1),id(2),id(3),5) = fmJ1(id(1),id(2),id(3)-1,6)-fmJ0(id(1),id(2),id(3),5)
        end do
        end do
        end if
        ! z1
        if ( kk /= ncm(3) ) then
        id(3) = id0(3)+fcz
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ1(id(1),id(2),id(3)+1,5) = fmJ1(id(1),id(2),id(3),6)
            !fmJn(id(1),id(2),id(3),6) = fmJ1(id(1),id(2),id(3),6)-fmJ0(id(1),id(2)+1,id(3)+1,5)
        end do
        end do
        end if
    end do
    end do
    end do


!    write(8,1), fmJn(:,:,:,1)
!    write(8,*)
!    write(8,1), fmJn(:,:,:,2)
!    write(8,*)
!    write(8,1), fmJn(:,:,:,3)
!    write(8,*)
!    write(8,1), fmJn(:,:,:,4)
!    write(8,*)
!    1 format(20es15.7)
!    stop



end subroutine


! =============================================================================
! L_REFJ
! =============================================================================
subroutine L_REFJ(fmF,fmJ0,fmJ1,fmJn)
    implicit none
    real(8), intent(in), dimension(:,:,:,:):: fmF, fmJ0, fmJ1, fmJn
    real(8):: ssum(2)

!    ! surface average
!    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
!    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
!    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
!        ! x0
!        if ( ii /= 1 ) then
!        ssum = 0;       id(1) = id0(1)+1
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),1)
!            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),1)
!        end do
!        end do
!        cmJ0(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
!        cmF(ii,jj,kk,1)  = ssum(2) / (fcr*fcz)
!        end if
!        ! x1
!        if ( ii /= ncm(1) ) then
!        ssum = 0;       id(1) = id0(1)+fcr
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),2)
!            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),2)
!        end do
!        end do
!        cmJ1(ii,jj,kk,2) = ssum(1) / (fcr*fcz)
!        cmF(ii,jj,kk,2)  = ssum(2) / (fcr*fcz)
!        end if
!        ! y0
!        if ( jj /= 1 ) then
!        ssum = 0;       id(2) = id0(2)+1
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do mm = 1, fcr; id(1) = id0(1)+mm
!            ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),3)
!            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),3)
!        end do
!        end do
!        cmJ0(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
!        cmF(ii,jj,kk,3)  = ssum(2) / (fcr*fcz)
!        end if
!        ! y1
!        if ( jj /= ncm(2) ) then
!        ssum = 0;       id(2) = id0(2)+fcr
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do mm = 1, fcr; id(1) = id0(1)+mm
!            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),4)
!            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),4)
!        end do
!        end do
!        cmJ1(ii,jj,kk,4) = ssum(1) / (fcr*fcz)
!        cmF(ii,jj,kk,4)  = ssum(2) / (fcr*fcz)
!        end if
!        ! z0
!        if ( kk /= 1 ) then
!        ssum = 0;       id(3) = id0(3)+1
!        do mm = 1, fcr; id(1) = id0(1)+mm
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),5)
!            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),5)
!        end do
!        end do
!        cmJ0(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
!        cmF(ii,jj,kk,5)  = ssum(2) / (fcr*fcr)
!        end if
!        ! z1
!        if ( kk /= ncm(3) ) then
!        ssum = 0;       id(3) = id0(3)+fcz
!        do mm = 1, fcr; id(1) = id0(1)+mm
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),6)
!            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),6)
!        end do
!        end do
!        cmJ1(ii,jj,kk,6) = ssum(1) / (fcr*fcr)
!        cmF(ii,jj,kk,6)  = ssum(2) / (fcr*fcr)
!        end if
!    end do
!    end do
!    end do
!
!    ! interface surface
!    do ii = 1, ncm(1)
!    do jj = 1, ncm(2)
!    do kk = 1, ncm(3)
!        ! x-direction
!        if ( ii /= 1 ) then
!            cmJn(ii,jj,kk,1) = cmJ1(ii-1,jj,kk,2)-cmJ0(ii,jj,kk,1)
!            cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)
!            cmF(ii,jj,kk,1) = (cmF(ii,jj,kk,1)+cmF(ii-1,jj,kk,2))/2D0
!            cmF(ii-1,jj,kk,2) = cmF(ii,jj,kk,1)
!        end if
!        ! y-direction
!        if ( jj /= 1 ) then
!            cmJn(ii,jj,kk,3) = cmJ1(ii,jj-1,kk,4)-cmJ0(ii,jj,kk,3)
!            cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)
!            cmF(ii,jj,kk,3) = (cmF(ii,jj,kk,3)+cmF(ii,jj-1,kk,4))/2D0
!            cmF(ii,jj-1,kk,4) = cmF(ii,jj,kk,3)
!        end if
!        ! z-direction
!        if ( kk /= 1 ) then
!            cmJn(ii,jj,kk,5) = cmJ1(ii,jj,kk-1,6)-cmJ0(ii,jj,kk,5)
!            cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)
!            cmF(ii,jj,kk,5) = (cmF(ii,jj,kk,5)+cmF(ii,jj,kk-1,6))/2D0
!            cmF(ii,jj,kk-1,6) = cmF(ii,jj,kk,5)
!        end if
!    end do
!    end do
!    end do


    ! surface average
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        ! x0
        if ( ii /= 1 ) then
        ssum = 0;       id(1) = id0(1)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),1)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),1)
        end do
        end do
        cmJn(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,1)  = ssum(2) / (fcr*fcz)
        end if
        ! x1
        if ( ii /= ncm(1) ) then
        ssum = 0;       id(1) = id0(1)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),2)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),2)
        end do
        end do
        cmJn(ii,jj,kk,2) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,2)  = ssum(2) / (fcr*fcz)
        end if
        ! y0
        if ( jj /= 1 ) then
        ssum = 0;       id(2) = id0(2)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),3)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),3)
        end do
        end do
        cmJn(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,3)  = ssum(2) / (fcr*fcz)
        end if
        ! y1
        if ( jj /= ncm(2) ) then
        ssum = 0;       id(2) = id0(2)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),4)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),4)
        end do
        end do
        cmJn(ii,jj,kk,4) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,4)  = ssum(2) / (fcr*fcz)
        end if
        ! z0
        if ( kk /= 1 ) then
        ssum = 0;       id(3) = id0(3)+1
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),5)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),5)
        end do
        end do
        cmJn(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,5)  = ssum(2) / (fcr*fcr)
        end if
        ! z1
        if ( kk /= ncm(3) ) then
        ssum = 0;       id(3) = id0(3)+fcz
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),6)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),6)
        end do
        end do
        cmJn(ii,jj,kk,6) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,6)  = ssum(2) / (fcr*fcr)
        end if
    end do
    end do
    end do

    ! interface surface
    do ii = 1, ncm(1)
    do jj = 1, ncm(2)
    do kk = 1, ncm(3)
        ! x-direction
        if ( ii /= 1 ) then
            cmJn(ii,jj,kk,1) = (cmJn(ii,jj,kk,1)+cmJn(ii-1,jj,kk,2))/2D0
            cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)
            cmF(ii,jj,kk,1) = (cmF(ii,jj,kk,1)+cmF(ii-1,jj,kk,2))/2D0
            cmF(ii-1,jj,kk,2) = cmF(ii,jj,kk,1)
        end if
        ! y-direction
        if ( jj /= 1 ) then
            cmJn(ii,jj,kk,3) = (cmJn(ii,jj,kk,3)+cmJn(ii,jj-1,kk,4))/2D0
            cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)
            cmF(ii,jj,kk,3) = (cmF(ii,jj,kk,3)+cmF(ii,jj-1,kk,4))/2D0
            cmF(ii,jj-1,kk,4) = cmF(ii,jj,kk,3)
        end if
        ! z-direction
        if ( kk /= 1 ) then
            cmJn(ii,jj,kk,5) = (cmJn(ii,jj,kk,5)+cmJn(ii,jj,kk-1,6))/2D0
            cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)
            cmF(ii,jj,kk,5) = (cmF(ii,jj,kk,5)+cmF(ii,jj,kk-1,6))/2D0
            cmF(ii,jj,kk-1,6) = cmF(ii,jj,kk,5)
        end if
    end do
    end do
    end do


    ! boundary surface
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    !   x0
    ii = 1; id(1) = 1
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),1)
    end do
    end do
    cmJn(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
    end do
    !   x1
    ii = ncm(1); id(1) = nfm(1)
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),2)
    end do
    end do
    cmJn(ii,jj,kk,2) = ssum(1) / (fcr*fcz)
    end do
    !   y0
    jj = 1; id(2) = 1
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum(1) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),3)
    end do
    end do
    cmJn(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
    end do
    !   y1
    jj = ncm(2); id(2) = nfm(2)
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum(1) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),4)
    end do
    end do
    cmJn(ii,jj,kk,4) = ssum(1) / (fcr*fcz)
    end do
    end do
    !   z0
    kk = 1; id(3) = 1
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1) = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),5)
    end do
    end do
    cmJn(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
    end do
    end do
    !   z1
    kk = ncm(3); id(3) = nfm(3)
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1) = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),6)
    end do
    end do
    cmJn(ii,jj,kk,6) = ssum(1) / (fcr*fcr)
    end do
    end do

end subroutine
    
    
! =============================================================================
! G_XS produces the flux-volume-weight group constants
! ============================================================================= 
subroutine G_XS(fm_t,fm_a,fm_nf,phi)
    implicit none
    real(8), intent(in), dimension(:,:,:):: fm_t, fm_a, fm_nf, phi

    ! homogenization
    do ii = 1, ncm(1); id(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id(2) = (jj-1)*fcr
    do kk = 1, ncm(3); id(3) = (kk-1)*fcz
        cm_phi1(ii,jj,kk) = sum(phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))
        cm_t(ii,jj,kk) = sum(fm_t(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz)*phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))/cm_phi1(ii,jj,kk)
        cm_a(ii,jj,kk) = sum(fm_a(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz)*phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))/cm_phi1(ii,jj,kk)
        cm_nf(ii,jj,kk) = sum(fm_nf(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz)*phi(id(1)+1:id(1)+fcr, &
            id(2)+1:id(2)+fcr,id(3)+1:id(3)+fcz))/cm_phi1(ii,jj,kk)
    end do
    end do
    end do
    cmD = 1D0 / (3D0 * cm_t)
    where ( cm_t == 0 ) cmD = 0
    cm_phi1 = cm_phi1 / (fcr*fcr*fcz)

    ! interface diffusion coefficient
    do ii = 1, ncm(1)
    do jj = 1, ncm(2)
    do kk = 1, ncm(3)
        cmDt(ii,jj,kk,1) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(1))
        cmDt(ii,jj,kk,2) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(1))
        cmDt(ii,jj,kk,3) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(2))
        cmDt(ii,jj,kk,4) = 2D0*cmD(ii,jj,kk)/(fcr*dfm(2))
        cmDt(ii,jj,kk,5) = 2D0*cmD(ii,jj,kk)/(fcz*dfm(3))
        cmDt(ii,jj,kk,6) = 2D0*cmD(ii,jj,kk)/(fcz*dfm(3))
    end do
    end do
    end do

end subroutine

end module

