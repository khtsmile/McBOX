module CMFD
    use FMFD_HEADER
    implicit none


    contains


! =============================================================================
! ONE_NODE_CMFD
! =============================================================================
subroutine ONE_NODE_CMFD(keff,fm_t,fm_a,fm_nf,fmD,fm_phi1,fmJ0,fmJ1,fmJn,fmF)
    use SOLVERS, only: BICG_G, SORL, BiCG_L, CG1
    implicit none
    real(8), intent(in):: keff
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

    k_fmfd = keff
    call L2G(fm_phi1,fm_t,fm_a,fm_nf,fmJ0,fmJ1,fmF)
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
    

    k_pre = k_fmfd
    do global = 1, 5
    cm_phi0 = cm_phi1
    !call G_SOURCE
    cm_s = cm_nf*cm_phi0/k_fmfd
    cm_phi1 = BiCG_G(Mcm,cm_s)
    !call CG1(Mcm,cm_s,cm_phi1)
    !call G_POWER
    k_fmfd = k_fmfd*sum(cm_nf*cm_phi1*cm_nf*cm_phi1) &
           / sum(cm_nf*cm_phi0*cm_nf*cm_phi1)
    end do
    error = abs(k_fmfd-k_pre)/k_fmfd
    !print*, k_fmfd, error
    if ( error < 1D-8 .or. isnan(k_fmfd) ) exit
    ! ------------------------------- LOCAL
    call G_INJ(cmDt,cmDh,cm_phi1)


    do local = 1, 2
    call G2L(fm_phi0,fm_phi1,fmJ0,fmJ1)
    call L_SOURCE(fm_phi0,fm_phi1,keff,fm_nf,fm_s,fmJ0,fmJ1)

!    write(8,1), Mfm(:,:,:,1)
!    write(8,*)
!    write(8,1), Mfm(:,:,:,2)
!    write(8,*)
!    write(8,1), Mfm(:,:,:,3)
!    write(8,*)
!    write(8,1), Mfm(:,:,:,4)
!    write(8,*)
!    write(8,1), Mfm(:,:,:,5)
!    write(8,*)
!    write(8,1), Mfm(:,:,:,6)
!    write(8,*)
!    write(8,1), Mfm(:,:,:,7)
!    write(8,*)
!    1 format(20es15.7)
!    stop


!    fm_phi1(:,:,:) = BICG_L(Mfm(:,:,:,:),fm_s(:,:,:))
    call SORL(Mfm,fm_s,fm_phi1)
    
!    print*, fm_phi1
!    stop

    call L_OUTJ(fm_phi0,fm_phi1,fmF,fmJ0,fmJ1,fmJn)
    call L_REFJ(fmF,fmJ0,fmJ1,fmJn)
    call G_XS(fm_t,fm_a,fm_nf,fm_phi1)
    end do
    end do



end subroutine

! =============================================================================
! L2G homogenizes the reactor parameters from local to global
! =============================================================================
subroutine L2G(phi,fm_tot,fm_abso,fm_nufiss,fmJ0,fmJ1,fmF)
    implicit none
    real(8), intent(in), dimension(:,:,:):: phi, fm_tot, fm_abso, fm_nufiss
    real(8), intent(in), dimension(:,:,:,:):: fmJ0, fmJ1, fmF
    real(8):: ssum(0:2)

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
            ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),1)
            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),1)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),1)
        end do
        end do
        cmJ0(ii,jj,kk,1) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,1)  = ssum(2) / (fcr*fcz)
        if ( ii /= 1 ) then
        cmJ0(ii-1,jj,kk,2) = cmJ0(ii,jj,kk,1)
        cmJ1(ii-1,jj,kk,2) = cmJ1(ii,jj,kk,1)
        cmF(ii-1,jj,kk,2)  = cmF(ii,jj,kk,1)
        end if
        ! y-direction
        ssum = 0;       id(2) = id0(2)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),3)
            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),3)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),3)
        end do
        end do
        cmJ0(ii,jj,kk,3) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,3)  = ssum(2) / (fcr*fcz)
        if ( jj /= 1 ) then
        cmJ0(ii,jj-1,kk,4) = cmJ0(ii,jj,kk,3)
        cmJ1(ii,jj-1,kk,4) = cmJ1(ii,jj,kk,3)
        cmF(ii,jj-1,kk,4)  = cmF(ii,jj,kk,3)
        end if
        ! z-direction
        ssum = 0;       id(3) = id0(3)+1
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),5)
            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),5)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),5)
        end do
        end do
        cmJ0(ii,jj,kk,5) = ssum(0) / (fcr*fcr)
        cmJ1(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,5)  = ssum(2) / (fcr*fcr)
        if ( kk /= 1 ) then
        cmJ0(ii,jj,kk-1,6) = cmJ0(ii,jj,kk,5)
        cmJ1(ii,jj,kk-1,6) = cmJ1(ii,jj,kk,5)
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
        ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),2)
        ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),2)
        ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),2)
    end do
    end do
    cmJ0(ii,jj,kk,2) = ssum(0) / (fcr*fcz)
    cmJ1(ii,jj,kk,2) = ssum(1) / (fcr*fcz)
    cmF(ii,jj,kk,2)  = ssum(2) / (fcr*fcz)
    end do
    end do
    !   y-direction
    jj = ncm(2);       id(2) = nfm(2)
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),4)
        ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),4)
        ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),4)
    end do
    end do
    cmJ0(ii,jj,kk,4) = ssum(0) / (fcr*fcz)
    cmJ1(ii,jj,kk,4) = ssum(1) / (fcr*fcz)
    cmF(ii,jj,kk,4)  = ssum(2) / (fcr*fcz)
    end do
    end do
    !   z-direction
    kk = ncm(3);       id(3) = nfm(3)
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),6)
        ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),6)
        ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),6)
    end do
    end do
    cmJ0(ii,jj,kk,6) = ssum(0) / (fcr*fcr)
    cmJ1(ii,jj,kk,6) = ssum(1) / (fcr*fcr)
    cmF(ii,jj,kk,6)  = ssum(2) / (fcr*fcr)
    end do
    end do
    cmJn = cmJ1 - cmJ0

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
        ! x-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            id(1) = id0(1)+1
            Dt(id(1),id(2),id(3),1) = 2D0*D(id(1),id(2),id(3))/dfm(1)
            deltf0(id(1),id(2),id(3),1) = 0D0
            id(1) = id0(1)+fcr
            Dt(id(1),id(2),id(3),2) = 2D0*D(id(1),id(2),id(3))/dfm(1)
            deltf0(id(1),id(2),id(3),2) = 0D0
        end do
        ! y-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
            id(2) = id0(2)+1
            Dt(id(1),id(2),id(3),3) = 2D0*D(id(1),id(2),id(3))/dfm(2)
            deltf0(id(1),id(2),id(3),3) = 0D0
            id(2) = id0(2)+fcr
            Dt(id(1),id(2),id(3),4) = 2D0*D(id(1),id(2),id(3))/dfm(2)
            deltf0(id(1),id(2),id(3),4) = 0D0
        end do
        end do
        ! z-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            id(3) = id0(3)+1
            Dt(id(1),id(2),id(3),5) = 2D0*D(id(1),id(2),id(3))/dfm(3)
            deltf0(id(1),id(2),id(3),5) = 0D0
            id(3) = id0(3)+fcz
            Dt(id(1),id(2),id(3),6) = 2D0*D(id(1),id(2),id(3))/dfm(3)
            deltf0(id(1),id(2),id(3),6) = 0D0
        end do
        end do
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

    ! inner boundary
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr

        ! x-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ! x0 --------------------------------------------------------------
            if ( ii /= 1 ) then
            id(1) = id0(1)+1
            deltf1(id(1),id(2),id(3),1) = -(Dt(id(1),id(2),id(3),1) &
                -Dh(id(1),id(2),id(3),1))/(1D0+2D0*Dt(id(1),id(2),id(3),1))
            end if
            ! x1 --------------------------------------------------------------
            if ( ii /= ncm(1) ) then
            id(1) = id0(1)+fcr
            deltf1(id(1),id(2),id(3),2) = (Dt(id(1),id(2),id(3),2) &
                +Dh(id(1),id(2),id(3),2))/(1D0+2D0*Dt(id(1),id(2),id(3),2))
            end if
        end do
        end do
        ! y-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ! y0 --------------------------------------------------------------
            if ( jj /= 1 ) then
            id(2) = id0(2)+1
            deltf1(id(1),id(2),id(3),3) = -(Dt(id(1),id(2),id(3),3) &
                -Dh(id(1),id(2),id(3),3))/(1D0+2D0*Dt(id(1),id(2),id(3),3))
            end if
            ! y1 --------------------------------------------------------------
            if ( jj /= ncm(2) ) then
            id(2) = id0(2)+fcr
            deltf1(id(1),id(2),id(3),4) = (Dt(id(1),id(2),id(3),4) &
                +Dh(id(1),id(2),id(3),4))/(1D0+2D0*Dt(id(1),id(2),id(3),4))
            end if
        end do
        end do
        ! z-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ! z0 --------------------------------------------------------------
            if ( kk /= 1 ) then
            id(3) = id0(3)+1
            deltf1(id(1),id(2),id(3),5) = -(Dt(id(1),id(2),id(3),5) &
                -Dh(id(1),id(2),id(3),5))/(1D0+2D0*Dt(id(1),id(2),id(3),5))
            end if
            ! z1 --------------------------------------------------------------
            if ( kk /= ncm(3) ) then
            id(3) = id0(3)+fcz
            deltf1(id(1),id(2),id(3),6) = (Dt(id(1),id(2),id(3),6) &
                +Dh(id(1),id(2),id(3),6))/(1D0+2D0*Dt(id(1),id(2),id(3),6))
            end if
        end do
        end do
    end do
    end do
    end do

!    !!! zigzag !!! BC
!    if ( zigzag ) then
!        ! quarter core
!        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
!            ! -----
!            do kk = fm0(3), fm1(3)
!            ii = fm1(1)-zz(2)
!            do jj = fm1(2)-zz(1)+1, fm1(2)
!                deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
!            end do
!            ! -----
!            ii = fm1(1)-zz(1)
!            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
!                deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
!            end do
!            ! -----
!            jj = fm1(2)-zz(2)
!            do ii = fm1(1)-zz(1)+1, fm1(1)
!                deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
!            end do
!            ! -----
!            jj = fm1(2)-zz(1)
!            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
!                deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
!            end do
!            end do
!
!        ! whole core
!        else
!            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
!            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
!            
!            if ( id(1) == afm(1)/2-zz(2) .and. &
!                 afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
!                if ( nint(ii-mp(1)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,1) = fm_dhat(ii,jj,kk,:,1)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(1) == afm(1)/2-zz(1) .and. &
!                 afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
!                if ( nint(ii-mp(1)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,1) = fm_dhat(ii,jj,kk,:,1)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == afm(2)/2-zz(2) .and. &
!                 afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
!                if ( nint(jj-mp(2)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,3) = fm_dhat(ii,jj,kk,:,3)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == afm(2)/2-zz(1) .and. &
!                 afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
!                if ( nint(jj-mp(2)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,3) = fm_dhat(ii,jj,kk,:,3)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
!                    end do
!                end if
!            end if
!
!            end do
!            end do
!        end if
!    end if

end subroutine

! =============================================================================
! L_MATRIX
! =============================================================================
subroutine L_MATRIX(Dt,Dh,abso)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:), abso(:,:,:)
    real(8):: deno  ! denominator

    ! Matrix formulation
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    
    ! -------------------------------------------------------------------------
    ! migration term
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
    do oo = 1, fcz; id(3) = id0(3)+oo
    
    if ( oo /= 1   ) Mfm(id(1),id(2),id(3),1) = &
      -(deltf0(id(1),id(2),id(3),5)+deltf1(id(1),id(2),id(3),5))/dfm(3)
    if ( nn /= 1   ) Mfm(id(1),id(2),id(3),2) = &
      -(deltf0(id(1),id(2),id(3),3)+deltf1(id(1),id(2),id(3),3))/dfm(2)
    if ( mm /= 1   ) Mfm(id(1),id(2),id(3),3) = &
      -(deltf0(id(1),id(2),id(3),1)+deltf1(id(1),id(2),id(3),1))/dfm(1)
    if ( mm /= fcr ) Mfm(id(1),id(2),id(3),5) = &
      -(deltf0(id(1),id(2),id(3),2)-deltf1(id(1),id(2),id(3),2))/dfm(1)
    if ( nn /= fcr ) Mfm(id(1),id(2),id(3),6) = &
      -(deltf0(id(1),id(2),id(3),4)-deltf1(id(1),id(2),id(3),4))/dfm(2)
    if ( oo /= fcz ) Mfm(id(1),id(2),id(3),7) = &
      -(deltf0(id(1),id(2),id(3),6)-deltf1(id(1),id(2),id(3),6))/dfm(3)
    
    Mfm(id(1),id(2),id(3),4) = &
     +(deltf0(id(1),id(2),id(3),1)-deltf1(id(1),id(2),id(3),1))/dfm(1) &
     +(deltf0(id(1),id(2),id(3),2)+deltf1(id(1),id(2),id(3),2))/dfm(1) &
     +(deltf0(id(1),id(2),id(3),3)-deltf1(id(1),id(2),id(3),3))/dfm(2) &
     +(deltf0(id(1),id(2),id(3),4)+deltf1(id(1),id(2),id(3),4))/dfm(2) &
     +(deltf0(id(1),id(2),id(3),5)-deltf1(id(1),id(2),id(3),5))/dfm(3) &
     +(deltf0(id(1),id(2),id(3),6)+deltf1(id(1),id(2),id(3),6))/dfm(3) &
     +abso(id(1),id(2),id(3))

    end do
    end do
    end do

    ! -------------------------------------------------------------------------
    ! source term
    ! x-direction
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ! x0 --------------------------------------------------------------
        if ( ii /= 1 ) then
        id(1) = id0(1)+1
        deno = 1D0+2D0*Dt(id(1),id(2),id(3),1)
        jsrc(id(1),id(2),id(3),1) = 4D0*Dt(id(1),id(2),id(3),1)/deno
        fsrc(id(1),id(2),id(3),1) = Dh(id(1),id(2),id(3),1)/deno
        end if
        ! x1 --------------------------------------------------------------
        if ( ii /= ncm(1) ) then
        id(1) = id0(1)+fcr
        deno = 1D0+2D0*Dt(id(1),id(2),id(3),2)
        jsrc(id(1),id(2),id(3),2) = 4D0*Dt(id(1),id(2),id(3),2)/deno
        fsrc(id(1),id(2),id(3),2) = Dh(id(1),id(2),id(3),2)/deno
        end if
    end do
    end do
    ! y-direction
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ! y0 --------------------------------------------------------------
        if ( jj /= 1 ) then
        id(2) = id0(2)+1
        deno = 1D0+2D0*Dt(id(1),id(2),id(3),3)
        jsrc(id(1),id(2),id(3),3) = 4D0*Dt(id(1),id(2),id(3),3)/deno
        fsrc(id(1),id(2),id(3),3) = Dh(id(1),id(2),id(3),3)/deno
        end if
        ! y1 --------------------------------------------------------------
        if ( jj /= ncm(2) ) then
        id(2) = id0(2)+fcr
        deno = 1D0+2D0*Dt(id(1),id(2),id(3),4)
        jsrc(id(1),id(2),id(3),4) = 4D0*Dt(id(1),id(2),id(3),4)/deno
        fsrc(id(1),id(2),id(3),4) = Dh(id(1),id(2),id(3),4)/deno
        end if
    end do
    end do
    ! z-direction
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ! z0 --------------------------------------------------------------
        if ( kk /= 1 ) then
        id(3) = id0(3)+1
        deno = 1D0+2D0*Dt(id(1),id(2),id(3),5)
        jsrc(id(1),id(2),id(3),5) = 4D0*Dt(id(1),id(2),id(3),5)/deno
        fsrc(id(1),id(2),id(3),5) = Dh(id(1),id(2),id(3),5)/deno
        end if
        ! z1 --------------------------------------------------------------
        if ( kk /= ncm(3) ) then
        id(3) = id0(3)+fcz
        deno = 1D0+2D0*Dt(id(1),id(2),id(3),6)
        jsrc(id(1),id(2),id(3),6) = 4D0*Dt(id(1),id(2),id(3),6)/deno
        fsrc(id(1),id(2),id(3),6) = Dh(id(1),id(2),id(3),6)/deno
        end if
    end do
    end do

    end do
    end do
    end do

!    !!! zigzag !!!
!    if ( zigzag ) then
!        ! quarter core
!        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
!            ! -----
!            do kk = fm0(3), fm1(3)
!            ii = fm1(1)-zz(2)
!            do jj = fm1(2)-zz(1)+1, fm1(2)
!                jsrc(ii,jj,kk,:,2) = 0D0
!                fsrc(ii,jj,kk,:,2) = 0D0
!            end do
!            ! -----
!            ii = fm1(1)-zz(1)
!            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
!                jsrc(ii,jj,kk,:,2) = 0D0
!                fsrc(ii,jj,kk,:,2) = 0D0
!            end do
!            ! -----
!            jj = fm1(2)-zz(2)
!            do ii = fm1(1)-zz(1)+1, fm1(1)
!                jsrc(ii,jj,kk,:,4) = 0D0
!                fsrc(ii,jj,kk,:,4) = 0D0
!            end do
!            ! -----
!            jj = fm1(2)-zz(1)
!            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
!                jsrc(ii,jj,kk,:,4) = 0D0
!                fsrc(ii,jj,kk,:,4) = 0D0
!            end do
!            end do
!
!        ! whole core
!        else
!            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
!            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
!
!            if ( id(1) == afm(1)/2-zz(2) .and. &
!                 afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
!                if ( nint(ii-mp(1)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,1) = 0D0
!                    fsrc(ii,jj,kk,:,1) = 0D0
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,2) = 0D0
!                    fsrc(ii,jj,kk,:,2) = 0D0
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(1) == afm(1)/2-zz(1) .and. &
!                 afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
!                if ( nint(ii-mp(1)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,1) = 0D0
!                    fsrc(ii,jj,kk,:,1) = 0D0
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,2) = 0D0
!                    fsrc(ii,jj,kk,:,2) = 0D0
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == afm(2)/2-zz(2) .and. &
!                 afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
!                if ( nint(jj-mp(2)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,3) = 0D0
!                    fsrc(ii,jj,kk,:,3) = 0D0
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,4) = 0D0
!                    fsrc(ii,jj,kk,:,4) = 0D0
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == afm(2)/2-zz(1) .and. &
!                 afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
!                if ( nint(jj-mp(2)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,3) = 0D0
!                    fsrc(ii,jj,kk,:,3) = 0D0
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                    jsrc(ii,jj,kk,:,4) = 0D0
!                    fsrc(ii,jj,kk,:,4) = 0D0
!                    end do
!                end if
!            end if
!
!            end do
!            end do
!        end if
!    end if

end subroutine

! =============================================================================
! G_DHAT
! =============================================================================
subroutine G_DHAT(Jn,Dt,vphi,sphi,Dh)
    implicit none
    real(8), intent(in) :: Jn(:,:,:,:), Dt(:,:,:,:), vphi(:,:,:), sphi(:,:,:,:)
    real(8), intent(out):: Dh(:,:,:,:)

    do kk = 1, ncm(3)
    do jj = 1, ncm(2)
    do ii = 1, ncm(1)
        ! x0 +
        Dh(ii,jj,kk,1) = (Jn(ii,jj,kk,1)+Dt(ii,jj,kk,1) &
            *(vphi(ii,jj,kk)-sphi(ii,jj,kk,1))) &
            /(vphi(ii,jj,kk)+sphi(ii,jj,kk,1))
        ! x1 -
        Dh(ii,jj,kk,2) = (Jn(ii,jj,kk,2)+Dt(ii,jj,kk,2) &
            *(sphi(ii,jj,kk,2)-vphi(ii,jj,kk))) &
            /(sphi(ii,jj,kk,2)+vphi(ii,jj,kk))
        ! y0 +
        Dh(ii,jj,kk,3) = (Jn(ii,jj,kk,3)+Dt(ii,jj,kk,3) &
            *(vphi(ii,jj,kk)-sphi(ii,jj,kk,3))) &
            /(vphi(ii,jj,kk)+sphi(ii,jj,kk,3))
        ! y1 -
        Dh(ii,jj,kk,4) = (Jn(ii,jj,kk,4)+Dt(ii,jj,kk,4) &
            *(sphi(ii,jj,kk,4)-vphi(ii,jj,kk))) &
            /(sphi(ii,jj,kk,4)+vphi(ii,jj,kk))
        ! z0 +
        Dh(ii,jj,kk,5) = (Jn(ii,jj,kk,5)+Dt(ii,jj,kk,5) &
            *(vphi(ii,jj,kk)-sphi(ii,jj,kk,5))) &
            /(vphi(ii,jj,kk)+sphi(ii,jj,kk,5))
        ! z1 -
        Dh(ii,jj,kk,6) = (Jn(ii,jj,kk,6)+Dt(ii,jj,kk,6) &
            *(sphi(ii,jj,kk,6)-vphi(ii,jj,kk))) &
            /(sphi(ii,jj,kk,6)+vphi(ii,jj,kk))
    end do
    end do
    end do

    ! boundary condition
    ii = 1;      Dh(ii,:,:,1) = Jn(ii,:,:,1) / vphi(ii,:,:)
    ii = ncm(1); Dh(ii,:,:,2) = Jn(ii,:,:,2) / vphi(ii,:,:)
    jj = 1;      Dh(:,jj,:,3) = Jn(:,jj,:,3) / vphi(:,jj,:)
    jj = ncm(2); Dh(:,jj,:,4) = Jn(:,jj,:,4) / vphi(:,jj,:)
    kk = 1;      Dh(:,:,kk,5) = Jn(:,:,kk,5) / vphi(:,:,kk)
    kk = ncm(3); Dh(:,:,kk,6) = Jn(:,:,kk,6) / vphi(:,:,kk)

!    !!! zigzag !!!
!    if  ( zigzag ) then
!        ! quarter core
!        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
!            do kk = 1, cm(3)
!            ! -----
!            ii = cm(1)-gzz(2)
!            do jj = cm(2)-gzz(1)+1, cm(2)
!                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
!            end do
!            ! -----
!            ii = cm(1)-gzz(1)
!            do jj = cm(2)-gzz(2)+1, cm(2)-gzz(1)
!                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
!            end do
!            ! -----
!            jj = cm(2)-gzz(2)
!            do ii = cm(1)-gzz(1)+1, cm(1)
!                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
!            end do
!            ! -----
!            jj = cm(2)-gzz(1)
!            do ii = cm(1)-gzz(2)+1, cm(1)-gzz(1)
!                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
!            end do
!            end do
!    
!        ! whole core
!        else
!            do ii = 1, cm(1); id(1) = abs(nint(ii-gmp(1)))
!            do jj = 1, cm(2); id(2) = abs(nint(jj-gmp(2)))
!        
!            ! boundary
!            if ( id(1) == cm(1)/2-gzz(2) .and. &
!                cm(2)/2-gzz(1) < id(2) .and. id(2) <= cm(2)/2 ) then
!                if ( nint(ii-gmp(1)) <= 0 ) then
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,1) = cmJn(ii,jj,kk,:,1)/cm_phi2(ii,jj,kk,:)
!                end do
!                else
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
!                end do
!                end if
!            end if
!            ! -----
!            if ( id(1) == cm(1)/2-gzz(1) .and. &
!                cm(2)/2-gzz(2) < id(2) .and. id(2) <= cm(2)/2-gzz(1) ) then
!                if ( nint(ii-gmp(1)) <= 0 ) then
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,1) = cmJn(ii,jj,kk,:,1)/cm_phi2(ii,jj,kk,:)
!                end do
!                else
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
!                end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == cm(2)/2-gzz(2) .and. &
!                cm(1)/2-gzz(1) < id(1) .and. id(1) <= cm(1)/2 ) then
!                if ( nint(jj-gmp(2)) <= 0 ) then
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,3) = cmJn(ii,jj,kk,:,3)/cm_phi2(ii,jj,kk,:)
!                end do
!                else
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
!                end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == cm(2)/2-gzz(1) .and. &
!                cm(1)/2-gzz(2) < id(1) .and. id(1) <= cm(1)/2-gzz(1) ) then
!                if ( nint(jj-gmp(2)) <= 0 ) then
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,3) = cmJn(ii,jj,kk,:,3)/cm_phi2(ii,jj,kk,:)
!                end do
!                else
!                do kk = 1, cm(3)
!                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
!                end do
!                end if
!            end if
!    
!            end do
!            end do
!        end if
!    end if

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

!    !!! zigzag !!!
!    if  ( zigzag ) then
!        ! quarter core
!        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
!            do kk = 1, cm(3)
!            ! -----
!            ii = cm(1)-gzz(2)
!            do jj = cm(2)-gzz(1)+1, cm(2)
!                deltc0(ii,jj,kk,:,2) = 0
!                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
!            end do
!            ! -----
!            ii = cm(1)-gzz(1)
!            do jj = cm(2)-gzz(2)+1, cm(2)-gzz(1)
!                deltc0(ii,jj,kk,:,2) = 0
!                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
!            end do
!            ! -----
!            jj = cm(2)-gzz(2)
!            do ii = cm(1)-gzz(1)+1, cm(1)
!                deltc0(ii,jj,kk,:,4) = 0
!                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
!            end do
!            ! -----
!            jj = cm(2)-gzz(1)
!            do ii = cm(1)-gzz(2)+1, cm(1)-gzz(1)
!                deltc0(ii,jj,kk,:,4) = 0
!                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
!            end do
!            end do
!    
!        ! whole core
!        else
!            do ii = 1, cm(1); id(1) = abs(nint(ii-gmp(1)))
!            do jj = 1, cm(2); id(2) = abs(nint(jj-gmp(2)))
!        
!            ! boundary
!            if ( id(1) == cm(1)/2-gzz(2) .and. &
!                cm(2)/2-gzz(1) < id(2) .and. id(2) <= cm(2)/2 ) then
!                if ( nint(ii-gmp(1)) <= 0 ) then
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,1) = 0
!                deltc1(ii,jj,kk,:,1) = cm_dhat(ii,jj,kk,:,1)
!                end do
!                else
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,2) = 0
!                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
!                end do
!                end if
!            end if
!            ! -----
!            if ( id(1) == cm(1)/2-gzz(1) .and. &
!                cm(2)/2-gzz(2) < id(2) .and. id(2) <= cm(2)/2-gzz(1) ) then
!                if ( nint(ii-gmp(1)) <= 0 ) then
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,1) = 0
!                deltc1(ii,jj,kk,:,1) = cm_dhat(ii,jj,kk,:,1)
!                end do
!                else
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,2) = 0
!                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
!                end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == cm(2)/2-gzz(2) .and. &
!                cm(1)/2-gzz(1) < id(1) .and. id(1) <= cm(1)/2 ) then
!                if ( nint(jj-gmp(2)) <= 0 ) then
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,3) = 0
!                deltc1(ii,jj,kk,:,3) = cm_dhat(ii,jj,kk,:,3)
!                end do
!                else
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,4) = 0
!                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
!                end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == cm(2)/2-gzz(1) .and. &
!                cm(1)/2-gzz(2) < id(1) .and. id(1) <= cm(1)/2-gzz(1) ) then
!                if ( nint(jj-gmp(2)) <= 0 ) then
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,3) = 0
!                deltc1(ii,jj,kk,:,3) = cm_dhat(ii,jj,kk,:,3)
!                end do
!                else
!                do kk = 1, cm(3)
!                deltc0(ii,jj,kk,:,4) = 0
!                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
!                end do
!                end if
!            end if
!    
!            end do
!            end do
!        end if
!    end if


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

!    !!! zigzag !!!
!    if ( zigzag ) then
!    ! out of domain
!    where ( cm_phi2(:,:,:,:) == 0 )
!        gm1(:,:,:,:) = 0D0
!        gm2(:,:,:,:) = 0D0
!        gm3(:,:,:,:) = 0D0
!        gm4(:,:,:,:) = 1D0
!        gm5(:,:,:,:) = 0D0
!        gm6(:,:,:,:) = 0D0
!        gm7(:,:,:,:) = 0D0
!    end where
!    end if
!
!    !!! zigzag !!!
!    if ( zigzag ) then
!        ! quarter core
!        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
!            do kk = 1, cm(3)
!            ! -----
!            ii = cm(1)-gzz(2)
!            do jj = cm(2)-gzz(1)+1, cm(2)
!                gm5(ii,jj,kk,:) = 0
!            end do
!            ! -----
!            ii = cm(1)-gzz(1)
!            do jj = cm(2)-gzz(2)+1, cm(2)-gzz(1)
!                gm5(ii,jj,kk,:) = 0
!            end do
!            ! -----
!            jj = cm(2)-gzz(2)
!            do ii = cm(1)-gzz(1)+1, cm(1)
!                gm6(ii,jj,kk,:) = 0
!            end do
!            ! -----
!            jj = cm(2)-gzz(1)
!            do ii = cm(1)-gzz(2)+1, cm(1)-gzz(1)
!                gm6(ii,jj,kk,:) = 0
!            end do
!            end do
!    
!        ! whole core
!        else
!            do ii = 1, cm(1); id(1) = abs(nint(ii-gmp(1)))
!            do jj = 1, cm(2); id(2) = abs(nint(jj-gmp(2)))
!        
!            ! boundary
!            if ( id(1) == cm(1)/2-gzz(2) .and. &
!                cm(2)/2-gzz(1) < id(2) .and. id(2) <= cm(2)/2 ) then
!                if ( nint(ii-gmp(1)) <= 0 ) then
!                    do kk = 1, cm(3)
!                    gm3(ii,jj,kk,:) = 0D0
!                    end do
!                else
!                    do kk = 1, cm(3)
!                    gm5(ii,jj,kk,:) = 0D0
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(1) == cm(1)/2-gzz(1) .and. &
!                cm(2)/2-gzz(2) < id(2) .and. id(2) <= cm(2)/2-gzz(1) ) then
!                if ( nint(ii-gmp(1)) <= 0 ) then
!                    do kk = 1, cm(3)
!                    gm3(ii,jj,kk,:) = 0D0
!                    end do
!                else
!                    do kk = 1, cm(3)
!                    gm5(ii,jj,kk,:) = 0D0
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == cm(2)/2-gzz(2) .and. &
!                cm(1)/2-gzz(1) < id(1) .and. id(1) <= cm(1)/2 ) then
!                if ( nint(jj-gmp(2)) <= 0 ) then
!                    do kk = 1, cm(3)
!                    gm2(ii,jj,kk,:) = 0D0
!                    end do
!                else
!                    do kk = 1, cm(3)
!                    gm6(ii,jj,kk,:) = 0D0
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == cm(2)/2-gzz(1) .and. &
!                cm(1)/2-gzz(2) < id(1) .and. id(1) <= cm(1)/2-gzz(1) ) then
!                if ( nint(jj-gmp(2)) <= 0 ) then
!                    do kk = 1, cm(3)
!                    gm2(ii,jj,kk,:) = 0D0
!                    end do
!                else
!                    do kk = 1, cm(3)
!                    gm6(ii,jj,kk,:) = 0D0
!                    end do
!                end if
!            end if
!    
!            end do
!            end do
!        end if
!    end if

end subroutine

!! =============================================================================
!! ALBEDO calculates an albedo parameter
!! =============================================================================
!subroutine ALBEDO(aa,jj1,jj0)
!    real(8), intent(in)::  jj1(:), jj0(:)
!    real(8), intent(out):: aa(:)
!    real(8):: bb(cm_eng)
!
!    do ee = 1, cm_eng
!    if ( jj0(ee) == 0 ) then
!        if ( jj1(ee) == 0 ) then
!            aa(ee) = 0D0
!        else
!            aa(ee) = 2D0
!        end if
!    else if ( jj0(ee) == jj1(ee) ) then
!        aa(ee) = 0D0
!    else
!        bb(ee) = jj1(ee)/jj0(ee)
!        aa(ee) = 2D0*(bb(ee)+1D0)/(bb(ee)-1D0)
!    end if
!    end do
!
!end subroutine

!! =============================================================================
!! G_SOURCE
!! =============================================================================
!subroutine G_SOURCE
!    implicit none
!
!    ! fission source
!    cm_s(:,:,:) = cm_nf(:,:,:)*cm_phi1(:,:,:)/k_fmfd
!
!!    write(*,2), cm_s(:,:,:)
!!    write(*,*)
!!    write(*,2), cm_nf(:,:,:)
!!    write(*,*)
!!    write(*,2), cm_phi1(:,:,:)
!!    write(*,*)
!!    stop
!!    write(*,2), Mcm(:,:,:,2)
!!    write(*,*)
!!    write(*,2), Mcm(:,:,:,3)
!!    write(*,*)
!!    write(*,2), Mcm(:,:,:,4)
!!    write(*,*)
!!    write(*,2), Mcm(:,:,:,5)
!!    write(*,*)
!!    write(*,2), Mcm(:,:,:,6)
!!    write(*,*)
!!    write(*,2), Mcm(:,:,:,7)
!!    write(*,*)
!!    2 format(2es15.7)
!!    stop
!
!end subroutine

! =============================================================================
! G_POWER
! =============================================================================
subroutine G_POWER
    implicit none

    ! k update
    k_fmfd = k_fmfd*sum(cm_nf*cm_phi1*cm_nf*cm_phi1) &
           / sum(cm_nf*cm_phi0*cm_nf*cm_phi1)

end subroutine

! =============================================================================
! G_INJ
! =============================================================================
subroutine G_INJ(Dt,Dh,phi)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:), phi(:,:,:)
    real(8):: deno

!    do kk = 1, ncm(3)
!    do jj = 1, ncm(2)
!    do ii = 1, ncm(1)
!        if ( ii /= 1 ) then ! x-direction
!        deno = Dt(ii-1,jj,kk,2)+Dt(ii,jj,kk,1)+Dh(ii,jj,kk,1)-Dh(ii-1,jj,kk,2)
!        cmF(ii,jj,kk,1) = ((Dt(ii,jj,kk,1)-Dh(ii,jj,kk,1)) &
!            *phi(ii,jj,kk)+(Dt(ii-1,jj,kk,2) &
!            +Dh(ii-1,jj,kk,2))*phi(ii-1,jj,kk))/deno
!        cmJn(ii,jj,kk,1) = -Dt(ii,jj,kk,1)*(phi(ii,jj,kk)-cmF(ii,jj,kk,1)) &
!                           +Dh(ii,jj,kk,1)*(phi(ii,jj,kk)+cmF(ii,jj,kk,1))
!        cmJ0(ii,jj,kk,1) = 25D-2*cmF(ii,jj,kk,1)-5D-1*cmJn(ii,jj,kk,1)
!        cmJ1(ii,jj,kk,1) = cmJ0(ii,jj,kk,1)+cmJn(ii,jj,kk,1)
!
!        cmF(ii-1,jj,kk,2)  = cmF(ii,jj,kk,1)
!        cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)
!        cmJ0(ii-1,jj,kk,2) = cmJ0(ii,jj,kk,1)
!        cmJ1(ii-1,jj,kk,2) = cmJ1(ii,jj,kk,1)
!        end if
!        if ( jj /= 1 ) then ! y-direction
!        deno = Dt(ii,jj-1,kk,4)+Dt(ii,jj,kk,3)+Dh(ii,jj,kk,3)-Dh(ii,jj-1,kk,4)
!        cmF(ii,jj,kk,3) = ((Dt(ii,jj,kk,3)-Dh(ii,jj,kk,3)) &
!            *phi(ii,jj,kk)+(Dt(ii,jj-1,kk,4) &
!            +Dh(ii,jj-1,kk,4))*phi(ii,jj-1,kk))/deno
!        cmJn(ii,jj,kk,3) = -Dt(ii,jj,kk,3)*(phi(ii,jj,kk)-cmF(ii,jj,kk,3)) &
!                           +Dh(ii,jj,kk,3)*(phi(ii,jj,kk)+cmF(ii,jj,kk,3))
!        cmJ0(ii,jj,kk,3) = 25D-2*cmF(ii,jj,kk,3)-5D-1*cmJn(ii,jj,kk,3)
!        cmJ1(ii,jj,kk,3) = cmJ0(ii,jj,kk,3)+cmJn(ii,jj,kk,3)
!
!        cmF(ii,jj-1,kk,4)  = cmF(ii,jj,kk,3)
!        cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)
!        cmJ0(ii,jj-1,kk,4) = cmJ0(ii,jj,kk,3)
!        cmJ1(ii,jj-1,kk,4) = cmJ1(ii,jj,kk,3)
!        end if
!        if ( kk /= 1 ) then ! z-direction
!        deno = Dt(ii,jj,kk-1,6)+Dt(ii,jj,kk,5)+Dh(ii,jj,kk,5)-Dh(ii,jj,kk-1,6)
!        cmF(ii,jj,kk,5) = ((Dt(ii,jj,kk,5)-Dh(ii,jj,kk,5)) &
!            *phi(ii,jj,kk)+(Dt(ii,jj,kk-1,6) &
!            +Dh(ii,jj,kk-1,6))*phi(ii,jj,kk-1))/deno
!        cmJn(ii,jj,kk,5) = -Dt(ii,jj,kk,5)*(phi(ii,jj,kk)-cmF(ii,jj,kk,5)) &
!                           +Dh(ii,jj,kk,5)*(phi(ii,jj,kk)+cmF(ii,jj,kk,5))
!        cmJ0(ii,jj,kk,5) = 25D-2*cmF(ii,jj,kk,5)-5D-1*cmJn(ii,jj,kk,5)
!        cmJ1(ii,jj,kk,5) = cmJ0(ii,jj,kk,5)+cmJn(ii,jj,kk,5)
!
!        cmF(ii,jj,kk-1,6)  = cmF(ii,jj,kk,5)
!        cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)
!        cmJ0(ii,jj,kk-1,6) = cmJ0(ii,jj,kk,5)
!        cmJ1(ii,jj,kk-1,6) = cmJ1(ii,jj,kk,5)
!        end if
!    end do
!    end do
!    end do


    ! incoming partial current for FMFD boundary condition
    do kk = 1, ncm(3)
    do jj = 1, ncm(2)
    do ii = 1, ncm(1)
        if ( ii /= 1 ) then ! x-direction
        deno = Dt(ii-1,jj,kk,2)+Dt(ii,jj,kk,1)+Dh(ii,jj,kk,1)-Dh(ii-1,jj,kk,2)
        cmF(ii,jj,kk,1) = ((Dt(ii,jj,kk,1)-Dh(ii,jj,kk,1)) &
            *phi(ii,jj,kk)+(Dt(ii-1,jj,kk,2) &
            +Dh(ii-1,jj,kk,2))*phi(ii-1,jj,kk))/deno
        cmF(ii-1,jj,kk,2) = cmF(ii,jj,kk,1)

        cmJn(ii,jj,kk,1) = -Dt(ii,jj,kk,1)*(phi(ii,jj,kk)-cmF(ii,jj,kk,1)) &
                           +Dh(ii,jj,kk,1)*(phi(ii,jj,kk)+cmF(ii,jj,kk,1))
        cmJn(ii-1,jj,kk,2) = -Dt(ii-1,jj,kk,2)*(cmF(ii-1,jj,kk,2)-phi(ii-1,jj,kk)) &
                             +Dh(ii-1,jj,kk,2)*(cmF(ii-1,jj,kk,2)+phi(ii-1,jj,kk))
!        cmJn(ii,jj,kk,1) = (cmJn(ii,jj,kk,1)+cmJn(ii-1,jj,kk,2))/2D0
!        cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)

        cmJ0(ii-1,jj,kk,2) = 25D-2*cmF(ii-1,jj,kk,2)-5D-1*cmJn(ii-1,jj,kk,2)
        cmJ1(ii,jj,kk,1)   = 25D-2*cmF(ii  ,jj,kk,1)+5D-1*cmJn(ii,jj,kk,1)
!        cmJ0(ii,jj,kk,1)   = cmJ0(ii-1,jj,kk,2)
!        cmJ1(ii-1,jj,kk,2) = cmJ1(ii,jj,kk,1)
        end if
        if ( jj /= 1 ) then ! y-direction
        deno = Dt(ii,jj-1,kk,4)+Dt(ii,jj,kk,3)+Dh(ii,jj,kk,3)-Dh(ii,jj-1,kk,4)
        cmF(ii,jj,kk,3) = ((Dt(ii,jj,kk,3)-Dh(ii,jj,kk,3)) &
            *phi(ii,jj,kk)+(Dt(ii,jj-1,kk,4) &
            +Dh(ii,jj-1,kk,4))*phi(ii,jj-1,kk))/deno
        cmF(ii,jj-1,kk,4)  = cmF(ii,jj,kk,3)

        cmJn(ii,jj,kk,3) = -Dt(ii,jj,kk,3)*(phi(ii,jj,kk)-cmF(ii,jj,kk,3)) &
                           +Dh(ii,jj,kk,3)*(phi(ii,jj,kk)+cmF(ii,jj,kk,3))
        cmJn(ii,jj-1,kk,4) = -Dt(ii,jj-1,kk,4)*(cmF(ii,jj-1,kk,4)-phi(ii,jj-1,kk)) &
                             +Dh(ii,jj-1,kk,4)*(cmF(ii,jj-1,kk,4)-phi(ii,jj-1,kk))
!        cmJn(ii,jj,kk,3) = (cmJn(ii,jj,kk,3)+cmJn(ii,jj-1,kk,4))/2D0
!        cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)

        cmJ0(ii,jj-1,kk,4) = 25D-2*cmF(ii,jj-1,kk,4)-5D-1*cmJn(ii,jj-1,kk,4)
        cmJ1(ii,jj,kk,3)   = 25D-2*cmF(ii,jj,kk,3)+5D-1*cmJn(ii,jj,kk,3)
!        cmJ0(ii,jj,kk,3)   = cmJ0(ii,jj-1,kk,4)
!        cmJ1(ii,jj-1,kk,4) = cmJ1(ii,jj,kk,3)
        end if
        if ( kk /= 1 ) then ! z-direction
        deno = Dt(ii,jj,kk-1,6)+Dt(ii,jj,kk,5)+Dh(ii,jj,kk,5)-Dh(ii,jj,kk-1,6)
        cmF(ii,jj,kk,5) = ((Dt(ii,jj,kk,5)-Dh(ii,jj,kk,5)) &
            *phi(ii,jj,kk)+(Dt(ii,jj,kk-1,6) &
            +Dh(ii,jj,kk-1,6))*phi(ii,jj,kk-1))/deno
        cmF(ii,jj,kk-1,6)  = cmF(ii,jj,kk,5)

        cmJn(ii,jj,kk,5) = -Dt(ii,jj,kk,5)*(phi(ii,jj,kk)-cmF(ii,jj,kk,5)) &
                           +Dh(ii,jj,kk,5)*(phi(ii,jj,kk)+cmF(ii,jj,kk,5))
        cmJn(ii,jj,kk-1,6) = -Dt(ii,jj,kk-1,6)*(cmF(ii,jj,kk-1,6)-phi(ii,jj,kk-1)) &
                             +Dh(ii,jj,kk-1,6)*(cmF(ii,jj,kk-1,6)+phi(ii,jj,kk-1))
!        cmJn(ii,jj,kk,5) = (cmJn(ii,jj,kk,5)+cmJn(ii,jj,kk-1,6))/2D0
!        cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)

        cmJ0(ii,jj,kk-1,6) = 25D-2*cmF(ii,jj,kk-1,6)-5D-1*cmJn(ii,jj,kk-1,6)
        cmJ1(ii,jj,kk,5)   = 25D-2*cmF(ii,jj,kk,5)+5D-1*cmJn(ii,jj,kk,5)
!        cmJ0(ii,jj,kk,5)   = cmJ0(ii,jj,kk-1,6)
!        cmJ1(ii,jj,kk-1,6) = cmJ1(ii,jj,kk,5)
        end if
    end do
    end do
    end do



!    write(8,1), cmJ0(:,:,:,1)
!    write(8,*)
!    write(8,1), cmJ0(:,:,:,2)
!    write(8,*)
!    write(8,1), cmJ0(:,:,:,3)
!    write(8,*)
!    write(8,1), cmJ0(:,:,:,4)
!    write(8,*)
!    write(8,1), cmJ0(:,:,:,5)
!    write(8,*)
!    write(8,1), cmJ0(:,:,:,6)
!    write(8,*)
!    1 format(2es15.7)
!    stop





end subroutine

! =============================================================================
! G2L carries out the flux and current modulation (from GLOBAL to LOCAL)
! =============================================================================
subroutine G2L(phi0,phi1,fmJ0,fmJ1)
    implicit none
    real(8), intent(inout):: phi1(:,:,:), fmJ0(:,:,:,:), fmJ1(:,:,:,:)
    real(8), intent(out)::   phi0(:,:,:)
    real(8):: ssum(1:2)

    phi0(:,:,:) = phi1(:,:,:)

    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr

    ! flux modulation
    phi1(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz) = &
    phi1(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz) &
    /sum(phi1(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz)) &
    *fcr*fcr*fcz*cm_phi1(ii,jj,kk)

    ! partial current modulation
    ! x0
    ssum = 0;       id(1) = id0(1)+1
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
!        ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),1)
        ssum(2) = ssum(2) + fmJ1(id(1),id(2),id(3),1)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
!        if ( ssum(1) /= 0 ) fmJ0(id(1),id(2),id(3),1) = &
!            fmJ0(id(1),id(2),id(3),1)/ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,1)
        if ( ssum(2) /= 0 ) fmJ1(id(1),id(2),id(3),1) = &
            fmJ1(id(1),id(2),id(3),1)/ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,1)
    end do
    end do
    ! x1
    ssum = 0;       id(1) = id0(1)+fcr
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),2)
!        ssum(2) = ssum(2) + fmJ1(id(1),id(2),id(3),2)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum(1) /= 0 ) fmJ0(id(1),id(2),id(3),2) = &
            fmJ0(id(1),id(2),id(3),2)/ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,2)
!        if ( ssum(2) /= 0 ) fmJ1(id(1),id(2),id(3),2) = &
!            fmJ1(id(1),id(2),id(3),2)/ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,2)
    end do
    end do
    ! y0
    ssum = 0;       id(2) = id0(2)+1
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
!        ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),3)
        ssum(2) = ssum(2) + fmJ1(id(1),id(2),id(3),3)
    end do 
    end do 
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
!        if ( ssum(1) /= 0 ) fmJ0(id(1),id(2),id(3),3) = &
!            fmJ0(id(1),id(2),id(3),3)/ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,3)
        if ( ssum(2) /= 0 ) fmJ1(id(1),id(2),id(3),3) = &
            fmJ1(id(1),id(2),id(3),3)/ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,3)
    end do
    end do 
    ! y1
    ssum = 0;       id(2) = id0(2)+fcr
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),4)
!        ssum(2) = ssum(2) + fmJ1(id(1),id(2),id(3),4)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        if ( ssum(1) /= 0 ) fmJ0(id(1),id(2),id(3),4) = &
            fmJ0(id(1),id(2),id(3),4)/ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,4)
!        if ( ssum(2) /= 0 ) fmJ1(id(1),id(2),id(3),4) = &
!            fmJ1(id(1),id(2),id(3),4)/ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,4)
    end do
    end do
    ! z0
    ssum = 0;       id(3) = id0(3)+1
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
!        ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),5)
        ssum(2) = ssum(2) + fmJ1(id(1),id(2),id(3),5)
    end do
    end do
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
!        if ( ssum(1) /= 0 ) fmJ0(id(1),id(2),id(3),5) = &
!            fmJ0(id(1),id(2),id(3),5)/ssum(1)*fcr*fcr*cmJ0(ii,jj,kk,5)
        if ( ssum(2) /= 0 ) fmJ1(id(1),id(2),id(3),5) = &
            fmJ1(id(1),id(2),id(3),5)/ssum(2)*fcr*fcr*cmJ1(ii,jj,kk,5)
    end do
    end do
    ! z1
    ssum = 0;       id(3) = id0(3)+fcz
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJ0(id(1),id(2),id(3),6)
!        ssum(2) = ssum(2) + fmJ1(id(1),id(2),id(3),6)
    end do
    end do
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
!        if ( ssum(1) /= 0 ) fmJ0(id(1),id(2),id(3),6) = &
!            fmJ0(id(1),id(2),id(3),6)/ssum(1)*fcr*fcr*cmJ0(ii,jj,kk,6)
        if ( ssum(2) /= 0 ) fmJ1(id(1),id(2),id(3),6) = &
            fmJ1(id(1),id(2),id(3),6)/ssum(2)*fcr*fcr*cmJ1(ii,jj,kk,6)
    end do
    end do

    end do
    end do
    end do

!    write(8,1), fmJ1(:,:,:,1)
!    write(8,*)
!    write(8,1), fmJ1(:,:,:,2)
!    write(8,*)
!    write(8,1), fmJ1(:,:,:,3)
!    write(8,*)
!    write(8,1), fmJ1(:,:,:,4)
!    write(8,*)
!    write(8,1), fmJ1(:,:,:,5)
!    write(8,*)
!    write(8,1), fmJ1(:,:,:,6)
!    write(8,*)
!    write(8,2), cmJ1(:,:,:,1)
!    write(8,*)
!    write(8,2), cmJ1(:,:,:,2)
!    write(8,*)
!    write(8,2), cmJ1(:,:,:,3)
!    write(8,*)
!    write(8,2), cmJ1(:,:,:,4)
!    write(8,*)
!    write(8,2), cmJ1(:,:,:,5)
!    write(8,*)
!    write(8,2), cmJ1(:,:,:,6)
!    write(8,*)
!    1 format(20es15.6)
!    2 format(2es15.6)
!    stop


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
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        ! x-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            if ( ii /= 1 ) then
            id(1) = id0(1)+1
            fm_s(id(1),id(2),id(3)) = fm_s(id(1),id(2),id(3)) &
                +(jsrc(id(1),id(2),id(3),1)*fmJ1(id(1),id(2),id(3),1) &
                +fsrc(id(1),id(2),id(3),1)*phi0(id(1)-1,id(2),id(3)))/dfm(1)
            end if
            if ( ii /= ncm(1) ) then
            id(1) = id0(1)+fcr
            fm_s(id(1),id(2),id(3)) = fm_s(id(1),id(2),id(3)) &
                +(jsrc(id(1),id(2),id(3),2)*fmJ0(id(1),id(2),id(3),2) &
                -fsrc(id(1),id(2),id(3),2)*phi0(id(1)+1,id(2),id(3)))/dfm(1)
            end if
        end do
        ! y-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
            if ( jj /= 1 ) then
            id(2) = id0(2)+1
            fm_s(id(1),id(2),id(3)) = fm_s(id(1),id(2),id(3)) &
                +(jsrc(id(1),id(2),id(3),3)*fmJ1(id(1),id(2),id(3),3) &
                +fsrc(id(1),id(2),id(3),3)*phi0(id(1),id(2)-1,id(3)))/dfm(2)
            end if
            if ( jj /= ncm(2) ) then
            id(2) = id0(2)+fcr
            fm_s(id(1),id(2),id(3)) = fm_s(id(1),id(2),id(3)) &
                +(jsrc(id(1),id(2),id(3),4)*fmJ0(id(1),id(2),id(3),4) &
                -fsrc(id(1),id(2),id(3),4)*phi0(id(1),id(2)+1,id(3)))/dfm(2)
            end if
        end do
        end do
        ! z-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            if ( kk /= 1 ) then
            id(3) = id0(3)+1
            fm_s(id(1),id(2),id(3)) = fm_s(id(1),id(2),id(3)) &
                +(jsrc(id(1),id(2),id(3),5)*fmJ1(id(1),id(2),id(3),5) &
                +fsrc(id(1),id(2),id(3),5)*phi0(id(1),id(2),id(3)-1))/dfm(3)
            end if
            if ( kk /= ncm(3) ) then
            id(3) = id0(3)+fcz
            fm_s(id(1),id(2),id(3)) = fm_s(id(1),id(2),id(3)) &
                +(jsrc(id(1),id(2),id(3),6)*fmJ0(id(1),id(2),id(3),6) &
                -fsrc(id(1),id(2),id(3),6)*phi0(id(1),id(2),id(3)+1))/dfm(3)
            end if
        end do
        end do
    end do
    end do
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
    real(8):: netJ

    ! outgoing partial current
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        ! x0
        id(1) = id0(1)+1
        if ( ii == 1 ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),1) = deltf1(id(1),id(2),id(3),1) &
                *phi1(id(1),id(2),id(3))
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = +jsrc(id(1),id(2),id(3),1)*fmJ1(id(1),id(2),id(3),1) &
                +deltf1(id(1),id(2),id(3),1)*phi1(id(1),id(2),id(3)) &
                +fsrc(id(1),id(2),id(3),1)*phi0(id(1)-1,id(2),id(3))
            fmF(id(1),id(2),id(3),1) = &
                4D0*fmJ1(id(1),id(2),id(3),1)-2D0*netJ
            fmJ0(id(1),id(2),id(3),1) = &
                25D-2*fmF(id(1),id(2),id(3),1)-5D-1*netJ
        end do
        end do
        end if
        ! x1
        id(1) = id0(1)+fcr
        if ( ii == ncm(1) ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),2) = deltf1(id(1),id(2),id(3),2) &
                *phi1(id(1),id(2),id(3))
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = -jsrc(id(1),id(2),id(3),2)*fmJ0(id(1),id(2),id(3),2) &
                +deltf1(id(1),id(2),id(3),2)*phi1(id(1),id(2),id(3)) &
                +fsrc(id(1),id(2),id(3),2)*phi0(id(1)+1,id(2),id(3))
            fmF(id(1),id(2),id(3),2) = &
                4D0*fmJ0(id(1),id(2),id(3),2)+2D0*netJ
            fmJ1(id(1),id(2),id(3),2) = &
                25D-2*fmF(id(1),id(2),id(3),2)+5D-1*netJ
        end do
        end do
        end if
        ! y0
        id(2) = id0(2)+1
        if ( jj == 1 ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJn(id(1),id(2),id(3),3) = deltf1(id(1),id(2),id(3),3) &
                *phi1(id(1),id(2),id(3))
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            netJ = +jsrc(id(1),id(2),id(3),3)*fmJ1(id(1),id(2),id(3),3) &
                +deltf1(id(1),id(2),id(3),3)*phi1(id(1),id(2),id(3)) &
                +fsrc(id(1),id(2),id(3),3)*phi0(id(1),id(2)-1,id(3))
            fmF(id(1),id(2),id(3),3) = &
                4D0*fmJ1(id(1),id(2),id(3),3)-2D0*netJ
            fmJ0(id(1),id(2),id(3),3) = &
                25D-2*fmF(id(1),id(2),id(3),3)-5D-1*netJ
        end do
        end do
        end if
        ! y1
        id(2) = id0(2)+fcr
        if ( jj == ncm(2) ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJn(id(1),id(2),id(3),4) = deltf1(id(1),id(2),id(3),4) &
                *phi1(id(1),id(2),id(3))
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            netJ = -jsrc(id(1),id(2),id(3),4)*fmJ0(id(1),id(2),id(3),4) &
                +deltf1(id(1),id(2),id(3),4)*phi1(id(1),id(2),id(3)) &
                +fsrc(id(1),id(2),id(3),4)*phi0(id(1),id(2)+1,id(3))
            fmF(id(1),id(2),id(3),4) = &
                4D0*fmJ0(id(1),id(2),id(3),4)+2D0*netJ
            fmJ1(id(1),id(2),id(3),4) = &
                25D-2*fmF(id(1),id(2),id(3),4)+5D-1*netJ
        end do
        end do
        end if
        ! z0
        id(3) = id0(3)+1
        if ( kk == 1 ) then
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),5) = deltf1(id(1),id(2),id(3),5) &
                *phi1(id(1),id(2),id(3))
        end do
        end do
        else
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = +jsrc(id(1),id(2),id(3),5)*fmJ1(id(1),id(2),id(3),5) &
                +deltf1(id(1),id(2),id(3),5)*phi1(id(1),id(2),id(3)) &
                +fsrc(id(1),id(2),id(3),5)*phi0(id(1),id(2),id(3)-1)
            fmF(id(1),id(2),id(3),5) = &
                4D0*fmJ1(id(1),id(2),id(3),5)-2D0*netJ
            fmJ0(id(1),id(2),id(3),5) = &
                25D-2*fmF(id(1),id(2),id(3),5)-5D-1*netJ
        end do
        end do
        end if
        ! z1
        id(3) = id0(3)+fcz
        if ( kk == ncm(3) ) then
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),6) = deltf1(id(1),id(2),id(3),6) &
                *phi1(id(1),id(2),id(3))
        end do
        end do
        else
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = -jsrc(id(1),id(2),id(3),6)*fmJ0(id(1),id(2),id(3),6) &
                +deltf1(id(1),id(2),id(3),6)*phi1(id(1),id(2),id(3)) &
                +fsrc(id(1),id(2),id(3),6)*phi0(id(1),id(2),id(3)+1)
            fmF(id(1),id(2),id(3),6) = &
                4D0*fmJ0(id(1),id(2),id(3),6)+2D0*netJ
            fmJ1(id(1),id(2),id(3),6) = &
                25D-2*fmF(id(1),id(2),id(3),6)+5D-1*netJ
        end do
        end do
        end if
    end do
    end do
    end do

!    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
!    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
!    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
!        ! x0
!        if ( ii /= 1 ) then
!        id(1) = id0(1)+1
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            fmJ0(id(1)-1,id(2),id(3),2) = fmJ0(id(1),id(2),id(3),1)
!        end do
!        end do
!        end if
!        ! x1
!        if ( ii /= ncm(1) ) then
!        id(1) = id0(1)+fcr
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            fmJ1(id(1)+1,id(2),id(3),1) = fmJ1(id(1),id(2),id(3),2)
!        end do
!        end do
!        end if
!        ! y0
!        if ( jj /= 1 ) then
!        id(2) = id0(2)+1
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do mm = 1, fcr; id(1) = id0(1)+mm
!            fmJ0(id(1),id(2)-1,id(3),4) = fmJ0(id(1),id(2),id(3),3)
!        end do
!        end do
!        end if
!        ! y1
!        if ( jj /= ncm(2) ) then
!        id(2) = id0(2)+fcr
!        do oo = 1, fcz; id(3) = id0(3)+oo
!        do mm = 1, fcr; id(1) = id0(1)+mm
!            fmJ1(id(1),id(2)+1,id(3),3) = fmJ1(id(1),id(2),id(3),4)
!        end do
!        end do
!        end if
!        ! z0
!        if ( kk /= 1 ) then
!        id(3) = id0(3)+1
!        do mm = 1, fcr; id(1) = id0(1)+mm
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            fmJ0(id(1),id(2),id(3)-1,6) = fmJ0(id(1),id(2),id(3),5)
!        end do
!        end do
!        end if
!        ! z1
!        if ( kk /= ncm(3) ) then
!        id(3) = id0(3)+fcz
!        do mm = 1, fcr; id(1) = id0(1)+mm
!        do nn = 1, fcr; id(2) = id0(2)+nn
!            fmJ1(id(1),id(2),id(3)+1,5) = fmJ1(id(1),id(2),id(3),6)
!        end do
!        end do
!        end if
!    end do
!    end do
!    end do


!    !!! zigzag !!!
!    if ( zigzag ) then
!        ! quarter core
!        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
!            do kk = fm0(3), fm1(3)
!            ! -----
!            ii = fm1(1)-zz(2)
!            do jj = fm1(2)-zz(1)+1, fm1(2)
!                fmJn(ii,jj,kk,:,2) = deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
!            end do
!            ! -----
!            ii = fm1(1)-zz(1)
!            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
!                fmJn(ii,jj,kk,:,2) = deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
!            end do
!            ! -----
!            jj = fm1(2)-zz(2)
!            do ii = fm1(1)-zz(1)+1, fm1(1)
!                fmJn(ii,jj,kk,:,4) = deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
!            end do
!            ! -----
!            jj = fm1(2)-zz(1)
!            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
!                fmJn(ii,jj,kk,:,4) = deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
!            end do
!            end do
!
!        ! whole core
!        else
!            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
!            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
!
!            if ( id(1) == afm(1)/2-zz(2) .and. &
!                 afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
!                if ( nint(ii-mp(1)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,1) = &
!                            deltf1(ii,jj,kk,:,1)*fm_phi2(ii,jj,kk,:)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,2) = &
!                            deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(1) == afm(1)/2-zz(1) .and. &
!                 afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
!                if ( nint(ii-mp(1)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,1) = &
!                            deltf1(ii,jj,kk,:,1)*fm_phi2(ii,jj,kk,:)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,2) = &
!                            deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == afm(2)/2-zz(2) .and. &
!                 afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
!                if ( nint(jj-mp(2)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,3) = &
!                            deltf1(ii,jj,kk,:,3)*fm_phi2(ii,jj,kk,:)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,4) = &
!                            deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
!                    end do
!                end if
!            end if
!            ! -----
!            if ( id(2) == afm(2)/2-zz(1) .and. &
!                 afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
!                if ( nint(jj-mp(2)) <= 0 ) then
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,3) = &
!                            deltf1(ii,jj,kk,:,3)*fm_phi2(ii,jj,kk,:)
!                    end do
!                else
!                    do kk = fm0(3), fm1(3)
!                        fmJn(ii,jj,kk,:,4) = &
!                            deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
!                    end do
!                end if
!            end if
!
!            end do
!            end do
!        end if
!    end if

end subroutine


! =============================================================================
! L_REFJ
! =============================================================================
subroutine L_REFJ(fmF,fmJ0,fmJ1,fmJn)
    implicit none
    real(8), intent(in), dimension(:,:,:,:):: fmF, fmJ0, fmJ1, fmJn
    real(8):: ssum(0:2)

    ! surface average
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        ! x0
        if ( ii /= 1 ) then
        ssum = 0;       id(1) = id0(1)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),1)
            !ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),1)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),1)
        end do
        end do
        cmJ0(ii,jj,kk,1) = ssum(0) / (fcr*fcz)
        !cmJ1(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,1)  = ssum(2) / (fcr*fcz)
        end if
        ! x1
        if ( ii /= ncm(1) ) then
        ssum = 0;       id(1) = id0(1)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            !ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),2)
            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),2)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),2)
        end do
        end do
        !cmJ0(ii,jj,kk,2) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,2) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,2)  = ssum(2) / (fcr*fcz)
        end if
        ! y0
        if ( jj /= 1 ) then
        ssum = 0;       id(2) = id0(2)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),3)
            !ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),3)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),3)
        end do
        end do
        cmJ0(ii,jj,kk,3) = ssum(0) / (fcr*fcz)
        !cmJ1(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,3)  = ssum(2) / (fcr*fcz)
        end if
        ! y1
        if ( jj /= ncm(2) ) then
        ssum = 0;       id(2) = id0(2)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            !ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),4)
            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),4)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),4)
        end do
        end do
        !cmJ0(ii,jj,kk,4) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,4) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,4)  = ssum(2) / (fcr*fcz)
        end if
        ! z0
        if ( kk /= 1 ) then
        ssum = 0;       id(3) = id0(3)+1
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),5)
            !ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),5)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),5)
        end do
        end do
        cmJ0(ii,jj,kk,5) = ssum(0) / (fcr*fcr)
        !cmJ1(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,5)  = ssum(2) / (fcr*fcr)
        end if
        ! z1
        if ( kk /= ncm(3) ) then
        ssum = 0;       id(3) = id0(3)+fcz
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            !ssum(0) = ssum(0) + fmJ0(id(1),id(2),id(3),6)
            ssum(1) = ssum(1) + fmJ1(id(1),id(2),id(3),6)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),6)
        end do
        end do
        !cmJ0(ii,jj,kk,6) = ssum(0) / (fcr*fcr)
        cmJ1(ii,jj,kk,6) = ssum(1) / (fcr*fcr)
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
            cmJn(ii,jj,kk,1) = cmJ1(ii-1,jj,kk,2)-cmJ0(ii,jj,kk,1)
            cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)
            cmF(ii,jj,kk,1) = (cmF(ii,jj,kk,1)+cmF(ii-1,jj,kk,2))/2D0
            cmF(ii-1,jj,kk,2) = cmF(ii,jj,kk,1)
        end if
        ! y-direction
        if ( jj /= 1 ) then
            cmJn(ii,jj,kk,3) = cmJ1(ii,jj-1,kk,4)-cmJ0(ii,jj,kk,3)
            cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)
            cmF(ii,jj,kk,3) = (cmF(ii,jj,kk,3)+cmF(ii,jj-1,kk,4))/2D0
            cmF(ii,jj-1,kk,4) = cmF(ii,jj,kk,3)
        end if
        ! z-direction
        if ( kk /= 1 ) then
            cmJn(ii,jj,kk,5) = cmJ1(ii,jj,kk-1,6)-cmJ0(ii,jj,kk,5)
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
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),1)
    end do
    end do
    cmJn(ii,jj,kk,1) = ssum(0) / (fcr*fcz)
    end do
    !   x1
    ii = ncm(1); id(1) = nfm(1)
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),2)
    end do
    end do
    cmJn(ii,jj,kk,2) = ssum(0) / (fcr*fcz)
    end do
    !   y0
    jj = 1; id(2) = 1
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),3)
    end do
    end do
    cmJn(ii,jj,kk,3) = ssum(0) / (fcr*fcz)
    end do
    !   y1
    jj = ncm(2); id(2) = nfm(2)
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),4)
    end do
    end do
    cmJn(ii,jj,kk,4) = ssum(0) / (fcr*fcz)
    end do
    end do
    !   z0
    kk = 1; id(3) = 1
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(0) = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),5)
    end do
    end do
    cmJn(ii,jj,kk,5) = ssum(0) / (fcr*fcr)
    end do
    end do
    !   z1
    kk = ncm(3); id(3) = nfm(3)
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(0) = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),6)
    end do
    end do
    cmJn(ii,jj,kk,6) = ssum(0) / (fcr*fcr)
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
!    write(8,1), fmJn(:,:,:,5)
!    write(8,*)
!    write(8,1), fmJn(:,:,:,6)
!    write(8,*)
!    write(8,2), cmJn(:,:,:,1)
!    write(8,*)
!    write(8,2), cmJn(:,:,:,2)
!    write(8,*)
!    write(8,2), cmJn(:,:,:,3)
!    write(8,*)
!    write(8,2), cmJn(:,:,:,4)
!    write(8,*)
!    write(8,2), cmJn(:,:,:,5)
!    write(8,*)
!    write(8,2), cmJn(:,:,:,6)
!    write(8,*)
!    1 format(20es15.7)
!    2 format(2es15.7)
!    stop


!    !!! zigzag !!!
!    if ( zigzag ) then
!        ! quarter core
!        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
!            do ee = 1, cm_eng
!            do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
!            ! -----
!            ii = cm(1)-gzz(2); id(1)  = fm1(1)-zz(2)
!            jj = cm(2);        id0(2) = fm1(2)-zz(1)
!            ssum(0) = 0
!            do oo = 1, fcz; id(3) = id0(3)+oo
!            do nn = 1, fcr; id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
!            ! -----
!            ii = cm(1)-gzz(1); id(1)  = fm1(1)-zz(1)
!            jj = cm(2)-gzz(1); id0(2) = fm1(2)-zz(2)
!            ssum(0) = 0
!            do oo = 1, fcz; id(3) = id0(3)+oo
!            do nn = 1, fcr; id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
!            ! -----
!            ii = cm(1);        id0(1) = fm1(1)-zz(1)
!            jj = cm(2)-gzz(2); id(2)  = fm1(2)-zz(2)
!            ssum(0) = 0
!            do oo = 1, fcz; id(3) = id0(3)+oo
!            do mm = 1, fcr; id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
!            ! -----
!            ii = cm(1)-gzz(1); id0(1) = fm1(1)-zz(2)
!            jj = cm(2)-gzz(1); id(2)  = fm1(2)-zz(1)
!            ssum(0) = 0
!            do oo = 1, fcz; id(3) = id0(3)+oo
!            do mm = 1, fcr; id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
!            end do
!            end do
!
!        ! whole core
!        else
!            do ee = 1, cm_eng
!            ! -----
!            ii = 1+gzz(1);     id(1)  = fm0(1)+zz(1)
!            jj = 1+gzz(1);     id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
!            end do
!            jj = cm(2)-gzz(1); id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
!            end do
!            ! -----
!            ii = 1+gzz(2);     id(1)  = fm0(1)+zz(2)
!            jj = 1;            id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
!            end do
!            jj = cm(2);        id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
!            end do
!            ! -----
!            ii = cm(1)-gzz(1); id(1)  = fm1(1)-zz(1)
!            jj = 1+gzz(1);     id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
!            end do
!            jj = cm(2)-gzz(1); id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
!            end do
!            ! -----
!            ii = cm(1)-gzz(2); id(1)  = fm1(1)-zz(2)
!            jj = 1;            id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
!            end do
!            jj = cm(2);        id0(2) = (jj-1)*fcr+fm0(2)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do nn = 1, fcr;    id(2) = id0(2)+nn
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
!            end do
!
!            ! -----
!            jj = 1+gzz(1);     id(2)  = fm0(2)+zz(1)
!            ii = 1+gzz(1);     id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
!            end do
!            ii = cm(1)-gzz(1); id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
!            end do
!            ! -----
!            jj = 1+gzz(2);     id(2)  = fm0(2)+zz(2)
!            ii = 1;            id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
!            end do
!            ii = cm(2);        id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
!            end do
!            ! -----
!            jj = cm(2)-gzz(1); id(2)  = fm1(2)-zz(1)
!            ii = 1+gzz(1);     id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
!            end do
!            ii = cm(1)-gzz(1); id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
!            end do
!            ! -----
!            jj = cm(2)-gzz(2); id(2)  = fm1(2)-zz(2)
!            ii = 1;            id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
!            end do
!            ii = cm(1);        id0(1) = (ii-1)*fcr+fm0(1)-1
!            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
!            do oo = 1, fcz;    id(3) = id0(3)+oo
!            do mm = 1, fcr;    id(1) = id0(1)+mm
!                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
!            end do
!            end do
!            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
!            end do
!            end do
!        end if
!    end if

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

    ! diffusion coefficient
    !   interface diffusion coefficient
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

