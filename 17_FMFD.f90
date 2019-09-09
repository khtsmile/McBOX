module FMFD
    use GEOMETRY_PARA,  only: bc_x0, bc_y0
    use CMFD_PARA
    use CMFD_TALLY,     only: L_ACCUM
    implicit none
    integer:: ii, jj, kk, ee
    integer:: id(3)

    contains
    
! =============================================================================
! FMFD_COMPT groups subroutines for FMFD computation
! =============================================================================
subroutine FMFD_COMPT(i)
    use TALLY_PARA, only: k_com
    use GEOMETRY_PARA
    use SOLVER, only: CG, CG1, CMFD_CG
    use CMFD
    use NEU_PARA, only: fmfdup, mprup
    implicit none
    integer, intent(in):: i
    real(8):: error
    integer:: skip = 1000
    integer:: local
    integer:: global
    real(8):: k_pre

!    cmt1(i,:,:) = cm_phi2(:,:,1,1)
!    cmt2(i,:,:) = cm_tot(:,:,1,1)
!    cmt3(i,:,:) = cm_remv(:,:,1,1)
!    cmt4(i,:,:) = cm_nufiss(:,:,1,1)
!    cmt5(i,:,:,:) = cmJn(:,:,1,1,:)
!    cmt6(i,:,:,:) = cmJ1(:,:,1,1,:)
!    cmt7(i,:,:,:) = cmJ0(:,:,1,1,:)

    !!! Conventional CMFD !!!
    if ( ccmfd ) then
    call CMFD_DHAT; call MMATRIX3
    error=1D0
    do while( error > 1D-8 )
    cm_phi1 = cm_phi2
    call CMFD_SOURCE
    call CMFD_CG(mm1(1:cm(1),1:cm(2),1:cm(3),:), &
                 mm2(1:cm(1),1:cm(2),1:cm(3),:), &
                 mm3(1:cm(1),1:cm(2),1:cm(3),:), &
                 mm4(1:cm(1),1:cm(2),1:cm(3),:), &
                 mm5(1:cm(1),1:cm(2),1:cm(3),:), &
                 mm6(1:cm(1),1:cm(2),1:cm(3),:), &
                 mm7(1:cm(1),1:cm(2),1:cm(3),:), &
                 cm_phi2(1:cm(1),1:cm(2),1:cm(3),:), &
                 cm_src(1:cm(1),1:cm(2),1:cm(3),:))
    call CMFD_POWER
    error = maxval(abs(cm_phi2-cm_phi1)/cm_phi2)
    end do

    else
    !!! FMFD !!!
    if ( .not. cmfdon ) then
    if ( pfmfdon ) then
    call FMFD_DHAT2(i); call MMATRIX2
    else
    call FMFD_DHAT1(i); call MMATRIX1
    end if
    ! -----    
    error=1D0
    do while( error > 1D-8 .and. .not. isnan(k_cmfd) )
    fm_phi1 = fm_phi2
    call SOURCE
    call CG(mm1,mm2,mm3,mm4,mm5,mm6,mm7,fm_phi2,fm_src)
    call POWER(i)
    error = maxval(abs(fm_phi2-fm_phi1)/fm_phi2)
    end do
    else

    !!! CMFD !!! (global-local)
    call FMFD_DHAT1(i)
    call L_BC
    call L_MATRIX
    do
    ! ------------------------------- GLOBAL
    error = 1D0
    call G_DHAT
    call G_MATRIX
    k_pre = k_cmfd
    do global = 1, 5
    cm_phi1 = cm_phi2
    call G_SOURCE
    call CG1(gm1,gm2,gm3,gm4,gm5,gm6,gm7,cm_phi2,cm_src)
    call G_POWER
    end do
    error = abs(k_cmfd-k_pre)/k_cmfd
    if ( error < 1D-8 .or. isnan(k_cmfd) ) exit
    ! ------------------------------- LOCAL
    call G_INJ
    do local = 1, 2
    call G_TO_L
    call L_SOURCE
    call L_SOLVER
    call L_OUTJ
    call L_REFJ
    call G_XS
    end do
    end do
    end if
    end if

    print*, "k_cmfd", k_cmfd

    ! divergence check
    if ( isnan(k_cmfd) ) then
        k_cmfd = 1D0
        p_fsd = 1D0
        if ( mprup ) then
            fmfdup = .true.
            return
        end if
    end if

    if ( i > skip ) then
        p_fsd = 1D0
    else
        if ( .not. ccmfd ) then
        call FISSPROB
        else
        call FISSPROB2
        end if
    end if
    !stop
    
end subroutine

! =============================================================================
! CMFD_DHAT
! ============================================================================= 
subroutine CMFD_DHAT
    implicit none

    do kk = 1, cm(3)
    do jj = 1, cm(2)
    do ii = 1, cm(1)
        ! x0
        if ( ii /= 1 ) then
        cm_dhat(ii,jj,kk,:,1) = &
            (cmJn(ii,jj,kk,:,1)+cm_dtild(ii,jj,kk,:,1) &
            *(cm_phi2(ii,jj,kk,:)-cm_phi2(ii-1,jj,kk,:))) &
            /(cm_phi2(ii,jj,kk,:)+cm_phi2(ii-1,jj,kk,:))
        end if
        ! x1
        if ( ii /= cm(1) ) then
        cm_dhat(ii,jj,kk,:,2) = &
            (cmJn(ii,jj,kk,:,2)+cm_dtild(ii,jj,kk,:,2) &
            *(cm_phi2(ii+1,jj,kk,:)-cm_phi2(ii,jj,kk,:))) &
            /(cm_phi2(ii+1,jj,kk,:)+cm_phi2(ii,jj,kk,:))
        end if
        ! y0
        if ( jj /= 1 ) then
        cm_dhat(ii,jj,kk,:,3) = &
            (cmJn(ii,jj,kk,:,3)+cm_dtild(ii,jj,kk,:,3) &
            *(cm_phi2(ii,jj,kk,:)-cm_phi2(ii,jj-1,kk,:))) &
            /(cm_phi2(ii,jj,kk,:)+cm_phi2(ii,jj-1,kk,:))
        end if
        ! y1
        if ( jj /= cm(2) ) then
        cm_dhat(ii,jj,kk,:,4) = &
            (cmJn(ii,jj,kk,:,4)+cm_dtild(ii,jj,kk,:,4) &
            *(cm_phi2(ii,jj+1,kk,:)-cm_phi2(ii,jj,kk,:))) &
            /(cm_phi2(ii,jj+1,kk,:)+cm_phi2(ii,jj,kk,:))
        end if
        ! z0
        if ( kk /= 1 ) then
        cm_dhat(ii,jj,kk,:,5) = &
            (cmJn(ii,jj,kk,:,5)+cm_dtild(ii,jj,kk,:,5) &
            *(cm_phi2(ii,jj,kk,:)-cm_phi2(ii,jj,kk-1,:))) &
            /(cm_phi2(ii,jj,kk,:)+cm_phi2(ii,jj,kk-1,:))
        end if
        ! z1
        if ( kk /= cm(3) ) then
        cm_dhat(ii,jj,kk,:,6) = &
            (cmJn(ii,jj,kk,:,6)+cm_dtild(ii,jj,kk,:,6) &
            *(cm_phi2(ii,jj,kk+1,:)-cm_phi2(ii,jj,kk,:))) &
            /(cm_phi2(ii,jj,kk+1,:)+cm_phi2(ii,jj,kk,:))
        end if
    end do
    end do
    end do

    ! Boundary condition
    do kk = 1, cm(3)
    do jj = 1, cm(2)
        ii = 1
        cm_dhat(ii,jj,kk,:,1) = cmJn(ii,jj,kk,:,1)/cm_phi2(ii,jj,kk,:)
        ii = cm(1)
        cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
    end do
    end do
    do kk = 1, cm(3)
    do ii = 1, cm(1)
        jj = 1
        cm_dhat(ii,jj,kk,:,3) = cmJn(ii,jj,kk,:,3)/cm_phi2(ii,jj,kk,:)
        jj = cm(2)
        cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
    end do
    end do
    do ii = 1, cm(1)
    do jj = 1, cm(2)
        kk = 1
        cm_dhat(ii,jj,kk,:,5) = cmJn(ii,jj,kk,:,5)/cm_phi2(ii,jj,kk,:)
        kk = cm(3)
        cm_dhat(ii,jj,kk,:,6) = cmJn(ii,jj,kk,:,6)/cm_phi2(ii,jj,kk,:)
    end do
    end do

end subroutine

! =============================================================================
! MMATRIX3
! =============================================================================    
subroutine MMATRIX3
    implicit none
    real(8):: hcm

    ! cell components
    do kk = 1, cm(3); hcm = sum(hfm((kk-1)*fcz+fm0(3):kk*fcz+fm0(3)-1))
    do jj = 1, cm(2)
    do ii = 1, cm(1)
        ! conventional FDM
        if ( kk /= 1 ) &
        mm1(ii,jj,kk,:) = -(cm_dtild(ii,jj,kk,:,5)+cm_dhat(ii,jj,kk,:,5))/hcm
        if ( jj /= 1 ) &
        mm2(ii,jj,kk,:) = -(cm_dtild(ii,jj,kk,:,3)+cm_dhat(ii,jj,kk,:,3))/(pfm*fcr)
        if ( ii /= 1 ) &
        mm3(ii,jj,kk,:) = -(cm_dtild(ii,jj,kk,:,1)+cm_dhat(ii,jj,kk,:,1))/(pfm*fcr)
        if ( ii /= cm(1) ) &
        mm5(ii,jj,kk,:) = -(cm_dtild(ii,jj,kk,:,2)-cm_dhat(ii,jj,kk,:,2))/(pfm*fcr)
        if ( jj /= cm(2) ) &
        mm6(ii,jj,kk,:) = -(cm_dtild(ii,jj,kk,:,4)-cm_dhat(ii,jj,kk,:,4))/(pfm*fcr)
        if ( kk /= cm(3) ) &   
        mm7(ii,jj,kk,:) = -(cm_dtild(ii,jj,kk,:,6)-cm_dhat(ii,jj,kk,:,6))/hcm
        
        mm4(ii,jj,kk,:)= &
            +(cm_dtild(ii,jj,kk,:,1)-cm_dhat(ii,jj,kk,:,1))/(pfm*fcr) &
            +(cm_dtild(ii,jj,kk,:,2)+cm_dhat(ii,jj,kk,:,2))/(pfm*fcr) &
            +(cm_dtild(ii,jj,kk,:,3)-cm_dhat(ii,jj,kk,:,3))/(pfm*fcr) &
            +(cm_dtild(ii,jj,kk,:,4)+cm_dhat(ii,jj,kk,:,4))/(pfm*fcr) &
            +(cm_dtild(ii,jj,kk,:,5)-cm_dhat(ii,jj,kk,:,5))/hcm &
            +(cm_dtild(ii,jj,kk,:,6)+cm_dhat(ii,jj,kk,:,6))/hcm &
            +cm_remv(ii,jj,kk,:)
    end do
    end do
    end do

    !!! zigzag !!!
    if ( zigzag ) then
    ! out of domain
    where ( cm_phi2(1:cm(1),1:cm(1),1:cm(3),:) == 0 )
        mm1(1:cm(1),1:cm(2),1:cm(3),:) = 0D0
        mm2(1:cm(1),1:cm(2),1:cm(3),:) = 0D0
        mm3(1:cm(1),1:cm(2),1:cm(3),:) = 0D0
        mm4(1:cm(1),1:cm(2),1:cm(3),:) = 1D0
        mm5(1:cm(1),1:cm(2),1:cm(3),:) = 0D0
        mm6(1:cm(1),1:cm(2),1:cm(3),:) = 0D0
        mm7(1:cm(1),1:cm(2),1:cm(3),:) = 0D0
    end where

    ! boundary surface
    if ( gmp(1) == 0 .and. gmp(2) == 0 ) then

        do kk = 1, cm(3)
        ! -----
        ii = cm(1)-gzz(2)
        do jj = cm(2)-gzz(1)+1, cm(2)
            mm5(ii,jj,kk,:) = 0D0
        end do
        ! -----
        ii = cm(1)-gzz(1)
        do jj = cm(2)-gzz(2)+1, cm(2)-gzz(1)
            mm5(ii,jj,kk,:) = 0D0
        end do
        ! -----
        jj = cm(2)-gzz(2)
        do ii = cm(1)-gzz(1)+1, cm(1)
            mm6(ii,jj,kk,:) = 0D0
        end do
        ! -----
        jj = cm(1)-gzz(1)
        do ii = cm(1)-gzz(2)+1, cm(1)-gzz(1)
            mm6(ii,jj,kk,:) = 0D0
        end do
        end do

    ! whole core
    else
        do ii = 1, cm(1); id(1) = abs(nint(ii-gmp(1)))
        do jj = 1, cm(2); id(2) = abs(nint(jj-gmp(2)))
        ! -----
        if ( id(1) == cm(1)/2-gzz(2) .and. &
            cm(2)/2-gzz(1) < id(2) .and. id(2) <= cm(2)/2 ) then
            if ( nint(ii-gmp(1)) < 0 ) then
                do kk = 1, cm(3)
                mm3(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = 1, cm(3)
                mm5(ii,jj,kk,:) = 0D0
                end do
            end if
        end if
        ! -----
        if ( id(1) == cm(1)/2-gzz(1) .and. &
            cm(2)/2-gzz(2) < id(2) .and. id(2) <= cm(2)/2-gzz(1) ) then
            if ( nint(ii-gmp(1)) < 0 ) then
                do kk = 1, cm(3)
                mm3(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = 1, cm(3)
                mm5(ii,jj,kk,:) = 0D0
                end do
            end if
        end if
        ! -----
        if ( id(2) == cm(2)/2-gzz(2) .and. &
            cm(1)/2-gzz(1) < id(1) .and. id(1) <= cm(1)/2 ) then
            if ( nint(jj-gmp(2)) < 0 ) then
                do kk = 1, cm(3)
                mm2(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = 1, cm(3)
                mm6(ii,jj,kk,:) = 0D0
                end do
            end if
        end if
        ! -----
        if ( id(2) == cm(2)/2-gzz(1) .and. &
            cm(1)/2-gzz(2) < id(1) .and. id(1) <= cm(1)/2-gzz(1) ) then
            if ( nint(jj-gmp(2)) < 0 ) then
                do kk = 1, cm(3)
                mm2(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = 1, cm(3)
                mm6(ii,jj,kk,:) = 0D0
                end do
            end if
        end if

        end do
        end do
    end if
    end if

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine CMFD_SOURCE
    implicit none
    real(8):: fsource

    ! fis1sion source
    do ii = 1, cm(1)
    do jj = 1, cm(2)
    do kk = 1, cm(3)
        fsource = sum(cm_nufiss(ii,jj,kk,:)*cm_phi1(ii,jj,kk,:))
    do ee=1, cm_eng
        cm_src(ii,jj,kk,ee) = cm_chi(ii,jj,kk,ee)*fsource/k_cmfd
    end do
    end do
    end do
    end do

    ! scattering source
    do ii=1, cm_eng
    do jj=1, cm_eng
    if ( ii /= jj ) then
        cm_src(:,:,:,ii) = &
        cm_src(:,:,:,ii) + cm_scat(:,:,:,jj,ii)*cm_phi1(:,:,:,jj)
    end if
    end do
    end do

end subroutine
    
! =============================================================================
! POWER updates multiplication factor by the power method
! =============================================================================
subroutine CMFD_POWER
    implicit none
    real(8):: sum1, sum2, vol
    
    ! k update
    sum1 = 0; sum2 = 0;
    do kk = 1, cm(3)
        vol = (pfm*fcr)*(pfm*fcr)*sum(hfm((kk-1)*fcz+fm0(3):kk*fcz+fm0(3)-1))
        sum1 = sum1 + sum(cm_nufiss(:,:,kk,:)*cm_phi1(:,:,kk,:))*vol
        sum2 = sum2 + sum(cm_nufiss(:,:,kk,:)*cm_phi2(:,:,kk,:))*vol
    end do
    k_cmfd = k_cmfd*sum2/sum1

end subroutine

! =============================================================================
! FMFD_DHAT1 calculates correction factors for the FMFD
! ============================================================================= 
subroutine FMFD_DHAT1(i)
    use CMFD,          only: D_BC
    implicit none
    integer, intent(in):: i

    do kk = fm0(3), fm1(3)
    do jj = fm0(2), fm1(2)
    do ii = fm0(1), fm1(1)
        ! x0
        if ( ii /= fm0(1) ) then
        fm_dhat(ii,jj,kk,:,1) = &
            (fmJn(ii,jj,kk,:,1)+dtild(ii,jj,kk,:,1) &
            *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii-1,jj,kk,:))) &
            /(fm_phi2(ii,jj,kk,:)+fm_phi2(ii-1,jj,kk,:))
        end if
        ! x1
        if ( ii /= fm1(1) ) then
        fm_dhat(ii,jj,kk,:,2) = &
            (fmJn(ii,jj,kk,:,2)+dtild(ii,jj,kk,:,2) &
            *(fm_phi2(ii+1,jj,kk,:)-fm_phi2(ii,jj,kk,:))) &
            /(fm_phi2(ii+1,jj,kk,:)+fm_phi2(ii,jj,kk,:))
        end if
        ! y0
        if ( jj /= fm0(2) ) then
        fm_dhat(ii,jj,kk,:,3) = &
            (fmJn(ii,jj,kk,:,3)+dtild(ii,jj,kk,:,3) &
            *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj-1,kk,:))) &
            /(fm_phi2(ii,jj,kk,:)+fm_phi2(ii,jj-1,kk,:))
        end if
        ! y1
        if ( jj /= fm1(2) ) then
        fm_dhat(ii,jj,kk,:,4) = &
            (fmJn(ii,jj,kk,:,4)+dtild(ii,jj,kk,:,4) &
            *(fm_phi2(ii,jj+1,kk,:)-fm_phi2(ii,jj,kk,:))) &
            /(fm_phi2(ii,jj+1,kk,:)+fm_phi2(ii,jj,kk,:))
        end if
        ! z0
        if ( kk /= fm0(3) ) then
        fm_dhat(ii,jj,kk,:,5) = &
            (fmJn(ii,jj,kk,:,5)+dtild(ii,jj,kk,:,5) &
            *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj,kk-1,:))) &
            /(fm_phi2(ii,jj,kk,:)+fm_phi2(ii,jj,kk-1,:))
        end if
        ! z1
        if ( kk /= fm1(3) ) then
        fm_dhat(ii,jj,kk,:,6) = &
            (fmJn(ii,jj,kk,:,6)+dtild(ii,jj,kk,:,6) &
            *(fm_phi2(ii,jj,kk+1,:)-fm_phi2(ii,jj,kk,:))) &
            /(fm_phi2(ii,jj,kk+1,:)+fm_phi2(ii,jj,kk,:))
        end if
    end do
    end do
    end do

    if ( unbiased ) call DHAT_UNBIASED1(i)
    if ( cmfdon ) call D_BC

    ! Boundary condition
    do kk=fm0(3), fm1(3)
    do jj=fm0(2), fm1(2)
        ii = fm0(1)
        fm_dhat(ii,jj,kk,:,1) = fmJn(ii,jj,kk,:,1)/fm_phi2(ii,jj,kk,:)
        ii = fm1(1)
        fm_dhat(ii,jj,kk,:,2) = fmJn(ii,jj,kk,:,2)/fm_phi2(ii,jj,kk,:)
    end do
    end do
    do kk=fm0(3), fm1(3)
    do ii=fm0(1), fm1(1)
        jj = fm0(2)
        fm_dhat(ii,jj,kk,:,3) = fmJn(ii,jj,kk,:,3)/fm_phi2(ii,jj,kk,:)
        jj = fm1(2)
        fm_dhat(ii,jj,kk,:,4) = fmJn(ii,jj,kk,:,4)/fm_phi2(ii,jj,kk,:)
    end do
    end do
    do ii=fm0(1), fm1(1)
    do jj=fm0(2), fm1(2)
        kk = fm0(3)
        fm_dhat(ii,jj,kk,:,5) = fmJn(ii,jj,kk,:,5)/fm_phi2(ii,jj,kk,:)
        kk = fm1(3)
        fm_dhat(ii,jj,kk,:,6) = fmJn(ii,jj,kk,:,6)/fm_phi2(ii,jj,kk,:)
    end do
    end do

    if ( unbiased ) call ALBEDO_UNBIASED1(i)

end subroutine


! =============================================================================    
! DHAT_UNBIASED1
! =============================================================================    
subroutine DHAT_UNBIASED1(i)
    use STATISTICS,     only: BIAS
    implicit none
    integer, intent(in):: i
    integer:: i1, i2

    call L_ACCUM(i,i1,i2)
    call DHAT_INOUT1(i1,i2)

!    if ( i2 == 5 ) then
!    kk = fm0(3)
!    do jj = fm1(2), fm0(2), -1
!        write(8,1), (fm_dhat(ii,jj,kk,ee,2), ii = fm0(1), fm1(1))
!    end do
!    write(8,*)
!    end if
!    1 format(1000es15.7)

    ! unbiased estimator
    do kk = fm0(3), fm1(3)
    do jj = fm0(2), fm1(2)
    do ii = fm0(1), fm1(1)
    do ee = 1, cm_eng
        ! x0
        if ( ii /= fm0(1) ) then
        if ( fm_phi2(ii-1,jj,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
                - BIAS(ure0(:,ii,jj,kk,ee,1),ure1(:,ii,jj,kk,ee,1))
        end if
        end if
        ! y0
        if ( jj /= fm0(2) ) then
        if ( fm_phi2(ii,jj-1,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
                - BIAS(ure0(:,ii,jj,kk,ee,3),ure1(:,ii,jj,kk,ee,3))
        end if
        end if
        ! z0
        if ( kk /= fm0(3) ) then
        if ( fm_phi2(ii,jj,kk-1,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,5) = fm_dhat(ii,jj,kk,ee,5) &
                - BIAS(ure0(:,ii,jj,kk,ee,5),ure1(:,ii,jj,kk,ee,5))
        end if
        end if
        ! x1
        if ( ii /= fm1(1) ) then
        if ( fm_phi2(ii+1,jj,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                - BIAS(ure0(:,ii,jj,kk,ee,2),ure1(:,ii,jj,kk,ee,2))
        end if
        end if
        ! y1
        if ( jj /= fm1(2) ) then
        if ( fm_phi2(ii,jj+1,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                - BIAS(ure0(:,ii,jj,kk,ee,4),ure1(:,ii,jj,kk,ee,4))
        end if
        end if
        ! z1
        if ( kk /= fm1(3) ) then
        if ( fm_phi2(ii,jj,kk+1,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,6) = fm_dhat(ii,jj,kk,ee,6) &
                - BIAS(ure0(:,ii,jj,kk,ee,6),ure1(:,ii,jj,kk,ee,6))
        end if
        end if
    end do
    end do
    end do
    end do

!    if ( i2 == 5 ) then
!    kk = fm0(3)
!    do jj = fm1(2), fm0(2), -1
!        write(8,1), (fm_dhat(ii,jj,kk,ee,2), ii = fm0(1), fm1(1))
!    end do
!    stop
!    end if

end subroutine

! =============================================================================
! DHAT_INOUT1
! =============================================================================
subroutine DHAT_INOUT1(i1,i2)
    integer, intent(in):: i1, i2
    integer:: cc    ! cycle index

    ! cycle length
    if ( i2 < accum ) then
        cc = i2
    else
        cc = accum
    end if
    
    ! update
    do ee = 1, cm_eng
    do ii = fm0(1), fm1(1)
    do jj = fm0(2), fm1(2)
    do kk = fm0(3), fm1(3)
        ! x0
        if ( ii /= fm0(1) ) then
            ure0(1:cc,ii,jj,kk,ee,1) = (fmJ1(i1:i2,ii,jj,kk,ee,1) &
                -fmJ0(i1:i2,ii,jj,kk,ee,1))+dtild(ii,jj,kk,ee,1) &
                *(fm_phi(i1:i2,ii,jj,kk,ee)-fm_phi(i1:i2,ii-1,jj,kk,ee))
            ure1(1:cc,ii,jj,kk,ee,1) = &
                fm_phi(i1:i2,ii,jj,kk,ee)+fm_phi(i1:i2,ii-1,jj,kk,ee)
        end if
        ! x1
        if ( ii /= fm1(1) ) then
            ure0(1:cc,ii,jj,kk,ee,2) = (fmJ1(i1:i2,ii,jj,kk,ee,2) &
                -fmJ0(i1:i2,ii,jj,kk,ee,2))+dtild(ii,jj,kk,ee,2) &
                *(fm_phi(i1:i2,ii+1,jj,kk,ee)-fm_phi(i1:i2,ii,jj,kk,ee))
            ure1(1:cc,ii,jj,kk,ee,2) = &
                fm_phi(i1:i2,ii+1,jj,kk,ee)+fm_phi(i1:i2,ii,jj,kk,ee)
        end if
        ! y0
        if ( jj /= fm0(2) ) then
            ure0(1:cc,ii,jj,kk,ee,3) = (fmJ1(i1:i2,ii,jj,kk,ee,3) &
                -fmJ0(i1:i2,ii,jj,kk,ee,3))+dtild(ii,jj,kk,ee,3) &
                *(fm_phi(i1:i2,ii,jj,kk,ee)-fm_phi(i1:i2,ii,jj-1,kk,ee))
            ure1(1:cc,ii,jj,kk,ee,3) = &
                fm_phi(i1:i2,ii,jj,kk,ee)+fm_phi(i1:i2,ii,jj-1,kk,ee)
        end if
        ! y1
        if ( jj /= fm1(2) ) then
            ure0(1:cc,ii,jj,kk,ee,4) = (fmJ1(i1:i2,ii,jj,kk,ee,4) &
                -fmJ0(i1:i2,ii,jj,kk,ee,4))+dtild(ii,jj,kk,ee,4) &
                *(fm_phi(i1:i2,ii,jj+1,kk,ee)-fm_phi(i1:i2,ii,jj,kk,ee))
            ure1(1:cc,ii,jj,kk,ee,4) = &
                fm_phi(i1:i2,ii,jj+1,kk,ee)+fm_phi(i1:i2,ii,jj,kk,ee)
        end if
        ! z0
        if ( kk /= fm0(3) ) then
            ure0(1:cc,ii,jj,kk,ee,5) = (fmJ1(i1:i2,ii,jj,kk,ee,5) &
                -fmJ0(i1:i2,ii,jj,kk,ee,5))+dtild(ii,jj,kk,ee,5) &
                *(fm_phi(i1:i2,ii,jj,kk,ee)-fm_phi(i1:i2,ii,jj,kk-1,ee))
            ure1(1:cc,ii,jj,kk,ee,5) = &
                fm_phi(i1:i2,ii,jj,kk,ee)+fm_phi(i1:i2,ii,jj,kk-1,ee)
        end if
        ! z1
        if ( kk /= fm1(3) ) then
            ure0(1:cc,ii,jj,kk,ee,6) = (fmJ1(i1:i2,ii,jj,kk,ee,6) &
                -fmJ0(i1:i2,ii,jj,kk,ee,6))+dtild(ii,jj,kk,ee,6) &
                *(fm_phi(i1:i2,ii,jj,kk+1,ee)-fm_phi(i1:i2,ii,jj,kk,ee))
            ure1(1:cc,ii,jj,kk,ee,6) = &
                fm_phi(i1:i2,ii,jj,kk+1,ee)+fm_phi(i1:i2,ii,jj,kk,ee)
        end if

    end do
    end do
    end do
    end do

!    ! data out
!    do ii = accum-1, 1, -1
!        ure0(ii+1,:,:,:,:,:) = ure0(ii,:,:,:,:,:)
!        ure1(ii+1,:,:,:,:,:) = ure1(ii,:,:,:,:,:)
!    end do
!
!    ! data in (update)
!    do ii = fm0(1), fm1(1)
!    do jj = fm0(2), fm1(2)
!    do kk = fm0(3), fm1(3)
!        ! x0
!        if ( ii /= fm0(1) ) then
!            ure0(1,ii,jj,kk,:,1) = fmJn(ii,jj,kk,:,1)+dtild(ii,jj,kk,:,1) &
!                *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii-1,jj,kk,:))
!            ure1(1,ii,jj,kk,:,1) = fm_phi2(ii,jj,kk,:)+fm_phi2(ii-1,jj,kk,:)
!        end if
!        ! x1
!        if ( ii /= fm1(1) ) then
!            ure0(1,ii,jj,kk,:,2) = fmJn(ii,jj,kk,:,2)+dtild(ii,jj,kk,:,2) &
!                *(fm_phi2(ii+1,jj,kk,:)-fm_phi2(ii,jj,kk,:))
!            ure1(1,ii,jj,kk,:,2) = fm_phi2(ii+1,jj,kk,:)+fm_phi2(ii,jj,kk,:)
!        end if
!        ! y0
!        if ( jj /= fm0(2) ) then
!            ure0(1,ii,jj,kk,:,3) = fmJn(ii,jj,kk,:,3)+dtild(ii,jj,kk,:,3) &
!                *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj-1,kk,:))
!            ure1(1,ii,jj,kk,:,3) = fm_phi2(ii,jj,kk,:)+fm_phi2(ii,jj-1,kk,:)
!        end if
!        ! y1
!        if ( jj /= fm1(2) ) then
!            ure0(1,ii,jj,kk,:,4) = fmJn(ii,jj,kk,:,4)+dtild(ii,jj,kk,:,4) &
!                *(fm_phi2(ii,jj+1,kk,:)-fm_phi2(ii,jj,kk,:))
!            ure1(1,ii,jj,kk,:,4) = fm_phi2(ii,jj+1,kk,:)+fm_phi2(ii,jj,kk,:)
!        end if
!        ! z0
!        if ( kk /= fm0(3) ) then
!            ure0(1,ii,jj,kk,:,5) = fmJn(ii,jj,kk,:,5)+dtild(ii,jj,kk,:,5) &
!                *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj,kk-1,:))
!            ure1(1,ii,jj,kk,:,5) = fm_phi2(ii,jj,kk,:)+fm_phi2(ii,jj,kk-1,:)
!        end if
!        ! z1
!        if ( kk /= fm1(3) ) then
!            ure0(1,ii,jj,kk,:,6) = fmJn(ii,jj,kk,:,6)+dtild(ii,jj,kk,:,6) &
!                *(fm_phi2(ii,jj,kk+1,:)-fm_phi2(ii,jj,kk,:))
!            ure1(1,ii,jj,kk,:,6) = fm_phi2(ii,jj,kk+1,:)+fm_phi2(ii,jj,kk,:)
!        end if
!
!    end do
!    end do
!    end do

end subroutine


! =============================================================================    
! ALBEDO_UNBIASED1
! =============================================================================    
subroutine ALBEDO_UNBIASED1(i)
    use STATISTICS,     only: BIAS
    use CMFD_PARA,      only: pfm
    implicit none
    integer, intent(in):: i
    integer:: i1, i2

    call L_ACCUM(i,i1,i2)

    do ee = 1, cm_eng
    do kk=fm0(3), fm1(3)
    do jj=fm0(2), fm1(2)
        ii = fm0(1)
        fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,1)-fmJ0(i1:i2,ii,jj,kk,ee,1) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
        ii = fm1(1)
        fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
    end do
    end do
    do kk=fm0(3), fm1(3)
    do ii=fm0(1), fm1(1)
        jj = fm0(2)
        fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,3)-fmJ0(i1:i2,ii,jj,kk,ee,3) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
        jj = fm1(2)
        fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
    end do
    end do
    do ii=fm0(1), fm1(1)
    do jj=fm0(2), fm1(2)
        kk = fm0(3)
        fm_dhat(ii,jj,kk,ee,5) = fm_dhat(ii,jj,kk,ee,5) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,5)-fmJ0(i1:i2,ii,jj,kk,ee,5) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/hfm(kk))
        kk = fm1(3)
        fm_dhat(ii,jj,kk,ee,6) = fm_dhat(ii,jj,kk,ee,6) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,6)-fmJ0(i1:i2,ii,jj,kk,ee,6) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/hfm(kk))
    end do
    end do
    end do

    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            do ee = 1, cm_eng
            do kk = fm0(3), fm1(3)
            ! -----
            ii = fm1(1)-zz(2)
            do jj = fm1(2)-zz(1)+1, fm1(2)
            fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            ! -----
            ii = fm1(1)-zz(1)
            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
            fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            ! -----
            jj = fm1(2)-zz(2)
            do ii = fm1(1)-zz(1)+1, fm1(1)
            fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            ! -----
            jj = fm1(2)-zz(1)
            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
            fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            end do
            end do


        ! whole core
        else
            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
            do ee = 1, cm_eng

            if ( id(1) == afm(1)/2-zz(2) .and. &
                afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,1)-fmJ0(i1:i2,ii,jj,kk,ee,1) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            ! -----
            if ( id(1) == afm(1)/2-zz(1) .and. &
                afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,1)-fmJ0(i1:i2,ii,jj,kk,ee,1) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(2) .and. &
                afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,3)-fmJ0(i1:i2,ii,jj,kk,ee,3) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(1) .and. &
                afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,3)-fmJ0(i1:i2,ii,jj,kk,ee,3) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            
            end do
            end do
            end do

        end if
    end if

end subroutine



! =============================================================================    
! FMFD_DHAT2 calculates correction factors based on partial currents
! 1/3/5 : negative correction factor
! 2/4/6 : positive correction factor
! =============================================================================    
subroutine FMFD_DHAT2(i)
    implicit none
    integer, intent(in):: i

    do ii = fm0(1), fm1(1)
    do jj = fm0(2), fm1(2)
    do kk = fm0(3), fm1(3)
        ! x0
        if ( ii /= fm0(1) ) &
        fm_dhat(ii,jj,kk,:,1) = &
            (fmJ0(0,ii,jj,kk,:,1)-5D-1*dtild(ii,jj,kk,:,1) &
            *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii-1,jj,kk,:))) &
            /fm_phi2(ii,jj,kk,:)
        ! x1
        if ( ii /= fm1(1) ) &   
        fm_dhat(ii,jj,kk,:,2) = &
            (fmJ1(0,ii,jj,kk,:,2)+5D-1*dtild(ii,jj,kk,:,2) &
            *(fm_phi2(ii+1,jj,kk,:)-fm_phi2(ii,jj,kk,:))) &
            /fm_phi2(ii,jj,kk,:)
        ! y0
        if ( jj /= fm0(2) ) &   
        fm_dhat(ii,jj,kk,:,3) = &
            (fmJ0(0,ii,jj,kk,:,3)-5D-1*dtild(ii,jj,kk,:,3) &
            *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj-1,kk,:))) &
            /fm_phi2(ii,jj,kk,:)
        ! y1
        if ( jj /= fm1(2) ) &   
        fm_dhat(ii,jj,kk,:,4) = &
            (fmJ1(0,ii,jj,kk,:,4)+5D-1*dtild(ii,jj,kk,:,4) &
            *(fm_phi2(ii,jj+1,kk,:)-fm_phi2(ii,jj,kk,:))) &
            /fm_phi2(ii,jj,kk,:)
        ! z0
        if ( kk /= fm0(3) ) &   
        fm_dhat(ii,jj,kk,:,5) = &
            (fmJ0(0,ii,jj,kk,:,5)-5D-1*dtild(ii,jj,kk,:,5) &
            *(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj,kk-1,:))) &
            /fm_phi2(ii,jj,kk,:)
        ! z1
        if ( kk /= fm1(3) ) &   
        fm_dhat(ii,jj,kk,:,6) = &
            (fmJ1(0,ii,jj,kk,:,6)+5D-1*dtild(ii,jj,kk,:,6) &
            *(fm_phi2(ii,jj,kk+1,:)-fm_phi2(ii,jj,kk,:))) &
            /fm_phi2(ii,jj,kk,:)
    end do
    end do
    end do

    if ( unbiased ) call DHAT_UNBIASED2(i)
    if ( zigzag ) call pDHATBC

    ! Boundary condition
    do kk=fm0(3), fm1(3)
    do jj=fm0(2), fm1(2)
        ii = fm0(1)
        fm_dhat(ii,jj,kk,:,1) = -fmJn(ii,jj,kk,:,1)/fm_phi2(ii,jj,kk,:)
        ii = fm1(1)
        fm_dhat(ii,jj,kk,:,2) = +fmJn(ii,jj,kk,:,2)/fm_phi2(ii,jj,kk,:)
    end do
    end do
    do kk=fm0(3), fm1(3)
    do ii=fm0(1), fm1(1)
        jj = fm0(2)
        fm_dhat(ii,jj,kk,:,3) = -fmJn(ii,jj,kk,:,3)/fm_phi2(ii,jj,kk,:)
        jj = fm1(2)
        fm_dhat(ii,jj,kk,:,4) = +fmJn(ii,jj,kk,:,4)/fm_phi2(ii,jj,kk,:)
    end do
    end do
    do ii=fm0(1), fm1(1)
    do jj=fm0(2), fm1(2)
        kk = fm0(3)
        fm_dhat(ii,jj,kk,:,5) = -fmJn(ii,jj,kk,:,5)/fm_phi2(ii,jj,kk,:)
        kk = fm1(3)
        fm_dhat(ii,jj,kk,:,6) = +fmJn(ii,jj,kk,:,6)/fm_phi2(ii,jj,kk,:)
    end do
    end do

    if ( unbiased ) call ALBEDO_UNBIASED2(i)

end subroutine

! =============================================================================    
! DHAT_UNBIASED2
! =============================================================================    
subroutine DHAT_UNBIASED2(i)
    use STATISTICS,     only: BIAS
    implicit none
    integer, intent(in):: i
    integer:: i1, i2

    call L_ACCUM(i,i1,i2)
    call DHAT_INOUT2(i)

    ! unbiased estimator
    do kk = fm0(3), fm1(3)
    do jj = fm0(2), fm1(2)
    do ii = fm0(1), fm1(1)
    do ee = 1, cm_eng
        ! x0
        if ( ii /= fm0(1) ) then
        if ( fm_phi2(ii-1,jj,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
                - BIAS(ure0(:,ii,jj,kk,ee,1),ure1(:,ii,jj,kk,ee,1))
        end if
        end if
        ! y0
        if ( jj /= fm0(2) ) then
        if ( fm_phi2(ii,jj-1,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
                - BIAS(ure0(:,ii,jj,kk,ee,3),ure1(:,ii,jj,kk,ee,3))
        end if
        end if
        ! z0
        if ( kk /= fm0(3) ) then
        if ( fm_phi2(ii,jj,kk-1,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,5) = fm_dhat(ii,jj,kk,ee,5) &
                - BIAS(ure0(:,ii,jj,kk,ee,5),ure1(:,ii,jj,kk,ee,5))
        end if
        end if
        ! x1
        if ( ii /= fm1(1) ) then
        if ( fm_phi2(ii+1,jj,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                - BIAS(ure0(:,ii,jj,kk,ee,2),ure1(:,ii,jj,kk,ee,2))
        end if
        end if
        ! y1
        if ( jj /= fm1(2) ) then
        if ( fm_phi2(ii,jj+1,kk,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                - BIAS(ure0(:,ii,jj,kk,ee,4),ure1(:,ii,jj,kk,ee,4))
        end if
        end if
        ! z1
        if ( kk /= fm1(3) ) then
        if ( fm_phi2(ii,jj,kk+1,ee) /= 0 ) then
            fm_dhat(ii,jj,kk,ee,6) = fm_dhat(ii,jj,kk,ee,6) &
                - BIAS(ure0(:,ii,jj,kk,ee,6),ure1(:,ii,jj,kk,ee,6))
        end if
        end if
    end do
    end do
    end do
    end do

end subroutine

! =============================================================================
! DHAT_INOUT2
! =============================================================================
subroutine DHAT_INOUT2(i)
    implicit none
    integer, intent(in):: i
    
    ! data out
    do ii = accum-1, 1, -1
        ure0(ii+1,:,:,:,:,:) = ure0(ii,:,:,:,:,:)
        ure1(ii+1,:,:,:,:,:) = ure1(ii,:,:,:,:,:)
    end do

    ! data in (update)
    do ii = fm0(1), fm1(1)
    do jj = fm0(2), fm1(2)
    do kk = fm0(3), fm1(3)
        ! x0
        if ( ii /= fm0(1) ) then
            ure0(1,ii,jj,kk,:,1) = fmJ0(i,ii,jj,kk,:,1)-5D-1 &
                *dtild(ii,jj,kk,:,1)*(fm_phi(i,ii,jj,kk,:)-fm_phi(i,ii-1,jj,kk,:))
            ure1(1,ii,jj,kk,:,1) = fm_phi(i,ii,jj,kk,:)
        end if
        ! x1
        if ( ii /= fm1(1) ) then
            ure0(1,ii,jj,kk,:,2) = fmJ1(i,ii,jj,kk,:,2)+5D-1 &
                *dtild(ii,jj,kk,:,2)*(fm_phi(i,ii+1,jj,kk,:)-fm_phi(i,ii,jj,kk,:))
            ure1(1,ii,jj,kk,:,2) = fm_phi(i,ii,jj,kk,:)
        end if
        ! y0
        if ( jj /= fm0(2) ) then
            ure0(1,ii,jj,kk,:,3) = fmJ0(i,ii,jj,kk,:,3)-5D-1 &
                *dtild(ii,jj,kk,:,3)*(fm_phi(i,ii,jj,kk,:)-fm_phi(i,ii,jj-1,kk,:))
            ure1(1,ii,jj,kk,:,3) = fm_phi(i,ii,jj,kk,:)
        end if
        ! y1
        if ( jj /= fm1(2) ) then
            ure0(1,ii,jj,kk,:,4) = fmJ1(i,ii,jj,kk,:,4)+5D-1 &
                *dtild(ii,jj,kk,:,4)*(fm_phi(i,ii,jj+1,kk,:)-fm_phi(i,ii,jj,kk,:))
            ure1(1,ii,jj,kk,:,4) = fm_phi(i,ii,jj,kk,:)
        end if
        ! z0
        if ( kk /= fm0(3) ) then
            ure0(1,ii,jj,kk,:,5) = fmJ0(i,ii,jj,kk,:,5)-5D-1 &
                *dtild(ii,jj,kk,:,5)*(fm_phi(i,ii,jj,kk,:)-fm_phi(i,ii,jj,kk-1,:))
            ure1(1,ii,jj,kk,:,5) = fm_phi(i,ii,jj,kk,:)
        end if
        ! z1
        if ( kk /= fm1(3) ) then
            ure0(1,ii,jj,kk,:,6) = fmJ1(i,ii,jj,kk,:,6)+5D-1 &
                *dtild(ii,jj,kk,:,6)*(fm_phi(i,ii,jj,kk+1,:)-fm_phi(i,ii,jj,kk,:))
            ure1(1,ii,jj,kk,:,6) = fm_phi(i,ii,jj,kk,:)
        end if

    end do
    end do
    end do

!    ! data in (update)
!    do ii = fm0(1), fm1(1)
!    do jj = fm0(2), fm1(2)
!    do kk = fm0(3), fm1(3)
!        ! x0
!        if ( ii /= fm0(1) ) then
!            ure0(1,ii,jj,kk,:,1) = fmJ0(0,ii,jj,kk,:,1)-5D-1 &
!                *dtild(ii,jj,kk,:,1)*(fm_phi2(ii,jj,kk,:)-fm_phi2(ii-1,jj,kk,:))
!            ure1(1,ii,jj,kk,:,1) = fm_phi2(ii,jj,kk,:)
!        end if
!        ! x1
!        if ( ii /= fm1(1) ) then
!            ure0(1,ii,jj,kk,:,2) = fmJ1(0,ii,jj,kk,:,2)+5D-1 &
!                *dtild(ii,jj,kk,:,2)*(fm_phi2(ii+1,jj,kk,:)-fm_phi2(ii,jj,kk,:))
!            ure1(1,ii,jj,kk,:,2) = fm_phi2(ii,jj,kk,:)
!        end if
!        ! y0
!        if ( jj /= fm0(2) ) then
!            ure0(1,ii,jj,kk,:,3) = fmJ0(0,ii,jj,kk,:,3)-5D-1 &
!                *dtild(ii,jj,kk,:,3)*(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj-1,kk,:))
!            ure1(1,ii,jj,kk,:,3) = fm_phi2(ii,jj,kk,:)
!        end if
!        ! y1
!        if ( jj /= fm1(2) ) then
!            ure0(1,ii,jj,kk,:,4) = fmJ1(0,ii,jj,kk,:,4)+5D-1 &
!                *dtild(ii,jj,kk,:,4)*(fm_phi2(ii,jj+1,kk,:)-fm_phi2(ii,jj,kk,:))
!            ure1(1,ii,jj,kk,:,4) = fm_phi2(ii,jj,kk,:)
!        end if
!        ! z0
!        if ( kk /= fm0(3) ) then
!            ure0(1,ii,jj,kk,:,5) = fmJ0(0,ii,jj,kk,:,5)-5D-1 &
!                *dtild(ii,jj,kk,:,5)*(fm_phi2(ii,jj,kk,:)-fm_phi2(ii,jj,kk-1,:))
!            ure1(1,ii,jj,kk,:,5) = fm_phi2(ii,jj,kk,:)
!        end if
!        ! z1
!        if ( kk /= fm1(3) ) then
!            ure0(1,ii,jj,kk,:,6) = fmJ1(0,ii,jj,kk,:,6)+5D-1 &
!                *dtild(ii,jj,kk,:,6)*(fm_phi2(ii,jj,kk+1,:)-fm_phi2(ii,jj,kk,:))
!            ure1(1,ii,jj,kk,:,6) = fm_phi2(ii,jj,kk,:)
!        end if
!
!    end do
!    end do
!    end do

end subroutine


! =============================================================================    
! ALBEDO_UNBIASED2
! =============================================================================    
subroutine ALBEDO_UNBIASED2(i)
    use STATISTICS,     only: BIAS
    use CMFD_PARA,      only: pfm
    implicit none
    integer, intent(in):: i
    integer:: i1, i2

    call L_ACCUM(i,i1,i2)

    do ee = 1, cm_eng
    do kk=fm0(3), fm1(3)
    do jj=fm0(2), fm1(2)
        ii = fm0(1)
        fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
            -BIAS(fmJ0(i1:i2,ii,jj,kk,ee,1)-fmJ1(i1:i2,ii,jj,kk,ee,1) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
        ii = fm1(1)
        fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
    end do
    end do
    do kk=fm0(3), fm1(3)
    do ii=fm0(1), fm1(1)
        jj = fm0(2)
        fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
            -BIAS(fmJ0(i1:i2,ii,jj,kk,ee,3)-fmJ1(i1:i2,ii,jj,kk,ee,3) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
        jj = fm1(2)
        fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
    end do
    end do
    do ii=fm0(1), fm1(1)
    do jj=fm0(2), fm1(2)
        kk = fm0(3)
        fm_dhat(ii,jj,kk,ee,5) = fm_dhat(ii,jj,kk,ee,5) &
            -BIAS(fmJ0(i1:i2,ii,jj,kk,ee,5)-fmJ1(i1:i2,ii,jj,kk,ee,5) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/hfm(kk))
        kk = fm1(3)
        fm_dhat(ii,jj,kk,ee,6) = fm_dhat(ii,jj,kk,ee,6) &
            -BIAS(fmJ1(i1:i2,ii,jj,kk,ee,6)-fmJ0(i1:i2,ii,jj,kk,ee,6) &
            ,fm_phi(i1:i2,ii,jj,kk,ee)/hfm(kk))
    end do
    end do
    end do

    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            do ee = 1, cm_eng
            do kk = fm0(3), fm1(3)
            ! -----
            ii = fm1(1)-zz(2)
            do jj = fm1(2)-zz(1)+1, fm1(2)
            fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            ! -----
            ii = fm1(1)-zz(1)
            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
            fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            ! -----
            jj = fm1(2)-zz(2)
            do ii = fm1(1)-zz(1)+1, fm1(1)
            fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            ! -----
            jj = fm1(2)-zz(1)
            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
            fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
            end do
            end do
            end do


        ! whole core
        else
            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
            do ee = 1, cm_eng

            if ( id(1) == afm(1)/2-zz(2) .and. &
                afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
                   - BIAS(fmJ0(i1:i2,ii,jj,kk,ee,1)-fmJ1(i1:i2,ii,jj,kk,ee,1) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            ! -----
            if ( id(1) == afm(1)/2-zz(1) .and. &
                afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,1) = fm_dhat(ii,jj,kk,ee,1) &
                   - BIAS(fmJ0(i1:i2,ii,jj,kk,ee,1)-fmJ1(i1:i2,ii,jj,kk,ee,1) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,2) = fm_dhat(ii,jj,kk,ee,2) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,2)-fmJ0(i1:i2,ii,jj,kk,ee,2) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(2) .and. &
                afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
                   - BIAS(fmJ0(i1:i2,ii,jj,kk,ee,3)-fmJ1(i1:i2,ii,jj,kk,ee,3) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(1) .and. &
                afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,3) = fm_dhat(ii,jj,kk,ee,3) &
                   - BIAS(fmJ0(i1:i2,ii,jj,kk,ee,3)-fmJ1(i1:i2,ii,jj,kk,ee,3) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,ee,4) = fm_dhat(ii,jj,kk,ee,4) &
                   - BIAS(fmJ1(i1:i2,ii,jj,kk,ee,4)-fmJ0(i1:i2,ii,jj,kk,ee,4) &
                   ,fm_phi(i1:i2,ii,jj,kk,ee)/pfm)
                end do
                end if
            end if
            
            end do
            end do
            end do

        end if
    end if

end subroutine

! =============================================================================
! pDHATBC determines correction factors of p-FMFD for BC
! =============================================================================    
subroutine pDHATBC
    implicit none

    ! quarter core
    if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
        do kk = fm0(3), fm1(3)
        ! -----
        ii = fm1(1)-zz(2)
        do jj = fm1(2)-zz(1)+1, fm1(2)
            fm_dhat(ii,jj,kk,:,2) = +fmJn(ii,jj,kk,:,2)/fm_phi2(ii,jj,kk,:)
        end do
        ! -----
        ii = fm1(1)-zz(1)
        do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
            fm_dhat(ii,jj,kk,:,2) = +fmJn(ii,jj,kk,:,2)/fm_phi2(ii,jj,kk,:)
        end do
        ! -----
        jj = fm1(2)-zz(2)
        do ii = fm1(1)-zz(1)+1, fm1(1)
            fm_dhat(ii,jj,kk,:,4) = +fmJn(ii,jj,kk,:,4)/fm_phi2(ii,jj,kk,:)
        end do
        ! -----
        jj = fm1(1)-zz(1)
        do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
            fm_dhat(ii,jj,kk,:,4) = +fmJn(ii,jj,kk,:,4)/fm_phi2(ii,jj,kk,:)
        end do
        end do

    ! whole core
    else
        do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
        do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
        ! -----
        if ( id(1) == afm(1)/2-zz(2) .and. &
            afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
            if ( nint(ii-mp(1)) < 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,1) = -fmJn(ii,jj,kk,:,1)/fm_phi2(ii,jj,kk,:)
                end do
            else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,2) = +fmJn(ii,jj,kk,:,2)/fm_phi2(ii,jj,kk,:)
                end do
            end if
        end if
        ! -----
        if ( id(1) == afm(1)/2-zz(1) .and. &
            afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
            if ( nint(ii-mp(1)) < 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,1) = -fmJn(ii,jj,kk,:,1)/fm_phi2(ii,jj,kk,:)
                end do
            else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,2) = +fmJn(ii,jj,kk,:,2)/fm_phi2(ii,jj,kk,:)
                end do
            end if
        end if
        ! -----
        if ( id(2) == afm(2)/2-zz(2) .and. &
            afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
            if ( nint(jj-mp(2)) < 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,3) = -fmJn(ii,jj,kk,:,3)/fm_phi2(ii,jj,kk,:)
                end do
            else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,4) = +fmJn(ii,jj,kk,:,4)/fm_phi2(ii,jj,kk,:)
                end do
            end if
        end if
        ! -----
        if ( id(2) == afm(2)/2-zz(1) .and. &
            afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
            if ( nint(jj-mp(2)) < 0 ) then
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,3) = -fmJn(ii,jj,kk,:,3)/fm_phi2(ii,jj,kk,:)
                end do
            else
                do kk = fm0(3), fm1(3)
                fm_dhat(ii,jj,kk,:,4) = +fmJn(ii,jj,kk,:,4)/fm_phi2(ii,jj,kk,:)
                end do
            end if
        end if

        end do
        end do
    end if

end subroutine


! =============================================================================
! MMATRIX1 constructs a migration matrix of coarse mesh structure
! =============================================================================    
subroutine MMATRIX1
    implicit none

    ! cell components
    do ii=fm0(1), fm1(1); id(1)=ii
    do jj=fm0(2), fm1(2); id(2)=jj
    do kk=fm0(3), fm1(3); id(3)=kk
        ! conventional FDM
        if ( kk /= fm0(3) ) &
        mm1(ii,jj,kk,:) = -(dtild(ii,jj,kk,:,5)+fm_dhat(ii,jj,kk,:,5))/hfm(kk)
        if ( jj /= fm0(2) ) &
        mm2(ii,jj,kk,:) = -(dtild(ii,jj,kk,:,3)+fm_dhat(ii,jj,kk,:,3))/pfm
        if ( ii /= fm0(1) ) &
        mm3(ii,jj,kk,:) = -(dtild(ii,jj,kk,:,1)+fm_dhat(ii,jj,kk,:,1))/pfm
        if ( ii /= fm1(1) ) &
        mm5(ii,jj,kk,:) = -(dtild(ii,jj,kk,:,2)-fm_dhat(ii,jj,kk,:,2))/pfm
        if ( jj /= fm1(2) ) &
        mm6(ii,jj,kk,:) = -(dtild(ii,jj,kk,:,4)-fm_dhat(ii,jj,kk,:,4))/pfm
        if ( kk /= fm1(3) ) &   
        mm7(ii,jj,kk,:) = -(dtild(ii,jj,kk,:,6)-fm_dhat(ii,jj,kk,:,6))/hfm(kk)
        
        mm4(ii,jj,kk,:)= &
            +(dtild(ii,jj,kk,:,1)-fm_dhat(ii,jj,kk,:,1))/pfm &
            +(dtild(ii,jj,kk,:,2)+fm_dhat(ii,jj,kk,:,2))/pfm &
            +(dtild(ii,jj,kk,:,3)-fm_dhat(ii,jj,kk,:,3))/pfm &
            +(dtild(ii,jj,kk,:,4)+fm_dhat(ii,jj,kk,:,4))/pfm &
            +(dtild(ii,jj,kk,:,5)-fm_dhat(ii,jj,kk,:,5))/hfm(kk) &
            +(dtild(ii,jj,kk,:,6)+fm_dhat(ii,jj,kk,:,6))/hfm(kk) &
            +fm_remv(ii,jj,kk,:)
    end do
    end do
    end do

    !!! zigzag !!!
    if ( zigzag ) call ZZ_OUT_MATRIX

end subroutine


! =============================================================================
! ZZ_OUT_MATRIX
! =============================================================================
subroutine ZZ_OUT_MATRIX
    implicit none

    ! out of domain
    where ( fm_phi2(:,:,:,:) == 0 )
        mm1(:,:,:,:) = 0D0
        mm2(:,:,:,:) = 0D0
        mm3(:,:,:,:) = 0D0
        mm4(:,:,:,:) = 1D0
        mm5(:,:,:,:) = 0D0
        mm6(:,:,:,:) = 0D0
        mm7(:,:,:,:) = 0D0
    end where

    ! boundary surface
    if ( mp(1) == 0 .and. mp(2) == 0 ) then

        do kk = fm0(3), fm1(3)
        ! -----
        ii = fm1(1)-zz(2)
        do jj = fm1(2)-zz(1)+1, fm1(2)
            mm5(ii,jj,kk,:) = 0D0
        end do
        ! -----
        ii = fm1(1)-zz(1)
        do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
            mm5(ii,jj,kk,:) = 0D0
        end do
        ! -----
        jj = fm1(2)-zz(2)
        do ii = fm1(1)-zz(1)+1, fm1(1)
            mm6(ii,jj,kk,:) = 0D0
        end do
        ! -----
        jj = fm1(1)-zz(1)
        do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
            mm6(ii,jj,kk,:) = 0D0
        end do
        end do

    ! whole core
    else
        do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
        do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
        ! -----
        if ( id(1) == afm(1)/2-zz(2) .and. &
            afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
            if ( nint(ii-mp(1)) < 0 ) then
                do kk = fm0(3), fm1(3)
                mm3(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = fm0(3), fm1(3)
                mm5(ii,jj,kk,:) = 0D0
                end do
            end if
        end if
        ! -----
        if ( id(1) == afm(1)/2-zz(1) .and. &
            afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
            if ( nint(ii-mp(1)) < 0 ) then
                do kk = fm0(3), fm1(3)
                mm3(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = fm0(3), fm1(3)
                mm5(ii,jj,kk,:) = 0D0
                end do
            end if
        end if
        ! -----
        if ( id(2) == afm(2)/2-zz(2) .and. &
            afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
            if ( nint(jj-mp(2)) < 0 ) then
                do kk = fm0(3), fm1(3)
                mm2(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = fm0(3), fm1(3)
                mm6(ii,jj,kk,:) = 0D0
                end do
            end if
        end if
        ! -----
        if ( id(2) == afm(2)/2-zz(1) .and. &
            afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
            if ( nint(jj-mp(2)) < 0 ) then
                do kk = fm0(3), fm1(3)
                mm2(ii,jj,kk,:) = 0D0
                end do
            else
                do kk = fm0(3), fm1(3)
                mm6(ii,jj,kk,:) = 0D0
                end do
            end if
        end if

        end do
        end do
    end if

end subroutine

! =============================================================================
! MMATRIX2 construcs a migration matrix for partial current based coarse mesh
! =============================================================================
subroutine MMATRIX2
    implicit none

    ! Cell components
    do ii=fm0(1), fm1(1); id(1)=ii
    do jj=fm0(2), fm1(2); id(2)=jj
    do kk=fm0(3), fm1(3); id(3)=kk
        ! Conventional FDM
        if ( kk /= fm0(3) ) &
        mm1(ii,jj,kk,:)= &
            -(dtild(ii,jj,kk,:,5)+fm_dhat(ii,jj,kk-1,:,6))/hfm(kk)
        if ( jj /= fm0(2) ) &
        mm2(ii,jj,kk,:)= &
            -(dtild(ii,jj,kk,:,3)+fm_dhat(ii,jj-1,kk,:,4))/pfm
        if ( ii /= fm0(1) ) &
        mm3(ii,jj,kk,:)= &
            -(dtild(ii,jj,kk,:,1)+fm_dhat(ii-1,jj,kk,:,2))/pfm
        if ( ii /= fm1(1) ) &
        mm5(ii,jj,kk,:)= &
            -(dtild(ii,jj,kk,:,2)+fm_dhat(ii+1,jj,kk,:,1))/pfm
        if ( jj /= fm1(2) ) &
        mm6(ii,jj,kk,:)= &
            -(dtild(ii,jj,kk,:,4)+fm_dhat(ii,jj+1,kk,:,3))/pfm
        if ( kk /= fm1(3) ) &
        mm7(ii,jj,kk,:)= &
            -(dtild(ii,jj,kk,:,6)+fm_dhat(ii,jj,kk+1,:,5))/hfm(kk)

        mm4(ii,jj,kk,:)= &
            +(dtild(ii,jj,kk,:,1)+fm_dhat(ii,jj,kk,:,1))/pfm &
            +(dtild(ii,jj,kk,:,2)+fm_dhat(ii,jj,kk,:,2))/pfm &
            +(dtild(ii,jj,kk,:,3)+fm_dhat(ii,jj,kk,:,3))/pfm &
            +(dtild(ii,jj,kk,:,4)+fm_dhat(ii,jj,kk,:,4))/pfm &
            +(dtild(ii,jj,kk,:,5)+fm_dhat(ii,jj,kk,:,5))/hfm(kk) &
            +(dtild(ii,jj,kk,:,6)+fm_dhat(ii,jj,kk,:,6))/hfm(kk) &
            +fm_remv(ii,jj,kk,:)
    end do
    end do
    end do

    !!! zigzag !!!
    if ( zigzag ) call ZZ_OUT_MATRIX

end subroutine

! =============================================================================
! SOURCE calculates the fission source for CMFD calculation
! =============================================================================
subroutine SOURCE
    implicit none
    real(8):: fsource

    ! fission source
    do ii=fm0(1), fm1(1)
    do jj=fm0(2), fm1(2)
    do kk=fm0(3), fm1(3)
        fsource = sum(fm_nufiss(0,ii,jj,kk,:)*fm_phi1(ii,jj,kk,:))
    do ee=1, cm_eng
        fm_src(ii,jj,kk,ee) = fm_chi(0,ii,jj,kk,ee)*fsource/k_cmfd
    end do
    end do
    end do
    end do

    ! scattering source
    do ii=1, cm_eng
    do jj=1, cm_eng
    if ( ii /= jj ) then
        fm_src(:,:,:,ii) = &
        fm_src(:,:,:,ii) + fm_scat(0,:,:,:,jj,ii)*fm_phi1(:,:,:,jj)
    end if
    end do
    end do

end subroutine
    
! =============================================================================
! POWER updates multiplication factor by the power method
! =============================================================================
subroutine POWER(i)
    implicit none
    integer, intent(in):: i
    real(8):: sum1, sum2, vol

    ! k update
    sum1 = 0; sum2 = 0;
    do kk=fm0(3), fm1(3)
        vol = pfm*pfm*hfm(kk)
        sum1 = sum1 + sum(fm_nufiss(0,:,:,kk,:)*fm_phi1(:,:,kk,:))*vol
        sum2 = sum2 + sum(fm_nufiss(0,:,:,kk,:)*fm_phi2(:,:,kk,:))*vol
    end do
    k_cmfd = k_cmfd*sum2/sum1

end subroutine
    
! =============================================================================
! FISSPROB calculates the probability of neutron production in a mesh cell
! =============================================================================    
subroutine FISSPROB
    implicit none
    real(8):: prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng)
    real(8):: prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng)
    real(8):: fsource   ! fission source
    real(8):: vol

    ! fission source distribution (CMFD)
    do ii = fm0(1), fm1(1)
    do jj = fm0(2), fm1(2)
    do kk = fm0(3), fm1(3); vol = pfm*pfm*hfm(kk)
        fsource = sum(fm_nufiss(0,ii,jj,kk,:)*fm_phi2(ii,jj,kk,:))
    do ee=1, cm_eng
        prob2(ii,jj,kk,ee) = fm_chi(0,ii,jj,kk,ee)*fsource*vol
    end do
    end do
    end do
    end do

    ! probability (normalization)
    prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) = &
        mc_fsd(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) &
        /sum(mc_fsd(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng))
    prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) = &
        prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) &
        /sum(prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng))
    where( prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) /= 0 ) &
    p_fsd(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) = &
        prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) &
        /prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng)

end subroutine

! =============================================================================
! FISSPROB calculates the probability of neutron production in a mesh cell
! =============================================================================    
subroutine FISSPROB2
    implicit none
    real(8):: prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng)
    real(8):: prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng)
    real(8):: fsource   ! fission source
    real(8):: vol

    ! flux reading
    do kk = 1, cm(3); id(3) = (kk-1)*fcz+fm0(3)
        vol = (pfm*fcr)*(pfm*fcr)*sum(hfm(id(3):id(3)+fcz-1))
    do jj = 1, cm(2); id(2) = (jj-1)*fcr+fm0(2)
    do ii = 1, cm(1); id(1) = (ii-1)*fcr+fm0(1)
        fsource = sum(cm_nufiss(ii,jj,kk,:)*cm_phi2(ii,jj,kk,:))
    do ee = 1, cm_eng
        prob2(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,ee) = &
        cm_chi(ii,jj,kk,ee)*fsource*vol
    end do
    end do
    end do
    end do

    ! average
    do ee = 1, cm_eng
    do ii = 1, cm(1); id(1) = (ii-1)*fcr+fm0(1)
    do jj = 1, cm(2); id(2) = (jj-1)*fcr+fm0(2)
    do kk = 1, cm(3); id(3) = (kk-1)*fcz+fm0(3)
        mc_fsd(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,ee) = &
        sum(mc_fsd(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,ee)) &
        /(fcr*fcr*fcz)
    end do
    end do
    end do
    end do

    ! probability (normalization)
    prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) = &
        mc_fsd(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) &
        /sum(mc_fsd(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng))
    prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) = &
        prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) &
        /sum(prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng))
    where( prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) /= 0 ) &
    p_fsd(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) = &
        prob2(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng) &
        /prob1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),1:cm_eng)

end subroutine

end module
