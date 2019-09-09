module CMFD
    use CMFD_PARA
    use GEOMETRY_PARA,  only: bc_x0, bc_y0
    implicit none
    integer:: ii, jj, kk, ee, mm, nn, oo
    integer:: id(3)
    integer:: id0(3)
    real(8):: vol
    logical:: zzout
    
    contains

subroutine DP1(var1)
    implicit none
    real(8), intent(in):: var1(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3))

    do kk = fm1(3), fm0(3), -1
    do jj = fm1(2), fm0(2), -1
        write(8,1), (var1(ii,jj,kk), ii = fm0(1), fm1(1))
    end do
    end do
    write(8,*)
    1 format(1000es15.7)

end subroutine
subroutine DP2(var1)
    implicit none
    real(8), intent(in):: var1(1:cm(1),1:cm(2),1:cm(3))

    do kk = cm(3), 1, -1
    do jj = cm(2), 1, -1
        write(8,1), (var1(ii,jj,kk), ii = 1, cm(1))
    end do
    end do
    write(8,*)
    1 format(1000es15.7)

end subroutine


! =============================================================================
! ZIGZAGOUT zigzag
! =============================================================================
subroutine ZIGZAGOUT(ii,jj)
    implicit none
    integer, intent(in):: ii, jj
    integer:: idz(2)

    zzout = .false.

    ! quarter core
    if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
        if ( ii > cm(1)-gzz(2) .and. jj > cm(2)-gzz(2) ) then
        if ( ii > cm(1)-gzz(1) .or.  jj > cm(2)-gzz(1) ) then
            zzout = .true.
        end if
        end if
    
    ! whole core
    else
        idz(1) = abs(nint(ii-gmp(1)))
        idz(2) = abs(nint(jj-gmp(2)))
        if ( idz(1) > cm(1)/2-gzz(2) .and. idz(2) > cm(2)/2-gzz(2) ) then
        if ( idz(1) > cm(1)/2-gzz(1) .or.  idz(2) > cm(2)/2-gzz(1) ) then
            zzout = .true.
        end if
        end if
    end if

end subroutine


! =============================================================================
! D_BC
! =============================================================================
subroutine D_BC
    implicit none

    ! diffusion coefficient at boundary
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
        ! x-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            id(1) = id0(1)+1
            dtild(id(1),id(2),id(3),:,1) = 2D0*fm_ddi(id(1),id(2),id(3),:)
            deltf0(id(1),id(2),id(3),:,1) = 0D0
            id(1) = id0(1)+fcr
            dtild(id(1),id(2),id(3),:,2) = 2D0*fm_ddi(id(1),id(2),id(3),:)
            deltf0(id(1),id(2),id(3),:,2) = 0D0
        end do
        ! y-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
            id(2) = id0(2)+1
            dtild(id(1),id(2),id(3),:,3) = 2D0*fm_ddj(id(1),id(2),id(3),:)
            deltf0(id(1),id(2),id(3),:,3) = 0D0
            id(2) = id0(2)+fcr
            dtild(id(1),id(2),id(3),:,4) = 2D0*fm_ddj(id(1),id(2),id(3),:)
            deltf0(id(1),id(2),id(3),:,4) = 0D0
        end do
        end do
        ! z-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            id(3) = id0(3)+1
            dtild(id(1),id(2),id(3),:,5) = 2D0*fm_ddk(id(1),id(2),id(3),:)
            deltf0(id(1),id(2),id(3),:,5) = 0D0
            id(3) = id0(3)+fcz
            dtild(id(1),id(2),id(3),:,6) = 2D0*fm_ddk(id(1),id(2),id(3),:)
            deltf0(id(1),id(2),id(3),:,6) = 0D0
        end do
        end do
    end do
    end do
    end do

end subroutine


! =============================================================================
! L_BC
! =============================================================================
subroutine L_BC
    implicit none
    real(8):: aa(cm_eng)

    deltf1 = fm_dhat

    ! inner boundary
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1

        if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle

        ! x-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ! x0 --------------------------------------------------------------
            if ( ii /= 1 ) then
            id(1) = id0(1)+1
            deltf1(id(1),id(2),id(3),:,1) = &
                -(dtild(id(1),id(2),id(3),:,1) &
                -fm_dhat(id(1),id(2),id(3),:,1)) &
                /(1D0+2D0*dtild(id(1),id(2),id(3),:,1))
            end if

            ! x1 --------------------------------------------------------------
            if ( ii /= cm(1) ) then
            id(1) = id0(1)+fcr
            deltf1(id(1),id(2),id(3),:,2) = &
                (dtild(id(1),id(2),id(3),:,2) &
                +fm_dhat(id(1),id(2),id(3),:,2)) &
                /(1D0+2D0*dtild(id(1),id(2),id(3),:,2))
            end if
        end do
        end do
        ! y-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ! y0 --------------------------------------------------------------
            if ( jj /= 1 ) then
            id(2) = id0(2)+1
            deltf1(id(1),id(2),id(3),:,3) = &
                -(dtild(id(1),id(2),id(3),:,3) &
                -fm_dhat(id(1),id(2),id(3),:,3)) &
                /(1D0+2D0*dtild(id(1),id(2),id(3),:,3))
            end if

            ! y1 --------------------------------------------------------------
            if ( jj /= cm(2) ) then
            id(2) = id0(2)+fcr
            deltf1(id(1),id(2),id(3),:,4) = &
                (dtild(id(1),id(2),id(3),:,4) &
                +fm_dhat(id(1),id(2),id(3),:,4)) &
                /(1D0+2D0*dtild(id(1),id(2),id(3),:,4))
            end if
        end do
        end do
        ! z-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ! z0 --------------------------------------------------------------
            if ( kk /= 1 ) then
            id(3) = id0(3)+1
            deltf1(id(1),id(2),id(3),:,5) = &
                -(dtild(id(1),id(2),id(3),:,5) &
                -fm_dhat(id(1),id(2),id(3),:,5)) &
                /(1D0+2D0*dtild(id(1),id(2),id(3),:,5))
            end if

            ! z1 --------------------------------------------------------------
            if ( kk /= cm(3) ) then
            id(3) = id0(3)+fcz
            deltf1(id(1),id(2),id(3),:,6) = &
                (dtild(id(1),id(2),id(3),:,6) &
                +fm_dhat(id(1),id(2),id(3),:,6)) &
                /(1D0+2D0*dtild(id(1),id(2),id(3),:,6))
            end if
        end do
        end do
    end do
    end do
    end do

    !!! zigzag !!! BC
    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            ! -----
            do kk = fm0(3), fm1(3)
            ii = fm1(1)-zz(2)
            do jj = fm1(2)-zz(1)+1, fm1(2)
                deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
            end do
            ! -----
            ii = fm1(1)-zz(1)
            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
                deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
            end do
            ! -----
            jj = fm1(2)-zz(2)
            do ii = fm1(1)-zz(1)+1, fm1(1)
                deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
            end do
            ! -----
            jj = fm1(2)-zz(1)
            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
                deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
            end do
            end do

        ! whole core
        else
            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))
            
            if ( id(1) == afm(1)/2-zz(2) .and. &
                 afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,1) = fm_dhat(ii,jj,kk,:,1)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
                    end do
                end if
            end if
            ! -----
            if ( id(1) == afm(1)/2-zz(1) .and. &
                 afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,1) = fm_dhat(ii,jj,kk,:,1)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,2) = fm_dhat(ii,jj,kk,:,2)
                    end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(2) .and. &
                 afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,3) = fm_dhat(ii,jj,kk,:,3)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
                    end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(1) .and. &
                 afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,3) = fm_dhat(ii,jj,kk,:,3)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        deltf1(ii,jj,kk,:,4) = fm_dhat(ii,jj,kk,:,4)
                    end do
                end if
            end if

            end do
            end do
        end if
    end if

end subroutine

! =============================================================================
! L_MATRIX
! =============================================================================
subroutine L_MATRIX
    use SOLVER,        only: CG2
    implicit none
    real(8):: deno(cm_eng)

    ! Matrix formulation
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1

    !!! zigzag !!! : out of domain
    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle
    
    ! -------------------------------------------------------------------------
    ! migration term
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
    do oo = 1, fcz; id(3) = id0(3)+oo
    
    if ( oo /= 1   ) mm1(id(1),id(2),id(3),:) = &
      -(deltf0(id(1),id(2),id(3),:,5)+deltf1(id(1),id(2),id(3),:,5))/hfm(id(3))
    if ( nn /= 1   ) mm2(id(1),id(2),id(3),:) = &
      -(deltf0(id(1),id(2),id(3),:,3)+deltf1(id(1),id(2),id(3),:,3))/pfm
    if ( mm /= 1   ) mm3(id(1),id(2),id(3),:) = &
      -(deltf0(id(1),id(2),id(3),:,1)+deltf1(id(1),id(2),id(3),:,1))/pfm
    if ( mm /= fcr ) mm5(id(1),id(2),id(3),:) = &
      -(deltf0(id(1),id(2),id(3),:,2)-deltf1(id(1),id(2),id(3),:,2))/pfm
    if ( nn /= fcr ) mm6(id(1),id(2),id(3),:) = &
      -(deltf0(id(1),id(2),id(3),:,4)-deltf1(id(1),id(2),id(3),:,4))/pfm
    if ( oo /= fcz ) mm7(id(1),id(2),id(3),:) = &
      -(deltf0(id(1),id(2),id(3),:,6)-deltf1(id(1),id(2),id(3),:,6))/hfm(id(3))
    
    mm4(id(1),id(2),id(3),:) = &
     +(deltf0(id(1),id(2),id(3),:,1)-deltf1(id(1),id(2),id(3),:,1))/pfm & 
     +(deltf0(id(1),id(2),id(3),:,2)+deltf1(id(1),id(2),id(3),:,2))/pfm & 
     +(deltf0(id(1),id(2),id(3),:,3)-deltf1(id(1),id(2),id(3),:,3))/pfm & 
     +(deltf0(id(1),id(2),id(3),:,4)+deltf1(id(1),id(2),id(3),:,4))/pfm & 
     +(deltf0(id(1),id(2),id(3),:,5)-deltf1(id(1),id(2),id(3),:,5))/hfm(id(3)) &
     +(deltf0(id(1),id(2),id(3),:,6)+deltf1(id(1),id(2),id(3),:,6))/hfm(id(3)) &
     +fm_remv(id(1),id(2),id(3),:)

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
        deno(:) = 1D0+2D0*dtild(id(1),id(2),id(3),:,1)
        jsrc(id(1),id(2),id(3),:,1) = 4D0*dtild(id(1),id(2),id(3),:,1)/deno(:)
        fsrc(id(1),id(2),id(3),:,1) = fm_dhat(id(1),id(2),id(3),:,1)/deno(:)
        end if
        ! x1 --------------------------------------------------------------
        if ( ii /= cm(1) ) then
        id(1) = id0(1)+fcr
        deno(:) = 1D0+2D0*dtild(id(1),id(2),id(3),:,2)
        jsrc(id(1),id(2),id(3),:,2) = 4D0*dtild(id(1),id(2),id(3),:,2)/deno(:)
        fsrc(id(1),id(2),id(3),:,2) = fm_dhat(id(1),id(2),id(3),:,2)/deno(:)
        end if
    end do
    end do
    ! y-direction
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ! y0 --------------------------------------------------------------
        if ( jj /= 1 ) then
        id(2) = id0(2)+1
        deno(:) = 1D0+2D0*dtild(id(1),id(2),id(3),:,3)
        jsrc(id(1),id(2),id(3),:,3) = 4D0*dtild(id(1),id(2),id(3),:,3)/deno(:)
        fsrc(id(1),id(2),id(3),:,3) = fm_dhat(id(1),id(2),id(3),:,3)/deno(:)
        end if
        ! y1 --------------------------------------------------------------
        if ( jj /= cm(2) ) then
        id(2) = id0(2)+fcr
        deno(:) = 1D0+2D0*dtild(id(1),id(2),id(3),:,4)
        jsrc(id(1),id(2),id(3),:,4) = 4D0*dtild(id(1),id(2),id(3),:,4)/deno(:)
        fsrc(id(1),id(2),id(3),:,4) = fm_dhat(id(1),id(2),id(3),:,4)/deno(:)
        end if
    end do
    end do
    ! z-direction
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ! z0 --------------------------------------------------------------
        if ( kk /= 1 ) then
        id(3) = id0(3)+1
        deno(:) = 1D0+2D0*dtild(id(1),id(2),id(3),:,5)
        jsrc(id(1),id(2),id(3),:,5) = 4D0*dtild(id(1),id(2),id(3),:,5)/deno(:)
        fsrc(id(1),id(2),id(3),:,5) = fm_dhat(id(1),id(2),id(3),:,5)/deno(:)
        end if
        ! z1 --------------------------------------------------------------
        if ( kk /= cm(3) ) then
        id(3) = id0(3)+fcz
        deno(:) = 1D0+2D0*dtild(id(1),id(2),id(3),:,6)
        jsrc(id(1),id(2),id(3),:,6) = 4D0*dtild(id(1),id(2),id(3),:,6)/deno(:)
        fsrc(id(1),id(2),id(3),:,6) = fm_dhat(id(1),id(2),id(3),:,6)/deno(:)
        end if
    end do
    end do

    end do
    end do
    end do

    !!! zigzag !!!
    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            ! -----
            do kk = fm0(3), fm1(3)
            ii = fm1(1)-zz(2)
            do jj = fm1(2)-zz(1)+1, fm1(2)
                jsrc(ii,jj,kk,:,2) = 0D0
                fsrc(ii,jj,kk,:,2) = 0D0
            end do
            ! -----
            ii = fm1(1)-zz(1)
            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
                jsrc(ii,jj,kk,:,2) = 0D0
                fsrc(ii,jj,kk,:,2) = 0D0
            end do
            ! -----
            jj = fm1(2)-zz(2)
            do ii = fm1(1)-zz(1)+1, fm1(1)
                jsrc(ii,jj,kk,:,4) = 0D0
                fsrc(ii,jj,kk,:,4) = 0D0
            end do
            ! -----
            jj = fm1(2)-zz(1)
            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
                jsrc(ii,jj,kk,:,4) = 0D0
                fsrc(ii,jj,kk,:,4) = 0D0
            end do
            end do

        ! whole core
        else
            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))

            if ( id(1) == afm(1)/2-zz(2) .and. &
                 afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,1) = 0D0
                    fsrc(ii,jj,kk,:,1) = 0D0
                    end do
                else
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,2) = 0D0
                    fsrc(ii,jj,kk,:,2) = 0D0
                    end do
                end if
            end if
            ! -----
            if ( id(1) == afm(1)/2-zz(1) .and. &
                 afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,1) = 0D0
                    fsrc(ii,jj,kk,:,1) = 0D0
                    end do
                else
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,2) = 0D0
                    fsrc(ii,jj,kk,:,2) = 0D0
                    end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(2) .and. &
                 afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,3) = 0D0
                    fsrc(ii,jj,kk,:,3) = 0D0
                    end do
                else
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,4) = 0D0
                    fsrc(ii,jj,kk,:,4) = 0D0
                    end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(1) .and. &
                 afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,3) = 0D0
                    fsrc(ii,jj,kk,:,3) = 0D0
                    end do
                else
                    do kk = fm0(3), fm1(3)
                    jsrc(ii,jj,kk,:,4) = 0D0
                    fsrc(ii,jj,kk,:,4) = 0D0
                    end do
                end if
            end if

            end do
            end do
        end if
    end if

end subroutine

! =============================================================================
! G_DHAT
! =============================================================================
subroutine G_DHAT
    implicit none

    do kk=1, cm(3)
    do jj=1, cm(2)
    do ii=1, cm(1)

    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle

        ! x0
        cm_dhat(ii,jj,kk,:,1) = &
            (cmJn(ii,jj,kk,:,1)+cm_dtild(ii,jj,kk,:,1) &
            *(cm_phi2(ii,jj,kk,:)-cmF(ii,jj,kk,:,1))) &
            /(cm_phi2(ii,jj,kk,:)+cmF(ii,jj,kk,:,1))
        ! x1
        cm_dhat(ii,jj,kk,:,2) = &
            (cmJn(ii,jj,kk,:,2)+cm_dtild(ii,jj,kk,:,2) &
            *(cmF(ii,jj,kk,:,2)-cm_phi2(ii,jj,kk,:))) &
            /(cmF(ii,jj,kk,:,2)+cm_phi2(ii,jj,kk,:))
        ! y0
        cm_dhat(ii,jj,kk,:,3) = &
            (cmJn(ii,jj,kk,:,3)+cm_dtild(ii,jj,kk,:,3) &
            *(cm_phi2(ii,jj,kk,:)-cmF(ii,jj,kk,:,3))) &
            /(cm_phi2(ii,jj,kk,:)+cmF(ii,jj,kk,:,3))
        ! y1
        cm_dhat(ii,jj,kk,:,4) = &
            (cmJn(ii,jj,kk,:,4)+cm_dtild(ii,jj,kk,:,4) &
            *(cmF(ii,jj,kk,:,4)-cm_phi2(ii,jj,kk,:))) &
            /(cmF(ii,jj,kk,:,4)+cm_phi2(ii,jj,kk,:))
        ! z0
        cm_dhat(ii,jj,kk,:,5) = &
            (cmJn(ii,jj,kk,:,5)+cm_dtild(ii,jj,kk,:,5) &
            *(cm_phi2(ii,jj,kk,:)-cmF(ii,jj,kk,:,5))) &
            /(cm_phi2(ii,jj,kk,:)+cmF(ii,jj,kk,:,5))
        ! z1
        cm_dhat(ii,jj,kk,:,6) = &
            (cmJn(ii,jj,kk,:,6)+cm_dtild(ii,jj,kk,:,6) &
            *(cmF(ii,jj,kk,:,6)-cm_phi2(ii,jj,kk,:))) &
            /(cmF(ii,jj,kk,:,6)+cm_phi2(ii,jj,kk,:))
        
    end do
    end do
    end do

    ! boundary condition
    do kk = 1, cm(3)
    do jj = 1, cm(2)
        ii = 1
        cm_dhat(ii,jj,kk,:,1) = cmJn(ii,jj,kk,:,1)/cm_phi2(ii,jj,kk,:)
        ii = cm(1)
        cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
    end do
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


    !!! zigzag !!!
    if  ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            do kk = 1, cm(3)
            ! -----
            ii = cm(1)-gzz(2)
            do jj = cm(2)-gzz(1)+1, cm(2)
                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
            end do
            ! -----
            ii = cm(1)-gzz(1)
            do jj = cm(2)-gzz(2)+1, cm(2)-gzz(1)
                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
            end do
            ! -----
            jj = cm(2)-gzz(2)
            do ii = cm(1)-gzz(1)+1, cm(1)
                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
            end do
            ! -----
            jj = cm(2)-gzz(1)
            do ii = cm(1)-gzz(2)+1, cm(1)-gzz(1)
                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
            end do
            end do
    
        ! whole core
        else
            do ii = 1, cm(1); id(1) = abs(nint(ii-gmp(1)))
            do jj = 1, cm(2); id(2) = abs(nint(jj-gmp(2)))
        
            ! boundary
            if ( id(1) == cm(1)/2-gzz(2) .and. &
                cm(2)/2-gzz(1) < id(2) .and. id(2) <= cm(2)/2 ) then
                if ( nint(ii-gmp(1)) <= 0 ) then
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,1) = cmJn(ii,jj,kk,:,1)/cm_phi2(ii,jj,kk,:)
                end do
                else
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
                end do
                end if
            end if
            ! -----
            if ( id(1) == cm(1)/2-gzz(1) .and. &
                cm(2)/2-gzz(2) < id(2) .and. id(2) <= cm(2)/2-gzz(1) ) then
                if ( nint(ii-gmp(1)) <= 0 ) then
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,1) = cmJn(ii,jj,kk,:,1)/cm_phi2(ii,jj,kk,:)
                end do
                else
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,2) = cmJn(ii,jj,kk,:,2)/cm_phi2(ii,jj,kk,:)
                end do
                end if
            end if
            ! -----
            if ( id(2) == cm(2)/2-gzz(2) .and. &
                cm(1)/2-gzz(1) < id(1) .and. id(1) <= cm(1)/2 ) then
                if ( nint(jj-gmp(2)) <= 0 ) then
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,3) = cmJn(ii,jj,kk,:,3)/cm_phi2(ii,jj,kk,:)
                end do
                else
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
                end do
                end if
            end if
            ! -----
            if ( id(2) == cm(2)/2-gzz(1) .and. &
                cm(1)/2-gzz(2) < id(1) .and. id(1) <= cm(1)/2-gzz(1) ) then
                if ( nint(jj-gmp(2)) <= 0 ) then
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,3) = cmJn(ii,jj,kk,:,3)/cm_phi2(ii,jj,kk,:)
                end do
                else
                do kk = 1, cm(3)
                cm_dhat(ii,jj,kk,:,4) = cmJn(ii,jj,kk,:,4)/cm_phi2(ii,jj,kk,:)
                end do
                end if
            end if
    
            end do
            end do
        end if
    end if

end subroutine


! =============================================================================
! G_MATRIX
! =============================================================================
subroutine G_MATRIX
    implicit none
    real(8):: deno(cm_eng)  ! denominator of the parameter
    real(8):: hcm           ! height of a coarse mesh

    ! diffusion coefficient
    do ii = 1, cm(1)
    do jj = 1, cm(2)
    do kk = 1, cm(3)

    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle

        ! x0
        if ( ii /= 1 ) then
        deno = cm_dtild(ii,jj,kk,:,1)+cm_dtild(ii-1,jj,kk,:,2) &
             + cm_dhat(ii,jj,kk,:,1)-cm_dhat(ii-1,jj,kk,:,2)
        deltc0(ii,jj,kk,:,1) = &
            (cm_dtild(ii,jj,kk,:,1)*cm_dtild(ii-1,jj,kk,:,2) &
            +cm_dhat(ii,jj,kk,:,1)*cm_dhat(ii-1,jj,kk,:,2))/deno
        deltc1(ii,jj,kk,:,1) = &
            (cm_dtild(ii,jj,kk,:,1)*cm_dhat(ii-1,jj,kk,:,2) &
            +cm_dtild(ii-1,jj,kk,:,2)*cm_dhat(ii,jj,kk,:,1))/deno
        end if
        ! x1
        if ( ii /= cm(1) ) then
        deno = cm_dtild(ii+1,jj,kk,:,1)+cm_dtild(ii,jj,kk,:,2) &
             + cm_dhat(ii+1,jj,kk,:,1)-cm_dhat(ii,jj,kk,:,2)
        deltc0(ii,jj,kk,:,2) = &
            (cm_dtild(ii+1,jj,kk,:,1)*cm_dtild(ii,jj,kk,:,2) &
            +cm_dhat(ii+1,jj,kk,:,1)*cm_dhat(ii,jj,kk,:,2))/deno
        deltc1(ii,jj,kk,:,2) = &
            (cm_dtild(ii+1,jj,kk,:,1)*cm_dhat(ii,jj,kk,:,2) &
            +cm_dtild(ii,jj,kk,:,2)*cm_dhat(ii+1,jj,kk,:,1))/deno
        end if
        ! y0
        if ( jj /= 1 ) then
        deno = cm_dtild(ii,jj,kk,:,3)+cm_dtild(ii,jj-1,kk,:,4) &
             + cm_dhat(ii,jj,kk,:,3)-cm_dhat(ii,jj-1,kk,:,4)
        deltc0(ii,jj,kk,:,3) = &
            (cm_dtild(ii,jj,kk,:,3)*cm_dtild(ii,jj-1,kk,:,4) &
            +cm_dhat(ii,jj,kk,:,3)*cm_dhat(ii,jj-1,kk,:,4))/deno
        deltc1(ii,jj,kk,:,3) = &
            (cm_dtild(ii,jj,kk,:,3)*cm_dhat(ii,jj-1,kk,:,4) &
            +cm_dtild(ii,jj-1,kk,:,4)*cm_dhat(ii,jj,kk,:,3))/deno
        end if
        ! y1
        if ( jj /= cm(2) ) then
        deno = cm_dtild(ii,jj+1,kk,:,3)+cm_dtild(ii,jj,kk,:,4) &
             + cm_dhat(ii,jj+1,kk,:,3)-cm_dhat(ii,jj,kk,:,4)
        deltc0(ii,jj,kk,:,4) = &
            (cm_dtild(ii,jj+1,kk,:,3)*cm_dtild(ii,jj,kk,:,4) &
            +cm_dhat(ii,jj+1,kk,:,3)*cm_dhat(ii,jj,kk,:,4))/deno
        deltc1(ii,jj,kk,:,4) = &
            (cm_dtild(ii,jj+1,kk,:,3)*cm_dhat(ii,jj,kk,:,4) &
            +cm_dtild(ii,jj,kk,:,4)*cm_dhat(ii,jj+1,kk,:,3))/deno
        end if
        ! z0
        if ( kk /= 1 ) then
        deno = cm_dtild(ii,jj,kk,:,5)+cm_dtild(ii,jj,kk-1,:,6) &
             + cm_dhat(ii,jj,kk,:,5)-cm_dhat(ii,jj,kk-1,:,6)
        deltc0(ii,jj,kk,:,5) = &
            (cm_dtild(ii,jj,kk,:,5)*cm_dtild(ii,jj,kk-1,:,6) &
            +cm_dhat(ii,jj,kk,:,5)*cm_dhat(ii,jj,kk-1,:,6))/deno
        deltc1(ii,jj,kk,:,5) = &
            (cm_dtild(ii,jj,kk,:,5)*cm_dhat(ii,jj,kk-1,:,6) &
            +cm_dtild(ii,jj,kk-1,:,6)*cm_dhat(ii,jj,kk,:,5))/deno
        end if
        ! z1
        if ( kk /= cm(3) ) then
        deno = cm_dtild(ii,jj,kk+1,:,5)+cm_dtild(ii,jj,kk,:,6) &
             + cm_dhat(ii,jj,kk+1,:,5)-cm_dhat(ii,jj,kk,:,6)
        deltc0(ii,jj,kk,:,6) = &
            (cm_dtild(ii,jj,kk+1,:,5)*cm_dtild(ii,jj,kk,:,6) &
            +cm_dhat(ii,jj,kk+1,:,5)*cm_dhat(ii,jj,kk,:,6))/deno
        deltc1(ii,jj,kk,:,6) = &
            (cm_dtild(ii,jj,kk+1,:,5)*cm_dhat(ii,jj,kk,:,6) &
            +cm_dtild(ii,jj,kk,:,6)*cm_dhat(ii,jj,kk+1,:,5))/deno
        end if
    end do
    end do
    end do
    ! boundary condition (J/phi)
    do kk = 1, cm(3)
    do jj = 1, cm(2)
        ii = 1;     deltc1(ii,jj,kk,:,1) = cm_dhat(ii,jj,kk,:,1)
        ii = cm(1); deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
    end do
    do ii = 1, cm(1)
        jj = 1;     deltc1(ii,jj,kk,:,3) = cm_dhat(ii,jj,kk,:,3)
        jj = cm(2); deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
    end do
    end do
    do ii = 1, cm(1)
    do jj = 1, cm(2)
        kk = 1;     deltc1(ii,jj,kk,:,5) = cm_dhat(ii,jj,kk,:,5)
        kk = cm(3); deltc1(ii,jj,kk,:,6) = cm_dhat(ii,jj,kk,:,6)
    end do
    end do


    !!! zigzag !!!
    if  ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            do kk = 1, cm(3)
            ! -----
            ii = cm(1)-gzz(2)
            do jj = cm(2)-gzz(1)+1, cm(2)
                deltc0(ii,jj,kk,:,2) = 0
                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
            end do
            ! -----
            ii = cm(1)-gzz(1)
            do jj = cm(2)-gzz(2)+1, cm(2)-gzz(1)
                deltc0(ii,jj,kk,:,2) = 0
                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
            end do
            ! -----
            jj = cm(2)-gzz(2)
            do ii = cm(1)-gzz(1)+1, cm(1)
                deltc0(ii,jj,kk,:,4) = 0
                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
            end do
            ! -----
            jj = cm(2)-gzz(1)
            do ii = cm(1)-gzz(2)+1, cm(1)-gzz(1)
                deltc0(ii,jj,kk,:,4) = 0
                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
            end do
            end do
    
        ! whole core
        else
            do ii = 1, cm(1); id(1) = abs(nint(ii-gmp(1)))
            do jj = 1, cm(2); id(2) = abs(nint(jj-gmp(2)))
        
            ! boundary
            if ( id(1) == cm(1)/2-gzz(2) .and. &
                cm(2)/2-gzz(1) < id(2) .and. id(2) <= cm(2)/2 ) then
                if ( nint(ii-gmp(1)) <= 0 ) then
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,1) = 0
                deltc1(ii,jj,kk,:,1) = cm_dhat(ii,jj,kk,:,1)
                end do
                else
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,2) = 0
                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
                end do
                end if
            end if
            ! -----
            if ( id(1) == cm(1)/2-gzz(1) .and. &
                cm(2)/2-gzz(2) < id(2) .and. id(2) <= cm(2)/2-gzz(1) ) then
                if ( nint(ii-gmp(1)) <= 0 ) then
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,1) = 0
                deltc1(ii,jj,kk,:,1) = cm_dhat(ii,jj,kk,:,1)
                end do
                else
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,2) = 0
                deltc1(ii,jj,kk,:,2) = cm_dhat(ii,jj,kk,:,2)
                end do
                end if
            end if
            ! -----
            if ( id(2) == cm(2)/2-gzz(2) .and. &
                cm(1)/2-gzz(1) < id(1) .and. id(1) <= cm(1)/2 ) then
                if ( nint(jj-gmp(2)) <= 0 ) then
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,3) = 0
                deltc1(ii,jj,kk,:,3) = cm_dhat(ii,jj,kk,:,3)
                end do
                else
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,4) = 0
                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
                end do
                end if
            end if
            ! -----
            if ( id(2) == cm(2)/2-gzz(1) .and. &
                cm(1)/2-gzz(2) < id(1) .and. id(1) <= cm(1)/2-gzz(1) ) then
                if ( nint(jj-gmp(2)) <= 0 ) then
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,3) = 0
                deltc1(ii,jj,kk,:,3) = cm_dhat(ii,jj,kk,:,3)
                end do
                else
                do kk = 1, cm(3)
                deltc0(ii,jj,kk,:,4) = 0
                deltc1(ii,jj,kk,:,4) = cm_dhat(ii,jj,kk,:,4)
                end do
                end if
            end if
    
            end do
            end do
        end if
    end if


    ! cell components
    do kk=1, cm(3); hcm = sum(hfm((kk-1)*fcz+fm0(3):kk*fcz+fm0(3)-1))
    do jj=1, cm(2)
    do ii=1, cm(1)
        ! conventional FDM
        if ( kk /= 1 )     gm1(ii,jj,kk,:) = &
            -(deltc0(ii,jj,kk,:,5)+deltc1(ii,jj,kk,:,5))/hcm
        if ( jj /= 1 )     gm2(ii,jj,kk,:) = &
            -(deltc0(ii,jj,kk,:,3)+deltc1(ii,jj,kk,:,3))/(pfm*fcr)
        if ( ii /= 1 )     gm3(ii,jj,kk,:) = &
            -(deltc0(ii,jj,kk,:,1)+deltc1(ii,jj,kk,:,1))/(pfm*fcr)
        if ( ii /= cm(1) ) gm5(ii,jj,kk,:) = &
            -(deltc0(ii,jj,kk,:,2)-deltc1(ii,jj,kk,:,2))/(pfm*fcr)
        if ( jj /= cm(2) ) gm6(ii,jj,kk,:) = &
            -(deltc0(ii,jj,kk,:,4)-deltc1(ii,jj,kk,:,4))/(pfm*fcr)
        if ( kk /= cm(3) ) gm7(ii,jj,kk,:) = &
            -(deltc0(ii,jj,kk,:,6)-deltc1(ii,jj,kk,:,6))/hcm
        
        gm4(ii,jj,kk,:)= &
            +(deltc0(ii,jj,kk,:,1)-deltc1(ii,jj,kk,:,1))/(pfm*fcr) &
            +(deltc0(ii,jj,kk,:,2)+deltc1(ii,jj,kk,:,2))/(pfm*fcr) &
            +(deltc0(ii,jj,kk,:,3)-deltc1(ii,jj,kk,:,3))/(pfm*fcr) &
            +(deltc0(ii,jj,kk,:,4)+deltc1(ii,jj,kk,:,4))/(pfm*fcr) &
            +(deltc0(ii,jj,kk,:,5)-deltc1(ii,jj,kk,:,5))/hcm &
            +(deltc0(ii,jj,kk,:,6)+deltc1(ii,jj,kk,:,6))/hcm &
            +cm_remv(ii,jj,kk,:)

    end do
    end do
    end do

    !!! zigzag !!!
    if ( zigzag ) then
    ! out of domain
    where ( cm_phi2(:,:,:,:) == 0 )
        gm1(:,:,:,:) = 0D0
        gm2(:,:,:,:) = 0D0
        gm3(:,:,:,:) = 0D0
        gm4(:,:,:,:) = 1D0
        gm5(:,:,:,:) = 0D0
        gm6(:,:,:,:) = 0D0
        gm7(:,:,:,:) = 0D0
    end where
    end if

    !!! zigzag !!!
    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            do kk = 1, cm(3)
            ! -----
            ii = cm(1)-gzz(2)
            do jj = cm(2)-gzz(1)+1, cm(2)
                gm5(ii,jj,kk,:) = 0
            end do
            ! -----
            ii = cm(1)-gzz(1)
            do jj = cm(2)-gzz(2)+1, cm(2)-gzz(1)
                gm5(ii,jj,kk,:) = 0
            end do
            ! -----
            jj = cm(2)-gzz(2)
            do ii = cm(1)-gzz(1)+1, cm(1)
                gm6(ii,jj,kk,:) = 0
            end do
            ! -----
            jj = cm(2)-gzz(1)
            do ii = cm(1)-gzz(2)+1, cm(1)-gzz(1)
                gm6(ii,jj,kk,:) = 0
            end do
            end do
    
        ! whole core
        else
            do ii = 1, cm(1); id(1) = abs(nint(ii-gmp(1)))
            do jj = 1, cm(2); id(2) = abs(nint(jj-gmp(2)))
        
            ! boundary
            if ( id(1) == cm(1)/2-gzz(2) .and. &
                cm(2)/2-gzz(1) < id(2) .and. id(2) <= cm(2)/2 ) then
                if ( nint(ii-gmp(1)) <= 0 ) then
                    do kk = 1, cm(3)
                    gm3(ii,jj,kk,:) = 0D0
                    end do
                else
                    do kk = 1, cm(3)
                    gm5(ii,jj,kk,:) = 0D0
                    end do
                end if
            end if
            ! -----
            if ( id(1) == cm(1)/2-gzz(1) .and. &
                cm(2)/2-gzz(2) < id(2) .and. id(2) <= cm(2)/2-gzz(1) ) then
                if ( nint(ii-gmp(1)) <= 0 ) then
                    do kk = 1, cm(3)
                    gm3(ii,jj,kk,:) = 0D0
                    end do
                else
                    do kk = 1, cm(3)
                    gm5(ii,jj,kk,:) = 0D0
                    end do
                end if
            end if
            ! -----
            if ( id(2) == cm(2)/2-gzz(2) .and. &
                cm(1)/2-gzz(1) < id(1) .and. id(1) <= cm(1)/2 ) then
                if ( nint(jj-gmp(2)) <= 0 ) then
                    do kk = 1, cm(3)
                    gm2(ii,jj,kk,:) = 0D0
                    end do
                else
                    do kk = 1, cm(3)
                    gm6(ii,jj,kk,:) = 0D0
                    end do
                end if
            end if
            ! -----
            if ( id(2) == cm(2)/2-gzz(1) .and. &
                cm(1)/2-gzz(2) < id(1) .and. id(1) <= cm(1)/2-gzz(1) ) then
                if ( nint(jj-gmp(2)) <= 0 ) then
                    do kk = 1, cm(3)
                    gm2(ii,jj,kk,:) = 0D0
                    end do
                else
                    do kk = 1, cm(3)
                    gm6(ii,jj,kk,:) = 0D0
                    end do
                end if
            end if
    
            end do
            end do
        end if
    end if

end subroutine

! =============================================================================
! ALBEDO calculates an albedo parameter
! =============================================================================
subroutine ALBEDO(aa,jj1,jj0)
    real(8), intent(in)::  jj1(:), jj0(:)
    real(8), intent(out):: aa(:)
    real(8):: bb(cm_eng)

    do ee = 1, cm_eng
    if ( jj0(ee) == 0 ) then
        if ( jj1(ee) == 0 ) then
            aa(ee) = 0D0
        else
            aa(ee) = 2D0
        end if
    else if ( jj0(ee) == jj1(ee) ) then
        aa(ee) = 0D0
    else
        bb(ee) = jj1(ee)/jj0(ee)
        aa(ee) = 2D0*(bb(ee)+1D0)/(bb(ee)-1D0)
    end if
    end do

end subroutine

! =============================================================================
! G_SOURCE
! =============================================================================
subroutine G_SOURCE
    implicit none
    real(8):: fsource

    ! fission source
    do ii = 1, cm(1)
    do jj = 1, cm(2)
    do kk = 1, cm(3)
        fsource = sum(cm_nufiss(ii,jj,kk,:)*cm_phi1(ii,jj,kk,:))
    do ee = 1, cm_eng
        cm_src(ii,jj,kk,ee) = cm_chi(ii,jj,kk,ee)*fsource/k_cmfd
    end do
    end do
    end do
    end do

    ! scattering source
    do ii = 1, cm_eng
    do jj = 1, cm_eng
    if ( ii /= jj ) then
        cm_src(:,:,:,ii) = &
        cm_src(:,:,:,ii) + cm_scat(:,:,:,jj,ii)*cm_phi1(:,:,:,jj)
    end if
    end do
    end do

end subroutine

! =============================================================================
! G_POWER
! =============================================================================
subroutine G_POWER
    implicit none
    real(8):: sum1, sum2, vol

    ! k update
    sum1 = 0; sum2 = 0;
    do kk = 1, cm(3); id(3) = (kk-1)*fcz+fm0(3)
        vol  = pfm*fcr
        vol  = vol*vol*sum(hfm(id(3):id(3)+fcz-1))
        sum1 = sum1 + sum(cm_nufiss(:,:,kk,:)*cm_phi1(:,:,kk,:))*vol
        sum2 = sum2 + sum(cm_nufiss(:,:,kk,:)*cm_phi2(:,:,kk,:))*vol
    end do
    k_cmfd = k_cmfd*sum2/sum1

end subroutine

! =============================================================================
! G_INJ
! =============================================================================
subroutine G_INJ
    implicit none
    real(8):: deno(cm_eng)  ! denominator of a parameter
    real(8):: surF(cm_eng)

    do kk = 1, cm(3)
    do jj = 1, cm(2)
    do ii = 1, cm(1)

    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle

        if ( ii /= 1 ) then ! x-direction
        deno(:) = cm_dtild(ii,jj,kk,:,1)+cm_dtild(ii-1,jj,kk,:,2) &
                + cm_dhat(ii,jj,kk,:,1)-cm_dhat(ii-1,jj,kk,:,2)
        cmF(ii,jj,kk,:,1) = ((cm_dtild(ii,jj,kk,:,1)-cm_dhat(ii,jj,kk,:,1)) &
            *cm_phi2(ii,jj,kk,:)+(cm_dtild(ii-1,jj,kk,:,2) &
            +cm_dhat(ii-1,jj,kk,:,2))*cm_phi2(ii-1,jj,kk,:))/deno(:)
        cmJn(ii,jj,kk,:,1) = -cm_dtild(ii,jj,kk,:,1) &
                            *(cm_phi2(ii,jj,kk,:)-cmF(ii,jj,kk,:,1)) &
                            +cm_dhat(ii,jj,kk,:,1) &
                            *(cm_phi2(ii,jj,kk,:)+cmF(ii,jj,kk,:,1))
        cmJ0(ii,jj,kk,:,1) = 25D-2*cmF(ii,jj,kk,:,1)-5D-1*cmJn(ii,jj,kk,:,1)
        cmJ1(ii,jj,kk,:,1) = cmJ0(ii,jj,kk,:,1)+cmJn(ii,jj,kk,:,1)

        cmF(ii-1,jj,kk,:,2)  = cmF(ii,jj,kk,:,1)
        cmJn(ii-1,jj,kk,:,2) = cmJn(ii,jj,kk,:,1)
        cmJ0(ii-1,jj,kk,:,2) = cmJ0(ii,jj,kk,:,1)
        cmJ1(ii-1,jj,kk,:,2) = cmJ1(ii,jj,kk,:,1)
        end if
        if ( jj /= 1 ) then ! y-direction
        deno(:) = cm_dtild(ii,jj,kk,:,3)+cm_dtild(ii,jj-1,kk,:,4) &
                + cm_dhat(ii,jj,kk,:,3)-cm_dhat(ii,jj-1,kk,:,4)
        cmF(ii,jj,kk,:,3) = ((cm_dtild(ii,jj,kk,:,3)-cm_dhat(ii,jj,kk,:,3)) &
            *cm_phi2(ii,jj,kk,:)+(cm_dtild(ii,jj-1,kk,:,4) &
            +cm_dhat(ii,jj-1,kk,:,4))*cm_phi2(ii,jj-1,kk,:))/deno(:)
        cmJn(ii,jj,kk,:,3) = -cm_dtild(ii,jj,kk,:,3) &
                            *(cm_phi2(ii,jj,kk,:)-cmF(ii,jj,kk,:,3)) &
                            +cm_dhat(ii,jj,kk,:,3) &
                            *(cm_phi2(ii,jj,kk,:)+cmF(ii,jj,kk,:,3))
        cmJ0(ii,jj,kk,:,3) = 25D-2*cmF(ii,jj,kk,:,3)-5D-1*cmJn(ii,jj,kk,:,3)
        cmJ1(ii,jj,kk,:,3) = cmJ0(ii,jj,kk,:,3)+cmJn(ii,jj,kk,:,3)

        cmF(ii,jj-1,kk,:,4)  = cmF(ii,jj,kk,:,3)
        cmJn(ii,jj-1,kk,:,4) = cmJn(ii,jj,kk,:,3)
        cmJ0(ii,jj-1,kk,:,4) = cmJ0(ii,jj,kk,:,3)
        cmJ1(ii,jj-1,kk,:,4) = cmJ1(ii,jj,kk,:,3)
        end if
        if ( kk /= 1 ) then ! z-direction
        deno(:) = cm_dtild(ii,jj,kk,:,5)+cm_dtild(ii,jj,kk-1,:,6) &
                + cm_dhat(ii,jj,kk,:,5)-cm_dhat(ii,jj,kk-1,:,6)
        cmF(ii,jj,kk,:,5) = ((cm_dtild(ii,jj,kk,:,5)-cm_dhat(ii,jj,kk,:,5)) &
            *cm_phi2(ii,jj,kk,:)+(cm_dtild(ii,jj,kk-1,:,6) &
            +cm_dhat(ii,jj,kk-1,:,6))*cm_phi2(ii,jj,kk-1,:))/deno(:)
        cmJn(ii,jj,kk,:,5) = -cm_dtild(ii,jj,kk,:,5) &
                            *(cm_phi2(ii,jj,kk,:)-cmF(ii,jj,kk,:,5)) &
                            +cm_dhat(ii,jj,kk,:,5) &
                            *(cm_phi2(ii,jj,kk,:)+cmF(ii,jj,kk,:,5))
        cmJ0(ii,jj,kk,:,5) = 25D-2*cmF(ii,jj,kk,:,5)-5D-1*cmJn(ii,jj,kk,:,5)
        cmJ1(ii,jj,kk,:,5) = cmJ0(ii,jj,kk,:,5)+cmJn(ii,jj,kk,:,5)

        cmF(ii,jj,kk-1,:,6)  = cmF(ii,jj,kk,:,5)
        cmJn(ii,jj,kk-1,:,6) = cmJn(ii,jj,kk,:,5)
        cmJ0(ii,jj,kk-1,:,6) = cmJ0(ii,jj,kk,:,5)
        cmJ1(ii,jj,kk-1,:,6) = cmJ1(ii,jj,kk,:,5)
        end if
    end do
    end do
    end do

end subroutine

! =============================================================================
! G_TO_L carries out the flux and current modulation (from GLOBAL to LOCAL)
! =============================================================================
subroutine G_TO_L
    implicit none
    integer:: ssum(1:2)

    fm_phi1(:,:,:,:) = fm_phi2(:,:,:,:)

    do ee = 1, cm_eng
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1

    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle

    ! flux modulation
    fm_phi2(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz,ee) = &
    fm_phi2(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz,ee) &
    /sum(fm_phi2(id0(1)+1:id0(1)+fcr,id0(2)+1:id0(2)+fcr,id0(3)+1:id0(3)+fcz,ee)) &
    *fcr*fcr*fcz*cm_phi2(ii,jj,kk,ee)

    ! partial current modulation
    ! x0
    ssum = 0;       id(1) = id0(1)+1
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJ0(0,id(1),id(2),id(3),ee,1)
        ssum(2) = ssum(2) + fmJ1(0,id(1),id(2),id(3),ee,1)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum(1) /= 0 ) fmJ0(0,id(1),id(2),id(3),ee,1) = &
            fmJ0(0,id(1),id(2),id(3),ee,1) &
            /ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,ee,1)
        if ( ssum(2) /= 0 ) fmJ1(0,id(1),id(2),id(3),ee,1) = &
            fmJ1(0,id(1),id(2),id(3),ee,1) &
            /ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,ee,1)
    end do
    end do
    ! x1
    ssum = 0;       id(1) = id0(1)+fcr
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJ0(0,id(1),id(2),id(3),ee,2)
        ssum(2) = ssum(2) + fmJ1(0,id(1),id(2),id(3),ee,2)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum(1) /= 0 ) fmJ0(0,id(1),id(2),id(3),ee,2) = &
            fmJ0(0,id(1),id(2),id(3),ee,2) &
            /ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,ee,2)
        if ( ssum(2) /= 0 ) fmJ1(0,id(1),id(2),id(3),ee,2) = &
            fmJ1(0,id(1),id(2),id(3),ee,2) &
            /ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,ee,2)
    end do
    end do
    ! y0
    ssum = 0;       id(2) = id0(2)+1
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(1) = ssum(1) + fmJ0(0,id(1),id(2),id(3),ee,3)
        ssum(2) = ssum(2) + fmJ1(0,id(1),id(2),id(3),ee,3)
    end do 
    end do 
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        if ( ssum(1) /= 0 ) fmJ0(0,id(1),id(2),id(3),ee,3) = &
            fmJ0(0,id(1),id(2),id(3),ee,3) &
            /ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,ee,3)
        if ( ssum(2) /= 0 ) fmJ1(0,id(1),id(2),id(3),ee,3) = &
            fmJ1(0,id(1),id(2),id(3),ee,3) &
            /ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,ee,3)
    end do
    end do 
    ! y1
    ssum = 0;       id(2) = id0(2)+fcr
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(1) = ssum(1) + fmJ0(0,id(1),id(2),id(3),ee,4)
        ssum(2) = ssum(2) + fmJ1(0,id(1),id(2),id(3),ee,4)
    end do
    end do
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        if ( ssum(1) /= 0 ) fmJ0(0,id(1),id(2),id(3),ee,4) = &
            fmJ0(0,id(1),id(2),id(3),ee,4) &
            /ssum(1)*fcr*fcz*cmJ0(ii,jj,kk,ee,4)
        if ( ssum(2) /= 0 ) fmJ1(0,id(1),id(2),id(3),ee,4) = &
            fmJ1(0,id(1),id(2),id(3),ee,4) &
            /ssum(2)*fcr*fcz*cmJ1(ii,jj,kk,ee,4)
    end do
    end do
    ! z0
    ssum = 0;       id(3) = id0(3)+1
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJ0(0,id(1),id(2),id(3),ee,5)
        ssum(2) = ssum(2) + fmJ1(0,id(1),id(2),id(3),ee,5)
    end do
    end do
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum(1) /= 0 ) fmJ0(0,id(1),id(2),id(3),ee,5) = &
            fmJ0(0,id(1),id(2),id(3),ee,5) &
            /ssum(1)*fcr*fcr*cmJ0(ii,jj,kk,ee,5)
        if ( ssum(2) /= 0 ) fmJ1(0,id(1),id(2),id(3),ee,5) = &
            fmJ1(0,id(1),id(2),id(3),ee,5) &
            /ssum(2)*fcr*fcr*cmJ1(ii,jj,kk,ee,5)
    end do
    end do
    ! z1
    ssum = 0;       id(3) = id0(3)+fcz
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJ0(0,id(1),id(2),id(3),ee,6)
        ssum(2) = ssum(2) + fmJ1(0,id(1),id(2),id(3),ee,6)
    end do
    end do
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        if ( ssum(1) /= 0 ) fmJ0(0,id(1),id(2),id(3),ee,6) = &
            fmJ0(0,id(1),id(2),id(3),ee,6) &
            /ssum(1)*fcr*fcr*cmJ0(ii,jj,kk,ee,6)
        if ( ssum(2) /= 0 ) fmJ1(0,id(1),id(2),id(3),ee,6) = &
            fmJ1(0,id(1),id(2),id(3),ee,6) &
            /ssum(2)*fcr*fcr*cmJ1(ii,jj,kk,ee,6)
    end do
    end do

    end do
    end do
    end do
    end do

end subroutine


! =============================================================================
! L_SOURCE
! =============================================================================
subroutine L_SOURCE
    implicit none
    real(8):: fsource

    ! neutron source (fission + scattering)
    do kk = fm0(3), fm1(3)
    do jj = fm0(2), fm1(2)
    do ii = fm0(1), fm1(1)

    !!! zigzag !!!
    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            if ( ii > fm1(1)-zz(2) .and. jj > fm1(2)-zz(2) ) then
            if ( ii > fm1(1)-zz(1) .or.  jj > fm1(2)-zz(1) ) then
                cycle
            end if
            end if
        
        ! whole core
        else
            id(1) = abs(nint(ii-mp(1)))
            id(2) = abs(nint(jj-mp(2)))
            if ( id(1) > afm(1)/2-zz(2) .and. id(2) > afm(2)/2-zz(2) ) then
            if ( id(1) > afm(1)/2-zz(1) .or.  id(2) > afm(2)/2-zz(1) ) then
                cycle
            end if
            end if
        end if
    end if

    fsource = sum(fm_nufiss(0,ii,jj,kk,:)*fm_phi2(ii,jj,kk,:))
    do mm = 1, cm_eng
        fm_src(ii,jj,kk,mm) = fm_chi(0,ii,jj,kk,mm)*fsource/k_cmfd
    do nn = 1, cm_eng
        if ( mm /= nn ) fm_src(ii,jj,kk,mm) = fm_src(ii,jj,kk,mm) &
                + fm_scat(0,ii,jj,kk,nn,mm)*fm_phi2(ii,jj,kk,nn)
    end do
    end do
    end do
    end do
    end do

    ! interface BC
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1

    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle

        ! x-direction
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            if ( ii /= 1 ) then
            id(1) = id0(1)+1
            fm_src(id(1),id(2),id(3),:) = fm_src(id(1),id(2),id(3),:) &
                +(jsrc(id(1),id(2),id(3),:,1)*fmJ1(0,id(1),id(2),id(3),:,1) &
                +fsrc(id(1),id(2),id(3),:,1) &
                *fm_phi1(id(1)-1,id(2),id(3),:))/pfm
            end if
            if ( ii /= cm(1) ) then
            id(1) = id0(1)+fcr
            fm_src(id(1),id(2),id(3),:) = fm_src(id(1),id(2),id(3),:) &
                +(jsrc(id(1),id(2),id(3),:,2)*fmJ0(0,id(1),id(2),id(3),:,2) &
                -fsrc(id(1),id(2),id(3),:,2) &
                *fm_phi1(id(1)+1,id(2),id(3),:))/pfm
            end if
        end do
        ! y-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
            if ( jj /= 1 ) &
            id(2) = id0(2)+1
            fm_src(id(1),id(2),id(3),:) = fm_src(id(1),id(2),id(3),:) &
                +(jsrc(id(1),id(2),id(3),:,3)*fmJ1(0,id(1),id(2),id(3),:,3) &
                +fsrc(id(1),id(2),id(3),:,3) &
                *fm_phi1(id(1),id(2)-1,id(3),:))/pfm
            if ( jj /= cm(2) ) &
            id(2) = id0(2)+fcr
            fm_src(id(1),id(2),id(3),:) = fm_src(id(1),id(2),id(3),:) &
                +(jsrc(id(1),id(2),id(3),:,4)*fmJ0(0,id(1),id(2),id(3),:,4) &
                -fsrc(id(1),id(2),id(3),:,4) &
                *fm_phi1(id(1),id(2)+1,id(3),:))/pfm
        end do
        end do
        ! z-direction
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            if ( kk /= 1 ) &
            id(3) = id0(3)+1
            fm_src(id(1),id(2),id(3),:) = fm_src(id(1),id(2),id(3),:) &
                +(jsrc(id(1),id(2),id(3),:,5)*fmJ1(0,id(1),id(2),id(3),:,5) &
                +fsrc(id(1),id(2),id(3),:,5) &
                *fm_phi1(id(1),id(2),id(3)-1,:))/hfm(id(3))
            if ( kk /= cm(3) ) &
            id(3) = id0(3)+fcz
            fm_src(id(1),id(2),id(3),:) = fm_src(id(1),id(2),id(3),:) &
                +(jsrc(id(1),id(2),id(3),:,6)*fmJ0(0,id(1),id(2),id(3),:,6) &
                -fsrc(id(1),id(2),id(3),:,6) &
                *fm_phi1(id(1),id(2),id(3)+1,:))/hfm(id(3))
        end do
        end do
    end do
    end do
    end do

end subroutine

! =============================================================================
! L_SOLVER
! =============================================================================
subroutine L_SOLVER
    use SOLVER, only: CG2
    implicit none

    do ee = 1, 2
    ! Matrix formulation
    do kk = 1, cm(3); id(3) = (kk-1)*fcz+fm0(3)
    do jj = 1, cm(2); id(2) = (jj-1)*fcr+fm0(2)
    do ii = 1, cm(1); id(1) = (ii-1)*fcr+fm0(1)

    ! zigzag : out of domain
    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle

    ! Solver
    call CG2(mm1(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             mm2(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             mm3(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             mm4(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             mm5(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             mm6(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             mm7(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             fm_phi2(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:), &
             fm_src(id(1):id(1)+fcr-1,id(2):id(2)+fcr-1,id(3):id(3)+fcz-1,:),ii,jj,kk)


    end do
    end do
    end do
    end do

end subroutine

! =============================================================================
! L_OUTJ
! =============================================================================
subroutine L_OUTJ
    implicit none
    real(8):: netJ(cm_eng)

    ! outgoing partial current
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
        ! x0
        id(1) = id0(1)+1
        if ( ii == 1 ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),:,1) = deltf1(id(1),id(2),id(3),:,1) &
                *fm_phi2(id(1),id(2),id(3),:)
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = +jsrc(id(1),id(2),id(3),:,1)*fmJ1(0,id(1),id(2),id(3),:,1) &
                +deltf1(id(1),id(2),id(3),:,1)*fm_phi2(id(1),id(2),id(3),:) &
                +fsrc(id(1),id(2),id(3),:,1)*fm_phi2(id(1)-1,id(2),id(3),:)
            fmF(0,id(1),id(2),id(3),:,1) = &
                4D0*fmJ1(0,id(1),id(2),id(3),:,1)-2D0*netJ
            fmJ0(0,id(1),id(2),id(3),:,1) = &
                25D-2*fmF(0,id(1),id(2),id(3),:,1)-5D-1*netJ
        end do
        end do
        end if
        ! x1
        id(1) = id0(1)+fcr
        if ( ii == cm(1) ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),:,2) = deltf1(id(1),id(2),id(3),:,2) &
                *fm_phi2(id(1),id(2),id(3),:)
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = -jsrc(id(1),id(2),id(3),:,2)*fmJ0(0,id(1),id(2),id(3),:,2) &
                +deltf1(id(1),id(2),id(3),:,2)*fm_phi2(id(1),id(2),id(3),:) &
                +fsrc(id(1),id(2),id(3),:,2)*fm_phi2(id(1)+1,id(2),id(3),:)
            fmF(0,id(1),id(2),id(3),:,2) = &
                4D0*fmJ0(0,id(1),id(2),id(3),:,2)+2D0*netJ
            fmJ1(0,id(1),id(2),id(3),:,2) = &
                25D-2*fmF(0,id(1),id(2),id(3),:,2)+5D-1*netJ
        end do
        end do
        end if
        ! y0
        id(2) = id0(2)+1
        if ( jj == 1 ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJn(id(1),id(2),id(3),:,3) = deltf1(id(1),id(2),id(3),:,3) &
                *fm_phi2(id(1),id(2),id(3),:)
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            netJ = +jsrc(id(1),id(2),id(3),:,3)*fmJ1(0,id(1),id(2),id(3),:,3) &
                +deltf1(id(1),id(2),id(3),:,3)*fm_phi2(id(1),id(2),id(3),:) &
                +fsrc(id(1),id(2),id(3),:,3)*fm_phi2(id(1),id(2)-1,id(3),:)
            fmF(0,id(1),id(2),id(3),:,3) = &
                4D0*fmJ1(0,id(1),id(2),id(3),:,3)-2D0*netJ
            fmJ0(0,id(1),id(2),id(3),:,3) = &
                25D-2*fmF(0,id(1),id(2),id(3),:,3)-5D-1*netJ
        end do
        end do
        end if
        ! y1
        id(2) = id0(2)+fcr
        if ( jj == cm(2) ) then
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJn(id(1),id(2),id(3),:,4) = deltf1(id(1),id(2),id(3),:,4) &
                *fm_phi2(id(1),id(2),id(3),:)
        end do
        end do
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            netJ = -jsrc(id(1),id(2),id(3),:,4)*fmJ0(0,id(1),id(2),id(3),:,4) &
                +deltf1(id(1),id(2),id(3),:,4)*fm_phi2(id(1),id(2),id(3),:) &
                +fsrc(id(1),id(2),id(3),:,4)*fm_phi2(id(1),id(2)+1,id(3),:)
            fmF(0,id(1),id(2),id(3),:,4) = &
                4D0*fmJ0(0,id(1),id(2),id(3),:,4)+2D0*netJ
            fmJ1(0,id(1),id(2),id(3),:,4) = &
                25D-2*fmF(0,id(1),id(2),id(3),:,4)+5D-1*netJ
        end do
        end do
        end if
        ! z0
        id(3) = id0(3)+1
        if ( kk == 1 ) then
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),:,5) = deltf1(id(1),id(2),id(3),:,5) &
                *fm_phi2(id(1),id(2),id(3),:)
        end do
        end do
        else
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = +jsrc(id(1),id(2),id(3),:,5)*fmJ1(0,id(1),id(2),id(3),:,5) &
                +deltf1(id(1),id(2),id(3),:,5)*fm_phi2(id(1),id(2),id(3),:) &
                +fsrc(id(1),id(2),id(3),:,5)*fm_phi2(id(1),id(2),id(3)-1,:)
            fmF(0,id(1),id(2),id(3),:,5) = &
                4D0*fmJ1(0,id(1),id(2),id(3),:,5)-2D0*netJ
            fmJ0(0,id(1),id(2),id(3),:,5) = &
                25D-2*fmF(0,id(1),id(2),id(3),:,5)-5D-1*netJ
        end do
        end do
        end if
        ! z1
        id(3) = id0(3)+fcz
        if ( kk == cm(3) ) then
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJn(id(1),id(2),id(3),:,6) = deltf1(id(1),id(2),id(3),:,6) &
                *fm_phi2(id(1),id(2),id(3),:)
        end do
        end do
        else
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            netJ = -jsrc(id(1),id(2),id(3),:,6)*fmJ0(0,id(1),id(2),id(3),:,6) &
                +deltf1(id(1),id(2),id(3),:,6)*fm_phi2(id(1),id(2),id(3),:) &
                +fsrc(id(1),id(2),id(3),:,6)*fm_phi2(id(1),id(2),id(3)+1,:)
            fmF(0,id(1),id(2),id(3),:,6) = &
                4D0*fmJ0(0,id(1),id(2),id(3),:,6)+2D0*netJ
            fmJ1(0,id(1),id(2),id(3),:,6) = &
                25D-2*fmF(0,id(1),id(2),id(3),:,6)+5D-1*netJ
        end do
        end do
        end if
    end do
    end do
    end do

    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
        ! x0
        id(1) = id0(1)+1
        if ( ii == 1 ) then
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ0(0,id(1)-1,id(2),id(3),:,2) = fmJ0(0,id(1),id(2),id(3),:,1)
        end do
        end do
        end if
        ! x1
        id(1) = id0(1)+fcr
        if ( ii == cm(1) ) then
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ1(0,id(1)+1,id(2),id(3),:,1) = fmJ1(0,id(1),id(2),id(3),:,2)
        end do
        end do
        end if
        ! y0
        id(2) = id0(2)+1
        if ( jj == 1 ) then
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJ0(0,id(1),id(2)-1,id(3),:,4) = fmJ0(0,id(1),id(2),id(3),:,3)
        end do
        end do
        end if
        ! y1
        id(2) = id0(2)+fcr
        if ( jj == cm(2) ) then
        else
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            fmJ1(0,id(1),id(2)+1,id(3),:,3) = fmJ1(0,id(1),id(2),id(3),:,4)
        end do
        end do
        end if
        ! z0
        id(3) = id0(3)+1
        if ( kk == 1 ) then
        else
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ0(0,id(1),id(2),id(3)-1,:,6) = fmJ0(0,id(1),id(2),id(3),:,5)
        end do
        end do
        end if
        ! z1
        id(3) = id0(3)+fcz
        if ( kk == cm(3) ) then
        else
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            fmJ1(0,id(1),id(2),id(3)+1,:,5) = fmJ1(0,id(1),id(2),id(3),:,6)
        end do
        end do
        end if
    end do
    end do
    end do


    !!! zigzag !!!
    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            do kk = fm0(3), fm1(3)
            ! -----
            ii = fm1(1)-zz(2)
            do jj = fm1(2)-zz(1)+1, fm1(2)
                fmJn(ii,jj,kk,:,2) = deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
            end do
            ! -----
            ii = fm1(1)-zz(1)
            do jj = fm1(2)-zz(2)+1, fm1(2)-zz(1)
                fmJn(ii,jj,kk,:,2) = deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
            end do
            ! -----
            jj = fm1(2)-zz(2)
            do ii = fm1(1)-zz(1)+1, fm1(1)
                fmJn(ii,jj,kk,:,4) = deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
            end do
            ! -----
            jj = fm1(2)-zz(1)
            do ii = fm1(1)-zz(2)+1, fm1(1)-zz(1)
                fmJn(ii,jj,kk,:,4) = deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
            end do
            end do

        ! whole core
        else
            do ii = fm0(1), fm1(1); id(1) = abs(nint(ii-mp(1)))
            do jj = fm0(2), fm1(2); id(2) = abs(nint(jj-mp(2)))

            if ( id(1) == afm(1)/2-zz(2) .and. &
                 afm(2)/2-zz(1) < id(2) .and. id(2) <= afm(2)/2 ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,1) = &
                            deltf1(ii,jj,kk,:,1)*fm_phi2(ii,jj,kk,:)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,2) = &
                            deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
                    end do
                end if
            end if
            ! -----
            if ( id(1) == afm(1)/2-zz(1) .and. &
                 afm(2)/2-zz(2) < id(2) .and. id(2) <= afm(2)/2-zz(1) ) then
                if ( nint(ii-mp(1)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,1) = &
                            deltf1(ii,jj,kk,:,1)*fm_phi2(ii,jj,kk,:)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,2) = &
                            deltf1(ii,jj,kk,:,2)*fm_phi2(ii,jj,kk,:)
                    end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(2) .and. &
                 afm(1)/2-zz(1) < id(1) .and. id(1) <= afm(1)/2 ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,3) = &
                            deltf1(ii,jj,kk,:,3)*fm_phi2(ii,jj,kk,:)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,4) = &
                            deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
                    end do
                end if
            end if
            ! -----
            if ( id(2) == afm(2)/2-zz(1) .and. &
                 afm(1)/2-zz(2) < id(1) .and. id(1) <= afm(1)/2-zz(1) ) then
                if ( nint(jj-mp(2)) <= 0 ) then
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,3) = &
                            deltf1(ii,jj,kk,:,3)*fm_phi2(ii,jj,kk,:)
                    end do
                else
                    do kk = fm0(3), fm1(3)
                        fmJn(ii,jj,kk,:,4) = &
                            deltf1(ii,jj,kk,:,4)*fm_phi2(ii,jj,kk,:)
                    end do
                end if
            end if

            end do
            end do
        end if
    end if

end subroutine


! =============================================================================
! L_REFJ
! =============================================================================
subroutine L_REFJ
    implicit none
    real(8):: ssum(0:2)
    real(8):: fsum(cm_eng)
    real(8):: surF(fm0(1):fm1(1),fm0(2):fm1(2),fm0(3):fm1(3),cm_eng,6)

    ! surface average
    do ee = 1, cm_eng
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
        ! x0
        if ( ii /= 1 ) then
        ssum = 0;       id(1) = id0(1)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(0) = ssum(0) + fmJ0(0,id(1),id(2),id(3),ee,1)
            ssum(1) = ssum(1) + fmJ1(0,id(1),id(2),id(3),ee,1)
            ssum(2) = ssum(2) + fmF(0,id(1),id(2),id(3),ee,1)
        end do
        end do
        cmJ0(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,ee,1) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,ee,1)  = ssum(2) / (fcr*fcz)
        end if
        ! x1
        if ( ii /= cm(1) ) then
        ssum = 0;       id(1) = id0(1)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(0) = ssum(0) + fmJ0(0,id(1),id(2),id(3),ee,2)
            ssum(1) = ssum(1) + fmJ1(0,id(1),id(2),id(3),ee,2)
            ssum(2) = ssum(2) + fmF(0,id(1),id(2),id(3),ee,2)
        end do
        end do
        cmJ0(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,ee,2) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,ee,2)  = ssum(2) / (fcr*fcz)
        end if
        ! y0
        if ( jj /= 1 ) then
        ssum = 0;       id(2) = id0(2)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(0) = ssum(0) + fmJ0(0,id(1),id(2),id(3),ee,3)
            ssum(1) = ssum(1) + fmJ1(0,id(1),id(2),id(3),ee,3)
            ssum(2) = ssum(2) + fmF(0,id(1),id(2),id(3),ee,3)
        end do
        end do
        cmJ0(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,ee,3) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,ee,3)  = ssum(2) / (fcr*fcz)
        end if
        ! y1
        if ( jj /= cm(2) ) then
        ssum = 0;       id(2) = id0(2)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(0) = ssum(0) + fmJ0(0,id(1),id(2),id(3),ee,4)
            ssum(1) = ssum(1) + fmJ1(0,id(1),id(2),id(3),ee,4)
            ssum(2) = ssum(2) + fmF(0,id(1),id(2),id(3),ee,4)
        end do
        end do
        cmJ0(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
        cmJ1(ii,jj,kk,ee,4) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,ee,4)  = ssum(2) / (fcr*fcz)
        end if
        ! z0
        if ( kk /= 1 ) then
        ssum = 0;       id(3) = id0(3)+1
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(0) = ssum(0) + fmJ0(0,id(1),id(2),id(3),ee,5)
            ssum(1) = ssum(1) + fmJ1(0,id(1),id(2),id(3),ee,5)
            ssum(2) = ssum(2) + fmF(0,id(1),id(2),id(3),ee,5)
        end do
        end do
        cmJ0(ii,jj,kk,ee,5) = ssum(0) / (fcr*fcr)
        cmJ1(ii,jj,kk,ee,5) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,ee,5)  = ssum(2) / (fcr*fcr)
        end if
        ! z1
        if ( kk /= cm(3) ) then
        ssum = 0;       id(3) = id0(3)+fcz
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(0) = ssum(0) + fmJ0(0,id(1),id(2),id(3),ee,6)
            ssum(1) = ssum(1) + fmJ1(0,id(1),id(2),id(3),ee,6)
            ssum(2) = ssum(2) + fmF(0,id(1),id(2),id(3),ee,6)
        end do
        end do
        cmJ0(ii,jj,kk,ee,6) = ssum(0) / (fcr*fcr)
        cmJ1(ii,jj,kk,ee,6) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,ee,6)  = ssum(2) / (fcr*fcr)
        end if
    end do
    end do
    end do
    end do


    ! interface surface
    do ii = 1, cm(1)
    do jj = 1, cm(2)
    do kk = 1, cm(3)
        ! x-direction
        if ( ii /= 1 ) then
            cmJn(ii,jj,kk,:,1) = cmJ1(ii-1,jj,kk,:,2)-cmJ0(ii,jj,kk,:,1)
            cmJn(ii-1,jj,kk,:,2) = cmJn(ii,jj,kk,:,1)
            fsum(:) = (cmF(ii,jj,kk,:,1)+cmF(ii-1,jj,kk,:,2))/2D0
            cmF(ii  ,jj,kk,:,1) = fsum(:)
            cmF(ii-1,jj,kk,:,2) = fsum(:)
        end if
        ! y-direction
        if ( jj /= 1 ) then
            cmJn(ii,jj,kk,:,3) = cmJ1(ii,jj-1,kk,:,4)-cmJ0(ii,jj,kk,:,3)
            cmJn(ii,jj-1,kk,:,4) = cmJn(ii,jj,kk,:,3)
            fsum(:) = (cmF(ii,jj,kk,:,3)+cmF(ii,jj-1,kk,:,4))/2D0
            cmF(ii,jj  ,kk,:,3) = fsum(:)
            cmF(ii,jj-1,kk,:,4) = fsum(:)
        end if
        ! z-direction
        if ( kk /= 1 ) then
            cmJn(ii,jj,kk,:,5) = cmJ1(ii,jj,kk-1,:,6)-cmJ0(ii,jj,kk,:,5)
            cmJn(ii,jj,kk-1,:,6) = cmJn(ii,jj,kk,:,5)
            fsum(:) = (cmF(ii,jj,kk,:,5)+cmF(ii,jj,kk-1,:,6))/2D0
            cmF(ii,jj,kk  ,:,5) = fsum(:)
            cmF(ii,jj,kk-1,:,6) = fsum(:)
        end if
    end do
    end do
    end do


    ! boundary surface
    do ee = 1, cm_eng
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    !   x0
    ii = 1; id(1) = fm0(1)
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
    end do
    end do
    cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
    end do
    !   x1
    ii = cm(1); id(1) = fm1(1)
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
    end do
    end do
    cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
    end do
    !   y0
    jj = 1; id(2) = fm0(2)
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
    end do
    end do
    cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
    end do
    !   y1
    jj = cm(2); id(2) = fm1(2)
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1; ssum(0) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
    end do
    end do
    cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
    end do
    end do
    !   z0
    kk = 1; id(3) = fm0(3)
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1; ssum(0) = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,5)
    end do
    end do
    cmJn(ii,jj,kk,ee,5) = ssum(0) / (fcr*fcr)
    end do
    end do
    !   z1
    kk = cm(3); id(3) = fm1(3)
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1; ssum(0) = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,6)
    end do
    end do
    cmJn(ii,jj,kk,ee,6) = ssum(0) / (fcr*fcr)
    end do
    end do
    end do


    !!! zigzag !!!
    if ( zigzag ) then
        ! quarter core
        if ( bc_x0 == -1 .and. bc_y0 == -1 ) then
            do ee = 1, cm_eng
            do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
            ! -----
            ii = cm(1)-gzz(2); id(1)  = fm1(1)-zz(2)
            jj = cm(2);        id0(2) = fm1(2)-zz(1)
            ssum(0) = 0
            do oo = 1, fcz; id(3) = id0(3)+oo
            do nn = 1, fcr; id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
            end do
            end do
            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
            ! -----
            ii = cm(1)-gzz(1); id(1)  = fm1(1)-zz(1)
            jj = cm(2)-gzz(1); id0(2) = fm1(2)-zz(2)
            ssum(0) = 0
            do oo = 1, fcz; id(3) = id0(3)+oo
            do nn = 1, fcr; id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
            end do
            end do
            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
            ! -----
            ii = cm(1);        id0(1) = fm1(1)-zz(1)
            jj = cm(2)-gzz(2); id(2)  = fm1(2)-zz(2)
            ssum(0) = 0
            do oo = 1, fcz; id(3) = id0(3)+oo
            do mm = 1, fcr; id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
            end do
            end do
            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
            ! -----
            ii = cm(1)-gzz(1); id0(1) = fm1(1)-zz(2)
            jj = cm(2)-gzz(1); id(2)  = fm1(2)-zz(1)
            ssum(0) = 0
            do oo = 1, fcz; id(3) = id0(3)+oo
            do mm = 1, fcr; id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
            end do
            end do
            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
            end do
            end do

        ! whole core
        else
            do ee = 1, cm_eng
            ! -----
            ii = 1+gzz(1);     id(1)  = fm0(1)+zz(1)
            jj = 1+gzz(1);     id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
            end do
            end do
            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
            end do
            jj = cm(2)-gzz(1); id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
            end do
            end do
            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
            end do
            ! -----
            ii = 1+gzz(2);     id(1)  = fm0(1)+zz(2)
            jj = 1;            id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
            end do
            end do
            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
            end do
            jj = cm(2);        id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,1)
            end do
            end do
            cmJn(ii,jj,kk,ee,1) = ssum(0) / (fcr*fcz)
            end do
            ! -----
            ii = cm(1)-gzz(1); id(1)  = fm1(1)-zz(1)
            jj = 1+gzz(1);     id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
            end do
            end do
            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
            end do
            jj = cm(2)-gzz(1); id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
            end do
            end do
            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
            end do
            ! -----
            ii = cm(1)-gzz(2); id(1)  = fm1(1)-zz(2)
            jj = 1;            id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
            end do
            end do
            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
            end do
            jj = cm(2);        id0(2) = (jj-1)*fcr+fm0(2)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do nn = 1, fcr;    id(2) = id0(2)+nn
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,2)
            end do
            end do
            cmJn(ii,jj,kk,ee,2) = ssum(0) / (fcr*fcz)
            end do

            ! -----
            jj = 1+gzz(1);     id(2)  = fm0(2)+zz(1)
            ii = 1+gzz(1);     id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
            end do
            end do
            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
            end do
            ii = cm(1)-gzz(1); id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
            end do
            end do
            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
            end do
            ! -----
            jj = 1+gzz(2);     id(2)  = fm0(2)+zz(2)
            ii = 1;            id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
            end do
            end do
            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
            end do
            ii = cm(2);        id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,3)
            end do
            end do
            cmJn(ii,jj,kk,ee,3) = ssum(0) / (fcr*fcz)
            end do
            ! -----
            jj = cm(2)-gzz(1); id(2)  = fm1(2)-zz(1)
            ii = 1+gzz(1);     id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
            end do
            end do
            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
            end do
            ii = cm(1)-gzz(1); id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
            end do
            end do
            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
            end do
            ! -----
            jj = cm(2)-gzz(2); id(2)  = fm1(2)-zz(2)
            ii = 1;            id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
            end do
            end do
            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
            end do
            ii = cm(1);        id0(1) = (ii-1)*fcr+fm0(1)-1
            do kk = 1, cm(3);  id0(3) = (kk-1)*fcz+fm0(3)-1; ssum(0) = 0
            do oo = 1, fcz;    id(3) = id0(3)+oo
            do mm = 1, fcr;    id(1) = id0(1)+mm
                ssum(0) = ssum(0) + fmJn(id(1),id(2),id(3),ee,4)
            end do
            end do
            cmJn(ii,jj,kk,ee,4) = ssum(0) / (fcr*fcz)
            end do
            end do
        end if
    end if

end subroutine
    
    
! =============================================================================
! G_XS produces the flux-volume-weight group constants
! ============================================================================= 
subroutine G_XS
    implicit none
    integer:: id1(3), id2(3)
    real(8):: phisum

    ! cross section
    do ee = 1, cm_eng
    do ii = 1, cm(1); id1(1) = (ii-1)*fcr+fm0(1); id2(1) = id1(1)+fcr-1
    do jj = 1, cm(2); id1(2) = (jj-1)*fcr+fm0(2); id2(2) = id1(2)+fcr-1
    do kk = 1, cm(3); id1(3) = (kk-1)*fcz+fm0(3); id2(3) = id1(3)+fcz-1

    if ( zigzag ) call ZIGZAGOUT(ii,jj); if ( zzout ) cycle
        phisum = sum(fm_phi2(id1(1):id2(1),id1(2):id2(2),id1(3):id2(3),ee))
        cm_tot(ii,jj,kk,ee) = &
        sum(fm_tot(0,id1(1):id2(1),id1(2):id2(2),id1(3):id2(3),ee) &
        *fm_phi2(id1(1):id2(1),id1(2):id2(2),id1(3):id2(3),ee))/phisum
        cm_nufiss(ii,jj,kk,ee) = &
        sum(fm_nufiss(0,id1(1):id2(1),id1(2):id2(2),id1(3):id2(3),ee) &
        *fm_phi2(id1(1):id2(1),id1(2):id2(2),id1(3):id2(3),ee))/phisum
    do mm = 1, cm_eng
        cm_scat(ii,jj,kk,ee,mm) = &
        sum(fm_scat(0,id1(1):id2(1),id1(2):id2(2),id1(3):id2(3),ee,mm) &
        *fm_phi2(id1(1):id2(1),id1(2):id2(2),id1(3):id2(3),ee))/phisum
    end do
    cm_remv(ii,jj,kk,ee) = cm_tot(ii,jj,kk,ee) - cm_scat(ii,jj,kk,ee,ee)
    cm_phi2(ii,jj,kk,ee) = phisum/(fcr*fcr*fcz)
    end do
    end do
    end do
    end do


    ! diffusion coefficient
    !   unit diffusion coefficient
    cm_diff = 1D0/(3D0*cm_tot)
    cm_ddi(:,:,:,:) = cm_diff(:,:,:,:)/(pfm*fcr)
    cm_ddj(:,:,:,:) = cm_diff(:,:,:,:)/(pfm*fcr)
    do kk = 1, cm(3); id(3) = (kk-1)*fcz+fm0(3)
    cm_ddk(:,:,kk,:) = cm_diff(:,:,kk,:)/sum(hfm(id(3):id(3)+fcz-1))
    end do
    !   interface diffusion coefficient
    do ii = 1, cm(1)
    do jj = 1, cm(2)
    do kk = 1, cm(3)
        cm_dtild(ii,jj,kk,:,1) = 2D0*cm_ddi(ii,jj,kk,:)
        cm_dtild(ii,jj,kk,:,2) = 2D0*cm_ddi(ii,jj,kk,:)
        cm_dtild(ii,jj,kk,:,3) = 2D0*cm_ddj(ii,jj,kk,:)
        cm_dtild(ii,jj,kk,:,4) = 2D0*cm_ddj(ii,jj,kk,:)
        cm_dtild(ii,jj,kk,:,5) = 2D0*cm_ddk(ii,jj,kk,:)
        cm_dtild(ii,jj,kk,:,6) = 2D0*cm_ddk(ii,jj,kk,:)
    end do
    end do
    end do

end subroutine

end module

