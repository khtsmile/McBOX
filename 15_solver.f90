module SOLVER
    implicit none
    
    contains
    
! =============================================================================
! CG is a matrix solver by the SOR
! =============================================================================
subroutine CG(m1,m2,m3,m4,m5,m6,m7,ff,ss)
    use CMFD_PARA, only: relax, n_inner, cm_eng, &
                        fm_dhat, dtild, fm_remv, fm0, fm1, fm_phi2
    implicit none
    real(8), intent(in   ):: m1(fm0(1):,fm0(2):,fm0(3):,:), &
                             m2(fm0(1):,fm0(2):,fm0(3):,:), &
                             m3(fm0(1):,fm0(2):,fm0(3):,:), &
                             m4(fm0(1):,fm0(2):,fm0(3):,:), &
                             m5(fm0(1):,fm0(2):,fm0(3):,:), &
                             m6(fm0(1):,fm0(2):,fm0(3):,:), &
                             m7(fm0(1):,fm0(2):,fm0(3):,:)
    real(8), intent(in   ):: ss(fm0(1):,fm0(2):,fm0(3):,:)
    real(8), intent(inout):: ff(fm0(1):,fm0(2):,fm0(3):,:)
    integer:: ii, jj, kk, ee, mm
    real(8):: temp
    integer:: id(3)

    do mm=1, n_inner
    do ee=1, cm_eng
    do kk=fm0(3), fm1(3)
    do jj=fm0(2), fm1(2)
    do ii=fm0(1), fm1(1)
       temp = ss(ii,jj,kk,ee)
       if ( ii /= fm0(1) ) temp = temp - m3(ii,jj,kk,ee)*ff(ii-1,jj,kk,ee)
       if ( ii /= fm1(1) ) temp = temp - m5(ii,jj,kk,ee)*ff(ii+1,jj,kk,ee)
       if ( jj /= fm0(2) ) temp = temp - m2(ii,jj,kk,ee)*ff(ii,jj-1,kk,ee)
       if ( jj /= fm1(2) ) temp = temp - m6(ii,jj,kk,ee)*ff(ii,jj+1,kk,ee)
       if ( kk /= fm0(3) ) temp = temp - m1(ii,jj,kk,ee)*ff(ii,jj,kk-1,ee)
       if ( kk /= fm1(3) ) temp = temp - m7(ii,jj,kk,ee)*ff(ii,jj,kk+1,ee)
       ff(ii,jj,kk,ee) = (1D0-relax)*ff(ii,jj,kk,ee)+relax*temp/m4(ii,jj,kk,ee)
!       if ( isnan(ff(ii,jj,kk,ee)) ) then
!           id(1)=ii;id(2)=jj;id(3)=kk
!           print*, ii, jj, kk
!           print*, ff(ii,jj,kk,ee)
!           print*, ss(ii,jj,kk,ee)
!           print*, temp
!           print*, "1", m1(ii,jj,kk,ee), dtild(ii,jj,kk,ee,5), fm_dhat(ii,jj,kk,ee,5)
!           print*, "2", m2(ii,jj,kk,ee), dtild(ii,jj,kk,ee,3), fm_dhat(ii,jj,kk,ee,3)
!           print*, "3", m3(ii,jj,kk,ee), dtild(ii,jj,kk,ee,1), fm_dhat(ii,jj,kk,ee,1)
!           print*, "4", m4(ii,jj,kk,ee)
!           print*, "5", m5(ii,jj,kk,ee), dtild(ii,jj,kk,ee,2), fm_dhat(ii,jj,kk,ee,2)
!           print*, "6", m6(ii,jj,kk,ee), dtild(ii,jj,kk,ee,4), fm_dhat(ii,jj,kk,ee,4)
!           print*, "7", m7(ii,jj,kk,ee), dtild(ii,jj,kk,ee,6), fm_dhat(ii,jj,kk,ee,6)
!           print*, ff(ii-1,jj,kk,ee)
!           print*, ff(ii+1,jj,kk,ee)
!           print*, ff(ii,jj-1,kk,ee)
!           print*, ff(ii,jj+1,kk,ee)
!           print*, ff(ii,jj,kk-1,ee)
!           print*, ff(ii,jj,kk+1,ee)
!           pause
!       end if
    end do
    end do
    end do
    end do
    end do

end subroutine

! =============================================================================
! CG1 for global calculation
! =============================================================================
subroutine CG1(m1,m2,m3,m4,m5,m6,m7,ff,ss)
    use CMFD_PARA, only: relax, n_inner, cm_eng, cm, cm_dhat, cm_dtild
    implicit none
    real(8), intent(in   ):: m1(1:,1:,1:,1:), &
                             m2(1:,1:,1:,1:), &
                             m3(1:,1:,1:,1:), &
                             m4(1:,1:,1:,1:), &
                             m5(1:,1:,1:,1:), &
                             m6(1:,1:,1:,1:), &
                             m7(1:,1:,1:,1:)
    real(8), intent(in   ):: ss(1:,1:,1:,1:)
    real(8), intent(inout):: ff(1:,1:,1:,1:)
    integer:: ii, jj, kk, ee, mm
    real(8):: temp
    integer:: id(3)

    do mm=1, n_inner
    do ee=1, cm_eng
    do kk=1, cm(3)
    do jj=1, cm(2)
    do ii=1, cm(1)
       temp = ss(ii,jj,kk,ee)
       if ( ii /= 1 )     temp = temp - m3(ii,jj,kk,ee)*ff(ii-1,jj,kk,ee)
       if ( ii /= cm(1) ) temp = temp - m5(ii,jj,kk,ee)*ff(ii+1,jj,kk,ee)
       if ( jj /= 1 )     temp = temp - m2(ii,jj,kk,ee)*ff(ii,jj-1,kk,ee)
       if ( jj /= cm(2) ) temp = temp - m6(ii,jj,kk,ee)*ff(ii,jj+1,kk,ee)
       if ( kk /= 1 )     temp = temp - m1(ii,jj,kk,ee)*ff(ii,jj,kk-1,ee)
       if ( kk /= cm(3) ) temp = temp - m7(ii,jj,kk,ee)*ff(ii,jj,kk+1,ee)
       ff(ii,jj,kk,ee) = (1D0-relax)*ff(ii,jj,kk,ee)+relax*temp/m4(ii,jj,kk,ee)
!       if ( ff(ii,jj,kk,ee) > huge(0Q0) ) then
!           print*, "INF"
!           print*, ii, jj, kk
!           print*, temp, m4(ii,jj,kk,ee)
!           print*, "1", m1(ii,jj,kk,ee), ff(ii,jj,kk-1,ee)
!           print*, "2", m2(ii,jj,kk,ee), ff(ii,jj-1,kk,ee)
!           print*, "3", m3(ii,jj,kk,ee), ff(ii-1,jj,kk,ee)
!           print*, "4", m4(ii,jj,kk,ee), ff(ii,jj,kk,ee)
!           print*, "5", m5(ii,jj,kk,ee), ff(ii+1,jj,kk,ee)
!           print*, "6", m6(ii,jj,kk,ee), ff(ii,jj+1,kk,ee)
!           print*, "7", m7(ii,jj,kk,ee), ff(ii,jj,kk+1,ee)
!           pause
!       end if
!       if ( isnan(ff(ii,jj,kk,ee)) ) then
!           print*, "GLOBAL"
!           id(1)=ii; id(2)=jj; id(3)=kk
!           print*, ii, jj, kk
!           print*, ff(ii,jj,kk,ee)
!           print*, ss(ii,jj,kk,ee)
!           print*, temp
!           print*, "1", m1(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,5), cm_dhat(ii,jj,kk,ee,5)
!           print*, "2", m2(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,3), cm_dhat(ii,jj,kk,ee,3)
!           print*, "3", m3(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,1), cm_dhat(ii,jj,kk,ee,1)
!           print*, "4", m4(ii,jj,kk,ee)                        
!           print*, "5", m5(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,2), cm_dhat(ii,jj,kk,ee,2)
!           print*, "6", m6(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,4), cm_dhat(ii,jj,kk,ee,4)
!           print*, "7", m7(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,6), cm_dhat(ii,jj,kk,ee,6)
!           print*, ff(ii-1,jj,kk,ee)
!           print*, ff(ii+1,jj,kk,ee)
!           print*, ff(ii,jj-1,kk,ee)
!           print*, ff(ii,jj+1,kk,ee)
!           print*, ff(ii,jj,kk-1,ee)
!           print*, ff(ii,jj,kk+1,ee)
!           pause
!       end if
    end do
    end do
    end do
    end do
    end do

end subroutine

! =============================================================================
! CG2 for local calculation
! =============================================================================
subroutine CG2(m1,m2,m3,m4,m5,m6,m7,ff,ss,iii,jjj,kkk)
    use CMFD_PARA, only: relax, n_inner, cm_eng, cm, &
                        fm_dhat, dtild, fm_remv, fm0, fm1, fm_phi2, fcr, fcz
    use GEOMETRY, only: SFM
    implicit none
    real(8), intent(in   ):: m1(1:,1:,1:,1:), &
                             m2(1:,1:,1:,1:), &
                             m3(1:,1:,1:,1:), &
                             m4(1:,1:,1:,1:), &
                             m5(1:,1:,1:,1:), &
                             m6(1:,1:,1:,1:), &
                             m7(1:,1:,1:,1:)
    real(8), intent(in   ):: ss(1:,1:,1:,1:)
    real(8), intent(inout):: ff(1:,1:,1:,1:)
    integer, intent(in)   :: iii, jjj, kkk
    integer:: ii, jj, kk, ee, mm
    real(8):: temp
    integer:: id(3)

    do mm = 1, n_inner
    do ee = 1, cm_eng
    do ii = 1, fcr
    do jj = 1, fcr
    do kk = 1, fcz
       temp = ss(ii,jj,kk,ee)
       if ( ii /= 1 )   temp = temp - m3(ii,jj,kk,ee)*ff(ii-1,jj,kk,ee)
       if ( ii /= fcr ) temp = temp - m5(ii,jj,kk,ee)*ff(ii+1,jj,kk,ee)
       if ( jj /= 1 )   temp = temp - m2(ii,jj,kk,ee)*ff(ii,jj-1,kk,ee)
       if ( jj /= fcr ) temp = temp - m6(ii,jj,kk,ee)*ff(ii,jj+1,kk,ee)
       if ( kk /= 1 )   temp = temp - m1(ii,jj,kk,ee)*ff(ii,jj,kk-1,ee)
       if ( kk /= fcz ) temp = temp - m7(ii,jj,kk,ee)*ff(ii,jj,kk+1,ee)
       ff(ii,jj,kk,ee) = (1D0-relax)*ff(ii,jj,kk,ee)+relax*temp/m4(ii,jj,kk,ee)
!       if ( isnan(ff(ii,jj,kk,ee)) ) then
!           print*, "LOCAL"
!           print*, "C", iii, jjj, kkk
!           print*, "F", ii, jj, kk
!           print*, ff(ii,jj,kk,ee)
!           print*, ss(ii,jj,kk,ee)
!           print*, temp
!           print*, "1", m1(ii,jj,kk,ee), dtild(ii,jj,kk,ee,5), fm_dhat(ii,jj,kk,ee,5)
!           print*, "2", m2(ii,jj,kk,ee), dtild(ii,jj,kk,ee,3), fm_dhat(ii,jj,kk,ee,3)
!           print*, "3", m3(ii,jj,kk,ee), dtild(ii,jj,kk,ee,1), fm_dhat(ii,jj,kk,ee,1)
!           print*, "4", m4(ii,jj,kk,ee)
!           print*, "5", m5(ii,jj,kk,ee), dtild(ii,jj,kk,ee,2), fm_dhat(ii,jj,kk,ee,2)
!           print*, "6", m6(ii,jj,kk,ee), dtild(ii,jj,kk,ee,4), fm_dhat(ii,jj,kk,ee,4)
!           print*, "7", m7(ii,jj,kk,ee), dtild(ii,jj,kk,ee,6), fm_dhat(ii,jj,kk,ee,6)
!           print*, ff(ii-1,jj,kk,ee)
!           print*, ff(ii+1,jj,kk,ee)
!           print*, ff(ii,jj-1,kk,ee)
!           print*, ff(ii,jj+1,kk,ee)
!           pause
!       end if
    end do
    end do
    end do
    end do
    end do

end subroutine

! =============================================================================
! CMFD_CG
! =============================================================================
subroutine CMFD_CG(m1,m2,m3,m4,m5,m6,m7,ff,ss)
    use CMFD_PARA, only: relax, n_inner, cm_eng, cm, cm_dtild, cm_dhat, &
                        fm_dhat, dtild, fm_remv, fm0, fm1, fm_phi2, fcr, fcz
    use GEOMETRY, only: SFM
    implicit none
    real(8), intent(in   ):: m1(1:cm(1),1:cm(2),1:cm(3),1:cm_eng), &
                             m2(1:cm(1),1:cm(2),1:cm(3),1:cm_eng), &
                             m3(1:cm(1),1:cm(2),1:cm(3),1:cm_eng), &
                             m4(1:cm(1),1:cm(2),1:cm(3),1:cm_eng), &
                             m5(1:cm(1),1:cm(2),1:cm(3),1:cm_eng), &
                             m6(1:cm(1),1:cm(2),1:cm(3),1:cm_eng), &
                             m7(1:cm(1),1:cm(2),1:cm(3),1:cm_eng)
    real(8), intent(in   ):: ss(1:cm(1),1:cm(2),1:cm(3),1:cm_eng)
    real(8), intent(inout):: ff(1:cm(1),1:cm(2),1:cm(3),1:cm_eng)
    integer:: ii, jj, kk, ee, mm
    real(8):: temp
    integer:: id(3)

    do mm = 1, n_inner
    do ee = 1, cm_eng
    do ii = 1, cm(1)
    do jj = 1, cm(2)
    do kk = 1, cm(3)
       temp = ss(ii,jj,kk,ee)
       if ( ii /= 1     ) temp = temp - m3(ii,jj,kk,ee)*ff(ii-1,jj,kk,ee)
       if ( ii /= cm(1) ) temp = temp - m5(ii,jj,kk,ee)*ff(ii+1,jj,kk,ee)
       if ( jj /= 1     ) temp = temp - m2(ii,jj,kk,ee)*ff(ii,jj-1,kk,ee)
       if ( jj /= cm(2) ) temp = temp - m6(ii,jj,kk,ee)*ff(ii,jj+1,kk,ee)
       if ( kk /= 1     ) temp = temp - m1(ii,jj,kk,ee)*ff(ii,jj,kk-1,ee)
       if ( kk /= cm(3) ) temp = temp - m7(ii,jj,kk,ee)*ff(ii,jj,kk+1,ee)
       ff(ii,jj,kk,ee) = (1D0-relax)*ff(ii,jj,kk,ee)+relax*temp/m4(ii,jj,kk,ee)
!       if ( isnan(ff(ii,jj,kk,ee)) ) then
!           print*, "C", ii, jj, kk
!           print*, ff(ii,jj,kk,ee)
!           print*, ss(ii,jj,kk,ee)
!           print*, temp
!           print*, "1", m1(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,5), cm_dhat(ii,jj,kk,ee,5)
!           print*, "2", m2(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,3), cm_dhat(ii,jj,kk,ee,3)
!           print*, "3", m3(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,1), cm_dhat(ii,jj,kk,ee,1)
!           print*, "4", m4(ii,jj,kk,ee)
!           print*, "5", m5(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,2), cm_dhat(ii,jj,kk,ee,2)
!           print*, "6", m6(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,4), cm_dhat(ii,jj,kk,ee,4)
!           print*, "7", m7(ii,jj,kk,ee), cm_dtild(ii,jj,kk,ee,6), cm_dhat(ii,jj,kk,ee,6)
!           print*, ff(ii-1,jj,kk,ee)
!           print*, ff(ii+1,jj,kk,ee)
!           print*, ff(ii,jj-1,kk,ee)
!           print*, ff(ii,jj+1,kk,ee)
!           pause
!       end if
    end do
    end do
    end do
    end do
    end do

end subroutine

end module
