module CMFD_TALLY
    use CMFD_PARA
    use GEOMETRY, only: CCM, SI
    implicit none
    integer:: ii, jj, kk, ee, mm, nn, oo
    integer:: id(3)
    
    contains
    
! =============================================================================
! CMFD_J calculates the partial current at the coarse mesh surfaces
! =============================================================================
subroutine CMFD_J(wgt,jj)
    implicit none
    real(8), intent(in   ):: wgt
    real(8), intent(inout):: jj

    jj = jj + wgt

end subroutine

! =============================================================================
! CMFD_SFLUX calculates the reference surface flux
! =============================================================================
subroutine CMFD_SFLUX(uvw,wgt,flux)
    implicit none
    real(8), intent(in   ):: uvw
    real(8), intent(in   ):: wgt
    real(8), intent(inout):: flux

    flux = flux + wgt/abs(uvw)

end subroutine

! =============================================================================
! CMFD_GC_TRK calculates the group constants for FDM
! =============================================================================    
subroutine CMFD_GC_TRK(iter,eng,wgt,mat,dist)
    use MATERIAL_PARA, only: xs, n_eng
    implicit none
    integer, intent(in):: iter
    integer, intent(in):: eng
    real(8), intent(in):: wgt
    integer, intent(in):: mat
    real(8), intent(in):: dist
    logical:: outside   ! is the particle outside the CMFD domain?

    ! CMFD cycle
    if ( iter < cycut ) return
    ! coordinate
    call CM_DOMAIN(fmijk,outside); if ( outside ) return

    ! flux
    fm_phi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_phi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*dist

    ! transport
    fm_tot(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_tot(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*dist*xs(mat)%tot(eng)
    
    ! chi
    fm_chi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_chi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*dist*xs(mat)%chi(eng)
    
    ! fission
    fm_fiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_fiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*dist*xs(mat)%fiss(eng)
    
    ! nufission
    fm_nufiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_nufiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*dist*xs(mat)%nufiss(eng)
    
    ! scattering
    do ee=1, n_eng
    fm_scat(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),cmE(ee)) = &
    fm_scat(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),cmE(ee)) + &
    wgt*dist*xs(mat)%scat(eng,ee)
    end do
    
end subroutine

! =============================================================================
! CMFD_GC_COL calculates the group constants for FDM
! =============================================================================    
subroutine CMFD_GC_COL(iter,eng,wgt,mat)
    use MATERIAL_PARA, only: xs, n_eng
    implicit none
    integer, intent(in):: iter
    integer, intent(in):: eng
    real(8), intent(in):: wgt
    integer, intent(in):: mat
    logical:: outside    ! is the particle outside the CMFD domain?

    ! CMFD cycle
    if ( iter < cycut ) return
    ! coordinate
    call CM_DOMAIN(fmijk(:),outside); if ( outside ) return

    ! flux
    fm_phi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_phi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt/xs(mat)%tot(eng)

    ! transport
    fm_tot(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_tot(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt
    
    ! chi
    fm_chi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_chi(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*xs(mat)%chi(eng)/xs(mat)%tot(eng)
    
    ! fission
    fm_fiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_fiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*xs(mat)%fiss(eng)/xs(mat)%tot(eng)
    
    ! nufission
    fm_nufiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) = &
    fm_nufiss(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng)) + &
    wgt*xs(mat)%nufiss(eng)/xs(mat)%tot(eng)
    
    ! scattering
    do ee=1, n_eng
    fm_scat(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),cmE(ee)) = &
    fm_scat(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),cmE(ee)) + &
    wgt*xs(mat)%scat(eng,ee)/xs(mat)%tot(eng)
    end do
    
end subroutine

! =============================================================================
!
! =============================================================================
subroutine CM_DOMAIN(fmijk,outside)
    integer, intent(in) :: fmijk(:)
    logical, intent(out):: outside

    !!! standard !!!
    outside = .false.
    do ii=1, 3
    if ( fmijk(ii) < fm0(ii) .or. fm1(ii) < fmijk(ii) ) then
        outside = .true.
        exit
    end if
    end do

!    !!! zigzag !!!
!    if ( zigzag ) then
!        ! quarter core
!        if ( mp(1) == 0 .and. mp(2) == 0 ) then
!        if ( fmijk(1) > afm(1)-zz(2) .and. fmijk(2) > afm(2)-zz(2) ) then
!        if ( fmijk(1) > afm(1)-zz(1) .or.  fmijk(2) > afm(2)-zz(1) ) then
!            outside = .true.
!        end if
!        end if
!        
!        ! whole core
!        else
!        id(1:2) = abs(nint(fmijk(1:2)-mp(1:2)))
!        if ( id(1) > afm(1)/2-zz(2) .and. id(2) > afm(2)/2-zz(2) ) then
!        if ( id(1) > afm(1)/2-zz(1) .or.  id(2) > afm(2)/2-zz(1) ) then
!            outside = .true.
!        end if
!        end if
!        end if
!    end if


end subroutine

! =============================================================================
! CMSURFTALLY calculates the surface quantity in the coarse mesh grid
! =============================================================================    
subroutine CMSURFTALLY(iter,eng,uvw,wgt)
    use GEOMETRY_PARA, only: pcross, ce1, ce2, su, ijk
    use GEOMETRY, only: CCM, CI
    implicit none
    integer, intent(in):: iter
    integer, intent(in):: eng
    real(8), intent(in):: uvw(:)
    real(8), intent(in):: wgt
    logical:: outside
    
    if ( .not. pcross ) return
    if ( iter < cycut ) return
    call CM_SURF(su,outside); if ( outside ) return


    select case(su)
    case(1)
    if ( fmijk(1) == fm1(1)+1 ) then
    call CMFD_J(wgt,fmJ0(iter,fmijk(1)-1,fmijk(2),fmijk(3),cmE(eng),2))
    call CMFD_SFLUX(uvw(1),wgt,fmF(iter,fmijk(1)-1,fmijk(2),fmijk(3),cmE(eng),2))
    else
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(1),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    case(2)
    if ( fmijk(1) == fm0(1)-1 ) then
    call CMFD_J(wgt,fmJ1(iter,fmijk(1)+1,fmijk(2),fmijk(3),cmE(eng),1))
    call CMFD_SFLUX(uvw(1),wgt,fmF(iter,fmijk(1)+1,fmijk(2),fmijk(3),cmE(eng),1))
    else
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(1),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    case(3)
    if ( fmijk(2) == fm1(2)+1 ) then
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2)-1,fmijk(3),cmE(eng),4))
    call CMFD_SFLUX(uvw(2),wgt,fmF(iter,fmijk(1),fmijk(2)-1,fmijk(3),cmE(eng),4))
    else
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(2),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    case(4)
    if ( fmijk(2) == fm0(2)-1 ) then
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2)+1,fmijk(3),cmE(eng),3))
    call CMFD_SFLUX(uvw(2),wgt,fmF(iter,fmijk(1),fmijk(2)+1,fmijk(3),cmE(eng),3))
    else
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(2),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    case(5)
    if ( fmijk(3) == fm1(3)+1 ) then
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2),fmijk(3)-1,cmE(eng),6))
    call CMFD_SFLUX(uvw(3),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3)-1,cmE(eng),6))
    else
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(3),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    case(6)
    if ( fmijk(3) == fm0(3)-1 ) then
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2),fmijk(3)+1,cmE(eng),5))
    call CMFD_SFLUX(uvw(3),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3)+1,cmE(eng),5))
    else
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(3),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    end select

end subroutine

! =============================================================================
! CM_SURF determines whether the particle crosses the CM surface
! =============================================================================
subroutine CM_SURF(surf,outside)
    use GEOMETRY_PARA, only: ijk
    implicit none
    integer, intent(in) :: surf
    logical, intent(out):: outside

    !!! standard !!!
    outside = .true. 
    select case(surf)
    case(1)
    if ( fmijk(1)>=fm0(1) .and. fmijk(2)>=fm0(2) .and. fmijk(3)>=fm0(3) ) then
    if ( fmijk(1)<=fm1(1)+1 .and. fmijk(2)<=fm1(2) .and. fmijk(3)<=fm1(3) )then
    if ( mod(ijk(1)-fm0(1),nfm(1)) == 0  ) outside = .false.
    end if
    end if
    case(2)
    if ( fmijk(1)>=fm0(1)-1 .and. fmijk(2)>=fm0(2) .and. fmijk(3)>=fm0(3) ) then
    if ( fmijk(1)<=fm1(1) .and. fmijk(2)<=fm1(2) .and. fmijk(3)<=fm1(3) )then
    if ( mod(ijk(1)-fm0(1)+1,nfm(1)) == 0 ) outside = .false.
    end if
    end if
    case(3)
    if ( fmijk(1)>=fm0(1) .and. fmijk(2)>=fm0(2) .and. fmijk(3)>=fm0(3) ) then
    if ( fmijk(1)<=fm1(1) .and. fmijk(2)<=fm1(2)+1 .and. fmijk(3)<=fm1(3) )then
    if ( mod(ijk(2)-fm0(2),nfm(2)) == 0  ) outside = .false.
    end if
    end if
    case(4)
    if ( fmijk(1)>=fm0(1) .and. fmijk(2)>=fm0(2)-1 .and. fmijk(3)>=fm0(3) ) then
    if ( fmijk(1)<=fm1(1) .and. fmijk(2)<=fm1(2) .and. fmijk(3)<=fm1(3) )then
    if ( mod(ijk(2)-fm0(2)+1,nfm(2)) == 0 ) outside = .false.
    end if
    end if
    case(5)
    if ( fmijk(1)>=fm0(1) .and. fmijk(2)>=fm0(2) .and. fmijk(3)>=fm0(3) ) then
    if ( fmijk(1)<=fm1(1) .and. fmijk(2)<=fm1(2) .and. fmijk(3)<=fm1(3)+1 )then
    if ( mod(ijk(3)-fm0(3),nfm(3)) == 0  ) outside = .false.
    end if
    end if
    case(6)
    if ( fmijk(1)>=fm0(1) .and. fmijk(2)>=fm0(2) .and. fmijk(3)>=fm0(3)-1 ) then
    if ( fmijk(1)<=fm1(1) .and. fmijk(2)<=fm1(2) .and. fmijk(3)<=fm1(3) )then
    if ( mod(ijk(3)-fm0(3)+1,nfm(3)) == 0 ) outside = .false.
    end if
    end if
    end select

!    !!! zigzag !!!
!    if ( zigzag ) then
!        ! quarter core
!        if ( mp(1) == 0 .and. mp(2) == 0 ) then
!        if ( fmijk(1) > afm(1)-zz(2) .and. fmijk(2) > afm(2)-zz(2) ) then
!        if ( fmijk(1) > afm(1)-zz(1) .or.  fmijk(2) > afm(2)-zz(1) ) then
!            outside = .true.
!        end if
!        end if
!
!        select case(surf)
!        case(1)
!        if ( fmijk(1) == afm(1)-zz(2)+1 .and. afm(2)-zz(1) < fmijk(2) .and. &
!             fmijk(2) <= afm(2) ) outside = .false.
!        if ( fmijk(1) == afm(1)-zz(1)+1 .and. afm(2)-zz(2) < fmijk(2) .and. &
!             fmijk(2) <= afm(2)-zz(1) ) outside = .false.
!
!        case(3)
!        if ( fmijk(2) == afm(2)-zz(2)+1 .and. afm(2)-zz(1) < fmijk(1) .and. &
!             fmijk(1) <= afm(1) ) outside = .false.
!        if ( fmijk(2) == afm(2)-zz(1)+1 .and. afm(2)-zz(2) < fmijk(1) .and. &
!             fmijk(1) <= afm(1)-zz(1) ) outside = .false.
!
!        end select
!    
!        ! whole core
!        else
!        id(1:2) = abs(nint(fmijk(1:2)-mp(1:2)))
!        if ( id(1) > afm(1)/2-zz(2) .and. id(2) > afm(2)/2-zz(2) ) then
!        if ( id(1) > afm(1)/2-zz(1) .or.  id(2) > afm(2)/2-zz(1) ) then
!            outside = .true.
!        end if
!        end if
!    
!        print*, "see 16_CMFD"
!        print*, "not developed yet"
!        stop
!    
!        end if
!    end if

end subroutine


! =============================================================================
! CMBOUNDTALLY calculates the group constants for FDM
! =============================================================================    
subroutine CMBOUNDTALLY(iter,eng,uvw,wgt)
    use GEOMETRY_PARA, only: bc_x0, bc_x1, bc_y0, bc_y1, bc_z0, bc_z1, &
                           pcross, ce1, ce2, su, ijk, bound, &
                           xcore, ycore, zcore, n_pin
    implicit none
    integer, intent(in):: iter
    integer, intent(in):: eng
    real(8), intent(in):: uvw(:)
    real(8), intent(in):: wgt

    if ( .not. bound  ) return
    if ( iter < cycut ) return

    select case(su)
    case(1)
    if ( bc_x0 == -1 .and. fm0(1) == 1 ) then
    if ( fm0(2) <= fmijk(2) .and. fmijk(2) <= fm1(2) ) then
    if ( fm0(3) <= fmijk(3) .and. fmijk(3) <= fm1(3) ) then
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(1),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    end if
    end if

    case(2)
    if ( bc_x1 == -1 .and. fm0(1) == xcore*n_pin ) then
    if ( fm0(2) <= fmijk(2) .and. fmijk(2) <= fm1(2) ) then
    if ( fm0(3) <= fmijk(3) .and. fmijk(3) <= fm1(3) ) then
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(1),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    end if
    end if

    case(3)
    if ( bc_y0 == -1 .and. fm0(2) == 1 ) then
    if ( fm0(1) <= fmijk(1) .and. fmijk(1) <= fm1(1) ) then
    if ( fm0(3) <= fmijk(3) .and. fmijk(3) <= fm1(3) ) then
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(2),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    end if
    end if
    
    case(4)
    if ( bc_y1 == -1 .and. fm0(2) == ycore*n_pin ) then
    if ( fm0(1) <= fmijk(1) .and. fmijk(1) <= fm1(1) ) then
    if ( fm0(3) <= fmijk(3) .and. fmijk(3) <= fm1(3) ) then
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(2),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    end if
    end if
    
    case(5)
    if ( bc_z0 == -1 .and. fm0(3) == 1 ) then
    if ( fm0(1) <= fmijk(1) .and. fmijk(1) <= fm1(1) ) then
    if ( fm0(2) <= fmijk(2) .and. fmijk(2) <= fm1(2) ) then
    call CMFD_J(wgt,fmJ1(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(3),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    end if
    end if
    
    case(6)
    if ( bc_z1 == -1 .and. fm0(3) == zcore ) then
    if ( fm0(1) <= fmijk(1) .and. fmijk(1) <= fm1(1) ) then
    if ( fm0(2) <= fmijk(2) .and. fmijk(2) <= fm1(2) ) then
    call CMFD_J(wgt,fmJ0(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    call CMFD_SFLUX(uvw(3),wgt,fmF(iter,fmijk(1),fmijk(2),fmijk(3),cmE(eng),su))
    end if
    end if
    end if
    end select

end subroutine


! =============================================================================
! CMFD_NORMAL
! =============================================================================
subroutine FMFD_NORMAL(i)
    use GLOBAL_PARA, only: n_hist
    use GEOMETRY_PARA, only: xcore, ycore, zcore
    use STATISTICS, only: AVG
    use GEOMETRY_PARA, only: pitch, height, n_pin
    implicit none
    integer, intent(in):: i     ! cycle index
    integer:: i1, i2            ! cycle index
    real(8):: surf
    real(8):: vol
    real(8):: summ0, summ1
    integer:: mm

    call SURFCOMB(i)
    call L_ACCUM(i,i1,i2)

    ! -------------------------------------------------------------------------
    ! cycle average
    do ii = fm0(1), fm1(1)
    do jj = fm0(2), fm1(2)
    do kk = fm0(3), fm1(3)

    !   volume quantity
    if ( sum(fm_chi(i,ii,jj,kk,:)) /= 0 ) &
    fm_chi(i,ii,jj,kk,:) = fm_chi(i,ii,jj,kk,:)/sum(fm_chi(i,ii,jj,kk,:))

    do ee=1, cm_eng
        fm_tot(0,ii,jj,kk,ee)     = AVG(fm_tot(i1:i2,ii,jj,kk,ee))
        fm_fiss(0,ii,jj,kk,ee)    = AVG(fm_fiss(i1:i2,ii,jj,kk,ee))
        fm_chi(0,ii,jj,kk,ee)     = AVG(fm_chi(i1:i2,ii,jj,kk,ee))
        fm_nufiss(0,ii,jj,kk,ee)  = AVG(fm_nufiss(i1:i2,ii,jj,kk,ee))
        fm_phi(0,ii,jj,kk,ee)     = AVG(fm_phi(i1:i2,ii,jj,kk,ee))
    do mm=1, cm_eng
        fm_scat(0,ii,jj,kk,ee,mm) = AVG(fm_scat(i1:i2,ii,jj,kk,ee,mm))
    end do
        fm_remv(ii,jj,kk,ee)      = &
        fm_tot(0,ii,jj,kk,ee)     - fm_scat(0,ii,jj,kk,ee,ee)

    ! surface average quantity
    do mm = 1, 6
        fmF(0,ii,jj,kk,ee,mm)  = AVG(fmF(i1:i2,ii,jj,kk,ee,mm))
        fmJ0(0,ii,jj,kk,ee,mm) = AVG(fmJ0(i1:i2,ii,jj,kk,ee,mm))
        fmJ1(0,ii,jj,kk,ee,mm) = AVG(fmJ1(i1:i2,ii,jj,kk,ee,mm))
    end do
    end do

    end do
    end do
    end do
    if ( cmfdon .or. ccmfd ) call LTOG1

    ! -------------------------------------------------------------------------
    ! normalization
    !   volume quantity
    fm_tot(0,:,:,:,:)     = fm_tot(0,:,:,:,:)/fm_phi(0,:,:,:,:)
    fm_fiss(0,:,:,:,:)    = fm_fiss(0,:,:,:,:)/fm_phi(0,:,:,:,:)
    fm_nufiss(0,:,:,:,:)  = fm_nufiss(0,:,:,:,:)/fm_phi(0,:,:,:,:)
    fm_remv(:,:,:,:)      = fm_remv(:,:,:,:)/fm_phi(0,:,:,:,:)
    do ee=1, cm_eng
    fm_scat(0,:,:,:,:,ee) = fm_scat(0,:,:,:,:,ee)/fm_phi(0,:,:,:,:)
    end do
    do kk=fm0(3), fm1(3); vol  = pfm*pfm*hfm(kk)
    fm_phi(0,:,:,kk,:)    = fm_phi(0,:,:,kk,:)/(vol*n_hist)
    end do
    if ( unbiased ) call XS_UNBIASED(i1,i2)
    !   surface quantity
    do kk = fm0(3), fm1(3); surf = pfm*hfm(kk)
    fmF(0,:,:,kk,:,1:4)  = fmF(0,:,:,kk,:,1:4) / (surf*n_hist)
    fmJ0(0,:,:,kk,:,1:4) = fmJ0(0,:,:,kk,:,1:4) / (surf*n_hist)
    fmJ1(0,:,:,kk,:,1:4) = fmJ1(0,:,:,kk,:,1:4) / (surf*n_hist)
    end do
    surf = pfm*pfm
    fmF(0,:,:,:,:,5:6)  = fmF(0,:,:,:,:,5:6) / (surf*n_hist)
    fmJ0(0,:,:,:,:,5:6) = fmJ0(0,:,:,:,:,5:6) / (surf*n_hist)
    fmJ1(0,:,:,:,:,5:6) = fmJ1(0,:,:,:,:,5:6) / (surf*n_hist)
    fmJn(:,:,:,:,:) = fmJ1(0,:,:,:,:,:) - fmJ0(0,:,:,:,:,:)
    if ( cmfdon .or. ccmfd ) call LTOG2

    where ( fm_phi(0,:,:,:,:) == 0 ) fm_tot(0,:,:,:,:)     = 0
    where ( fm_phi(0,:,:,:,:) == 0 ) fm_fiss(0,:,:,:,:)    = 0
    where ( fm_phi(0,:,:,:,:) == 0 ) fm_nufiss(0,:,:,:,:)  = 0
    do ee=1, cm_eng
    where ( fm_phi(0,:,:,:,:) == 0 ) fm_scat(0,:,:,:,:,ee) = 0
    end do

end subroutine


! =============================================================================
! L_ACCUM
! =============================================================================
subroutine L_ACCUM(i,i1,i2)
    integer, intent(in):: i         ! cycle number
    integer, intent(out):: i1, i2   ! initial cycle, closing cycle

    if (      i <  cycut+accut ) then;        i1 = i;           i2 = i
    else if ( i >= cycut+accut .and. &
              i <  cycut+accut+accum ) then;  i1 = cycut+accut; i2 = i
    else if ( i >= cycut+accut+accum ) then;  i1 = i-accum+1;   i2 = i
    end if

end subroutine

! =============================================================================
! XS_UNBIASED
! =============================================================================
subroutine XS_UNBIASED(i1,i2)
    use STATISTICS, only: BIAS
    use GLOBAL_PARA, only: filename
    implicit none
    integer, intent(in):: i1, i2

    do ee = 1, cm_eng
    do ii = fm0(1), fm1(1)
    do jj = fm0(2), fm1(2)
    do kk = fm0(3), fm1(3)
        if ( fm_phi(0,ii,jj,kk,ee) /= 0 ) then
        fm_tot(0,ii,jj,kk,ee) = fm_tot(0,ii,jj,kk,ee) &
            -BIAS(fm_tot(i1:i2,ii,jj,kk,ee),fm_phi(i1:i2,ii,jj,kk,ee))
        fm_fiss(0,ii,jj,kk,ee) = fm_fiss(0,ii,jj,kk,ee) &
            -BIAS(fm_fiss(i1:i2,ii,jj,kk,ee),fm_phi(i1:i2,ii,jj,kk,ee))
        fm_nufiss(0,ii,jj,kk,ee) = fm_nufiss(0,ii,jj,kk,ee) &
            -BIAS(fm_nufiss(i1:i2,ii,jj,kk,ee),fm_phi(i1:i2,ii,jj,kk,ee))
        do mm = 1, cm_eng
        fm_scat(0,ii,jj,kk,ee,mm) = fm_scat(0,ii,jj,kk,ee,mm) &
            -BIAS(fm_scat(i1:i2,ii,jj,kk,ee,mm),fm_phi(i1:i2,ii,jj,kk,ee))
        end do
        fm_remv(ii,jj,kk,ee) = fm_tot(0,ii,jj,kk,ee)-fm_scat(0,ii,jj,kk,ee,ee)
        end if
    end do
    end do
    end do
    end do

end subroutine

! =============================================================================
! SURFCOMB
! =============================================================================
subroutine SURFCOMB(iter)
    integer:: iter
    real(8):: ssum(cm_eng)
   
    do ii = fm0(1), fm1(1)
    do jj = fm0(2), fm1(2)
    do kk = fm0(3), fm1(3)
        
        ! x-direction
        if ( ii /= fm1(1) ) then
        ssum = fmJ0(iter,ii+1,jj,kk,:,1)+fmJ0(iter,ii,jj,kk,:,2)
        fmJ0(iter,ii+1,jj,kk,:,1) = ssum
        fmJ0(iter,ii,jj,kk,:,2)   = ssum
        ssum = fmJ1(iter,ii+1,jj,kk,:,1)+fmJ1(iter,ii,jj,kk,:,2)
        fmJ1(iter,ii+1,jj,kk,:,1) = ssum
        fmJ1(iter,ii,jj,kk,:,2)   = ssum
        ssum = fmF(iter,ii+1,jj,kk,:,1)+fmF(iter,ii,jj,kk,:,2)
        fmF(iter,ii+1,jj,kk,:,1)  = ssum
        fmF(iter,ii,jj,kk,:,2)    = ssum
        end if

        ! y-direction
        if ( jj /= fm1(2) ) then
        ssum = fmJ0(iter,ii,jj+1,kk,:,3)+fmJ0(iter,ii,jj,kk,:,4)
        fmJ0(iter,ii,jj+1,kk,:,3) = ssum
        fmJ0(iter,ii,jj,kk,:,4)   = ssum
        ssum = fmJ1(iter,ii,jj+1,kk,:,3)+fmJ1(iter,ii,jj,kk,:,4)
        fmJ1(iter,ii,jj+1,kk,:,3) = ssum
        fmJ1(iter,ii,jj,kk,:,4)   = ssum
        ssum = fmF(iter,ii,jj+1,kk,:,3)+fmF(iter,ii,jj,kk,:,4)
        fmF(iter,ii,jj+1,kk,:,3)  = ssum
        fmF(iter,ii,jj,kk,:,4)    = ssum
        end if

        ! z-direction
        if ( kk /= fm1(3) ) then
        ssum = fmJ0(iter,ii,jj,kk+1,:,5)+fmJ0(iter,ii,jj,kk,:,6)
        fmJ0(iter,ii,jj,kk+1,:,5) = ssum
        fmJ0(iter,ii,jj,kk,:,6)   = ssum
        ssum = fmJ1(iter,ii,jj,kk+1,:,5)+fmJ1(iter,ii,jj,kk,:,6)
        fmJ1(iter,ii,jj,kk+1,:,5) = ssum
        fmJ1(iter,ii,jj,kk,:,6)   = ssum
        ssum = fmF(iter,ii,jj,kk+1,:,5)+fmF(iter,ii,jj,kk,:,6)
        fmF(iter,ii,jj,kk+1,:,5)  = ssum
        fmF(iter,ii,jj,kk,:,6)    = ssum
        end if

    end do
    end do
    end do

end subroutine


! =============================================================================
! FMFD_PARAMTER
! =============================================================================
subroutine FMFD_PARAMETER
    use GEOMETRY_PARA, only: n_pin
    use STATISTICS, only: AVG
    implicit none
    integer:: id0(3)
    real(8):: summ(0:2)         ! parameters

    fm_phi2(:,:,:,:) = fm_phi(0,:,:,:,:)            ! flux
    fm_diff(:,:,:,:) = 1D0/(3D0*fm_tot(0,:,:,:,:))  ! diffusion coefficient

    !!! zigzag !!! zero diffusion coefficient
    if ( zigzag ) where ( fm_tot(0,:,:,:,:) == 0 ) fm_diff(:,:,:,:) = 0D0

    !!! standard !!!
    !   unit diffusion coefficient
    fm_ddi(:,:,:,:)  = fm_diff(:,:,:,:)/pfm
    fm_ddj(:,:,:,:)  = fm_diff(:,:,:,:)/pfm
    do kk=fm0(3), fm1(3)
    fm_ddk(:,:,kk,:) = fm_diff(:,:,kk,:)/hfm(kk)
    end do
    !   interface diffusion coefficient
    do kk = fm0(3), fm1(3)
    do jj = fm0(2), fm1(2)
    do ii = fm0(1), fm1(1)
        ! x0
        if ( ii /= fm0(1) ) then
        dtild(ii,jj,kk,:,1)= &
            2D0*fm_ddi(ii,jj,kk,:)*fm_ddi(ii-1,jj,kk,:) &
            /(fm_ddi(ii,jj,kk,:)+fm_ddi(ii-1,jj,kk,:))
        deltf0(ii,jj,kk,:,1) = dtild(ii,jj,kk,:,1)
        end if
        ! x1
        if ( ii /= fm1(1) ) then
        dtild(ii,jj,kk,:,2)= &
            2D0*fm_ddi(ii+1,jj,kk,:)*fm_ddi(ii,jj,kk,:) &
            /(fm_ddi(ii+1,jj,kk,:)+fm_ddi(ii,jj,kk,:))
        deltf0(ii,jj,kk,:,2) = dtild(ii,jj,kk,:,2)
        end if
        ! y0
        if ( jj /= fm0(2) ) then
        dtild(ii,jj,kk,:,3)= &
            2D0*fm_ddj(ii,jj,kk,:)*fm_ddj(ii,jj-1,kk,:) &
            /(fm_ddj(ii,jj,kk,:)+fm_ddj(ii,jj-1,kk,:))
        deltf0(ii,jj,kk,:,3) = dtild(ii,jj,kk,:,3)
        end if
        ! y1
        if ( jj /= fm1(2) ) then
        dtild(ii,jj,kk,:,4)= &
            2D0*fm_ddj(ii,jj+1,kk,:)*fm_ddj(ii,jj,kk,:) &
            /(fm_ddj(ii,jj+1,kk,:)+fm_ddj(ii,jj,kk,:))
        deltf0(ii,jj,kk,:,4) = dtild(ii,jj,kk,:,4)
        end if
        ! z0
        if ( kk /= fm0(3) ) then
        dtild(ii,jj,kk,:,5)= &
            2D0*fm_ddk(ii,jj,kk,:)*fm_ddk(ii,jj,kk-1,:) &
            /(fm_ddk(ii,jj,kk,:)+fm_ddk(ii,jj,kk-1,:))
        deltf0(ii,jj,kk,:,5) = dtild(ii,jj,kk,:,5)
        end if
        ! z1
        if ( kk /= fm1(3) ) then
        dtild(ii,jj,kk,:,6)= &
            2D0*fm_ddk(ii,jj,kk+1,:)*fm_ddk(ii,jj,kk,:) &
            /(fm_ddk(ii,jj,kk+1,:)+fm_ddk(ii,jj,kk,:))
        deltf0(ii,jj,kk,:,6) = dtild(ii,jj,kk,:,6)
        end if
    end do
    end do
    end do

    !!! CMFD !!!
    if ( cmfdon .or. ccmfd ) then
    !   unit diffusion coefficient
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
    end if

    if ( ccmfd ) then
    cm_dtild = 0D0
    do kk = 1, cm(3)
    do jj = 1, cm(2)
    do ii = 1, cm(1)
        ! x0
        if ( ii /= 1 ) then
        cm_dtild(ii,jj,kk,:,1)= &
            2D0*cm_ddi(ii,jj,kk,:)*cm_ddi(ii-1,jj,kk,:) &
            /(cm_ddi(ii,jj,kk,:)+cm_ddi(ii-1,jj,kk,:))
        end if
        ! x1
        if ( ii /= cm(1) ) then
        cm_dtild(ii,jj,kk,:,2)= &
            2D0*cm_ddi(ii+1,jj,kk,:)*cm_ddi(ii,jj,kk,:) &
            /(cm_ddi(ii+1,jj,kk,:)+cm_ddi(ii,jj,kk,:))
        end if
        ! y0
        if ( jj /= 1 ) then
        cm_dtild(ii,jj,kk,:,3)= &
            2D0*cm_ddj(ii,jj,kk,:)*cm_ddj(ii,jj-1,kk,:) &
            /(cm_ddj(ii,jj,kk,:)+cm_ddj(ii,jj-1,kk,:))
        end if
        ! y1
        if ( jj /= cm(2) ) then
        cm_dtild(ii,jj,kk,:,4)= &
            2D0*cm_ddj(ii,jj+1,kk,:)*cm_ddj(ii,jj,kk,:) &
            /(cm_ddj(ii,jj+1,kk,:)+cm_ddj(ii,jj,kk,:))
        end if
        ! z0
        if ( kk /= 1 ) then
        cm_dtild(ii,jj,kk,:,5)= &
            2D0*cm_ddk(ii,jj,kk,:)*cm_ddk(ii,jj,kk-1,:) &
            /(cm_ddk(ii,jj,kk,:)+cm_ddk(ii,jj,kk-1,:))
        end if
        ! z1
        if ( kk /= cm(3) ) then
        cm_dtild(ii,jj,kk,:,6)= &
            2D0*cm_ddk(ii,jj,kk+1,:)*cm_ddk(ii,jj,kk,:) &
            /(cm_ddk(ii,jj,kk+1,:)+cm_ddk(ii,jj,kk,:))
        end if
    end do
    end do
    end do
    end if

end subroutine

! =============================================================================
! FMFD_RESULT
! =============================================================================
subroutine FMFD_RESULT(i)
    use GLOBAL_PARA, only: n_icycle
    use TALLY_PARA, only: tally5
    implicit none
    integer, intent(in):: i
    integer:: iter

    if ( i <= n_icycle ) return
    if ( .not. tally5 )  return

    iter = i - n_icycle 
    
    ! multiplication factor
    fm_keff(iter) = k_cmfd

    ! power distribution
    fm_power(iter,:,:,:,:) = fm_fiss(0,:,:,:,:)*fm_phi2(:,:,:,:)

end subroutine

! =============================================================================
! cmE determines the energy in terms of CMFD
! =============================================================================
function cmE(eng)
    integer:: cmE
    integer, intent(in):: eng

    cmE = 1
    do ii=cm_eng, 2, -1
    if ( eng >= cmeng(ii) ) then
        cmE = ii
        exit
    end if
    end do

end function

! =============================================================================
! LTOG1 converts node parameters from LOCAL to GLOBAL (average)
! =============================================================================
subroutine LTOG1
    use GLOBAL_PARA, only: n_hist
    implicit none
    integer:: id1(3)
    integer:: id2(3)
    real(8):: vol

    ! homogenization
    do ee = 1, cm_eng
    do ii = 1, cm(1); id2(1) = (ii-1)*fcr+fm0(1)
    do jj = 1, cm(2); id2(2) = (jj-1)*fcr+fm0(2)
    do kk = 1, cm(3); id2(3) = (kk-1)*fcz+fm0(3)
    cm_tot(ii,jj,kk,ee) = sum(fm_tot(0,id2(1):id2(1)+fcr-1, &
                            id2(2):id2(2)+fcr-1,id2(3):id2(3)+fcz-1,ee))
    cm_chi(ii,jj,kk,ee) = sum(fm_chi(0,id2(1):id2(1)+fcr-1, &
                            id2(2):id2(2)+fcr-1,id2(3):id2(3)+fcz-1,ee)* &
                            fm_nufiss(0,id2(1):id2(1)+fcr-1, &
                            id2(2):id2(2)+fcr-1,id2(3):id2(3)+fcz-1,ee))
    cm_fiss(ii,jj,kk,ee) = sum(fm_fiss(0,id2(1):id2(1)+fcr-1, &
                            id2(2):id2(2)+fcr-1,id2(3):id2(3)+fcz-1,ee))
    cm_nufiss(ii,jj,kk,ee) = sum(fm_nufiss(0,id2(1):id2(1)+fcr-1, &
                            id2(2):id2(2)+fcr-1,id2(3):id2(3)+fcz-1,ee))
    cm_phi2(ii,jj,kk,ee) = sum(fm_phi(0,id2(1):id2(1)+fcr-1, &
                            id2(2):id2(2)+fcr-1,id2(3):id2(3)+fcz-1,ee))
    do mm = 1, cm_eng
    cm_scat(ii,jj,kk,ee,mm) = sum(fm_scat(0,id2(1):id2(1)+fcr-1, &
                         id2(2):id2(2)+fcr-1,id2(3):id2(3)+fcz-1,ee,mm))
    end do
    cm_remv(ii,jj,kk,ee) = cm_tot(ii,jj,kk,ee) - cm_scat(ii,jj,kk,ee,ee)
    end do
    end do
    end do
    end do

    ! cross-section
    cm_tot = cm_tot / cm_phi2
    cm_chi = cm_chi / cm_nufiss
    cm_fiss = cm_fiss / cm_phi2
    cm_nufiss = cm_nufiss / cm_phi2
    cm_remv = cm_remv / cm_phi2
    do mm = 1, cm_eng
    cm_scat(:,:,:,:,mm) = cm_scat(:,:,:,:,mm) / cm_phi2(:,:,:,:)
    end do
    do kk = 1, cm(3); id2(3) = (kk-1)*fcz+fm0(3)
    vol = pfm*fcr; vol = vol*vol*sum(hfm(id2(3):id2(3)+fcz-1))
    cm_phi2(:,:,kk,:) = cm_phi2(:,:,kk,:) / (vol*n_hist)
    end do

    ! diffusion coefficient
    cm_diff = 1D0/(3D0*cm_tot)
    where ( cm_tot == 0 ) cm_diff = 0D0

    !!! zigzag !!!
    if ( zigzag ) then
        where ( cm_phi2 == 0 ) cm_tot    = 0
        where ( cm_phi2 == 0 ) cm_nufiss = 0
        where ( cm_phi2 == 0 ) cm_diff   = 0
        where ( cm_phi2 == 0 ) cm_chi    = 0
        do ee = 1, cm_eng
        where ( cm_phi2(:,:,:,:) == 0 ) cm_scat(:,:,:,:,ee) = 0
        end do
    end if

end subroutine

! =============================================================================
! LTOG2 converts surface parameters from LOCAL to GLOBAL (average)
! =============================================================================
subroutine LTOG2
    implicit none
    integer:: id0(3)
    real(8):: ssum(0:2)

    do ee = 1, cm_eng
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
        ! x-direction
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
        if ( ii /= 1 ) then
        cmJ0(ii-1,jj,kk,ee,2) = cmJ0(ii,jj,kk,ee,1)
        cmJ1(ii-1,jj,kk,ee,2) = cmJ1(ii,jj,kk,ee,1)
        cmF(ii-1,jj,kk,ee,2)  = cmF(ii,jj,kk,ee,1)
        end if
        ! y-direction
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
        if ( jj /= 1 ) then
        cmJ0(ii,jj-1,kk,ee,4) = cmJ0(ii,jj,kk,ee,3)
        cmJ1(ii,jj-1,kk,ee,4) = cmJ1(ii,jj,kk,ee,3)
        cmF(ii,jj-1,kk,ee,4)  = cmF(ii,jj,kk,ee,3)
        end if
        ! z-direction
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
        if ( kk /= 1 ) then
        cmJ0(ii,jj,kk-1,ee,6) = cmJ0(ii,jj,kk,ee,5)
        cmJ1(ii,jj,kk-1,ee,6) = cmJ1(ii,jj,kk,ee,5)
        cmF(ii,jj,kk-1,ee,6)  = cmF(ii,jj,kk,ee,5)
        end if
    end do
    end do
    end do
    ! Closure
    !   x-direction
    ii = cm(1);       id(1) = fm1(1)
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1; ssum = 0
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
    end do
    end do
    !   y-direction
    jj = cm(2);       id(2) = fm1(2)
    do kk = 1, cm(3); id0(3) = (kk-1)*fcz+fm0(3)-1
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1; ssum = 0
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
    end do
    end do
    !   z-direction
    kk = cm(3);       id(3) = fm1(3)
    do ii = 1, cm(1); id0(1) = (ii-1)*fcr+fm0(1)-1
    do jj = 1, cm(2); id0(2) = (jj-1)*fcr+fm0(2)-1; ssum = 0
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
    end do
    end do
    end do
    cmJn = cmJ1 - cmJ0

end subroutine

end module



