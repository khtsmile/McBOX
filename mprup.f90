module MPRUP
    use VARIABLES, only: icore, score, ngen, n_totcyc, n_inact, n_act
    use ENTROPY
    use FMFD_HEADER, only: fmfdon, k_fmfd
    implicit none

    ! * Convergence criteria
    !   1 convergence : crt1
    !    (1) relative error
    !    (2) length of cycles
    !   2 convergence : crt2
    !    (1) relative error
    !    (2) final generation size

    contains

! =============================================================================
! 
! =============================================================================
subroutine GENSIZE(cyc)
    use FMFD_HEADER, only: acc, n_acc
    implicit none
    integer, intent(in) :: cyc
    integer:: ii, jj 

    ! FMFD divergence
    if ( cyc > 1 .and. fmfdon ) then
    if ( isnan(k_fmfd(cyc-1)) .or. &
        ( k_fmfd(cyc-1) < 0D0 .or. k_fmfd(cyc-1) > 2D0 ) ) then
        print*, "1"
        do ii = 2, n_acc
        acc(ii)%fm(:,:,:)%phi     = 0
        acc(ii)%fm(:,:,:)%sig_t   = 0
        acc(ii)%fm(:,:,:)%sig_a   = 0
        acc(ii)%fm(:,:,:)%nusig_f = 0
        do jj = 1, 6
        acc(ii)%fm(:,:,:)%Jn(jj)  = 0
        acc(ii)%fm(:,:,:)%J0(jj)  = 0
        acc(ii)%fm(:,:,:)%J1(jj)  = 0
        end do
        end do
        call GENSIZEUP
        return
    end if
    end if

    if ( genup ) then
    ! 1 convergence test
    !print*, "dsh", crt1, dshannon
    select case(ccrt)
    case(1); if ( dshannon > crt1 ) return
    case(2); if ( cyc == 1 .or. mod(cyc-1,int(crt1)) /= 0 ) return
    end select

    ! 2 stationary point
    select case(scrt)
    case(1)
        entrp2 = sum(shannon(1:elength))/dble(elength)
        dentrp = abs(entrp2-entrp1)/entrp2
        !print*, "den", crt2, dentrp, "*"
        if ( dentrp < crt2 ) then
            genup = .false.
            return
        end if
        entrp1 = entrp2

        call GENSIZEUP

    case(2)
        if ( ngen >= crt2 ) then
            genup = .false.
            return
        end if
        call GENSIZEUP

    end select

    ! 3 stopping test
    else
        if ( dshannon < crt3 ) then
            up_sign = .true.
            mprupon = .false.
            genup = .true.
            !call CYCLECHANGE(cyc)
            return
        end if

    end if


end subroutine


! =============================================================================
! GENSIZEUP increases the generation size and reallocate the fission neutron
! =============================================================================
subroutine GENSIZEUP
    use BANK_HEADER,    only: source_bank
    use FMFD_HEADER,    only: fmfdon
    implicit none

!    source_bank(:)%wgt  = ngen/size(source_bank)
    source_bank(:)%wgt = (ngen+rampup)/dble(ngen)
    ngen = ngen + rampup
    up_sign = .true.

    ! update criteria
    if ( fmfdon ) then
        if ( ccrt == 1 ) then
            crt1 = 2D-1/sqrt(dble(ngen))
        end if
        if ( scrt == 1 ) then
            crt2 = 9D-2/sqrt(dble(ngen))
        end if
        crt3 = 9D-2/sqrt(dble(ngen))
    else
        if ( ccrt == 1 ) then
            crt1 = 9D-2/sqrt(dble(ngen))
        end if
        if ( scrt == 1 ) then
            crt2 = 3D-2/sqrt(dble(ngen))
        end if
        crt3 = 3D-2/sqrt(dble(ngen))
    end if

end subroutine

! =============================================================================
! CYCLECHANGE
! =============================================================================
subroutine CYCLECHANGE(cyc)
    integer, intent(in):: cyc

    n_inact  = cyc + 1
    n_totcyc = n_inact + n_act
    genup = .false.
    !if ( icore == score ) print*, "end", n_inact, ngen, n_totcyc

end subroutine

! =============================================================================
! MPRUP_DIST distributes the necessary parameters to the nodes
! =============================================================================
subroutine MPRUP_DIST(sz,source_bank)
    use MPI,         only: MPI_REAL8, MPI_COMM_WORLD
    use BANK_HEADER, only: Bank
    implicit none
    integer, intent(in):: sz
    type(Bank):: source_bank(:)
    integer:: TP
    integer:: WOR

    TP = MPI_REAL8
    WOR = MPI_COMM_WORLD

    do ii = 1, sz
    call MPI_BCAST(source_bank(ii)%wgt,1,TP,score,WOR,ierr)
    end do
    call MPI_BCAST(ngen,1,TP,score,WOR,ierr)
    call MPI_BCAST(genup,1,TP,score,WOR,ierr)
    call MPI_BCAST(mprupon,1,TP,score,WOR,ierr)
   
end subroutine

end module
