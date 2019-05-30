module MPRUP
    use VARIABLES, only: icore, score
    use ENTROPY
    use VARIABLES, only: ngen, n_totcyc, n_inact, n_act
    implicit none

    contains

! =============================================================================
! 
! =============================================================================
subroutine GENSIZE(cyc)
    implicit none
    integer, intent(in) :: cyc

!    ! FMFD divergence
!    if ( fmfdup ) then
!        call GENSIZEUP(wgt)
!        fmfdup = .false.
!        return
!    end if

    if ( genup ) then
    ! 1 convergence test
    !print*, "dsh", crt1, dshannon
    select case(ccrt)
    case(1); if ( dshannon > crt1 )     return
    case(2); if ( mod(cyc,int(crt1)) /= 0 ) return
    end select

    ! 2 stationary point
    select case(scrt)
    case(1)
        entrp2 = sum(shannon(1:elength))/dble(elength)
        dentrp = abs(entrp2-entrp1)/entrp2
        !print*, "den", crt2, dentrp, "*"
        if ( dentrp < crt2 ) then
            genup = .false.
            !mprupon = .false.
            !call CYCLECHANGE(cyc)
            return
        end if
        entrp1 = entrp2

        call GENSIZEUP
        !print*, "up", ngen

    case(2)
        if ( ngen >= crt2 ) then
            genup = .false.
            return
        end if
        call GENSIZEUP
        !print*, "up", ngen

    end select

    ! 3 stopping test
    else
        if ( dshannon > crt3 ) then
            mprupon = .false.
            call CYCLECHANGE(cyc)
            return
        end if

    end if


end subroutine


! =============================================================================
! GENSIZEUP increases the generation size and reallocate the fission neutron
! =============================================================================
subroutine GENSIZEUP
    use BANK_HEADER,    only: source_bank
    implicit none

    ngen = ngen + rampup
    source_bank(:)%wgt  = ngen/size(source_bank)

    ! update criteria
!    if ( fmfdon ) then
!        if ( ccrt == 1 ) then
!            crt1 = 2D-1/sqrt(dble(ngen))
!        end if
!        if ( scrt == 1 ) then
!            crt2 = 9D-2/sqrt(dble(ngen))
!        end if
!        crt3 = 9D-2/sqrt(dble(ngen))
!    else
        if ( ccrt == 1 ) then
            crt1 = 9D-2/sqrt(dble(ngen))
        end if
        if ( scrt == 1 ) then
            crt2 = 3D-2/sqrt(dble(ngen))
        end if
        crt3 = 3D-2/sqrt(dble(ngen))
!    end if

end subroutine

! =============================================================================
! CYCLECHANGE
! =============================================================================
subroutine CYCLECHANGE(cyc)
    integer, intent(in):: cyc

    n_inact  = cyc + 1
    n_totcyc = n_inact + n_act
    !print*, "end", n_inact, ngen

end subroutine

end module
