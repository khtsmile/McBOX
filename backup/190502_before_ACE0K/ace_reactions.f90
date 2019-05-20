module ace_reactions
use omp_lib
use constants,         only : PI, wgt_min, barn, tiny_bit
use variables,         only : E_mode, keff, k_col, k_tl
use material_header 
use ace_header
use ace_xs
use randoms,          only : rang, rand_vec
use scattering_laws 
use particle_header
use bank_header,     only: fission_bank, thread_bank, bank_idx

implicit none 

contains

! ================================================== !
!    
! ================================================== !
subroutine collision_CE (p) 
    type(particle), intent(inout) :: p
    integer :: iso, i, i_iso
    real(8) :: rn, el, noel, r, sigt_sum, temp, sum1, sum2
    real(8) :: micro_xs(6), macro_xs(4)
    ! * microscopic cross section
    ! 1 : total
    ! 2 : elastic
    ! 3 : absorption
    ! 4 : fission
    ! 5 : nufission
    ! 6 : thermal elastic
    
    p%n_collision = p%n_collision + 1
    p % n_coord = 1
    !===============================================
    ! Sample a target isotope in the mixture
    macro_xs = getMacroXS(materials(p%material), p%E)
    rn = rang(); temp = 0 
    do i = 1, materials(p%material)%n_iso
        micro_xs = getMicroXS( materials(p%material)%ace_idx(i), p%E)
        ! S(a,b)
        call GET_SAB_MIC(materials(p%material),i,p%E,micro_xs)
        temp = temp + micro_xs(1)*materials(p%material)%numden(i)*barn
        if ( rn < temp/macro_xs(1) ) then
            iso = materials(p%material)%ace_idx(i)
            i_iso = i
            if ( materials(p%material)%sab  .and. ace(iso)%sab_iso /= 0 &
                .and. p%E < 4D-6 ) p%yes_sab = .true.
            exit
        endif
    enddo
    
    !> Collision estimator 
    !$omp atomic
    k_col = k_col + p%wgt * macro_xs(4)/macro_xs(1)
    
    call fissionSite_CE(p, iso, micro_xs)
    !!===============================================
    !Sampling reaction: elastic vs.non-elastic
    el   = micro_xs(2)
    noel = micro_xs(1)-micro_xs(3)-el
    if (abs(noel) < 1.d-5) noel = 0 
    r = rang()*(noel+el)-el

    if ( ace(iso)%nxs(5) == 0 .or. r <= 0D0 ) then
        if ( p%yes_sab ) then
        call SAB_CE(p,iso,micro_xs(2),micro_xs(6))
        else
        call elastic_CE (p, iso)
        end if
    else
        call notElastic_CE (p, iso)
    end if
    
    p%wgt = p%wgt * ((el+noel)/micro_xs(1)) !(1 - micro_xs(3)/micro_xs(1))
    
    
    if (p%E /= p%E) then 
        print *, ace(iso)%library, p%E
        print*, "energy"
        stop
    endif 
    
    call absorption_CE(p)

    
    
    !print '(I3,F10.5,A4,2F10.5)', p%n_collision, p%last_E, '->', p%E!, 100*(p%last_E-p%E)/p%last_E
    !print *, p%n_collision, p%last_E, '->', p%E
    
end subroutine


! =============================================================================
! SAB_CE
! =============================================================================
subroutine SAB_CE(p,iso,th,thel)
    type(particle), intent(inout):: p
    integer, intent(in):: iso
    real(8), intent(in):: th
    real(8), intent(in):: thel
    real(8):: re

    re = rang()*th-thel
    if ( re <= 0D0 ) then
        call SAB_EL_CE(p,iso)    ! elastic scattering
    else
        call SAB_IN_CE(p,iso)    ! inelastic scattering
    end if
    p%yes_sab = .false.

end subroutine


! =============================================================================
! SAB_EL_CE
! =============================================================================
subroutine SAB_EL_CE(p,iso)
    type(particle), intent(inout):: p
    integer, intent(in):: iso
    type(SAB_EL_ANG), pointer:: ab1
    type(SAB_EL_XS), pointer:: ab2
    integer:: isab, ierg    ! index of S(a,b) and energy
    integer:: iang          ! index of angle
    real(8):: ssum          ! summation (CDF)
    real(8):: rn            ! random number
    real(8):: mu            ! cosine angle
    real(8):: ipfac         ! interpolation factor
    real(8):: awr           ! atomic weight ratio
    real(8):: aa            ! parameter

    isab = ace(iso)%sab_iso
    ab1 => sab(isab)%itca
    ab2 => sab(isab)%itce
    
    ! energy index
    call GET_IERG_SABE(isab,ierg,p%e)
    
    ! angle index
    iang = int(rang()*(sab(isab)%nxs(6)+1))+1
    
    ! interpolation factor
    ipfac=max(0D0,min(1D0,(p%e-ab2%erg(ierg))/(ab2%erg(ierg+1)-ab2%erg(ierg))))

    ! outgoing angle
    mu = ab1%ang(ierg,iang) + ipfac*(ab1%ang(ierg+1,iang)-ab1%ang(ierg,iang))
    
    ! coordinate change
    awr = ace(iso)%atn
    aa = 1D0+awr*(awr+2D0*mu)
    p%E = p%E*aa / ((1D0+awr)*(1D0+awr))
    mu = (1D0+mu*awr)/sqrt(aa)
    p%coord(1)%uvw = rotate_angle(p%coord(1)%uvw,mu)
    
    if ( associated(ab1) ) nullify(ab1)
    if ( associated(ab2) ) nullify(ab2)

end subroutine


! =============================================================================
! SAB_IN_CE
! =============================================================================
subroutine SAB_IN_CE(p,iso)
    type(particle), intent(inout):: p
    integer, intent(in):: iso
    type(SAB_INEL_E), pointer:: ab1
    type(SAB_INEL_XS), pointer:: ab2
    integer:: isab, ierg, ierg2, iang   ! index of S(a,b), energy, angle
    real(8):: mu    ! cosine angle
    real(8):: Ecm   ! energy (COM)
    real(8):: Eout  ! outgoing energy
    real(8):: ipfac ! interpolation factor
    real(8):: awr   ! atomic weight ratio

    isab = ace(iso)%sab_iso
    ab1 => sab(isab)%itxe
    ab2 => sab(isab)%itie

    ! energy & angle indices
    call GET_IERG_SABI(isab,ierg,p%e)
    ierg2 = int(rang()*(sab(isab)%nxs(4)))+1
    call SKEWED_SECOND(sab(isab)%nxs(4),ierg2)
    iang = int(rang()*(sab(isab)%nxs(3)+1))+1

    ! interpolation factor
    ipfac=max(0D0,min(1D0,(p%e-ab2%erg(ierg))/(ab2%erg(ierg+1)-ab2%erg(ierg))))

    ! outgoing energy
    p%E = ab1%erg(ierg,ierg2)+ipfac*(ab1%erg(ierg+1,ierg2)-ab1%erg(ierg,ierg2))

    ! outgoing angle
    mu = ab1%ang(ierg,ierg2,iang) + &
        ipfac*(ab1%ang(ierg+1,ierg2,iang)-ab1%ang(ierg,ierg2,iang))

    ! coordinate change
    p%coord(1)%uvw = rotate_angle(p%coord(1)%uvw,mu)
    
    if ( associated(ab1) ) nullify(ab1)
    if ( associated(ab2) ) nullify(ab2)

end subroutine


! =============================================================================
! SKEWED_SECOND
! =============================================================================
subroutine SKEWED_SECOND(ne,ierg2)
    integer, intent(in) :: ne
    integer, intent(out):: ierg2
    real(8):: pp    ! proability

    pp = rang()*(ne-3)
    if ( pp > 1D0 ) then
        ierg2 = int(pp) + 2
    elseif ( pp > 6D-1 ) then
        ierg2 = ne - 1
    elseif ( pp > 5D-1 ) then
        ierg2 = ne
    elseif ( pp > 1D-1 ) then
        ierg2 = 2
    else
        ierg2 = 1
    end if

end subroutine



! ================================================== !
!    acecos() samples the scattering cosine of the 
!    secondary neutron in the CM frame.
! ================================================== !
subroutine acecos (erg, iso, mu, iMT)
    !type(particle), intent(inout) :: p
    real(8), intent(in) :: erg 
    integer, intent(in) :: iso, iMT
    real(8), intent(inout):: mu
    real(8) :: LOCB, LC
    real(8) :: rn1, rn2 
    real(8) :: ipfac
    real(8) :: temp
    real(8) :: CSOUT1, CSOUT2, PDF1, PDF2, CDF1, CDF2
    integer :: NE, pt1, pt2, pt3, IE, i, JJ, ilaw
    type (AngularDist), pointer :: an
    real(8), allocatable :: P_tbl(:), CSOUT(:), PDF(:), CDF(:)
    integer :: K, NP 

    an => ace(iso)%ang(iMT)
    LOCB = ace(iso) % ang_flag(iMT) 
    if (LOCB == 0)  then
        !> isotropic scattering (CM system)
        mu = 2.0d0*rang()-1.0d0

    elseif (LOCB == -1) then 
        ! No angular distribution data are given for this reaction in the AND Block. 
        ! Angular distribution data are specified through LAW_i=44 in the DLW Blaock.
        !print *, 'law 44 called'!, p%last_E , '->', p%E, 100*(p%last_E - p%E)/p%last_E        
        do i = 1, ace(iso)%pneg(1)%nlaw
            if (ace(iso)%pneg(1)%dist(i)%law == 44) then 
                ilaw = i
                exit 
            endif
        enddo 
        call LAW44_ANG(erg, iso, mu, ace(iso)%pneg(1)%dist(ilaw), iMT)
        
    elseif (LOCB > 0) then 
        !print *, 'tabular or 32 equiprobable', LOCB
        !> Find energy index
        NE  = an%NE 
        pt1 = 1; pt2 = NE
        do 
            if (pt2-pt1 ==1) exit 
            pt3 = (pt1+pt2)/2
            if(erg >= an%E(pt3)) then 
                pt1 = pt3 
            else 
                pt2 = pt3
            endif 
        enddo 
        
        ipfac = max(0.d0, min(1.d0,(erg-an%E(pt1))/(an%E(pt2)-an%E(pt1))))
        IE = pt1
        if (rang()  < ipfac) IE = pt2
        !if (erg < an%E(1)) IE = 1
        
        
        LC = an%dist_flag(IE) 
        if (LC == 0) then   
            !> isotropic scattering (CM system)
            mu = 2.0d0*rang()-1.0d0
        elseif (LC > 0) then  !> 32 Equiprobable bin distribution 
            rn1 = rang()
            i  = int(32.0*rn1)
            allocate(P_tbl(33)) 
            !if( rang()*(an%E(pt2)-an%E(pt1)) < p%E-an%E(pt1) )  IE=pt2
            P_tbl(1:33) = an % dist(IE) % LDAT(1:33)
            mu = P_tbl(i) + (32.*rn1-i)*(P_tbl(i+1) - P_tbl(i))
            !mu = an%dist(IE)%P(i) + (32.*rn1-i)*(an%dist(IE)%P(i+1) - an%dist(IE)%P(i))
        else  !> Tabular probability angular distribution        
            NP = an % dist(IE) % LDAT(2)  !(size(an % dist(IE) % LDAT)-2)/3
            allocate(CSOUT(1:NP)); CSOUT(1:NP) = an % dist(IE) % LDAT(3:3+NP-1)
            allocate(PDF(1:NP)); PDF(1:NP) = an % dist(IE) % LDAT(3+NP:3+2*NP-1)
            allocate(CDF(1:NP)); CDF(1:NP) = an % dist(IE) % LDAT(3+2*NP:3+3*NP-1)
                        
            rn1 = rang() 
            pt1 = 1; pt2 = NP
            do 
                if (pt2-pt1 ==1) exit 
                pt3 = (pt1+pt2)/2 
                if(rn1 >= CDF(pt3)) then 
                    pt1 = pt3
                else 
                    pt2 = pt3
                endif
            enddo 
            
            CSOUT1    = CSOUT(pt1)
            CSOUT2    = CSOUT(pt2)
            PDF1    = PDF(pt1)
            PDF2    = PDF(pt2)
            CDF1    = CDF(pt1)
            CDF2    = CDF(pt2)
            
            JJ = an % dist(IE) % LDAT(1) !an % dist(IE) % JJ
            temp = (PDF2-PDF1)/(CSOUT2-CSOUT1)
            if (JJ == 1 .or. temp == 0) then !> Histogram 
                mu = CSOUT1 + (rn1-CDF1)/PDF1
            elseif (JJ == 2) then !> Linear-Linear 
                mu = CSOUT1+(sqrt(max(0.0d0,PDF1**2 +2.0d0*temp*(rn1-CDF1)))-PDF1)/temp
            else
                print *, "ERROR :: UNKNOWN INTERPOLATION TYPE IN SCATTERING ANGLE DISTRIBUTIOIN"
            endif
        endif         
        
    else 
        print *, "ERROR: NO SCATTERING ANGLE DISTRIBUTION", ace(iso)%library, ace(iso)%MT(iMT)
        stop
    endif
    
    
end subroutine

! ================================================== !
!    elastic_CE() handles elastic scattering of the 
!    particle 
!    Output :: Angle & Energy of p object
! ================================================== !
subroutine elastic_CE (p, iso)
    type(particle), intent(inout) :: p
    integer, intent(in) :: iso
    real(8) :: mu
    integer :: iMT, i 
    
    iMT = 0
    p%last_E = p%E
    !sample the neutron output direction and calculate its energy.
    !> elastic scattering is always considered as CM frame
    call acecos(p%E, iso, mu, iMT) !scattering angle in COM    
    if (abs(mu) > 1. ) then 
        !print *, "cos larger than 1"
        mu = sign(1.0d0, mu)
    endif
    
    call directionEnergy (p, mu, iMT, iso) 
    
end subroutine


! ================================================== !
!    notElastic_CE() handles all scattering event 
!    excluding elastic scattering. 
!    Output :: Angle & Energy of p object
! ================================================== !
subroutine notElastic_CE (p,iso)
    type(particle), intent(inout) :: p
    type (AngularDist), pointer :: an
    type (EnergyDist),  pointer :: eg
    type (CrossSectionDataForm), pointer :: sigmt
    integer, intent(in) :: iso
    integer :: i, j
    integer :: iMT, MT 
    integer :: pt1, pt2, pt3
    integer :: ierg
    integer :: law, ilaw        ! collision law
    real(8) :: ipfac    ! interpolation factor
    real(8) :: F         ! collision probability
    real(8) :: rn 
    real(8) :: mu 
    real(8) :: erg
    real(8) :: sig_arr(1:ace(iso)%NXS(5)), sig_sum, temp
    real(8) :: u, v, w, phi
    
    !real(8) :: micro(5)
    
    
    
    erg = p%E
    call getierg(iso,ierg,erg)
    ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
    
    sig_arr(:) = 0
    do i = 1, ace(iso)%NXS(5) !> through the reaction types...
    !> Determine reaction type number from 1 to NXS(5)
        !print *, int(ace(iso)%MT(i)), int(ace(iso)%ty(i)), ace(iso)%sig_MT(i)%cx(ierg), ipfac
    
        if (int(ace(iso)%ty(i))==19) cycle
        !> 1. locate energy grid in SIG(i) block 
        sigmt => ace(iso)%sig_MT(i)
        ! if p%E is outside the energy grid :: cycle
        !print *, sigmt%IE, ierg, (sigmt%IE+sigmt%NE-1)
        if (ierg >= (sigmt%IE+sigmt%NE-1) .or. ierg < sigmt%IE  ) cycle 
        
        !> 2. calculate XS for the reaction type
        sig_arr(i) = sigmt%cx(ierg) + ipfac*(sigmt%cx(ierg+1)-sigmt%cx(ierg))
        
        !print *, i, sigmt%cx(ierg), ipfac
        !print *, i, int(ace(iso)%MT(i)), sig_arr(i) 
    enddo 
    
    
    rn = rang()
    sig_sum = sum(sig_arr(:)); temp = 0; iMT = -10
    do i = 1, ace(iso)%NXS(5) !> through the reaction types...
        !print *, i, sig_arr(i)
        temp = temp + sig_arr(i)
        if (rn < temp/sig_sum) then
            MT = ace(iso)%MT(i)
            iMT = i
            exit
        endif 
    enddo 
    
    rn = rang()
    !> determine collision law number
    eg => ace(iso) % pneg(iMT)
    law = -1
    !print *, iMT, 'MT num', ace(iso)%MT(iMT), 'nlaw', eg%nlaw
    F = 0.
    law_search: do ilaw = 1, eg%nlaw
        !> Binary search to corresponding energy grid index
        pt1 = 1; pt2 =  eg % dist(ilaw) % NE
        do 
            if(pt2-pt1 == 1) exit 
            pt3 = (pt2+pt1)/2 
            if (p%E >= eg%dist(ilaw)%E(pt3) ) then 
                pt1 = pt3 
            else 
                pt2 = pt3
            endif
        enddo
        !print *, 'inelastic' , eg%nlaw, pt1, pt2 , eg % dist(ilaw) % NE
        !> Interpolate F(E) (P(E) in ENDF Manual...)
        ipfac = max(0.d0, min(1.d0,(p%E-eg%dist(ilaw)%E(pt1))/(eg%dist(ilaw)%E(pt1+1)-eg%dist(ilaw)%E(pt1))))
        F     = F + eg%dist(ilaw)%F(pt1) + ipfac*(eg%dist(ilaw)%F(pt1+1)-eg%dist(ilaw)%F(pt1))
        
        if (rn < F) then 
            law = eg % dist(ilaw) % law
            exit law_search
        endif
    enddo law_search
    
    if (law < 0) then 
        print *, 'ERROR :: law not selected'
        print *, F, eg%dist(ilaw)%F(pt1), ipfac
        print *, ace(iso)%library, iMT, ace(iso)%TY(iMT), eg%nlaw, law
        stop
    endif 
    
    !print *, ace(iso)%library, 'MT num', ace(iso)%MT(iMT), ace(iso)%TY(iMT), 'law', law
    
    
    
    p%last_E = p%E
    !> Sample the scattering cosine
    if (law /= 44 .and. law /= 61) then 
        !print *, 'law', law, p%E
        call acecos (p%E, iso, mu, iMT)
    endif 
    !p%E = p%last_E
    
    
    !> Sample the outgoing energy by the specific scattering Law
    !> call the corresponding law subroutines
    call law_selector (p%E, iso, mu, eg%dist(ilaw), iMT, law)
    
    if (abs(mu) > 1. ) then 
        !print *, 'WARNING :: abs cosine larger than 1', mu 
        mu = sign(1.0d0, mu)
    endif
    
    call directionEnergy (p, mu, iMT, iso) 
    
    p % wgt = p%wgt * dble(abs(ace(iso)%TY(iMT)))

    
end subroutine


! ================================================== !
!    fissionSite_CE() samples the potential fission 
!    source through the implicit capture process. 
! ================================================== !
subroutine fissionSite_CE (p, iso, micro_xs)
    type(particle), intent(in) :: p
    real(8), intent(in) :: micro_xs(5) 
    integer, intent(in) :: iso
    real(8) :: sig_arr(1:ace(iso)%NXS(5)), sig_sum, temp
    type (CrossSectionDataForm), pointer :: sigmt
    type (EnergyDist),  pointer :: eg
    real(8) :: rn
    integer :: i_source, i, n, ierg
    integer :: iMT, MT
    integer :: pt1, pt2, pt3
    integer :: law, ilaw        ! collision law
    real(8) :: ipfac    ! interpolation factor
    real(8) :: F         ! collision probability
    real(8) :: erg_out, mu = 1
    
    
    n = int(p%wgt*(micro_xs(5)/micro_xs(1))*(1.0/keff) + rang())
    do i_source = 1, n
        
        !print *, 'fission site **********'
        bank_idx = bank_idx + 1
        thread_bank(bank_idx)%xyz = p%coord(1)%xyz
        thread_bank(bank_idx)%uvw = rand_vec()
        
        !> Sample fission neutron energy 
        call getierg(iso,ierg,p%E)
        ipfac = max(0.d0, min(1.d0,(p%E-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
        sig_arr(:) = 0
        do i = 1, ace(iso)%NXS(5) !> through the reaction types...
        !> Determine reaction type number from 1 to NXS(5)
            !> 1. locate energy grid in SIG(i) block 
            if (int(ace(iso)%ty(i))/=19) cycle  ! consider fission reaction only
            sigmt => ace(iso)%sig_MT(i)
            if (ierg >= (sigmt%IE+sigmt%NE-1) .or. ierg < sigmt%IE  ) cycle 
            !> 2. calculate XS for the reaction type
            sig_arr(i) = sigmt%cx(ierg) + ipfac*(sigmt%cx(ierg+1)-sigmt%cx(ierg))
        enddo 
        
        rn = rang()
        sig_sum = sum(sig_arr(:)); temp = 0; 
        do i = 1, ace(iso)%NXS(5) !> through the reaction types...
            temp = temp + sig_arr(i)
            if (rn < temp/sig_sum) then
                MT = ace(iso)%MT(i)
                iMT = i 
                exit
            endif 
        enddo
        
        rn = rang()
        !> determine collision law number
        eg => ace(iso) % pneg(iMT)
        law = -1
        !print *, eg%nlaw
        law_search: do ilaw = 1, eg%nlaw
            !> Binary search to corresponding energy grid index
            pt1 = 1; pt2 =  eg % dist(ilaw) % NE
            do 
                if(pt2-pt1 == 1) exit 
                pt3 = (pt2+pt1)/2 
                if (p%E >= eg%dist(ilaw)%E(pt3) ) then 
                    pt1 = pt3 
                else 
                    pt2 = pt3
                endif 
            enddo 
            
            !> Interpolate F(E) (P(E) in ENDF Manual...)
            ipfac = max(0.d0, min(1.d0,(p%E-eg%dist(ilaw)%E(pt1))/(eg%dist(ilaw)%E(pt1+1)-eg%dist(ilaw)%E(pt1))))
            F     = eg%dist(ilaw)%F(pt1) + ipfac*(eg%dist(ilaw)%F(pt1+1)-eg%dist(ilaw)%F(pt1))
            if (rn < F) then 
                law = eg % dist(ilaw) % law
                exit law_search
            endif 
        enddo law_search
        
        if (law < 0) then 
            print *, '**************************   law not selected', F, eg%dist(ilaw)%F(pt1), ipfac
            stop
        endif 
        erg_out = p%E
        !print *, 'input E ',erg_out
        !> call the corresponding law subroutines
        call law_selector (erg_out, iso, mu, eg%dist(ilaw), iMT, law)
        thread_bank(bank_idx)%E = erg_out
        
    enddo 
        
        
end subroutine

! ================================================== !
!    absorption_CE() kills the particle if the weight
!    is too small than the cutoff + Russian Roulette
! ================================================== !
subroutine absorption_CE (p)
    type(particle), intent(inout) :: p
    real(8) :: wgt_s
    
    !if (p%n_collision > 1000) then 
    !    p%wgt = 0 
    !    p%alive = .false.
    !    print *, 'killed for too many collisions'
    !endif
    
    !if (p%E < 1.0d-11) then 
    !    p%wgt = 0 
    !    p%alive = .false.
    !    !print *, 'killed for too small energy'
    !endif 
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

! ================================================== !
!    directionEnergy() calculates mu and E in the 
!    target-at-rest (lab) frame 
! ================================================== !
subroutine directionEnergy (p, mu, iMT, iso) 

    type(particle), intent(inout) :: p 
    real(8), intent(inout) :: mu
    integer, intent(in) :: iMT, iso 
    real(8) :: mu_CM, Eout_CM, Eout_lab, Ein, A
    real(8) :: u,v,w, u0, v0, w0, phi
    real(8) :: temp, val
    logical :: found = .false. 
    integer :: j 
    real(8) :: kT        ! temperature (MeV)
    real(8):: uvw(3)
    
    
    !> Elastic scattering 
    if (iMT == 0) then
        !print *, 'elastic'
        !> Elastic Scattering Energy calculation 
        
        mu_CM = mu
        Ein = p%E
        A = ace(iso)%atn
        kT = ace(iso)%temp
        
        ! target at rest
        if ( p%E > 4D2*kT .and. A > 1D0 .or. kT == 0D0 ) then
            temp = 1.+A*(A+2.*mu_CM)
            Eout_lab = Ein*temp/(1.+A)**2
            mu     = (1.+mu_CM*A)/sqrt(temp)
        
            p%E = Eout_lab 
            p%coord(1)%uvw = rotate_angle(p%coord(1)%uvw, mu)
        
        ! target in motion (free gas model)
        else
            call FREE_GAS(p,mu,kT,A)
        end if
        
        
    !> Inelastic scattering 
    elseif (ace(iso)%TY(iMT) < 0) then !> CM frame
        !> mu and p%E are in CM frame
        !> convert to lab frame (existing code )
        A = ace(iso)%atn
        Eout_CM = p%E
        Ein = p%last_E
        mu_CM = mu
        Eout_lab = Eout_cm + (Ein + 2.*mu_cm*(A+1.)*sqrt(Ein*Eout_cm))/(A+1.)**2
        mu         = mu_CM*sqrt(Eout_CM/Eout_lab) + sqrt(Ein/Eout_lab)/(A+1.)
        
        p%E = Eout_lab
        p%coord(1)%uvw = rotate_angle(p%coord(1)%uvw, mu)
        !do j = 1, p%n_coord
        !    p%coord(j)%uvw = rotate_angle (p%coord(j)%uvw, mu)
        !enddo 
        
    else!if(ace(iso)%TY(iMT)>0) then  !> Reaction was in TAR frame 
        !> mu and p%E are in lab frame
        p%coord(1)%uvw = rotate_angle(p%coord(1)%uvw, mu)

        !do j = 1, p%n_coord
        !    p%coord(j)%uvw = rotate_angle (p%coord(j)%uvw, mu)
        !enddo 
    endif 

    !do j = 1, p%n_coord
    !    if (p % coord(j) % uvw(1) /= p % coord(j) % uvw(1)) then 
    !        print *, j
    !        print *, p % coord(j) % uvw(:)
    !        print *, u, v, w
    !        print *, mu, sqrt(1-mu**2)
    !        print *, sqrt(1-w0**2), (1-w0**2)
    !        
    !        
    !        print *, iMT, ace(iso)%library
    !        stop
    !    endif 
    !enddo
    
end subroutine 

! =============================================================================
! FREE_GAS
! =============================================================================
subroutine FREE_GAS(p,mu,kT,A)
    type(Particle), intent(inout):: p
    real(8), intent(inout):: mu
    real(8), intent(in):: kT
    real(8), intent(in):: A
    real(8):: speedn                    ! neutron speed
    real(8):: uvw(3)                    ! incoming & outgoing direction
    real(8):: uvw_cm(3)                 ! direction of COM
    real(8):: v_t(3), v_n(3), v_cm(3)   ! velocity of target, neutron, COM
    real(8):: mu_cm                     ! cosine angle
    
    real(8), parameter:: m_u = 1.660540D-27 ! amu to Kg
    real(8), parameter:: m_n = 1.008664     ! neutron mass (amu)
    real(8), parameter:: mevj = 1.6022D-13  ! MeV to Joule

    ! -------------------------------------------------------------------------
    ! LAB system
    ! - speed of neutron
!    speedn = sqrt(p%E)
    speedn = sqrt(2D0*p%E*mevj/(m_u*m_n))   ! m/s
    ! - neutron incoming direction
    uvw = p%coord(1)%uvw
    ! - velocity of neutron
    v_n = speedn * uvw
    ! - target velocity
    call TARGETV(p%E,uvw,A,kT,v_t)
!    v_t = 0D0


    ! -------------------------------------------------------------------------
    ! COM system
    ! - velocity of COM
    v_cm = (v_n + A * v_t)/(A+1D0)
    ! - transform to COM
    v_n = v_n - v_cm
    ! - speed of neutron
    speedn = sqrt(dot_product(v_n,v_n))
    ! - scattering angle
    mu_cm = mu
    ! - direction cosine
    uvw_cm = v_n / speedn
    ! - neutron velocity
    v_n = speedn * rotate_angle(uvw_cm,mu_cm)


    ! -------------------------------------------------------------------------
    ! LAB system
    ! - transform to LAB
    v_n = v_n + v_cm
    ! - neutron speed
    speedn = sqrt(dot_product(v_n,v_n))
    ! - neutron outgoing energy
!    p%E = dot_product(v_n,v_n)
    p%E = (m_u*m_n)*speedn*speedn/(mevj*2D0)
    ! - angle test
    mu = dot_product(uvw,v_n) / speedn
    if ( abs(mu) > 1D0 ) mu = sign(1D0,mu)
    ! - neutron outgoing direction ( p%E = speed of neutron )
    p%coord(1)%uvw = v_n / speedn
!    p%coord(1)%uvw = uvw

!    print*, "E2", p%E

end subroutine

! =============================================================================
! TARGETV
! =============================================================================
subroutine TARGETV(E0,uvw,A,kT,v_t)
    real(8), intent(in):: E0        ! incident energy (MeV)
    real(8), intent(in):: uvw(3)    ! incident direction
    real(8), intent(in):: A         ! atomic weight ratio
    real(8), intent(in):: kT        ! temperature (MeV)
    real(8), intent(inout):: v_t(3) ! target velocity
    real(8):: urn1, urn2            ! uniform random number
    real(8):: accept                ! acceptance parameter
    real(8):: speedt                ! target speed
    real(8):: mut                   ! cosine angle
    real(8):: xx, xx2, yy           ! parameters
    real(8):: aa, bb, cc, ss

    real(8), parameter:: m_u = 1.660540D-27 ! AMU to Kg
    real(8), parameter:: m_n = 1.008664     ! neutron mass (amu)
    real(8), parameter:: mevj = 1.6022D-13  ! MeV to Joule
    integer:: ii

    ! parameters
    bb = sqrt(A*m_u*m_n/(2D0*kT*mevj))
    yy = sqrt(A*E0/kT)        ! (beta) X (v_n)
    aa = 2D0/(sqrt(pi)*yy+2D0)

    ! sampling target speed from Maxwellian distribution
    do
        ! 1) target speed
        urn1 = rang()
        urn2 = rang()
        if ( rang() < aa ) then
            xx2 = -log(urn1*urn2)
        else
            cc = cos(pi/2D0*rang())
            xx2 = -log(urn1)-log(urn2)*cc*cc
        end if
        xx = sqrt(xx2)  ! (beta) X (v_t)

        ! 2) cosine angle (isotropic)
        mut = 2D0*rang()-1D0

        ! 3) rejection method
        accept = sqrt(yy*yy+xx2-2D0*xx*yy*mut)/(xx+yy)
        if ( rang() < accept ) exit
    end do

!    do
!        if ( rang()*(yy+1.12837917D0) > yy ) then
!            xx2 = -log(rang()*rang())
!        else
!            do
!                urn1 = rang()
!                urn2 = rang()
!                ss = urn1*urn1+urn2*urn2
!                if ( ss <= 1D0 ) exit
!            end do
!            xx2 = -urn1*urn1*log(ss)/ss-log(rang())
!        end if
!
!        xx = sqrt(xx2)
!        mut = 2D0*rang()-1D0
!        accept = yy*yy+xx2-2D0*xx*yy*mut
!        if ( (rang()*(yy+xx))**2 <= accept ) exit
!
!    end do

    ! target speed
    speedt = xx / bb
    !speedt = xx*sqrt(kT*mevj/(A*m_u*m_n))
    !print*, "v_t", speedt
!    do ii = 1, 500
!    if ( derg*(ii-1) <= speedt .and. speedt < derg*ii ) cnt(ii) = cnt(ii)+1
!    end do

    ! target velocity
    v_t = speedt * rotate_angle(uvw,mut)

end subroutine


function rotate_angle (uvw0, mu) result (uvw)
    implicit none 
    
    real(8), intent(in) :: uvw0(3) 
    real(8), intent(in) :: mu 
    real(8) :: uvw(3)
    
    real(8) :: phi
    real(8) :: u0, v0, w0 
    real(8) :: a, b, sinphi, cosphi, temp
    integer :: i 
    
    u0 = uvw0(1); v0 = uvw0(2); w0 = uvw0(3)
    phi = 2*PI*rang()
    sinphi = sin(phi) 
    cosphi = cos(phi)
    a = sqrt(max(0.0d0,1-mu**2))
    b = sqrt(max(0.0d0, 1-w0**2))
    
    if (b > 1.0d-10) then 
        uvw(1) = mu*u0 + a*(u0*w0*cosphi-v0*sinphi)/b
        uvw(2) = mu*v0 + a*(v0*w0*cosphi+u0*sinphi)/b
        uvw(3) = mu*w0 - a*b*cosphi
    else 
        b = sqrt(max(0.0d0, 1-v0**2))
        uvw(1) = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
        uvw(2) = mu*v0 - a*b*cosphi
        uvw(3) = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    endif

    do i = 1, 3 
        if (uvw(i) /= uvw(i)) then 
            print*, "rotation error"
            print *, a, b
            print *, uvw(:) 
            stop
        endif 
    enddo 
    !temp = sqrt(uvw(1)**2 + uvw(2)**2 + uvw(3)**2) 
    !uvw(:) = uvw(:)/temp
    
    !print *, sqrt(uvw(1)**2 + uvw(2)**2 + uvw(3)**2) 
end function 

end module
