module ace_reactions
use omp_lib
use constants,         only : PI, wgt_min, barn, tiny_bit
use variables,         only : E_mode, keff, k_col, icore, score
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
    integer :: iso, i, i_iso, xn
    real(8) :: rn, el, noel, r, sigt_sum, temp, sum1, sum2
    real(8) :: micro_xs(5), macro_xs(5)
    real(8) :: ipfac
    integer :: ierg
    
    p%n_collision = p%n_collision + 1
    p % n_coord = 1
    xn = 1
    !===============================================
    ! Sample a target isotope in the mixture
    macro_xs = getMacroXS(materials(p%material), p%E)
    rn = rang(); temp = 0 
    do i = 1, materials(p%material)%n_iso
        micro_xs = getMicroXS( materials(p%material)%ace_idx(i), p%E)
        temp = temp + micro_xs(1)*materials(p%material)%numden(i)*barn
        if ( rn < temp/macro_xs(1) ) then
            iso = materials(p%material)%ace_idx(i)
            i_iso = i
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
    !noel = micro_xs(1)-micro_xs(3)-el
    !if (abs(noel) < 1.d-5) noel = 0 

    noel = 0 
    call getierg(iso,ierg,p%E)
    ipfac = max(0.d0, min(1.d0,(p%E-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
    do i = 1, ace(iso)%NXS(5) !> through the reaction types...
        if (abs(ace(iso)%TY(i)) == 19) cycle 
        noel = noel + ace(iso)%sig_MT(i)%cx(ierg) & 
                    + ipfac*(ace(iso)%sig_MT(i)%cx(ierg+1) - ace(iso)%sig_MT(i)%cx(ierg))
    enddo 

    r = rang()*(noel+el)-el
    if( ace(iso)%nxs(5) == 0 .or. r <= 0.0d0 ) then 
        call elastic_CE (p, iso)
    else
        call notElastic_CE (p, iso, xn)
    endif 
    
    
    p%wgt = p%wgt * ((el+noel)/micro_xs(1)) !(1 - micro_xs(3)/micro_xs(1))
    
    
    !> (n, xn) reaction
    p%wgt = p%wgt * dble(xn)
    
    call absorption_CE(p)
    

    !print '(I3,F10.5,A4,2F10.5)', p%n_collision, p%last_E, '->', p%E!, 100*(p%last_E-p%E)/p%last_E
    !print *, p%n_collision, p%last_E, '->', p%E
    
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
            
            CSOUT1 = CSOUT(pt1)
            CSOUT2 = CSOUT(pt2)
            PDF1   = PDF(pt1)
            PDF2   = PDF(pt2)
            CDF1   = CDF(pt1)
            CDF2   = CDF(pt2)
            
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
subroutine notElastic_CE (p,iso,xn)
    type(particle), intent(inout) :: p
    type (AngularDist), pointer :: an
    type (EnergyDist),  pointer :: eg
    type (CrossSectionDataForm), pointer :: sigmt
    integer, intent(in) :: iso
    integer, intent(inout) :: xn
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
    real(8) :: micro(5)
    
    
    
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
        temp = temp + sig_arr(i)
        if (rn < temp/sig_sum) then
            MT = ace(iso)%MT(i)
            iMT = i
            exit
        endif 
    enddo 
    
    if (iMT < 0) then 
        print *, '*********** warning :: iMT not selected ************'
        print *, 'rn',rn
        print *, 'NXS(5)', ace(iso)%NXS(5)
        print *, 'sig_arr', sig_arr(:)
        
        micro = getMicroXS( iso, p%E)
        
        print *, '  elastic', micro(2)
        print *, 'inelastic', micro(1)-micro(3)-micro(2)

        stop 
    endif
    
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
    
    xn = abs(ace(iso)%TY(iMT))
    if(xn > 4 .or. xn == 0) xn = 1
    
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
    if (n > 10) then 
        print *, 'fission site number ', n
        print *, 'particle wgt', p%wgt 
        print *, 'nusig_f', micro_xs(5) 
        print *, 'sig_t', micro_xs(1) 
        
        stop 
    endif 
    

    do i_source = 1, n
        !if (icore /= score) print *, bank_idx, n
        !print *, bank_idx, n , p%wgt
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
    !    !print *, 'killed for too many collisions'
    !endif
    !if (p%E < 1.0d-13) then 
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
    
    real(8) :: vel       ! magnitude of velocity
    real(8) :: v_n(3)    ! velocity of neutron
    real(8) :: v_cm(3)   ! velocity of center-of-mass
    real(8) :: v_t(3)    ! velocity of target nucleus
    real(8) :: uvw_cm(3) ! directional cosines in center-of-mass
    real(8) :: uvw(3)
    
    
    !> Elastic scattering 
    if (iMT == 0) then
        !print *, 'elastic'
        !> Elastic Scattering Energy calculation 
        
        mu_CM = mu
        Ein = p%E
        A = ace(iso)%atn
        
        temp = 1.+A*(A+2.*mu_CM)
        Eout_lab = Ein*temp/(1.+A)**2
        mu     = (1.+mu_CM*A)/sqrt(temp)
        p%E = Eout_lab 
        
        p%coord(1)%uvw = rotate_angle(p%coord(1)%uvw, mu)
        
        !uvw = p%coord(1)%uvw
        !vel = sqrt(p%E)
        !A = ace(iso)%atn
        !
        !! Neutron velocity in LAB
        !v_n = vel * uvw
        !
        !v_t = 0.
        !
        !! Velocity of center-of-mass
        !v_cm = (v_n + A*v_t)/(A + 1.)
        !
        !! Transform to CM frame
        !v_n = v_n - v_cm
        !
        !! Find speed of neutron in CM
        !vel = sqrt(dot_product(v_n, v_n))
        !! Determine direction cosines in CM
        !uvw_cm = v_n/vel
        !
        !! Rotate neutron velocity vector to new angle -- note that the speed of the
        !! neutron in CM does not change in elastic scattering. However, the speed
        !! will change when we convert back to LAB
        !v_n = vel * rotate_angle(uvw_cm, mu_cm)
        !
        !! Transform back to LAB frame
        !v_n = v_n + v_cm
        !
        !p%E = dot_product(v_n, v_n)
        !vel = sqrt(p%E)
        !
        !! compute cosine of scattering angle in LAB frame by taking dot product of
        !! neutron's pre- and post-collision angle
        !!mu_lab = dot_product(uvw, v_n) / vel
        !
        !! Set energy and direction of particle in LAB frame
        !uvw = v_n / vel
        !    
        !do j = 1, p%n_coord
        !    p%coord(j)%uvw = uvw
        !enddo 
        
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
            print *, b
            print *, uvw(:) 
            stop
        endif 
    enddo 
    !temp = sqrt(uvw(1)**2 + uvw(2)**2 + uvw(3)**2) 
    !uvw(:) = uvw(:)/temp
    
    !print *, sqrt(uvw(1)**2 + uvw(2)**2 + uvw(3)**2) 
end function 

end module
