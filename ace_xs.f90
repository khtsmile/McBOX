module ace_xs

use constants, only : barn
use variables, only : E_mode
use material_header 
use ace_header 
use ace_module 


implicit none 

contains

function getMacroXS (mat, erg) result (macro_xs)
	type(Material_CE), intent(in) :: mat
	real(8), intent(in) :: erg
	real(8) :: macro_xs(5)
	
	integer :: i 
	integer :: i_iso, iso_, ierg_
	integer :: pt1, pt2, pt3, pt4
	real(8) :: ipfac
	real(8) :: micro_t, micro_d, micro_f, micro_nuf, micro_a, micro_el, micro_xn
	real(8) :: macro_t, macro_f, macro_nuf, macro_a, macro_qf
	real(8) :: xn_xs(4) 
	
	macro_t   = 0.0d0
	macro_a   = 0.0d0
	macro_f   = 0.0d0
	macro_nuf = 0.0d0
	macro_qf  = 0.0d0
	xn_xs(:)  = 0.0d0
	
	!print *, mat%mat_name, mat%n_iso
	do i_iso = 1, mat%n_iso
	
		iso_ = mat%ace_idx(i_iso)
	
		call getierg(iso_,ierg_,erg)
		
		ipfac = max(0.d0, min(1.d0,(erg-ace(iso_)%E(ierg_))/(ace(iso_)%E(ierg_+1)-ace(iso_)%E(ierg_))))
		
		!==============================================================
		! Microscopic XS
		micro_t   = ace(iso_)%sigt(ierg_) + ipfac*(ace(iso_)%sigt(ierg_+1)-ace(iso_)%sigt(ierg_))
		micro_d   = ace(iso_)%sigd(ierg_) + ipfac*(ace(iso_)%sigd(ierg_+1)-ace(iso_)%sigd(ierg_))
		micro_f   = 0.d0
		micro_nuf = 0.d0
		micro_a   = micro_d
		
		
		! Fissionable Material
		!if(ace(iso_)%jxs(21)/=0) then
		!if(allocated(ace(iso_)%sigf)) then 
		if(ace(iso_)%jxs(21)/=0 .or. allocated(ace(iso_)%sigf)) then
			micro_f   = ace(iso_)%sigf(ierg_) + ipfac*(ace(iso_)%sigf(ierg_+1)-ace(iso_)%sigf(ierg_))
			micro_nuf = getnu(iso_,erg)*micro_f
			micro_a   = micro_d + micro_f
		endif
		!micro_el = ace(iso_)%sigel(ierg_) + ipfac*(ace(iso_)%sigel(ierg_+1)-ace(iso_)%sigel(ierg_))
		 
		!>Summation for macroscopic cross sections
		macro_t   = macro_t   + mat%numden(i_iso) * micro_t   * barn
		macro_a   = macro_a   + mat%numden(i_iso) * micro_a   * barn
		macro_f   = macro_f   + mat%numden(i_iso) * micro_f   * barn
		macro_nuf = macro_nuf + mat%numden(i_iso) * micro_nuf * barn
		macro_qf  = macro_qf  + mat%numden(i_iso) * micro_f   * barn * ace(iso_)%qval
		
		!> Macro_xs of Sig_abs is only used for CMFD, (n,xn) XS is subtracted. 
		do i = 1, ace(iso_)%NXS(5) !> through the reaction types...
			pt1 = abs(ace(iso_)%TY(i))
			if (pt1 > 1 .and. pt1 < 5) then 
				micro_xn   = ace(iso_)%sig_MT(i)%cx(ierg_) & 
							+ ipfac*(ace(iso_)%sig_MT(i)%cx(ierg_+1) - ace(iso_)%sig_MT(i)%cx(ierg_))
				
				xn_xs(pt1) = xn_xs(pt1) + mat%numden(i_iso) * micro_xn * barn
			endif
		enddo
	enddo 

	do i = 2, 4 
		macro_a = macro_a - (dble(i)-1.0d0)*xn_xs(i)
	enddo 	
	
	macro_xs(1) = macro_t  
	macro_xs(2) = macro_a  
	macro_xs(3) = macro_f  
	macro_xs(4) = macro_nuf
	macro_xs(5) = macro_qf
	
	
	
	
end function

function getmicroXS (iso, erg) result (micro_xs)
	integer, intent(in) :: iso
	real(8), intent(in) :: erg
	real(8) :: micro_xs(5)
	
	integer :: ierg_, i 
	integer :: pt1, pt2, pt3, pt4
	real(8) :: ipfac
	real(8) :: micro_t, micro_d, micro_f, micro_nuf, micro_a, micro_el
	
	
	call getierg(iso,ierg_,erg)
	
	ipfac = max(0.d0, min(1.d0,(erg-ace(iso)%E(ierg_))/(ace(iso)%E(ierg_+1)-ace(iso)%E(ierg_))))
	!==============================================================
	! Microscopic XS
	micro_t   = ace(iso)%sigt(ierg_) + ipfac*(ace(iso)%sigt(ierg_+1)-ace(iso)%sigt(ierg_))
	micro_d   = ace(iso)%sigd(ierg_) + ipfac*(ace(iso)%sigd(ierg_+1)-ace(iso)%sigd(ierg_))
	micro_f   = 0.d0
	micro_nuf = 0.d0
	micro_a   = micro_d
	
	! Fissionable Material
	if(ace(iso)%jxs(21)/=0 .or. allocated(ace(iso)%sigf)) then
		micro_f   = ace(iso)%sigf(ierg_) + ipfac*(ace(iso)%sigf(ierg_+1)-ace(iso)%sigf(ierg_))
		micro_nuf = getnu(iso,erg)*micro_f
		micro_a   = micro_d + micro_f
	endif
	micro_el = ace(iso)%sigel(ierg_) + ipfac*(ace(iso)%sigel(ierg_+1)-ace(iso)%sigel(ierg_))
	
	
	micro_xs(1) = micro_t
	micro_xs(2) = micro_el
	micro_xs(3) = micro_a
	micro_xs(4) = micro_f
	micro_xs(5) = micro_nuf
	
	!print '(5F10.4)', micro_t, micro_el, micro_a, micro_f
	
	
end function


function getnu (iso_,erg0) result (nu)
	integer, intent(in) :: iso_
	real(8), intent(in) :: erg0
	real(8) :: nu
	integer :: nublock, n, i, NC, ierg, pt1, pt2, pt3 
	real(8) :: ipfac
	type (AceFormat), pointer :: ac
	
	
	ac => ace(iso_)
	
	if(ac % nu_block_exist == .false.) return
	
	select case( ac % nu_tot_flag )
	case(1) !> polynomial function form
		nu = 0.0d0
		do i = 1, ac % nu_tot % NC
			nu = nu + ac%nu_tot%C(i)*erg0**(i-1)
		enddo
		
	case(2) !> tabular data form
		!ac%nu_tot%E(:) !> nu energy grid 
		!ac%nu_tot%F(:) !> corresponding nu
		
		!if (ac%nu_tot%NR /= 0) print *, "WARNING :: nu is not lin-lin for", ace(iso_)%library
		! 1. binary search to find ierg of erg in E(:) 
		pt1 = 1
		pt2 = ac % nu_tot % NE
		Do
			if(pt2 - pt1 == 1) exit
			pt3 = (pt2 + pt1)/2
			if(erg0 >= ac%nu_tot%E(pt3)) then
			  pt1 = pt3
			else
			  pt2 = pt3
			endif
		Enddo
		ierg = pt1 !store low bound energy index 
		
		!if (ac % nu_tot % NR /= 0) print *, "ENDF interpolation required for Nu "
		! 2. calculate interpolation factor
		ipfac = max(0.d0, min(1.d0,(erg0-ac%nu_tot%E(ierg))/(ac%nu_tot%E(ierg+1)-ac%nu_tot%E(ierg))))
		
		! 3. linear-linear interpolation
		nu = ac%nu_tot%F(ierg) + ipfac*(ac%nu_tot%F(ierg+1)-ac%nu_tot%F(ierg))

	end select
	
end function


subroutine getierg(iso_,ierg_,erg0)
	implicit none
	real(8), intent(in) :: erg0
	integer, intent(in) :: iso_
	integer, intent(out) :: ierg_

	integer :: pt1, pt2, pt3, uidx


	!Energy grid search algorithm
	!erg0 = incidient neutron energy in lab system
	!iso_ = collision isotope index
	!ierg_ = low bound energy grid index
	
	
	uidx = 1 + int(log10(erg0/Emin)/udelta)
	if(uidx<1) uidx = 1
	pt1 = ugrid(uidx-1,iso_)
	pt2 = min(ugrid(uidx,iso_)+1,ace(iso_)%nxs(3))
	
	!pt1 = 1 
	!pt2 = size(ace(iso_)%E)
	
	if(pt1==pt2) then
	  pt1 = pt1 - 1
	else
	  Do
		if(pt2 - pt1 == 1) exit
		pt3 = (pt2 + pt1)/2
		if(erg0 >= ace(iso_)%E(pt3)) then
		  pt1 = pt3
		else
		  pt2 = pt3
		endif
	  Enddo
	endif
	ierg_ = pt1 !store low bound energy index 
	
	!if(erg0 < 1.d-11) print *, erg0, ierg_ !, ace(iso_)%E(ierg_), ace(iso_)%E(ierg_+1)
end subroutine 


subroutine setugrid
	implicit none
	real(8) :: Etmp
	integer :: i, j, k, iso_, idx


	if(E_mode==0) return
	!Set ugrid to accelerate energy-grid search
	Emin = 1.d-11
	allocate(ugrid(0:nugrid,1:num_iso))
	!print *, 'ugrid size', nugrid, num_iso
	do iso_ = 1, num_iso
	  ugrid(0,iso_) = 1
	  ugrid(nugrid,iso_) = ace(iso_)%nxs(3)
	end do
	udelta = log10(Emax/Emin)/dble(nugrid)

	do iso_ = 1, num_iso
	  idx = 1
	  do i=1, nugrid-1
		Etmp = Emin*10.d0**(dble(i)*udelta)
		if(Etmp > ace(iso_)%E( ace(iso_)%NXS(3) )) then 
		  idx = ace(iso_)%nxs(3)
		  go to 10  
		end if
		do
		  if(Etmp < ace(iso_)%E(idx)) go to 10
		  idx = idx + 1
		end do
10	    ugrid(i,iso_) = idx - 1
	  end do
	enddo

!> Only if OTF-DB is Enabled

!	if(num_iso0K) then
!	  allocate(ugrid0K(0:nugrid,1:num_iso0K))
!	  do iso_ = 1, num_iso0K
!		ugrid0K(0,iso_) = 1
!		ugrid0K(nugrid,iso_) = ace0K(iso_)%nxs(3)
!	  end do
!	  
!	  do iso_ = 1, num_iso0K
!		idx = 1
!		do i=1, nugrid-1
!		  Etmp = Emin*10.d0**(dble(i)*udelta)
!		  if(Etmp > ace0K(iso_)%xss(ace0K(iso_)%nxs(3))) then
!			idx = ace0K(iso_)%nxs(3)
!			go to 20  
!		  end if
!		  do
!			if(Etmp < ace0K(iso_)%xss(idx)) go to 20
!			idx = idx + 1
!		  end do
!20	    ugrid0K(i,iso_) = idx - 1
!		end do
!	  end do
!	end if


	if(icore==score) print *, "   Setting ugrid..."
	
	
	!do i=1, nugrid-1
	!	write(99, *) i, Emin*10.d0**(dble(i)*udelta)
	!enddo 
	!stop
	
end subroutine 



end module 