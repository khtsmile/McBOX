module ace_header

implicit none
character(100):: library_path    !> Path of ACE library


type TabularDataForm
  integer :: NR        !> number of interpolation regions
  real(8), allocatable :: NBT(:) !> ENDF interpolation parameters, NBT(I), I=1,NR
  real(8), allocatable :: INT(:) !> ENDF interpolation parameters, INT(I), I=1,NR
  integer :: NE        !> number of energies
  real(8), allocatable :: E(:)   !> number of energies, E(I), I=1,NE
  real(8), allocatable :: F(:)   !> energy dependent data, F(I), 1=1,NE
end type 



type, extends (TabularDataForm) :: NuTotDataForm 
  integer :: NC                !> number of coefficients
  real(8), allocatable :: C(:) !> Coefficients C(I), I=1,NC
end type


type CrossSectionDataForm
  integer :: IE        !> first energy grid index correspond to E(:) in ESZ Block
  integer :: NE        !> number of consecutive entries
  real(8), allocatable :: cx(:) !> cross sections, sig(I), I=IE,IE+NE-1
end type



type AngularDistDataForm
  !> case( dist_flag(IE) = 0 ) :: isotropic distribution, no data is needed. 
  !> case( dist_flag(IE) < 0 ) :: tabular angular distribution
  !>   LDAT( 1 : 2 + 3*NP )
  !>   LDAT(1) = JJ, interpolation flag
  !>     LDAT(1)=1 :: histogram
  !>     LDAT(1)=2 :: lin-lin
  !>   LDAT(2) = NP, number of points in the distribution
  !>   LDAT( 3      : 3+NP-1   ) = CSOUT(1:NP), cosine scattering angular grid
  !>   LDAT( 3+NP   : 3+2*NP-1 ) = PDF(1:NP), probability density function 
  !>   LDAT( 3+2*NP : 3+3*NP-1 ) = CDF(1:NP), cumulative density function
  !> case( dist_flag(IE) > 0 ) :: 32 equiprobable cosine bins for scattering 
  !real(8), allocatable :: LDAT(:) 
  
  !> 32 equiprobable distribution 
  !real(8), allocatable :: P(:) 
  
  !> tabular angular distribution
  !integer :: JJ = 0, NP = 0
  !real(8), allocatable :: CSOUT(:) 
  !real(8), allocatable :: PDF(:) 
  !real(8), allocatable :: CDF(:) 
  
  real(8), allocatable :: LDAT(:)
end type


type AngularDist
  integer :: flag !> Reaction MT dependent flag !> flag = 0  :: no angular distribution data (isotropic distribution)
                                                !> flag = -1 :: angular distribution data are given in DLW Block (LAW=44)
                                                !> flag > 0  :: angular distribution data are given
  integer :: NE     !> number of energies at which angular distributions are tabulated
  real(8), allocatable :: E(:)         !> energy grid, E(IE), IE=1,NE
  integer, allocatable :: dist_flag(:) !> distribution flag, dist_flag(IE), IE=1,NE
                                       !> dist_flag(IE) = 1  :: 32 equiprobable bin distribution
                                       !> dist_flag(IE) = -1 :: tabular angular distribution  
                                       !> dist_flag(IE) = 0  :: isotropic
  type (AngularDistDataForm), allocatable :: dist(:) !> angular distribution array, dist(IE), IE=1,NE
end type


type, extends (TabularDataForm) :: EnergyDistDataForm
  integer :: law
  integer :: IDAT
  real(8), allocatable :: LDAT(:)
end type


type EnergyDist
  integer :: nlaw
  type (EnergyDistDataForm), allocatable :: dist(:) !> energy distribution array, dist(I), I=1,nlaw
end type


type, extends (TabularDataForm) :: PrecursorDataForm
  real(8) :: decay_const
end type

type UNRtype 
	integer :: N		!> Number of incident energies where there is a probability table
	integer :: M		!> Length of table; i.e., number of probabilities, typically 20
	integer :: INT 		!> Interpolation parameter between tables =2 lin-lin; =5 log-log
	integer :: ILF		!> Inelastic competition flag
	integer :: IOA		!> Other absorption flag
	integer :: IFF		!> Factors flag
	integer, allocatable :: E(:), P(:,:,:)
endtype 


!Nuclear data library in ace format
type AceFormat
  character(20) :: library       !> name of library for each isotope
  integer :: ZAID                !> ZAID number
  integer :: NXS(1:16)           !> number array in ace format
  integer :: JXS(1:32)           !> pointer array in ace format
  real(8) :: kT                  !> temperature in [MeV]
  real(8) :: atn                 !> ratio of atomic mass to neutron mass

  !Data blocks in ace format
  !> ESZ_Block // FIS_Block
  real(8), allocatable :: E(:)     !> energies, E(I), I=1,NXS(3)
  real(8), allocatable :: sigt(:)  !> total cross sections, E(I), I=1,NXS(3)
  real(8), allocatable :: sigd(:)  !> disappearance cross sections, sigd(I), I=1,NXS(3)
  real(8), allocatable :: sigel(:) !> elastic cross sections, sigel(I), I=1,NXS(3)
  real(8), allocatable :: sigf(:)  !> fission cross sections, sigf(I), I=IE, IE+NE-1
  real(8), allocatable :: H(:)     !> average heating numbers, H(I), I=1,NXS(3)
  !real(8), allocatable :: signuf(:)!> nu*sigf, identical size with sigf(:) // nu is generated in NU block

  real(8) :: qval
  
  !> NU_block
  logical :: nu_block_exist = .true.      !> existence of NU block
  integer :: nu_tot_flag         !> 1 = polynomial function flag, 2 = tabular data flag
  integer :: nu_del_flag         !> 1 = polynomial function flag, 2 = tabular data flag
  type(NuTotDataForm) :: nu_tot 
  type(TabularDataForm) :: nu_del                       !> exist if JXS(24) > 0 
  type(PrecursorDataForm), allocatable :: prcr(:)       !> exist if JXS(24) > 0

  !> MTR Block
  integer, allocatable :: MT(:) !> ENDF MT numbers, MT(I), I=1,NXS(4)


  !> LQR Block
  real(8), allocatable :: Q(:) !> Q-value of reaction MT, Q(I), I=1,NXS(4)


  !> TYR Block
  integer, allocatable :: TY(:) !> Neutron release for reaction MT, TY(I), I=1,NXS(4)


  !> SIG Block
  type (CrossSectionDataForm), allocatable :: sig_MT(:) !> cross section arrays for reaction MT(I), sig_MT(I), I=1,NXS(4)


  !> AND Block
  integer, allocatable :: ang_flag(:)       !> flags for angular distribution arrays, for reaction MT(I), ang_flag(I), I=0,NXS(5)
  type (AngularDist), allocatable :: ang(:) !> angular distribution arrays, for reaction MT(I), ang_flag(I), I=0,NXS(5)


  !> DLW Block // Delayed Neutron Energy Distribution // DLWP Block
  type (EnergyDist), allocatable :: pneg(:)  !> prompt neutron energy distribution arrays, for reaction MT(I), pneg(I), I=1,NXS(5)
  type (EnergyDist), allocatable :: dneg(:)  !> delayed neutron energy distribution arrays, for delayed neutron group I, dneg(I), I=1,NXS(8)
  type (EnergyDist), allocatable :: ppeg(:)  !> prompt photon energy distribution arrays, for reaction MT(I), ppeg(I), I=1,NXS(8)


  !> Energy-Dependent Neutron Yields
  type(TabularDataForm), allocatable :: nyd(:)  !> neutron yield data, nyd(I), I=1,NXS(4) 


  !> UNR Block
  type(UNRtype) :: UNR
  

  integer :: isab
  integer :: iso0K
end type
type (AceFormat), allocatable, target :: ace(:)
integer :: num_iso              !> total number of isotopes


! Hash-based Energy Search algorithm 
real(8) :: Emax, Emin 
integer, allocatable :: ugrid(:,:),&   !Lethargy-grid for hash-based search
					  & ugrid0K(:,:)   !Lethargy-grid for hash-based search
real(8) :: udelta 
integer :: nugrid




!On-the-fly Doppler broadening
logical :: do_OTFDB     !> on-the-fly Doppler broadening in resolved resonance region via Gauss Hermite Quadrature method (Y.G. Jo, KNS 2017)
integer :: scat_kernel  !> scattering kernel : cons = 0, wcm = 1, dbrc = 2
real(8) :: DBRC_min_E,  DBRC_max_E !>[MeV], For U-238, DBRC_min_E = 0.4eV, DBRC_max_E = 210eV
integer :: num_iso0K    !> number of isotopes treated by exact scattering kernel 
type AceFormat0K
  character(20) :: library      !> name of library for each isotope
  integer :: ZAID !> ZAID number
  integer :: NXS(1:16)  !> number array in ace format
  integer :: JXS(1:32)  !> pointer array in ace format
  real(8) :: kT !> temperature in [MeV]
  real(8) :: atn  !> ratio of atomic mass to neutron mass

  !Data block to read
  !> ESZ Block
  real(8), allocatable :: E(:)
  real(8), allocatable :: sigel(:)

  integer :: iso  !> mapping index to original ace library 
end type
type (AceFormat), pointer :: ace0Kptr
type (AceFormat0K), allocatable, target :: ace0K(:)


!S(alpha,beta) thermal scattering treatment
!* Sab tables for different temperatures have same structures
integer :: nsab
logical :: do_sab
type sab_AceFormat
  integer :: iso
  integer :: ntemp
  real(8), allocatable :: kT(:)
  character(20), allocatable :: library(:)
  integer :: NXS(1:16)
  integer :: JXS(1:32)

  !data blocks in ace format S(alfa,beta)


  real(8) :: esab  !> energy upper bound of S(alfa,beta) table
  integer :: ierg_el, ierg_inel, itemp
  real(8) :: inp_el, inp_inel
end type
type (sab_AceFormat), allocatable :: sab(:)


real(8), allocatable :: XSS(:)        !> temporary XSS array for currently reading isotope
real(8), allocatable :: sab_XSS(:,:)  !> temporary XSS array for currently reading isotope for S(alfa,beta)


!Fission energy spectrum
real(8), parameter :: fes(33) = &  
& (/  0.0d0,       .162524d0,  .266043d0,  .358425d0,  .445672d0, &
&      .530293d0,  .613680d0,  .696783d0,  .780264d0,  .864702d0, &
&      .950540d0, 1.038286d0, 1.128326d0, 1.221111d0, 1.317206d0, &
&     1.417070d0, 1.521302d0, 1.630646d0, 1.745929d0, 1.868073d0, &
&     1.998282d0, 2.138046d0, 2.289259d0, 2.454356d0, 2.636707d0, &
&     2.840830d0, 3.073518d0, 3.344965d0, 3.672134d0, 4.086420d0, &
&     4.656234d0, 5.588725d0, 9.0d0  /)


!Fission ZAIDS for fission Q-values.
integer, parameter :: mfiss(22) =  & 
& (/     90232,     91233,     92233,     92234,     92235, &
&        92236,     92237,     92238,     92239,     92240, &
&        93237,     94238,     94239,     94240,     94241, &
&        94242,     94243,     95241,     95242,     95243, &
&        96242,     96244 /)

!Fission Q-values. [MeV]
real(8), parameter :: qfiss(23) =  & 
& (/ 171.91d+0, 175.57d+0, 180.84d+0, 179.45d+0, 180.88d+0, &
&    179.50d+0, 180.40d+0, 181.31d+0, 180.40d+0, 180.40d+0, &
&    183.67d+0, 186.65d+0, 189.44d+0, 186.36d+0, 188.99d+0, &
&    185.98d+0, 187.48d+0, 190.83d+0, 190.54d+0, 190.25d+0, &
&    190.49d+0, 190.49d+0, 180.00d+0 /)

!MT numbers for reactions in ENDF Library
integer, parameter :: &
& ENDF_TOT = 1,   ENDF_ELASTIC = 2,     ENDF_INELASTIC = 4, &
& ENDF_N2N = 16,  ENDF_N3N = 17, &
& ENDF_FISS = 19, ENDF_DISAPPEAR = 101, ENDF_NG = 102,  ENDF_NP = 103, &
& ENDF_NALFA = 107


!Gauss-Hermite Quadratures for on-the-fly Doppler broadening
real(8),parameter :: node_ghq(1:16) = (/-4.688738939305818364688, -3.869447904860122698719, -3.176999161979956026814,&
 -2.546202157847481362159, -1.951787990916253977435, -1.380258539198880796372, -0.8229514491446558925825, &
 -0.2734810461381524521583, 0.2734810461381524521583, 0.8229514491446558925825, 1.380258539198880796372, &
 1.951787990916253977435, 2.546202157847481362159, 3.176999161979956026814, 3.869447904860122698719, &
 4.688738939305818364688/)
real(8),parameter :: weight_ghq(1:16) = (/2.65480747401118224471d-10, 2.32098084486521065339d-7, 2.71186009253788151202d-5, &
 9.32284008624180529914d-4, 0.01288031153550997368346d0, 0.0838100413989858294154d0, 0.2806474585285336753695d0, &
 0.5079294790166137419135d0, 0.5079294790166137419135d0, 0.2806474585285336753695d0, 0.0838100413989858294154d0, &
 0.01288031153550997368346d0, 9.32284008624180529914d-4, 2.71186009253788151202d-5, 2.32098084486521065339d-7, &
 2.65480747401118224471d-10/)
real(8),parameter :: node_ghq2(1:16) = (/1.94840741569E-01, 5.84978765436E-01, 9.76500463590E-01, 1.37037641095E+00, &
 1.76765410946E+00, 2.16949918361E+00, 2.57724953773E+00, 2.99249082500E+00, 3.41716749282E+00, 3.85375548547E+00, &
 4.30554795335E+00, 4.77716450350E+00, 5.27555098652E+00, 5.81222594952E+00, 6.40949814927E+00, 7.12581390983E+00/)
real(8),parameter :: x2_ghq2(1:16) = (/3.7962914575E-02, 3.4220015601E-01, 9.5355315539E-01, 1.8779315077E+00, &
 3.1246010507E+00, 4.7067267077E+00, 6.6422151797E+00, 8.9550013377E+00, 1.1677033674E+01, 1.4851431342E+01, &
 1.8537743179E+01, 2.2821300694E+01, 2.7831438211E+01, 3.3781970488E+01, 4.1081666525E+01, 5.0777223878E+01/)
real(8),parameter :: wx2_ghq2(1:16) = (/1.42451415249E-02, 9.49462195824E-02, 1.44243732244E-01, 1.13536229019E-01, &
 5.48474621706E-02, 1.72025699141E-02, 3.56200987793E-03, 4.85055175195E-04, 4.26280054876E-05, 2.33786448915E-06, &
 7.59830980027E-08, 1.35405428589E-09, 1.17309796257E-11, 4.04486402497E-14, 3.79255121844E-17, 3.71215853650E-21/)

 
 
	contains 
	function find_ACE_iso_idx (this, iso_id) result (idx) 
		type(AceFormat) :: this(:) 
		character(*) :: iso_id 
		integer :: i, idx
		
		do i = 1, size(this) 
			if (trim(this(i)%library) == iso_id) then 
				idx = i 
				return 
			endif
		enddo 
		print *, "no such isotope id : ", iso_id 
		stop 
		
	end function 


end module 