MODULE constants_m
    integer    , parameter :: ABOVE=+100000000
    integer    , parameter :: BELOW=-100000000
    real*8     , parameter :: PI=3.141592653589793238462643383279502884197d0 
    real*8     , parameter :: PIO2=1.57079632679489661923132169163975144209858d0 
    real*8     , parameter :: TWOPI=6.283185307179586476925286766559005768394d0 
    real*8     , parameter :: FOURPI=12.56637061435917295385d0
    real*8     , parameter :: PI4 = FOURPI
    real*8     , parameter :: SQRT2=1.41421356237309504880168872420969807856967d0 
    real*8     , parameter :: SQRT2PI=2.50662827463100050242d0
    real*8     , parameter :: EULER=0.5772156649015328606065120900824024310422d0 
    real*8     , parameter :: HALF   = 5.d-1
    real*8     , parameter :: THIRD  = 1.d0/3.d0
    real*8     , parameter :: TWO    = 2.d0
    real*8     , parameter :: THREE  = 3.d0
    real*8     , parameter :: FOUR   = 4.d0
    real*8     , parameter :: FIVE   = 5.d0
    real*8     , parameter :: SIX    = 6.d0
    real*8     , parameter :: TWELVE = 12.d0
    real*8     , parameter :: ThreeQuarters = 3.d0/4.d0
    real*8     , parameter :: a_Bohr = 0.52917720859d0
    real*8     , parameter :: Hartree_2_eV=27.21138386d0        
    real*8     , parameter :: h_bar=6.58264d-4                      ! <== Planck's  constant  (eV * ps)
    real*8     , parameter :: debye_inv  = 2.0819436d-1             ! <== e[C]*d[Angs] = p[Debye] * 0.20819436
    real*8     , parameter :: debye_2_au = 0.393430294d0            ! <== converts Debye units to atomic units (a.u.)
    real*8     , parameter :: low_prec   = 1.1d-7
    real*8     , parameter :: mid_prec   = 1.d-10
    real*8     , parameter :: high_prec  = 1.d-14
    real*8     , parameter :: large      = 1.d+6
    real*8     , parameter :: real_large = 1.d+10
    real*8     , parameter :: infty      = 1.d+30
    real*8     , parameter :: deg_2_rad  = PI / 180.0d0
    real*8     , parameter :: rad_2_deg  = 180.0d0 / PI
    real*8     , parameter :: cutoff_Angs = 12.d0                   ! <== cutoff radius for overlap related calculations  

    complex*16 , parameter :: zi = (0.d0,1.d0)

    complex*16 , parameter :: C_zero = (0.d0,0.d0)
    real*8     , parameter :: D_zero =  0.d0
    integer    , parameter :: I_zero =  0

    complex*16 , parameter :: C_one = (1.d0,0.d0)
    complex*16 , parameter :: C_two = (2.d0,0.d0)
    real*8     , parameter :: D_one =  1.d0
    integer    , parameter :: I_one =  1

!   Molecular dynamics constants ...
    real*8  , parameter :: ee          = 1.60219d-19                ! in Coulomb
    real*8  , parameter :: rsqpi       = 0.56418958354d0            ! sqrt(pi)
    real*8  , parameter :: mol         = 6.02214129d26              ! mol X 1000
    real*8  , parameter :: imol        = 0.166057788d-26            ! 1/(mol x 1000)
    real*8  , parameter :: boltz       = 1.3806488d-23              ! kB in J/K
    real*8  , parameter :: iboltz      = 2.41432176d22              ! 1/(3kB) in K/J
    real*8  , parameter :: coulomb     = 230.7113d0                 ! e^2/(4.pi.epsilon_0) x 1.E+20 N.m^2
    real*8  , parameter :: kJmol_2_eV  = 0.010364272301331d0        ! converts kJ/mol to eV units 
    real*8  , parameter :: eV_2_kJmol  = 96.485307499258d0          ! converts eV to kJ/mol units 
    real*8  , parameter :: cal_2_J     = 4.184d0                    ! converts cal to J units 
    real*8  , parameter :: milli       = 1.0d-3                     ! milli unit
    real*8  , parameter :: micro       = 1.0d-6                     ! micro unit
    real*8  , parameter :: pico_2_sec  = 1.0d-12                    ! converts picosecond units to second units
    real*8  , parameter :: nano_2_angs = 1.0d1                      ! converts nanometer units to angstron units
    real*8  , parameter :: mts_2_nano  = 1.0d9                      ! converts meter units to nanometer units
    real*8  , parameter :: Angs_2_mts  = 1.0d-10                    ! converts Angs units to meter units
    real*8  , parameter :: mts_2_Angs  = 1.0d+10                    ! converts meter units to Angs units
    real*8  , parameter :: cm_2_Angs   = 1.0d+8                     ! converts centimeter units to Angs units
    real*8  , parameter :: amber_LJ    = 1.12246204830937d0         ! converts sigma LJ (amber) to OPLS 
    real*8  , parameter :: factor1     = 1.0d26                     ! factor used to not work with small numbers
    real*8  , parameter :: factor2     = 1.0d24                     ! factor used to not work with small numbers
    real*8  , parameter :: factor3     = 1.0d-20                    ! factor used to compensate the factor1 and factor2 factors

END MODULE constants_m
