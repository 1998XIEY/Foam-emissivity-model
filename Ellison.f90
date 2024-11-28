!
! Ellison Ocean Permittivity module.
!
! Module containing routines to compute the complex permittivities for
! sea water based on
!
!   Ellison, W.J. et al. (2003) A comparison of ocean emissivity models
!     using the Advanced Microwave Sounding Unit, the Special Sensor
!     Microwave Imager, the TRMM Microwave Imager, and airborne radiometer
!     observations. Journal of Geophysical Research, v108, D21, Pages ACL 1,1-14
!     doi:10.1029/2002JD0032132
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 11-Apr-2007
!                       paul.vandelst@noaa.gov
!

MODULE Ellison

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
!  USE CSEM_Type_Kinds, ONLY: fp

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! ... Datatypes
  PUBLIC :: iVar_type
  ! ... Procedures
  PUBLIC :: Ellison_Ocean_Permittivity
  INTEGER,  PARAMETER :: fp        = SELECTED_REAL_KIND(15)


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Ellison.f90 17503 2012-01-31 22:54:09Z paul.vandelst@noaa.gov $'
  REAL(fp), PARAMETER :: PI = 3.141592653589793238462643_fp  
  REAL(fp), PARAMETER :: ZERO   = 0.0_fp
  REAL(fp), PARAMETER :: POINT5 = 0.5_fp
  REAL(fp), PARAMETER :: ONE    = 1.0_fp
  REAL(fp), PARAMETER :: TWO    = 2.0_fp
  REAL(fp), PARAMETER :: THREE  = 3.0_fp
  REAL(fp), PARAMETER :: FOUR   = 4.0_fp
  REAL(fp), PARAMETER :: FIVE   = 5.0_fp
  REAL(fp), PARAMETER :: TWOPI  = TWO*PI

  REAL(fp), PARAMETER :: K_TO_C = 273.15_fp
  
  ! Permittivity of vacuum (F/m)
  REAL(fp), PARAMETER :: E0 = 8.854187817620389E-012 ! Permittivity of vacuum (F/m)

  ! Scaling factors (used here for documenting the conversion
  ! to SI units of the double Debye model denominator)
  ! ---------------------------------------------------------
  REAL(fp), PARAMETER :: PS_TO_S   = 1.0e-12_fp ! Picoseconds -> Seconds
  REAL(fp), PARAMETER :: GHZ_TO_HZ = 1.0e+09_fp ! Gigahertz   -> Hertz
  REAL(fp), PARAMETER :: SCALE_FACTOR = PS_TO_S * GHZ_TO_HZ


  ! Parameters for the Ellison et al (2003) permittivity model
  ! ----------------------------------------------------------
  ! The coefficients used to fit the Double Debye model
  REAL(fp), PARAMETER :: TAU1_COEFF(0:2)   = (/    17.535_fp, &
                                                 -0.61767_fp, &
                                                0.0089481_fp /)
  REAL(fp), PARAMETER :: TAU2_COEFF(0:3)   = (/    3.1842_fp, &
                                                 0.019189_fp, &
                                                -0.010873_fp, &
                                               0.00025818_fp /)
  REAL(fp), PARAMETER :: DELTA1_COEFF(0:3) = (/    68.396_fp, &
                                                 -0.40643_fp, &
                                                 0.022832_fp, &
                                              -0.00053061_fp /)
  REAL(fp), PARAMETER :: DELTA2_COEFF(0:3) = (/    4.7629_fp, &
                                                   0.1541_fp, &
                                                -0.033717_fp, &
                                               0.00084428_fp /)
  REAL(fp), PARAMETER :: EINF_COEFF(0:1)   = (/   5.31250_fp, &
                                               -0.0114770_fp /)
  REAL(fp), PARAMETER :: SIGMA_COEFF(0:1) = (/      2.906_fp, &
                                                  0.09437_fp /)


  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    REAL(fp) :: t=ZERO                    ! Temperature in degC
    REAL(fp) :: f=ZERO, f2=ZERO, f0=ZERO  ! Frequency terms
    REAL(fp) :: tau1=ZERO  , tau2=ZERO    ! Relaxation frequencies
    REAL(fp) :: delta1=ZERO, delta2=ZERO  ! Delta terms
    REAL(fp) :: d1=ZERO    , d2=ZERO      ! Denominator terms
  END TYPE iVar_type


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Ellison_Ocean_Permittivity
!
! PURPOSE:
!       Subroutine to compute ocean permittivity according to the reference,
!         Ellison, W.J. et al. (2003) A comparison of ocean emissivity models
!           using the Advanced Microwave Sounding Unit, the Special Sensor
!           Microwave Imager, the TRMM Microwave Imager, and airborne radiometer
!           observations. Journal of Geophysical Research, v108, D21, ACL 1,1-14
!           doi:10.1029/2002JD0032132
!
! CALLING SEQUENCE:
!       CALL Ellison_Ocean_Permittivity( Temperature , & ! Input
!                                        Frequency   , & ! Input
!                                        Permittivity, & ! Output
!                                        iVar          ) ! Internal variable output
!
! INPUTS:
!       Temperature:   Sea surface temperature
!                      UNITS:      Kelvin (K)
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!       Frequency:     Frequency
!                      UNITS:      GHz
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Permittivity:  Ocean permittivity
!                      UNITS:      N/A
!                      TYPE:       COMPLEX(fp)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!       iVar:          Structure containing internal variables required for
!                      subsequent tangent-linear or adjoint model calls.
!                      The contents of this structure are NOT accessible
!                      outside of this module.
!                      UNITS:      N/A
!                      TYPE:       TYPE(iVar_type)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
! COMMENTS:
!       There is currently no salinity dependence.
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE Ellison_Ocean_Permittivity( &
    Temperature , & ! Input
    Frequency   , & ! Input
    Permittivity  ) ! Output

    ! Arguments
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Frequency
    COMPLEX(fp),     INTENT(OUT)    :: Permittivity
    TYPE(iVar_type)  :: iVar
    ! Local variables
    REAL(fp) :: einf
    REAL(fp) :: re1, re2
    REAL(fp) :: ie1, ie2
    REAL(fp) :: sigma, iesigma
    REAL(fp) :: re, ie


    ! Compute the various polynomial components of the double Debye model
    ! -------------------------------------------------------------------
    ! Compute temperature value for polynomials
    iVar%t = Temperature - K_TO_C

    ! Compute the Debye model relaxation frequencies
    ! (eqn on pg ACL 1-4 of Ellison et al. 2003)
    iVar%tau1 = TAU1_COEFF(0) + iVar%t*(TAU1_COEFF(1) + &
                                  iVar%t*TAU1_COEFF(2))
    iVar%tau2 = TAU2_COEFF(0) + iVar%t*(TAU2_COEFF(1)   + &
                                  iVar%t*(TAU2_COEFF(2) + &
                                    iVar%t*TAU2_COEFF(3)))

    ! Compute the delta terms
    ! (eqn on pg ACL 1-4 of Ellison et al. 2003)
    iVar%delta1 = DELTA1_COEFF(0) + iVar%t*(DELTA1_COEFF(1) + &
                                      iVar%t*(DELTA1_COEFF(2) + &
                                        iVar%t*DELTA1_COEFF(3)))
    iVar%delta2 = DELTA2_COEFF(0) + iVar%t*(DELTA2_COEFF(1) + &
                                      iVar%t*(DELTA2_COEFF(2) + &
                                        iVar%t*DELTA2_COEFF(3)))

    ! Compute the "infinite" permittivity term
    ! (No coeffs provided in ref. Taken from existing code)
    einf = EINF_COEFF(0) + iVar%t*EINF_COEFF(1)


    ! Compute the permittivities using the double Debye model
    ! (eqn on pg ACL 1-3 of Ellison et al. 2003)
    ! -------------------------------------------------------
    ! The common frequency terms
    iVar%f  = TWOPI * Frequency * SCALE_FACTOR
    iVar%f2 = iVar%f**2
    iVar%f0 = TWOPI * Frequency * GHZ_TO_HZ * E0

    ! The denominators of the double Debye model
    iVar%d1 = ONE + iVar%f2*iVar%tau1**2
    iVar%d2 = ONE + iVar%f2*iVar%tau2**2

    ! The real parts of the "delta" terms
    re1 = iVar%delta1 / iVar%d1
    re2 = iVar%delta2 / iVar%d2

    ! The imaginary parts of the "delta" terms
    ie1 = iVar%delta1 * iVar%f * iVar%tau1 / iVar%d1
    ie2 = iVar%delta2 * iVar%f * iVar%tau2 / iVar%d2

    ! The conductivity term
    sigma   = SIGMA_COEFF(0) + SIGMA_COEFF(1)*iVar%t
    iesigma = sigma / iVar%f0

    ! Construct the complex permittivity, e = e' - j.e"
    re = re1 + re2 + einf
    ie = ie1 + ie2 + iesigma
    Permittivity = CMPLX(re, -ie, fp)

  END SUBROUTINE Ellison_Ocean_Permittivity


END MODULE Ellison
