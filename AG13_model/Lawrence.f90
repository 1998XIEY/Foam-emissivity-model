!
! Lawrence Ocean Permittivity module.
!
! Module containing routines to compute the complex permittivities for
! sea water based on
!
!   Lawrence H , Bormann N , English S J .
!   Uncertainties in the permittivity model for seawater in FASTEM and implications for the calibration/validation of microwave imagers[J]. 
!   Journal of Quantitative Spectroscopy and Radiative Transfer, 2019, 243:106813.
!
!
! CREATION HISTORY:
!       Written by:     Lingli He, JCSRA 5-Nov-2020
!                       hllsat@163.com

MODULE Lawrence

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
  PUBLIC :: Lawrence_Ocean_Permittivity
  !PUBLIC :: Lawrence_Ocean_Permittivity_TL
  !PUBLIC :: Lawrence_Ocean_Permittivity_AD
  INTEGER,  PARAMETER :: fp        = SELECTED_REAL_KIND(15)

  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Liu.f90 12622 2011-03-03 23:45:25Z paul.vandelst@noaa.gov $'

  ! Literal constants
  ! -----------------
  REAL(fp), PARAMETER :: ZERO   = 0.0_fp
  REAL(fp), PARAMETER :: ONE    = 1.0_fp
  REAL(fp), PARAMETER :: TWO    = 2.0_fp
  REAL(fp), PARAMETER :: THREE  = 3.0_fp
  REAL(fp), PARAMETER :: FOUR   = 4.0_fp
  
  REAL(fp), PARAMETER :: PI = 3.141592653589793238462643_fp
  REAL(fp), PARAMETER :: K_TO_C = 273.15_fp
  
  ! Permittivity of vacuum (F/m)
  REAL(fp), PARAMETER :: E0 = 8.854187817620389E-012 ! Permittivity of vacuum (F/m)
 

  ! Scaling factors
  ! ---------------
  REAL(fp), PARAMETER :: GHZ_TO_HZ = 1.0e+09_fp ! Gigahertz   -> Hertz

  ! Fixed value for ionic conductivity denominator term, scaled
  ! for frequency. 
  ! -----------------------------------------------------------
  REAL(fp), PARAMETER :: TAU0 = TWO*PI*E0*GHZ_TO_HZ
  
  ! Parameters for the Lawrence et al (2019) permittivity model
  ! ------------------------------------------------------
  ! The coefficients for the high-frequency permittivity temperature
  ! polynomial.
  REAL(fp), PARAMETER :: EINF_COEFF(0:1) = (/ 6.94420520_fp, &
                                              -4.85259069e-02_fp /)
  
  ! The coefficients for the static permittivity temperature
  REAL(fp), PARAMETER :: ES_T_COEFF(0:3) = (/ 87.8010451_fp      , &
                                              -2.68769795e-01_fp, &
                                              -4.74768472e-03_fp, &
                                              5.68884987E-05_fp /)
  REAL(fp), PARAMETER :: ES_S_COEFF(0:2) = (/-3.85142359e-03_fp, &
                                             1.73819707E-05_fp  , &
                                             -8.86453805E-06_fp   /)

  ! The coefficients for the intermediate frequency permittivity
  REAL(fp), PARAMETER :: E1_T_COEFF(0:3) = (/ 11.3724126_fp     , &
                                              -7.02072218e-03_fp     , &
                                              9.45363190e-02_fp, &
                                             -1.34396224e-03_fp /)
  
  ! The coefficients for the relaxation time temperature and
  ! salinity polynomials
  REAL(fp), PARAMETER :: TAU1_T_COEFF(0:3) = (/ 1.11539833e-01_fp , &
                                               -3.13471814e-03_fp, &
                                                8.74997935e-05_fp , &
                                               -1.19366748e-06_fp /)
  REAL(fp), PARAMETER :: TAU1_S_COEFF(0:2) = (/-5.33189243E-04_fp, &  
                                                -1.50409346E-04_fp , &  
                                               7.62214936E-06_fp /)  

  REAL(fp), PARAMETER :: TAU2_T_COEFF(0:3) = (/ 2.42020866e-02_fp, &
                                               1.67067466E-03_fp, &
                                               -4.28919588E-05_fp, &
                                               2.25254269E-07_fp /)
  

  ! The coefficients for the ionic conductivity exponential.
  REAL(fp), PARAMETER :: ALPHA_COEFF = -4.259775841E-08_fp

  ! The coefficients for the ionic conductivity exponent term
  REAL(fp), PARAMETER :: BETA_COEFF(0:5) = (/ 2.59590677E-02_fp, &
                                              -1.44059387E-03_fp, &
                                              5.90688298E-05_fp, &
                                             -1.37008575E-04_fp, &
                                              4.69808667E-05_fp, &
                                             -1.70044096E-06_fp /)
     
  ! The coefficients for the ionic conductivity at 25C polynomial.
  REAL(fp), PARAMETER :: ALPHA25_COEFF(0:3) = (/ 1.65346382e-01_fp, &
                                                -3.58691908E-05_fp, &
                                                 -4.68000709E-06_fp, &
                                                3.89879336E-08_fp /)


  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    REAL(fp) :: t=ZERO                             ! Temperature in deg.C
    REAL(fp) :: s=ZERO                             ! Salinity
    REAL(fp) :: delta=ZERO, beta=ZERO              ! Ionic conductivity components
    REAL(fp) :: alpha25=ZERO, alpha=ZERO           ! Ionic conductivity terms
    REAL(fp) :: es_t=ZERO, es_s=ZERO               ! The temperature and salinity es terms
    REAL(fp) :: e1_t=ZERO, e1_s=ZERO               ! The temperature and salinity e1 terms
    REAL(fp) :: tau1_t=ZERO, tau1_s=ZERO           ! The temperature and salinity tau1 terms
    REAL(fp) :: tau2_t=ZERO, tau2_s=ZERO           ! The temperature and salinity tau2 terms
    REAL(fp) :: f1=ZERO, f2=ZERO                   ! The relaxation compound terms, f.tau
    REAL(fp) :: del1=ZERO, del2=ZERO               ! The permittivity differences
  END TYPE iVar_type


CONTAINS


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Lawrence_Ocean_Permittivity
!
! PURPOSE:
!       Subroutine to compute ocean permittivity according to the reference,
!       Lawrence H , Bormann N , English S J .
!       Uncertainties in the permittivity model for seawater in FASTEM and implications for the calibration/validation of microwave imagers[J]. 
!       Journal of Quantitative Spectroscopy and Radiative Transfer, 2019, 243:106813.
!
! CALLING SEQUENCE:
!       CALL Lawrence_Ocean_Permittivity( Temperature , & ! Input
!                                    Salinity    , & ! Input
!                                    Frequency   , & ! Input
!                                    Permittivity, & ! Output
!                                    iVar          ) ! Internal variable output
!
! INPUTS:
!       Temperature:   Sea surface temperature
!                      UNITS:      Kelvin (K)
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!       Salinity:      Water salinity
!                      UNITS:      ppt (parts per thousand)
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
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE Lawrence_Ocean_Permittivity( &
    Temperature , & ! Input
    Salinity    , & ! Input
    Frequency   , & ! Input
    Permittivity  ) ! Output
    ! Arguments
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Salinity
    REAL(fp),        INTENT(IN)     :: Frequency
    COMPLEX(fp),     INTENT(OUT)    :: Permittivity
    TYPE(iVar_type) :: iVar
    ! Local variables
    REAL(fp) :: einf
    REAL(fp) :: tau1, tau2, es, e1
    REAL(fp) :: re, ie

    ! Setup
    ! -----
    ! ...Initialise imaginary component of result
    ie = ZERO
    ! ...Save the inputs
    iVar%t = Temperature - K_TO_C
    iVar%s = Salinity


    ! Compute the TEMPERATURE polynomial parameterisations
    ! ----------------------------------------------------
    ! ...The high-frequency permittivity temperature polynomial 
    einf = EINF_COEFF(0) + iVar%t*EINF_COEFF(1)
    ! ...The static permittivity temperature polynomial (eqn.4b)
    es = ES_T_COEFF(0) + iVar%t*(ES_T_COEFF(1) + &
                           iVar%t*(ES_T_COEFF(2) + &
                             iVar%t*ES_T_COEFF(3)))
    iVar%es_t = es  ! Save it
    ! ...The intermediate frequency permittivity temperature polynomial (eqn.4c)
    e1 = E1_T_COEFF(0) + iVar%t*(E1_T_COEFF(1) + &
                        iVar%t*(E1_T_COEFF(2) + &
                        iVar%t*E1_T_COEFF(3)))
    iVar%e1_t = e1  ! Save it
    ! ...The Debye relaxation time constants temperature polynomials
    ! ...Units of tau: nanoseconds (for use with GHz frequencies)
    tau1 = TAU1_T_COEFF(0) + iVar%t*(TAU1_T_COEFF(1) + &
                               iVar%t*(TAU1_T_COEFF(2) + &
                                 iVar%t*TAU1_T_COEFF(3)))
    iVar%tau1_t = tau1  ! Save it
    tau2 = TAU2_T_COEFF(0) + iVar%t*(TAU2_T_COEFF(1) + &
                               iVar%t*(TAU2_T_COEFF(2) + &
                                 iVar%t*TAU2_T_COEFF(3)))
    iVar%tau2_t = tau2  ! Save it 


    ! Compute the SALINITY polynomial parameterisations
    ! -------------------------------------------------
    IF ( iVar%s > ZERO ) THEN
      ! ...The temperature difference from 25C used to compute ionic conductivity.
      iVar%delta = 25.0_fp - iVar%t
      ! ...The beta term (eqn.4i) used to compute ionic conductivity
      iVar%beta = BETA_COEFF(0) + iVar%delta*(BETA_COEFF(1) + &
                                    iVar%delta*BETA_COEFF(2)) + &
                  (BETA_COEFF(3) + iVar%delta*(BETA_COEFF(4) + &
                                     iVar%delta*BETA_COEFF(5)))*iVar%S
      ! ...The ionic conductivity at 25C
      iVar%alpha25 = iVar%s*(ALPHA25_COEFF(0) + &
                       iVar%s*(ALPHA25_COEFF(1) + &
                         iVar%s*(ALPHA25_COEFF(2) + &
                           iVar%s*ALPHA25_COEFF(3))))
      ! ...The ionic conductivity 
      iVar%alpha = iVar%alpha25*EXP(-iVar%delta*iVar%beta)
      !  ...The imaginary component dependent on ionic conductivity 
      ie = -iVar%alpha/(Frequency*TAU0)


      ! ...The static permittivity salinity polynomial
      iVar%es_s = ONE + iVar%s*(ES_S_COEFF(0) + iVar%s*ES_S_COEFF(1) + iVar%t*ES_S_COEFF(2))
      es = es * iVar%es_s 


      ! ...The intermediate frequency permittivity have no salinity polynomial 
      iVar%e1_s = ONE
      e1 = e1 * iVar%e1_s 


      ! ...The Debye relaxation time constants salinity polynomials
      ! ...Units of tau: nanoseconds (for use with GHz frequencies)
      iVar%tau1_s = ONE + iVar%s*(TAU1_S_COEFF(0) + iVar%t*(TAU1_S_COEFF(1) + &
                                                      iVar%t*TAU1_S_COEFF(2)))
      tau1 = tau1 * iVar%tau1_s
      
      !tau2 have no salinity polynomial
      iVar%tau2_s = ONE
      tau2 = tau2 * iVar%tau2_s

    END IF
    
    
    ! Compute the complex permittivity
    ! --------------------------------
    ! ...The compound terms
    iVar%f1 = Frequency*tau1  ! Note there is no GHz->Hz conversion.
    iVar%f2 = Frequency*tau2  ! That is embedded in the tau values.
    iVar%del1 = es - e1
    iVar%del2 = e1 - einf
    ! ...The real part
    re = einf + iVar%del1/(ONE + iVar%f1**2) + &
                iVar%del2/(ONE + iVar%f2**2)
    ! ...The imaginary part
    ie = -ie + iVar%del1*iVar%f1/(ONE + iVar%f1**2) + &
               iVar%del2*iVar%f2/(ONE + iVar%f2**2)
    ! ...Combine them (e = re - j.ie)
    Permittivity = CMPLX(re,-ie,fp)                                  

  END SUBROUTINE Lawrence_Ocean_Permittivity
END MODULE Lawrence