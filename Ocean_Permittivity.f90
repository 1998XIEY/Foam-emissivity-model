!
! Ocean_Permittivity
!
! Container for modules to compute complex permittivities for
! sea water.
!
! Three models are included in this module; that of
!
!   Guillou, C. et al. (1998) Impact of new permittivity measurements
!      on sea surface emissivity modeling in microwaves.
!      Radio Science, Volume 33, Number 3, Pages 649-667
!
! and of
!
!   Ellison, W.J. et al. (2003) A comparison of ocean emissivity models
!     using the Advanced Microwave Sounding Unit, the Special Sensor
!     Microwave Imager, the TRMM Microwave Imager, and airborne radiometer
!     observations. Journal of Geophysical Research, v108, D21, Pages ACL 1,1-14
!     doi:10.1029/2002JD0032132
!
! and of
!
!   Liu, Q. et al. (2010) An improved fast microwave water emissivity model.
!      IEEE Trans. Geosci. Remote Sensing, accepted June 25, 2010
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 11-Apr-2007
!                       paul.vandelst@noaa.gov
!      Modified by:     ming Chen, 01-Apr-2016
!                       ming.chen@noaa.gov
!

MODULE Ocean_Permittivity

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
!  USE CSEM_Type_Kinds, ONLY: fp
  
  USE Guillou, ONLY: GuillouVar_type => iVar_type , &
                     Guillou_Ocean_Permittivity
  USE Ellison, ONLY: EllisonVar_type => iVar_type , &
                     Ellison_Ocean_Permittivity
  USE Liu    , ONLY: LiuVar_type => iVar_type , &
                     Liu_Ocean_Permittivity
  
  USE Lawrence    , ONLY: LawrenceVar_type => iVar_type , &
                     Lawrence_Ocean_Permittivity
  
  USE Meissner, ONLY:Meissner_Ocean_Permittivity
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! For Guillou model
  PUBLIC :: GuillouVar_type
  PUBLIC :: Guillou_Ocean_Permittivity
  ! For Ellison model
  PUBLIC :: EllisonVar_type
  PUBLIC :: Ellison_Ocean_Permittivity
  ! For Liu model
  PUBLIC :: LiuVar_type
  PUBLIC :: Liu_Ocean_Permittivity
  ! For Lawrence model
  PUBLIC :: LawrenceVar_type
  PUBLIC :: Lawrence_Ocean_Permittivity
  ! For Meissner and Wentz model
  PUBLIC :: Meissner_Ocean_Permittivity
  
  PUBLIC :: CSEM_Ocean_Permittivity
  
  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Ocean_Permittivity.f90 12622 2011-03-03 23:45:25Z paul.vandelst@noaa.gov $'
  INTEGER,  PARAMETER :: fp        = SELECTED_REAL_KIND(15)
  
CONTAINS
!===========================================================================
! This subroutine provide a general interface to the optional  permittivity 
! models of ocean. The default is set the one by A. Klein and Calvin Swift 
! was originally implemented in SPM model by Simon H. Yueh 
! Double Debye models by Ellison and Liu are the same as those currently
! implemented in NOAA FASTEM
!
! Input parameters
!	Frequency:   frequency in GHz
!	Temperature: Temperature K
! 	Salinity:    Salinity in parts per thousand
! Output parameters
!	Permittivity: complex permittivity
!
! Author:Ming Chen, Dec-15-2015
!       ming.chen@noaa.gov
!=============================================================================

  SUBROUTINE CSEM_Ocean_Permittivity(&
    Temperature , & ! Input
    Salinity    , & ! Input
    Frequency   , & ! Input
    Permittivity, & ! Output
    Alg_Name      )  
  
    ! Arguments
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Salinity
    REAL(fp),        INTENT(IN)     :: Frequency
    COMPLEX(fp),     INTENT(OUT)    :: Permittivity
    CHARACTER(*), OPTIONAL ::  Alg_Name 
    
!    print*,TRIM(Alg_Name)
    SELECT CASE (TRIM(Alg_Name))
      CASE('Liu')
        CALL Liu_Ocean_Permittivity(     &
             Temperature,  & ! Input
             Salinity,     & ! Input
             Frequency,    & ! Input
             Permittivity  )
     
      CASE('Guillou')
        CALL Guillou_Ocean_Permittivity( &
             Temperature,  & ! Input
             Salinity,     & ! Input
             Frequency,    & ! Input
             Permittivity  )
      
      CASE('Ellison')
        CALL Ellison_Ocean_Permittivity( &
             Temperature,  & ! Input
             Frequency,    & ! Input
             Permittivity  ) 
        
      CASE('Lawrence')
       CALL Lawrence_Ocean_Permittivity(     &
             Temperature,  & ! Input
             Salinity,     & ! Input
             Frequency,    & ! Input
             Permittivity  )
      CASE('Meissner')
       CALL Meissner_Ocean_Permittivity(     &
             Temperature,  & ! Input
             Salinity,     & ! Input
             Frequency,    & ! Input
             Permittivity  ) 
      CASE Default
        CALL Klein_Ocean_Permittivity(   &
             Temperature,  & ! Input
             Salinity,     & ! Input
             Frequency,    & ! Input
             Permittivity  ) 
    END SELECT

    Permittivity = CONJG(Permittivity)

  END SUBROUTINE CSEM_Ocean_Permittivity

!===========================================================================
! This subroutine calculates the permittivity
! of ocean. The formula is based on the Debye equation
! given in
!`An Improved model for the Dielectric Constant
! of Sea Water at Microwave Frequencies,''
! Lawrence A. Klein and Calvin Swift,
!       IEEE Trans. AP-25, No. 1, 104-111, 1977.
!
! Input parameters
!	Frequency:   frequency in GHz
!	Temperature: Temperature K
! 	Salinity:    Salinity in parts per thousand
! Output parameters
!	Permittivity: complex permittivity
!
! Author:Ming Chen, Dec-15-2015
!       ming.chen@noaa.gov
!=============================================================================
 SUBROUTINE Klein_Ocean_Permittivity( &
    Temperature , & ! Input
    Salinity    , & ! Input
    Frequency   , & ! Input
    Permittivity  ) ! Output
    ! Arguments
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Salinity
    REAL(fp),        INTENT(IN)     :: Frequency
    COMPLEX(fp),     INTENT(OUT)    :: Permittivity

    complex(fp) :: cj=cmplx(0.0,1.0)

    REAL(fp) :: T,S,f
    REAL(fp) :: a,b,alpha,beta,delta,omega,epsinf, eps0
    REAL(fp) :: tau, tau0, es, es0, sigma, sigma25
  
    f = Frequency
    T = Temperature - 273.15
    s = Salinity
    
    alpha=0
    es0=87.134-0.1949*T-1.276E-2*T**2+2.491E-4*T**3
    a=1+1.613E-5*S*T-3.656E-3*S+3.210E-5*S**2-4.232E-7*S**3
    es=es0*a

    tau0=1.768E-11-6.086E-13*T+1.104E-14*T**2-8.111E-17*T**3
    b=1+2.282E-5*S*T-7.638E-4*S-7.760E-6*S**2+1.105E-8*S**3
    tau=tau0*b

    sigma25 = S*(0.182521-1.46192E-3*S+2.09324E-5*S**2-1.28205E-7*S**3)
    delta = 25-T
    beta = 2.033E-2+1.266E-4*delta+2.464E-6*delta**2 &
          -S*(1.849E-5-2.551E-7*delta+2.551E-8*delta**2)
    sigma = sigma25*exp(-delta*beta)

    epsinf=4.9

    omega=f*2*3.1415926535*1.E9

    eps0=8.854E-12

    ! Debye Equation
    Permittivity=epsinf+(es-epsinf)/(1-(cj*omega*tau)**(1-alpha)) &
         +cj*sigma/(omega*eps0)

    Permittivity = CONJG(Permittivity)
    RETURN
  END SUBROUTINE Klein_Ocean_Permittivity
  
  

END MODULE Ocean_Permittivity
