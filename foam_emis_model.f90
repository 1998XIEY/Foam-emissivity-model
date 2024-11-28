Module foam_emis_model
  use Ocean_Permittivity

implicit none

contains

subroutine foam_emis_xie(theta,frequency,temperature,Salinity,z_max,vf_top,e_v,e_h)

  !theta,入射角度, degree
  !Frequency, 频率, GHz
  !Temperature, 温度, K
  !Salinity, 盐度, PSU
  !z_max, foam thickness, cm
  !vf_top, void fraction at the atmosphere-foam interface

  !---海沫光学厚度定义
  integer,parameter :: n_layer=10
  real(8) vf_top,vf_bot,a,b,m,z_max
  real(8) vf(n_layer),z(n_layer),z_step
  real(8) Temperature,Salinity,Frequency,k0,incident_angle
  real(8) temp(n_layer)  !温度
  real(8) ke(n_layer),albedo(n_layer)
  complex(8) Permittivity,foam_permittivity(n_layer),sea_permittivity,eair
  CHARACTER(LEN=20)  ::  Permittivity_Model

  real(8) alpha(n_layer),beta(n_layer),p(n_layer),q(n_layer)
  real(8) sigma(n_layer),t1(n_layer),tau(n_layer)

  real(8) e_v,e_h
  real(8) theta,theta_i,theta_t,mu,pi
  real(8) ssalb_h,ssalb_v,gh,gv
  real(8) r21_v, r21_h,t21_v, t21_h,r23_v, r23_h
  integer i

  pi=3.1415926536
  theta= theta*pi/180.0
  Permittivity_Model = TRIM('Liu')

  vf_bot=0.01
  m=0.1
  a=vf_top+m
  b=log((a-vf_bot)/m)/z_max  !z_max用于归一化
  z_step=z_max/(n_layer-1)
  do i=1,n_layer
     z(i)=0.0+z_step*(i-1)
  enddo
  vf(:)=a-m*exp(b*z(:))
  vf(n_layer)=vf_bot

  CALL CSEM_Ocean_Permittivity( &
     Temperature,             &
     Salinity,                & 
     Frequency,               &
     Permittivity,            &
     Permittivity_Model) 

  sea_permittivity=Permittivity

  !
  Permittivity=Permittivity**(1.0/2.0)
  foam_permittivity(:)=(vf+(1-vf)*Permittivity)**2.0

  !
  k0=Frequency/30.0   !1/cm
  alpha=k0*abs(aimag(sqrt(foam_permittivity)))
  beta=k0*real(sqrt(foam_permittivity))
  p=2.0*alpha*beta
  q=beta**2-alpha**2-k0**2*sin(theta)**2
  t1=sqrt(p**2+q**2)+q
  sigma=atan(sqrt(2.0)*k0*sin(theta)/sqrt(t1))
  tau=2.0*alpha/cos(sigma)*z

  ssalb_h=0.05
  ssalb_v=0.05
  gh=0.0   ! 非对称因子，大于0为前向散射，等于0为各项同性散射
  gv=0.0
  ! 进行二流近似辐射传输
  eair = CMPLX(1.0,-0.0)
  mu      = COS(theta)
  theta_i = ASIN(REAL(SIN(theta)*SQRT(eair)/SQRT(foam_permittivity(1))))
  CALL Reflectance(foam_permittivity(1), eair, theta_i,  theta, r21_v, r21_h)
  CALL Transmittance(foam_permittivity(1), eair, theta_i, theta, t21_v, t21_h)
  theta_t = ASIN(REAL(SIN(theta_i)*SQRT(foam_permittivity(1))/SQRT(sea_permittivity)))
  CALL Reflectance(foam_permittivity(1), sea_permittivity, theta_i, theta_t, r23_v, r23_h)
  CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau(n_layer),tau(n_layer), &
                          r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,e_v,e_h, &
                          frequency, Temperature, Temperature)
  
  theta=theta/pi*180

end subroutine



subroutine Reflectance(em1, em2, theta_i, theta_t, rv, rh)
    REAL(8) :: theta_i, theta_t
    REAL(8) :: rh, rv,cos_i,cos_t
    COMPLEX(8) :: em1, em2, m1, m2, angle_i, angle_t
    ! compute the refractive index ratio between medium 2 and 1
    ! using dielectric constant (n = SQRT(e))
    cos_i = COS(theta_i)
    cos_t = COS(theta_t)
  
    angle_i = CMPLX(cos_i, 0.0, 8)
    angle_t = CMPLX(cos_t, 0.0, 8)
  
    m1 = SQRT(em1)
    m2 = SQRT(em2)
  
    rv = (ABS((m1*angle_t-m2*angle_i)/(m1*angle_t+m2*angle_i)))**2
    rh = (ABS((m1*angle_i-m2*angle_t)/(m1*angle_i+m2*angle_t)))**2
  
  end subroutine Reflectance

  
subroutine Transmittance(em1,em2,theta_i,theta_t,tv,th)
  REAL(8) :: theta_i, theta_t
  REAL(8) :: th, tv, rr, cos_i,cos_t
  COMPLEX(8) :: em1, em2, m1, m2, angle_i, angle_t
  ! compute the refractive index ratio between medium 2 and 1
  ! using dielectric constant (n = SQRT(e))
  cos_i = COS(theta_i)
  cos_t = COS(theta_t)
  angle_i = CMPLX(cos_i, 0.0, 8)
  angle_t = CMPLX(cos_t, 0.0, 8)
  m1 = SQRT(em1)
  m2 = SQRT(em2)  
  rr = ABS(m2/m1)*cos_t/cos_i
  tv = rr*(ABS(2.0*m1*angle_i/(m1*angle_t + m2*angle_i)))**2
  th = rr*(ABS(2.0*m1*angle_i/(m1*angle_i + m2*angle_t)))**2
end subroutine Transmittance


!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Two_Stream_Solution
! PURPOSE:
!       Two stream solution 
!       Updated with the more accurate formula of total upwelling radiance emanating from the surface.
! INPUT:
!      b:              Scattering layer temperature (k)         (gdas)   (not used here)
!      mu:             cos(theta)
!      gv:             Asymmetry factor for v pol
!      gh:             Asymmetry factor for h pol
!      ssalb_v:        Single scattering albedo at v. polarization
!      ssalb_h:        Single scattering albedo at h. polarization
!      tau_v:          Optical depth at v. polarization
!      tau_h:          Optical depth at h. polarization
!      r12_v:          Reflectivity at vertical polarization   (not used here)
!      r12_h:          Reflectivity at horizontal polarization (not used here)
!      r21_v:          Reflectivity at vertical polarization
!      r21_h:          Reflectivity at horizontal polarization
!      r23_v:          Reflectivity at vertical polarization
!      r23_h:          Reflectivity at horizontal polarization
!      t21_v:          Transmisivity at vertical polarization
!      t21_h:          Transmisivity at horizontal polarization
!      t12_v:          Transmisivity at vertical polarization   (not used here)
!      t12_h:          Transmisivity at horizontal polarization (not used here)
!      Frequency:      Frequency
!      t_soil:         Soil temperature
!      t_skin:         Land surface temperature
!
! OUTPUT:
!      esv:             Emissivity at vertical polarization
!      esh:             Emissivity at horizontal polarization
! REFERENCES:
!    Weng, F., B. Yan, and N. Grody, 2001: "A microwave land emissivity model",
!     J. Geophys. Res., 106, 20, 115-20, 123
!
!----------------------------------------------------------------------------------

subroutine Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
      r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,esv,esh,frequency,t_soil,t_skin)

  REAL(8) :: mu, gv, gh, ssalb_h, ssalb_v, tau_h,tau_v,                 &
              r21_h, r21_v, r23_h, r23_v, t21_v, t21_h, esv, esh
  REAL(8) :: alfa_v, alfa_h, kk_h, kk_v, gamma_h, gamma_v, beta_v, beta_h
  REAL(8) :: fact1,fact2
  REAL(8) :: frequency, t_soil, t_skin
  REAL(8) :: gsect0, gsect1_h, gsect1_v, gsect2_h, gsect2_v

  REAL(8) :: C_2

  C_2 = 6.62606876e-34 * 2.99792458e+08 /(8.314472/6.02214199e+23)

  alfa_h  = SQRT((1.0 - ssalb_h)/(1.0 - gh*ssalb_h))
  kk_h    = SQRT((1.0 - ssalb_h)*(1.0 -  gh*ssalb_h))/mu
  beta_h  = (1.0 - alfa_h)/(1.0 + alfa_h)
  gamma_h = (beta_h -r23_h)/(1.0-beta_h*r23_h)

  alfa_v  = SQRT((1.0-ssalb_v)/(1.0 - gv*ssalb_v))
  kk_v    = SQRT((1.0-ssalb_v)*(1.0 - gv*ssalb_v))/mu
  beta_v  = (1.0 - alfa_v)/(1.0 + alfa_v)
  gamma_v = (beta_v -r23_v)/(1.0-beta_v*r23_v)


  fact1=gamma_h*EXP(-2.0*kk_h*tau_h)
  fact2=gamma_v*EXP(-2.0*kk_v*tau_v)

  gsect0  =(EXP(C_2*frequency/t_skin) -1.0)/(EXP(C_2*frequency/t_soil) -1.0)

  
  gsect1_h=(1.0-r23_h)*(gsect0-1.0)
  gsect2_h=((1.0-beta_h*beta_h)/(1.0-beta_h*r23_h))*EXP(-kk_h*tau_h)
  gsect1_v=(1.0-r23_v)*(gsect0-1.0)
  gsect2_v=((1.0-beta_v*beta_v)/(1.0-beta_v*r23_v))*EXP(-kk_h*tau_v)

  esh  = t21_h*((1.0 - beta_h)*(1.0 + fact1)+gsect1_h*gsect2_h) /(1.0-beta_h*r21_h-(beta_h-r21_h)*fact1)
  esv  = t21_v*((1.0 - beta_v)*(1.0 + fact2)+gsect1_v*gsect2_v) /(1.0-beta_v*r21_v-(beta_v-r21_v)*fact2)
  if (esh > 1.0) esh = 1.0
  if (esv > 1.0) esv = 1.0

end subroutine Two_Stream_Solution


end module