
program foam_emis
   use Ocean_Permittivity
   use two_streams
   implicit none


    !---海沫光学厚度定义
   integer,parameter :: n_layer=10
   real(8) vf_top,vf_bot,a,b,m,z_max
   real(8) vf(n_layer),z(n_layer),z_step
   real(8) Temperature,Salinity,Frequency,k0,incident_angle
   real(8) temp(n_layer)  !温度
   real(8) ke(n_layer),albedo(n_layer)
   complex(8) Permittivity,foam_permittivity(n_layer),ocean_permittivity,eair
   CHARACTER(LEN=20)  ::  Permittivity_Model

   real(8) alpha(n_layer),beta(n_layer),p(n_layer),q(n_layer)
   real(8) sigma(n_layer),t1(n_layer),tau(n_layer)

   !
   real(8) theta,theta_i,theta_t,mu,pi
   real(8) ssalb_h,ssalb_v,gh,gv
   real(8) r21_v, r21_h,t21_v, t21_h,r23_v, r23_h,e_v,e_h

   real(8) e_v0,e_h0
   integer count


   integer i_freq,i,i_angle,i_z

   pi=3.1415926536
   e_v0=0.0
   e_h0=0.0


   open(2001,file='/mnt/data/home/xieyc/work/Two-streams/results/TS_Z_saturation.txt')

   do i_freq=1,37
      do i_angle=1,80     !7
         count=0
         do i_z=1,200  !
            Permittivity_Model = TRIM('Liu')
            
            theta=i_angle
            !theta = 15.0+i_angle*5.0
            !theta = 25.0+i_angle*5.0
            theta= theta*pi/180.0

            Frequency=i_freq
            Temperature = 290 
            Salinity =  25.0
            !z_max=1.099  !cm
            z_max=i_z*0.1
            vf_top=0.90
            vf_bot=0.01

            !-----------------------------------------------------------------
            ! Camps 1
            ! Frequency=1.4  !GHz
            ! Temperature = 273.15+20.6
            ! Salinity = 10.49
            ! z_max=1.099  !cm
            ! vf_top=0.2227        
            ! vf_bot=0.01

            
            ! Rose 1
            ! Frequency=10.8 !GHz
            ! Temperature = 273.15+19.0
            ! Salinity = 10.0
            ! z_max=2.8  !cm
            ! vf_top=0.93      
            ! vf_bot=0.01     
            
            ! rose2
            ! Frequency=36.5 !GHz
            ! Temperature = 273.15+19.0
            ! Salinity = 10.0
            ! z_max=2.8  !cm
            ! vf_top=0.91      
            ! vf_bot=0.01  
            !-------------------------------------------------------------
            


            m=0.1
            a=vf_top+m
            b=log((a-vf_bot)/m)/z_max  !z_max用于归一化
            z_step=z_max/(n_layer-1)
            do i=1,n_layer
               z(i)=0.0+z_step*(i-1)
            enddo
            vf(:)=a-m*exp(b*z(:))
            vf(n_layer)=vf_bot

            !-------------------------------------------------
            !接入一个海洋介电常数模型
            CALL CSEM_Ocean_Permittivity( &
               Temperature,             &
               Salinity,                & 
               Frequency,               &
               Permittivity,            &
               Permittivity_Model) 

            ocean_permittivity=Permittivity
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
            theta_t = ASIN(REAL(SIN(theta_i)*SQRT(foam_permittivity(1))/SQRT(ocean_permittivity)))
            CALL Reflectance(foam_permittivity(1), ocean_permittivity, theta_i, theta_t, r23_v, r23_h)
            CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau(n_layer),tau(n_layer), &
                                    r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,e_v,e_h, &
                                    frequency, Temperature, Temperature)

                                 
            ! 计算饱和厚度
            if (abs(e_v0-e_v)<0.00005) then
               if (count==0) then
                  print*,z_max,e_v,e_h
                  write(2001,'(5e13.5)') Frequency,z_max,e_v,e_h,theta/pi*180
               endif
               count=1
            endif
            e_v0=e_v
            !----------------------结束



            ! print*,theta/pi*180,e_v,e_h
            ! write(2001,'(3e13.5)') theta/pi*180,e_v,e_h

            ! print*,z_max,e_v,e_h
            ! write(2001,'(3e13.5)') z_max,e_v,e_h

            !write(2001,'(3e13.5)') Frequency,e_v,e_h
            !输出foam 介电常数,100层
            ! do i=1,100
            !    write(2001,'(3e13.5)') vf(i),real(foam_permittivity(i)),aimag(foam_permittivity(i))
            ! enddo

         enddo
      enddo
   enddo
   close(2001)

end program foam_emis
