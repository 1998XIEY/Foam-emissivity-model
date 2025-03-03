!
!
!
! Run AG13 foam emissivity model
!
program AG131

    use Ocean_Permittivity
    use mod_foam_emiss

    implicit none
    real(8) Temperature,Salinity,Frequency,incident_angle
    real(8) z_max, vf_top,vf_bot

    complex(8) Permittivity
    CHARACTER(LEN=20)  ::  Permittivity_Model
    real(8) e_v,e_h

    integer i


    open(2001,file='/mnt/data/home/xieyc/work/Foam_emissivity_model/results/AG_Z_1.4.txt')

    do i=1,100

        incident_angle=20.0

        !incident_angle=15.0+i*5.0
        !incident_angle=25.0+i*5.0
        Permittivity_Model = TRIM('Liu')

        Frequency=36.5  !GHz
        Temperature = 290 
        Salinity =  25.0
        !z_max=1.099  !cm
        z_max=i*0.1
        vf_top=0.90
        vf_bot=0.01


        ! camps1
        ! Frequency=1.4  !GHz
        ! Temperature = 273.15+20.6
        ! Salinity = 10.49
        ! z_max=1.099  !cm
        ! vf_top=0.2227        
        ! vf_bot=0.01
        
        ! rose1
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
        ! vf_bot=0.005     

        CALL CSEM_Ocean_Permittivity( &
            Temperature,             &
            Salinity,                & 
            Frequency,               &
            Permittivity,            &
            Permittivity_Model) 
        !!!

        CALL esf_anguelova(incident_angle,Frequency,Permittivity,z_max,vf_top,vf_bot,e_v,e_h)

        print*,z_max,e_v,e_h
        write(2001,'(3e13.5)') z_max,e_v,e_h

    enddo

    close(2001)

end program

