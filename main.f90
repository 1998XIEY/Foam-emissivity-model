
program foam_emis
    use foam_emis_model
    implicit none
    real(8) theta,e_v,e_h
    real(8) vf_top,z_max
    real(8) Temperature,Salinity,Frequency
    integer i_freq,i_angle,i_z,i_vf

    open(2001,file='TS_test.txt')
 
    do i_freq=1,1
       do i_angle=1,1     !7
          do i_z=1,100  !

            !theta=15.0+i_angle*5.0 
            theta=20.0
            Frequency=1.4
            Temperature = 290.0
            Salinity = 35.0
            z_max=i_z*0.1
            vf_top=0.90  

            call foam_emis_xie(theta,Frequency,Temperature,Salinity,z_max,vf_top,e_v,e_h)

            !  print*,Frequency,vf_top,e_v,e_h
            !  write(2001,'(4e13.5)') Frequency,vf_top,e_v,e_h
            print*,z_max,e_v,e_h
            write(2001,'(3e13.5)') z_max,e_v,e_h             

          enddo
       enddo
    enddo
    close(2001)
 
 end program foam_emis
 
