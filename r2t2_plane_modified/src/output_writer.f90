module output
!!RT-Engine
!!Copyright (C) 2016 Karri Muinonen, Timo Väisänen and University of Helsinki
!!All rights reserved.
!!The new BSD License is applied to this software, see LICENSE.txt
    use constants
    use typedefinitions
    use scatterer
    use geometry
    use grid_holder
    implicit none

    character(len=*), parameter  :: fmt1 = "(A20,3X,I2,'=',E12.5)"

    private fmt1
contains

    subroutine print_arr_with_message(arr,r_arr,stream,msg,fmt0)
        integer, intent(in) :: stream
        character(*), intent(in) :: msg
        character(*), intent(in) :: fmt0
        real(kind=rk), intent(in) :: arr(:)
        integer, intent(in) :: r_arr(:)
        integer :: j1
        do j1=1,size(arr)
            write(stream,*) msg,j1,r_arr(j1)
        enddo
    end subroutine

    subroutine print_single_val_with_message(arr,r_arr,stream,msg,fmt0,rays)
        integer, intent(in) :: stream
        character(*), intent(in) :: msg
        character(*), intent(in) :: fmt0
        real(kind=rk), intent(in) :: arr(:)
        integer, intent(in) :: r_arr(:)
        integer, intent(in) :: rays(:)
        write(stream,*) msg,sum(arr(:))/(sum(rays))
    end subroutine

    subroutine printOutput(aD,dB)
        type(assembledData),intent(inout) :: aD
        type(dataBlock), intent(in) :: dB
        character(len=120) :: addString
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        real(kind=rk) :: renorm,Atmp,Atra,Adt,Asca
        real(kind=rk) :: norm,rd,phiout,M(4,4),the
        integer :: j0,j1,j2,j3,j4
        real(kind=rk), allocatable :: MBS(:,:,:,:,:,:),MRT(:,:,:,:,:)
        integer :: ps_count,stream,k2   

        !Atot=sum(Aref)+sum(Asca)+sum(Aspec)+sum(Aabs)+sum(Astop)

        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)


        write(addString,*) "#RT-Engine: ",trim(adjustl(geometry_output_data1())), &
            & ",date,",trim(adjustl(date)),",time,",trim(adjustl(time)) 

        ps_count=size(dB%run_table)



        rd=pi/180.0_rk

        !parameter output
        open(unit=1,file=trim(adjustl(dB%output_details)))
        stream = 1
        write(stream,*) trim(adjustl(addString))        
        if(dB%addInfo) write(1,*) "#",trim(adjustl(dB%outputExt))
        call writeInfoScreen(dB,-1,stream,.true.)
        write(stream,'(A)') "run_table"
        write(stream,*) dB%run_table 
        !write(stream,'(A)') ""
        !write(stream,'(A)') "Raw albedos"  
        !call print_arr_with_message(sum(aD%Atra)/ps_count,db%ntot,stream,"Atra",fmt1)
        !call print_arr_with_message(sum(aD%Adt)/ps_count,db%ntot,stream,"Adt",fmt1)
        !call print_arr_with_message(sum(aD%Asca)/ps_count,db%ntot,stream,"Asca",fmt1)
        !call print_arr_with_message(sum(aD%Atra+aD%Adt)/ps_count,db%ntot,stream,"Atra+Adt",fmt1)
        !call print_arr_with_message(sum(aD%Aabs)/ps_count,db%ntot,stream,"Aabs",fmt1)
        !call print_arr_with_message(sum(aD%Asca+Aabs)/ps_count,db%ntot,stream,"Asca+Aabs",fmt1)
        !call print_arr_with_message(sum(aD%Astop)/ps_count,db%ntot,stream,"Astop",fmt1)
        !write(1,'(A20,"=",3X,E12.5)') "Atot",Atot
        !write(1,'(A20,"=",3X,E12.5)') "1 - Atot",1.0_rk - Atot

        write(stream,'(A)') str_line
        call print_single_val_with_message(aD%ATra,dB%run_table,stream,"Atra",fmt1,aD%rays)
        call print_single_val_with_message(aD%Adt,dB%run_table,stream,"Adt",fmt1,aD%rays)
        call print_single_val_with_message(aD%Aref,dB%run_table,stream,"Aref",fmt1,aD%rays)
        call print_single_val_with_message(aD%Aabs,dB%run_table,stream,"Aabs",fmt1,aD%rays)
        call print_single_val_with_message(aD%Astop,dB%run_table,stream,"Astop",fmt1,aD%rays)
        write(stream,'(A)') str_line
        call print_arr_with_message(aD%ATra,dB%run_table,stream,"Atra",fmt1)
        call print_arr_with_message(aD%Adt,dB%run_table,stream,"Adt",fmt1)
        call print_arr_with_message(aD%Aref,dB%run_table,stream,"Aref",fmt1)
        call print_arr_with_message(aD%Aabs,dB%run_table,stream,"Aabs",fmt1)
        call print_arr_with_message(aD%Astop,dB%run_table,stream,"Astop",fmt1)
        write(stream,*) "nrays",aD%rays


        !MRT = aD%MRT  
        Asca = sum(aD%Aref)/(sum(aD%rays))
        call move_alloc(aD%MRT , MRT)
        renorm=1.0_rk/(Asca)
        do j1=1,nthe
            norm=renorm/INORM(j1)
            MRT(:,:,j1,:,:)=norm*MRT(:,:,j1,:,:)
        enddo

        norm=1.0_rk/(Asca)  
        !!normalization


        

        write(stream,'(A)') "polarization states: Q+=2,Q-=3,U+=4,U-=5,V+=6,V+=7"
        write(stream,'(A)') "Electric field: EA=1,EB=2,EAB=3"
        write(stream,'(A,A)') dB%output_rt,"contains: phi, theta, polarization state,M11,M12, ... ,M43,M44"
        write(stream,'(A,A)') dB%output_rt,"contains: phi, theta, polarization state, Electric field ,M11,M12, ... ,M43,M44"
        

        call printScattererData(dB,1,stream)
        close(stream)


        
        !call outputReflectedRT(aD,dB)
        open(unit=1234,file=trim(adjustl(dB%output_rt)))
        if(dB%addInfo) write(1234,*) "#",trim(adjustl(dB%outputExt))
        if(dB%addInfo) write(1234,*) trim(adjustl(addString))
        write(1234,*) "nrays",aD%rays
        write(1234,*) ps_count,nthe,nphi 
        do j3=1,ps_count
            do j1=1,nthe
                do j2=1,nphi
                    write(1234,'(F7.2,F7.2,I3,4(1X,E30.16))') theGrid(j1)/rd,    &
         &          phiGrid(j2)/rd, dB%run_table(j3),                                   &
         &          MRT(1,:,j1,j2,j3)/aD%rays(j3)
                enddo
            enddo
        enddo
        close(1234)
        if(dB%cb_on) then
            !!Automatic allocation
            !MBS = aD%MBS
            call move_alloc(aD%MBS , MBS)
            !4x4x 3(EA,EB,EAB) x ntheb x nphib x pol_states
            !Complete output of the scattering radiation:
            MBS=norm*MBS
            open(unit=1234,file=trim(adjustl(dB%output_cb)))
            if(dB%addInfo) write(1234,*) "#",trim(adjustl(dB%outputExt))
            if(dB%addInfo) write(1234,*) trim(adjustl(addString))
            write(1234,*) ps_count,3,ntheb,nphib  
            do j4=1,ps_count
                do j1=1,3
                    do j2=1,ntheb
                        do j3=1,nphib
                            write(1234,'(F7.2,F7.2,I3,I3,4(1X,E30.16))')       & 
                                &   theGrid_cb(j2)/rd,                                  & 
                                &   phiGrid_cb(j3)/rd,dB%run_table(j4),j1,               &
                                &   MBS(1,:,j1,j2,j3,j4)/aD%rays(j4)
                        enddo
                    enddo
                enddo
            enddo
            close(1234)
        endif
       
    end subroutine

    !For test purposes
    subroutine writeInfoScreen(dB,nodes,stream,print_header)
        integer, intent(in) :: nodes,stream

        logical, intent(in) :: print_header
        type(dataBlock), intent(in) :: dB
        write(stream,*) "////////////////////////////"
        if(print_header) then
            write(stream,*) "RT-Engine"
            write(stream,*) trim(adjustl(geometry_output_data1()))
            if(nodes/=-1) write(stream,*) "MPI PROCESSES",nodes
            write(stream,*) "////////////////////////////"
        endif
!#OUTPUT_GENERATOR_SCRIPT1
        write(stream,'(A30,20X,E12.5)') "wavelength", dB%wavel
        write(stream,'(A30,20X,E12.5)') "mean_free_path", dB%ell
        write(stream,'(A30,20X,E12.5)') "cutoff_intensity", dB%Fstop
        write(stream,'(A30,20X,E12.5)') "incident_theta", dB%the0
        write(stream,'(A30,20X,E12.5)') "incident_phi", dB%phi0
        write(stream,'(A30,20X,E12.5)') "phi_output", dB%phiout
        write(stream,'(A30,20X,E12.5)') "tau_threshold", dB%tau_c
        write(stream,'(A30,20X,E12.5)') "volume_element_radius", dB%volrad
        write(stream,'(A30,20X,E12.5)') "medium_thickness_or_radius", dB%hr
        write(stream,'(A30,20X,E12.5)') "default_run_time", dB%runtime
        write(stream,'(A30,20X,I12)') "seed", dB%seed
        write(stream,'(A30,20X,I12)') "number_of_rays", dB%ntot
        write(stream,'(A30,20X,I12)') "max_sca_procs", dB%nsca
        write(stream,'(A30,20X,I12)') "rotations_every_n_path", dB%nrot
        write(stream,'(A30,20X,L12)') "additional_info", dB%addInfo
        write(stream,'(A30,20X,L12)') "estimateTime", dB%estimateTime
        write(stream,'(A30,20X,L12)') "force_to_stay_in", dB%forced_vm
        write(stream,'(A30,20X,L12)') "surface_attenuation", dB%surface_att
        write(stream,'(A30,20X,L12)') "overlap_probability", dB%overlap_prob
        write(stream,'(A30,20X,L12)') "override_albedo_to_1", dB%albedo_override
        write(stream,'(A30,20X,I12)') "cb_theta_angles", dB%ntheb
        write(stream,'(A30,20X,I12)') "cb_phi_angles", dB%nphib
        write(stream,'(A30,20X,I12)') "rt_the_angles", dB%nthe
        write(stream,'(A30,20X,I12)') "rt_phi_angles", dB%nphi
        write(stream,'(A30,20X,E12.5)') "distance_to_far_field", dB%ffboundary
        write(stream,'(A30,20X,L12)') "pol_state_Q_plus", dB%polQp
        write(stream,'(A30,20X,L12)') "pol_state_Q_minus", dB%polQm
        write(stream,'(A30,20X,L12)') "pol_state_U_plus", dB%polUp
        write(stream,'(A30,20X,L12)') "pol_state_U_minus", dB%polUm
        write(stream,'(A30,20X,L12)') "pol_state_V_plus", dB%polVp
        write(stream,'(A30,20X,L12)') "pol_state_V_minus", dB%polVm
        write(stream,'(A30,20X,L12)') "compute_cb", dB%cb_on
        write(stream,'(A30,2X,A30)') "output_details",adjustl(trim( dB%output_details))
        write(stream,'(A30,2X,A30)') "output_rt",adjustl(trim( dB%output_rt))
        write(stream,'(A30,2X,A30)') "output_cb",adjustl(trim( dB%output_cb))
        write(stream,'(A30,2X,A30)') "own_text_to_output",adjustl(trim( dB%outputExt))
        write(stream,'(A30,2X,A30)') "read_theta_grid_points_cb",adjustl(trim( dB%cb_the_points))
        write(stream,'(A30,2X,A30)') "tmatrix_location",adjustl(trim( dB%t_matrix_location))
!#OUTPUT_GENERATOR_SCRIPT2
        write(stream,*) str_line
    end subroutine

end module


