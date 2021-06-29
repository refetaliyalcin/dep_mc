#include "macros.inc"
module input_reader
!!Input reader module for rtcb
!!
!!Uses hashtable and InputParser by Antti Penttilä, 2012, 2015
!!
!!There are two possible ways to give input
!!1: program reads input from "input.dat"
!!2: program reads input from the file given
!!   through the command line
!!
!!SYNTAX for the input is:
!!wavelength=0.71123123
!!
!!Possible parameters are listed below and in the
!!documentation coming with the program
!!
!!If the program can't find "input.dat", it will 
!!use default values which are written below.
!!If the program fails to read or find the given file,
!!the program will be terminated
!!
!!Copyright (C) 2016 Timo Väisänen and University of Helsinki
!!All rights reserved.
!!The new BSD License is applied to this software, see LICENSE.txt

    use InputParser
    !!Input parser for input value handling
    use typedefinitions
    !!get definition for dataBlock
    use error_handler
    implicit none

    private readInput2,filldataBlock,setDefaultValues
contains

!////////////////////////////////////////////////////////
    subroutine readInput(dB)
        type(dataBlock),intent(inout) :: dB
        real(kind=rk) :: time
        call setDefaultValues()
        call readInput2(time)
        if(check_if_errors()) return
        call filldataBlock(dB) 
        if(dB%cb_on) call read_cb_grid(dB)
        if(time>0.0_rk) dB%runtime = time
    end subroutine

!////////////////////////////////////////////////////////
    subroutine setDefaultValues()
    !!Read input to dataBlock
    !!
    !!INOUT: dB
    !!
    !!Fills inputParser with default values and 
    !!inputParser will overwrite them depending on the choice of input.
    !!After it, these values are read to dataBlock
        character(len=64) :: tmp    
        !init inputParser

!#INPUT_READER_SCRIPT1
        call init_input(42)
        call setup_value("wavelength",6.283185307179586_rk)      !!wavelength
        call setup_value("mean_free_path",112.028_rk)            !!incoherent mean free path
        call setup_value("cutoff_intensity",0.0001_rk)           !!cutoff intensity
        call setup_value("incident_theta",180.0_rk)              !!incident theta angle
        call setup_value("incident_phi",0.0_rk)                  !!incident phi angle
        call setup_value("phi_output",0.0_rk)                    !!phi output angle
        call setup_value("tau_threshold",50.0_rk)                !!threshold for optical depth
        call setup_value("medium_thickness_or_radius",86.7481_rk)!!thickness/radius
        call setup_value("volume_element_radius",7.5_rk)         !!volume element radius
        call setup_value("default_run_time",0.1_rk)              !!default run time in hours
        call setup_value("seed",0)                               !!seed
        call setup_value("number_of_rays",1000)                  !!number of rays
        call setup_value("max_sca_procs",100)                    !!max num of scatt. processes
        call setup_value("rotations_every_n_path",1)             !!asd
        call setup_value("additional_info",.true.)               !!add additional info to output
        call setup_value("estimateTime",.false.)                 !!estimate time
        call setup_value("force_to_stay_in",.false.)             !!force the volume element to stay inside the medium
        call setup_value("surface_attenuation",.false.)          !!attenuation starts from the surface
        call setup_value("overlap_probability",.true.)           !!if
        call setup_value("override_albedo_to_1",.false.)         !!override albedo to be 1
        call setup_value("cb_phi_angles",8)                      !!number of phi angles(cb-solution)
        call setup_value("rt_the_angles",64)                     !!number of theta angles(rt-solution)
        call setup_value("rt_phi_angles",48)                     !!number of theta angles(rt-solution)
        call setup_value("distance_to_far_field",1000000.0_rk)   !!distance to far field
        call setup_value("pol_state_Q_plus",.true.)              !!compute RT with polarization state Q+
        call setup_value("pol_state_Q_minus",.true.)             !!compute RT with polarization state Q-
        call setup_value("pol_state_U_plus",.true.)              !!compute RT with polarization state U+
        call setup_value("pol_state_U_minus",.true.)             !!compute RT with polarization state U-
        call setup_value("pol_state_V_plus",.true.)              !!compute RT with polarization state V+
        call setup_value("pol_state_V_minus",.true.)             !!compute RT with polarization state V-
        call setup_value("compute_cb",.true.)                    !!Enable/Disable coherent backscattering
        write(tmp,"(A)") "x_details.out"
        call setup_value("output_details",tmp,64)                !!Print input/output parameters
        write(tmp,"(A)") "x_output_rt.out"
        call setup_value("output_rt",tmp,64)                     !!RT-solution output filename
        write(tmp,"(A)") "x_output_cb.out"
        call setup_value("output_cb",tmp,64)                     !!CB-solution output filename
        write(tmp,"(A)") "dear diary vol 1"
        call setup_value("own_text_to_output",tmp,64)            !!extra information to output
        write(tmp,"(A)") "defaults/cbTheGrid.dat"
        call setup_value("read_theta_grid_points_cb",tmp,64)     !!read theta grid points from this file
        write(tmp,"(A)") "defaults/T_multi.h5"
        call setup_value("tmatrix_location",tmp,64)              !!The file which contains incoherent T-matrice
!#INPUT_READER_SCRIPT2
        !----------------------------------------
    end subroutine



!////////////////////////////////////////////////////////
    subroutine readInput2(time)
        !!readInput2()
        !!
        !!check which kind of input is given and
        !!overwrite default values
        !!
        !!The subroutine will check if there are
        !!any command line arguments. If there are none, 
        !!it will try to find "input.dat". Otherwise
        !!default values are used.
        !file where input is located
        character(len=128) :: filename
        !integer for various tasks
        integer :: i,ierr
        real(kind=rk), intent(out) :: time
        logical :: fileExists   
        !default input location
        filename="input.dat"
        
        !check whether user has given filename
        i = command_argument_count()      
        if(i==0) then
            write(6,*) "RT-engine will use default values. Give input file as an argument."
            write(6,*) str_line 
            return
            !Leave, use default values
        else
            !get filename from command line
            call get_command_argument(1,filename)
            !check that file exists
            inquire(file=filename,exist=fileExists)
            if(.not. fileExists) then
                !Critical error so terminate program
                call addError("input_reader: file """//trim(adjustl(filename))//""" does not exists")
                write(6,*) str_line 
                return
            endif
        endif
        write(6,*) "Reading input from: ",adjustl(trim(filename))
        write(6,*) str_line 
        call read_input(filename, .false. )

        if(i==2) then
            !get time
            call get_command_argument(2,filename)
            !read time
            read(filename,*,iostat=ierr) time  
            if(ierr/=0) then 
                call addError("error")
                return
            endif    
            write(6,*) "Time was given, overwrite default time"
        else 
            write(6,*) "Use default time"
            time = -1.0_rk
        endif
    end subroutine




    subroutine read_cb_grid(dB)
        type(dataBlock), intent(inout) :: dB
        logical :: fileExists
        integer :: count0,iostatus
        real(kind=rk), allocatable :: buffer(:)
        inquire(file=trim(adjustl(dB%cb_the_points)),exist=fileExists)
        if(.not. fileExists) then
            !Critical error so terminate program
            call addError("input_reader: the file """//trim(adjustl(dB%cb_the_points))//""" does not exists.")
            return
        endif      
        !!I don't think anyone is going to read more than 400 theta points
        allocate(buffer(1000))
        count0=0
        iostatus = 0
        open(unit=1337,file=trim(adjustl(dB%cb_the_points)))
        do while(iostatus==0)
            count0=count0+1
            read(1337,*,IOSTAT=iostatus) buffer(count0)
            if(iostatus/=0) exit 
        enddo
        if(count0==1 .and. iostatus<0) call addError("file: " &
            &   //trim(adjustl(dB%cb_the_points)) // "is empty")
        if(iostatus>0) call addError("input_reader: File """//trim(adjustl(dB%cb_the_points)) &
            &   // """ contains non-valid values. Ensure that the values are real.")
        close(1337)
        count0=count0-1
        !!copy buffer
        allocate(dB%theGrid_cb(count0))
        dB%theGrid_cb = buffer(1:count0)
        dB%ntheb = count0
        deallocate(buffer)

        ASSERTI(dB%ntheb>0,"Input_reader: read_cb_grid: ntheb must be >0")
        ASSERTI(.not. contains_elem(-pi,dB%theGrid_cb),"Input_reader: read_cb_grid: -pi is invalid value")
        ASSERTI(.not. contains_elem(pi,dB%theGrid_cb),"Input_reader: read_cb_grid: pi is invalid value")
        if(debugging) then
            ASSERTC(allocated(dB%theGrid_cb),.eqv.,.true.)
            ASSERTC(size(dB%theGrid_cb),>,0)
        endif

    end subroutine


!////////////////////////////////////////////////////////
    subroutine filldataBlock(dB)
    !!Get data from inputParser to datablock
    !!
    !!INOUT: DB: input will be put into this

        type(dataBlock),intent(inout) :: dB
        integer :: tmp
        tmp=64
!#INPUT_READER_SCRIPT3
        call get_value("wavelength",dB%wavel)                    !!wavelength
        call get_value("mean_free_path",dB%ell)                  !!incoherent mean free path
        call get_value("cutoff_intensity",dB%Fstop)              !!cutoff intensity
        call get_value("incident_theta",dB%the0)                 !!incident theta angle
        call get_value("incident_phi",dB%phi0)                   !!incident phi angle
        call get_value("phi_output",dB%phiout)                   !!phi output angle
        call get_value("tau_threshold",dB%tau_c)                 !!threshold for optical depth
        call get_value("medium_thickness_or_radius",dB%hr)       !!thickness/radius
        call get_value("volume_element_radius",dB%volrad)        !!volume element radius
        call get_value("default_run_time",dB%runtime)            !!default run time in hours
        call get_value("seed",dB%seed)                           !!seed
        call get_value("number_of_rays",dB%ntot)                 !!number of rays
        call get_value("max_sca_procs",dB%nsca)                  !!max num of scatt. processes
        call get_value("rotations_every_n_path",dB%nrot)         !!asd
        call get_value("additional_info",dB%addInfo)             !!add additional info to output
        call get_value("estimateTime",dB%estimateTime)           !!estimate time
        call get_value("force_to_stay_in",dB%forced_vm)          !!force the volume element to stay inside the medium
        call get_value("surface_attenuation",dB%surface_att)     !!attenuation starts from the surface
        call get_value("overlap_probability",dB%overlap_prob)    !!if
        call get_value("override_albedo_to_1",dB%albedo_override)!!override albedo to be 1
        call get_value("cb_phi_angles",dB%nphib)                 !!number of phi angles(cb-solution)
        call get_value("rt_the_angles",dB%nthe)                  !!number of theta angles(rt-solution)
        call get_value("rt_phi_angles",dB%nphi)                  !!number of theta angles(rt-solution)
        call get_value("distance_to_far_field",dB%ffboundary)    !!distance to far field
        call get_value("pol_state_Q_plus",dB%polQp)              !!compute RT with polarization state Q+
        call get_value("pol_state_Q_minus",dB%polQm)             !!compute RT with polarization state Q-
        call get_value("pol_state_U_plus",dB%polUp)              !!compute RT with polarization state U+
        call get_value("pol_state_U_minus",dB%polUm)             !!compute RT with polarization state U-
        call get_value("pol_state_V_plus",dB%polVp)              !!compute RT with polarization state V+
        call get_value("pol_state_V_minus",dB%polVm)             !!compute RT with polarization state V-
        call get_value("compute_cb",dB%cb_on)                    !!Enable/Disable coherent backscattering
        call get_value("output_details",dB%output_details,tmp)   !!Print input/output parameters
        call get_value("output_rt",dB%output_rt,tmp)             !!RT-solution output filename
        call get_value("output_cb",dB%output_cb,tmp)             !!CB-solution output filename
        call get_value("own_text_to_output",dB%outputExt,tmp)    !!extra information to output
        call get_value("read_theta_grid_points_cb",dB%cb_the_points,tmp)!!read theta grid points from this file
        call get_value("tmatrix_location",dB%t_matrix_location,tmp)!!The file which contains incoherent T-matrice
!#INPUT_READER_SCRIPT4
    end subroutine


end module
