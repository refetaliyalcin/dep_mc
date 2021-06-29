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
        call init_input(80)

        call setup_value("wavelength", 6.283185307179586_rk)
        call setup_value("single_scattering_albedo", 1.0_rk)
        call setup_value("mean_free_path",86.028_rk)
        call setup_value("medium_thickness_or_radius",86.7481_rk)
        call setup_value("number_of_rays",1000)
        call setup_value("min_relative_flux",0.000001_rk)
        call setup_value("max_number_of_scattering",9999)
        call setup_value("seed",0)                    
        call setup_value("threshold_for_optical_depth",50.0_rk)
        call setup_value("number_of_theta_angles",64)
        call setup_value("number_of_phi_angles",48)
        
        !PLANE/SLAB ONLY
        call setup_value("finite",.true.)
        call setup_value("phi_output_angle",0.0_rk)
        call setup_value("theta_angle_of_incidence",120.0_rk)
        call setup_value("phi_angle_of_incidence",0.0_rk)

        !FOR SPLINE PRESENTATION
        call setup_value("nrn",1000)
        call setup_value("ncm",32)
        
        !HM MEDIUM/surrounding environment of a spherical scatterer 
        call setup_value("host_medium_refractive_index_real", 1.0000_rk)
        call setup_value("host_medium_refractive_index_imag", 0.0_rk)        

        !parameters which require strings
        !Notice syntax
        write(tmp,"(A)") "rayleigh"
        call setup_value("scatterer_type",tmp,64)
        write(tmp,"(A)") "scatterer.in"
        call setup_value("load_scatterer_data_from",tmp,64)
        
        !OUTPUT
        write(tmp,"(A)") "details.out"
        call setup_value("details_output",tmp,64)
        write(tmp,"(A)") "rt1.out"
        call setup_value("rt_solution",tmp,64)
        write(tmp,"(A)") "rtcb1.out"
        call setup_value("cb_solution",tmp,64)
        call setup_value("output_add_additional_info",.false.)     
        call setup_value("estimate_time",.true.)
        call setup_value("I21_safety_check",.false.)
        write(tmp,"(A)") "simulation num. X"
        call setup_value("output_extra_info",tmp,64)
        
!MIE
        call setup_value("generate_mie_scatterer",.true.)
        call setup_value("scatterer_radius",0.5_rk)
        call setup_value("scatterer_real_ref",1.3_rk)
        call setup_value("scatterer_imag_ref",0.0_rk)
        call setup_value("points_in_spline_presentation",180)
        call setup_value("volume_fraction",0.1_rk)
        write(tmp,"(A)") "scatterer.out"
        call setup_value("save_scatterer_data_to",tmp,64)
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
        call get_value("wavelength", dB%wavel)
        call get_value("single_scattering_albedo", dB%ssalbedo)
        call get_value("mean_free_path",dB%ell1)
        call get_value("medium_thickness_or_radius",dB%hr)
        call get_value("number_of_rays",dB%ntot)
        call get_value("min_relative_flux",dB%Fstop)
        call get_value("max_number_of_scattering",dB%nsca)
        call get_value("seed",dB%seed)
        call get_value("threshold_for_optical_depth",dB%tau_c)
        call get_value("number_of_theta_angles",dB%nthe)
        call get_value("number_of_phi_angles",dB%nphi)

        !PLANE/SLAB ONLY
        call get_value("finite",dB%finite)
        call get_value("phi_output_angle",dB%phiout)
        call get_value("theta_angle_of_incidence",dB%the0)
        call get_value("phi_angle_of_incidence",dB%phi0)

        
        !SPLINE PRESENTATION
        call get_value("nrn",dB%nrn)
        call get_value("ncm",dB%ncm)
 
        !HM MEDIUM
        call get_value("host_medium_refractive_index_real", dB%mre)
        call get_value("host_medium_refractive_index_imag", dB%mim)       
                  
        !Scatterer type
        call get_value("scatterer_type",dB%scattererType,tmp)
        call get_value("load_scatterer_data_from",dB%scattererData,tmp)
          
        !output filenames
        call get_value("details_output",dB%output_details,tmp)
        call get_value("rt_solution",dB%output_rt,tmp)
        call get_value("cb_solution",dB%output_cb,tmp)
        call get_value("output_add_additional_info",dB%addInfo)
        call get_value("estimate_time",dB%estimateTime)
        call get_value("I21_safety_check",dB%I21test)
        call get_value("output_extra_info",dB%outputExt,tmp)

        call get_value("generate_mie_scatterer",dB%generateMieScatterer)
        call get_value("scatterer_radius",dB%scRadius)
        call get_value("scatterer_real_ref",dB%scRealRef)
        call get_value("scatterer_imag_ref",dB%scImagRef)
        call get_value("points_in_spline_presentation",dB%np)
        call get_value("volume_fraction",dB%vf)
        call get_value("save_scatterer_data_to",dB%saveScatterer,tmp)
!#INPUT_READER_SCRIPT4
    end subroutine


end module
