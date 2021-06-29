#include "macros.inc"

#if defined (__INTEL_COMPILER)
#define INCLUDE_ME USE IFPORT
#elif (__GFORTRAN__)
#define INCLUDE_ME
#endif

!!Copyright (C) 2018 Timo Väisänen and University of Helsinki
!!All rights reserved.
!!The new BSD License is applied to this software, see LICENSE.txt
module prepareData 
    use constants
    use error_handler
    !dBlock type
    use typedefinitions
    use geometry
    use scatterer
    use grid_holder
    use input_reader
    use algorithm
    use mpi_tools
    use rng
    INCLUDE_ME 

    implicit none


contains

    subroutine prepare(dB)
    !!Initialize data which only master process requires to prepare
        type(dataBlock), intent(inout) :: dB
        real :: arg
        integer :: j1,clock
        integer, allocatable :: seed(:)
        integer :: size0,id
        !!init rng if seed = 0
        if(db%seed<=0) then
            call random_seed(size=size0)
            allocate(seed(size0))
            call system_clock(count=clock)
            id = GETPID()
            seed = clock+id*(/(j1-1,j1=1,size0)/)
            call random_seed(put=seed)
            call random_number(arg)
            dB%seed=int(arg*100000)
        endif
    end subroutine 


    subroutine prepare2(dB,aD,id)
    !!Distribute data and do final preparations
        type(dataBlock), intent(inout) :: dB    !dataBlock
        type(assembledData), intent(inout) :: aD
        integer, intent(in) :: id               !id of the process
        integer :: mySeed     
        mySeed=dB%seed
        !!initRng so that the every process has different seed. Do this first
        call init_rng(mySeed,id)
        if(debugging) write(6,*) "Seed of ",id," => ", mySeed,", jumps =>",id

        call fill_run_table(dB)
        call init_grid_holder(dB)
        call init_scatterer(dB)
        call init_geometry(dB)
        call init_algorithm(dB)
        call initAssembledData(dB,aD)
    end subroutine


    subroutine initAssembledData(dB,aD)
    !Initialize assembled Data
        type(dataBlock), intent(in) :: dB
        type(assembledData), intent(inout) :: aD
        integer :: ps_count   
        !MuellerMatrices,theta Angle, phiAngle
        if(.not. allocated(dB%run_table)) write(6,*) "You have to allocate run_table first"
        ps_count=size(dB%run_table)
        if(dB%cb_on) allocate(aD%MBS(4,4,5,ntheb,nphib,ps_count))
        if(dB%cb_on)  aD%MBS=0.0_rk
        allocate(aD%MRT(4,4,nthe,nphi,ps_count))        
        aD%MRT=0.0_rk    
        allocate(aD%Aref(ps_count))
        allocate(aD%Atra(ps_count))
        allocate(aD%Adt(ps_count))
        allocate(aD%Aabs(ps_count))
        allocate(aD%Astop(ps_count))
        allocate(aD%Aspec(ps_count))
        allocate(aD%rays(ps_count))
        aD%Aspec = 0.0_rk
        aD%Astop = 0.0_rk
        aD%Aabs = 0.0_rk
        aD%Adt = 0.0_rk
        aD%Atra = 0.0_rk
        aD%Aref = 0.0_rk
        aD%rays = 0
    end subroutine

    subroutine fill_run_table(dB)
        !!Fill run table which contains the info about init Stoke's configuration
        type(dataBlock), intent(inout) :: dB
        integer :: ps_count
        integer :: pol_st(6)
        ps_count=0
        if(dB%polQp) then
            ps_count=ps_count+1
            pol_st( ps_count) = 1
        endif        
        if(dB%polQm) then
            ps_count=ps_count+1
            pol_st( ps_count) = 2
        endif        
        if(dB%polUp) then
            ps_count=ps_count+1
            pol_st( ps_count) = 3
        endif        
        if(dB%polUm) then
            ps_count=ps_count+1
            pol_st( ps_count) = 4
        endif      
        !!right hand     
        if(dB%polVp) then
            ps_count=ps_count+1
            pol_st( ps_count) = 5
        endif        
        !!left hand
        if(dB%polVm) then
            ps_count=ps_count+1
            pol_st( ps_count) = 6
        endif     
        allocate(dB%run_table(ps_count))
        db%run_table(:) = pol_st(1:ps_count) 
        if(ps_count == 0) write(6,*) "WARNING! No polarization states selected"       
    end subroutine

end module
