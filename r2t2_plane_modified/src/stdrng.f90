!------------------------------------------------
!Get random numbers from fortran's intrinsic PRNG
!Don't use in parallel work!
!1. Ask module to generate seed from the given input
!(initSeed())
!2. Use seed to initialise PRNG
!3. Get numbers by calling randNum()
!
!Copyright (C) 2016,2017 Timo Väisänen and University of Helsinki
!All rights reserved.
!The new BSD License is applied to this software, see LICENSE.txt
!------------------------------------------------
module rng
    implicit none
    integer, parameter :: rk=selected_real_kind(15)
    private :: serialSeed,rk
contains

    !------------------------------------------------
    !Initialise PRNG with the given seed
    !
    !IN: GivenSeed: seed for the RNG
    !------------------------------------------------
    subroutine initRng(givenSeed,steps)
        integer, intent(in) :: givenSeed
        integer, intent(in) :: steps
        call serialSeed(givenSeed)
    end subroutine

    !------------------------------------------------
    !Initialise seed    
    !Seed will be the current time from OS's clock
    !if givenSeed<=0. Otherwise seed=givenSeed
    !
    !IN: givenSeed: input for seed
    !OUT: seed: seed for RNG
    !------------------------------------------------
    function initSeed(givenSeed) result(seed)
        integer, intent(in) :: givenSeed
        integer :: seed
        if(givenSeed>0) then
            seed=givenSeed
        else
            call system_clock(count=seed)
        endif
    end function

    !------------------------------------------------
    !Get a random number using Fortran's intrinsic
    !random_number()
    !
    !OUT: retVal: Random number [0,1) 
    !------------------------------------------------    
    function randNum() result(retVal)
        real(kind=rk) :: retVal
        call random_number(retVal)
    end function

    !------------------------------------------------
    !Private function which creates a seed which can
    !be used in the RGN. 
    !Single scalar value is not applicable and that's
    !why further processing is needed.
    !
    !IN: givenSeed: seed for the RNG
    !------------------------------------------------
    subroutine serialSeed(givenSeed)
        integer, intent(in) :: givenSeed
        integer :: j1, n
        integer, allocatable :: seed(:)
        call random_seed(size = n)
        allocate(seed(n))
        seed = givenSeed + 37*(/(j1-1,j1=1,n)/)
        call random_seed(put = seed)
        deallocate(seed)
    end

end module
