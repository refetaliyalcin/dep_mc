!------------------------------------------------
!Get random number's from 
!SIMD-oriented Fast Mersenne Twister (dSFMT) by
!Mutsuo Saito, Makoto Matsumoto(2007,2008)
!
!Uses Fortran interface for dSFMT by James Spencer(2012)
!
!2. Use the seed to initialise the RNG
!3. Get numbers by calling randNum() 
!
!Each module has own generator(rng0) which are 
!initialised with the givenSeed. The generator stores
!a certain amount of random numbers(RANDOMNUMBERSSTORED)
!from where the numbers can be accessed fast
!
!Copyright (C) 2016 Timo Väisänen and University of Helsinki
!All rights reserved.
!The new BSD License is applied to this software, see LICENSE.txt
!------------------------------------------------

!Define how many random numbers generator should store
#define     RANDOMNUMBERSSTORED 65536
module rng
    use dSFMT_interface
    implicit none
    integer, parameter :: rk=selected_real_kind(15)
    real(kind=rk), parameter :: pi = 4.0_rk*atan(1.0_rk)
    integer :: idum
    !module's private rng0
    type(dSFMT_t) :: rng0
    !$OMP threadprivate(rng0)
    private :: rng0,rk,pi
contains


    !------------------------------------------------
    !Initialise PRNG with the given seed
    !Calls dSFMT interface's dSFMT_init()
    !
    !IN: givenSeed: seed for the RNG
    !------------------------------------------------
    subroutine init_rng(givenSeed,steps)
        integer, intent(in) :: givenSeed
        integer, intent(in) :: steps
        !Parameter: seed, number of randoms, generator to operate
        !write(6,*) givenSeed,tid
        idum = 2
        call dSFMT_init(givenSeed, RANDOMNUMBERSSTORED, rng0) 
        if(steps/=0) then
            call dSFMT_jumpf(rng0,steps) 
        endif
    end subroutine

    !------------------------------------------------
    !Get random number from the generator's array
    !
    !OUT: retVal: Random number [0,1) 
    !------------------------------------------------
    function randNum() result(retVal)
        real(kind=rk) :: retVal
        retVal = get_rand_close_open(rng0)
    end function

    function rnd_normal(mean,stdev) RESULT(c)
        real(kind=rk) :: mean,stdev,c,temp(2), theta,r 
        IF(stdev <= 0.0_rk) THEN
            c = mean
        ELSE
            c = -1.0_rk
            do while(c<0)
                r=(-2.0_rk*log(randNum()))**0.5
                theta = 2.0_rk*PI*randNum()
                c= mean+stdev*r*sin(theta)
            end do
        endif
    end function



end module
