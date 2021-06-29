
!-------------------------------------------------------------------
! Interface for random number generator so that the switching
! of random number generators would be easy. 
! If openMP is used, check that the used PRNG is thread safe
!
! Copyright (C) 2018 Timo Väisänen, Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! **.*.2017 AP: Normal and poisson distribution 
! 19.3.2018 TV: Maintenance, cleaning, commenting
!
! List of functions:
! randNum() result(retVal)
!
! TODO
! *Remove copy-pasta
! *Add compression
! *Better error messaging
!-----------------------------


module rng
use common
!use iso_fortran_env
#ifdef __INTEL_COMPILER
use ifport
#endif
implicit none

!integer, parameter :: dp = real64
!real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)
contains

!!Generates one uniformly distributed random variable
function rand_num() result(retVal)
  real(kind=dp) :: retVal
  call random_number(retVal)
end function

!!Generates 3 uniformly distributed random variables
function rand_num_vec() result(vec)
  real(kind=dp) :: vec(3)
  call random_number(vec)
end function

!!Initialize random number generator. 
!!Requires getpid(), because 
!!each MPI process is seeded using
!!the their pid value
subroutine init_random_seed()
  INTEGER :: i, n, clock
  integer :: seed1
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed1 = abs( mod((clock*181)*((getpid()-83)*359), 104729)) 
  seed = seed1 + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  DEALLOCATE(seed)
end subroutine init_random_seed

!!Generates poisson distributed integers
!!takes lambda 
FUNCTION rnd_poisson(lambda) RESULT(k)
  INTEGER, INTENT(IN) :: lambda
  INTEGER :: k
  REAL(KIND=dp) :: u, p, L
  k = 0
  do while(k==0)
    p = 1.0_dp
    L = EXP(-1.0_dp * lambda)
    DO WHILE(p > L)
      k = k+1
      CALL RANDOM_NUMBER(u)
      p = p*u
    END DO
    k = k-1
  end do
END FUNCTION rnd_poisson


!!Generate random number with normal distribution
!!Uses Box Muller-algorithm
!!Give: mean and stdev
!!if stdev<=0 returns mean
function rnd_normal(mean,stdev) result(c)
  real(dp), intent(in) :: mean,stdev
  real(dp) :: c,temp(2), theta,r 
  if(stdev <= 0.0_dp) then
    c = mean
  else
    c = -1.0_dp
    do while(c<0.0_dp)
      call random_number(temp)
      r=(-2.0_dp*log(temp(1)))**0.5_dp
      theta = 2.0_dp*pi*temp(2)
      c= mean+stdev*r*sin(theta)
    enddo
  endif
end function

!!Generate random number with power law distribution Cx^alpha
!!n is the power index (must be positive)
!!lowlim: lower limit
!!maxlim: upper limit
function rnd_power_law(alpha,lowlim,maxlim) result(c)
  real(kind=dp), intent(in) :: alpha,lowlim,maxlim
  real(kind=dp) :: c,tmp
  call random_number(tmp)
  c = ((maxlim**(-alpha+1.0_dp)-lowlim**(-alpha+1.0_dp))*tmp+lowlim**(-alpha+1.0_dp))**(1.0_dp/(-alpha+1.0_dp))
end function

end module
