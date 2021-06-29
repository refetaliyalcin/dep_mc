!-------------------------------------------------------------------
! Common mathematical functions and defitions between different
! modules
!
! Copyright (C) 2018 Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 19.3.2018 TV: Maintenance, cleaning and commenting
!
! List of types:
! *data_struct
! *Tmatrix
!
! List of functions:
! *factorial(a) result(c)
! *cart2sph(x) result(vec)
! *sph2cart(r,theta,phi) result(x)
! *sph2cart_vec(theta, phi, vec) result(vec2)
! *norm(r, rp) result(c)
! *binomial(n,k) result(c)
! *rotation_angles(x) result(vec)
! *rotation_matrix(a,b,g) result(rot)
! *crossCC(a,b) result(c)
! *truncation_order(ka) 
!-----------------------------

module common 
  use iso_fortran_env, only : real64
  implicit none 

  integer, parameter :: dp = real64
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: epsilon_C = 8.854187817_dp*(10.0_dp**(-12.0_dp))
  real(dp), parameter :: mu = 4.0_dp*pi*10.0_dp**(-7.0_dp)

  type data_struct
    integer :: Nmax, ind_a1, ind_a2, ind_b1, ind_b2
    integer :: Tmat_ind, ifT
    real(dp) :: euler_angles(3)
    real(dp) :: cp(3), r
    complex(dp) :: eps_r
    complex(dp), dimension(:,:), allocatable :: A, B
  end type data_struct

  type Tmatrix
    complex(dp), dimension(:,:), allocatable :: Taa, Tbb, Tab, Tba
  end type Tmatrix

contains

!Computes factorial a!
!Returns double precision real c
elemental function factorial(a) result(c)
  integer, intent(in) :: a
  integer :: i1
  real(kind=dp) :: c

  c=1
  do i1 = 1,a
    c = c * i1
  end do
end function factorial

!Converts cartesian vector (x) coordinates
!to spherical coordinates
!Returns vec=[r,theta,phi]
pure function cart2sph(x) result(vec)
  real(dp), intent(in) :: x(3)
  real(dp) :: vec(3)

  vec(1) = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! r
  vec(2) = acos(x(3)/vec(1)) ! theta
  vec(3) = atan2(x(2),x(1))   

end function cart2sph

!Converts cartesian vector (x) coordinates
!to spherical coordinates
!Returns vec=[r,theta,phi]
pure function sph2cart(r,theta,phi) result(x)
  real(dp), intent(in) :: r, theta, phi
  real(dp) :: x(3)

  x(1) = r*sin(theta)*cos(phi)
  x(2) = r*sin(theta)*sin(phi)
  x(3) = r*cos(theta)

end function sph2cart

!I think this is named wrongly because 
!it also rotates
pure function sph2cart_vec(theta, phi, vec) result(vec2)
  complex(dp), intent(in) :: vec(3)
  real(dp),    intent(in) :: theta, phi
  complex(dp) :: vec2(3)
  real(dp) :: H(3,3)

  H(:,1) = [sin(theta)*cos(phi), sin(theta)*sin(phi),cos(theta)]
  H(:,2) = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]
  H(:,3) = [-sin(phi), cos(phi), dble(0.0)];
  
  !H(1,:) = [sin(theta)*cos(phi), cos(theta)*cos(phi), -sin(phi)]
  !H(2,:) = [sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi)]
  !H(3,:) = [cos(theta), -sin(theta), dble(0.0)];
  vec2 = matmul(H,vec)
end function sph2cart_vec

!Returns the distance (c) between two points (r and rp)
pure function norm(r, rp) result(c)
  real(dp), dimension(3), intent(in) :: r, rp
  real(dp) :: c
  c = sqrt((r(1)-rp(1))**2 + (r(2)-rp(2))**2 + (r(3)-rp(3))**2)
end function norm

!Computes binomial coeffients (n,k)
pure function binomial(n,k) result(c)
  integer, intent(in) :: n, k
  integer :: i1
  real(dp) :: c

  c = 1.0_dp
  do i1 = 1,k
    c = c * (n + 1.0_dp - i1)/i1
  end do
end function binomial

!Generate rotation angles 
!CHECK! Also cart2sph, compatibility reasons?
pure function rotation_angles(x) result(vec)
  real(dp), intent(in) :: x(3)
  real(dp) :: vec(3)

  vec(1) = sqrt(x(1)**2 + x(2)**2 + x(3)**2) 
  vec(3) = atan2(x(2),x(1))     ! phi      
  vec(2) = atan2(x(1),x(3))   

end function rotation_angles


!Generate rotation matrix ()
pure function rotation_matrix(a,b,g) result(rot)
  real(dp), intent(in) :: a, b, g
  real(dp) :: rot(3,3)

  rot(1,1) = cos(b)*cos(g);
  rot(1,2) = cos(a)*sin(g) + sin(a)*sin(b)*cos(g);
  rot(1,3) = sin(a)*sin(g) - cos(a)*sin(b)*cos(g);

  rot(2,1) = -cos(b)*sin(g);
  rot(2,2) = cos(a)*cos(g) - sin(a)*sin(b)*sin(g);
  rot(2,3) = sin(a)*cos(g) + cos(a)*sin(b)*sin(g);

  rot(3,1) = sin(b);
  rot(3,2) = -sin(a)*cos(b);
  rot(3,3) = cos(a)*cos(b);
end function rotation_matrix

!Complex cross product
pure function crossCC(a,b) result(c)
  complex(dp), intent(in) :: a(3), b(3)
  complex(dp) :: c(3)
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = -a(1)*b(3) + a(3)*b(1)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function crossCC

!Generate truncation order 
!See Wiscombe’s criterion (x+cx**(1/3)+b) 
! W. J. Wiscombe, Improved Mie scattering algorithms, Appl. Opt. 19 (1980) 1505–150
!
!NOTICE: Wiscombe suggests to use c=4 and b=1. Should be tested?
pure function truncation_order(ka) result(Nmax)
  real(dp), intent(in) :: ka 
  integer :: Nmax

  if(ka > 1.0_dp) then
    Nmax = floor(ka + 3.0_dp* (ka)**(1.0_dp/3.0_dp))
  elseif(ka>=0.1_dp) then
    Nmax =  4
  elseif(ka>=0.01_dp) then
    Nmax = 3
  elseif(ka>=0.0001_dp) then
    Nmax = 2
  else
    Nmax = 1
  end if

end function truncation_order
end module common 
