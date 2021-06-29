!-------------------------------------------------------------------
! Copyright (C) 2018 Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Also some code from Bohren and Huffman
!
! Modifications:
! 26.3.2018 TV: Maintenance, cleaning, commenting
! removed R2T2 and unused code: find them from
! the commit:f2d3ed19add39416b29ccef6a3e379a3bba21bc4
!
! TODO: 
! Comments
!-----------------------------

module mie
use sfunctions
use common
use translations


implicit none

contains


subroutine sphere_absorption2(a_in, b_in, k, a, epsr, Nmax, absb)
  complex(dp), dimension(:), intent(in) :: a_in, b_in
  real(dp), intent(in)  :: k, a
  real(dp), intent(out) :: absb
  integer, intent(in)  :: Nmax
  complex(dp), intent(in)  :: epsr

  complex(dp), dimension((Nmax+1)**2-1) ::a_n, b_n, c_n, d_n
  complex(dp) :: dn, cn
  real(dp) ::  int, x
  integer :: n, m, las

  x = k*a

  call mie_coeff_nm(Nmax, x, sqrt(epsr), a_n, b_n, c_n, d_n)

  las = 1
  int = 0.0
  do n = 1, Nmax
    do m = -n, n
      cn = -dble(1.0_dp / a_n(las) + 1.0_dp)
      dn = -dble(1.0_dp / b_n(las) + 1.0_dp)

      int = int + 1/k**2 * real(dn*abs(a_in(las))**2 + cn*abs(b_in(las))**2)
      
      las = las + 1
    end do

  end do

  absb =  int 
end subroutine sphere_absorption2


subroutine cross_sections(a_out, b_out, a_in, b_in, k, Nmax, Cext, Csca, Cabs)
  real(dp), intent(out) :: Cext, Cabs, Csca
  complex(dp), intent(in) :: k
  complex(dp), dimension(:), intent(in) :: a_out, a_in, b_out, b_in
  integer, intent(in) :: Nmax

  integer :: las, n, m
  real(dp) :: ee

  Cext=0.0
  Csca=0.0
  las = 1
  do n = 1, Nmax
    do m = -n,n
      ee = n*(n+1) / (2*n+1) * factorial(n+m) / factorial(n-m)
      Csca = Csca + dble(1/k**2 * (a_out(las)*conjg(a_out(las)) + b_out(las)*conjg(b_out(las))))
      Cext = Cext - real(1/k**2 * (a_out(las)*conjg(a_in(las)) + b_out(las)*conjg(b_in(las))))

      las = las + 1
    end do
  end do

  Cabs = Cext - Csca

end subroutine cross_sections
!******************************************************

subroutine scattering_cross_section(a_out, b_out, k, Nmax, Csca)
  real(dp), intent(out) :: Csca
  complex(dp), intent(in) :: k
  complex(dp), dimension(:), intent(in) :: a_out, b_out
  integer :: las, n, m
  integer, intent(in) :: Nmax

  Csca=0.0
  las = 1
  do n = 1, Nmax
      do m = -n,n
        Csca = Csca + dble(1/k**2 * (a_out(las)*conjg(a_out(las)) + &
                                      & b_out(las)*conjg(b_out(las))))
        las = las + 1
      end do
  end do
end subroutine scattering_cross_section



subroutine mie_coeff_nm(N, x, mr, a_nm, b_nm, c_nm, d_nm)
  integer, intent(in) :: N 
  real(dp), intent(in)  :: x
  complex(dp), intent(in) :: mr 
  complex(dp), dimension((N+1)**2-1), intent(out)  :: a_nm, b_nm, c_nm, d_nm

  complex(dp), dimension(N) :: a_n, b_n
  integer :: las, i1, m

  call BHMIE(N+1, x, mr, a_n, b_n)

  las = 0;
  do i1 = 1,N
    do  m = -i1,i1
        las = las + 1;
        a_nm(las) = a_n(i1);
        b_nm(las) = b_n(i1);
        c_nm(las) = 0.0_dp!c_n(i1);
        d_nm(las) = 0.0_dp!d_n(i1);

    
    end do
  end do
end subroutine mie_coeff_nm

!______________________________________________________!!
!
! The routine computes the SVWF expansion coefficients 
! for a time harmonic x-polarized planewave
! propagating +z-direction  with the wave number k
! 
!_____________________________________________________!!
subroutine planewave(Nmax, k, a_nm, b_nm)
  complex(dp), dimension((Nmax+1)**2-1), intent(out)  :: a_nm, b_nm
  integer, intent(in) :: Nmax
  real(dp), intent(in)  :: k

  integer :: ind, n, m, mm
  real(dp) :: scale, C, E0, omega
  complex(dp) :: q

  E0 = 1
  omega = k*299792458.0
  ind = 0

  do n = 1,Nmax

    scale = sqrt(dble(n*(n+1)))

    do m = -n, n 
        ind = ind + 1
        mm = abs(m)
        
        C = scale*E0*sqrt(pi*(2*n+1.0_dp)) / (n*(n+1.0_dp)) * sqrt(factorial(n+1)/factorial(n-1))

        q = -(cmplx(0.0_dp,1.0_dp,kind=dp)**n*k)/(omega*mu)

        if(mm == 1) then 
          a_nm(ind) = cmplx(0.0_dp,1.0_dp,kind=dp)**(n-1.0_dp) * C
          b_nm(ind) = -cmplx(0.0_dp,1.0_dp,kind=dp)**(n+1.0_dp) * C

          if(m == -1) then 
              a_nm(ind) = -a_nm(ind)
          
          end if

        else 
          a_nm(ind) = cmplx(0.0_dp,0.0_dp,kind=dp)
          b_nm(ind) = cmplx(0.0_dp,0.0_dp,kind=dp)
        end if
        a_nm(ind) = -a_nm(ind) 
        b_nm(ind) = -b_nm(ind) 
    end do
  end do
end subroutine planewave




!*****************************************************************
!
! Calculates fields from coefficients
!
! F = electric field
! G = Magnetic field
!
!*******************************************************************
subroutine calc_fields(a_nm, b_nm, k, Po, F, G, inou)
  complex(dp), dimension(:), intent(in) :: a_nm, b_nm
  real(dp), intent(in)                  :: Po(3)
  complex(dp), intent(in)               :: k
  complex(dp), intent(out)              :: F(3), G(3) 
  integer, intent(in)                   :: inou
  integer :: Nmax, n, m, ind, mm
  real(dp) :: r, theta, phi, vec(3), q, omega
  complex(dp) :: kr, alpha, beta, gamma, cc
  complex(dp), dimension(:), allocatable :: sphj, sphy, sphh
  complex(dp) :: P(3), B(3), C(3), Y, Y1, Y2, M_nm(3), N_nm(3)
  real(dp), dimension(:), allocatable :: L, L1, L2

  Nmax = int(sqrt(dble(1+size(a_nm))))-1

  vec = cart2sph(Po)

  r = vec(1)
  theta = vec(2)
  phi = vec(3)

  kr = k*r

  omega = real(k)*299792458.0

  allocate(sphj(Nmax+2), sphy(Nmax+2), sphh(Nmax+2))


  ! spherical bessel functions at kr
  if(inou == 0) then
    call cspherebessel(Nmax+1,kr, sphj, sphy)
    sphh = sphj
  else 
    sphh = sphankel(Nmax+1, kr)
  end if


  ind = 0
  F(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
  G(:) = cmplx(0.0_dp,0.0_dp,kind=dp)

  do n = 1, Nmax
    alpha = sphh(n+1)
    beta = sqrt(dble(n*(n+1)))/kr * sphh(n+1)
    gamma = (n+1.0_dp)/kr * sphh(n+1) - sphh(n+2)

    allocate(L(n+1),L1(n+2),L2(n))

    call legendre2(n,cos(theta),L)  
    call legendre2(n+1,cos(theta),L1)
    call legendre2(n-1,cos(theta),L2)

    q=(sqrt(n*(n+1.0_dp)))/((n*2.0_dp+1.0_dp)*sin(theta));

    do m = -n, n
      ind = ind+1
      mm = abs(m)
      cc = sqrt((2.0_dp*n+1.0_dp)*factorial(n-mm)/factorial(n+mm)/(4.0_dp*pi));
    
      ! Unnormalized complex scalar spherical harmonics
      Y=L(mm+1)*exp(cmplx(0.0_dp, m*phi,kind=dp));
      Y1=L1(mm+1)*exp(cmplx(0.0_dp, m*phi,kind=dp));
    

      if(mm == n) then
        Y2 = cmplx(0.0_dp,0.0_dp,kind=dp)
      else 
        Y2 = L2(mm+1)*exp(cmplx(0.0_dp, m*phi,kind=dp)) 
      end if

      ! vector spherical harmonics
      P(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
      P(1) = Y

      Y1=Y1*((n-mm+1.0_dp)/(n+1.0_dp))

    
      Y2=Y2*(dble(n+mm)/dble(n))

      B(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
      B(2) = Y1-Y2
      B(3)=((cmplx(0.0_dp, m*(2*n+1.0_dp),kind=dp))/(n*(n+1.0_dp)))*Y

      B = B*q

      C(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
      C(2) = B(3) 
      C(3) = -B(2)

      ! Spherical vector wave functions
      M_nm = cc * alpha * C
      N_nm = cc * (beta*P + gamma * B)
      
      F = F + a_nm(ind) * M_nm + b_nm(ind) * N_nm 
      G = G + k * (a_nm(ind) * N_nm + b_nm(ind) * M_nm)
      
    end do
    deallocate(L,L1,L2)
  end do

  F = sph2cart_vec(theta, phi, F)
  G = sph2cart_vec(theta, phi, G)
  G = G/(cmplx(0.0_dp,1.0_dp,kind=dp) * omega * mu)

end subroutine calc_fields

subroutine compute_cluster_absorption(sphere, x, k, Cabs)
  type (data_struct), dimension(:), intent(in) :: sphere
  complex(dp), dimension(:), intent(in) :: x
  complex(dp), intent(in) :: k
  real(dp), intent(out) :: Cabs

  integer :: sph, a1, a2, b1, b2, Nmax
  real(dp) :: Cabs2, rad,cp(3)
  complex(dp) ::  epsr

  Cabs = 0.0
  do sph = 1, size(sphere)
      
    a1 = sphere(sph)%ind_a1
    a2 = sphere(sph)%ind_a2
    b1 = sphere(sph)%ind_b1
    b2 = sphere(sph)%ind_b2
    cp = sphere(sph)%cp
    rad = sphere(sph)%r
    epsr = sphere(sph)%eps_r

    Nmax = sphere(sph)%Nmax
    

    call sphere_absorption2(x(a1:a2),x(b1:b2), dble(k), sphere(sph)%r, sphere(sph)%eps_r, Nmax, Cabs2)

    Cabs = Cabs + Cabs2

  end do
end subroutine compute_cluster_absorption



subroutine mueller_matrix_coeff(a_nm, b_nm, a_nm2, b_nm2, k, N_theta, N_phi, S_out)

  complex(dp), dimension(:), intent(in) :: a_nm, b_nm, a_nm2, b_nm2
  complex(dp), intent(in) :: k
  integer, intent(in) :: N_phi, N_theta
  real(dp), dimension(:,:), allocatable, intent(out) :: S_out

  complex(dp), dimension(3) :: E_out, E_out2, H_out, H_out2
  integer :: i1, i2, las
  real(dp) :: theta, phi, abcd(2,2), r(3), RR, unit_th(3), unit_phi(3)
  complex(dp), dimension(N_phi*N_theta) :: f11,f12,f21,f22
  complex(dp) ::S1, S2, S3, S4
  complex(dp) :: i

  i = cmplx(0.0_dp,1.0_dp,kind=dp)
  RR = 1.0d6
  allocate(S_out(N_phi*N_theta,18))


  las = 0
  do i1 = 1, N_phi
    do i2 = 1, N_theta

      theta = pi * (i2-1) / (N_theta)  + pi/N_theta/2.0
      phi = 2*pi * (i1-1) / N_phi
        
      abcd(1,:) = [cos(phi), sin(phi)]
      abcd(2,:) = [sin(phi), -cos(phi)]

      r(1) = RR*sin(theta)*cos(phi)
      r(2) = RR*sin(theta)*sin(phi)
      r(3) = RR*cos(theta)

      call calc_fields(a_nm, b_nm, k, r, E_out, H_out, 1)
      call calc_fields(a_nm2, b_nm2, k, r, E_out2, H_out2, 1)


      unit_th = [cos(theta) * cos(phi), sin(phi) * cos(theta), -sin(theta)]
      unit_phi = [-sin(phi), cos(phi),0.0_dp];

      f11(las+1) = k*RR / exp(i*k*(RR)) * dot_product(E_out,unit_th);
      f21(las+1) = k*RR / exp(i*k*(RR)) * dot_product(E_out,unit_phi);    
      f12(las+1) = k*RR / exp(i*k*(RR)) * dot_product(E_out2,unit_th);
      f22(las+1) = k*RR / exp(i*k*(RR)) * dot_product(E_out2,unit_phi);

      S1 = -i * (f21(las+1)*abcd(1,2) + f22(las+1)*abcd(2,2))
      S2 = -i * (f11(las+1)*abcd(1,1) + f12(las+1)*abcd(2,1))
      S3 = i * (f11(las+1)*abcd(1,2) + f12(las+1)*abcd(2,2))
      S4 = i * (f21(las+1)*abcd(1,1) + f22(las+1)*abcd(2,1))

      S_out(las+1,1) = phi
      S_out(las+1,2) = theta

      ! Mueller matrix
      S_out(las+1,3)=(abs(S1)**2 + abs(S2)**2 + abs(S3)**2 + abs(S4)**2)/2.0 !S11
      S_out(las+1,4)=(abs(S2)**2 - abs(S1)**2 + abs(S4)**2 - abs(S3)**2)/2.0 !S12
      S_out(las+1,5) = -real(S2*conjg(S3) + S1*conjg(S4)) !S13
      S_out(las+1,6) = -imag(S2*conjg(S3) - S1*conjg(S4)) !S14
      S_out(las+1,7)=(abs(S2)**2 - abs(S1)**2 + abs(S3)**2 - abs(S4)**2)/2.0 !S21
      S_out(las+1,8)=(abs(S1)**2 - abs(S3)**2 - abs(S4)**2 + abs(S2)**2)/2.0 !S22
      S_out(las+1,9) = real(S2*conjg(S3) - S1*conjg(S4)) !S23
      S_out(las+1,10) = -imag(S2*conjg(S3) + S1*conjg(S4)) !S24
      S_out(las+1,11) = real(S2*conjg(S4) + S1*conjg(S3)) !S31
      S_out(las+1,12) = -real(S2*conjg(S4) - S1*conjg(S3)) !S32
      S_out(las+1,13) = -real(S1*conjg(S2) + S3*conjg(S4)) !S33
      S_out(las+1,14) = imag(S2*conjg(S1) + S4*conjg(S3)) ! S34 
      S_out(las+1,15) = imag(S4*conjg(S2) + S1*conjg(S3)) ! S41 
      S_out(las+1,16) = imag(S4*conjg(S2) - S1*conjg(S3)) ! S42 
      S_out(las+1,17) = imag(S1*conjg(S2) - S3*conjg(S4)) ! S43 
      S_out(las+1,18) = -real(S1*conjg(S2) - S3*conjg(S4)) ! S44

      las = las + 1 
    end do
  end do
end subroutine mueller_matrix_coeff


pure subroutine tr_T(Taa, Tbb, Tab, Tba, k, crs)
  complex(dp), dimension(:,:), intent(in) :: Taa, Tbb, Tab, Tba
  real(dp), intent(out)     :: crs(4)
  real(dp), intent(in)      :: k
  integer :: nm, mu, nu
  real(dp) :: T1, T2

  nm = size(Taa,1)

  T1 = 0.0_dp
  T2 = 0.0_dp

  do mu = 1, nm

    T1 = T1 + (real(Taa(mu, mu)) + real(Tbb(mu, mu)) + real(Tab(mu, mu)) + real(Tba(mu, mu))) 

    do nu = 1,nm

      T2 = T2 + (abs(Taa(mu, nu))**2 + abs(Tbb(mu, nu))**2 + &
            abs(Tab(mu, nu))**2 + abs(Tba(mu, nu))**2)
    
    end do
  end do

  crs(1) = -2*pi*T1/k**2
  crs(2) = 2*pi*T2/k**2
  crs(3) = (-2*pi*T1/k**2-2*pi*T2/k**2)
  crs(4) = T1 + T2

  !print*,'Ave. Extinction cross section:', -2*pi*T1/k**2 
  !print*,'Ave. Scattering cross section:', 2*pi*T2/k**2
  !print*,'Ave. absorption cross section:', -2*pi*T1/k**2 -2*pi*T2/k**2
  !print*,'Tr_(real(T) + T adj(T)):', T1 + T2

end subroutine tr_T


include "BHMIE.f90"

end module mie

