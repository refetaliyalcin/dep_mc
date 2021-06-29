!-------------------------------------------------------------------
! Compute interaction between volume elements
!
! Copyright (C) 2018 Johannes Markkanen, Timo Väisänen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 4.4.2018 TV: Cleaning
!
!-----------------------------

#include "../macros.inc"
module mie_extension
    use sfunctions
    use common
    use mie
    implicit none

    integer, parameter :: ERRORCODE = 999
    integer, parameter :: rk = dp;
    complex(kind=rk), allocatable :: vecMArr(:,:,:,:)
    complex(kind=rk), allocatable :: vecNArr(:,:,:,:)
    complex(kind=rk), allocatable :: vecMArr_cb(:,:,:,:)
    complex(kind=rk), allocatable :: vecNArr_cb(:,:,:,:)
    complex(kind=rk), allocatable :: vecMarr_tmp(:,:)
    complex(kind=rk), allocatable :: vecNarr_tmp(:,:)  
    !!Base vectors which are heavily used
    private :: vecMArr, vecNArr, vecMArr_cb, vecNArr_cb,vecMarr_tmp,vecNarr_tmp
    private :: rk,ERRORCODE,fillBaseVectors,init_base_vectors
    private :: calc_e_field
contains



!!Initialize mie_extension
subroutine init_mie(cb_on,ffBoundary,Nmax,theGrid,phiGrid,theGrid_cb,phiGrid_cb,k)
    logical, intent(in)                 :: cb_on                        !!check whether cb in on
    real(kind=rk), intent(in)           :: ffBoundary                   !!Distance to far field
    real(kind=rk), intent(in)           :: theGrid(:),phiGrid(:)        !!the and phi angles
    real(kind=rk), intent(in)           :: theGrid_cb(:),phiGrid_cb(:)  !!the and phi angles for cb
    complex(kind=rk), intent(in)        :: k                            !!wavenumber
    integer, intent(in)                 :: Nmax

    !!Make checks
    ASSERTI(Nmax>0,"NMAX has to be over 0")
    ASSERTI(size(theGrid)>0,"Module mie: size(theGrid)>0")
    ASSERTI(size(phiGrid)>0,"Module mie: size(phiGrid)>0")
    ASSERTI(ffboundary>0.0_rk,"Module mie: distance to far field boundary > 0")

    !!Initialize base vectors
    call init_base_vectors(Nmax,theGrid,phiGrid,ffboundary,k,vecMArr,vecNArr)
    ASSERTC2(allocated(vecMArr),.eqv.,.true.)
    ASSERTC2(allocated(vecNArr),.eqv.,.true.)

    !!And base vectors for cb
    if(cb_on) then
        ASSERTI(size(theGrid_cb)>0,"Module mie:: theGrid_cb<=0")
        ASSERTI(size(phiGrid_cb)>0,"Module mie:: phiGrid_cb<=0")
        call init_base_vectors(Nmax,theGrid_cb,phiGrid_cb, &
            & ffboundary,k,vecMArr_cb,vecNArr_cb)
            ASSERTC2(allocated(vecMArr_cb),.eqv.,.true.)
            ASSERTC2(allocated(vecNArr_cb),.eqv.,.true.)
    endif
end subroutine




!!Create base vectors
subroutine init_base_vectors(Nmax,theGrid,phiGrid,distToFF,k,M_arr,N_arr)
    integer, intent(in)                         :: Nmax
    complex(kind=rk), intent(in)                :: k                 !!wavenumber
    real(kind=rk), intent(in)                   :: phiGrid(:)        !!phi angles
    real(kind=rk), intent(in)                   :: theGrid(:)        !!the angles
    real(kind=rk), intent(in)                   :: distToFF          !!Distance to far field
    complex(kind=rk), allocatable, intent(out)  :: M_arr(:,:,:,:)    !!M base vectors
    complex(kind=rk), allocatable, intent(out)  :: N_arr(:,:,:,:)    !!N base vectors

    real(kind=rk) :: X(3)
    integer :: j1,j2,nn

    nn = (((Nmax*2+1)+3)*Nmax)/2
    ASSERTC2(nn,>,0)

    allocate(M_Arr(3,nn,size(theGrid),size(phiGrid)),N_Arr(3,nn,size(theGrid),size(phiGrid)))
    do j1=1,size(theGrid)
        do j2=1,size(phiGrid)
            x(1)=distToFF
            x(2)=theGrid(j1)
            x(3)=phiGrid(j2)
            if(theGrid(j1)<=0.0_rk) write(6,*) &
                & "WARNING! Module mie: Theta angles must be greater than 0"
            call fillBaseVectors(Nmax,x,M_Arr(:,:,j1,j2),N_Arr(:,:,j1,j2),k)
        enddo
    enddo
end subroutine



!!Compute electric field using base vectors M and N
subroutine calc_e_field(F,a_nm, b_nm,base_M,base_N)
    complex(kind=rk), intent(out)       :: F(3)
    complex(kind=rk), intent(in)        :: base_M(:,:),base_N(:,:)
    complex(kind=rk), intent(in)        :: a_nm(:), b_nm(:)

    integer :: j1

    F=cmplx(0.0,0.0,kind=rk)
    do j1 = 1, size(a_nm)
        F = F + a_nm(j1) * base_M(:,j1) + b_nm(j1) * base_N(:,j1)
    end do
end subroutine

!!Get the electric field at the direction the(j1) and phi(j2) using
!!standard base vectors
subroutine get_e_field(E, a_nm, b_nm,j1,j2)
    complex(kind=rk), intent(in)            :: a_nm(:), b_nm(:)
    complex(kind=rk), intent(out)           :: E(3)
    integer, intent(in)                     :: j1,j2

    call calc_e_field(E,a_nm, b_nm,vecMArr(:,:,j1,j2),vecNArr(:,:,j1,j2))
end subroutine

!!Get the electric field at the direction the(j1) and phi(j2) using
!!cb base vectors
subroutine get_e_field_cb(E,a_nm, b_nm,j1,j2)
    complex(kind=rk), intent(in)                :: a_nm(:), b_nm(:)
    integer, intent(in)                         :: j1,j2
    complex(kind=rk), intent(out)               :: E(3)

    call calc_e_field(E,a_nm, b_nm,vecMArr_cb(:,:,j1,j2),vecNArr_cb(:,:,j1,j2))
end subroutine

!!Get the electric field at the direction x(r,the,phi) using
!!cb base vectors
subroutine get_e_field2(E,x,a_nm, b_nm,  k)
    real(kind=rk), intent(in)       :: x(3)
    complex(kind=rk), intent(in)    :: a_nm(:), b_nm(:)
    complex(kind=rk), intent(in)    :: k
    complex(kind=rk), intent(out)   :: E(3)

    integer :: nn,Nmax

    nn = (((size(a_nm)*2+1)+3)*size(a_nm))/2
    if(.not. allocated(vecMarr_tmp)) then
        allocate(vecMarr_tmp(3,nn),vecNarr_tmp(3,nn))
    endif
    Nmax = int(sqrt(size(a_nm)+1.0_rk)) - 1  
    call fillBaseVectors(Nmax,x,vecMArr_tmp,vecNArr_tmp,k)
    call calc_e_field(E,a_nm, b_nm,vecMArr_tmp(:,:),vecNArr_tmp(:,:))
end subroutine



!!______________________________________________________!!
!!
!! The routine computes the SVWF expansion coefficients 
!! for time harmonic linearly and circularly polarized planewaves
!! propagating +z-direction  with the wave number k
!! 
!!@Note: TV: 23.12 Changed the order of polarizations to match IQUV
!!_____________________________________________________!!
subroutine inc_waves(Nmax, k, a_nm2, b_nm2)
    complex(kind=rk), dimension(:,:), intent(out)   :: a_nm2, b_nm2
    integer, intent(in)                             :: Nmax
    real(kind=rk), intent(in)                       :: k
    
    integer :: ind, n, m, mm, las, nm_in
    integer, dimension(:,:), allocatable :: indD
    real(kind=rk) :: scale, C, E0, omega
    complex(kind=rk) :: q
    complex(kind=rk), dimension(:), allocatable :: rotD
    complex(kind=rk), dimension((Nmax+1)**2-1) :: a_nm, b_nm
    
    E0 = 1.0_rk
    omega = k*299792458.0_rk
    ind = 0
    
    nm_in = (Nmax+1)**2-1
    
    do n = 1,Nmax
    
        scale = sqrt(dble(n*(n+1)))
    
        do m = -n, n 
            ind = ind + 1
            mm = abs(m)
            
            C = scale*E0*sqrt(pi*(2*n+1.0_rk)) / (n*(n+1.0_rk)) * sqrt(factorial(n+1)/factorial(n-1))
    
            q = -(cmplx(0.0_rk,1.0_rk,kind=rk)**n*k)/(omega*mu)
    
            if(mm == 1) then 
                a_nm(ind) = cmplx(0.0_rk,1.0_rk,kind=rk)**(n-1.0_rk) * C
                b_nm(ind) = -cmplx(0.0_rk,1.0_rk,kind=rk)**(n+1.0_rk) * C
    
                if(m == -1) then 
                a_nm(ind) = -a_nm(ind)
                
                end if
    
            else 
                a_nm(ind) = cmplx(0.0_rk,0.0_rk,kind=rk)
                b_nm(ind) = cmplx(0.0_rk,0.0_rk,kind=rk)
            end if
            a_nm(ind) = -a_nm(ind) 
                b_nm(ind) = -b_nm(ind) 
        end do
    end do
    
    
    las = int(Nmax+2*Nmax*(Nmax+1)+(4*Nmax**3/3.0_dp+2*Nmax**2+2*Nmax/3.0_dp)+0.1)
    
    allocate(rotD(las))
    allocate(indD(las,2))
    
    ! x-polarized wave
    call sph_rotation_sparse(pi, 0.0_rk, Nmax, rotD, indD)
    a_nm2(:,1) = sparse_matmul(rotD,indD,a_nm,nm_in)
    b_nm2(:,1) = sparse_matmul(rotD,indD,b_nm,nm_in)
    
    ! y-polarized wave
    call sph_rotation_sparse(pi, -pi/2.0_rk, Nmax, rotD, indD)
    a_nm2(:,2) = sparse_matmul(rotD,indD,a_nm,nm_in)
    b_nm2(:,2) = sparse_matmul(rotD,indD,b_nm,nm_in)
    
    ! rotated (45 deg.) linearly polarized waves
    call sph_rotation_sparse(0.0_rk, -pi/4.0_rk, Nmax, rotD, indD)
    a_nm2(:,4) = sparse_matmul(rotD,indD,a_nm2(:,2),nm_in)
    b_nm2(:,4) = sparse_matmul(rotD,indD,b_nm2(:,2),nm_in)
    
    a_nm2(:,3) = sparse_matmul(rotD,indD,a_nm2(:,1),nm_in)
    b_nm2(:,3) = sparse_matmul(rotD,indD,b_nm2(:,1),nm_in)
    
    ! left handed wave
    a_nm2(:,6) = a_nm2(:,2) + cmplx(0.0_rk,1.0_rk,kind=rk)*a_nm2(:,1)
    b_nm2(:,6) = b_nm2(:,2) + cmplx(0.0_rk,1.0_rk,kind=rk)*b_nm2(:,1)
    
    ! right handed wave
    a_nm2(:,5) = a_nm2(:,2) - cmplx(0.0_rk,1.0_rk,kind=rk)*a_nm2(:,1)
    b_nm2(:,5) = b_nm2(:,2) - cmplx(0.0_rk,1.0_rk,kind=rk)*b_nm2(:,1)

    deallocate(rotD,indD)
end subroutine inc_waves


!!Create base vectors 
subroutine fillBaseVectors(Nmax,vec,vecM,vecN,k)
    complex(kind=rk), dimension(:,:), intent(inout) :: vecM,vecN
    real(kind=rk), intent(in)                       :: vec(3)
    complex(kind=rk), intent(in)                    :: k

    integer, parameter :: inou = 1
    integer :: Nmax, n, m, ind, mm
    real(kind=rk), dimension(:), allocatable :: L, L1, L2
    real(kind=rk) :: r, theta, phi, q, omega
    complex(kind=rk) :: kr, alpha, beta, gamma, cc
    complex(kind=rk), dimension(:), allocatable :: sphj, sphy, sphh
    complex(kind=rk) :: P(3), B(3), C(3), Y, Y1, Y2, M_nm(3), N_nm(3)
    complex(kind=rk) :: F(3), G(3)

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
    F(:) = cmplx(0.0,0.0,kind=rk)
    G(:) = cmplx(0.0,0.0,kind=rk)

    do n = 1, Nmax
        alpha = sphh(n+1)
        beta = sqrt(dble(n*(n+1)))/kr * sphh(n+1)
        gamma = (n+1.0_rk)/kr * sphh(n+1) - sphh(n+2)

        allocate(L(n+1),L1(n+2),L2(n))

        call legendre2(n,cos(theta),L)
        call legendre2(n+1,cos(theta),L1)
        call legendre2(n-1,cos(theta),L2)
        !call system_clock(t2)

        do m = -n, n
            ind = ind+1
            mm = abs(m)
            
            cc = sqrt((2_rk*n+1.0_rk)*factorial(n-mm)/factorial(n+mm)/(4_rk*pi));
            
            ! Unnormalized complex scalar spherical harmonics
            Y=L(mm+1)*exp(cmplx(0.0_rk, m*phi,kind=rk));
            Y1=L1(mm+1)*exp(cmplx(0.0_rk, m*phi,kind=rk));
            
        
            if(mm == n) then
                Y2 = cmplx(0.0_rk,0.0_rk,kind=rk)
            else 
                Y2 = L2(mm+1)*exp(cmplx(0.0, m*phi,kind=rk)) 
            end if

            ! vector spherical harmonics
            P(:) = cmplx(0.0_rk,0.0_rk,kind=rk)
            P(1) = Y

            Y1=Y1*((n-mm+1.0_rk)/(n+1.0_rk))

            
            Y2=Y2*(dble(n+mm)/dble(n))

            B(:) = cmplx(0.0_rk,0.0_rk,kind=rk)
            B(2) = Y1-Y2
            B(3)=((cmplx(0.0_rk, m*(2.0_rk*n+1.0_rk),kind=rk))/(n*(n+1.0_rk)))*Y

        
            q=(sqrt(n*(n+1.0_rk)))/((n*2.0_rk+1.0_rk)*sin(theta));
            B = B*q

            C(:) = cmplx(0.0_rk,0.0_rk,kind=rk)
            C(2) = B(3) 
            C(3) = -B(2)
            
            ! Spherical vector wave functions
            M_nm = cc * alpha * C
            N_nm = cc * (beta*P + gamma * B)
            vecM(:,ind)=M_nm
            vecN(:,ind)=N_nm

        end do
        deallocate(L,L1,L2)
    end do
    deallocate(sphj, sphy, sphh)
end subroutine


end module

