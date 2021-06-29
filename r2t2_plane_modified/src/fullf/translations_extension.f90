
!-------------------------------------------------------------------
!
! Copyright (C) 2018 Johannes Markkanen, Timo Väisänen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 26.03.18 TV: Cleaning and commenting
!-----------------------------

module translations_extension
    use common 
    use sfunctions 
    implicit none
    
    !!These are used to prevent multiple allocation/deallocation
    !!inside function
    complex(kind=dp), allocatable   :: DD(:,:), locD(:,:)
    real(kind=dp), allocatable      :: denom_mu(:)
    private :: DD,locD,denom_mu
    
    contains
    
    !!Initialize module
    subroutine init_translations(Nmax)
        integer, intent(in) :: Nmax
        if(.not. allocated(DD)) then
            allocate(DD(2*(Nmax+2)+1,2*(Nmax+2)+1))
            allocate(locD(2*Nmax+1,2*Nmax+1))
            allocate(denom_mu(2*Nmax-1))
        endif
    end subroutine
     
    subroutine sph_rotation_gen2(angles, Nmax, spD)
        real(dp), intent(in) :: angles(3)
        integer, intent(in) :: Nmax
        complex(dp), dimension(:,:),intent(out) :: spD
    
        integer :: n, m, mu, ind1, ind2, mm, x, y, en, las, m2, mu2, en2
        real(dp) :: R(3,3)
        complex(dp) :: C(3,3), invC(3,3), D1(7,7)
        complex(dp) :: a1,a2,a3
        real(dp) ::  a, b, b2, cm, dm, dmm, delta, nom_a, nom_b, nom_b2
    
        ! Note Change of signs alpha and gamma in the rotation matrix 
        R = rotation_matrix(-angles(1), angles(2), -angles(3))
    
        spD(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    
        a1 = 1.0_dp/sqrt(2.0_dp)
        a2 = dcmplx(0.0,1.0)*(1.0_dp)/sqrt(2.0_dp)
        a3 = 1.0_dp
    
        C(1,:) = [a1, cmplx(0.0_dp,0.0_dp,kind=dp), -a1]
        C(2,:) = [-a2, cmplx(0.0_dp,0.0_dp,kind=dp), -a2]
        C(3,:) = [0.0_dp, 1.0_dp, 0.0_dp]
    
        invC(1,:) = [a1, a2, 0.0*a3]
        invC(2,:) = [0.0*a3, 0.0*a3, a3]
        invC(3,:) = [-a1, a2, 0.0*a3]
    
        D1(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        D1(3:5,3:5) = (transpose(matmul(invC,matmul(R,C))))
    
        DD = cmplx(0.0_dp,0.0_dp,kind=dp)
        DD(1:7,1:7) = D1
    
        las = 0 
        !********** n = 1 *************************************!!
        n = 1 
        do m2 = -n,n
            do mu2 = -n,n
                las = las + 1
                
                ind1 = n*n + n + m2
                ind2 = n*n + n + mu2
                    
                if(m2 >= 0 .and. mu2 >= 0) then
                    delta = 1
                end if
                
                if(m2 >= 0 .and. mu2 < 0) then
                    delta = (-1)**mu2
                end if
                
                if(m2 < 0 .and. mu2 >= 0) then
                    delta = (-1)**m2        
                end if
                
                if(m2 < 0 .and. mu2 < 0) then
                    delta = (-1)**(m2+mu2)   
                end if
                        
                !delta=1
                spD(ind1,ind2) = delta * D1(2+m2+n+1,2+mu2+n+1)
            
                
            end do
        end do
    
        !****************************************************
    
        do n = 2,Nmax
    
            locD(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    
            do mu = -n+1, n-1
                denom_mu(mu+n) = sqrt(dble((n+mu)*(n-mu)))
            end do
    
            do m = -n, n
                mm = -m
    
                nom_a = sqrt(dble((n+m)*(n-m)))
                nom_b = sqrt(dble((n+m)*(n+m-1)))        
                nom_b2 = sqrt(dble((n+mm)*(n+mm-1)))
    
    
                do mu = -n+1, n-1
                
                    a = nom_a / denom_mu(mu+n)
                    b = nom_b / (sqrt(2.0_dp)*denom_mu(mu+n))
                    b2 = nom_b2 / (sqrt(2.0_dp)*denom_mu(mu+n))
                        
                    x = mu+n+1
                    y = m+n+1
                    
                    locD(x,y) = D1(4,4) * a * DD(x+1,y+1) &
                    + b * D1(4,5) * DD(x+1,y) &
                    + b2 * D1(4,3) * DD(x+1,y+2)
                end do
                
                !___________________________________________________!
    
                mu = -n
    
                cm = sqrt(dble(n+m)*(n-m)/dble(n*(2*n-1)) );
                dm = sqrt( dble(n+m)*(n+m-1)/dble(2*n*(2*n-1)) );  
                dmm = sqrt( dble(n-m)*(n-m-1)/dble(2*n*(2*n-1)) );
                
                x = mu+n+1;
                y = m+n+1;
            
    
                locD(x,y) = (D1(3,4) * cm * DD(x+2,y+1) &
                + dm * D1(3,5) * DD(x+2,y) &
                + dmm * D1(3,3) * DD(x+2,y+2))
    
                !_________________________________________________!
    
                mu = n;
    
                x = mu+n+1;
                y = m+n+1;
    
                locD(x,y) = (D1(5,4) * cm * DD(x,y+1) &
                + dm * D1(5,5) * DD(x,y) &
                + dmm * D1(5,3) * DD(x,y+2));
    
            end do
            
            DD(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
            en = 2*(n+2)-1
            en2 = 2*n+1
            DD(3:en,3:en) = locD(1:en2,1:en2)
    
            ind1 = n**2;
            ind2 = (n+1)**2 - 1; 
    
            do m2 = -n,n
                do mu2 = -n,n
                    las = las + 1
    
                    ind1 = n*n + n + m2
                    ind2 = n*n + n + mu2
    
                    if(m2 >= 0 .and. mu2 >= 0) then
                        delta = 1
                    end if
                    
                    if(m2 >= 0 .and. mu2 < 0) then
                        delta = (-1)**mu2
                    end if
                    
                    if(m2 < 0 .and. mu2 >= 0) then
                        delta = (-1)**m2
                        
                    end if
                    
                    if(m2 < 0 .and. mu2 < 0) then
                        delta = (-1)**(m2+mu2)
                        
                    end if
    
                    !delta = 1
            
                    spD(ind1,ind2) = delta * locD(m2+n+1,mu2+n+1)
                
    
    
                end do
            end do
    
        end do
    
    end subroutine
    
    
    
    
    end module 
    