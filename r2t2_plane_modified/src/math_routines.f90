module math_routines
!!Coordinate system transformations and mathematical operations
!!
!!Copyright (C) 2016 Karri Muinonen, Timo Väisänen and University of Helsinki
!!All rights reserved.
!!The new BSD License is applied to this software, see LICENSE.txt

    use constants
    implicit none
    contains

!///////////////////////////////////////////////////////// 
    pure subroutine rayToNorm(K,cthe,sthe,cphi,sphi)
    !!Rotate direction vector (DV) from the ray to 
    !!the normal coordinate system.
        real(kind=rk), intent(inout) :: K(3)             !!IN: DV in ray coord. system
                                                         !!OUT: DV in norm coord. system
        
        real(kind=rk), intent(in) :: cthe,sthe,cphi,sphi !angles in norm coord. system
        real(kind=rk) :: q1,q2,q3                        !temporary values
        q1=K(1)
        q2=K(2)
        q3=K(3)
        K(1)= q1*cthe*cphi-q2*sphi+q3*sthe*cphi
        K(2)= q1*cthe*sphi+q2*cphi+q3*sthe*sphi
        K(3)=-q1*sthe+q3*cthe
    end subroutine

!///////////////////////////////////////////////////////// 
    pure subroutine normToRay(K,cthe,sthe,cphi,sphi)
    !!Rotate direction vector (DV) from the normal to 
    !!the ray coordinate system. 
        real(kind=rk), intent(inout) :: K(3)            !!IN: DV in ray coord. system
                                                        !!OUT: DV in norm coord. system
        real(kind=rk), intent(in)::cthe,sthe,cphi,sphi  !!angles in norm. coord. system
        real(kind=rk) :: q1,q2,q3                       !!temporary values
        q1=K(1)
        q2=K(2)
        q3=K(3)
        K(1)= q1*cthe*cphi+q2*cthe*sphi-q3*sthe
        K(2)=-q1*sphi+q2*cphi
        K(3)= q1*sthe*cphi+q2*sthe*sphi+q3*cthe    
    end subroutine    


!///////////////////////////////////////////////////////// 
    pure subroutine getRotationMatrix(R,c2psi,s2psi)
    !!Creates a [4,4] rotation matrice R from cos(2psi), sin(2psi).
        real(kind=rk),intent(out) :: R(4,4)             !!rotation matrice
        real(kind=rk), intent(in) :: c2psi,s2psi        !!cos(2psi),sin(2psi)
        R(1,1)=1.0_rk
        R(2,2)=c2psi
        R(3,3)=R(2,2)
        R(2,3)=s2psi
        R(3,2)=-R(2,3)
        R(4,4)=1.0_rk    
        R(1,2:4)=0.0_rk
        R(2,1)=0.0_rk
        R(2,4)=0.0_rk
        R(3,1)=0.0_rk
        R(3,4)=0.0_rk  
        R(4,1:3)=0.0_rk
    end subroutine


!///////////////////////////////////////////////////////// 
    subroutine getSphericalCoordinates(X,r,cthe,phi)
    !!Computes the spherical coordinates (r,cthe,phi) 
    !!from the Cartesian coordinates X=(x,y,z).
        real(kind=rk),intent(in) :: X(3)
        real(kind=rk),intent(out) :: r,cthe,phi
        real(kind=rk) :: cphi,sphi,sthe
        real(kind=rk),parameter :: tol = (10.0_rk**(-14))
        r=sqrt(X(1)**2+X(2)**2+X(3)**2)
        if(r<tol) then
            cthe=0.0_rk
            phi=0.0_rk
            write(6,*) "Warning in cartesianToSpherical: null vector."
            return
        endif
        cthe=X(3)/r
        if(abs(cthe)>1.0_rk-tol) then
            phi=0.0_rk
            return
        endif

        sthe=sqrt(1.0_rk-cthe**2)
        if (r*sthe<tol) then
            phi=0.0_rk
            return
        endif

        cphi=X(1)/(r*sthe)
        sphi=X(2)/(r*sthe)
        if (cphi>1.0_rk) then
            phi=0.0_rk
        elseif (cphi<-1.0_rk) then
            phi=pi
        else
            phi=acos(cphi)
        if (sphi<0.0_rk) phi=2.0_rk*pi-phi
        endif
    end subroutine

!///////////////////////////////////////////////////////// 
    pure subroutine getSphericalCoordinates2(X,r,cthe,phi)
    !!Computes the spherical coordinates (r,mu,phi) 
    !!from the Cartesian  coordinates X=(x,y,z).
    !!Expects that the length of vector X is=> abs(X)=1
        real(kind=rk),intent(in) :: X(3)
        real(kind=rk),intent(out) :: r,cthe,phi
        real(kind=rk) :: cphi,sphi,sthe
        real(kind=rk),parameter :: tol = (10.0_rk**(-12))
        r=1.0_rk

        cthe=X(3)
        if(abs(cthe)>1.0_rk-tol) then
            phi=0.0_rk
            return
        endif

        sthe=sqrt(1.0_rk-cthe**2)
        if (sthe<tol) then
            phi=0.0_rk
            return
        endif

        cphi=X(1)/sthe
        sphi=X(2)/sthe
        if (cphi>1.0_rk) then
            phi=0.0_rk
        elseif (cphi<-1.0_rk) then
            phi=pi
        else
            phi=acos(cphi)
        if (sphi<0.0_rk) phi=2.0_rk*pi-phi
        endif
    end subroutine

!///////////////////////////////////////////////////////// 
    pure function matXvec(M,Z) result(A)
    !A=MZ, A(4)=M(4,4)xZ(4)
        real(kind=rk) :: A(4)
        real(kind=rk), intent(in) :: M(4,4),Z(4)
        integer :: i,j
        A=0.0_rk
        do j=1,4            
            do i=1,4
                A(j)=A(j)+M(j,i)*Z(i)
            enddo
        enddo
    end function
    

!///////////////////////////////////////////////////////// 
    pure function dotProduct(X,Y) result(Z)
    !!Dot product XY=dot(X,Y)
        real(kind=rk),intent(in) :: X(3),Y(3)
        real(kind=rk) :: Z
        Z=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)    
    end function

!/////////////////////////////////////////////////////////
    pure function crossProduct(X,Y) result(Z)
    !Cross product. XY
        real(kind=rk),intent(in) :: X(3),Y(3)
        real(kind=rk) :: Z(3)
        Z(1)=X(2)*Y(3)-X(3)*Y(2)
        Z(2)=X(3)*Y(1)-X(1)*Y(3)
        Z(3)=X(1)*Y(2)-X(2)*Y(1)    
    end function



!/////////////////////////////////////////////////////////
    pure function combinedRotationFunc(I1,c2psi,s2psi) result(I2)
    !!Create rotation matrix getRotationMatrix() and do I2=R*I1
        real(kind=rk), intent(in) :: c2psi,s2psi,I1(4)
        real(kind=rk) :: I2(4)
        I2(1)=I1(1)
        I2(2)=I1(2)*c2psi+I1(3)*s2psi
        I2(3)=-I1(2)*s2psi+I1(3)*c2psi
        I2(4)=I1(4)
    end function
 

!/////////////////////////////////////////////////////////
    pure subroutine gaussLQuads(a,b,xn,wt)
    !!Gauss-Legendre Quadratures between [a,b]
    !!Coefficients are used to calculate integral
    !!/int_b^a(f(x))=dot_product(wt,f(x))
    !!
    !!Input:     a: lower limit
    !!           b: upper limit
    !!
    !!Output:    xn: abscissas
    !!           wt: weights
    !!
    !!Important phases
    !!#1:    Recursion for Legendre polynomials
    !!       also derivative
    !!#2:    Find roots using newton's method
    !!#3:    Root is now abscissas and derivative(pd)
    !!       can be used calculate wt
    !!#4     expand xn,wt to fit interval
        real(kind=rk),intent(in) :: a,b             !!lower and upper bound
        real(kind=rk), intent(out) :: xn(:),wt(:)   !!abscissas, weights
        real(kind=qp) :: dx,x,pd,x0                 
        real(kind=qp), allocatable :: p(:)
        integer :: n,j1,j2,j3
        n = size(xn)
        allocate(p(0:n))
        p(0) = 1.0_qp
        do j1=1,n
            x = cos(pi*(j1-0.25_qp)/(n+0.5_qp))  
            do j3=1,12
                p(1)=x
                do j2=1,n-1
                    p(j2+1)=((2.0_qp*j2+1)*x*p(j2)-j2*p(j2-1))/(j2*1.0_qp+1)
                enddo
                pd=n*(x*p(n)-p(n-1))/(x**2-1.0_qp)
                x0=x            
                x=x-p(n)/pd
                x0=x-x0
                if(abs(x0)<10*epsilon(x0)) exit
            enddo
            xn(n+1-j1) = x*(b-a)/2.0_qp+(a+b)/2.0_qp
            wt(n+1-j1) = 2.0_qp/((1-x**2)*pd**2)*((b-a)/2.0_qp)
        enddo
        deallocate(p)
        if(mod(n,2)==1) then
            xn((n+1)/2)=(a+b)/2.0_rk   
        endif
    end subroutine


!/////////////////////////////////////////////////////////
subroutine binarySearch(val,retVal,arr,nthe,nphi)
    integer,intent(out) :: retVal 
    real(kind=rk), intent(in) :: val
    real(kind=rk), intent(in) :: arr(:)
    integer :: imin,imid,imax
    integer, intent(in) :: nphi,nthe
    imin = 1
    imax = size(arr)
    if(val<=arr(1)) then
        retVal=1
        return
    endif
    if(val>arr(nphi*nthe-1)) then
        retVal=imax
        return
    endif
    do while(imax-imin>1)
        imid=ceiling((imax+imin)/2.0)
        if(arr(imid)>val) then
            imax=imid
        else if(arr(imid)<val) then
            imin=imid
        else
            retVal = imid
            return
        endif
        
    end do
    
    retVal=imax
end subroutine



end module

