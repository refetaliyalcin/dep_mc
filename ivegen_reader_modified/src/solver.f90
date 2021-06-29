!-------------------------------------------------------------------

! Copyright (C) 2018 Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 23.3.2018 TV: Maintenance, cleaning, commenting
!
! TODO
! *Commenting
!
! TODO!
! check intents!!!!! intent(out) can be messy and ignore data which can be 
! relevant to the algorithm
!-----------------------------
module solver
use common
use translations
use possu
use octtree
use mie
implicit none


type data
   complex(dp), dimension(:,:), allocatable :: a_out, b_out
end type data


contains 

!****************************************************************************


pure subroutine rot_Tmatmul(Taa, Tab, Tba, Tbb, a_in, b_in, rotD, nm, a_out, b_out)
    complex(dp), dimension(nm,nm), intent(in)   :: Taa, Tab, Tba, Tbb
    complex(dp), dimension(nm,2), intent(in)    :: a_in, b_in
    complex(dp), dimension(nm,2), intent(out)   :: a_out, b_out
    complex(dp), dimension(:,:), intent(in)     :: rotD
    integer, intent(in)                         :: nm

    complex(dp), dimension(nm,nm) :: Taa_r, Tab_r, Tba_r, Tbb_r
    integer :: Nmax, i1

    Nmax = int(sqrt(nm + 1.0_dp) - 1.0_dp)

    !____rotate T-matrix_______

    do i1 = 1, nm

        !Taa_r(:,i1) = sparse_matmul(rotD,indD,Taa(:,i1),nm)
        !Tab_r(:,i1) = sparse_matmul(rotD,indD,Tab(:,i1),nm)
        !Tba_r(:,i1) = sparse_matmul(rotD,indD,Tba(:,i1),nm)
        !Tbb_r(:,i1) = sparse_matmul(rotD,indD,Tbb(:,i1),nm)

        
    end do

    Taa_r = matmul(transpose(conjg(rotD)),matmul(Taa,rotD))
    Tab_r = matmul(transpose(conjg(rotD)),matmul(Tab,rotD))
    Tba_r = matmul(transpose(conjg(rotD)),matmul(Tba,rotD))
    Tbb_r = matmul(transpose(conjg(rotD)),matmul(Tbb,rotD))

    a_out(:,1) = matmul(Taa_r,a_in(:,1)) + matmul(Tab_r,b_in(:,1)) 
    b_out(:,1) = matmul(Tba_r,a_in(:,1)) + matmul(Tbb_r,b_in(:,1)) 

    a_out(:,2) = matmul(Taa_r,a_in(:,2)) + matmul(Tab_r,b_in(:,2)) 
    b_out(:,2) = matmul(Tba_r,a_in(:,2)) + matmul(Tbb_r,b_in(:,2)) 
end subroutine rot_Tmatmul


subroutine vinc_xy(sphere, b_vecx, b_vecy, k)
    type (data_struct), dimension(:), intent(inout) :: sphere
    complex(dp), dimension(:), allocatable, intent(out) :: b_vecx, b_vecy
    complex(dp), intent(in) :: k

    complex(dp), dimension(:,:), allocatable :: a_out, b_out
    integer ::  Nspheres, n, loc, sph, bsize, las, nm_in, N_in, Nmax2
    complex(dp) :: phase_shift
    complex(dp), dimension(:), allocatable :: rotD
    integer, dimension(:,:), allocatable :: indD
    real(dp) :: delta, cp(3)

    Nspheres = size(sphere)

    bsize = 0
    do sph = 1, Nspheres
        Nmax2 = sphere(sph)%Nmax
        loc = (Nmax2+1)**2 - 1
        sphere(sph)%ind_a1 = bsize + 1 
        sphere(sph)%ind_a2 = bsize + loc
        sphere(sph)%ind_b1 = bsize + loc + 1 
        sphere(sph)%ind_b2 = bsize + 2*loc  
        bsize = bsize + 2*loc
    end do

    allocate(b_vecx(bsize), b_vecy(bsize))
    b_vecx(:) = dcmplx(0.0,0.0)
    b_vecy(:) = dcmplx(0.0,0.0)



    do sph = 1, Nspheres
        N_in = sphere(sph)%Nmax
        nm_in = (N_in+1)**2-1

        allocate(a_out(nm_in,2),b_out(nm_in,2))
        call planewave(N_in, real(k), a_out(:,1), b_out(:,1))

        las = 0
        do n = 1,N_in
            las = las + (2*n+1)**2
        end do

        allocate(rotD(las))
        allocate(indD(las,2))

        call sph_rotation_sparse(0.0d0, -pi/2.0d0, N_in, rotD, indD)
        a_out(:,2) = sparse_matmul(rotD,indD,a_out(:,1),nm_in)
        b_out(:,2) = sparse_matmul(rotD,indD,b_out(:,1),nm_in)


        cp = sphere(sph)%cp
        delta = cp(3)
        phase_shift = cdexp(dcmplx(0.0, dble(k)*delta))

        a_out = a_out * phase_shift
        b_out = b_out * phase_shift
        

        b_vecx(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) =   a_out(:,1)
        b_vecx(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) =   b_out(:,1)

        b_vecy(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) =   a_out(:,2)
        b_vecy(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) =   b_out(:,2)

        deallocate(a_out, b_out)
        deallocate(rotD,indD)
    end do
end subroutine vinc_xy

!*********************************

subroutine cluster_rhs(Nmax, a_in, b_in, sphere, b_vec, k)
    type (data_struct), dimension(:), intent(inout) :: sphere
    complex(dp), dimension(:), intent(in) :: a_in, b_in
    complex(dp), dimension(:), allocatable, intent(out) :: b_vec

    complex(dp), dimension(:), allocatable :: a_out, b_out, c_nm1, d_nm1,a_nm1,b_nm1
    complex(dp) :: k
    integer ::  Nspheres, loc, sph, bsize, nm, nm_in, N_in, Nmax2
    integer, intent(in) :: Nmax

    nm = size(a_in)


    Nspheres = size(sphere)

    bsize = 0
    do sph = 1, Nspheres
        Nmax2 = sphere(sph)%Nmax
        loc = (Nmax2+1)**2 - 1
        sphere(sph)%ind_a1 = bsize + 1 
        sphere(sph)%ind_a2 = bsize + loc
        sphere(sph)%ind_b1 = bsize + loc + 1 
        sphere(sph)%ind_b2 = bsize + 2*loc  

        bsize = bsize + 2*loc
        
    end do

    allocate(b_vec(bsize))

    b_vec(:) = dcmplx(0.0,0.0)


    do sph = 1, Nspheres
        N_in = sphere(sph)%Nmax
        nm_in = (N_in+1)**2-1
        
        allocate(a_nm1(nm_in),b_nm1(nm_in),c_nm1(nm_in),d_nm1(nm_in))
        call mie_coeff_nm(N_in, real(k)*sphere(sph)%r, sqrt(sphere(sph)%eps_r), a_nm1, b_nm1, c_nm1, d_nm1)
        
        ! Translation (in to in)
        allocate(a_out(nm_in),b_out(nm_in))
        call translate(sphere(sph)%cp, Nmax, N_in, k, a_in, b_in, a_out, b_out,0)

        b_vec(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) =  a_out * b_nm1 
        b_vec(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) =  b_out * a_nm1


        deallocate(a_nm1,b_nm1,c_nm1,d_nm1)
        deallocate(a_out, b_out)
        
    end do

end subroutine cluster_rhs

!****************************************************************************

subroutine matvec(sphere, k, x, x_out)
    type (data_struct), dimension(:), intent(in) :: sphere
    complex(dp), dimension(:), intent(in) :: x 
    complex(dp), dimension(:), intent(out) :: x_out
    complex(dp), intent(in) :: k

    integer :: sph1, sph2, Nspheres, Nmax, nm, s1, Nmax2, nm2
    integer :: i1_a1, i1_a2, i1_b1, i1_b2, i2_a1, i2_a2, i2_b1, i2_b2
    complex(dp), dimension(:), allocatable :: a_loc, b_loc, a_nm, b_nm, a_out, b_out, c_nm, d_nm

    Nspheres = size(sphere)

    do sph1 = 1, Nspheres

        Nmax = sphere(sph1)%Nmax
        nm = (Nmax+1)**2-1

        allocate(a_nm(nm),b_nm(nm), c_nm(nm), d_nm(nm))
        call mie_coeff_nm(Nmax, real(k)*sphere(sph1)%r, sqrt(sphere(sph1)%eps_r), a_nm, b_nm, c_nm, d_nm)

        
        s1 = sphere(sph1)%ind_a2 - sphere(sph1)%ind_a1 + 1
        allocate(a_loc(s1), b_loc(s1))
        allocate(a_out(s1), b_out(s1))
        a_loc(:) = dcmplx(0.0,0.0)
        b_loc(:) = dcmplx(0.0,0.0)

        i1_a1 = sphere(sph1)%ind_a1
        i1_a2 = sphere(sph1)%ind_a2
        i1_b1 = sphere(sph1)%ind_b1
        i1_b2 = sphere(sph1)%ind_b2


        do sph2 = 1, Nspheres 

            i2_a1 = sphere(sph2)%ind_a1
            i2_a2 = sphere(sph2)%ind_a2
            i2_b1 = sphere(sph2)%ind_b1
            i2_b2 = sphere(sph2)%ind_b2

            Nmax2 = sphere(sph2)%Nmax
            nm2 = (Nmax2+1)**2-1

            if(sph1 == sph2) then
                !a_loc = a_loc + x(i2_a1:i2_a2)
                !b_loc = b_loc + x(i2_b1:i2_b2)
            else 

                call translate(sphere(sph1)%cp-sphere(sph2)%cp, Nmax2, Nmax, k, &
                    x(i2_a1:i2_a2), x(i2_b1:i2_b2), a_out, b_out,1)

                a_loc = a_loc + b_nm * a_out  
                b_loc = b_loc + a_nm * b_out
            end if

        end do

        
        x_out(i1_a1:i1_a2) = a_loc
        x_out(i1_b1:i1_b2) = b_loc


        deallocate(a_out, b_out)
        deallocate(a_loc, b_loc)
        deallocate(a_nm, b_nm,c_nm,d_nm)
    end do

end subroutine matvec

subroutine matvecI(sphere, k, x, x_out)
    type (data_struct), dimension(:), intent(inout) :: sphere
    complex(dp), dimension(:), intent(in)   :: x 
    complex(dp), dimension(:), intent(out)  :: x_out
    complex(dp), intent(in)                 :: k
    call matvec(sphere, k, x, x_out)
    x_out =  x - x_out 
end subroutine matvecI


subroutine gmres(sphere, kk, b, x, tol, restart, maxit)
    type (data_struct), dimension(:) :: sphere
    complex(dp), dimension(:) :: b, x
    real(dp) :: tol
    complex(dp) :: kk
    integer :: restart, maxit

    complex(dp), dimension(:), allocatable :: r, w, cs, sn, g, y, Ax, mx 
    complex(dp), dimension(:,:), allocatable :: v, h

    integer :: N, max_iter, k, j, i, iter, m, ite, T1, T2, rate
    real(dp) :: err_tol, b_norm, error, nu, normav, normav2, res_norm
    complex(dp) :: temp(2), tmp, hr

    N = size(b)

    err_tol = tol
    max_iter = restart ! restart number
    m = maxit ! number of iterations / restarts

    allocate(Ax(N),mx(N))


    mx = b ! Initial guess
    x=mx

    b_norm = dble(sqrt(dot_product(b,b)))

    allocate(y(m))
    
    allocate(r(N), w(N))
    allocate(v(N,m+1))
    allocate(h(m+1,m))
    allocate(cs(m+1), sn(m+1), g(m+1))



    v(:,:) = dcmplx(0.0, 0.0)
    h(:,:) = dcmplx(0.0, 0.0)

    cs(:) = dcmplx(0.0, 0.0)
    sn(:) = dcmplx(0.0, 0.0)

    w(:) = dcmplx(0.0, 0.0)
    ! GMRES ITERATIONS
    ite = 0


    do iter = 1,max_iter
    

        mx = x
        call matvecI(sphere, kk, mx, Ax) 
        r = b-Ax

        res_norm = dble(sqrt(dot_product(r,r)))
        if (res_norm < err_tol*1e-8) then
            exit
        end if

        v(:,1) = r / res_norm
        g(:) = dcmplx(0.0,0.0)
        g(1) = res_norm
        
        do i = 1, m
            call system_clock(T1,rate)
            
            mx = v(:,i)
            
            !call matvecI(sphere, kk, mx, Ax)
            call matvecI(sphere, kk, mx, Ax)
            w = Ax
            
            normav = dble(sqrt(dot_product(Ax, Ax)))

            
            !_______Modified Gram-Schmidt____________________________
            do k = 1,i
                h(k,i) = dot_product(v(:,k),w)
                w = w - h(k,i)*v(:,k)
            end do
            
            h(i+1,i) = sqrt(dot_product(w, w))
            normav2 = dble(h(i+1,i))
            v(:,i+1) = w
            
            !_____________Reorthogonalize?________________________________
            if(normav + 0.001*normav2 == normav) then
                do j = 1,i
                    hr = dot_product(v(:,j), v(:,i+1))
                    h(j,i) = h(j,i) + hr
                    v(:,i+1) = v(:,i+1) - hr*v(:,j)
                end do
                
                h(i+1,i) = sqrt(dot_product(v(:,i+1), v(:,i+1)))
                print*, 'reorthogonalize'
            end if
            !______________________________________________________

            if(h(i+1,i) .ne. 0.0) then
                v(:,i+1) = v(:,i+1) / h(i+1,i)
            end if
            

            !_____ apply Givens rotations_________________________________
            if(i>1) then        
                do k = 1,i-1                
                    tmp = cs(k)*h(k,i) - sn(k)*h(k+1,i)
                    h(k+1,i) = sn(k)*h(k,i) + conjg(cs(k))*h(k+1,i)  
                    h(k,i) = tmp
                end do
            end if
            !________________________________________________
            nu = dble(sqrt(dot_product(H(i:i+1,i), H(i:i+1,i))))
            
            if(nu .ne. 0.0) then 
                
                cs(i) = conjg(h(i,i)/nu) 
                sn(i) = -h(i+1,i)/nu  
                H(i,i) = cs(i)*H(i,i) - sn(i)*H(i+1,i);
                H(i+1,i) = 0.0;
                temp(1:2) = g(i:i+1)
                g(i) = cs(i)*temp(1) - sn(i)*temp(2)
                g(i+1) = sn(i)*temp(1) + conjg(cs(i))*temp(2)
                
            end if
            
            error  = abs(g(i+1)) / b_norm;
            
            
            if(error < err_tol) then
                y = matmul(Cinv(H(1:i,1:i)) , g(1:i));        
                x = x + matmul(V(:,1:i),y) 
                exit       
            end if
        
            call system_clock(T2)
            
            
            ! print *,'RE (',ite+1,')','=', real(error), 'time/iter =',real(T2-T1) / real(rate)
            
            ite = ite + 1
            
        end do

        if (error < err_tol) then
            exit
        end if
        
        y = matmul(Cinv(H(1:m,1:m)) , g(1:m));
        x = x + matmul(V(:,1:m),y)
        mx = x

        call matvecI(sphere, kk, mx, Ax)
        r = b - Ax
        
        !  print *, 'Restart GMRES',abs(sqrt(dot_product(r,r))) / b_norm

        if (error < err_tol) then
            exit
        end if
    
    end do

    !print*, 'GMRES converged in', ite, 'iterations'
end subroutine gmres


subroutine cluster_coeff(sphere, x, k, a_nm, b_nm, N_out)
    complex(dp), dimension((N_out+1)**2-1), intent(out) :: a_nm, b_nm
    type (data_struct), dimension(:), intent(in)    :: sphere
    complex(dp), dimension(:), intent(in)           :: x
    complex(dp), intent(in)                         :: k
    integer, intent(in)                             :: N_out

    complex(dp), dimension(:), allocatable :: a_out, b_out
    integer :: N_in, nm_in, nm_out, Nspheres, a1,a2,b1,b2, sph


    Nspheres = size(sphere)

    nm_out = (N_out+1)**2-1

    allocate(a_out(nm_out),b_out(nm_out))

    a_nm(:) = dcmplx(0.0,0.0)
    b_nm(:) = dcmplx(0.0,0.0)


    do sph = 1, Nspheres
        a1 = sphere(sph)%ind_a1
        a2 = sphere(sph)%ind_a2
        b1 = sphere(sph)%ind_b1
        b2 = sphere(sph)%ind_b2

        N_in = sphere(sph)%Nmax
        nm_in = (N_in+1)**2-1

        call translate(-sphere(sph)%cp, N_in, N_out, k, x(a1:a2),x(b1:b2), a_out, b_out,0)

        a_nm = a_nm + a_out
        b_nm = b_nm + b_out

    end do

end subroutine cluster_coeff


!**********************************


!****************************************************************************
! Compute right hand side for a planewave incident

subroutine inc_xy(sphere, b_vecx, b_vecy, k)
    type (data_struct), dimension(:), intent(inout)        :: sphere
    complex(dp), dimension(:), allocatable, intent(out) :: b_vecx, b_vecy
    complex(dp), intent(in) :: k

    complex(dp), dimension(:,:), allocatable :: a_out, b_out
    complex(dp) :: phase_shift
    integer ::  Nspheres, n, loc, sph, bsize, las, nm_in, N_in, Nmax2
    complex(dp), dimension(:), allocatable :: rotD
    integer, dimension(:,:), allocatable :: indD
    real(dp) :: delta, cp(3)

    Nspheres = size(sphere)

    bsize = 0
    do sph = 1, Nspheres
        Nmax2 = sphere(sph)%Nmax
        loc = (Nmax2+1)**2 - 1
        sphere(sph)%ind_a1 = bsize + 1 
        sphere(sph)%ind_a2 = bsize + loc
        sphere(sph)%ind_b1 = bsize + loc + 1 
        sphere(sph)%ind_b2 = bsize + 2*loc  
        bsize = bsize + 2*loc
    end do

    allocate(b_vecx(bsize), b_vecy(bsize))
    b_vecx(:) = dcmplx(0.0,0.0)
    b_vecy(:) = dcmplx(0.0,0.0)



    !$omp parallel default(private) &
    !$omp firstprivate(k, Nspheres, Nmax, a_in,b_in) &
    !$omp shared(b_vecx,b_vecy,sphere)
    !$omp do
    do sph = 1, Nspheres
        N_in = sphere(sph)%Nmax
        nm_in = (N_in+1)**2-1

        allocate(a_out(nm_in,2),b_out(nm_in,2))
        call planewave(N_in, real(k), a_out(:,1), b_out(:,1))

        las = 0
        do n = 1,N_in
            las = las + (2*n+1)**2
        end do

        allocate(rotD(las))
        allocate(indD(las,2))

        call sph_rotation_sparse(0.0d0, -pi/2.0d0, N_in, rotD, indD)
        a_out(:,2) = sparse_matmul(rotD,indD,a_out(:,1),nm_in)
        b_out(:,2) = sparse_matmul(rotD,indD,b_out(:,1),nm_in)


        cp = sphere(sph)%cp
        delta = cp(3)
        phase_shift = cdexp(dcmplx(0.0, dble(k)*delta))

        a_out = a_out * phase_shift
        b_out = b_out * phase_shift
        

        b_vecx(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) =   a_out(:,1)
        b_vecx(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) =   b_out(:,1)

        b_vecy(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) =   a_out(:,2)
        b_vecy(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) =   b_out(:,2)

        
        deallocate(a_out, b_out)
        deallocate(rotD,indD)

    end do
    !$omp end do
    !$omp end parallel

end subroutine inc_xy


!****************************************************************************

subroutine rhs2_xy(sphere, Tmat, b_vecx, b_vecy, k)
    type (data_struct), dimension(:), intent(inout) :: sphere
    type (Tmatrix), dimension(:), intent(in) :: Tmat
    complex(dp), intent(in) :: k
    complex(dp), dimension(:), allocatable, intent(out) :: b_vecx, b_vecy

    complex(dp), dimension(:), allocatable :: a_nm1, b_nm1
    complex(dp), dimension(:,:), allocatable :: a_out, b_out, a_out2, b_out2
    complex(dp), dimension(:), allocatable :: c_nm1, d_nm1
    complex(dp) :: phase_shift
    integer ::  Nspheres, n, loc, sph, bsize, las, nm_in, N_in, Nmax2, ind
    complex(dp), dimension(:,:), allocatable :: rotD2
    complex(dp), dimension(:), allocatable :: rotD
    integer, dimension(:,:), allocatable :: indD
    real(dp) :: delta, cp(3)

    Nspheres = size(sphere)

    bsize = 0
    do sph = 1, Nspheres
        Nmax2 = sphere(sph)%Nmax
        loc = (Nmax2+1)**2 - 1
        sphere(sph)%ind_a1 = bsize + 1 
        sphere(sph)%ind_a2 = bsize + loc
        sphere(sph)%ind_b1 = bsize + loc + 1 
        sphere(sph)%ind_b2 = bsize + 2*loc  
        bsize = bsize + 2*loc
    end do

    allocate(b_vecx(bsize), b_vecy(bsize))
    b_vecx(:) = dcmplx(0.0,0.0)
    b_vecy(:) = dcmplx(0.0,0.0)



    !$omp parallel default(private) &
    !$omp firstprivate(k, Nspheres, Nmax, a_in,b_in) &
    !$omp shared(b_vecx,b_vecy,sphere, Tmat)
    !$omp do
    do sph = 1, Nspheres
        N_in = sphere(sph)%Nmax
        nm_in = (N_in+1)**2-1

        allocate(a_nm1(nm_in),b_nm1(nm_in),c_nm1(nm_in),d_nm1(nm_in))
        call mie_coeff_nm(N_in, real(k)*sphere(sph)%r, sqrt(sphere(sph)%eps_r), a_nm1, b_nm1, c_nm1, d_nm1)



        allocate(a_out(nm_in,2),b_out(nm_in,2))

        call planewave(N_in, real(k), a_out(:,1), b_out(:,1))

        las = 0
        do n = 1,N_in
            las = las + (2*n+1)**2
        end do

        allocate(rotD(las))
        allocate(indD(las,2))

        call sph_rotation_sparse(0.0d0, -pi/2.0d0, N_in, rotD, indD)
        a_out(:,2) = sparse_matmul(rotD,indD,a_out(:,1),nm_in)
        b_out(:,2) = sparse_matmul(rotD,indD,b_out(:,1),nm_in)


        cp = sphere(sph)%cp
        delta = cp(3)
        phase_shift = cdexp(dcmplx(0.0, dble(k)*delta))

        a_out = a_out * phase_shift
        b_out = b_out * phase_shift
        
        ind = sphere(sph)%Tmat_ind

        if(size(Tmat) == 0) then
            b_vecx(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) =  b_nm1 * a_out(:,1)
            b_vecx(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) =  a_nm1 * b_out(:,1)

            b_vecy(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) =  b_nm1 * a_out(:,2)
            b_vecy(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) =  a_nm1 * b_out(:,2)


        else 

            allocate(a_out2(nm_in,2),b_out2(nm_in,2))


            allocate(rotD2(nm_in,nm_in))
            call sph_rotation_gen(sphere(sph)%euler_angles, N_in, rotD2)

            call rot_Tmatmul(Tmat(ind)%Taa, Tmat(ind)%Tab, Tmat(ind)%Tba, &
                Tmat(ind)%Tbb, a_out, b_out, rotD2, nm_in, a_out2, b_out2)
            
            b_vecx(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) = a_out2(:,1)
            b_vecx(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) = b_out2(:,1)

            b_vecy(sphere(sph)%ind_a1 : sphere(sph)%ind_a2) = a_out2(:,2)          
            b_vecy(sphere(sph)%ind_b1 : sphere(sph)%ind_b2) = b_out2(:,2)

            deallocate(a_out2, b_out2)
            deallocate(rotD2)

        end if



        deallocate(a_nm1,b_nm1,c_nm1,d_nm1)
        deallocate(a_out, b_out)
        deallocate(rotD,indD)

    end do
    !$omp end do
    !$omp end parallel

end subroutine rhs2_xy

!************************************************************************
!
!  MLFMM accelerated matvec multiplication (2 rhs)
!
!*************************************************************************

subroutine matvec_mlfmm2(sphere, otree, k, x, x_out, x2, x_out2)
    use iso_fortran_env, only : int64
    type (data_struct), dimension(:), intent(in) :: sphere
    type (level_struct), dimension(:), intent(in) :: otree
    complex(dp), dimension(:), intent(in) :: x , x2
    complex(dp), dimension(:), intent(out) :: x_out, x_out2
    complex(dp), intent(in) :: k

    type (data), dimension(:), allocatable :: aggregates, aggregates2

    integer :: sph1, sph2, Nspheres, Nmax, nm, s1, Nmax2, nm2
    integer :: i1_a1, i1_a2, i1_b1, i1_b2, i2_a1, i2_a2, i2_b1, i2_b2, box, box2
    integer :: sph1_loc, sph2_loc, box2_loc, i1, near_box, Nmax_box
    integer :: Nmax_box1, Nmax_box2
    complex(dp), dimension(:,:), allocatable :: a_loc, b_loc, a_out, b_out
    complex(dp), dimension(:,:), allocatable :: a_box2, b_box2, a_box2_2, b_box2_2
    complex(dp), dimension(:,:), allocatable :: a_box_loc, b_box_loc, a_box2_loc, b_box2_loc
    real(dp) :: cp_box(3), cp_box2(3)
    integer :: nm_box, nm_box2, max_level, level_ind, level, parent, neighbour, i2, i3
    integer :: T1, T2, rate
    integer(int64) :: mem
    call system_clock(T1,rate)

    max_level = size(otree) - 1
    level_ind = max_level + 1
    Nspheres = size(sphere)

    !*** Allocate aggregates *******

    allocate(aggregates(max_level))
    allocate(aggregates2(max_level))


    Nmax_box = otree(max_level+1)%tree(1)%Nmax   
    nm_box = (Nmax_box+1)**2-1

    allocate(aggregates(max_level)%a_out(nm_box,size(otree(max_level+1)%tree)))
    allocate(aggregates(max_level)%b_out(nm_box,size(otree(max_level+1)%tree)))

    allocate(aggregates2(max_level)%a_out(nm_box,size(otree(max_level+1)%tree)))
    allocate(aggregates2(max_level)%b_out(nm_box,size(otree(max_level+1)%tree)))

    mem = sizeof(aggregates(max_level)%a_out)

    if(max_level >= 3) then
    
        do level = max_level, 3, -1
            Nmax_box2 = otree(level)%tree(1)%Nmax  ! low
            nm_box = (Nmax_box2+1)**2-1 !low

            allocate(aggregates(level-1)%a_out(nm_box,size(otree(level)%tree)))  
            allocate(aggregates(level-1)%b_out(nm_box,size(otree(level)%tree))) 

            allocate(aggregates2(level-1)%a_out(nm_box,size(otree(level)%tree)))  
            allocate(aggregates2(level-1)%b_out(nm_box,size(otree(level)%tree))) 

            mem = mem + sizeof(aggregates(level-1)%a_out)

        end do
    end if



    !*************** Near-zone interactions *******************************

    !$omp parallel default(private) &
    !$omp firstprivate(k, level_ind) &
    !$omp shared(otree,sphere, x_out,x_out2,x,x2) &
    !$omp shared(aggregates,aggregates2)
    !$omp do schedule(dynamic,1)
    do box = 1, size(otree(level_ind)%tree)
        
        do sph1_loc = 1,  otree(level_ind)%tree(box)%N_source

            sph1 =  otree(level_ind)%tree(box)%sources(sph1_loc)
            Nmax = sphere(sph1)%Nmax
            nm = (Nmax+1)**2-1

            !allocate(a_nm(nm),b_nm(nm), c_nm(nm), d_nm(nm))
            !call mie_coeff_nm(Nmax, real(k)*sphere(sph1)%r, sqrt(sphere(sph1)%eps_r), a_nm, b_nm, c_nm, d_nm)

        
            s1 = sphere(sph1)%ind_a2 - sphere(sph1)%ind_a1 + 1
            allocate(a_loc(s1,2), b_loc(s1,2))
            allocate(a_out(s1,2), b_out(s1,2))
            a_loc(:,:) = dcmplx(0.0,0.0)
            b_loc(:,:) = dcmplx(0.0,0.0)

            i1_a1 = sphere(sph1)%ind_a1
            i1_a2 = sphere(sph1)%ind_a2
            i1_b1 = sphere(sph1)%ind_b1
            i1_b2 = sphere(sph1)%ind_b2


            do box2_loc = 1,27

                box2 =   otree(level_ind)%tree(box)%near_neighbours(box2_loc)

                
                if(box2 == 0) exit
            
                do sph2_loc = 1,   otree(level_ind)%tree(box2)%N_source

                    sph2 =  otree(level_ind)%tree(box2)%sources(sph2_loc)

                    i2_a1 = sphere(sph2)%ind_a1
                    i2_a2 = sphere(sph2)%ind_a2
                    i2_b1 = sphere(sph2)%ind_b1
                    i2_b2 = sphere(sph2)%ind_b2

                    Nmax2 = sphere(sph2)%Nmax
                    nm2 = (Nmax2+1)**2-1

                    if(sph1 .ne. sph2) then
                

                    call translate_xy(sphere(sph1)%cp-sphere(sph2)%cp, Nmax2, Nmax, k, &
                            [x(i2_a1:i2_a2),x2(i2_a1:i2_a2)], &
                            [x(i2_b1:i2_b2),x2(i2_b1:i2_b2)], a_out, b_out,1)

                    !a_loc(:,1) = a_loc(:,1) + b_nm * a_out(:,1)  
                    !b_loc(:,1) = b_loc(:,1) + a_nm * b_out(:,1)

                    !a_loc(:,2) = a_loc(:,2) + b_nm * a_out(:,2)  
                    !b_loc(:,2) = b_loc(:,2) + a_nm * b_out(:,2)

                    a_loc(:,1) = a_loc(:,1) + a_out(:,1)  
                    b_loc(:,1) = b_loc(:,1) + b_out(:,1)

                    a_loc(:,2) = a_loc(:,2) + a_out(:,2)  
                    b_loc(:,2) = b_loc(:,2) + b_out(:,2)

                    end if

                end do

            end do

            x_out(i1_a1:i1_a2) = a_loc(:,1)
            x_out(i1_b1:i1_b2) = b_loc(:,1)

            x_out2(i1_a1:i1_a2) = a_loc(:,2)
            x_out2(i1_b1:i1_b2) = b_loc(:,2)


            deallocate(a_out, b_out)
            deallocate(a_loc, b_loc)
            !deallocate(a_nm, b_nm,c_nm,d_nm)
        end do
        
    end do

    !$omp end do
    !$omp end parallel

    call system_clock(T2)

    if(sphere(1)%ifT == 0) print*, 'Near-zone done in ', real(T2-T1) / real(rate)

    ! *******************************************************
    !
    ! mlfmm part
    !
    !*******************************************************
    call system_clock(T1,rate)

    !**** Aggregation at the finest level ****************



    Nmax_box = otree(max_level+1)%tree(1)%Nmax   
    nm_box = (Nmax_box+1)**2-1

    allocate(a_box_loc(nm_box,2), b_box_loc(nm_box,2))

    aggregates(max_level)%a_out(:,:) = dcmplx(0.0,0.0)
    aggregates(max_level)%b_out(:,:) = dcmplx(0.0,0.0)

    aggregates2(max_level)%a_out(:,:) = dcmplx(0.0,0.0)
    aggregates2(max_level)%b_out(:,:) = dcmplx(0.0,0.0)



    !$omp parallel default(private) &
    !$omp firstprivate(k, max_level, Nmax_box) &
    !$omp shared(otree, aggregates, aggregates2, sphere, x, x2)
    !$omp do schedule(dynamic,1)
    do box = 1, size(otree(max_level+1)%tree)

        cp_box = otree(max_level+1)%tree(box)%cp

        
        do sph1_loc = 1, otree(max_level+1)%tree(box)%N_source

            sph1 = otree(max_level+1)%tree(box)%sources(sph1_loc)
        
            Nmax = sphere(sph1)%Nmax
            nm = (Nmax+1)**2-1

            i1_a1 = sphere(sph1)%ind_a1
            i1_a2 = sphere(sph1)%ind_a2
            i1_b1 = sphere(sph1)%ind_b1
            i1_b2 = sphere(sph1)%ind_b2
            
            call translate_xy(cp_box-sphere(sph1)%cp, Nmax, Nmax_box, k, &
                    [x(i1_a1:i1_a2),x2(i1_a1:i1_a2)], &
                    [x(i1_b1:i1_b2),x2(i1_b1:i1_b2)], a_box_loc, b_box_loc,0)
                
            aggregates(max_level)%a_out(:,box) = aggregates(max_level)%a_out(:,box) + a_box_loc(:,1)
            aggregates(max_level)%b_out(:,box) = aggregates(max_level)%b_out(:,box) + b_box_loc(:,1)
            
            aggregates2(max_level)%a_out(:,box) = aggregates2(max_level)%a_out(:,box) + a_box_loc(:,2)
            aggregates2(max_level)%b_out(:,box) = aggregates2(max_level)%b_out(:,box) + b_box_loc(:,2)
            


        end do

    end do
    !$omp end do
    !$omp end parallel

    deallocate(a_box_loc,b_box_loc)

    !******** Aggregation to lower levels ****************!
    ! high to low 
    !
    if(max_level >= 3) then
    
        do level = max_level, 3, -1

            Nmax_box1 = otree(level)%tree(1)%Nmax ! parent
            Nmax_box2 = otree(level+1)%tree(1)%Nmax ! child 

            nm_box = (Nmax_box1+1)**2-1 

            allocate(a_box_loc(nm_box,2),b_box_loc(nm_box,2)) !low
            
            
            aggregates(level-1)%a_out(:,:) = dcmplx(0.0,0.0) ! low
            aggregates(level-1)%b_out(:,:) = dcmplx(0.0,0.0) ! low

            aggregates2(level-1)%a_out(:,:) = dcmplx(0.0,0.0) ! low
            aggregates2(level-1)%b_out(:,:) = dcmplx(0.0,0.0) ! low

            
            !$omp parallel default(private) &
            !$omp firstprivate(Nmax_box2, k, level) &
            !$omp shared(otree, aggregates, aggregates2)
            !$omp do schedule(dynamic,1) 
            
            do box = 1, size(otree(level)%tree) ! high parent box

                cp_box = otree(level)%tree(box)%cp 
                Nmax_box1 = otree(level)%tree(box)%Nmax

                do i1 = 1, 8 ! children

                    box2 = otree(level)%tree(box)%children(i1)   
                    if(box2 == 0) exit
                    cp_box2 = otree(level+1)%tree(box2)%cp 
                
                    
                    ! high to low translation 
                    call translate_xy(cp_box-cp_box2, Nmax_box2, Nmax_box1, k, &
                    [aggregates(level)%a_out(:,box2),aggregates2(level)%a_out(:,box2)], &
                    [aggregates(level)%b_out(:,box2),aggregates2(level)%b_out(:,box2)], &
                    a_box_loc, b_box_loc,0)
                    
                    aggregates(level-1)%a_out(:,box) = aggregates(level-1)%a_out(:,box) + a_box_loc(:,1)
                    aggregates(level-1)%b_out(:,box) = aggregates(level-1)%b_out(:,box) + b_box_loc(:,1)
        
                    aggregates2(level-1)%a_out(:,box) = aggregates2(level-1)%a_out(:,box) + a_box_loc(:,2)
                    aggregates2(level-1)%b_out(:,box) = aggregates2(level-1)%b_out(:,box) + b_box_loc(:,2)

                end do
            end do

            !$omp end do
            !$omp end parallel

            deallocate(a_box_loc, b_box_loc)
            
            
        end do
    end if

    call system_clock(T2)
    if(sphere(1)%ifT == 0) print*, 'Aggregation done in', real(T2-T1) / real(rate)
    !********** Translation **************************************
    do level = max_level,2,-1
        call system_clock(T1,rate)
        Nmax_box = otree(level+1)%tree(1)%Nmax   
        nm_box = (Nmax_box+1)**2-1
        allocate(a_box2_loc(nm_box,2), b_box2_loc(nm_box,2))

        allocate(a_box2(nm_box, size(otree(level+1)%tree)))
        allocate(b_box2(nm_box, size(otree(level+1)%tree)))
        
        allocate(a_box2_2(nm_box, size(otree(level+1)%tree)))
        allocate(b_box2_2(nm_box, size(otree(level+1)%tree)))


        a_box2(:,:) = dcmplx(0.0,0.0)
        b_box2(:,:) = dcmplx(0.0,0.0)

        a_box2_2(:,:) = dcmplx(0.0,0.0)
        b_box2_2(:,:) = dcmplx(0.0,0.0)

        !$omp parallel default(private) &
        !$omp firstprivate(level, Nmax_box, k) &
        !$omp shared(a_box2, b_box2, otree, aggregates) &
        !$omp shared(a_box2_2, b_box2_2, aggregates2)
        !$omp do schedule(dynamic,1) 
        do box = 1, size(otree(level+1)%tree) 
            parent = otree(level+1)%tree(box)%parent ! level - 1 
            
            cp_box =  otree(level+1)%tree(box)%cp ! from

            do i1 = 1, 27 ! near neighbours of parent cube
                neighbour = otree(level)%tree(parent)%near_neighbours(i1)
                if(neighbour == 0) exit
            
                do i2 = 1, 8 ! children of near neighbours
                
                    box2 = otree(level)%tree(neighbour)%children(i2) ! level 
                    
                    if(box2 == 0) exit
                
                    ! ***check if near box ****
                    near_box = 0
                    do i3 = 1,27
                    if(box2 == otree(level+1)%tree(box)%near_neighbours(i3)) then
                        near_box = 1
                    end if
                    end do
                
                    if(near_box == 0) then
            
                    cp_box2 =  otree(level+1)%tree(box2)%cp !to

                    !**** Translate box to accessible box
                    call translate_xy(cp_box-cp_box2, Nmax_box, Nmax_box, k, &
                            [aggregates(level)%a_out(:,box2),aggregates2(level)%a_out(:,box2)], &
                            [aggregates(level)%b_out(:,box2),aggregates2(level)%b_out(:,box2)], &
                            a_box2_loc, b_box2_loc,1)

                    a_box2(:,box) = a_box2(:,box) + a_box2_loc(:,1)
                    b_box2(:,box) = b_box2(:,box) + b_box2_loc(:,1)

                    a_box2_2(:,box) = a_box2_2(:,box) + a_box2_loc(:,2)
                    b_box2_2(:,box) = b_box2_2(:,box) + b_box2_loc(:,2)


                    end if

                end do

            end do

        end do
        
        !$omp end do
        !$omp end parallel

        aggregates(level)%a_out = a_box2
        aggregates(level)%b_out = b_box2

        aggregates2(level)%a_out = a_box2_2
        aggregates2(level)%b_out = b_box2_2

        
        deallocate(a_box2_loc,b_box2_loc)
        deallocate(a_box2,b_box2, a_box2_2, b_box2_2)

        call system_clock(T2)
        if(sphere(1)%ifT == 0) print*, 'Translation done in', real(T2-T1) / real(rate), 'level',level
    end do



    !***********  Disagregation  *********************************************
    ! 
    ! low to high 
    call system_clock(T1,rate)

    if(max_level >= 3) then
        do level = 2, max_level-1
        
            Nmax_box2 = otree(level+2)%tree(1)%Nmax ! high
            nm_box2 = (Nmax_box2+1)**2-1 ! high
            allocate(a_box_loc(nm_box2,2), b_box_loc(nm_box2,2)) ! high 

            !$omp parallel default(private) &
            !$omp firstprivate(Nmax_box2, k, level) &
            !$omp shared(otree, aggregates,aggregates2)
            !$omp do schedule(dynamic,1)
            do box = 1, size(otree(level+1)%tree) ! low
                Nmax_box1 = otree(level+1)%tree(box)%Nmax ! low
                cp_box = otree(level+1)%tree(box)%cp ! from low

                do i1 = 1, 8
                    box2 = otree(level+1)%tree(box)%children(i1) ! high
                    if(box2 == 0) exit
                    cp_box2 = otree(level+2)%tree(box2)%cp ! to high

                    ! translate low to high
                    call translate_xy(cp_box2-cp_box, Nmax_box1, Nmax_box2, k, &
                        [aggregates(level)%a_out(:,box),aggregates2(level)%a_out(:,box)], &
                        [aggregates(level)%b_out(:,box),aggregates2(level)%b_out(:,box)], &
                        a_box_loc, b_box_loc,0)
                
                    aggregates(level+1)%a_out(:,box2) = & 
                        aggregates(level+1)%a_out(:,box2) + a_box_loc(:,1)
                    
                    aggregates(level+1)%b_out(:,box2) = &
                        aggregates(level+1)%b_out(:,box2) + b_box_loc(:,1)
                    
                    aggregates2(level+1)%a_out(:,box2) = & 
                        aggregates2(level+1)%a_out(:,box2) + a_box_loc(:,2)
                    
                    aggregates2(level+1)%b_out(:,box2) = &
                        aggregates2(level+1)%b_out(:,box2) + b_box_loc(:,2)
                    
                    
                end do
                
            end do
            !$omp end do
            !$omp end parallel

            deallocate(a_box_loc,b_box_loc)
        end do
    end if

    !************   Dissaggregate   ******************

    !$omp parallel default(private) &
    !$omp firstprivate(k, max_level) &
    !$omp shared(otree, aggregates, aggregates2, sphere, x_out,x_out2)
    !$omp do schedule(dynamic,1)

    do box = 1, size(otree(max_level+1)%tree)

        Nmax_box = otree(max_level+1)%tree(box)%Nmax
        cp_box = otree(max_level+1)%tree(box)%cp

    
        do sph1_loc = 1, otree(max_level+1)%tree(box)%N_source

            sph1 = otree(max_level+1)%tree(box)%sources(sph1_loc)
        
            Nmax = sphere(sph1)%Nmax
            nm = (Nmax+1)**2-1

            i1_a1 = sphere(sph1)%ind_a1
            i1_a2 = sphere(sph1)%ind_a2
            i1_b1 = sphere(sph1)%ind_b1
            i1_b2 = sphere(sph1)%ind_b2
            
            
            !allocate(a_nm(nm),b_nm(nm), c_nm(nm), d_nm(nm))
            !call mie_coeff_nm(Nmax, real(k)*sphere(sph1)%r, sqrt(sphere(sph1)%eps_r), a_nm, b_nm, c_nm, d_nm)
                    
            s1 = sphere(sph1)%ind_a2 - sphere(sph1)%ind_a1 + 1           
            allocate(a_out(s1,2), b_out(s1,2))
                    
            call translate_xy(sphere(sph1)%cp-cp_box, Nmax_box, Nmax, k, &
                    [aggregates(max_level)%a_out(:,box),aggregates2(max_level)%a_out(:,box)], &
                    [aggregates(max_level)%b_out(:,box),aggregates2(max_level)%b_out(:,box)],&
                    a_out, b_out,0)
            
            !x_out(i1_a1:i1_a2) = x_out(i1_a1:i1_a2) + b_nm * a_out(:,1)
            !x_out(i1_b1:i1_b2) =  x_out(i1_b1:i1_b2) + a_nm * b_out(:,1)

            !x_out2(i1_a1:i1_a2) = x_out2(i1_a1:i1_a2) + b_nm * a_out(:,2)
            !x_out2(i1_b1:i1_b2) = x_out2(i1_b1:i1_b2) + a_nm * b_out(:,2)

            x_out(i1_a1:i1_a2) = x_out(i1_a1:i1_a2) + a_out(:,1)
            x_out(i1_b1:i1_b2) =  x_out(i1_b1:i1_b2) + b_out(:,1)

            x_out2(i1_a1:i1_a2) = x_out2(i1_a1:i1_a2) + a_out(:,2)
            x_out2(i1_b1:i1_b2) = x_out2(i1_b1:i1_b2) + b_out(:,2)
            

            !deallocate(a_nm,b_nm,c_nm,d_nm)
            deallocate(a_out, b_out) 

        end do

    end do

    !$omp end do
    !$omp end parallel

    call system_clock(T2)
    if(sphere(1)%ifT == 0) print*, 'Disaggregation done in', real(T2-T1) / real(rate)
    deallocate(aggregates,aggregates2)
end subroutine matvec_mlfmm2



!***********************************************************************
subroutine matvecI_mlfmm2(sphere, otree, Tmat, k, x, x_out, x2, x_out2)
    type (data_struct), dimension(:), intent(in) :: sphere
    type (level_struct), dimension(:), intent(in) :: otree
    type (Tmatrix), dimension(:), intent(in) :: Tmat
    complex(dp), dimension(:), intent(in) :: x , x2
    complex(dp), dimension(:), intent(out) :: x_out, x_out2
    complex(dp), intent(in) :: k



    integer :: sph, Nmax, nm, a1,a2, b1, b2, ind, las, n
    complex(dp), dimension(:), allocatable :: a_nm,b_nm,c_nm,d_nm
    complex(dp), dimension(:,:), allocatable :: a_in,b_in,a_out,b_out, rotD2
    complex(dp), dimension(:), allocatable :: rotD
    integer, dimension(:,:), allocatable :: indD 

    call matvec_mlfmm2(sphere,otree, k, x, x_out, x2, x_out2)

    if(size(Tmat) == 0) then

        do sph = 1,size(sphere)
            Nmax = sphere(sph)%Nmax
            nm = (Nmax+1)**2-1
            allocate(a_nm(nm),b_nm(nm), c_nm(nm), d_nm(nm))
            call mie_coeff_nm(Nmax, real(k)*sphere(sph)%r, sqrt(sphere(sph)%eps_r), a_nm, b_nm, c_nm, d_nm)

            a1 = sphere(sph)%ind_a1
            a2 = sphere(sph)%ind_a2
            b1 = sphere(sph)%ind_b1
            b2 = sphere(sph)%ind_b2
            
            x_out(a1:a2) =  x(a1:a2) - x_out(a1:a2) * b_nm 
            x_out(b1:b2) =  x(b1:b2) - x_out(b1:b2) * a_nm 
            
            x_out2(a1:a2) =  x2(a1:a2) - x_out2(a1:a2) * b_nm
            x_out2(b1:b2) =  x2(b1:b2) - x_out2(b1:b2) * a_nm 
            
            deallocate(a_nm,b_nm,c_nm,d_nm)
        end do   
    else
    
        do sph = 1,size(sphere)
            Nmax = sphere(sph)%Nmax
            nm = (Nmax+1)**2-1
                
            a1 = sphere(sph)%ind_a1
            a2 = sphere(sph)%ind_a2
            b1 = sphere(sph)%ind_b1
            b2 = sphere(sph)%ind_b2
            
            ind = sphere(sph)%Tmat_ind
            
            !x_out(a1:a2) =  x(a1:a2) - matmul(Tmat(ind)%Taa,x_out(a1:a2)) &
            !     - matmul(Tmat(ind)%Tab, x_out(b1:b2))

            !x_out(b1:b2) =  x(b1:b2) - matmul(Tmat(ind)%Tbb, x_out(b1:b2)) &
            !     - matmul(Tmat(ind)%Tba, x_out(a1:a2))
            
            !x_out2(a1:a2) =  x2(a1:a2) - matmul(Tmat(ind)%Taa, x_out2(a1:a2)) &
            !     - matmul(Tmat(ind)%Tab, x_out2(b1:b2))

            !x_out2(b1:b2) =  x2(b1:b2) - matmul(Tmat(ind)%Tbb, x_out2(b1:b2)) &
            !     - matmul(Tmat(ind)%Tba, x_out2(a1:a2))
        

            
            las = 0
            do n = 1,Nmax
                las = las + (2*n+1)**2
            end do

            allocate(rotD(las))
            allocate(indD(las,2))

            !call sph_rotation_sparse_gen(rot, Nmax, rotD, indD)

            allocate(rotD2(nm,nm))
            call sph_rotation_gen(sphere(sph)%euler_angles, Nmax, rotD2)


            allocate(a_in(nm,2),b_in(nm,2))
            allocate(a_out(nm,2),b_out(nm,2))
            
            a_in(:,1) = x_out(a1:a2)
            b_in(:,1) = x_out(b1:b2)

            a_in(:,2) = x_out2(a1:a2)
            b_in(:,2) = x_out2(b1:b2)

            call rot_Tmatmul(Tmat(ind)%Taa, Tmat(ind)%Tab, Tmat(ind)%Tba, &
                Tmat(ind)%Tbb, a_in, b_in, rotD2, nm, a_out, b_out)

            !call Tmatmul(Tmat(ind)%Taa, Tmat(ind)%Tab, Tmat(ind)%Tba, &
            !     Tmat(ind)%Tbb, a_in, b_in, nm, a_out, b_out)

            x_out(a1:a2) =  x(a1:a2) - a_out(:,1)
            x_out(b1:b2) =  x(b1:b2) - b_out(:,1)

            x_out2(a1:a2) =  x2(a1:a2) - a_out(:,2)
            x_out2(b1:b2) =  x2(b1:b2) - b_out(:,2)

            deallocate(a_in,b_in,a_out,b_out)
            deallocate(rotD,indD)
            deallocate(rotD2)

        end do

    end if


end subroutine matvecI_mlfmm2
!***********************************************************************

end module 
