!-------------------------------------------------------------------
! Generates incoherent T-matrices and also a coherent T-matrix.
!
! Copyright (C) 2018 Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 27.3.2018 TV: Maintenance, cleaning, commenting
!
! TODO
! *More comments and some refactoring
!-----------------------------

module T_matrix
use mpi
use common
use solver
use translations
use geometry_x
use mie

implicit none 
    
contains
    

    ! Routine computes the T-matrix (of order Nmax) for a cluster of spheres. 
    ! the cluster must be centered at the origin
     
    subroutine compute_T_matrix(sphere, k, Nmax, Taa, Tab, Tba, Tbb)
        type (data_struct), dimension(:), intent(inout) :: sphere
        complex(dp), intent(in)                         :: k
        complex(dp), intent(out)                        :: Taa((Nmax + 1)**2-1, (Nmax + 1)**2-1)
        complex(dp), intent(out)                        :: Tbb((Nmax + 1)**2-1, (Nmax + 1)**2-1)
        complex(dp), intent(out)                        :: Tab((Nmax + 1)**2-1, (Nmax + 1)**2-1)
        complex(dp), intent(out)                        :: Tba((Nmax + 1)**2-1, (Nmax + 1)**2-1)
        integer, intent(in)                             :: Nmax 

        integer :: i1
        complex(dp) :: a_in((Nmax + 1)**2-1), b_in((Nmax + 1)**2-1)
        complex(dp) :: a_nm((Nmax + 1)**2-1), b_nm((Nmax + 1)**2-1)

        complex(dp), dimension(:), allocatable :: rhs, x
        
        real(dp) :: orig(3)
        
        a_in(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        b_in(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        
        orig(:) = 0.0
        do i1 = 1, size(a_in) 
            a_in(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
            b_in(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
            
            a_in(i1) = dcmplx(1.0,0.0)
            
            call cluster_rhs(Nmax, a_in, b_in, sphere, rhs, k)
            allocate(x(size(rhs)))
            
            call gmres(sphere, k, rhs, x, dble(1e-3), 10, 50) 
            call cluster_coeff(sphere, x, k, a_nm, b_nm, Nmax)
            
            Taa(:,i1) = a_nm
            Tba(:,i1) = b_nm
            
            deallocate(x, rhs)
            
            a_in(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
            b_in(:) = cmplx(0.0_dp,0.0_dp,kind=dp)
            
            b_in(i1) = dcmplx(1.0,0.0)
            
            call cluster_rhs(Nmax, a_in, b_in, sphere, rhs, k)
            allocate(x(size(rhs)))
            
            call gmres(sphere, k, rhs, x, dble(1e-3), 10, 50)
            
            
            call cluster_coeff(sphere, x, k, a_nm, b_nm, Nmax)
            
            Tab(:,i1) = a_nm
            Tbb(:,i1) = b_nm
            
            deallocate(x, rhs)
        
        end do
        
    end subroutine compute_T_matrix
    
    
    

    !******************************************************************
    !*
    !* The function computes the realization avegared T-matrix that 
    !* defines a mapping from incoming coefficient to the scattered 
    !* coefficients of the coherent field   
    !*
    !*******************************************************************
    subroutine compute_T_matrix_ave(volrad, dens, epsr, k, Nave, Taa2, &
         &  Tab2, Tba2, Tbb2, N_theta, N_phi, S_ave, kappa_s, kappa_a, kappa_s_tot, &
         &  Taa_data,Tab_data, Tba_data, Tbb_data,n_samples,cexts)
        real(dp), intent(in)                        :: volrad       !volume element radius
        real(dp), intent(in)                        :: dens         !density
        complex(dp), intent(in)                     :: epsr, k      !permittivity and wavenumber
        complex(dp), dimension(:,:), allocatable, intent(out) :: Taa2, Tab2  !Coherent T-matrix
        complex(dp), dimension(:,:), allocatable, intent(out) :: Tba2, Tbb2  !Coherent T-matrix
        complex(dp), dimension(:,:,:), allocatable, intent(out) :: Taa_data, Tab_data   !List of incoherent T-matrices
        complex(dp), dimension(:,:,:), allocatable, intent(out) :: Tba_data, Tbb_data   !List of incoherent T-matrices
        integer, intent(in)                                     :: N_theta, N_phi       !How many angles should be covered
        integer, intent(in)  :: n_samples                       !how many samples of incoherent volume elements are stored
        real(dp), dimension(:,:), allocatable, intent(out)      :: S_ave                !Averaged incoherent S-matrix
        real(dp), intent(out) :: kappa_a, kappa_s_tot,kappa_s   !kappa (absorption, coherent and incoherent scattering)
        real(dp), allocatable, intent(out) :: cexts(:,:)        !cexts (1,:) contains extinction cs, and (2,:) albedos 

        type (Tmatrix), dimension(:), allocatable :: Tmat
        type (data_struct), dimension(:), allocatable :: sphere
        real(dp), dimension(:,:), allocatable     :: S_out
        real(dp) :: ka, maxrad, vol
        real(dp), pointer :: radius(:), coord(:,:)
        real(dp) :: Csca, Csca2, Csca_ave
        real(dp) :: Cabs, Cabs2, Cabs_ave, Csca_ave_tot, crs1(4),crs2(4)
        complex(dp), dimension(:,:), allocatable :: Taa, Tab, Tba, Tbb
        complex(dp), dimension(:), allocatable :: b_vecx, b_vecy, xx, xy  
        complex(dp), dimension(:), pointer :: param
        complex(dp), dimension(:), allocatable :: a_out, b_out, a_out2, b_out2 
        complex(dp), dimension(:), allocatable :: ac_out, bc_out, ac_out2, bc_out2 
        complex(dp), dimension(:), allocatable :: a_in, b_in, a_in2, b_in2
        complex(dp), dimension(:), allocatable :: rotD
        integer :: sph, i1, t1,t2,rate
        integer :: Nspheres, Nmax, Nave
        integer :: my_id, N_procs, ierr, nstart, nstop, block_size
        integer :: las, n
        integer, dimension(:,:), allocatable :: indD


        !!Allocate space for Cexts
        allocate(Cexts(2,n_samples))

        allocate(Tmat(0))
        Cexts(:,:)=0.0_dp
        call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
        call MPI_COMM_SIZE (MPI_COMM_WORLD, N_procs, ierr)
        
        maxrad = volrad
        vol = maxrad**3*4.0_dp/3.0_dp*pi
        
        ka =  dble(k)*maxrad
        Nmax = truncation_order(ka) 
        
        if(my_id == 0) then
            write(6,*) "#1: Compute coherent T-matrix"
            write(6,*) 'ka =', ka
            write(6,*) 'T-matrix order', Nmax
        end if
        
        allocate(Taa( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        allocate(Tab( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        allocate(Tba( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        allocate(Tbb( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        
        allocate(Taa2( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        allocate(Tab2( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        allocate(Tba2( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        allocate(Tbb2( (Nmax+1)**2-1,(Nmax+1)**2-1 ))
        
        allocate(Taa_data((Nmax+1)**2-1,(Nmax+1)**2-1,n_samples))
        allocate(Tab_data((Nmax+1)**2-1,(Nmax+1)**2-1,n_samples))
        allocate(Tba_data((Nmax+1)**2-1,(Nmax+1)**2-1,n_samples))
        allocate(Tbb_data((Nmax+1)**2-1,(Nmax+1)**2-1,n_samples))
        
        
        Taa2(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        Tab2(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        Tba2(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        Tbb2(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        Taa_data(:,:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        Tab_data(:,:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        Tba_data(:,:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        Tbb_data(:,:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
        block_size = int(ceiling(dble(Nave)/dble(N_procs)))
        nstart = 1 + my_id * block_size
        nstop =  (my_id + 1) * block_size
        
        if((my_id)*block_size > Nave) then
            nstart = Nave
        end if
        
        if((my_id+1)*block_size > Nave) then 
            nstop = Nave
        end if
        
        
        
        do i1 = nstart, nstop 
            call system_clock(t1,rate)
            !*** Build geometry ****************
            
            call pack_spheres(maxrad,dens, epsr,i1, Nspheres, coord,radius,param)
            allocate(sphere(Nspheres))
            do sph = 1, Nspheres   
                sphere(sph)%cp = coord(:,sph)
                sphere(sph)%r = radius(sph)
                sphere(sph)%eps_r = param(sph)
                
                ka = real(k) * radius(sph)
                sphere(sph)%Nmax = truncation_order(ka)
            end do
            !**************************
            
            call compute_T_matrix(sphere, k, Nmax, Taa, Tab, Tba, Tbb)
            !!get total extinctions
            call tr_T(Taa, Tbb, Tab, Tba, real(k), crs1)
            
            Taa2 = Taa2 + Taa/dble(Nave)
            Tab2 = Tab2 + Tab/dble(Nave)
            Tba2 = Tba2 + Tba/dble(Nave)
            Tbb2 = Tbb2 + Tbb/dble(Nave)
            
                if(i1<=n_samples) then
                    cexts(2,i1) = crs1(3) 
                    Taa_data(:,:,i1) = Taa
                    Tab_data(:,:,i1) = Tab
                    Tba_data(:,:,i1) = Tba
                    Tbb_data(:,:,i1) = Tbb
                endif
            
            
            
            deallocate(sphere)
            deallocate(coord,radius,param)
            call system_clock(t2)
            
            write(6,'(I4,2X,A,F10.5,2X,A,2X,I5,A,I5,2X,A,2X,I5,2X,A,2X,F5.1,A)') my_id,'time:',  &
                                        &   real(t2-t1)/real(rate),'T-matrix', i1, '/',      &
                                        &   Nave, "spheres" ,Nspheres, "work done",          &
                                        &   (i1-nstart)*1.0/(nstop-nstart)*100, "%"
        end do
        
        
        if(my_id == 0) then
            call MPI_REDUCE(MPI_IN_PLACE,Taa2,size(Taa2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Tab2,size(Tab2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Tba2,size(Tba2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Tbb2,size(Tbb2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            
        else
            call MPI_REDUCE(Taa2,Taa2,size(Taa2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Tab2,Tab2,size(Tab2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Tba2,Tba2,size(Tba2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Tbb2,Tbb2,size(Tbb2),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            
        end if
        
        call MPI_BCAST(Taa2, size(Taa2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Tab2, size(Tab2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Tba2, size(Tba2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Tbb2, size(Tbb2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD,ierr)
        
        !*************************************************************************
        !*
        !*
        !*********************************************************************
        
        allocate(a_out((Nmax+1)**2-1), b_out((Nmax+1)**2-1))
        allocate(a_out2((Nmax+1)**2-1), b_out2((Nmax+1)**2-1))
        
        allocate(ac_out((Nmax+1)**2-1), bc_out((Nmax+1)**2-1))
        allocate(ac_out2((Nmax+1)**2-1), bc_out2((Nmax+1)**2-1))
        
        allocate(a_in((Nmax+1)**2-1), b_in((Nmax+1)**2-1))
        allocate(a_in2((Nmax+1)**2-1), b_in2((Nmax+1)**2-1))
        
        allocate(S_ave(N_theta*N_phi,18))
        
        
        S_ave(:,:) = 0.0
        
        Csca_ave = 0.0
        Cabs_ave = 0.0
        Csca_ave_tot = 0.0
        
        call planewave(Nmax, dble(k), a_in, b_in)
        
        las = 0
        do n = 1,Nmax
            las = las + (2*n+1)**2
        end do
        
        allocate(rotD(las))
        allocate(indD(las,2))
        
        call sph_rotation_sparse(0.0_dp, -pi/2.0_dp, Nmax, rotD, indD)
        a_in2 = sparse_matmul(rotD,indD,a_in,(Nmax+1)**2-1)
        b_in2 = sparse_matmul(rotD,indD,b_in,(Nmax+1)**2-1)
        
        ! Coherent field coefficients
        ac_out = matmul(Taa2,a_in) + matmul(Tab2,b_in)
        bc_out = matmul(Tbb2,b_in) + matmul(Tba2,a_in)
        
        ac_out2 =  matmul(Taa2,a_in2) + matmul(Tab2,b_in2)
        bc_out2 =  matmul(Tbb2,b_in2) + matmul(Tba2,a_in2)
        
        
        if(my_id == 0 ) write(6,*) "#2: Compute incoherent T-matrices and Mueller matrix"
        !******* Compute incoherent mueller matrix
        do i1 = nstart, nstop 
            
            call system_clock(t1,rate)
            
            call pack_spheres(maxrad,dens,epsr,i1,Nspheres, coord,radius,param)
            allocate(sphere(Nspheres))
            
            do sph = 1, Nspheres   
                sphere(sph)%cp = coord(:,sph)
                sphere(sph)%r = radius(sph)
                sphere(sph)%eps_r = param(sph)
                
                ka = real(k) * radius(sph)
                sphere(sph)%Nmax = truncation_order(ka)
            end do
            !**************************
            
            call rhs2_xy(sphere, Tmat, b_vecx, b_vecy, k)
            allocate(xx(size(b_vecx)), xy(size(b_vecy)))
            
            
            call gmres(sphere, k, b_vecx, xx, dble(1e-3), 10, 50)
            call gmres(sphere, k, b_vecy, xy, dble(1e-3), 10, 50)
            
            call compute_cluster_absorption(sphere, xx, k, Cabs)
            call compute_cluster_absorption(sphere, xy, k, Cabs2)
            
            Cabs_ave = Cabs_ave + (Cabs+Cabs2)/2.0/Nave
            
            call cluster_coeff(sphere, xx, k, a_out, b_out, Nmax)
            call cluster_coeff(sphere, xy, k, a_out2, b_out2, Nmax)
            
            call scattering_cross_section(a_out, b_out, k, Nmax, Csca)
            call scattering_cross_section(a_out2, b_out2, k, Nmax, Csca2)
            
            Csca_ave_tot = Csca_ave_tot + (Csca + Csca2)/2.0/Nave
            
            
            !Incoherent field coeff.
            a_out = a_out - ac_out 
            b_out = b_out - bc_out
            
            a_out2 =  a_out2 - ac_out2 
            b_out2 =  b_out2 - bc_out2 
            
            
            call scattering_cross_section(a_out, b_out, k, Nmax, Csca)
            call scattering_cross_section(a_out2, b_out2, k, Nmax, Csca2)
            
            Csca_ave = Csca_ave + (Csca + Csca2)/2.0/Nave
            
            call mueller_matrix_coeff(a_out, b_out, a_out2, b_out2, k, N_theta, N_phi, S_out)
            
            S_ave = S_ave + S_out/dble(Nave)
            
            
            if(i1<=n_samples) then
                Taa_data(:,:,i1) = Taa_data(:,:,i1) - Taa2
                Tab_data(:,:,i1) = Tab_data(:,:,i1) - Tab2
                Tba_data(:,:,i1) = Tba_data(:,:,i1) - Tba2
                Tbb_data(:,:,i1) = Tbb_data(:,:,i1) - Tbb2
                call tr_T(Taa_data(:,:,i1),Tbb_data(:,:,i1),Tab_data(:,:,i1),Tba_data(:,:,i1), real(k), crs2)
                cexts(1,i1) = (cexts(2,i1)+crs2(2))
                cexts(2,i1) = crs2(2)/(cexts(2,i1)+crs2(2))
            endif
            
            
            
            deallocate(S_out)
            deallocate(xx,xy, b_vecx, b_vecy)
            deallocate(sphere)
            deallocate(coord,radius,param)
            
            call system_clock(t2)
            
            write(6,'(I4,2X,A,F10.5,2X,A,2X,I5,A,I5,2X,A,2X,I5,2X,A,2X,F5.1,A)') my_id,'time:',  &
                                        &   real(t2-t1)/real(rate),'T-matrix', i1, '/',      &
                                        &   Nave, "spheres" ,Nspheres, "work done",          &
                                        &   (i1-nstart)*1.0/(nstop-nstart)*100, "%"
        end do
        
        
        
        if(my_id == 0) then
            write(6,*) "Gathering information from other processes"
            call MPI_REDUCE(MPI_IN_PLACE,Cexts,size(Cexts),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            
            
            call MPI_REDUCE(MPI_IN_PLACE,S_ave,size(S_ave),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Csca_ave,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Cabs_ave,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Csca_ave_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            
            call MPI_REDUCE(MPI_IN_PLACE,Taa_data,size(Taa_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Tab_data,size(Tab_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Tba_data,size(Tba_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(MPI_IN_PLACE,Tbb_data,size(Tbb_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            write(6,*)  "Ready"
        else
            
            call MPI_REDUCE(Cexts,Cexts,size(Cexts),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            
            call MPI_REDUCE(S_ave,S_ave,size(S_ave),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Csca_ave,Csca_ave,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Cabs_ave,Cabs_ave,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Csca_ave_tot,Csca_ave_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            
            call MPI_REDUCE(Taa_data,Taa_data,size(Taa_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Tab_data,Tab_data,size(Tab_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Tba_data,Tba_data,size(Tba_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(Tbb_data,Tbb_data,size(Tbb_data),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)  
        end if
        
        kappa_s = Csca_ave / vol
        kappa_a = Cabs_ave / vol
        kappa_s_tot = Csca_ave_tot / vol
    end subroutine compute_T_matrix_ave
    
end module T_matrix
    
