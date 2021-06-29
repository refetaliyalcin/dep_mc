!-------------------------------------------------------------------
!The program computes incoherent scattering coefficients with the 
!Fast superposition T-matrix Method 
!
! Main program
!
! Copyright (C) 2018 Johannes Markkanen, University of Helsinki
! All rights reserved.
!
!
! @Note: TV: Some changes were made in 28.12.2016 
! @Note: AP: Changes to output files were made 27.1.2017
! @Note: TV: fixed: x1 was not initialized before used
! @Note: TV: 27.03.2018, Cleaning
!
! TODO!
! Refactoring
!-----------------------------

program main 
use Common
use IO
use T_matrix
use MPI
use RadiusGenerator
use RNG
implicit none

character(len=64) :: arg_name, arg, T_out, T_out2, lmk_file, S_out_prefix, kappa_out_prefix, &
     S_out, kappa_out

character(8)  :: date0
character(10) :: time0
character(5)  :: zone

complex(dp) :: k, epsr
complex(dp), dimension(:,:), allocatable :: Taa, Tab, Tba, Tbb
complex(dp), dimension(:,:,:), allocatable :: Taa_data, Tab_data, Tba_data, Tbb_data

real(dp), dimension(:,:), allocatable :: S_ave, S_ave2, data2
real(dp) :: k_in, packdens, eps_r, eps_i, lambda
real(dp) :: kappa_s, kappa_a, kappa_s_tot, rad, vol_rad

real(kind=dp),allocatable :: Cexts(:,:) 
integer :: Nspheres,  Nave, poi
integer :: ierr, my_id, rc,  read_lmk
integer :: i_arg, num_args, th, phi, N_phi, N_theta, ii, jj
integer,dimension(8) :: values
integer :: k1, col, row, T1, T2, rate
integer :: errorcode,n_samples,j2, world_size

character(len=8) :: fmt
character(5) :: x1


fmt = '(I5.5)'


call MPI_INIT(ierr)

if (ierr/=MPI_SUCCESS) then
   write(6,*) "Error starting MPI program. Terminating."
   call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   stop
end if

call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call init_random_seed()


if(my_id == 0) then
   ! Default args.
   n_samples = 128
   rad = 2.0_dp
   k_in = 1.0_dp
   packdens = 0.15d0
   eps_r = 1.7161d0
   eps_i = 0.01d0
   T_out = 'T.h5'
   T_out2 = 'T_multi.h5'
   S_out_prefix = 'mueller'
   Nave = 128
   N_theta = 180
   N_phi = 64
   kappa_out_prefix = 'kappa'
   read_lmk = 0 
   lmk_file = 'data.txt'
   vol_rad = 2*pi

    ! Read command line arguments 
    num_args = command_argument_count()
    do i_arg = 1,num_args,2
        call get_command_argument(i_arg,arg_name)
        select case(arg_name)
        case('-n_samples')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) n_samples
        case('-k')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) k_in 
        case('-vol_rad')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) vol_rad
        case('-N_ave')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) Nave 
        case('-read_lmk')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) read_lmk
        case('-read_file')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) lmk_file 
        case('-N_theta')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) N_theta
        case('-N_phi')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) N_phi
        case('-density')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) packdens 
        case('-eps_r')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) eps_r
        case('-eps_i')
            call get_command_argument(i_arg+1,arg)
            read(arg,*) eps_i 
        case('-T_out')
            call get_command_argument(i_arg+1,arg)
            T_out = arg
        case('-T_out2')
            call get_command_argument(i_arg+1,arg)
            T_out2 = arg
        case('-S_out_prefix')
            call get_command_argument(i_arg+1,arg)
            S_out_prefix = arg
        case('-kappa_out_prefix')
            call get_command_argument(i_arg+1,arg)
            kappa_out_prefix = arg
        case default 
            if(parse_arguments(arg_name,i_arg)) cycle
            print '(a,a,/)', 'Unrecognized command-line option: ', arg_name
            call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)  
        end select
    end do

    call MPI_Comm_size(MPI_COMM_WORLD, world_size,ierr)
    if(Nave< world_size) Nave =  world_size
    if(n_samples>Nave) n_samples=Nave

    k = dcmplx(k_in,0.0)

    if(read_lmk /= 0) then
        allocate(data2(4,read_lmk))
        open(unit=98,file = lmk_file ,status='old')
        do row = 1, read_lmk
            read(98,*) (data2(col,row),col=1, 4)
        end do
        close(98)
    end if

    if(read_lmk == 0) then
        read_lmk = 1
        allocate(data2(4,1))
        data2(1,1) = real(k)
        data2(2,1) = imag(k)
        data2(3,1) = eps_r
        data2(4,1) = eps_i
    end if

end if


call MPI_Bcast(n_samples,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(read_lmk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(poi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(N_phi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(N_theta,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(k,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Nave,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(packdens,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(eps_r,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(eps_i,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(vol_rad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!!Distribute radius generator information and print information
call distribute_args_MPI()
if(my_id == 0) call print_RG_information()  

if(my_id /= 0 ) allocate(data2(4,read_lmk))

call MPI_Bcast(data2,size(data2),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

if(my_id == 0) then
    write(6,*) 'n_samples', n_samples
    write(6,*) 'Volume element radius', vol_rad
end if


!**************************************************************************
!
!    T-matrix
!***********************************************
do k1 = 1, read_lmk
    if(my_id == 0 .and. read_lmk>1) then
        write(6,*) 'solving:', k1, "out of",read_lmk,"jobs"
    end if
    call system_clock(T1,rate)

    lambda = 2*pi/data2(1,k1)

    k =  cmplx(data2(1,k1),data2(2,k1),kind=dp)
    epsr = cmplx(data2(3,k1),data2(4,k1),kind=dp)

    if(my_id == 0) then
        write(6,*) 'wavelength', lambda
        write(6,*) 'eps_r =',  epsr
        write(6,*) 'k', k    
    end if

    call compute_T_matrix_ave(vol_rad, packdens, epsr, k, Nave, &
        &   Taa, Tab, Tba, Tbb, N_theta, N_phi, S_ave, kappa_s, kappa_a, &
        &   kappa_s_tot, Taa_data, Tab_data, Tba_data, Tbb_data,n_samples,cexts)
   
    if(my_id == 0) then
      
        allocate(S_ave2(N_theta,17))
        S_ave2(:,:) = 0.0_dp
        do th = 1,N_theta
            do phi = 1,N_phi
                S_ave2(th,1:17) = S_ave2(th,1:17) + S_ave(th+N_theta*(phi-1),2:18)/N_phi;           
            end do
        end do
    
        S_ave2(:,1) = S_ave2(:,1)*180/pi;

        write (x1,fmt) k1
        call T_write2file(Taa, Tab, Tba, Tbb, T_out)


        !If eps_i == 0, set albedo to 1
        if(eps_i==0.0_dp) then
            Cexts(2,:) = 1.0_dp
        endif

        call T_write2file2(Taa_data, Tab_data, Tba_data, Tbb_data, Cexts, T_out2)
        S_out = trim(S_out_prefix)//'_'//trim(x1)//'.h5'
        call real_write2file(S_ave2, S_out)

        S_out = trim(S_out_prefix)//'_'//trim(x1)//'.dat'
        kappa_out = trim(kappa_out_prefix)//'_'//trim(x1)//'.dat'

        OPEN(UNIT=12, FILE=S_out, ACTION="write", STATUS="replace")
        DO ii=1,N_theta
            do j2=1,17
                WRITE(12,'(E15.6)',advance='no') S_ave2(ii,j2)
            enddo
            write(12,*) 
        END DO
        CLOSE(12)

        call date_and_time(date0,time0,zone,values)  
        call print_output(6)            !!Print the output to console
        call print_output(123)       !!Print the output to a file

        deallocate(S_ave2)

    end if

    deallocate(Taa, Tbb, Tab, Tba)
    deallocate(Taa_data, Tbb_data, Tab_data, Tba_data)
    deallocate(S_ave)
    call system_clock(T2)
    if(my_id == 0) then
        write(6,*) 'Done in', real(T2-T1)/real(rate), 'seconds'
    end if
end do

call MPI_FINALIZE(ierr)
stop        !!Print ieee flags

contains

subroutine print_output(stream)
    integer, intent(in) :: stream
    if(stream>6) OPEN(UNIT=stream, FILE=kappa_out, ACTION="write", STATUS="replace")
    WRITE(stream,*) "time stamp: ", date0, time0
    WRITE(stream,*) 'n_samples', n_samples
    WRITE(stream,*) 'n_ave', Nave
    WRITE(stream,*) 'lambda_e = ', lambda
    WRITE(stream,*) 'eps_r =',  epsr
    WRITE(stream,*) 'n =',  sqrt(epsr)
    WRITE(stream,*) 'vol_rad =', vol_rad 
    WRITE(stream,*) 'density =', packdens 
    WRITE(stream,*) 'k', k    
    WRITE(stream,*) 'kappa_s_ic/k =', kappa_s/dble(k)
    WRITE(stream,*) 'kappa_a/k =', kappa_a/dble(k)
    WRITE(stream,*) 'kappa_s/k =', kappa_s_tot/dble(k)
    WRITE(stream,*) 'albedo =', kappa_s / (kappa_s + kappa_a)
    WRITE(stream,*) "incoherent mean free path", 1.0_dp/(kappa_s+kappa_a)
    WRITE(stream,*) "free space mean free path", 1.0_dp/(kappa_s_tot+kappa_a)
    call print_scatterer_information(stream)
    if(stream>6) close(stream)
end subroutine

end program
