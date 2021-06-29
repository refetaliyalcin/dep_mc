!-------------------------------------------------------------------
! HDF5 interfacing
!
! Copyright (C) 2018 Johannes Markkanen, Timo Väisänen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 19.3.2018 TV: Complete overhaul
!
! List of functions:
! *write_array2file_real2D(M,dims,file_id,dsetname,dspace_id) result(error)
! *write_array2file_real3D(M,dims,file_id,dsetname,dspace_id) result(error)
! *read_1real2D_array(M,file_id,dsetname,dspace_id) result(error)
!
! List of subroutines:
! *subroutine T_write2file(Taa, Tab, Tba, Tbb, fname)
! *subroutine read_T(fname, Taa, Tab, Tba, Tbb)
! *subroutine T_write2file2(Taa, Tab, Tba, Tbb, Arr2D, fname)
! *subroutine read_T2(fname, Taa, Tab, Tba, Tbb, Arr2D)
! *subroutine real_write2file(A, fname)
! *subroutine readTxx2D(Txx,Mr,dsetname1,dsetname2,file_id,dspace_id)
! *subroutine readTxx3D(Txx,Mr,dsetname1,dsetname2,file_id,dspace_id)
!
!
! TODO
! *Remove copy-pasta
! *Add compression
! *Better error messaging
!-----------------------------
module io
  use iso_fortran_env, only : real64
  use hdf5
  implicit none
  integer, parameter :: dp = real64
  integer, parameter :: ESTREAM = 6
  character(len=64), parameter :: ERRSTR= "ERROR: io module: " 
  private :: write_array2file_real2D,write_array2file_real3D,read_1real2D_array
  private :: readTxx2D,readTxx3D,dp
contains 

!****************************************************
!!Writes single T-matrix (Taa,Tab,Tba,Tbb) to the file (fname)
subroutine T_write2file(Taa, Tab, Tba, Tbb, fname)
  complex(dp),      intent(in) :: Taa(:,:)
  complex(dp),      intent(in) :: Tab(:,:)
  complex(dp),      intent(in) :: Tba(:,:)
  complex(dp),      intent(in) :: Tbb(:,:)
  character(LEN=38), intent(in):: fname

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dspace_id     ! Dataspace identifier
  integer(HSIZE_T)  :: dims(2)  ! Dataset dimensions
  integer     ::    rank = 2                       ! Dataset rank
  integer     ::   error ! Error flag

  dims = (/size(Taa,1),size(Taa,2)/) 
  call h5open_f(error)
  ! Create a new file using default properties.
  call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
  ! Create the dataspace.
  call h5screate_simple_f(rank, dims, dspace_id, error)
  !******************************************************************
  !!imag is not part of the fortran standard but bewary that aimag might
  !some issues with precisions when using different compiler
  if(write_array2file_real2D(real(Taa), dims,file_id,"Taa_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,1
  if(write_array2file_real2D(aimag(Taa),dims,file_id,"Taa_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,2
  if(write_array2file_real2D(real(Tab), dims,file_id,"Tab_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,3
  if(write_array2file_real2D(aimag(Tab),dims,file_id,"Tab_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,4
  if(write_array2file_real2D(real(Tba), dims,file_id,"Tba_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,5
  if(write_array2file_real2D(aimag(Tba),dims,file_id,"Tba_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,6
  if(write_array2file_real2D(real(Tbb), dims,file_id,"Tbb_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,7
  if(write_array2file_real2D(aimag(Tbb),dims,file_id,"Tbb_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,8
  !******************************************************************
  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, error)
  ! Close the file.
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)
end subroutine

!****************************************************
!!Read T-matrix (Taa,Tab,Tba,Tbb) from the file (fname)
subroutine read_T(fname, Taa, Tab, Tba, Tbb)
  complex(dp), allocatable, intent(out) :: Taa(:,:),Tab(:,:)
  complex(dp), allocatable, intent(out) :: Tba(:,:),Tbb(:,:)
  character(len=38), intent(in)       :: fname  

  integer(HID_T) :: file_id 
  real(dp), allocatable :: Mr(:,:), Mi(:,:)
  integer :: error

  call h5open_f(error) ! initialize interface
  call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error) ! open file

  !****************************************
  call readTxx2D(Taa,Mr,"Taa_r","Taa_i",file_id)
  call readTxx2D(Tab,Mr,"Tab_r","Tab_i",file_id)
  call readTxx2D(Tba,Mr,"Tba_r","Tba_i",file_id)
  call readTxx2D(Tbb,Mr,"Tbb_r","Tbb_i",file_id)
  !****************************************

  call h5fclose_f(file_id, error) ! close file
  call h5close_f(error)           ! close inteface

  if(allocated(Mr)) deallocate(Mr)
  if(allocated(Mi)) deallocate(Mi)
end subroutine



!****************************************************
!!New R2T2 compatible
!!Writes multiple T-matrices (Taa,Tab,Tba,Tbb) and one 2D array
!!to the file (fname)
subroutine T_write2file2(Taa, Tab, Tba, Tbb, Arr2D, fname)
  complex(dp),      intent(in) :: Taa(:,:,:),Tab(:,:,:)
  complex(dp),      intent(in) :: Tba(:,:,:),Tbb(:,:,:)
  real(dp),         intent(in) :: Arr2D(:,:)
  character(len=38), intent(in):: fname

  integer(HID_T)   :: file_id    ! File identifier
  integer(HID_T)   :: dspace_id  ! Dataspace identifier
  integer(HSIZE_T) :: dims(3)    ! Dataset dimensions
  integer(HSIZE_T) :: dims2(2)   ! Dataset dimensions
  integer          :: rank = 3   ! Dataset rank
  integer          :: rank2 = 2  ! Dataset rank
  integer          :: error      ! Error flag

  dims2 = (/size(Arr2D,1),size(Arr2D,2)/) 
  dims = (/size(Taa,1),size(Taa,2),size(Taa,3)/) 
  call h5open_f(error)
  ! Create a new file using default properties.
  call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
  ! Create the dataspace.
  call h5screate_simple_f(rank, dims, dspace_id, error)
  !******************************************************************
  !!imag is not part of the fortran standard but bewary that aimag might
  !some issues with precisions when using different compiler
  if(write_array2file_real3D(real(Taa), dims,file_id,"Taa_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,1
  if(write_array2file_real3D(aimag(Taa),dims,file_id,"Taa_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,2
  if(write_array2file_real3D(real(Tab), dims,file_id,"Tab_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,3
  if(write_array2file_real3D(aimag(Tab),dims,file_id,"Tab_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,4
  if(write_array2file_real3D(real(Tba), dims,file_id,"Tba_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,5
  if(write_array2file_real3D(aimag(Tba),dims,file_id,"Tba_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,6
  if(write_array2file_real3D(real(Tbb), dims,file_id,"Tbb_r",dspace_id)/=0) write(ESTREAM,*) ERRSTR,7
  if(write_array2file_real3D(aimag(Tbb),dims,file_id,"Tbb_i",dspace_id)/=0) write(ESTREAM,*) ERRSTR,8
  call h5sclose_f(dspace_id, error)



  call h5screate_simple_f(rank2, dims2, dspace_id, error)
  if(write_array2file_real2D(Arr2D,dims2,file_id,"Cexts",dspace_id)/=0) write(ESTREAM,*) ERRSTR,9
  call h5sclose_f(dspace_id, error)
  !******************************************************************

  ! Close the file.
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)
end subroutine


!****************************************************
!!New R2T2 compatible
!!Reads multiple T-matrices (Taa,Tab,Tba,Tbb) and one 2D array
!!from the file (fname)
subroutine read_T2(fname, Taa, Tab, Tba, Tbb, Arr2D)
  complex(dp), allocatable, intent(out) :: Taa(:,:,:), Tab(:,:,:)
  complex(dp), allocatable, intent(out) :: Tba(:,:,:), Tbb(:,:,:)
  real(dp), allocatable,    intent(out) :: Arr2D(:,:)
  character(len=38),        intent(in)  :: fname 

  integer(HID_T) :: file_id 
  real(dp), allocatable :: Mr(:,:,:)
  integer :: error

  call h5open_f(error) ! initialize interface
  call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error) ! open file

  !****************************************
  call readTxx3D(Taa,Mr,"Taa_r","Taa_i",file_id)
  call readTxx3D(Tba,Mr,"Tba_r","Tba_i",file_id)
  call readTxx3D(Tab,Mr,"Tab_r","Tab_i",file_id)
  call readTxx3D(Tbb,Mr,"Tbb_r","Tbb_i",file_id)
  if(read_1real2D_array(Arr2D,file_id,"Cexts")/=0) write(ESTREAM,*) ERRSTR,8
  !****************************************

  call h5fclose_f(file_id, error) ! close file
  call h5close_f(error)           ! close inteface

  if(allocated(Mr)) deallocate(Mr)
end subroutine


!****************************************************
!!Writes single 2D real array to the file (fname)
subroutine real_write2file(A, fname)
  real(dp),      intent(in) :: A(:,:)
  character(LEN=38), intent(in):: fname

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dspace_id     ! Dataspace identifier
  integer(HSIZE_T)  :: dims(2)  ! Dataset dimensions
  integer     ::    rank = 2                       ! Dataset rank
  integer     ::   error ! Error flag

  dims = (/size(A,1),size(A,2)/) 
  call h5open_f(error)
  ! Create a new file using default properties.
  call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
  ! Create the dataspace.
  call h5screate_simple_f(rank, dims, dspace_id, error)
  !******************************************************************
  !!imag is not part of the fortran standard but bewary that aimag might
  !some issues with precisions when using different compiler
  if(write_array2file_real2D(A, dims,file_id,"mueller",dspace_id)/=0) write(6,*) "asd"
  !******************************************************************
  ! Terminate access to the data space.
  call h5sclose_f(dspace_id, error)
  ! Close the file.
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)
end subroutine



!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

!****************************************************
!Writes single 2D array to a file (described by the file_id)
function write_array2file_real2D(M,dims,file_id,dsetname,dspace_id) result(error)
  real(dp),          intent(in) :: M(:,:)       !!dp 2D array
  integer(HSIZE_T),  intent(in) :: dims(:)      !!Dimensions of the array M
  integer(HID_T),    intent(in) :: file_id      !!File id
  integer(HID_T),    intent(in) :: dspace_id    !!dataspace identifier
  character(len=5),  intent(in) :: dsetname     !!The name of the dataset

  integer ::        error                       !!Returns error code
  integer(HID_T) :: dset_id

  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
        & dset_id, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, M(:,:), dims, error)
  call h5dclose_f(dset_id, error)
end function

!****************************************************
!Writes single 3D array to a file (described by the file_id)
function write_array2file_real3D(M,dims,file_id,dsetname,dspace_id) result(error)
  real(dp),          intent(in) :: M(:,:,:)     !!dp 3D array
  integer(HSIZE_T),  intent(in) :: dims(:)      !!Dimensions of the array M
  integer(HID_T),    intent(in) :: file_id      !!File id
  integer(HID_T),    intent(in) :: dspace_id    !!dataspace identifier
  character(len=5),  intent(in) :: dsetname     !!The name of the dataset

  integer ::        error                       !!Returns error code
  integer(HID_T) :: dset_id

  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
        & dset_id, error, H5P_DEFAULT_F, H5P_DEFAULT_F, H5P_DEFAULT_F)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, M(:,:,:), dims, error)
  call h5dclose_f(dset_id, error)
end function

!****************************************************
!Reads single 2D array from a file (described by the file_id)
function read_1real2D_array(M,file_id,dsetname) result(error)
  real(dp),allocatable, intent(inout) :: M(:,:)     !!dp 2D array
  integer(HID_T),       intent(in)    :: file_id    !!File id
  character(len=5),     intent(in)    :: dsetname   !!The name of the dataset

  integer ::            error                       !!Returns error code
  integer(HID_T)   ::   dset_id
  integer(HSIZE_T) ::   dims_out(2), dims(2)
  integer(HID_T)   ::   dspace_id  !!dataspace identifier
  logical :: reallocate

  reallocate = .false.
  call h5dopen_f(file_id, dsetname, dset_id, error)
  call h5dget_space_f(dset_id, dspace_id, error) 
  call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)

  !Do not reallocate the array if already sized correctly
  if(.not. allocated(M)) then
    reallocate =  .true.
  else
    if(size(M,dim=1)/=dims(1) .or. size(M,dim=2)/=dims(2)) then
      deallocate(M)
      reallocate =  .true.
    endif
  endif
  if(reallocate) allocate(M(dims(1), dims(2)))
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,M,dims,error)
  call h5dclose_f(dset_id,error) 
end function


!****************************************************
!Reads single 3D array from a file (described by the file_id)
function read_1real3D_array(M,file_id,dsetname) result(error)
  real(dp),allocatable, intent(inout) :: M(:,:,:)   !!dp 3D array
  integer(HID_T),       intent(in)    :: file_id    !!File id
  character(len=5),     intent(in)    :: dsetname   !!The name of the dataset

  integer ::            error,j1                    !!Returns error code
  integer(HID_T)   ::   dset_id
  integer(HSIZE_T) ::   dims_out(3), dims(3)
  integer(HID_T)   ::   dspace_id                   !!dataspace identifier
  logical :: reallocate

  reallocate = .false.
  call h5dopen_f(file_id, dsetname, dset_id, error)
  call h5dget_space_f(dset_id, dspace_id, error) 
  call H5sget_simple_extent_dims_f(dspace_id, dims_out, dims, error)
  !Do not reallocate the array if already sized correctly
  if(.not. allocated(M)) then
    reallocate =  .true.
  else
    do j1=1,3
      if(size(M,dim=j1)/=dims(j1)) then
        deallocate(M)
        reallocate =  .true.
        exit
      endif
    enddo
  endif
  if(reallocate) allocate(M(dims(1), dims(2), dims(3)))
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,M,dims,error)
  call h5dclose_f(dset_id,error) 
end function


!****************************************************
subroutine readTxx2D(Txx,Mr,dsetname1,dsetname2,file_id)
  complex(dp), allocatable, intent(inout) :: Txx(:,:)
  real(dp),allocatable, intent(inout)     :: Mr(:,:)
  integer(HID_T), intent(in)              :: file_id 
  character(len=5), intent(in)            :: dsetname1,dsetname2

  !if(read_1real2D_array(Mr,file_id,dsetname1)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  !if(read_1real2D_array(Mi,file_id,dsetname2)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  !allocate(Txx(size(Mr,dim=1),size(Mr,dim=2)))
  !Txx = cmplx(Mr, Mi,kind=dp)  

  if(read_1real2D_array(Mr,file_id,dsetname1)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  allocate(Txx(size(Mr,dim=1),size(Mr,dim=2)))
  Txx = cmplx(Mr, 0.0_dp,kind=dp)  
  if(read_1real2D_array(Mr,file_id,dsetname2)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  Txx = cmplx(real(Txx), Mr,kind=dp)  
end subroutine


!****************************************************
subroutine readTxx3D(Txx,Mr,dsetname1,dsetname2,file_id)
  complex(dp), allocatable, intent(inout) :: Txx(:,:,:)
  real(dp),allocatable, intent(inout)     :: Mr(:,:,:)
  integer(HID_T), intent(in)              :: file_id 
  character(len=5), intent(in)            :: dsetname1,dsetname2

  !if(read_1real3D_array(Mr,file_id,dsetname1)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  !if(read_1real3D_array(Mi,file_id,dsetname2)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  !allocate(Txx(size(Mr,dim=1),size(Mr,dim=2),size(Mr,dim=2)))
  !Txx = cmplx(Mr, Mi,kind=dp) 
  if(read_1real3D_array(Mr,file_id,dsetname1)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  allocate(Txx(size(Mr,dim=1),size(Mr,dim=2),size(Mr,dim=3)))
  Txx = cmplx(Mr, 0.0_dp,kind=dp)  
  if(read_1real3D_array(Mr,file_id,dsetname2)/=0) write(ESTREAM,*) ERRSTR,dsetname1
  Txx = cmplx(real(Txx), Mr,kind=dp)   
end subroutine


end module 

