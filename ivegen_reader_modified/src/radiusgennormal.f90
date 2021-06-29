
!-------------------------------------------------------------------
! Radius generator of the packing algorithm (single scatterers). 
! This version generates normal distributed radiuses (stdev>0) or
! gives the mean radius (stdev<=0). 
! Arguments are checked with parseArgument(), if the given
! argument is not recognized in the main.f90.
!
! Arguments given to a command line:
! -radius       :: mean radius          (Default: 2)
! -stdev        :: standard deviation   (Default: 0)
!
! Copyright (C) 2018 Timo Väisänen,  University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 19.3.2018 TV: Maintenance, cleaning, commenting
!
! List of functions/subroutine:
! *parseArgument(arg,i_arg) result(retVal)
! *rad_gen_get_radius() result(retVal)
! *distribute_args_MPI()
! *print_RG_information()
!-----------------------------

module RadiusGenerator
use common
use rng
use mpi
implicit none

!Default values
real(kind=dp) :: mean = 2.0
real(kind=dp) :: stdev = 0.0

private :: mean,stdev
contains

!!Parse command line arguments
!!Returns .true. if recognised otherwise .false. 
function parse_arguments(arg,i_arg) result(retVal)
    character(32), intent(in)   :: arg
    integer, intent(in)         :: i_arg
    logical :: retVal
    character(32) :: tmp
    retVal = .true.
    select case(arg)
        case('-mean')
            call get_command_argument(i_arg+1,tmp)
            read(tmp,*) mean
        case('-stdev')
            call get_command_argument(i_arg+1,tmp)
            read(tmp,*) stdev
        case default
            retVal = .false.
    end select
end function

!!Call this to generate radiuses
function rad_gen_get_radius() result(retVal)
    real(kind=dp) :: retVal
    retVal = rnd_normal(mean,stdev) 
end function 

!!Share arguments with other MPI processes
subroutine distribute_args_MPI()
    integer :: ierr
    call MPI_Bcast(stdev,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(mean,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end subroutine

!!Print initial information
subroutine print_RG_information()
    write(6,*) "====================================="
    write(6,*) "Radius generator: Normal distribution"
    write(6,*) "mean=",mean,", stdev=",stdev
    write(6,*) "====================================="
end subroutine

!!Print scatterer information to stream
subroutine print_scatterer_information(stream)
    integer, intent(in) :: stream
    write(stream,*) "Normal distribution"
    write(stream,*) "Mean",mean
    write(stream,*) "stdev",stdev
end subroutine

end module
