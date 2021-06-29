
!-------------------------------------------------------------------
! Radius generator of the packing algorithm (single scatterers). 
! This version generates power law distributed radiuses (Cx^-n)
! Arguments are checked with parseArgument(), if the given
! argument is not recognized in the main.f90.
!
! Arguments given to a command line:
! -n            :: power law index (must be positive)  (Default: 2.0)
! -lowB         :: lower bound    (Default: 1.0)
! -upB          :: upper bound    (Default: 2.0)
!
! Copyright (C) 2018 Timo Väisänen,  University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 19.3.2018 TV: Maintenance, cleaning, commenting
!
! List of functions/subroutine:
! *parse_arguments(arg,i_arg) result(retVal)
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
real(kind=dp) :: lowlim = 1.0
real(kind=dp) :: maxlim = 2.0
real(kind=dp) :: n = 2.0
private lowlim,maxlim,n

contains

!!Parse command line arguments
!!Returns .true. if recognised otherwise .false. 
function parse_arguments(arg,i_arg) result(retVal)
    character(32), intent(in) :: arg
    integer, intent(in) :: i_arg
    logical :: retVal
    character(32) :: tmp
    retVal = .true.
    select case(arg)
        case('-lowB')
            call get_command_argument(i_arg+1,tmp)
            read(tmp,*) lowlim
        case('-upB')
            call get_command_argument(i_arg+1,tmp)
            read(tmp,*) maxlim
        case('-npower')
            call get_command_argument(i_arg+1,tmp)
            read(tmp,*) n
            if(n<=0) then
                write(6,*) "Given power law index must be positive (Cx^-n)"
                stop
            endif
        case default
            retVal = .false.
    end select
end function

!!Call this to generate radiuses
function rad_gen_get_radius() result(retVal)
    real(kind=dp) :: retVal
    retVal = rnd_power_law(n,lowlim,maxlim)
end function 

!!Share arguments with other MPI processes
subroutine distribute_args_MPI()
    integer :: ierr
    call MPI_Bcast(lowlim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(maxlim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine

!!Print initial information
subroutine print_RG_information()
    write(6,*) "====================================="
    write(6,*) "Radius generator: Power law distribution"
    write(6,*) "lowB=",lowlim,",upB",maxlim
    write(6,*) "npower=",n
    write(6,*) "====================================="
end subroutine

!!Print scatterer information to stream
subroutine print_scatterer_information(stream)
    integer, intent(in) :: stream
    write(stream,*) "Power law distribution"
    write(stream,*) "lowB",lowlim
    write(stream,*) "upB",maxlim
    write(stream,*) "npower",n
end subroutine

end module