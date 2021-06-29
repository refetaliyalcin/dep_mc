#define SIZEOFDATABLOCK 39
#define NUMOFARRAYSINDB 15
#define SIZEOFBUFFER 20000
#define SIZEOFBUFFER2 1000


module mpi_tools
!!MPI tools which are used when MPI libraries are not present. 
!!
!!This needs update!
!!
!! Copyright (C) 2016 Timo Väisänen and University of Helsinki
!! All rights reserved.
!! The new BSD License is applied to this software, see LICENSE.txt

    use typeDefinitions
    use constants
    implicit none
    character :: buffer(SIZEOFBUFFER)
    real(kind=rk) :: buffer2(SIZEOFBUFFER2)
    integer :: listOfInds(NUMOFARRAYSINDB)
contains

    subroutine initializeMPI(rank,ierr,MPIprocs)
        integer, intent(out) :: ierr,rank,MPIprocs
        rank=0
        MPIprocs=1
        ierr=0
    end subroutine

    subroutine sendDataDB(dB)
        type(dataBlock), intent(in) :: dB
        integer :: lengths(SIZEOFDATABLOCK)
        integer :: pos,ierr,nodes,j1,w4root
    end subroutine



    subroutine receiveDataDB(dB)
        type(dataBlock), intent(inout) :: dB
        integer :: pos,ierr,w4Root
    end subroutine

    
#define sizeOfBufferMerge 50000
    subroutine mergeAssembledData(aD)
        type(assembledData), intent(inout) :: aD
        real(kind=rk) :: buffer(sizeOfBufferMerge)
        integer :: listOfInds(9),ierr
        integer :: rank,j1,j2,j3,j4,j5,x,nodes,finished
    end subroutine

    subroutine finishMPISession()
    end subroutine

    subroutine ABORT(msg)
        integer :: ierr,errorcode
        character(*), intent(in) :: msg
        write(6,*) trim(msg)
        stop
    end subroutine


    subroutine bcastThreadCount(threadCount)
        integer, intent(inout) :: threadCount
    end subroutine

end module



