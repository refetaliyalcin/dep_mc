#include "macros.inc"

#define MPIADDRESS(ARG) call MPI_Get_Address(ARG,disp(j1)); j1=j1+1

#define MPIUNPACKREAL(ARG) call MPI_Unpack(buffer,SIZEOFBUFFER,pos,ARG,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD)   
#define MPIUNPACKINTEGER(ARG) call MPI_Unpack(buffer,SIZEOFBUFFER,pos,ARG,1,MPI_INTEGER,MPI_COMM_WORLD)
#define MPIUNPACKLOGICAL(ARG) call MPI_Unpack(buffer,SIZEOFBUFFER,pos,ARG,1,MPI_LOGICAL,MPI_COMM_WORLD)
#define MPIUNPACKCMPLX(ARG) call MPI_Unpack(buffer,SIZEOFBUFFER,pos,ARG,1,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD)  

#define MPIUNPACKINTARR(ARG,SIZE) call MPI_Unpack(buffer,SIZEOFBUFFER,pos,ARG,SIZE,MPI_INTEGER,MPI_COMM_WORLD)
#define MPIUNPACKREALARR(ARG,SIZE) call MPI_Unpack(buffer,SIZEOFBUFFER,pos,ARG,SIZE,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD)
#define MPIUNPACKCHARARR(ARG,SIZE) call MPI_Unpack(buffer,SIZEOFBUFFER,pos,ARG,SIZE,MPI_CHARACTER,MPI_COMM_WORLD)  

#define ADDLENGTH(ARG)  length(ARG)=1; j1=j1+1
#define TYPELISTADDREAL typelist(j1)=MPI_DOUBLE_PRECISION;j1=j1+1
#define TYPELISTADDINTEGER typelist(j1)=MPI_INTEGER;j1=j1+1
#define TYPELISTADDCOMPLEX typelist(j1)=MPI_DOUBLE_COMPLEX;j1=j1+1
#define TYPELISTADDLOGICAL typelist(j1)=MPI_LOGICAL;j1=j1+1

#define SUBREADBUFFER(ARG) call subReadBuffer(ind,indPos,ARG)
#define COPYARRAY(ARG) call copyArray(ARG,listOfInds,initPos,indsPos)

#define SIZEOFDATABLOCK 200
#define NUMOFARRAYSINDB 15
#define SIZEOFBUFFER 50000
#define SIZEOFBUFFER2 4000

module mpi_tools
!! Copyright (C) 2016 Timo Väisänen and University of Helsinki
!! All rights reserved.
!! The new BSD License is applied to this software, see LICENSE.txt
    use typeDefinitions
    use constants
    use error_handler
    use mpi_f08
    implicit none
    character, allocatable :: buffer(:)
    real(kind=rk), allocatable :: buffer2(:)
    integer :: listOfInds(NUMOFARRAYSINDB)

    type(MPI_Datatype) :: newType
    private buffer,buffer2,fillTypeListAndLengths,fillBuffer2
    private newType,listOfInds,readBuffer,subReadBuffer
contains


    subroutine initializeMPI(rank,ierr,processes)
    !!The rank, error flag and the number of MPI processes
        integer, intent(out) :: ierr      !error flag
        integer, intent(out) :: rank      !rank(master must have always 0)
        integer, intent(out) :: processes !number of MPI processes
        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,processes,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    end subroutine



    subroutine send_data_DB(dB)
    !!Broadcast dataBlock to every MPI process
        type(dataBlock), intent(in) :: dB       !datablock
        integer :: lengths(SIZEOFDATABLOCK)
        type(MPI_Datatype) :: newType, typeList(SIZEOFDATABLOCK)
        integer(kind=MPI_ADDRESS_KIND) :: disp(SIZEOFDATABLOCK)  
        integer :: pos,ierr,nodes,j1,w4root
        w4root = 1
        allocate(buffer2(SIZEOFBUFFER2))
        allocate(buffer(SIZEOFBUFFER))
        call fillBuffer2(dB)
        call fillTypeListAndLengths(typeList,lengths)
        call fillDisplacements(dB,disp)
        call MPI_Type_create_struct(SIZEOFDATABLOCK,lengths,disp,typeList,newtype)
        call MPI_Type_commit(newType)
        pos = 0
        call MPI_Pack(MPI_BOTTOM, 1, newType, buffer, SIZEOFBUFFER, pos, MPI_COMM_WORLD,ierr)
        !call MPI_BCAST(w4Root,1,MPI_INTEGER,0,MPI_COMM_WORLD)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nodes, ierr)
        do j1=1,nodes-1
            call MPI_Send(buffer, pos, MPI_PACKED, j1, 123456, MPI_COMM_WORLD)
        enddo
        deallocate(buffer)
        deallocate(buffer2)
    end subroutine


    subroutine receive_data_DB(dB)
    !!Receive a package from the master process and unpack it.
        type(dataBlock), intent(inout) :: dB        !!datablock
        integer :: pos,ierr,w4Root
        pos = 0
        allocate(buffer2(SIZEOFBUFFER2))
        allocate(buffer(SIZEOFBUFFER))
        !call MPI_BCAST(w4Root,1,MPI_INTEGER,0,MPI_COMM_WORLD)
        call MPI_Recv(buffer,SIZEOFBUFFER,MPI_PACKED,0,123456,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        call openPackage(dB,pos)
        call readBuffer(dB)
        deallocate(buffer)
        deallocate(buffer2)
    end subroutine

    subroutine openPackage(dB,pos0)
    !!Unpack the package in the same order as it was packed
        integer, intent(in) :: pos0             !!position of the first variable
        type(dataBlock), intent(inout) :: dB    !!dataBlock
        integer :: ierr,pos
        pos = pos0
!#MPITOOLS_SCRIPT001
        MPIUNPACKREAL(dB%wavel)
        MPIUNPACKREAL(dB%mre)
        MPIUNPACKREAL(dB%mim)
        MPIUNPACKREAL(dB%ssalbedo)
        MPIUNPACKREAL(dB%ell1)
        MPIUNPACKREAL(dB%Fstop)
        MPIUNPACKREAL(dB%the0)
        MPIUNPACKREAL(dB%phi0)
        MPIUNPACKREAL(dB%phiout)
        MPIUNPACKREAL(dB%tauc)
        MPIUNPACKREAL(dB%h)
        MPIUNPACKREAL(dB%mu0)
        MPIUNPACKREAL(dB%nu0)
        MPIUNPACKREAL(dB%cphi0)
        MPIUNPACKREAL(dB%sphi0)
        MPIUNPACKREAL(dB%tauh)
        MPIUNPACKREAL(dB%albedo)
        MPIUNPACKREAL(dB%mre12)
        MPIUNPACKREAL(dB%mre21)
        MPIUNPACKREAL(dB%waven)
        MPIUNPACKREAL(dB%xell)
        MPIUNPACKREAL(dB%dthe)
        MPIUNPACKREAL(dB%nre)
        MPIUNPACKREAL(dB%nim)
        MPIUNPACKREAL(dB%rho)
        MPIUNPACKREAL(dB%predeflambda)
        MPIUNPACKREAL(dB%appFBoundary)
        MPIUNPACKINTEGER(dB%radPoints)
        MPIUNPACKINTEGER(dB%seed)
        MPIUNPACKINTEGER(dB%nrn)
        MPIUNPACKINTEGER(dB%ntot)
        MPIUNPACKINTEGER(dB%nsca)
        MPIUNPACKINTEGER(dB%np)
        MPIUNPACKINTEGER(dB%tauflg)
        MPIUNPACKINTEGER(dB%ncm)
        MPIUNPACKINTEGER(dB%ns)
        MPIUNPACKINTEGER(dB%ntheb)
        MPIUNPACKINTEGER(dB%nphib)
        MPIUNPACKINTEGER(dB%nthe)
        MPIUNPACKINTEGER(dB%nphi)
        MPIUNPACKCMPLX(dB%m12)
        MPIUNPACKCMPLX(dB%m21)
!#MPITOOLS_SCRIPT002
        call MPI_Unpack(buffer,SIZEOFBUFFER,pos,listOfInds,NUMOFARRAYSINDB,MPI_INTEGER,MPI_COMM_WORLD)
        call MPI_Unpack(buffer,SIZEOFBUFFER,pos,buffer2,SIZEOFBUFFER2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD)
    end subroutine


    !--------------------------------------------
    !subroutine subreadBuffer(ind,indPos,array)
    !Read arrays from the buffer
    !--------------------------------------------
    subroutine subreadBuffer(ind,indPos,array)
        integer, intent(inout) :: ind                           !!position of the size data in listOfInds
        integer, intent(inout) :: indPos                        !!position of the data in buffer
        real(kind=rk), allocatable, intent(inout) :: array(:)   !!array to be filled
        integer :: j1
        if(listOfInds(ind)>0) then
            allocate(array(listOfInds(ind)))
            do j1=1,size(array)
                array(j1) = buffer2(indPos+j1)
            enddo
            indPos=indPos+size(array)
           
        endif
        ind=ind+1
    end subroutine



    !--------------------------------------------
    !subroutine fillTypeListAndLengths(typelist,lengths)

    !--------------------------------------------
    subroutine fillTypeListAndLengths(typelist,lengths)
    !Fill information about types and size of the data
        type(MPI_Datatype), intent(out) :: typelist(:)  !!list of types
        integer, intent(out) :: lengths(:)              !!list of lengths
        integer :: n,j1
        j1=1
!#MPITOOLS_SCRIPT101
        do j1=1,22
            typelist(j1)=MPI_DOUBLE_PRECISION
        end do
        do j1=23,34
            typelist(j1)=MPI_INTEGER
        end do
        do j1=35,36
            typelist(j1)=MPI_DOUBLE_COMPLEX
        end do

!#MPITOOLS_SCRIPT102
        typelist(j1) = MPI_INTEGER
        j1=j1+1
        typelist(j1) = MPI_DOUBLE_PRECISION
        j1=0
!#MPITOOLS_SCRIPT201
        
!#MPITOOLS_SCRIPT202
        lengths(:)=1
        lengths(size(lengths)-1)=NUMOFARRAYSINDB
        lengths(size(lengths))=SIZEOFBUFFER2
    end subroutine


    subroutine fillDisplacements(dB,disp)
    !!Get the location of the data 
        type(dataBlock), intent(in) :: dB                           !!datablock
        integer(kind=MPI_ADDRESS_KIND),intent(inout) :: disp(:)     !!list of addresses
        integer :: j1
        j1=1
!#MPITOOLS_SCRIPT301
        MPIADDRESS(dB%wavel)
        MPIADDRESS(dB%mre)
        MPIADDRESS(dB%mim)
        MPIADDRESS(dB%ssalbedo)
        MPIADDRESS(dB%ell1)
        MPIADDRESS(dB%Fstop)
        MPIADDRESS(dB%the0)
        MPIADDRESS(dB%phi0)
        MPIADDRESS(dB%phiout)
        MPIADDRESS(dB%tauc)
        MPIADDRESS(dB%h)
        MPIADDRESS(dB%mu0)
        MPIADDRESS(dB%nu0)
        MPIADDRESS(dB%cphi0)
        MPIADDRESS(dB%sphi0)
        MPIADDRESS(dB%tauh)
        MPIADDRESS(dB%albedo)
        MPIADDRESS(dB%mre12)
        MPIADDRESS(dB%mre21)
        MPIADDRESS(dB%waven)
        MPIADDRESS(dB%xell)
        MPIADDRESS(dB%dthe)
        MPIADDRESS(dB%seed)
        MPIADDRESS(dB%nrn)
        MPIADDRESS(dB%ntot)
        MPIADDRESS(dB%nsca)
        MPIADDRESS(dB%np)
        MPIADDRESS(dB%tauflg)
        MPIADDRESS(dB%ncm)
        MPIADDRESS(dB%ns)
        MPIADDRESS(dB%ntheb)
        MPIADDRESS(dB%nphib)
        MPIADDRESS(dB%nthe)
        MPIADDRESS(dB%nphi)
        MPIADDRESS(dB%m12)
        MPIADDRESS(dB%m21)
!#MPITOOLS_SCRIPT302
        MPIADDRESS(listOfInds) 
        MPIADDRESS(buffer2)
    end subroutine



    subroutine readBuffer(dB)
    !!Read arrays from the buffer
        type(dataBlock), intent(inout) :: dB        !!dataBlock
        integer :: ind,indPos
        ind=1
        indPos=0
!#MPITOOLS_SCRIPT401
        call subReadBuffer(ind,indPos,dB%CTHEIF)
        call subReadBuffer(ind,indPos,dB%STHEIF)
        call subReadBuffer(ind,indPos,dB%INORM)
        call subReadBuffer(ind,indPos,dB%STHEI)
        call subReadBuffer(ind,indPos,dB%CTHEI)
        call subReadBuffer(ind,indPos,dB%CPHII)
        call subReadBuffer(ind,indPos,dB%SPHII)
        call subReadBuffer(ind,indPos,dB%THEBF)
        call subReadBuffer(ind,indPos,dB%CTHEBF)
        call subReadBuffer(ind,indPos,dB%STHEBF)
        call subReadBuffer(ind,indPos,dB%STHEB)
        call subReadBuffer(ind,indPos,dB%CTHEB)
        call subReadBuffer(ind,indPos,dB%CPHIB)
        call subReadBuffer(ind,indPos,dB%SPHIB)
        call subReadBuffer(ind,indPos,dB%PHIB)
!#MPITOOLS_SCRIPT402
    end subroutine



    subroutine fillBuffer2(dB)
    !!Fill buffer with arrays
        type(dataBlock), intent(in) :: dB       !!dataBlock
        integer :: n,j1,indsPos,initPos
        initPos = 0
        indsPos = 1
!#MPITOOLS_SCRIPT501
        call copyArray(dB%CTHEIF,listOfInds,initPos,indsPos)
        call copyArray(dB%STHEIF,listOfInds,initPos,indsPos)
        call copyArray(dB%INORM,listOfInds,initPos,indsPos)
        call copyArray(dB%STHEI,listOfInds,initPos,indsPos)
        call copyArray(dB%CTHEI,listOfInds,initPos,indsPos)
        call copyArray(dB%CPHII,listOfInds,initPos,indsPos)
        call copyArray(dB%SPHII,listOfInds,initPos,indsPos)
        call copyArray(dB%THEBF,listOfInds,initPos,indsPos)
        call copyArray(dB%CTHEBF,listOfInds,initPos,indsPos)
        call copyArray(dB%STHEBF,listOfInds,initPos,indsPos)
        call copyArray(dB%STHEB,listOfInds,initPos,indsPos)
        call copyArray(dB%CTHEB,listOfInds,initPos,indsPos)
        call copyArray(dB%CPHIB,listOfInds,initPos,indsPos)
        call copyArray(dB%SPHIB,listOfInds,initPos,indsPos)
        call copyArray(dB%PHIB,listOfInds,initPos,indsPos)
!#MPITOOLS_SCRIPT502
    end subroutine


    subroutine copyArray(array,listOfInds,initPos,indsPos)
        !!Copies array to buffer and fills data needed for 
        !!unpacking the them. 
        real(kind=rk), allocatable,intent(in) :: array(:)   !!array
        integer, intent(inout) :: listOfInds(:)             !!list of sizes
        integer, intent(inout) :: initPos                   !!position in buffer
        integer, intent(inout) :: indsPos                   !!position in listOfInds
        integer :: n,j1
        if(allocated(array)) then
            do j1=1,size(array)
                buffer2(initPos+j1)=array(j1)
            enddo
            initPos=initPos+size(array)
            listOfInds(indsPos) = size(array)
        else
            listOfInds(indsPos) = 0
        endif
        indsPos=indsPos+1
    end subroutine


    subroutine mergeAssembledData(aD)
    !!Merge assembledData from different MPI processes
        type(assembledData), intent(inout) :: aD        !!assembledData
        real(kind=rk), allocatable :: buffer(:)
        integer :: listOfInds(10),ierr
        integer :: rank,j1,j2,j3,j4,j5,j6,x,nodes,finished
        integer :: bufferSize
        type(MPI_STATUS) :: status0
        type(MPI_REQUEST) :: request
        request = MPI_REQUEST_NULL
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        listOfInds=(/size(aD%MRT,dim=1),size(aD%MRT,dim=2),size(aD%MRT,dim=3),size(aD%MRT,dim=4),size(aD%MRT,dim=5),&
                size(aD%MBS,dim=1),size(aD%MBS,dim=2),size(aD%MBS,dim=3),size(aD%MBS,dim=4),size(aD%MBS,dim=5)/)

        bufferSize = listOfInds(1)*listOfInds(2)*listOfInds(3)*listOfInds(4)*listOfInds(5)
        bufferSize = bufferSize + listOfInds(5)*listOfInds(6)*listOfInds(7)*listOfInds(8)*listOfInds(9)*listOfInds(10)
        bufferSize = bufferSize+7*listOfInds(5)+100 !!little extra
        !!compute buffer size
        
        allocate(buffer(bufferSize))
        !!Allocate buffer

        
        if(rank==0) then
        !!Master process collects data from other processes
            finished=0
            call MPI_COMM_SIZE(MPI_COMM_WORLD,nodes,ierr)
            do while(finished<nodes-1)
                !!Receive non blocking
                call MPI_Irecv(buffer,bufferSize,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,0, MPI_COMM_WORLD,request)               
                if(request/=MPI_REQUEST_NULL) then
                    call MPI_WAIT(request,status0)
                    x=1
                    do j5=1,listOfInds(5)
                        do j4=1,listOfInds(4)
                            do j3=1,listOfInds(3)
                                do j2=1,listOfInds(2)
                                    do j1=1,listOfInds(1)
                                        aD%MRT(j1,j2,j3,j4,j5)=aD%MRT(j1,j2,j3,j4,j5)+buffer(x)
                                        x=x+1
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                    if(allocated(ad%MBS)) then
                        do j6=1,listOfInds(5)
                            do j5=1,listOfInds(10)
                                do j4=1,listOfInds(9)
                                    do j3=1,listOfInds(8)
                                        do j2=1,listOfInds(7)
                                            do j1=1,listOfInds(6)
                                                aD%MBS(j1,j2,j3,j4,j5,j6)= &
                                                    & aD%MBS(j1,j2,j3,j4,j5,j6)+buffer(x)
                                                x=x+1
                                            enddo
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    endif
                    !!Then store other parameters
                    do j1=1,listOfInds(5)
                        aD%Aref(j1)=buffer(x)+aD%Aref(j1)
                        x=x+1
                        aD%Aspec(j1)=buffer(x)+aD%Aspec(j1)
                        x=x+1
                        aD%Aabs(j1)=buffer(x)+aD%Aabs(j1)
                        x=x+1
                        aD%Astop(j1)=buffer(x)+aD%Astop(j1)
                        x=x+1
                        aD%Adt(j1)=buffer(x)+aD%Adt(j1)
                        x=x+1
                        aD%Atra(j1)=buffer(x)+aD%Atra(j1)
                        x=x+1
                        aD%rays(j1)=buffer(x)+aD%rays(j1)
                        x=x+1
                    enddo 
                    finished=finished+1
                else
                    write(6,*) "Master Process: waiting for job"
                endif
       
            enddo
        else
            !!Send to master process
            x=1
            !!First pack RT solution
            do j5=1,listOfInds(5)
                do j4=1,listOfInds(4)
                    do j3=1,listOfInds(3)
                        do j2=1,listOfInds(2)
                            do j1=1,listOfInds(1)
                                buffer(x) = aD%MRT(j1,j2,j3,j4,j5)
                                x=x+1
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            if(allocated(aD%MBS)) then
                do j6=1,listOfInds(5)
                    do j5=1,listOfInds(10)
                        do j4=1,listOfInds(9)
                            do j3=1,listOfInds(8)
                                do j2=1,listOfInds(7)
                                    do j1=1,listOfInds(6)
                                        buffer(x) = aD%MBS(j1,j2,j3,j4,j5,j6)
                                        x=x+1
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            endif
            !!Then store other parameters
            do j1=1,listOfInds(5)
                buffer(x)=aD%Aref(j1)
                x=x+1
                buffer(x)=aD%Aspec(j1)
                x=x+1
                buffer(x)=aD%Aabs(j1)
                x=x+1
                buffer(x)=aD%Astop(j1)
                x=x+1
                buffer(x)=aD%Adt(j1)
                x=x+1
                buffer(x)=aD%Atra(j1)
                x=x+1
                buffer(x)=aD%rays(j1)
                x=x+1
            enddo
            !!Send
            call MPI_SEND(buffer,bufferSize,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD)     
        endif
    end subroutine

    function mpi_get_my_id() result(id)
        integer :: id
        integer :: ierr
        call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
    end function

    subroutine mpi_tools_bcast_Tmat(Taa,Tab,Tba,Tbb)
        complex(kind=rk), allocatable, intent(inout) :: Taa(:,:,:), Tab(:,:,:), Tba(:,:,:), Tbb(:,:,:)
        integer :: id,ierr,size1,size2,size3
        call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
        if(id==0) then
            ASSERTI(allocated(Taa),"mpi_tools: mpi_tools_bcast_Tmat: Taa must be allocated")
            ASSERTI(allocated(Tab),"mpi_tools: mpi_tools_bcast_Tmat: Tab must be allocated")
            ASSERTI(allocated(Tba),"mpi_tools: mpi_tools_bcast_Tmat: Tba must be allocated")
            ASSERTI(allocated(Tbb),"mpi_tools: mpi_tools_bcast_Tmat: Tbb must be allocated")
            ASSERTI(size(Taa,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Taa must be > 0")
            ASSERTI(size(Tab,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Tab must be > 0")
            ASSERTI(size(Tba,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Tba must be > 0")
            ASSERTI(size(Tbb,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Tbb must be > 0")
            size1=size(Taa,dim=1)
            size2=size(Taa,dim=2)
            size3=size(Taa,dim=3)
        endif
        call MPI_Bcast(size1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(size2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(size3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(id/=0) then
            allocate(Taa(size1,size2,size3), Tab(size1,size2,size3), Tba(size1,size2,size3), Tbb(size1,size2,size3))
        endif
        call MPI_Bcast(Taa,size1*size2*size3,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)    
        call MPI_Bcast(Tab,size1*size2*size3,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 
        call MPI_Bcast(Tba,size1*size2*size3,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 
        call MPI_Bcast(Tbb,size1*size2*size3,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 

        ASSERTI(allocated(Taa),"mpi_tools: mpi_tools_bcast_Tmat: Taa must be allocated")
        ASSERTI(allocated(Tab),"mpi_tools: mpi_tools_bcast_Tmat: Tab must be allocated")
        ASSERTI(allocated(Tba),"mpi_tools: mpi_tools_bcast_Tmat: Tba must be allocated")
        ASSERTI(allocated(Tbb),"mpi_tools: mpi_tools_bcast_Tmat: Tbb must be allocated")
        ASSERTI(size(Taa,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Taa must be > 0")
        ASSERTI(size(Tab,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Tab must be > 0")
        ASSERTI(size(Tba,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Tba must be > 0")
        ASSERTI(size(Tbb,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Tbb must be > 0")
    end subroutine




    subroutine mpi_tools_bcast_real_2darray(Cexts)
        real(kind=rk), allocatable, intent(inout) :: Cexts(:,:)
        integer :: id,ierr,size1
        call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
        if(id==0) then
            ASSERTI(allocated(Cexts),"mpi_tools: mpi_tools_bcast_Tmat: Taa must be allocated")
            ASSERTI(size(Cexts,dim=2)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Taa must be > 0")
            size1=size(Cexts,dim=2)
        endif
        call MPI_Bcast(size1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(id/=0) then
            allocate(Cexts(2,size1))
        endif
        call MPI_Bcast(Cexts,2*size1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)    
        ASSERTI(allocated(Cexts),"mpi_tools: mpi_tools_bcast_Tmat: Taa must be allocated")
        ASSERTI(size(Cexts,1)>0,"mpi_tools: mpi_tools_bcast_Tmat: size of Taa must be > 0")
    end subroutine




    function mpi_get_time() result(retVal)
        real(kind=rk) :: retVal
        retVal = MPI_Wtime()
    end function


    subroutine finishMPISession()
    !!Read arrays from the buffer
        integer :: ierr
        call MPI_FINALIZE(ierr)
    end subroutine


    subroutine ABORT0(msg)
    !!Make clean Abort
        integer :: ierr,errorcode
        character(*), intent(in) :: msg             !!output message
        write(6,*) trim(msg)
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)  
        stop
    end subroutine


        

end module



