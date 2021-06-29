program RTENGINE
    use input_reader
    use typedefinitions
    use error_handler
    use preparedata
    use algorithm
    use output
    use mpi_tools
    use misc
    use constants
    use grid_holder
    implicit none

    type(dataBlock) :: dB       !!Transfer variables into modules
    type(assembledData) :: aD   !!Store collected data
    integer :: id               !!id of the process
    integer :: ierr             !!integer of the error
    integer :: MPIprocs         !!number of MPI processess
    real(kind=rk) :: t1,t2      !!estimate time    
    integer :: myShare          !!get number of rays
    real(kind=rk) :: timeFactor !!ID 0 is under heavier load than other processes so let it
    real(kind=rk), parameter :: timeFactorID0 = 0.9 !!finish little earlier

    !!Initialize MPI    
    call initializeMPI(id,ierr,MPIprocs)
    

    !master process initializes, other processes wait
    if(id==0) then   
        write(6,*) str_line 
        write(6,*) "RT-Engine"
        write(6,*) "Start initialization"  
        call readInput(dB)
        if(check_if_errors()) call abort0("Abort: read_input")
        call prepare(dB)  
        if(check_if_errors()) call abort0("Abort: prepare")
        call send_data_DB(dB) 
        if(check_if_errors()) call abort0("Abort: send_data_DB")
    else
        call receive_data_DB(dB)
        if(check_if_errors()) call abort0("Abort: receive_data_DB")
    endif
    
    call prepare2(dB,aD,id)
    !!make final preparations
    if(check_if_errors()) call abort0("Abort: prepare2")
    if(id==0) write(6,*) "Finished initialization" 
    if(id==0) write(6,*) "RT-engine will run",dB%runtime," hours"
    if(id==0) write(6,*) "Execution will stop after",dB%ntot," rays"
    if(id==0) write(6,*) str_line 
    if(id==0) call writeInfoScreen(dB,MPIprocs,6,.false.)
    
    !!estimate time
    if(id==0 .and. (dB%estimateTime)) then
        write(6,*) "Estimating time using 10 rays (not accurate for large media)..."
        t1 = mpi_get_time();
        call start(dB%run_table,aD,10,5.0*60.0_rk,id)
        t2 = mpi_get_time();
        if(t2-t1>=5.0*60.0_rk) then 
            write(6,*) "WARNING! Time estimator couldn't finish 10 rays in 5 minutes"
            write(6,*) "Time estimate : at least", (t2-t1)/size(dB%run_table)/(aD%rays*dB%ntot*60.0_rk*MPIprocs),"minutes"
        else
            write(6,*) "Time estimate: ", (t2-t1)/dB%ntot*size(dB%run_table)/(10*60.0_rk*MPIprocs),"minutes"
        endif
        write(6,*) str_line
        call clean_assembled_data(aD)   
    endif   
    
    !!start algorithm and keep track of time usage
    if(id==0) t1 = mpi_get_time();
    myShare = getMyShare(dB%ntot,MPIprocs,id,dB%estimateTime)
    if(debugging) write(6,*) "ID:", id, ", my share:",myShare
    timeFactor = 1.0_rk
    if(id==0) timeFactor = timeFactorID0
    call start(dB%run_table,aD,myShare,dB%runtime*60.0_rk*60.0_rk*timeFactor,id)
    
    !!collect data
    call mergeAssembledData(aD)

    if(id==0) t2 = mpi_get_time()
    if(id==0) write(6,*) "Time: ", (t2-t1)/60.0_rk

    !!print output
    if(id==0) call printOutput(aD,dB)

    !!check if there were any errors
    if(check_if_errors()) call abort0("Abort: main algorithm. output was saved anyway")
    !!make clean finish
    call finishMPISession()
    stop

contains


function getMyShare(rays,nprocs,id,estimatedTime) result(myShare)
!!distributes the jobs to nodes equally
    integer, intent(in) :: rays,nprocs,id
    logical, intent(in) :: estimatedTime
    real(kind=rk), parameter :: id0load = 0.75_rk
    real(kind=rk) :: mult
    integer :: myshare

    mult=1.0_rk
    if(estimatedTime) mult=id0load

    myshare = int(rays/nprocs)
    if(id==0) then
        myshare=int(myshare*id0load*mult)
    else
        myshare=int(myshare*(1+(1.0_rk-id0load*mult)/nprocs))
    endif

    if(myshare<=0) myshare = 1 

end function

end program
