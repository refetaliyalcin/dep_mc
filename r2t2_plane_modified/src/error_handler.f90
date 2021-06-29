module error_handler
!!a simple error handler. The name is not self explanatory because
!!previously all the errors went through addError but I wanted to go with assertions
!!so the job of the function changed during the development
!!
!!Copyright (C) 2016 Timo Väisänen and University of Helsinki
!!All rights reserved.
!!The new BSD License is applied to this software, see LICENSE.txt
    use constants
    implicit none


    integer :: errorCount = 0                       !!number of failed assertion encountered
    real(kind=rk), parameter :: a_diff = 0.00001_rk !!variance that is allowed
    private errorCount

    interface toString
        !!Use this to convert int,real,cmplx,logical to string
        module procedure intToString, realToString, complexToString, logicalToString
    end interface

contains

!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine addError(msg)
        !!Old function which was used to increase error counter present in the module
        character(*), intent(in) :: msg     !message to be displayed
        errorCount=errorCount+1
        write(6,*) "Error in module: "//trim(msg)
    end subroutine

    !-----------------------------
    !function checkIfErrors() result(retVal))
    !
    !Check whether the error count is /= 0
    !-----------------------------
!///////////////////////////////////////////////////////////////////////////////////////////
    function check_if_errors() result(retVal)
       !!Returns true if error counter in the module is not zero
       logical :: retVal                    
       retVal=.false.
       if(errorCount/=0) retVal=.true.
    end function

!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine assert_c(cond,msg)
        !!One of the first versions of the assert
        logical, intent(in) :: cond
        character(*), intent(in) :: msg 
        if(.not. cond) call addError(msg)
    end subroutine


!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine assert(cond,msg)
        !!Use this assert in the routines which are not pure
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        real(kind=rk) :: r 
        if(.not. cond) write(6,*) adjustl(trim(msg))
        !if(.not. cond) r=r/0.0_rk
        if(.not. cond) stop
    end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine assert2(cond,msg,line_num,msg2)
        !!Use this with the macro defined in the "macros.inc"
        logical, intent(in) :: cond             !!condition
        character(*), intent(in) :: msg,msg2    !!name of the module and message
        integer, intent(in) :: line_num         !!line number of the assert
        real(kind=rk) :: r 
        if(.not. cond) write(6,*) "Error in: ",adjustl(trim(msg)),&
            &   " line "//trim(adjustl(toString(line_num)))//". Expected: ",adjustl(trim(msg2))
        !if(.not. cond) r=r/0.0_rk
        if(.not. cond) stop
    end subroutine


!///////////////////////////////////////////////////////////////////////////////////////////
    pure subroutine assert2_pure(cond,msg,line_num,msg2)
        !!when the function is pure
        logical, intent(in) :: cond
        character(*), intent(in) :: msg,msg2
        integer, intent(in) :: line_num
        real(kind=rk) :: r 
        ! if(.not. cond) write(6,*) "Error in: ",adjustl(trim(msg)),&
        !    &   " line "//trim(adjustl(toString(line_num)))//". Expected: ",adjustl(trim(msg2))
        if(.not. cond) r=r/0.0_rk
        !if(.not. cond) stop
    end subroutine



!///////////////////////////////////////////////////////////////////////////////////////////
    !!Different kind of conversion functions which are used with a
    function intToString(val) result(retVal)
        integer, intent(in) :: val
        character(len=100) :: retVal
        write(retVal,'(I5)') val 
    end function

    function realToString(val) result(retVal)
        real(kind=rk), intent(in) :: val
        character(len=100) :: retVal
        write(retVal,'(E14.5)') val 
    end function

    function complexToString(val) result(retVal)
        complex(kind=rk), intent(in) :: val
        character(len=100) :: retVal
        write(retVal,'(E14.5, E14.5)') val 
    end function

    function logicalToString(val) result(retVal)
        logical, intent(in) :: val
        character(len=100) :: retVal
        write(retVal,'(L1)') val 
    end function

!///////////////////////////////////////////////////////////////////////////////////////////
    function contains_elem(elem,arr) result(retVal) 
        !!Check whether the array(arr) contains the element(elem)
        real(kind=rk), intent(in) :: arr(:)
        real(kind=rk), intent(in) :: elem
        logical :: retVal
        integer :: j1
        retVal = .true.
        do j1=1,size(arr)
            if(arr(j1)==elem) return
        enddo
        retVal = .false.
    end function

!///////////////////////////////////////////////////////////////////////////////////////////
    function is_cumulative(arr) result(retVal)
        !!Check that the cdf is valid
        real(kind=rk), intent(in) :: arr(:)
        integer :: j1
        logical :: retVal
        if(size(arr)<2) return
        retVal = .false.
        do j1=2,size(arr)
            if(arr(j1-1)>arr(j1)) return
        enddo 
        retVal = .true.
    end function


!///////////////////////////////////////////////////////////////////////////////////////////
    pure subroutine assert_pure(cond,msg)
        !!Assertion for pure function
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        real(kind=rk) :: r 
        if(.not. cond) r=r/0.0_rk
    end subroutine


!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine stop_program_gracefully()
        !!Still I'm not sure should I use this kind of solution or not
        !!It is just a reminder here
	    write(6,*) "Exit ""Gracefully"""
        stop !:) "gracefully"

    end subroutine


end module
