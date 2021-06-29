#include "macros.inc"
module misc
!!Some miscellaneous init,reset,update functions.
!!
!!Copyright (C) 2016 Karri Muinonen, Timo Väisänen and University of Helsinki
!!All rights reserved.
!!The new BSD License is applied to this software, see LICENSE.txt

    use constants
    use typedefinitions
    use error_handler
    implicit none
contains

!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine collectData(aD, accD, pol_state, ind)
        !!collect data from accD to aD
        type(assembledData), intent(inout) :: aD    !!Collect temporary results to here
        type(accessibleData), intent(in) :: AccD    !!temporary results from the algorithm
        integer, intent(in) :: pol_state, ind       !!number of rays
        real(kind=rk) :: tmp
        logical, allocatable :: arr(:,:,:,:,:,:),arr2(:,:,:,:,:)
        integer :: j1
        real(kind=rk) :: tmp_val

        ASSERTI(ind>0,"Module misc: collectData: ind must be >0")
        ASSERTI(pol_state>0,"Module misc: collectData: pol_state must be >0")
        ASSERTI(pol_state<=6,"Module misc: collectData: pol_state must be <=6")
        ASSERTI(accD%rays>=0,"Module misc: collectData: rayCount > 0")
        tmp = (accD%Aref+accD%Aspec+accD%Aabs+accD%Astop+accD%Adt+accD%Atra)/accD%rays
        ASSERTI(tmp<=1.000001_rk,"Module misc: collectData: Conservation laws broken, sum(A)>1")
        ASSERTI(tmp>=0.999999_rk,"Module misc: collectData: Conservation laws broken, sum(A)<1")
        ASSERTI(accD%Aref>=0,"Module misc: collectData: Negative reflection")
        ASSERTI(accD%Aspec>=0,"Module misc: collectData: Negative specular reflection")
        ASSERTI(accD%Aabs>=0,"Module misc: collectData: Negative absorption")
        ASSERTI(accD%Astop>=0,"Module misc: collectData: Negative waste intensity")
        ASSERTI(accD%Adt>=0,"Module misc: collectData: Negative direct transmission")
        ASSERTI(accD%Atra>=0,"Module misc: collectData: Negative transmitted energy")
        if(debugging) then
            ASSERTC(size(aD%Aref),>,0)
            ASSERTC(size(aD%Aspec),>,0)
            ASSERTC(size(aD%Aabs),>,0)
            ASSERTC(size(aD%Astop),>,0)
            ASSERTC(size(aD%Adt),>,0)
            ASSERTC(size(aD%Atra),>,0)
        endif

        aD%Aref(ind)=aD%Aref(ind)+accD%Aref
        aD%Aspec(ind)=aD%Aspec(ind)+accD%Aspec
        aD%Aabs(ind)=aD%Aabs(ind)+accD%Aabs
        aD%Astop(ind)=aD%Astop(ind)+accD%Astop
        aD%Adt(ind)=aD%Adt(ind)+accD%Adt
        aD%Atra(ind)=aD%Atra(ind)+accD%Atra      
        aD%rays(ind) = accD%rays
        aD%MRT(:,:,:,:,ind)=aD%MRT(:,:,:,:,ind)+accD%MRT(:,:,:,:)

        if(allocated(aD%MBS)) aD%MBS(:,:,:,:,:,ind)=aD%MBS(:,:,:,:,:,ind)+accD%MBS(:,:,:,:,:)


        if (debugging) then
            if(allocated(aD%MBS)) then
                allocate(arr(size(aD%MBS,1),size(aD%MBS,2),&
                & size(aD%MBS,3),size(aD%MBS,4),size(aD%MBS,5),size(aD%MBS,6)))
                arr = .false.
                where(aD%MBS>HUGE(aD%MBS(1,1,1,1,1,1))) arr = .true. !!check infinity
                ASSERTC(count(arr),==,0)
                where(aD%MBS /= aD%MBS) arr = .true.  !!check NaN
                ASSERTC(count(arr),==,0)
                deallocate(arr)
            endif
            allocate(arr2(size(aD%MRT,1),size(aD%MRT,2),&
            & size(aD%MRT,3),size(aD%MRT,4),size(aD%MRT,5)))
            arr2=.false.
            where(aD%MRT>HUGE(aD%MRT(1,1,1,1,1))) arr2 = .true. !!check infinity
            ASSERTC(count(arr2),==,0)
            where(aD%MRT /= aD%MRT) arr2 = .true.  !!check NaN
            ASSERTC(count(arr2),==,0)
            deallocate(arr2)
        endif 
    end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
    pure subroutine setIKEHEV(I,K,EH,EV,I1,K1,EH1,EV1)
    !!Sets the Stokes vector, wave vector, and the transversal basis vectors.
        real(kind=rk), intent(out) :: I(4),K(3),EH(3),EV(3)
        real(kind=rk), intent(in) :: I1(4),K1(3),EH1(3),EV1(3)
        if(debugging) then
            ASSERTCP(dot_product(K1,K1),<=,1+a_diff)
            ASSERTCP(dot_product(K1,K1),>=,1-a_diff)
            ASSERTCP(I1(1),>=,0.0_rk)
        endif

        I(:)=I1(:)
        K(:)=K1(:)
        EH(:)=EH1(:)
        EV(:)=EV1(:)
    end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine clean_assembled_data(aD)
    !!Clean asssembledData between test run and the proper run
        type(assembledData), intent(inout) :: aD !!zero variables/arrays in this aD
        ASSERTI(allocated(aD%MRT),"MRT not allocated")
        aD%Aref = 0.0_rk
        aD%Atra = 0.0_rk
        aD%Adt = 0.0_rk
        aD%Aabs = 0.0_rk
        aD%Astop = 0.0_rk
        aD%Aspec = 0.0_rk
        aD%MRT = 0.0_rK
        ad%rays = 0
        if(allocated(aD%MBS)) aD%MBS = 0.0_rk
    end subroutine



!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine resetAccessibleData(accD)
    !!Zero data in accessibleData(accD)
        type(accessibleData), intent(inout) :: accD !!zero variables/arrays in this accD
        ASSERTI(allocated(accD%MRT),"MRT not allocated")
        accD%Aref = 0.0_rk
        accD%Atra = 0.0_rk
        accD%Adt = 0.0_rk
        accD%Aabs = 0.0_rk
        accD%Astop = 0.0_rk
        accD%Aspec = 0.0_rk
        accD%MRT = 0.0_rk
        if(allocated(accD%MBS)) accD%MBS = 0.0_rk
    end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
    subroutine initAccessibleData(accD,aD)
    !!Initialize accessibleData
        type(accessibleData), intent(inout) :: accD !!zero variables/arrays in this accD
        type(assembledData), intent(in) :: aD
        integer :: ntheb, nphib,nthe,nphi
        ASSERTI(size(aD%MRT,dim=3)>0,"Module misc: initAccessibleData: size of MRT > 0")
        ASSERTI(size(aD%MRT,dim=4)>0,"Module misc: initAccessibleData: size of MRT > 0")
        nthe=size(aD%MRT,dim=3)
        nphi=size(aD%MRT,dim=4)
        allocate(accD%MRT(4,4,nthe,nphi))

        if(allocated(aD%MBS)) then
            ASSERTI(size(aD%MBS,dim=4)>0,"Module misc: initAccessibleData: size of MBS > 0")
            ASSERTI(size(aD%MBS,dim=5)>0,"Module misc: initAccessibleData: size of MBS > 0")
            ntheb=size(aD%MBS,dim=4)
            nphib=size(aD%MBS,dim=5) 
            allocate(accD%MBS(4,4,5,ntheb,nphib))
        endif

        if(debugging) then
            if(allocated(aD%MBS)) then 
                ASSERTC(ntheb,>,0)
                ASSERTC(nphib,>,0)    
            endif
            ASSERTC(nthe,>,0)
            ASSERTC(nphi,>,0)
        endif
    end subroutine


end module
