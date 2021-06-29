#include "macros.inc"

module grid_holder
!!Contains grid which defines the location of the detectors
use constants
use typedefinitions
use error_handler
use rng
use math_routines
implicit none

real(kind=rk), allocatable :: KPATH(:,:)
!!Directions
real(kind=rk), allocatable :: XPATH(:,:)
!!Locations (no unit)
real(kind=rk), allocatable :: INORM(:)
!!normalization
integer :: nthe,nphi,ntheb,nphib
!!number of theta and phi angles
real(kind=rk), allocatable :: theGrid(:),theGrid_cb(:)
!!theta grid points
real(kind=rk), allocatable :: phiGrid(:),phiGrid_cb(:)
!!phi grid_points
real(kind=rk), allocatable :: CTHE_arr(:),CTHEB_arr(:)
real(kind=rk), allocatable :: STHE_arr(:),STHEB_arr(:)
!!theta grid points
real(kind=rk), allocatable :: SPHI_arr(:),SPHIB_arr(:)
real(kind=rk), allocatable :: CPHI_arr(:),CPHIB_arr(:)
!!phi grid_points
integer :: g_jsca
!!number of scattering processes
integer :: g_max_jsca
!!max number of scattering processes
real(kind=rk) :: far_field_boundary
!!distance to far field boundary(where detectors are located)

protected :: KPATH, XPATH,g_jsca,INORM,nthe,nphi,nphib,ntheb,far_field_boundary
protected :: theGrid,theGrid_cb,phiGrid,phiGrid_cb
!!Let others see these, but do not let them modify them

contains

subroutine init_grid_holder(db)
    !!Initialize module grid_holder
    type(dataBlock), intent(in) :: db

    real(kind=rk), allocatable :: X(:),W(:)
    real(kind=rk) :: phi,tmp1,tmp2
    integer :: j1
    
    ASSERTI(dB%nsca>0,"Grid_holder: Number of scattering processes has to be >0")
    ASSERTI(dB%nthe>0,"Grid_holder: Number of theta angles in detector grid has to be >0")
    ASSERTI(dB%nphi>0,"Grid_holder: Number of phi angles in detector grid has to be >0")
    if(dB%cb_on) then
        ASSERTI(dB%nphib>0,"Grid_holder: Number of phi angles in grid (cb) has to be >0")
        ASSERTI(dB%ntheb>0,"Grid_holder: Number of theta angles in grid (cb) has to be >0")
    endif
    ASSERTI(dB%ffboundary>0,"Grid_holder: Far field boundary has to be >0")
    allocate(KPATH(3,dB%nsca+1),XPATH(3,dB%nsca+1))
    nthe=db%nthe
    nphi=db%nphi
    g_max_jsca = dB%nsca
    far_field_boundary = dB%ffboundary
    allocate(theGrid(nthe),phiGrid(nphi))
    allocate(CTHE_arr(nthe),STHE_arr(nthe))
    allocate(CPHI_arr(nphi),SPHI_arr(nphi))
    allocate(X(nthe),W(nthe))
    allocate(INORM(nthe))

    !Gaussian Legengdre
    call gaussLQuads(-1.0_rk,1.0_rk,X,W)
    do j1 = 1,nthe
        INORM(j1)=W(j1)*(2.0_rk*pi/nphi)/(4.0_rk*pi)
        CTHE_arr(j1)=X(j1)
        theGrid(j1) = acos(CTHE_arr(j1))
        STHE_arr(j1)=sqrt(1.0_rk-CTHE_arr(j1)**2)    
    end do

    do j1 = 1,nphi
        phi=(j1-1)*2.0_rk*pi/nphi
        phiGrid(j1) = phi
        CPHI_arr(j1)=cos(phi)
        SPHI_arr(j1)=sin(phi)
    end do


    if(dB%cb_on) then
        ntheb=size(dB%theGrid_cb)
        nphib=db%nphib
        allocate(CTHEB_arr(ntheb),STHEB_arr(ntheb))
        allocate(SPHIB_arr(nphib),CPHIB_arr(nphib))
        allocate(theGrid_cb(ntheb),phiGrid_cb(nphib))
        !!cb related angles
        theGrid_cb = pi-dB%theGrid_cb*pi/180
        do j1 = 1,ntheb
            CTHEB_arr(j1)=cos(theGrid_cb(j1))
            STHEB_arr(j1)=sin(theGrid_cb(j1))
        end do

        do j1 = 1,nphib
            phi=(j1-1)*2.0_rk*pi/nphib
            phiGrid_cb(j1) = phi
            CPHIB_arr(j1)=cos(phi)
            SPHIB_arr(j1)=sin(phi)
        end do

    endif

    deallocate(X,W)

    if(debugging) then
        ASSERTC(allocated(XPATH),.eqv.,.true.)
        ASSERTC(allocated(KPATH),.eqv.,.true.)
        ASSERTC(allocated(INORM),.eqv.,.true.)
        ASSERTC(allocated(CPHI_arr),.eqv.,.true.)
        ASSERTC(allocated(theGrid),.eqv.,.true.)
        ASSERTC(allocated(phiGrid),.eqv.,.true.)
        ASSERTC(allocated(SPHI_arr),.eqv.,.true.)
        ASSERTC(sum(abs(theGrid)),==,sum(theGrid))
        if(dB%cb_on) then
            tmp1 = sum(abs(theGrid_cb))
            tmp2 = sum(theGrid_cb)
            ASSERTC(tmp1,==,tmp2)
            tmp1 = sum(abs(phiGrid_cb))
            tmp2 = sum(phiGrid_cb)
            ASSERTC(tmp1,==,tmp2)
            ASSERTC(allocated(CTHEB_arr),.eqv.,.true.)
            ASSERTC(allocated(STHEB_arr),.eqv.,.true.)
            ASSERTC(allocated(theGrid_cb),.eqv.,.true.)
            ASSERTC(allocated(phiGrid_cb),.eqv.,.true.)
            ASSERTC(nphib,>,0)
            ASSERTC(ntheb,>,0)
        endif
        ASSERTC(sum(abs(phiGrid)),==,sum(phiGrid))
        ASSERTC(sum(abs(INORM)),==,sum(INORM))
        ASSERTC(nthe,>,0)
        ASSERTC(far_field_boundary,>,0)
        ASSERTC(nphi,>,0)
        ASSERTC(g_max_jsca ,>,0)
    endif

end subroutine

subroutine add_location(X,K)
    !!Store location and direction of the ray
    real(kind=rk), intent(in) :: X(3) !!New location (no unit)
    real(kind=rk), intent(in) :: K(3) !!Direction vector
    if(debugging) then
        ASSERTC(g_jsca,<=,g_max_jsca)
        ASSERTC(allocated(XPATH),.eqv.,.true.)
        ASSERTC(allocated(KPATH),.eqv.,.true.)       
        ASSERTC(dot_product(K,K),>,1.0_rk-0.0001_rk)  
        ASSERTC(dot_product(K,K),<,1.0_rk+0.0001_rk) 
    endif

    g_jsca=g_jsca+1
    KPATH(:,g_jsca)=K
    XPATH(:,g_jsca)=X
end subroutine




subroutine clear_grid_holder()
    !!Clear scattering count
    g_jsca = 0
end subroutine




pure function get_direction_vector_cb(j1,j2) result(KF)
    !!Get direction vector to the detector at the grid point j1,j2 (coherent backscattering)
    integer, intent(in) :: j1,j2        !!j1 (theta), j2 (phi)
    real(kind=rk) :: KF(3)

    KF(1)=STHEB_arr(j1)*CPHIB_arr(j2)
    KF(2)=STHEB_arr(j1)*SPHIB_arr(j2)
    KF(3)=CTHEB_arr(j1)

    if(debugging) then   
        ASSERTCP(dot_product(KF,KF),>,1.0_rk-0.0001_rk)  
        ASSERTCP(dot_product(KF,KF),<,1.0_rk+0.0001_rk)    
    endif
end function

pure function get_direction_vector(j1,j2) result(KF)
    !!Get direction vector to the detector at the grid point j1,j2
    integer, intent(in) :: j1,j2
    real(kind=rk) :: KF(3)
    KF(1)=STHE_arr(j1)*CPHI_arr(j2)
    KF(2)=STHE_arr(j1)*SPHI_arr(j2)
    KF(3)=CTHE_arr(j1)
    if(debugging) then   
        ASSERTCP(dot_product(KF,KF),>,1.0_rk-0.0001_rk)  
        ASSERTCP(dot_product(KF,KF),<,1.0_rk+0.0001_rk)    
    endif
end function




end module


