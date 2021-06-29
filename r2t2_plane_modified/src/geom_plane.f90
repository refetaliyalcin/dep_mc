#include "macros.inc"

!-------------------------------------------------------------------
! Defines sphere geometry
!
! Copyright (C) 2018 Timo Väisänen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 4.4.2018 TV: Cleaning
!
! TODO:
! Apply support for different volume element size
!-----------------------------

module geometry
    use constants
    use typedefinitions
    use rng
    implicit none

    real(kind=rk) :: rad_ve         !!Radius of the VE      [optical depth]
    real(kind=rk) :: threshold      !!too deep inside the media

    real(kind=rk) :: ell            !!incoherent mean free path

    real(kind=rk) :: forced_rad     !!To keep the VE inside the medium
    logical       :: overlap_prob   !!overlap is allowed sometimes, otherwise always         
    real(kind=rk) :: surface_att    !!the attenuation starts from the surface
    real(kind=rk) :: planeN(3),R(3,3)

    private
    public :: geometry_output_data1,init_geometry
    public :: randomize_initial_location, not_too_far
    public :: propagate,propagate_in,distance_to_surface,distance_to_surface2
contains

!///////////////////////////////////////////////////////////////////////////////////////////
function geometry_output_data1() result(retVal)
    !!Return the scatterer type as string
    character(64) :: retVal
    retVal = "Sphere geometry"
end function


!///////////////////////////////////////////////////////////////////////////////////////////
!!Initialize geometry using data from dataBlock (db)
subroutine init_geometry(dB)
    type(dataBlock), intent(in) :: dB
    ASSERTI(dB%volrad>=0,"The radius of the VE must be > 0 ")
    ASSERTI(dB%ell>0,"Mean free path must be be > 0")
    ASSERTI(dB%tau_c>0,"The threshold must be over > 0")
    ASSERTI(dB%the0>90,"theta angle must be over 90%")

    rad_ve =  dB%volrad/dB%ell
    threshold = dB%tau_c

    ell = dB%ell
    forced_rad = 0.0_rk
    surface_att = 0.0_rk

    planeN = plane(pi-dB%the0*pi/180.0_rk)
    R = rotation_x(pi-dB%the0*pi/180.0_rk)

    overlap_prob = dB%overlap_prob
    if(dB%surface_att) surface_att = rad_ve
    if(db%forced_vm) forced_rad = rad_ve
    
end subroutine

pure function plane(the) result(K)
    real(kind=rk), intent(in) :: the
    real(kind=rk) :: K(3)
    K(1)=0.0_rk
    K(2)=sin(the)
    K(3)=cos(the)
end function


pure function distance_to_plane(X,K,N) result(dx)
    real(kind=rk), intent(in) :: X(3),K(3),N(3)
    real(kind=rk) :: dX,denom,p(3),t
    dX=huge(dX)
    denom = dot_product(K,N)
    if(denom > epsilon(denom)) then
        p = -X
        t = dot_product(p,N)/denom
        if(t < 0) return
        dX = t
    endif
end function


!////////////////////////////////////////////
!!Distance to surface from X at the direction of K
function distance_to_surface(X,K) result(dX)
    real(kind=rk), intent(in) :: X(3)   !!Current location      [optical depth]
    real(kind=rk), intent(in) :: K(3)   !!Propagation direction 
    real(kind=rk) :: dX                 !!distance to surface   [optical depth]
    real(kind=rk) :: kx, xx
    real(kind=rk) :: tmpX(3)
    dX = distance_to_plane(X,K,planeN)
    if(dX>surface_att) then
        dX = dX-surface_att
    endif

    !!Asserts, turn them off if not used
    if(debugging) then
        ASSERTC(dX,>=,0.0_rk)
    endif 
end function


!///////////////////////////////////////////////
!!Distance from the entry  to the exit point 
!!Use this with the direct transmission
function distance_to_surface2(X,K) result(dX)
    real(kind=rk), intent(in) :: X(3)   !!Current location      [optical depth]
    real(kind=rk), intent(in) :: K(3)   !!Propagation direction 
    real(kind=rk) :: dX                 !!distance to surface   [optical depth]
    real(kind=rk) :: kx, xx
    dX = 999999.0_rk    
    !!Asserts, turn them off if not used
    if(debugging) then
        ASSERTC(dX,>=,0.0_rk)
    endif 
end function

function rotation_x(the) result(R)
    real(kind=rk), intent(in) :: the
    real(kind=rk) :: R(3,3)
    R(:,1)=(/1.0_rk,0.0_rk,0.0_rk/)
    R(:,2)=(/0.0_rk,cos(the),sin(the)/)
    R(:,3)=(/0.0_rk,-sin(the),cos(the)/)
end function

!////////////////////////////////////////////
!!!Check whether light has gone too deep
function not_too_far(X) result(bool)
    real(kind=rk), intent(in) :: X(3)  !!Current location of the ray [optical depth]
    logical :: bool                    !!return value
    real(kind=rk) :: D
    D = dot_product(planeN,X)
    bool = (D<threshold)
end function





!///////////////////////////////////////////////////////////////////////////////////////////
subroutine propagate_in(X,K)
    !!Propagates into the spherical media
    real(kind=rk), intent(inout) :: X(3)   
    !!in: entry point, out: next scattering location (X unit)
    real(kind=rk), intent(in) :: K(3)      
    !!propagation direction (no unit)
    real(kind=rk) :: kx,xx,tauc,t,rn
    t=-log(randNum())
    X(:)=X(:)+K(:)*t
end subroutine


!//////////////////////////////////////////////////////////////
!!Propagates the ray within the spherical scattering 
!!and absorbing medium. The subroutine will tell if the
!!ray escaped the system through the variable "escapes"
subroutine propagate(X,K,escapes)
    real(kind=rk), intent(inout) :: X(3)   !!entry point, out: next scattering location [OPTICAL DEPTH]
    real(kind=rk), intent(in) :: K(3)      !!propagation direction [OPTICAL DEPTH]
    logical, intent(out) :: escapes        !!did the ray escape?
    real(kind=rk) :: t
    real(kind=rk) :: Y(3),yy
    real(kind=rk) :: ovolume,dX

    escapes=.false.

    t=-log(randNum())
    
    !!Checks whether the volume elements overlap
    if(overlap_prob) then
        if(t<2*rad_ve) then
            ovolume=1/12.0_rk*pi*(4*rad_ve+t)*(2*rad_ve-t)**2
            if(randNum()>(1-(ovolume/(4.0_rk/3.0_rk*pi*rad_ve**3)))) then
                escapes = .true.
                return
            endif
        endif
    endif

    dX = distance_to_plane(X,K,planeN)

    if(t<dX-forced_rad) then
        X(:)=X(:)+K(:)*t
    else
        escapes=.true.
    endif

end subroutine


!///////////////////////////////////////////////////////
!!Get the location where the ray enters the media
function randomize_initial_location() result(X)
    real(kind=rk) :: X(3)    !!The entering point (X unit)

    X(1)=0.0_rk
    X(2)=0.0_rk
    X(3)=0.0_rk-forced_rad

    X = matmul(R,x)

end function




end module 
