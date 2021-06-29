#include "macros.inc"
!-------------------------------------------------------------------
!! Overall algorithm, see PLOS ONE paper for longer explanation.
!!  Contains parts which are not dependent on the geometry or scatterer
!!
!! Copyright (C) 2018 Timo Väisänen, Karri Muinonen, University of Helsinki
!! All rights reserved.
!! The new BSD License is applied to this software, see LICENSE.txt
!!
!! Modifications:
!! 4.4.2018 TV: Maintenance, cleaning, commenting
!!-----------------------------
module algorithm

    use constants
    use typedefinitions
    use math_routines
    use scatterer
    use geometry
    use misc
    use rng
    use grid_holder
    use absorption_module
    use error_handler
    use mpi_tools

    implicit none

    real(kind=rk) :: ell            !!mean free path
    real(kind=rk) :: threshold      !!threshold for propagation
    logical :: cb_on                !!compute coherent backscattering
    real(kind=rk) :: K0(3)          !!Initial direction 
    real(kind=rk) :: Fstop          !!Stopping intensity
    real(kind=rk) :: waven          !!wavenumber
    integer :: nsca                 !!number of consequent scattering processes
    real(kind=rk) :: mu0,nu0,sphi0,cphi0

    private
    public :: init_algorithm,start

contains



!//////////////////////////////////////////
!!Initialize module
subroutine init_algorithm(dB)
    type(dataBlock), intent(in) :: dB               
    real(kind=rk) :: mu,nu,cphi,sphi    

    ASSERTI(dB%tau_c>0.0_rk,"Algorithm init: Threshold(threshold) has to be >0")
    ASSERTI(dB%ell>=0.0_rk,"Algorithm init: Mean free path has to be >0")
    ASSERTI(dB%nsca>=0,"Algorithm init: number of scattering processes has to be >0")
    ASSERTI(dB%wavel>0.0_rk,"Algorithm init: wavelength has to be >0")
    ASSERTI(dB%Fstop>0.0_rk,"Algorithm init: Fstop has to be over >0")

    threshold = dB%tau_c
    cb_on = dB%cb_on        
    ell = dB%ell
    nsca = dB%nsca
    waven = 2*pi/dB%wavel
    Fstop = dB%Fstop    

    mu0=-1.0_rk
    nu0=0.0_rk
    cphi0=cos(dB%phi0*pi/180)
    sphi0=sin(dB%phi0*pi/180)


    K0(:) = 0.0_rk
    K0(3) = 1.0_rk

    call rayToNorm(K0,mu0,nu0,cphi0,sphi0)

    !!prepare incident direction vector which stays constant
    ASSERTC(dot_product(k0,k0),==,1.0_rk)
    ASSERTC(threshold,>,0.0_rk)
    ASSERTC(ell,>,0.0_rk)
    ASSERTC(nsca,>,0)
    ASSERTC(waven,>,0)
    ASSERTC(Fstop,>,0)
    ASSERTC(abs(mu0),<=,1)
    ASSERTC(abs(nu0),<=,1)        
    ASSERTC(abs(cphi0),<=,1)
    ASSERTC(abs(sphi0),<=,1) 

end subroutine



!///////////////////////////////////////////////////////////////////////////////////////////
!!Start algorithm
subroutine start(run_table,aD,ray_max,runtime,id)
    type(assembledData), intent(inout) :: aD        !!Results are stored here
    integer, intent(in) :: id                       !!id of the process
    type(accessibleData) :: accD                    !!temporary store for results
    integer, intent(in) :: run_table(:)             !!Contains data of initial Stoke's configurations
    real(kind=rk), intent(in) :: runtime
    integer, intent(in) :: ray_max
    real(kind=rk) :: t1,t2
    real(kind=rk) :: IS(6,4)
    integer :: j1,ps_count

    ASSERTI(ray_max>0,"Algorithm start: Number of rays has to be >0")
    ASSERTC(size(run_table),>,0)
    ASSERTC(id,>=,0)

    !initial configurations
    IS(1,:) = (/1.0_rk,1.0_rk,0.0_rk,0.0_rk/)
    IS(2,:) = (/1.0_rk,-1.0_rk,0.0_rk,0.0_rk/)
    IS(3,:) = (/1.0_rk,0.0_rk,1.0_rk,0.0_rk/)
    IS(4,:) = (/1.0_rk,0.0_rk,-1.0_rk,0.0_rk/)
    IS(5,:) = (/1.0_rk,0.0_rk,0.0_rk,1.0_rk/)
    IS(6,:) = (/1.0_rk,0.0_rk,0.0_rk,-1.0_rk/)

    ps_count = 0
    call initAccessibleData(accD,aD)
    !!run the algorithm with all configurations
    do j1=1,size(run_table)
        ps_count=ps_count+1      
        call resetAccessibleData(accD)  !!init albedos and empty arrays"
        t1 = mpi_get_time() 
        call MS(accD,IS(run_table(j1),:),run_table(j1),runtime/size(run_table),ray_max)  
        t2 = mpi_get_time() 
        !run algorithm with selected Stoke's configuration        
        write(6,'(A,2X,I4,2X,A,2X,I1,2X,A,2X,F5.3,2X,A,2X,F5.1,2X,A,2X,I6)') "process:",id,"Stoke's:"&
        &   ,j1,"Aref",accD%Aref/accD%rays,"time(min)",(t2-t1)/60.0_rk,"rays",accD%rays
        !save results
        call collectData(aD,accD,j1,ps_count)            
    enddo
end subroutine


!///////////////////////////////////////////////////////////////////////////////////////////
!!First direct transmission is computed and after that the ray propagates into media.
!!Then compute first scattering, but do not change the ray's location.
!!Compute absorption, the scattering with coherent backscattering
!!take peel-off, change the place of the scattered ray and normalize the intensity
!!of the current ray. Go back to step where next location is generated and continue
!!until the intensity drops below Stopping intensity
subroutine MS(accD,I0,polarization,runtime,ray_max)
    real(kind=rk), intent(in) :: I0(4)          !!initial Stoke's vector
    integer, intent(in) :: ray_max
    type(accessibleData), intent(inout) :: accD !!collected results
    integer, intent(in) :: polarization         !!Polarization
    real(kind=rk), intent(in) :: runtime
    real(kind=rk) :: I(4),X(3),K(3)
    real(kind=rk) :: EH(3),EV(3),EH1(3),EV1(3)  !!obselete artifacts maybe
    real(kind=rk) :: I1(4),X1(3),K1(3)
    real(kind=rk) :: I00(4)
    real(kind=rk) :: nk,N(3),rtsum,Atmp,norm
    real(kind=rk) :: X0(3),t,next,albedo
    real(kind=rk) :: phi,r,tau,tau_c,exp_a,b
    real(kind=rk) :: start_time,end_time,time_limit_sec
    integer :: jsca,ind,totref,j1,iter,nphib,ntheb,nthe,nphi,tauflg,j2,best,num0,numI,ierr,ray
    logical :: escaped

    ASSERTI(ray_max>0,"Algorithm start: The number of rays has to be>0")
    ASSERTC(I0(1),>=,0.0_rk)
    ASSERTC(dot_product(k0,k0),==,1.0_rk)
    start_time = mpi_get_time()
    ray = 0
    do while(ray<ray_max)
        call scatterer_preparations(polarization,0,0,ray)
        ray = ray+1
        !!Start tracing
        if(debugging) write(6,*) "Ray:",ray
        call clear_grid_holder()
        !!Clear grid_holder data
        jsca=0
        
        !!set initial direction, inten        
        call setIKEHEV(I,K,EH,EV,I0,K0,EH1,EV1)
        X = randomize_initial_location()
 
        !Direct transmission
        t= distance_to_surface2(X,K)

        if (t<=threshold) then
            exp_a=exp(-t)
        else
            exp_a=0.0_rk
        endif 
        do j1=1,4
            I1(j1)=exp_a*I(j1)
            I(j1)=I(j1)-I1(j1)
        enddo
        accD%Adt=accD%Adt+I1(1)
        !initial propagation
        call propagate_in(X,K)
        
        !Propagate conditionally:
        do while(I(1)>=Fstop .and. jsca<nsca)
            call add_location(X,K)
            jsca=jsca+1   
            
            call scatter0(I,albedo)

            escaped=.true.
            
            !Dependent subsequent scattering and propagation:
            do while(escaped)
                I1(:)=I(:); X1(:)=X(:); K1(:)=K
                call generate_direction(I1,K1,EH1,EV1)
                call propagate(X1,K1,escaped)
            enddo

            !Absorption:
            call ABSORB(I,accD%Aabs,albedo)

            !Bypass due to large optical depth:
            if (not_too_far(X)) then

                !!Peel-off, Diffuse reflection:
                rtsum=0.0_rk
                call peel_off(I,X,K,Atmp,accD)
                !Coherent backscattering:
                if(cb_on) call coherent_backscattering(X,K,K0,I0,I,accD)

                accD%Aref=accD%Aref+Atmp   
                rtsum=rtsum+Atmp
                !!Maintain energy conservation:
                I(1)=I(1)-rtsum
            endif
            !!Set the new position and Stokes vector generated earlier:
            norm=I1(1)
            I1(:)=I(1)*I1(:)/norm
            if(debugging) then
                ASSERTC(I1(1),<=,I(1))
            endif
            call setIKEHEV(I,K,EH,EV,I1,K1,EH1,EV1)
            X(:)=X1(:)  
            
            if(mod(jsca,20)==0) then
                end_time = mpi_get_time()
                if((end_time-start_time)>runtime) exit
            endif

        enddo
        end_time = mpi_get_time()
        accD%Astop=accD%Astop+I(1)
        if((end_time-start_time)>runtime) exit
    end do
    accD%rays = ray
end subroutine



!///////////////////////////////////////
!!Handles peel off strategy
subroutine peel_off(I0,X,K,Atmp,accD)
    real(kind=rk), intent(in) :: I0(4)              !!current Stoke's vector
    real(kind=rk), intent(in) :: X(3)               !!current location of the ray
    real(kind=rk), intent(in) :: K(3)               !!incident direction
    real(kind=rk), intent(out) :: Atmp              !!collected intensity
    type(accessibleData), intent(inout) :: accD     !!collected data

    integer :: j1,j2
    real(kind=rk) :: next,dX,K1(3),I1(4),S(4,4)

    !!Asserts, turn them off if not used    
    if(debugging) then
        ASSERTC(I0(1),>=,0.0_rk)
        ASSERTC(dot_product(K,K),<=,(1.0_rk+a_diff))
        ASSERTC(dot_product(K,K),>=,(1.0_rk-a_diff))
    endif

    Atmp=0.0_rk
    do j2=1,nphi
        do j1=1,nthe
            K1(:) = get_direction_vector(j1,j2)
            
            dX = distance_to_surface(X,K1)

            if(dX<=threshold) then
                call get_SMAT(S,I0,K1,K,j1,j2)
                next=exp(-dX)
                S(:,:)=S(:,:)*next
                accD%MRT(:,:,j1,j2)=accD%MRT(:,:,j1,j2)+S(:,:) 
                Atmp=Atmp+S(1,1) 
            endif
        enddo
    enddo

    ASSERTC(Atmp,>=,0.0_rk)
end subroutine


!////////////////////////////////////
!!Handles peel off strategy
subroutine First_order_scattering(I0,X,K,accD)
    real(kind=rk), intent(in) :: I0(4)          !!current Stoke's vector
    real(kind=rk), intent(in) :: X(3)           !!current location of the ray
    real(kind=rk), intent(in) :: K(3)           !!incident direction
    type(accessibleData), intent(inout) :: accD !!collected data
    real(kind=rk) :: next,dX,K1(3),S(4,4)       
    integer :: j1,j2

    !!Asserts, turn them off if not used    
    if(debugging) then
        ASSERTC(I0(1),>=,0.0_rk)
        ASSERTC(dot_product(K,K),<=,(1.0_rk+a_diff))
        ASSERTC(dot_product(K,K),>=,(1.0_rk-a_diff))
    endif

    do j2=1,nphib
        do j1=1,ntheb
            K1(:) = get_direction_vector_cb(j1,j2)
            
            dX = distance_to_surface(X,K1)
            
            if(dX<=threshold) then
                call get_SMAT_cb(S,I0,K1,K,j1,j2)
                next=exp(-dX)
                S(:,:)=S(:,:)*next
                !!Save EA
                accD%MBS(:,:,1,j1,j2)=accD%MBS(:,:,1,j1,j2)+S(:,:)
                !!Save EB
                accD%MBS(:,:,2,j1,j2)=accD%MBS(:,:,2,j1,j2)+S(:,:)
                !!Save EA+EB
                accD%MBS(:,:,3,j1,j2)=accD%MBS(:,:,3,j1,j2)+S(:,:)
            endif
        enddo
    enddo

end subroutine


!/////////////////////////////////////
!!Computes coherent backscattering
subroutine coherent_backscattering(X,K,K0,I_i,I_cur,accD)
    real(kind=rk), intent(in) :: I_i(4)
    real(kind=rk), intent(in) :: I_cur(4)       !!current Stoke's vector
    real(kind=rk), intent(in) :: X(3)           !!current location of the ray
    real(kind=rk), intent(in) :: K(3),K0(3)     !!incident direction
    type(accessibleData), intent(inout) :: accD !!collected data
    real(kind=rk) :: taua,taub,dx0

    integer :: j1,j2
    real(kind=rk) :: KF(3)
    real(kind=rk) :: phd,X2(3),dx1,dx2,XF1(3),XF2(3)
    real(kind=rk) :: pos1(3),pos2(3),S1(4,4),S2(4,4),S3(4,4)
    complex(kind=rk) :: EA(3),EB(3),EAB(3)


    !!Asserts, turn them off if not used    
    if(debugging) then
        ASSERTC(I_cur(1),>=,0.0_rk)
        ASSERTC(dot_product(K0,K0),<=,(1.0_rk+a_diff))
        ASSERTC(dot_product(K0,K0),>=,(1.0_rk-a_diff))
        ASSERTC(dot_product(K,K),<=,(1.0_rk+a_diff))
        ASSERTC(dot_product(K,K),>=,(1.0_rk-a_diff))
    endif

    call prepare_scatterer_4_cb(I_cur)

    if(g_jsca==1) then
        call First_order_scattering(I_cur,X,K,accD)
        return
    endif

    taua=distance_to_surface(XPATH(:,g_jsca),-K0)
    taub=distance_to_surface(XPATH(:,1),-K0)

    call set_EA_EB(I_cur,I_i,mu0,nu0,sphi0,cphi0,taua,taub)

    do j1=1,ntheb
        do j2=1,nphib
            
            KF(:) = get_direction_vector_cb(j1,j2)
            
            phd = waven*ell*dot_product((XPATH(:,g_jsca)-XPATH(:,1)),(K0+KF))
            
            dx0 = taua-taub
            dx1 = distance_to_surface(XPATH(:,g_jsca),KF)
            dx2 = distance_to_surface(XPATH(:,1),KF)

            call cb_step1(KF)
            
            if(dx1<=threshold) then
                EA(:) = get_EA_field(j1,j2)
                EA(:)=EA(:)*exp(-0.5*abs(dx1))*exp(-cmplx(0.0_rk,phd,kind=rk))
            else
                EA(:)=cmplx(0.0_rk,0.0_rk,kind=rk)
            endif

            if(dx2<=threshold) then
                EB(:) = get_EB_field(j1,j2)
                EB(:)=EB(:)*exp(-0.5*(dx2+dx0))
            else
                EB(:)=cmplx(0.0_rk,0.0_rk,kind=rk)
            endif
            
            EAB=EA+EB
            
            call cb_step2(KF,j2)
        
            call stokes_from_E(S1,EA)
            call stokes_from_E(S2,EB)
            call stokes_from_EAB(S3,S1(1,1),S2(1,1),EAB)
            accD%MBS(:,:,1,j1,j2)=accD%MBS(:,:,1,j1,j2)+S1(:,:)
            accD%MBS(:,:,2,j1,j2)=accD%MBS(:,:,2,j1,j2)+S2(:,:)
            accD%MBS(:,:,3,j1,j2)=accD%MBS(:,:,3,j1,j2)+S3(:,:)
        enddo
        
    enddo

end subroutine


end module
