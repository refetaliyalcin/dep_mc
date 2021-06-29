#include "macros.inc"
!-------------------------------------------------------------------
! Scatterers are presented by incoherent T-matrices
!
! Copyright (C) 2018 Timo Väisänen, Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 4.4.2018 TV: Maintenance, cleaning, commenting
!
! TODO
! *Some refactoring and support for volume elements
! of different sizes
!-----------------------------

module scatterer
    use typedefinitions
    use constants, only : rk
    use rng
    use io
    use grid_holder
    use mpi_tools
    use mie_extension
    use T_matrix
    use math_routines
    implicit none

    integer, parameter      :: ERRORCODE = 1001
    integer                 :: Nmax                 !!truncation order
    integer, allocatable    :: cluster_num_store(:) !!Store index of volume elemnets

    real(kind=rk), allocatable    :: rotation_store(:,:) !!Store rotations



    real(kind=rk)           :: rad_ve               !!volume element radius     [PHYS UNIT]
    real(kind=rk)           :: ell                  !!mean free path
    real(kind=rk), allocatable :: cdf(:),cdf2(:)    !!Cumulative distribution for the probability of the scattering angle
    real(kind=rk), allocatable :: Cexts(:,:),Cexts2(:)!!Storage for extinction coefficients
    real(kind=rk), allocatable :: SMat(:,:,:,:)     !!Store S-matrices between the function calls 
    complex(kind=rk)        :: k                    !!wavenumber

    !!Storage for T-matrices
    complex(kind=rk), allocatable :: Taa(:,:,:), Tab(:,:,:), Tba(:,:,:), Tbb(:,:,:)

    !!Field coefficients
    complex(kind=rk), allocatable :: initWaveA_in(:),initWaveb_in(:)    !!initial waves
    complex(kind=rk), allocatable :: a_out(:),b_out(:),ac_out(:),bc_out(:)
    complex(kind=rk), allocatable :: a_in(:),b_in(:)
    complex(kind=rk), allocatable :: a_back_out(:),b_back_out(:),ac_back_out(:),bc_back_out(:)
    complex(kind=rk), allocatable :: a_back_in(:),b_back_in(:)
    complex(kind=rk), allocatable :: rotD(:,:)

    real(kind=rk) :: normalize_fields_coeff
    integer :: current_polarization,nrot
    
    private
    public :: init_scatterer, scatter0,stokes_from_E,stokes_from_EAB
    public :: get_EA_field,get_EB_field,scatterer_preparations
    public :: get_SMAT,get_SMAT_cb,cb_step1,cb_step2,generate_direction
    public :: set_EA_EB,prepare_scatterer_4_cb,printscattererdata
contains


!//////////////////////////////////////////////////////
!!Initialize T-matrix volume element scatterer
subroutine init_scatterer(dB)
    !!Init scatterer module
    type(dataBlock), intent(in) :: dB
    integer :: id,j1
    logical :: fileExists
    !!asserts
    ASSERTI(dB%nsca>0,"Init scatterer: number of sc4atterers has to be > 0")
    ASSERTI(dB%wavel>0,"Init scatterer: Wavelength has to be >0 ")
    ASSERTI(dB%volrad>=0.0_rk,"Init scatterer: volrad has to be > 0")
    ASSERTI(dB%ell>0,"Init scatterer: mean free path has to be > 0")

    allocate(cluster_num_store(dB%nsca))     
    allocate(rotation_store(3,dB%nsca))   

    !!Store values
    k = 2*pi/dB%wavel
    nrot = dB%nrot
    ell = dB%ell
    rad_ve = dB%volrad



    !!Read T-matrices
    id = mpi_get_my_id()
    if(id==0) then
        inquire(file=trim(adjustl(db%t_matrix_location)),exist=fileExists)
        if(.not. fileExists) then
            write(6,*) "File",trim(adjustl(db%t_matrix_location)),"does not exist"
            call ABORT0("scatterer r2t2 ")
        endif
        call read_T2(db%t_matrix_location, Taa, Tab, Tba, Tbb, Cexts)
    endif

    !!Broadcast data
    call mpi_tools_bcast_Tmat(Taa,Tab,Tba,Tbb)
    call mpi_tools_bcast_real_2darray(Cexts)
    
    !!force albedo to go 1
    if(dB%albedo_override) then
        Cexts(:,:) = 1.0_rk
    endif

    !!make Cexts2(:) with a cumulative distribution
    !!of the Cexts(1,:) (extinction coefficients)
    allocate(Cexts2(size(Cexts,dim=2)+1))
    Cexts2(1)=0.0_rk
    do j1=1,size(Cexts,dim=2)
        Cexts2(j1+1)=Cexts2(j1)+Cexts(1,j1)
    enddo

    !Compute truncation order and initialize mie_extension and translation
    Nmax = int(sqrt(size(Taa,1)+1.0_rk)) - 1  
    call init_mie(dB%cb_on,dB%ffboundary,Nmax,theGrid,phiGrid,theGrid_cb,phiGrid_cb,k)
    call init_translations(Nmax)

    !!Allocate field coefficietns
    allocate(a_in(size(Taa,1)), b_in(size(Taa,1)))
    allocate(a_out(size(Taa,1)), b_out(size(Taa,1)))
    allocate(ac_out(size(Taa,1)), bc_out(size(Taa,1)))
    allocate(a_back_in(size(Taa,1)), b_back_in(size(Taa,1)))
    allocate(a_back_out(size(Taa,1)), b_back_out(size(Taa,1)))
    allocate(ac_back_out(size(Taa,1)), bc_back_out(size(Taa,1)))
    allocate(cdf(size(phiGrid)*size(theGrid)),SMat(4,4,size(theGrid),size(phiGrid)))
    allocate(initWaveA_in(size(Taa,1)), initWaveB_in(size(Taa,1)))
    allocate(cdf2(size(phiGrid)*size(theGrid)))

    allocate(rotD(size(Taa,dim=1),size(Taa,dim=2)))

    ASSERTC2(real(k),>,0)       
    ASSERTC2(allocated(bc_out),.eqv.,.true.)
    ASSERTC2(allocated(ac_out),.eqv.,.true.)
    ASSERTC2(allocated(a_in),.eqv.,.true.)
    ASSERTC2(allocated(b_in),.eqv.,.true.)
    ASSERTC2(allocated(a_out),.eqv.,.true.)
    ASSERTC2(allocated(b_out),.eqv.,.true.)
    ASSERTC2(size(cluster_num_store),>,0)   
    ASSERTC2(Nmax,>,0)

end subroutine


!//////////////////////////////////////////////////////////
!!Set polarization state and rotate T-matrices
subroutine scatterer_preparations(polarization,flag2,flag3,ray)
    integer, intent(in) :: polarization,flag2,flag3 !!polarization and some flags
    complex(kind=rk), allocatable :: tmpA(:,:), tmpB(:,:)
    integer, intent(in) :: ray
    if(ray==0) then
        allocate(tmpA(size(initWaveA_in),6),tmpB(size(initWaveA_in),6)) 
        call inc_waves(Nmax, real(k), tmpA, tmpB)
        current_polarization=polarization
        initWaveA_in(:) = tmpA(:,polarization)
        initWaveB_in(:) = tmpB(:,polarization)
        deallocate(tmpA,tmpB)
    endif
    !if(nrot>=0 .and. ray==0) then
    !    call rotate_T_matrices()
    !elseif(nrot>0) then
    !    if(mod(ray,nrot)==0) call rotate_T_matrices2()
    !endif
    !call rotate_T_matrices2()
end subroutine


!////////////////////////////////////////////
!!Select T-matrix by using extinction cross sections
subroutine generate_cluster()
    integer :: j1
    real(kind=rk) :: cValue,phi,theta,psi

    cValue = Cexts2(size(Cexts2))*randNum()
    j1 =  bisectionSearch(Cexts2,cValue)-1
    if(j1==0) j1 = 1 
    if(j1>size(Taa,dim=3)) j1 = size(Taa,dim=3)
    cluster_num_store(g_jsca) = j1

    
    call generate_random_euler(phi,theta,psi)
    rotation_store(:,g_jsca) = (/phi,theta,psi/)
end subroutine

!////////////////////////////////////////
!!Fields hits the volume element and scatters
subroutine scatter0(I0, albedo)
    real(kind=rk), intent(out) :: albedo        !!Albedo of the VE

    real(kind=rk) :: cp1(3),cp2(3)
    real(kind=rk), intent(in) :: I0(4)
    real(kind=rk) :: aux_tmp
    complex(kind=rk) :: lambda
    
    !!Genererate volume element
    call generate_cluster()

    if(g_jsca==1) then
        !!This was the first scattering
        cp1(:)=0.0_rk
        cp2(:)=0.0_rk
        !!Scatter using incident planewave
        lambda = cmplx(0.0_rk,real(k)*ell*XPATH(3,1)*1.0_rk,kind=rk)
        call scatter_new(cluster_num_store(g_jsca), k, cp1, cp2, &
            & initWaveA_in*exp(lambda), initWaveB_in*exp(lambda), &
            & ac_out, bc_out, Taa, Tab, Tba, Tbb, Nmax,.true.,rad_ve,rad_ve, rotation_store(:,g_jsca),rotD) 
    else
        cp1(:) = XPATH(:,g_jsca-1)*ell
        cp2(:) = XPATH(:,g_jsca)*ell
  
        !!Scatter using incident planewave
        call scatter_new(cluster_num_store(g_jsca), k, cp1, cp2, a_in, &
            & b_in, ac_out, bc_out, Taa, Tab, Tba, Tbb, Nmax,.false.,rad_ve,rad_ve, rotation_store(:,g_jsca),rotD)
    endif
    
    !!Every 10th scattering we will need to normalize coefficients to prevent underflow
    if(mod(g_jsca,10)==0) then
        aux_tmp=sum(real(ac_out))
        a_in=ac_out/aux_tmp
        b_in=bc_out/aux_tmp
    else
        aux_tmp=1.0_rk
        a_in=ac_out
        b_in=bc_out
    endif        

    !!Generate intensity map
    call solveWholeSpace()

    !!Set albedo
    albedo = Cexts(2,cluster_num_store(g_jsca))
    if(albedo>1.0_rk) then
        albedo=1.0_rk
    elseif(albedo<0.0_rk) then
        albedo=0.0_rk
    endif 
    ASSERTC(normalize_fields_coeff,>,0.0_rk)
end subroutine






!///////////////////////////////////////////////////////////////////////////////////////////
subroutine solveWholeSpace()
    !!Solve scattering angle distribution
    integer :: j1,j2,ind
    real(kind=rk) :: K2(3)
    real(kind=rk) :: S(4,4)
    real(kind=rk) :: total
    complex(kind=rk) :: EA(3)
    real(kind=rk) :: X2(3)
    real(kind=rk) :: tmp2,tmp1,t

    if(debugging) then
        ASSERTC(far_field_boundary,>,1.0_rk)
        ASSERTC(size(theGrid),==,nthe)
        ASSERTC(size(phiGrid),==,nphi)
        ASSERTC(allocated(SMat),.eqv.,.true.)
        ASSERTC(size(cdf),==,nphi*nthe)
        ASSERTC(size(SMat,1),==,4)
        ASSERTC(size(SMat,2),==,4)
        ASSERTC(size(SMat,3),==,nthe)
        ASSERTC(size(SMat,4),==,nphi)
    endif

    ind = 0
    total = 0.0_rk
    t=0
    do j1=1,nthe
        do j2=1,nphi
            ind=ind+1
            K2(1)=sin(theGrid(j1))*cos(phiGrid(j2))
            K2(2)=sin(theGrid(j1))*sin(phiGrid(j2))
            K2(3)=cos(theGrid(j1))

            X2=K2*far_field_boundary
            call get_e_field(EA,ac_out, bc_out,j1,j2)
            call stokes_from_E(S,EA)

            if(j1==1) then
                tmp1=(theGrid(j1)+theGrid(j1+1))*0.5_rk
                tmp2=pi
            elseif(j1==nthe) then
                tmp1=0.0_rk
                tmp2=(theGrid(j1)+theGrid(j1-1))*0.5_rk
            else
                tmp1=(theGrid(j1)+theGrid(j1+1))*0.5_rk
                tmp2=(theGrid(j1)+theGrid(j1-1))*0.5_rk
            endif
            cdf2(ind)=t+(2*pi)/nphi*abs(cos(tmp1)-cos(tmp2))*S(1,1)
            t = cdf2(ind)

            cdf(ind)=total+S(1,1)*INORM(j1)
            SMat(:,:,j1,j2)=S(:,:)
            total=total+S(1,1)*INORM(j1)
        enddo
    enddo

    normalize_fields_coeff = 1.0_rk/cdf(size(cdf))
    if(debugging) then
        ASSERTC(is_cumulative(cdf),.eqv.,.true.)
    endif
end subroutine





!///////////////////////////////////////////////
!!Generate scattering direction from the intensity map
subroutine generate_direction(I1,K1,EH1,EV1)
    real(kind=rk), intent(inout) :: I1(4)    !!incident and then outgoing wave
    real(kind=rk), intent(inout) :: K1(3)    !!incident and then outgoing direction
    real(kind=rk), intent(out) :: EH1(3)    !!scattered electric fields, obsolete
    real(kind=rk), intent(out) :: EV1(3)    !!scattered electric fields, obsolete
    real(kind=rk) :: the,phi
    integer :: retInd
    real(kind=rk) :: rnd
    if(debugging) then
        ASSERTC(size(cdf),==,nphi*nthe)
    endif

    rnd = cdf2(size(cdf2))*randNum()
    call binarySearch(rnd,retInd,cdf2,nthe,nphi)
    call getDirection0(retInd,the,phi)
    K1(1)=sin(the)*cos(phi)
    K1(2)=sin(the)*sin(phi)
    K1(3)=cos(the)

    EH1=0.0_rk
    EV1=0.0_rk

    if(debugging) then
        !write(6,*) the,phi
        ASSERTC(dot_product(K1,K1),<=,1.0_rk+a_diff)
        ASSERTC(dot_product(K1,K1),>=,1.0_rk-a_diff)
    endif
end subroutine




!!----------------------------------------------------
!!Get Stokes generated from electric field E
!!----------------------------------------------------
pure subroutine SFromE(S,E)
    complex(kind=rk),intent(in) :: E(3)     !Electric field
    real(kind=rk),intent(out) :: S(4,4)       !Stokes
    real(kind=rk) :: qh,qv
    complex(kind=rk) :: img
    
    img=cmplx(0.0_rk,1.0_rk,kind=rk)
    
    qh=real(E(2)*conjg(E(2)))
    qv=real(E(3)*conjg(E(3)))
    
    S = 0.0_rk
    S(1,1)=qh+qv
    S(1,2)=qh-qv
    S(1,3)=real(-E(3)*conjg(E(2))+conjg(-E(3))*E(2))
    S(1,4)=real(img*(E(3)*conjg(E(2))-conjg(E(3))*E(2)))
end subroutine



!///////////////////////////////////////////////////////////////////////////////////////////
subroutine set_EA_EB(I_cur,I_i,mu0,nu0,sphi0,cphi0,taua,taub)
    !!Get reversed path EB
    integer :: j1
    real(kind=rk) :: cp1(3),cp2(3)
    logical :: last
    real(kind=rk) :: tmp_norm1,tmp_norm2
    real(kind=rk) :: KF(3),aux_tmp
    real(kind=rk) :: pos1(3)
    real(kind=rk), intent(in) :: I_cur(4),I_i(4),mu0,nu0,sphi0,cphi0,taua,taub

    complex(kind=rk) :: lambda
    complex(kind=rk) :: EA(3),EB(3),EAB(3)
    !!asserts turn them off, if not needed
    if (debugging) then
        ASSERTC(real(k),>,0.0_rk)
        ASSERTC(size(initWaveA_in),==,size(a_in))
        ASSERTC(size(initWaveB_in),==,size(b_in))
        ASSERTC(g_jsca,>,0)
    endif 

    cp1(:)=0.0_rk
    cp2(:)=0.0_rk

    lambda = cmplx(0.0_rk,real(k)*ell*XPATH(3,1)*1.0_rk,kind=rk)

    !!Reversed path
    call scatter_new(cluster_num_store(g_jsca), k, cp1, cp2,  &
        & initWaveA_in*exp(lambda),initWaveB_in*exp(lambda),  &
        & ac_back_out, bc_back_out,Taa, Tab, Tba, Tbb, Nmax,  &
        & .true.,rad_ve,rad_ve, rotation_store(:,g_jsca),rotD) 

    a_back_out=ac_back_out
    b_back_out=bc_back_out

    last = .false.
    do j1=g_jsca-1,1,-1
        a_back_in = a_back_out
        b_back_in = b_back_out
        cp1 = XPATH(:,j1+1)*ell
        cp2 = XPATH(:,j1)*ell
       ! write(6,*) cp1,cp2,sqrt(dot_product(cp1-cp2,cp1-cp2))
        if(j1==1) last = .true.
        call scatter_new(cluster_num_store(j1), k, cp1, cp2,  &
            & a_back_in, b_back_in, ac_back_out, bc_back_out, &
            & Taa, Tab, Tba, Tbb, Nmax,.false.,rad_ve,rad_ve, rotation_store(:,j1),rotD) !!!!!!
        if(mod(j1,10)==0) then
            aux_tmp=sum(real(ac_back_out))
        else
            aux_tmp=1.0_rk
        endif    
        a_back_out  = ac_back_out/aux_tmp
        b_back_out  = bc_back_out/aux_tmp
    enddo

    !!reciprocity normalization. Notice that the backscattering angle is shifted a little
    KF(1)=0.0_rk
    KF(2)=1.00_rk
    KF(3)=far_field_boundary
    call getSphericalCoordinates(KF,pos1(1),pos1(2),pos1(3))
    pos1(2)=acos(pos1(2))
  

    call get_e_field2(EA,pos1,ac_out, bc_out,  k)
    call get_e_field2(EB,pos1,ac_back_out, bc_back_out, k)

    tmp_norm1 = real(EA(2)*conjg(EA(2))+EA(3)*conjg(EA(3)))
    tmp_norm2 = real(EB(2)*conjg(EB(2))+EB(3)*conjg(EB(3)))

    if(tmp_norm2==0.0_rk) then
        tmp_norm2=1.0_rk
    endif
    ac_back_out = sqrt(tmp_norm1)/sqrt(tmp_norm2)*ac_back_out 
    bc_back_out = sqrt(tmp_norm1)/sqrt(tmp_norm2)*bc_back_out

end subroutine




!///////////////////////////////////////////////////////////////////////////////////////////
subroutine  get_SMAT(S,I0,K1,K0,j1,j2)
    !!Get scattering matrix
    real(kind=rk), intent(out) :: S(4,4)
    !!Output scattering matrix
    real(kind=rk), intent(in) :: I0(4)
    !!Stokes vector of the incident field
    real(kind=rk), intent(in) :: K1(3),K0(3)
    !!outgoing and incident direction
    integer, intent(in) :: j1,j2
    !!index location of the detector
    logical :: arr(4,4)
    if (debugging) then
        ASSERTC(size(cdf),==,nphi*nthe)
        ASSERTC(is_cumulative(cdf),.eqv.,.true.)
        ASSERTC(SMat(1,1,j1,j2),>=,0.0_rk)
        ASSERTC(I0(1),>=,0.0_rk)
        ASSERTC(dot_product(K1,K1),<=,1+a_diff)
        ASSERTC(dot_product(K1,K1),>=,1-a_diff)
        ASSERTC(dot_product(K0,K0),<=,1+a_diff)
        ASSERTC(dot_product(K0,K0),>=,1-a_diff)
        ASSERTC(j1,<=,nthe)
        ASSERTC(j1,>,0)
        ASSERTC(j2,<=,nphi) 
        ASSERTC(j2,>,0) 
        ASSERTC(INORM(j1),>,0.0_rk)
    endif 

    if(cdf(size(cdf))/=0) then
        S(:,:)=I0(1)*SMat(:,:,j1,j2)/cdf(size(cdf))*INORM(j1)
    else
        S(:,:) = 0.0_rk
    endif
    if (debugging) then
        arr = .false.
        where(S>HUGE(S(1,1))) arr = .true. !!check infinity
        ASSERTC(count(arr),==,0)
        where(S /= S) arr = .true.  !!check NaN
        ASSERTC(count(arr),==,0)
        ASSERTC(S(1,1),>=,0.0_rk)
    endif 
end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
subroutine  prepare_scatterer_4_cb(I0)
    real(kind=rk), intent(in) :: I0(4)

    !!Get field normalization constant which cb requires
    if (debugging) then
        ASSERTC(normalize_fields_coeff,>,0.0_rk)
    endif 
    ac_out = sqrt(I0(1)*normalize_fields_coeff)*ac_out
    bc_out = sqrt(I0(1)*normalize_fields_coeff)*bc_out
end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
subroutine  get_SMAT_cb(S,I0,K1,K0,j1,j2)
    !!Get scattering matrix
    real(kind=rk), intent(out) :: S(4,4)
    !!Output scattering matrix
    real(kind=rk), intent(in) :: I0(4)
    !!Stokes vector of the incident field
    real(kind=rk), intent(in) :: K1(3),K0(3)
    !!outgoing and incident direction
    integer, intent(in) :: j1,j2
    !!index location of the detector
    real(kind=rk) :: X2(3),S1(4,4)
    real(kind=rk) :: XF1(3)
    complex(kind=rk) :: EA(3)
    logical :: arr(4,4)
    if (debugging) then
        ASSERTC(size(cdf),==,nphi*nthe)
        ASSERTC(is_cumulative(cdf),.eqv.,.true.)
        ASSERTC(SMat(1,1,j1,j2),>=,0.0_rk)
        ASSERTC(I0(1),>=,0.0_rk)
        ASSERTC(dot_product(K1,K1),<=,1+a_diff)
        ASSERTC(dot_product(K1,K1),>=,1-a_diff)
        ASSERTC(dot_product(K0,K0),<=,1+a_diff)
        ASSERTC(dot_product(K0,K0),>=,1-a_diff)
        ASSERTC(j1,<=,ntheb)
        ASSERTC(j1,>,0)
        ASSERTC(j2,<=,nphib) 
        ASSERTC(j2,>,0) 
    endif

    !!Get detector point in the far field region
    !X2(:)=K1(:)*far_field_boundary
    !call cartesianToSpherical(X2,pos(1),pos(2),pos(3))
    !pos(2) = acos(pos(2))

    !!Compute electric field there
    call get_e_field_cb(EA,ac_out, bc_out,j1,j2)
    call get_S_from_E(EA,real(k),S1,0,0.0_rk,0.0_rk,X2,XF1,far_field_boundary,current_polarization)
    !Etmp(1)=EA(2)
    !Etmp(2)=EA(3)
    !call stokesFromE(II,Etmp)
    !S(1,1) = II(1)
    !S(:,:)=S1(:,:)/cdf(size(cdf))
    S(:,:) = S1(:,:)
    if (debugging) then
        arr = .false.
        where(S>HUGE(S(1,1))) arr = .true. !!check infinity
        ASSERTC(count(arr),==,0)
        where(S /= S) arr = .true.  !!check NaN
        ASSERTC(count(arr),==,0)
        ASSERTC(S(1,1),>=,0.0_rk)
    endif 
end subroutine

subroutine cb_step1(KF)
    real(kind=rk), intent(in) :: KF(3)
end subroutine

subroutine cb_step2(KF,j2)
    real(kind=rk), intent(in) :: KF(3)
    integer, intent(in) :: j2
end subroutine




!///////////////////////////////////////////////////////////////////////////////////////////
function get_EA_field(j1,j2) result(EA)
    !!get EA
    integer, intent(in) :: j1,j2            !!location in grid holder
    complex(kind=rk) :: EA(3)               !!EA field
    if (debugging) then
        ASSERTC(j1,<=,ntheb)
        ASSERTC(j1,>,0)
        ASSERTC(j2,<=,nphib) 
        ASSERTC(j2,>,0) 
    endif
    call get_e_field_cb(EA,ac_out, bc_out,j1,j2)
end function


!///////////////////////////////////////////////////////////////////////////////////////////
function get_EB_field(j1,j2) result(EB)
    !!Get EB
    integer, intent(in) :: j1,j2            !!location in grid holder
    complex(kind=rk) :: EB(3)               !!EA field
    if (debugging) then
        ASSERTC(j1,<=,ntheb)
        ASSERTC(j1,>,0)
        ASSERTC(j2,<=,nphib) 
        ASSERTC(j2,>,0) 
    endif
    call get_e_field_cb(EB,ac_back_out, bc_back_out,j1,j2)
end function

!///////////////////////////////////////////////////////////////////////////////////////////
subroutine stokes_from_E(S1,EA)
    real(kind=rk), intent(out) :: S1(4,4)
    complex(kind=rk), intent(in) :: EA(3)
    call SfromE(S1,EA)
end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
subroutine stokes_from_EAB(S1,IAF,IBF,EAB)
    !!Generate 2x2 scattering matrix from E
    real(kind=rk), intent(in) :: IAF,IBF
    complex(kind=rk), intent(in) :: EAB(3)
    real(kind=rk), intent(out) :: S1(4,4)
    call stokes_from_E(S1,EAB)
    if(IAF/=0.0_rk) then
        S1(:,:) = S1(:,:)*IAF/(IAF+IBF)
    else
        S1(:,:) = 0.0_rk
    endif
end subroutine



!///////////////////////////////////////////////////////
!!Generate 2x2 scattering matrix from E
subroutine get_S_from_E(F,k,S,polarization,the,phi,K1,K2,R,flag)
    
    real(kind=rk), intent(in) :: R
    real(kind=rk), intent(in) :: the,phi
    real(kind=rk), intent(in) :: k
    real(kind=rk), intent(in) :: K1(3),K2(3)
    real(kind=rk), intent(out) :: S(4,4)
    complex(kind=rk), intent(in) :: F(3)
    integer, intent(in) :: flag
    integer, intent(in) :: polarization
    call SfromE(S,F)
end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
subroutine printScattererData(dB,flag,stream)
    !!RT-CB will require this because
    type(dataBlock), intent(in)         :: dB
    integer, intent(in) :: stream
    integer, intent(in)                 :: flag
    write(stream,'(A)') "R2T2, reciprocal transactions"
end subroutine


!!Get direction near the grid point ind.
!!@note
!!\($1+1=1$ \),
!!$How this is handled needs to be done better.
!!$Maybe in the future.
!!$In order to get uniformly distributed points on the sphere, it is enough to 
!!divide phi angles equally but theta angles are distributed according to
!!the = acos(2v-1), where v is uniform random number.
!!First find the grid point which was picked from the 
!!cumulative distribution cdf.
!!After that pick a new random point near that grid point.
!!The new phi angle is phi_new = phi_grid-dPhi*0.5+dPhi*u
!! where u is uniform random number.
!!Theta angle is then picked in the following way
!!Get lower theta angle between detectors the_1
!!Get higher theta angle between detectors the_2
!!Get correct v from the theta's distribution shown above
!!v_1=(cos(the_1)+1)*0.5
!!v_2=(cos(the_2)+1)*0.5
!!Pick from the distribution a new v => v_3 = (v_2-v_1)*rand()
!!and get the corresponding theta angle with
!!the_3 = acos(2(v_1+v_3)-1)
!!tarkista gradusta
subroutine getDirection0(ind,thetaX,phiX)

    integer, intent(in) :: ind
    real(kind=rk), intent(out) :: thetaX,phiX

    integer :: ind1,ind2
    real(kind=rk) :: dPhi,dy,phi0,the0,the1,phi,the,y1,y2,tmp  
    integer :: a,b,c
    volatile c,a,b

    if(debugging) then   
        ASSERTCP(ind,>,0)  
        ASSERTCP(ind,<=,size(theGrid)*size(phiGrid))    
        ASSERTCP(size(theGrid),>,0)   
        ASSERTCP(size(phiGrid),>,0)   
    endif

    !!extract indices
    c=ind
    a = ceiling(c/(1.0_rk*size(phiGrid)))
    b = c-(a-1)*size(phiGrid)
    ind1=a
    ind2=b

    !!Generate Theta and Phi angle
    if(ind1==1) then
        the0=(theGrid(ind1)+theGrid(ind1+1))*0.5_rk
        the1=pi
    elseif(ind1==nthe) then
        the0=0.0_rk
        the1=(theGrid(ind1)+theGrid(ind1-1))*0.5_rk 
    else
        the0=(theGrid(ind1)+theGrid(ind1-1))*0.5_rk
        the1=(theGrid(ind1)+theGrid(ind1+1))*0.5_rk      
    endif 
    dPhi=2.0_rk*pi/size(phiGrid)
    phi0=phiGrid(ind2)-dPhi*0.5_rk
    phi=phi0+dPhi*randNum()
    y1=(cos(the0)+1.0_rk)*0.5_rk
    y2=(cos(the1)+1.0_rk)*0.5_rk
    dy=abs(y2-y1)
    if(y1>y2) then
        the=acos(2.0_rk*(y2+dy*randNum())-1)
    else
        the=acos(2.0_rk*(y1+dy*randNum())-1)
    endif
    thetaX=the
    phiX=phi

    if(debugging) then   
        ASSERTC(b+(a-1)*size(phiGrid),==,ind)
        if(a==nthe) then
            tmp = pi
        else
            tmp = theGrid(a)+(theGrid(a)+theGrid(a+1))*0.5_rk
        endif
        ASSERTC(thetaX,<=,tmp)   
        if(a==1) then
            tmp = -pi
        else
            tmp = theGrid(a)-(theGrid(a)+theGrid(a-1))*0.5_rk
        endif
        ASSERTC(thetaX,>=,tmp)   
        tmp = phiGrid(b)+pi/nphi
        ASSERTC(phiX,<=,tmp)   
        tmp = phiGrid(b)-pi/nphi
        ASSERTC(phiX,>=,tmp)   
    endif
end subroutine


!!---------------------------------------------
!!bisectionSearch(xp,the)
!!Search index for the value(lower) closest to
!!the value "the" using bisection search.
!!
!!IN: xp:    data   
!!   the:    the value to look for
!!---------------------------------------------
pure function bisectionSearch(xp,the) result(retVal)
    real(kind=rk), intent(in) :: the         !point to be evaluated
    real(kind=rk), intent(in) :: xp(:)       !sorted array
    integer :: retVal
    integer :: x1,x2,xmid 
    x1=1
    x2=size(xp)
    if(the<xp(1)) then
        retVal = 1
        return
    endif
    do while(x2-x1>1)
        xmid=(x1+x2)/2
        if(xp(xmid) > the) then
            x2=xmid
        else
            x1=xmid
        endif
    enddo
    retVal=x2
end function



!!Ken Shoemake. In "Graphics Gems IV", pp 175-192. 
!!Generate uniform random euler angles
subroutine generate_random_euler(phi,theta,psi)
    real(kind=rk), intent(out) :: phi,theta,psi
    real(kind=rk) :: s,sigma1,sigma2,theta1,theta2,q0,q1,q2,q3
    s=randNum()
    sigma1=sqrt(1-s)
    sigma2=sqrt(s)
    theta1 = 2*pi*randNum()
    theta2 = 2*pi*randNum()
    q0 = cos(theta2)*sigma2
    q1 = sin(theta1)*sigma1
    q2 = cos(theta1)*sigma1
    q3 = sin(theta2)*sigma2

    phi = atan2(2*(q0*q1+q2*q3),1-2*(q1**2+q2**2))
    theta = asin(2*(q0*q2-q3*q1))
    psi = atan2(2*(q0*q3+q1*q2),1-2*(q2**2+q3**2))
end subroutine


end module 
