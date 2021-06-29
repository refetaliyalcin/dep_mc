!-------------------------------------------------------------------
! Compute interaction between volume elements
!
! Copyright (C) 2018 Johannes Markkanen, Timo Väisänen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 4.4.2018 TV: Cleaning
!-----------------------------
module T_matrix
    use iso_fortran_env, only : real64 
    use common, only : dp
    use translations, only : translate
    use translations_extension
    implicit none 
contains


!********************************************************************
!* Computes incoherent scattering coefficients oa a cluster of spheres
!*
!* If volume elements overlap (cp1-cp2<ve_rad1+ve_rad2), the
!* translation distance will be ve_rad1+ve_rad2 instead of |cp1-cp2|.
!*  If(first_scatt) the translation check will be omitted
!*
!* The T-matrix information will be fetched from Txx(:,:,cluster_num) 
!*
!* ALL UNITS NEED TO BE IN PHYSICAL UNITS (cp1,cp2,ve_rad1,ve_rad2)
!* and the wavelength in k must be in the same units
!*
!* TODO!
!* This is heavily used and computationally heavy, if I remember correctly. 
!* Significant boost to performance will be appreciated  
!************************************************************************
subroutine scatter_new(cluster_num, k, cp1, cp2, a_in, b_in, ac_nm, bc_nm,&
     Taa, Tab, Tba, Tbb, Nmax, first_scatt,ve_rad1,ve_rad2,rotation,rotD)
    integer, intent(in) :: cluster_num                              !!The index of T-matrix 
    real(dp), intent(in) :: cp1(3), cp2(3)                          !!VE position from (cp1) to to (cp2)
    real(dp), intent(in) :: ve_rad1                                 !!volume element radius (1)
    real(dp), intent(in) :: ve_rad2                                 !!volume element radius (2)
    logical, intent(in) :: first_scatt                              !!Is this the first scattering event
    integer, intent(in) :: Nmax                                     !!Nmax 
    complex(dp), intent(in) :: k                                    !!wavenumber
    complex(dp), dimension(:,:,:),intent(in) :: Taa, Tab, Tba, Tbb  !!T-matrices
    complex(dp), dimension(:), intent(in) :: a_in, b_in             !!Field in
    complex(dp), dimension(:), intent(out) :: ac_nm, bc_nm          !!Field out
    complex(dp), dimension((Nmax + 1)**2-1) ::  a_in2, b_in2        !!temporary storage
    real(dp) :: new_loc(3)                                          !!Where to translate
    complex(dp), intent(inout) :: rotD(:,:)
    real(dp), intent(in) :: rotation(3)


    !real(kind=rk) :: phi,theta,psi
    !integer :: i
    
    
    


    !do i=1,size(Taa,dim=3)
    !    call generate_random_euler(phi,theta,psi)
    !    rotation_store(:,i) = (/phi,theta,psi/)
        !call sph_rotation_gen2([phi,theta,psi], Nmax, rotD)
        !Taa(:,:,i)= matmul(transpose(conjg(rotD)),matmul(Taa(:,:,i),rotD))
        !Tab(:,:,i)= matmul(transpose(conjg(rotD)),matmul(Tab(:,:,i),rotD))
        !Tba(:,:,i)= matmul(transpose(conjg(rotD)),matmul(Tba(:,:,i),rotD))
        !Tbb(:,:,i)= matmul(transpose(conjg(rotD)),matmul(Tbb(:,:,i),rotD))
    !enddo 
    !deallocate(rotD)



    if(dot_product(cp2-cp1,cp2-cp1)<(ve_rad1+ve_rad2)**2 .and. .not. first_scatt) then
        new_loc = (cp2-cp1)/sqrt(dot_product(cp2-cp1,cp2-cp1))*(ve_rad1+ve_rad2)
    else
        new_loc = cp2-cp1
    endif

    !! If this is the first scattering, do not translate fields
    if(first_scatt) then
        a_in2 = a_in
        b_in2 = b_in
    else 
        ! H(cp1) --> J(cp2)
        call translate(new_loc, Nmax, Nmax, k, a_in, b_in, a_in2, b_in2, 1)
    end if

    
    call sph_rotation_gen2(rotation, Nmax, rotD)

    ac_nm = matmul(transpose(conjg(rotD)),matmul(Taa(:,:,cluster_num),matmul(rotD,a_in2))) + &
       &    matmul(transpose(conjg(rotD)),matmul(Tab(:,:,cluster_num),matmul(rotD,b_in2))) 
    bc_nm = matmul(transpose(conjg(rotD)),matmul(Tbb(:,:,cluster_num),matmul(rotD,b_in2))) + &
       &    matmul(transpose(conjg(rotD)),matmul(Tba(:,:,cluster_num),matmul(rotD,a_in2))) 

end subroutine scatter_new


end module T_matrix
