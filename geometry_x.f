
!---------------------------------------------------------
! Generates volume elements made of spheres specified 
! by the RadiusGenerator module. 
! SHORT DESCRIPTION:
! Pack a large container periodically with scatterers and
! extract spherical volume element from it. The packing is accelerated
! by dividing the container into cells which are used to check
! whether a new scatterer will overlap with existing
! scatterers.
!
! LONG DESCRIPTION
! #1 Find the maximum amount of scatterers in a single cell
! #2 Generate a grid of cells. 
!       (1st column (pos,radius), 2nd column (nth sphere), 3-5th columns (pos in grid))
! #3 Pack periodically the grid with scatterers until the
!       given density of the container is the same as the given density
! #4 Extract a volume element from a container
! 
! NOTICE:
! The packed container can be used again by adjusting 
! N_GEN_NEW0, which tells how many times the same container 
! can be used to extract volume elements until a new container
! will be generated. (Default: N_GEN_NEW0=1, always generate a new)
! The average number of scatterers inside the container: 
! (Default: N_AVG = 10**5, the smaller, the faster)
! The algorithm adjusts the maximum number of scatterers
! by trying how many fits inside a cell. This trial is not
! enough for each cell and that is why the number of scatterers
! must be multiplied (Default: NSPH_MULTIPLY=4).
!
! MAX_FAIL_COUNT (10000), number of tries before the packer
! starts to signal that the container is quite crowded. Won't
! kill the program. 
!
! Copyright (C) 2018 Timo Väisänen, Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! List of public functions:
! subroutine pack_spheres(volrad, dens, eps, N, coords, radiuses, params)
!
! Modifications:
! 23.3.2018 TV: Maintenance, cleaning, commenting
!
! TODO:
! *The algorithm kills the whole program if the 
! number of spheres inside the cell is exceeded.
! A better scheme would be to restart the program.
! Easily done.
! *More explaining
! *Make it possible to force the whole scatterer
! stay inside the volume element
!-----------------------------
module Geometry_x
    use Rng
    use Common
    use RadiusGenerator
    implicit none

    real(kind=dp), parameter :: NSPH_MULTIPLY = 3
    real(kind=dp), parameter :: p_fix = 4*epsilon(p_fix)
    integer, parameter :: N_AVG = 10**6
    integer, parameter :: N_GEN_NEW = 1
    integer, parameter :: MAX_FAIL_COUNT = 10000

    real(kind=dp) :: l_edge
    real(kind=dp) :: l_box,dens
    integer :: n_edge
    integer :: extract_counter
    integer :: n_max 
    logical :: initialized = .false.   

    private 
    public :: pack_spheres

contains


!!Generates a volume element centered to 0. 
subroutine pack_spheres(volrad, dens, eps,i1 , N, coords, radiuses, params)
    real(kind=dp), intent(in)       :: volrad       !!Volume element (VE) radius
    real(kind=dp), intent(in)       :: dens         !!density
    complex(kind=dp), intent(in)    :: eps          !!permittivity of the scatterers
    integer, intent(in)                    :: i1            !!Number of ensemble
    integer, intent(out)                    :: N            !!Number of particles inside the VE
    real(kind=dp), pointer, intent(out)     :: coords(:,:)  !!Coordinates of the scatterer
    real(kind=dp), pointer, intent(out)     :: radiuses(:)  !!Radiuses of the scatterers
    complex(kind=dp), pointer, intent(out)  :: params(:)    !!permittivies of the scatterers



	real(kind=dp) , dimension(:) ,allocatable :: params1_temp, params2_temp
	integer:: i2,n1 ,io
	character(5) x1
	character(len=8) :: fmt ! format descriptor

	fmt = '(I0.0)' ! an integer of width 5 with zeros at the left
	write (x1,fmt) i1 ! converting integer to string using a 'internal file'

	open(2, file = './spheres/refet_'//trim(x1)//'.dat', status = 'old', action = 'read')

	n1 = 0 
	DO
		READ(2,*,iostat=io)
		IF (io/=0) EXIT
		n1 = n1 + 1
	END DO

	allocate( coords(3,n1), radiuses(n1), params(n1), params1_temp(n1) , params2_temp(n1))
	rewind(2)
	DO i2 =1,n1
		READ(2,*) coords(1,i2),coords(2,i2), coords(3,i2), radiuses(i2), params1_temp(i2) , params2_temp(i2) 
		params(i2)=complex(params1_temp(i2),params2_temp(i2))
		!!params(i2)=complex(2.25,0.0)
	END DO

	
	!!coords(:,:)  =coords_temp
	!!params(:) = params_temp
	!!radiuses(:) =radiuses_temp

    N=size(params)
	write(*,*) i1

	
end subroutine

end module
