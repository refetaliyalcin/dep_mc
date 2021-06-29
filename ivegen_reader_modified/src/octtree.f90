
!-------------------------------------------------------------------
! Copyright (C) 2018 Johannes Markkanen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 26.3.2018 TV: Cleaning
!-----------------------------

module octtree
use common
implicit none


type octtree_struct
    integer :: near_neighbours(27) 
    integer, dimension(:), allocatable :: sources
    real(dp):: cp(3), dl
    integer :: N_source, Nmax
    integer :: level, parent, children(8)
end type octtree_struct

type level_struct
   type (octtree_struct), dimension(:), allocatable :: tree
end type level_struct


contains


!*************************************************************
subroutine sources_in_cube(cubes, parent, sphere, cube_ind)
    type (octtree_struct), dimension(:), intent(inout) :: cubes
    type (octtree_struct), dimension(:), intent(in) :: parent
    type (data_struct), dimension(:), intent(in) :: sphere
    integer, intent(in) :: cube_ind

    integer :: i1, tri, parent_ind, las
    real(dp):: cp(3)
    real(dp):: xmax, xmin, ymax, ymin, zmax, zmin

    parent_ind = cubes(cube_ind)%parent
    las = 0

    do i1 = 1, size(parent(parent_ind)%sources)
        tri = parent(parent_ind)%sources(i1)
        cp = sphere(tri)%cp

        xmax = cubes(cube_ind)%cp(1) + cubes(cube_ind)%dl / 2.0
        xmin = cubes(cube_ind)%cp(1) - cubes(cube_ind)%dl / 2.0

        ymax = cubes(cube_ind)%cp(2) + cubes(cube_ind)%dl / 2.0
        ymin = cubes(cube_ind)%cp(2) - cubes(cube_ind)%dl / 2.0

        zmax = cubes(cube_ind)%cp(3) + cubes(cube_ind)%dl / 2.0
        zmin = cubes(cube_ind)%cp(3) - cubes(cube_ind)%dl / 2.0

        if(cp(1) <= xmax .and. cp(1) > xmin .and. &
                cp(2) <= ymax .and. cp(2) > ymin .and. &
                cp(3) <= zmax .and. cp(3) > zmin) then
            las = las + 1
        end if
    end do

    cubes(cube_ind)%N_source = las

    if(las>0) then
        allocate(cubes(cube_ind)%sources(las))
    else 
        allocate(cubes(cube_ind)%sources(1))   
        cubes(cube_ind)%sources(1) = 0
    end if

    las = 0

    do i1 = 1, size(parent(parent_ind)%sources)
        tri = parent(parent_ind)%sources(i1)
        
        cp = sphere(tri)%cp

        xmax = cubes(cube_ind)%cp(1) + cubes(cube_ind)%dl / 2.0
        xmin = cubes(cube_ind)%cp(1) - cubes(cube_ind)%dl / 2.0

        ymax = cubes(cube_ind)%cp(2) + cubes(cube_ind)%dl / 2.0
        ymin = cubes(cube_ind)%cp(2) - cubes(cube_ind)%dl / 2.0

        zmax = cubes(cube_ind)%cp(3) + cubes(cube_ind)%dl / 2.0
        zmin = cubes(cube_ind)%cp(3) - cubes(cube_ind)%dl / 2.0

        if(cp(1) <= xmax .and. cp(1) > xmin .and. &
                cp(2) <= ymax .and. cp(2) > ymin .and. &
                cp(3) <= zmax .and. cp(3) > zmin) then
            las = las + 1
            cubes(cube_ind)%sources(las) = tri
        end if
    end do
end subroutine sources_in_cube



!*************************************************************
subroutine divide_cube(cubes, parent, sphere, cube_ind, parent_ind, k)
    type (octtree_struct), dimension(:), intent(out) :: cubes
    type (octtree_struct), dimension(:), intent(in) :: parent
    type (data_struct), dimension(:), intent(in) :: sphere
    real(dp):: cp(3), dl
    real(dp):: cp2(3), dl2, ka
    integer :: Nmax
    real(dp), intent(in) :: k
    integer, intent(in) :: cube_ind,parent_ind


    cp = parent(parent_ind)%cp
    dl = parent(parent_ind)%dl

    dl2 = dl/2

    ka = dl2 * sqrt(3.0d0) / 2.0 * k

    Nmax = truncation_order(ka)

    cp2 = cp + [-dl2/2, -dl2/2, -dl2/2]
    cubes(cube_ind+1)%cp = cp2
    cubes(cube_ind+1)%dl = dl2
    cubes(cube_ind+1)%parent = parent_ind
    cubes(cube_ind+1)%Nmax = Nmax
    call sources_in_cube(cubes, parent, sphere, cube_ind+1)

    cp2 = cp + [dl2/2, -dl2/2, -dl2/2]
    cubes(cube_ind+2)%cp = cp2
    cubes(cube_ind+2)%dl = dl2
    cubes(cube_ind+2)%parent = parent_ind
    cubes(cube_ind+2)%Nmax = Nmax 
    call sources_in_cube(cubes, parent, sphere, cube_ind+2)

    cp2 = cp + [-dl2/2, dl2/2, -dl2/2]
    cubes(cube_ind+3)%cp = cp2
    cubes(cube_ind+3)%dl = dl2
    cubes(cube_ind+3)%parent = parent_ind
    cubes(cube_ind+3)%Nmax = Nmax 
    call sources_in_cube(cubes, parent, sphere, cube_ind+3)

    cp2 = cp + [dl2/2, dl2/2, -dl2/2]
    cubes(cube_ind+4)%cp = cp2
    cubes(cube_ind+4)%dl = dl2
    cubes(cube_ind+4)%parent = parent_ind
    cubes(cube_ind+4)%Nmax = Nmax 
    call sources_in_cube(cubes, parent, sphere, cube_ind+4)

    cp2 = cp + [-dl2/2, -dl2/2, dl2/2]
    cubes(cube_ind+5)%cp = cp2
    cubes(cube_ind+5)%dl = dl2
    cubes(cube_ind+5)%parent = parent_ind
    cubes(cube_ind+5)%Nmax = Nmax 
    call sources_in_cube(cubes, parent, sphere, cube_ind+5)

    cp2 = cp + [dl2/2, -dl2/2, dl2/2]
    cubes(cube_ind+6)%cp = cp2
    cubes(cube_ind+6)%dl = dl2
    cubes(cube_ind+6)%parent = parent_ind
    cubes(cube_ind+6)%Nmax = Nmax 
    call sources_in_cube(cubes, parent, sphere, cube_ind+6)

    cp2 = cp + [-dl2/2, dl2/2, dl2/2]
    cubes(cube_ind+7)%cp = cp2
    cubes(cube_ind+7)%dl = dl2
    cubes(cube_ind+7)%parent = parent_ind
    cubes(cube_ind+7)%Nmax = Nmax 
    call sources_in_cube(cubes, parent, sphere, cube_ind+7)

    cp2 = cp + [dl2/2, dl2/2, dl2/2]
    cubes(cube_ind+8)%cp = cp2
    cubes(cube_ind+8)%dl = dl2
    cubes(cube_ind+8)%parent = parent_ind
    cubes(cube_ind+8)%Nmax = Nmax 
    call sources_in_cube(cubes, parent, sphere, cube_ind+8)

end subroutine divide_cube




!****************************************************************
subroutine create_octtree(sphere, otree, k, max_level)
    type (data_struct), dimension(:), intent(in) :: sphere
    type (level_struct), dimension(:), allocatable, intent(out) :: otree
    type (octtree_struct), dimension(:), allocatable :: tree
    real(dp), intent(in) :: k
    integer, intent(out) :: max_level

    integer, allocatable, dimension(:) :: las_vec
    integer ::  sph, level, i1, i2
    integer :: parent_ind, cube_ind, las, Nmax
    real(dp):: rmax, cp(3), r, dl, cp1(3), cp2(3), a_max,a, ka
    real(dp):: xyz_max(3), xyz_min(3), min_dl

    rmax = 0.0
    a_max = 0.0

    xyz_max(:) = 0.0
    xyz_min(:) = 0.0

    do sph = 1, size(sphere)
        cp = sphere(sph)%cp
        r = sqrt(dot_product(cp,cp)) + sphere(sph)%r
        
        a = sphere(sph)%r
        if(r > rmax) then
            rmax = r
        end if
        if(a > a_max) then
            a_max = a
        end if


        do i1 = 1,3
            if(xyz_max(i1) < cp(i1) + sphere(sph)%r) then
                xyz_max(i1) = cp(i1) + sphere(sph)%r
            end if
            
            if(xyz_min(i1) > cp(i1) - sphere(sph)%r) then
                xyz_min(i1) = cp(i1) - sphere(sph)%r
            end if

        end do

    end do


    !********************************
    !dl = 2*rmax 
    !cp(:) = 0.0d0 
    ! define max level
    !dl = 4.0 * a_max
    !max_level = 0
    !do while (dl < 2*rmax)
    !   dl = dl*2.0d0  
    !   max_level = max_level + 1
    !end do
    !***********************************

    dl = maxval(xyz_max - xyz_min)
    cp = (xyz_min + [dl,dl,dl]/ 2.0)
    min_dl = dl
    max_level = 0
    do while (min_dl > 8.0 * a_max)
        min_dl = min_dl/2.0d0  
        max_level = max_level + 1
    
    end do
    !print*, 'origin:', real(cp)

    !********************************

    if(max_level < 2) max_level = 2
    
    allocate(otree(max_level+1))


    !********** level 0 **********************
    allocate(otree(1)%tree(1))
    otree(1)%tree(1)%cp = cp
    otree(1)%tree(1)%dl = dl
    otree(1)%tree(1)%N_source = size(sphere)

    ka = dl * sqrt(3.0d0) / 2.0 * k
    Nmax = ceiling(2+ka + 4.0*(ka)**(1.0/3.0))  
    otree(1)%tree(1)%Nmax =  Nmax

    allocate(otree(1)%tree(1)%sources(size(sphere)))

    do sph = 1, size(sphere)
        otree(1)%tree(1)%sources(sph) = sph
    end do

    !*****************************************

    do level = 1, max_level
        allocate(tree(8**level))
        cube_ind = 0

        do i1 = 1, size(tree)
            tree(i1)%N_source = 0 
        end do

        do parent_ind = 1, size(otree(level)%tree)
            otree(level)%tree(parent_ind)%children(:) = 0
            call divide_cube(tree, otree(level)%tree, sphere, cube_ind, parent_ind,k)
            cube_ind = cube_ind + 8
        end do

        
        ! remove empty boxes
        las = 0
        do i1 = 1, size(tree)
            if(tree(i1)%N_source > 0) then
                las = las + 1
            end if
        end do

        allocate(otree(level+1)%tree(las))

        las = 0
        do i1 = 1, size(tree)
            if(tree(i1)%N_source > 0) then
                las = las + 1
                otree(level+1)%tree(las) = tree(i1)
                
            end if
        end do

        deallocate(tree)

        !***** find near neighbours ********************
        do i1 = 1, size(otree(level+1)%tree)
            cp1 = otree(level+1)%tree(i1)%cp
            otree(level+1)%tree(i1)%near_neighbours(:) = 0
            las = 0

            do i2 = 1, size(otree(level+1)%tree)
                cp2 = otree(level+1)%tree(i2)%cp

                if(sqrt(dot_product(cp1-cp2,cp1-cp2)) < 1.9*otree(level+1)%tree(i1)%dl) then
                    las = las + 1      
                    otree(level+1)%tree(i1)%near_neighbours(las) = i2
                end if

            end do

        end do
        !*******************************************
        
        end do

        !***** Children cubes **********************!
        do level = 1,max_level 

        allocate(las_vec(size(otree(level)%tree)))
        las_vec(:) = 0
        do i1 = 1, size(otree(level+1)%tree)

            parent_ind = otree(level+1)%tree(i1)%parent ! level
        
            las_vec(parent_ind) = las_vec(parent_ind) + 1 
        
            otree(level)%tree(parent_ind)%children(las_vec(parent_ind)) = i1

        end do
        deallocate(las_vec)
    end do
    !**********************************
end subroutine create_octtree

!**************************************************************

end module octtree
