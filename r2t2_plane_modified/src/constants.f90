#ifndef DEBUG
#define    DEBUGGING .false.
#else
#define    DEBUGGING .true.
#endif

!-------------------------------------------------------------------
! Copyright (C) 2018 Timo Väisänen, University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt
!
! Modifications:
! 4.4.2018 TV: Cleaning
!-----------------------------
module constants
    use, intrinsic :: iso_fortran_env, only : REAL64,REAL128
    use common, only : pi
    implicit none
    integer, parameter :: rk = REAL64                       !!global precision(double)
    integer, parameter :: qp = REAL128                      !!global precision(quad)
    !real(kind=rk), parameter :: pi = 4.0_rk*atan(1.0_rk)    !!pi
    logical, parameter :: debugging = DEBUGGING             !!debugging
    character(*), parameter :: str_line ="//////////////////////////////////////////////////////////////////"
    !!String line in output
end module

