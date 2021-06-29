!!dataBlock: Holds variables
!!Copyright (C) 2018 Timo Väisänen and University of Helsinki
!!All rights reserved.
!!The new BSD License is applied to this software, see LICENSE.txt
module typedefinitions

    use constants
    implicit none

    type :: dataBlock 
    !!abbreviated int the code as dB
    !!Keeps track of input and precalculated tables.
    !!The policy was to create a type
    !!which can transfer variables easily between
    !!subroutines.


!#INPUT_GENERATOR_SCRIPT1
!#MPITOOLS_SCRIPT1
        real(kind=rk) :: wavel                           !!wavelength
        real(kind=rk) :: ell                             !!incoherent mean free path
        real(kind=rk) :: Fstop                           !!cutoff intensity
        real(kind=rk) :: the0                            !!incident theta angle
        real(kind=rk) :: phi0                            !!incident phi angle
        real(kind=rk) :: phiout                          !!phi output angle
        real(kind=rk) :: tau_c                           !!threshold for optical depth
        real(kind=rk) :: hr                              !!thickness/radius
        real(kind=rk) :: volrad                          !!volume element radius
        real(kind=rk) :: runtime                         !!default run time in hours
        integer :: seed                                  !!seed
        integer :: ntot                                  !!number of rays
        integer :: nsca                                  !!max num of scatt. processes
        integer :: nrot                                  !!asd
        logical :: addInfo                               !!add additional info to output
        logical :: estimateTime                          !!estimate time
        logical :: forced_vm                             !!force the volume element to stay inside the medium
        logical :: surface_att                           !!attenuation starts from the surface
        logical :: overlap_prob                          !!if
        logical :: albedo_override                       !!override albedo to be 1
        integer :: ntheb                                 !!number of theta angles(cb-solution)
        integer :: nphib                                 !!number of phi angles(cb-solution)
        integer :: nthe                                  !!number of theta angles(rt-solution)
        integer :: nphi                                  !!number of theta angles(rt-solution)
        real(kind=rk),allocatable :: CTHEIF(:)           !!cos(theta)
        real(kind=rk),allocatable :: STHEIF(:)           !!sin(theta)
        real(kind=rk),allocatable :: INORM(:)            !!gauss-legendre weighting values used in integration
        real(kind=rk),allocatable :: STHEI(:)            !!cos(theta)
        real(kind=rk),allocatable :: CTHEI(:)            !!sin(theta)
        real(kind=rk),allocatable :: CPHII(:)            !!cos(phi)
        real(kind=rk),allocatable :: SPHII(:)            !!sin(phi)
        real(kind=rk),allocatable :: THEBF(:)            !!theta
        real(kind=rk),allocatable :: CTHEBF(:)           !!cos(theta)
        real(kind=rk),allocatable :: STHEBF(:)           !!sin(theta)
        real(kind=rk),allocatable :: STHEB(:)            !!cos(theta)
        real(kind=rk),allocatable :: CTHEB(:)            !!sin(theta)
        real(kind=rk),allocatable :: CPHIB(:)            !!cos(phi)
        real(kind=rk),allocatable :: SPHIB(:)            !!sin(phi)
        real(kind=rk),allocatable :: PHIB(:)             !!phi
        real(kind=rk),allocatable :: theGrid_cb(:)       !!theta grid points
        real(kind=rk) :: ffboundary                      !!distance to far field
        logical :: polQp                                 !!compute RT with polarization state Q+
        logical :: polQm                                 !!compute RT with polarization state Q-
        logical :: polUp                                 !!compute RT with polarization state U+
        logical :: polUm                                 !!compute RT with polarization state U-
        logical :: polVp                                 !!compute RT with polarization state V+
        logical :: polVm                                 !!compute RT with polarization state V-
        logical :: cb_on                                 !!Enable/Disable coherent backscattering
!#MPITOOLS_SCRIPT2
        integer,allocatable :: run_table(:)              !!list of polarization states
        character(len=64) :: output_details              !!Print input/output parameters
        character(len=64) :: output_rt                   !!RT-solution output filename
        character(len=64) :: output_cb                   !!CB-solution output filename
        character(len=64) :: outputExt                   !!extra information to output
        character(len=64) :: cb_the_points               !!read theta grid points from this file
        character(len=64) :: t_matrix_location           !!The file which contains incoherent T-matrice
!#INPUT_GENERATOR_SCRIPT2
    end type


    type :: assembledData
        !!The store for final results.
        !!When a thread has finished a job with
        !!one Stokes parameter configuration, it updates 
        !!these values using the gathered data(accessibleData).
        !!assembledData has to be initialised before
        !!starting the simulation.
        !!abbreviated as aD
        !!MBS(1:2,1:4,1:4,ntheb,nphib): Mueller cb-enhanced
        !!MRT(1:4,1:4,nthe,nphi): Mueller rt-only
        integer, allocatable :: pol_state(:)            !!polarization states
        real(kind=rk), allocatable :: Aref(:)                  !!reflected/scattered
        real(kind=rk), allocatable :: Atra(:)                  !!transmitted
        real(kind=rk), allocatable :: Adt(:)                   !!direct transmission
        real(kind=rk), allocatable :: Aabs(:)                  !!absorption
        real(kind=rk), allocatable :: Astop(:)                 !!leftover intensity
        real(kind=rk), allocatable :: Aspec(:)                 !!specular
        integer,       allocatable :: rays(:)
        !!4x4x 3(EA,EB,EAB) x ntheb x nphib x pol_states
        real(kind=rk), allocatable :: MBS(:,:,:,:,:,:)  !!CB-solution. Accumulated Mueller matrices
        !!4x4 x nthe x nphi
        real(kind=rk), allocatable :: MRT(:,:,:,:,:)    !!RT-solution. Accumulated Mueller matrices
    end type


    type :: accessibleData
        !!This type holds variables which are
        !!updated frequently during the simulation.
        !!abbreviated as accD
        !!IBS(1:4,1:2,ntheb,nphib): Stokes cb-enhanced
        !!IRT(1:4,nthe,nphi): Mueller rt-only
        real(kind=rk) :: Aref = 0.0_rk                  !!reflected/scattered
        real(kind=rk) :: Atra = 0.0_rk                  !!transmitted
        real(kind=rk) :: Adt = 0.0_rk                   !!direct transmission
        real(kind=rk) :: Aabs = 0.0_rk                  !!absorption
        real(kind=rk) :: Astop = 0.0_rk                 !!leftover intensity
        real(kind=rk) :: Aspec = 0.0_rk                 !!specular
        integer       :: rays = 0
        !!4x4 x 3(EA,EB,EAB) x ntheb x nphib
        real(kind=rk), allocatable :: MBS(:,:,:,:,:)      !!CB-solution. Accumulated intensities
        !!4x4 x nthe x nphi
        real(kind=rk), allocatable :: MRT(:,:,:,:)       !!RT-solution. Accumulated intensities
    end type



end module


