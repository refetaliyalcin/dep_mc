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
        real(kind=rk) :: wavel        !wavelength
        real(kind=rk) :: mre          !medium refractive index(real)
        real(kind=rk) :: mim          !medium refractive index(imag)
        real(kind=rk) :: ssalbedo     !single-scattering albedo
        real(kind=rk) :: ell1         !mean free path
        real(kind=rk) :: Fstop        !minimum relative flux
        real(kind=rk) :: the0         !incident theta angle
        real(kind=rk) :: phi0         !incident phi angle
        real(kind=rk) :: phiout       !phi output angle
        real(kind=rk) :: tau_c        !threshold for optical depth
        real(kind=rk) :: hr           !thickness/radius

        integer :: seed               !seed
        integer :: nrn                !number of points in cdf
        integer :: ntot               !number of rays
        integer :: nsca               !max num of scatt. processes
        integer :: ncm                !integration points in spline presen.
        logical :: finite             !medium is finite
        logical :: addInfo            !add additional info to output
        logical :: estimateTime       !estimate time
        logical :: I21test            !before run test that 
                                      !the spline presention is valid.

       
        real(kind=rk) :: mu0                !incident mu-angle
        real(kind=rk) :: nu0                !incident mu-angle
        real(kind=rk) :: cphi0              !incident cos(phi)
        real(kind=rk) :: sphi0              !incident sin(phi)
        
        real(kind=rk) :: tau                !thickness/radius as optical depth
        real(kind=rk) :: albedo             !albedo

        real(kind=rk) :: mre12              !refrel    
        real(kind=rk) :: mre21              !1/refrel
        complex(kind=rk) :: m12             !complex refractive index
        complex(kind=rk) :: m21             !1/complex refrective index

        real(kind=rk) :: waven              !wavenumber
        real(kind=rk) :: xell               !mean free path in a host media

        real(kind=rk) :: dthe               !angle between theta bins
        integer :: ns                       !number of spline
        integer :: ntheb                    !number of theta angles(cb-solution)
        integer :: nphib                    !number of phi angles(cb-solution)
        integer :: nthe                     !number of theta angles(rt-solution)
        integer :: nphi                     !number of theta angles(rt-solution)


        !Detector bins(RT-solution)
        real(kind=rk),allocatable :: CTHEIF(:)  !cos(theta)
        real(kind=rk),allocatable :: STHEIF(:)  !sin(theta)
        real(kind=rk),allocatable :: INORM(:)   !gauss-legendre weighting values used in integration
        real(kind=rk),allocatable :: STHEI(:)   !cos(theta), takes account of host media
        real(kind=rk),allocatable :: CTHEI(:)   !sin(theta), takes account of host media
        real(kind=rk),allocatable :: CPHII(:)   !cos(phi), takes account of host media
        real(kind=rk),allocatable :: SPHII(:)   !sin(phi), takes account of host media


        !Detector bins(CB-solution)
        real(kind=rk),allocatable :: THEBF(:)   !theta
        real(kind=rk),allocatable :: CTHEBF(:)  !cos(theta)
        real(kind=rk),allocatable :: STHEBF(:)  !sin(theta)
        real(kind=rk),allocatable :: STHEB(:)   !cos(theta), takes account of host media
        real(kind=rk),allocatable :: CTHEB(:)   !sin(theta), takes account of host media
        real(kind=rk),allocatable :: CPHIB(:)   !cos(phi), takes account of host media
        real(kind=rk),allocatable :: SPHIB(:)   !sin(phi), takes account of host media
        real(kind=rk),allocatable :: PHIB(:)    !phi


        !Scatterer Specific (Spline)
        real(kind=rk) :: vf                         !volume fraction
        real(kind=rk) :: scRadius                   !scatterer radius
        real(kind=rk) :: scRealRef                  !scatterer real refractive index
        real(kind=rk) :: scImagRef                  !scatterer imag refractive index
        logical :: generateMieScatterer             !generate Mie scattering otherwise load from file
        integer :: np                               !points in spline presentation

!#MPITOOLS_SCRIPT2

        character(len=64) :: scattererData          !load scattering data
        character(len=64) :: saveScatterer          !save scatterer data

        character(len=64) :: output_details !Print input/output parameters
        character(len=64) :: output_rt      !RT-solution output filename
        character(len=64) :: output_cb      !CB-solution output filename
        character(len=64) :: scattererType  !Scatterer type
        character(len=64) :: outputExt      !extra information to output
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


