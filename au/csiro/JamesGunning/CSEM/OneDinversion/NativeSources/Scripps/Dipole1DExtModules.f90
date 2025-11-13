!==============================================================================!
!================================================================= dipole1d_io !
!==============================================================================!        
    module dipole1d_io
!   
! This is the input/output (io) module for the dipole1D code.
! External calling routines need to allocate and define the parameters 
! listed below. Electric and magnetic field values are passed out using 
! this module.  
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
! 
    implicit none
!
! Transmitter paramters
!
    real(8)    :: sdm1D             ! (Am)      Source dipole moment 
    real(8)    :: ftx1D             ! (Hz)      Source frequency 
    real(8)    :: azimuthTx1D       ! (degrees) Source azimuth from x axis (positive clockwise)
    real(8)    :: dipTx1D           ! (degrees) Vertical dip angle of source along azimuthTx1D, positive down 
    real(8)    :: xTx1D,yTx1D,zTx1D ! (m)   Transmitter x, y and z location   
!
! Model parameters
!
    integer                             :: nlay1D    ! Number of layers
    real(8), dimension(:), allocatable  :: sig1D     ! (S/m) Layer conductivities 
    real(8), dimension(:), allocatable  :: zlay1D    ! (m)   Depth to top of each layer, first layer ignored 

!
!
! Site locations (where E and B are computed)
!
    integer                             :: n1D          ! number of receiver sites
    real(8), dimension(:), allocatable  :: x1D, y1D, z1D   
!
! Option to only compute electric fields.  lbcomp = .false. or .true. 
!
    logical :: lbcomp  
!
! Fast Hankel Transform Filter method.  Choose the filter coefficients:
!
    character(32)    :: HTmethod1D   ! Filter options listed in FilterModules.f90

!
! Output domain (either spatial (x,y,z) or (kx,y,z) wavenumber domain for 2.5D modeling primary fields)
!
    character(32)    :: outputdomain1D     ! Use 'spatial' or 'kx'
!
! If using outputdomain1D='kx', then set kx1d, kxmode1D and CTmethod1D:
!
! kxmode1D is used to define the principle axis for the Tx
! 
    real(8)          :: kx1D     ! If outputdomain1D='kx', then the fields are 
                                 ! computed at at this wavenumber.  
                              
    integer          :: kxmode1D ! 1,2,3 for x,y,z dipoles.  Ignored if spatial computation.
!
! (Co)Sine transform using digital filter method. 
!
    character(32)    :: CTmethod1D  ! Filter options listed in FilterModules.f90

!
! Use Spline interpolation of coefficients for massive speed up when there are many receivers
!
    logical          :: lUseSpline1D 
    
!
!  Electric and magnetic field vectors.  Note that I output vertical current:
!
    complex(8), dimension(:), allocatable :: ex1D,ey1D,jz1D,bx1D,by1D,bz1D  ! units are V/m and T, e(-iwt) convention used

!
! Invesion sensitivity computations:
!

    logical :: linversion  ! set =.true. to output the field derivatives with respect to conductivity
!
!  Electric and magnetic field derivatives with respect to layer conductivities
!  This are computed when linversion = .true.
! 
    complex(8),dimension(:,:), allocatable :: dexdsig,deydsig,djzdsig,dbxdsig,dbydsig,dbzdsig   ! df_i / dsig_j array
!
! This module contains an initialization set of defaults:
!
    contains 
    
    subroutine init_defaults_Dipole1D 
    
!
! Specify some defaults for parameters required by Dipole1D:
!
    HTmethod1D      = 'kk_ht_201'    ! Use 201 point HT digital filters.
    outputdomain1D  = 'spatial'      ! Assume spatial domain comps
    lbcomp          = .false.        ! This is changed to true if magnetics in data file
    sdm1D           = 1.0            ! (Am), dipole moment. Normalize to unit source moment
    lUseSpline1D    = .true.         ! JG HACK SHOULD BE TRUE Use spline interpolation for faster 1D computations
    linversion      = .false.        ! Compute derivatives with respect to sigma(layers)
    
    end subroutine init_defaults_Dipole1D 
    
    
  end module dipole1d_io
