! 
! Modules for interfacing CSEM1D to java
!
! All allocatables from Kerry Key's original f90 
! that are performed at parse-input
! time are now declared as f90 pointers.
! Java routines do the parsing, allocations, and JNI
! routines set the fortran module members to point to the
! Java VM heap space allocated

! Revision  Date           Author           Notes
!
! 1         Oct 08         James Gunning
! 2         Nov 09         ''               Update to match Scripps release V3.9


!==============================================================================!    
!================================================================! CSEM1D_MOD  ! 
!==============================================================================!

    module csem1d_mod

!    use Occam, only : realprec    ! FIXME JG
     
    implicit none    
!
! Module to store various parameters for 1D CSEM inversion
! This is stuff that is not used in Occam.f90 and not Dipole1D.f90
! but is required for interfacing between them.

!   Magic codes and constants in this include
    include 'CSEM1D_parameters.inc'

! While the data parameter array dp is passed in/out of Occam, it is only used
! in the forward routines.  So we set cnDataParams in the forward routines
! and then allocate the dp array in readData.
!        
    integer :: cnDataParams
 
    ! JG real(8), dimension(:,:), allocatable :: layer_params
    real(8), pointer  :: layer_params(:,:) 
!
! layer_params is [ top_depth, rho, penalty, Preference, prej_penalty ]          
!
! top_depth :: Depth to the top of the layer (ignores first value since top of 1D model is infinity)
! rho       :: Either a fixed value or else -1, the flag for an inversion free parameter.
!              If the resistivity is fixed, then none of the remaining 3 layer parameters are used 

! penalty       :: The roughness penalty to apply across the layer top boundary. Use zero to allow for a jump
! Preference     :: A linear resistivity Preference value that you can use to push the inversion in a certain direction.  
!                  If you don't want a Preference, set this and the prej_penalty to 0.
! prej_penalty  :: The roughness penalty to apply to the difference of the current model parameter value and the Preference.
!                  If you really want to force the inversion towards the Preference value, make prej_penalty >> penalty for 
!                  the layer
!

!
! Store a copy of the d and sd vectors to use for the OPRA method:
!
    ! real(RealPrec), dimension(:), allocatable :: d_copy, sd_copy   ! using RealPrec from module Occam_realprecision
! JG d & sd come from java, so will be float[]s
    real,   pointer  :: d_copy(:), sd_copy(:)    

    !
    ! Data dependent Tx-Rx variables read from data file in subroutine readData:
    !
    integer                               :: nTxIn, nRxIn, nFreqIn
    !   real(8), dimension(:),  allocatable    :: azimuthTxIn,dipTxIn
    real(8),   pointer     :: azimuthTxIn(:),dipTxIn(:)
    !   real(8), dimension(:),  allocatable    :: xTxIn,yTxIn,zTxIn,xRxIn,yRxIn,zRxIn
    real(8),   pointer     :: xTxIn(:),yTxIn(:),zTxIn(:),xRxIn(:),yRxIn(:),zRxIn(:)
    !   real(8), dimension(:),  allocatable    :: ThetaRxIn,AlphaRxIn, BetaRxIn
    real(8),   pointer     :: ThetaRxIn(:),AlphaRxIn(:), BetaRxIn(:)
    !   real(8), dimension(:),  allocatable    :: fTxIn, sigsite 
    real(8),   pointer     :: fTxIn(:), sigsite(:)

    !
    ! Rx orientation solver variables:
    !
    integer                               :: nRxRotSolve ! Number of Rx to solve for orientations
    !logical, dimension(:), allocatable    :: lRxRotSolve ! [nRxIn] logical array of flags for solving for Rx rotations
    !JG use ints for portability
    integer,   pointer     :: lRxRotSolve(:) ! [nRxRotSolve] array of Rx indexes to solve for rotations


!
! Transformed Tx-Rx arrays used to speed up computations (reduces nTx when dip, azimuth and depth are constant): 
!
    integer                               :: nTxTfm
    !real(8), dimension(:), allocatable    :: xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm
    real(8),   pointer     :: xTxTfm(:), yTxTfm(:), zTxTfm(:), azimuthTxTfm(:), dipTxTfm(:)
    integer                               :: nRxTfm
    ! real(8), dimension(:), allocatable    :: xRxTfm, yRxTfm, zRxTfm
    real(8),   pointer     :: xRxTfm(:), yRxTfm(:), zRxTfm(:)
    !integer, dimension(:), allocatable    :: iRxTfm ! index pointing from Tfm site number to input site number (RxIn)
    integer,   pointer     :: iRxTfm(:) ! index pointing from Tfm site number to input site number (RxIn)
    !   integer, dimension(:),  allocatable    :: RxInd  ! index pointing from Tfm site number to site sent to Dipole1D
    integer,   pointer     :: RxInd(:)  ! index pointing from Tfm site number to site sent to Dipole1D
    !   integer, dimension(:),  allocatable    :: rx1DtoRxIn ! index pointing from Dipole1D site number to RxIn number, used for rotations
    integer,   pointer     :: rx1DtoRxIn(:) ! index pointing from Dipole1D site number to RxIn number, used for rotations
    !   integer, dimension(:,:),  allocatable  :: dpTfm  ! data parameter array for transformed Tx to Rx's
    integer,  pointer      :: dpTfm(:,:) ! data parameter array for transformed Tx to Rx's



    end module csem1d_mod


module Occam
! JG membes of Occam module needed
  integer, parameter, public :: RealPrec = kind(0e0)   ! Double precision:  RealPrec = kind(0d0) 
                                                         ! Single precision:  RealPrec = kind(0e0) 
    !
    ! This module contains objects that are used in the Forward
    ! and sensitivity computations, and within Occam.
    !
    ! *This is the only Occam module that Forward routines need to use.
    
    
    integer       :: nCurrentIter     ! Current Occam iteration number
    character(50) :: cRootName        ! The Occam iteration file root name. 
                                      ! Defaults to 'ITER' if not given
                                      ! as a command line argument

    
    ! Data and model file names:
    character(50) :: dataFileName, modelFileName

    ! Data to be inverted:
    integer                              :: nd
    ! real, dimension(:), allocatable      :: d, sd
	real, pointer :: d(:), sd(:)	
    ! integer, dimension(:,:), allocatable :: dp
    integer, pointer :: dp(:,:)

    ! nd       = number of data.  + 1 or 2 if have statics
    ! d(nd)    = data 
    ! sd(nd)   = data errors (one standard error) - usually read from data file
    ! dp(nd,:) = array of parameters associated with each datum (freq., etc)
    !            Note that dp() is not used in Occam, only in the forward routines
    !
    !-----------------------------------------------------------------------
    ! model parameters for the current iteration and its response
    integer                           :: nParams
    ! real, dimension(:), allocatable   :: pm, dm 
    real, pointer :: pm(:), dm(:)
    ! real, dimension(:,:), allocatable :: wj
    real, pointer :: wj(:,:)
    !   nParams        = number of free model parameters
    !   pm(nParams)    = array of model parameters
    !   dm(nd)         = array of model response (data domain)
    !   wj(nd,nParams) = model response Jacobian matrix
    !-----------------------------------------------------------------------
    ! associated data/model parameters for the current iteration 
    ! these are things that are derived from the current model and/or data
    ! but are not used for the model roughness or misfit computations.
    ! These are only along for the ride in the Occam routines and are enirely
    ! dealt with in the Forward code.  They need to be included here
    ! since Occam can keep track of pm_assoc for the best fitting model
    ! of each iteration.
    integer  :: npm_assoc
    ! real, dimension(:), allocatable     :: pm_assoc  ! associated model responses (not inverted for but important to Fwd responses)
	real, pointer :: pm_assoc(:)  ! associated model responses (not inverted for but important to Fwd responses)	
    
    !real, dimension(:), allocatable     :: prewts, premod   
    real, pointer     :: prewts(:), premod(:)   	

    !   prewts(nParams) = weights for model prejudice
    !   premod(nParams) = model prejudice


!
!------------------------
! Contained subroutines:
!------------------------
! Thes variables above have global scope within the contained subroutines below.
!
    contains
 
!-----------------------------------------------------------------------
! Public subroutines
!-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    subroutine Lower( s )
    !-----------------------------------------------------------------------
    ! OCCAM 3.0 Package
    ! David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine Revision 3.0, November 2006
    ! DGM Nov 2006 - convert string to lower case
     
    character(*), intent(out)  :: s
    integer i
    do  i=1,len_trim(s)
      if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
        s(i:i) = char(ichar(s(i:i)) + 32)
      endif
    enddo
    
    end subroutine Lower

    
end module Occam



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
!   JGHACK make pointers
!    real(8), dimension(:), allocatable  :: sig1D     ! (S/m) Layer conductivities 
!    real(8), dimension(:), allocatable  :: zlay1D    ! (m)   Depth to top of each layer, first layer ignored 
    real(8), pointer :: sig1D(:)                      ! (S/m) Layer conductivities 
    real(8), pointer :: zlay1D(:)                     ! (m)   Depth to top of each layer, first layer ignored 
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
    lUseSpline1D    = .true.         ! Use spline interpolation for faster 1D computations
    linversion      = .false.        ! Compute derivatives with respect to sigma(layers)
    
    end subroutine init_defaults_Dipole1D 
    
    
    end module dipole1d_io

!==============================================================================!    
!==================================================================! MT1D_MOD  ! 
!==============================================================================!
    module mt1d_mod
!
! Module for 1D MT response and layer deriviatives:
!
    real(8) :: zRx_MT
    real(8) :: freqMT
    
    complex(8) :: CResp
    complex(8), dimension(:), allocatable :: dCRespdRho  
     
    
    end module mt1d_mod

