! Standard modules for CSEM1D and MT1D with allocatable memory, for use
! in stand-alone f90 codes
! James Gunning Oct 08

    module csem1d_mod

!
! Module to store various parameters for 1D CSEM inversion
! This is stuff that is not used in Occam.f90 and not Dipole1D.f90
! but is required for interfacing between them.

!   Magic codes and constants in this include
    include 'CSEM1DParameters.inc'


! While the data parameter array dp is passed in/out of Occam, it is only used
! in the forward routines.  So we set cnDataParams in the forward routines
! and then allocate the dp array in readData.
!        
    integer :: cnDataParams
 
    real(8), dimension(:,:), allocatable :: layer_params 
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
!                  If you really want to force the inversion towards the Preference value, make prej_penalty >> penalty for the layer
!

!
! Store a copy of the d and sd vectors to use for the OPRA method:
!
    real, dimension(:), allocatable :: d_copy, sd_copy 
       
!
! Data dependent Tx-Rx variables read from data file in subroutine readData:
!
    integer                               :: nTxIn, nRxIn, nFreqIn
    real(8), dimension(:), allocatable    :: azimuthTxIn,dipTxIn
    real(8), dimension(:), allocatable    :: xTxIn,yTxIn,zTxIn,xRxIn,yRxIn,zRxIn
    real(8), dimension(:), allocatable    :: ThetaRxIn,AlphaRxIn, BetaRxIn
    real(8), dimension(:), allocatable    :: fTxIn, sigsite 


!
! Transformed Tx-Rx arrays used to speed up computations (reduces nTx when dip, azimuth and depth are constant): 
!
	integer                               :: nTxTfm
	real(8), dimension(:), allocatable    :: xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm
 	integer                               :: nRxTfm
 	real(8), dimension(:), allocatable    :: xRxTfm, yRxTfm, zRxTfm
 	integer, dimension(:), allocatable    :: iRxTfm ! index pointing from Tfm site number to input site number (RxIn)
    integer, dimension(:), allocatable    :: RxInd  ! index pointing from Tfm site number to site sent to Dipole1D
    integer, dimension(:), allocatable    :: rx1DtoRxIn ! index pointing from Dipole1D site number to RxIn number, used for rotations
 	integer, dimension(:,:), allocatable  :: dpTfm  ! data parameter array for transformed Tx to Rx's
    
!
! Rx orientation solver variables:
!
    integer :: nRxRotSolve  ! Number of Rx to solve for orientations
    integer, dimension(:), allocatable    :: IndRxRotSolve ! [nRxRotSolve] array of Rx indexes to solve for rotations


    end module 
    
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

