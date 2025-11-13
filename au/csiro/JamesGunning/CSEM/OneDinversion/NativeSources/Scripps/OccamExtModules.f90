
! Standard module for external driving of occam with explicitly allocatable memory
! James Gunning Oct 08

module Occam_interface_mod
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
    real, dimension(:), allocatable      :: d, sd
    integer, dimension(:,:), allocatable :: dp
    ! nd       = number of data.  + 1 or 2 if have statics
    ! d(nd)    = data 
    ! sd(nd)   = data errors (one standard error) - usually read from data file
    ! dp(nd,:) = array of parameters associated with each datum (freq., etc)
    !            Note that dp() is not used in Occam, only in the forward routines
    !
    !-----------------------------------------------------------------------
    ! model parameters for the current iteration and its response
    integer                           :: nParams
    real, dimension(:), allocatable   :: pm, dm 
    real, dimension(:,:), allocatable :: wj
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
    real, dimension(:), allocatable     :: pm_assoc  ! associated model responses (not inverted for but important to Fwd responses)
    
    real, dimension(:), allocatable     :: prewts, premod   
    !   prewts(nParams) = weights for model prejudice
    !   premod(nParams) = model prejudice
    
end module
