! Occam_Mod.f90 - module of shared memory objects.
!   Accumulated from a variety of old common blocks and made
!   so that the memory can be allocated at runtime instead of
!   at compile time.
!
! David Myer
! Oct 2006
! 
! Revision  Date           Author      Notes
!
! 2.0       Sept 4, 2008   Kerry Key   Split into Occam_mod into two modules to isolate 
!                                      critical working arrays from the Forward code.                                     
!                                      Added Occam_Timers module, nCurrentIter and cRootName.
!
! 1.1       May 2008       Kerry Key   Added in associated parameter array
!
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!
! Everything below here should only be used by the Occam code. 
!
! Forward codes stay outta here!
!
!----------------------------------------------------------------------
module Occam_internal_mod 
    !
    ! This module contains objects used only within the Occam routines:
    !
    logical                              :: gbAddDiags       ! DGM Oct 2006
    ! gbAddDiags = true if diagonal penalties should be added to horiz & vert.
    

    !-----------------------------------------------------------------------
    ! stuff to be carried into critical subroutines 
    real, dimension(:,:), allocatable   :: ptp, wjtwj
    real, dimension(:), allocatable     :: wjtwd, pwk1, pwk2, pwk3, dwk1, dwk3
    real                                :: frac
    integer                             :: iout, idebug, nfor, nFracTimes
    real, dimension(:), allocatable     :: dwk1_assoc, dwk2_assoc, dwk3_assoc ! working arrays for passing associated parameters

    !   ptp() = product of penalty matrix with its transpose
    !   wjtwj() = product of weighted jacobian matrix with its transpose
    !   wjtwd() = product of weighted jacobian matrix with weighted translated data
    !   pwk1,2,3() = vectors to carry around working copies of models
    !   frac = current reduction in step size
    !   nfor = tally for total number of forward calculations
    !   iout = unit number of output (usually LOGFIL)
    !   idebug = current debug level
    !   dwk1,3 = working copies of the model responses - 'dm' is used as dwk2
    !   dwk1_assoc,2,3= working copies of dm_assoc
    !   nFracTimes = # of times to cut stepsize before giving up (DGM Nov 2006 - was hardcoded)

    !-----------------------------------------------------------------------
    ! Working variables used by TOFMU - a function that is called many times
    !   per occam iteration and which forms 90% of the program run time.
    !   These vars were allocated automatically at every call, but it is more
    !   efficient to do it one time and let the memory stand.
    real, dimension(:,:), allocatable   :: aMat
    real, dimension(:), allocatable     :: pwk4, dwk4
    
end module Occam_internal_mod


!----------------------------------------------------------------------

module Occam_ModelLimits

! DGM Nov 2006 - Values read from the startup iteration file controlling
!   the span of values the model can take as well as the number of steps
!   per decade.  Values are read in ItrIn() and used in TOFMU().

    logical ::  gbModelLimited
    real    ::  gnModelMin, gnModelMax
    logical ::  gbModelStepped
    real    ::  gnModelSteps  
contains

    subroutine ModelLimits_Init
        implicit none
        gbModelLimited  = .false.
        gbModelStepped  = .false.
        gnModelMin      = 0
        gnModelMax      = 0
        gnModelSteps    = 0
    end subroutine ModelLimits_Init
end module Occam_ModelLimits

!----------------------------------------------------------------------

    module Occam_Timers
!
!  Module to store some timers for the Occam program and iterations
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    real(8) :: time0     ! Start time of Occam program (in seconds)
    real(8) :: timeIter0 ! Start time of an Occam iteration 
    
    contains 

    subroutine  get_time_offset(timein,timeout)
!    
! timein is the clock start time or 0.
! timeout is the time in seconds since the input time
!
! Version 2.0  February 25, 2008  Now uses date_and_time Fortran intrinsic
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu

   implicit none
   
   integer, dimension(8) :: values
   integer               :: i,j,k,mjd
   
   real(8) :: timein, timeout, fracday
    
 !
 ! Fortran90 Time function:
 !
    call date_and_time(values=values) !this ouputs only values, ignoring other optional arguments
    
 ! Convert year, month day to modified julian day:
      
    i = values(1)
    j = values(2)
    k = values(3)
    mjd = -678927 + k + 1461*(i+(j-14)/12)/4 + 367*(j-2-12 * &
          & ((j-14)/12))/12 + (24002-12*i-j)/1200
    
               
  ! Add on fractional day:
                ! hour            ! minute          ! sec       ! millisec
    fracday = ((values(5)*60.d0 + values(6))*60.d0 + values(7) + values(8)/1000.d0 )/86400.d0
    
    timeout = mjd + fracday
    
  ! Finally, convert timeout to time difference between it and timein:  
    timeout =  timeout*86400.d0  - timein
               
           
    end subroutine  get_time_offset
    
    end module Occam_Timers

