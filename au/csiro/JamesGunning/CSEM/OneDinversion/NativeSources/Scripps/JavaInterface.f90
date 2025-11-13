! 
! Interface of fortrn 90 codes to Java drivers
! Uses modules defined in JavaModules.f90
! James Gunning Oct 08

subroutine set_dipole_mod_pointers(zlay1D_ptr, sig1D_ptr)                              ! allocateL
     
  use dipole1D ! Model nLay1D 
  implicit none

  real(8), target :: zlay1D_ptr(nlay1D) !  double[nlay1D]
  real(8), target :: sig1D_ptr(nlay1D)  !  double[nlay1D]

  zlay1D => zlay1D_ptr
  sig1D => sig1D_ptr  

end subroutine set_dipole_mod_pointers


subroutine computeFwdWrapper( bDoPartialsInt, currentMod )
  use occam                                     ! for nParams

  integer, intent(in) :: bDoPartialsInt        ! int flag for computing Jacobians, c style. 0=false
  real, intent(in)    :: currentMod(nParams)   ! input model parameters

  if (bDoPartialsInt == 0) then
     call computeFwd(.false., currentMod )
  else
     call computeFwd(.true., currentMod )
  endif

end subroutine computeFwdWrapper


subroutine set_csem1d_mod_pointers( &
     nTxIn_v, azimuthTxIn_ptr, dipTxIn_ptr, xTxIn_ptr, yTxIn_ptr, zTxIn_ptr,              & ! allocateTX
     nFreqIn_v,fTxIn_ptr,                                                                 & ! allocateF  
     nRxIn_v, xRxIn_ptr, yRxIn_ptr, zRxIn_ptr, ThetaRxIn_ptr, AlphaRxIn_ptr, BetaRxIn_ptr,& ! allocateR,
     nRxTfm_v, xRxTfm_ptr, yRxTfm_ptr, zRxTfm_ptr, iRxTfm_ptr,                            & ! allocateRfm
     nd, dpTfm_ptr,                                                                       & ! allocate_dpTfm,
     nTxTfm_v, xTxTfm_ptr, yTxTfm_ptr, zTxTfm_ptr, azimuthTxTfm_ptr, dipTxTfm_ptr,        & ! allocateTfm 
     d_copy_ptr, sd_copy_ptr,                                                             & ! allocateDcopies 
     sigsite_ptr, RxInd_ptr, rx1DtoRxIn_ptr,                                              & ! allocateDM
     nlay1D_v, layer_params_ptr,                                                          & ! allocateL
     nRxRotSolve_v,lRxRotSolve_ptr                                                      & ! allocatelRxRotSolve
     )
  use csem1d_mod
  use dipole1d
 
  implicit none

  ! Arguments, in order corresponding to the class allocation calls
  
  ! allocateTX
  integer :: nTxIn_v
  real(8), target :: azimuthTxIn_ptr(nTxIn_v)
  real(8), target :: dipTxIn_ptr(nTxIn_v)
  real(8), target :: xTxIn_ptr(nTxIn_v)
  real(8), target :: yTxIn_ptr(nTxIn_v) 
  real(8), target :: zTxIn_ptr(nTxIn_v)

  !  allocateF  
  integer :: nFreqIn_v
  real(8), target :: fTxIn_ptr(nFreqIn_v)

  !allocateR
  integer :: nRxIn_v
  real(8), target ::   xRxIn_ptr(nRxIn)
  real(8), target ::   yRxIn_ptr(nRxIn)
  real(8), target ::   zRxIn_ptr(nRxIn)
  real(8), target ::   ThetaRxIn_ptr(nRxIn)
  real(8), target ::   AlphaRxIn_ptr(nRxIn)
  real(8), target ::   BetaRxIn_ptr(nRxIn)

  ! allocateRfm
  integer :: nRxTfm_v
  real(8), target :: xRxTfm_ptr(nRxTfm_v)
  real(8), target :: yRxTfm_ptr(nRxTfm_v)
  real(8), target :: zRxTfm_ptr(nRxTfm_v)
  integer, target :: iRxTfm_ptr(nRxTfm_v)

  !allocate_dpTfm
  integer :: nd
  integer, target :: dpTfm_ptr(nd,2)

  ! allocateTfm(int nTxTfm)
  integer :: nTxTfm_v
  real(8), target ::   xTxTfm_ptr(nTxTfm_v)
  real(8), target ::   yTxTfm_ptr(nTxTfm_v)
  real(8), target ::   zTxTfm_ptr(nTxTfm_v)
  real(8), target ::   azimuthTxTfm_ptr(nTxTfm_v)
  real(8), target ::   dipTxTfm_ptr(nTxTfm_v)
  
  ! allocateDcopies(int nd)
  real, target ::   d_copy_ptr(nd)
  real, target ::   sd_copy_ptr(nd)

  ! allocateDM(int nRxTfm)
  real(8), target ::  sigsite_ptr(nRxTfm_v)    !  double[nRxTfm];
  integer, target ::  RxInd_ptr(nRxTfm_v)      !  int[nRxTfm];
  integer, target :: rx1DtoRxIn_ptr(nRxTfm_v)  !  int[nRxTfm];


  !allocateL(int nlay1D)
  integer :: nlay1D_v
  real(8), target :: layer_params_ptr(nlay1D_v,5) !  double[nlay1D][5];

  !allocatelRxRotSolve
  integer :: nRxRotSolve_v
  integer, target :: lRxRotSolve_ptr(nRxRotSolve_v) !  int[]

  ! internal variables
  integer :: i

  ! Here, the main body.
  ! Set all the module pointers to the alloced memory passed in
  

  ! allocateTX
  nTxIn = nTxIn_v
  azimuthTxIn => azimuthTxIn_ptr
  dipTxIn => dipTxIn_ptr
  xTxIn => xTxIn_ptr
  yTxIn => yTxIn_ptr
  zTxIn => zTxIn_ptr




  ! allocateF
  nFreqIn = nFreqIn_v
  fTxIn => fTxIn_ptr

  ! allocateR
  nRxIn = nRxIn_v
  xRxIn => xRxIn_ptr 
  yRxIn => yRxIn_ptr 
  zRxIn => zRxIn_ptr 
  ThetaRxIn => ThetaRxIn_ptr 
  AlphaRxIn => AlphaRxIn_ptr 
  BetaRxIn => BetaRxIn_ptr 



  ! allocateRfm
  nRxTfm = nRxTfm_v  
  xRxTfm => xRxTfm_ptr
  yRxTfm => yRxTfm_ptr
  zRxTfm => zRxTfm_ptr
  iRxTfm => iRxTfm_ptr
 

  ! allocate_dpTfm(int nd)  
  dpTfm => dpTfm_ptr ! int[nd][2];


  ! allocateTfm(int nTxTfm)  
  nTxTfm = nTxTfm_v
  xTxTfm => xTxTfm_ptr !  double[nTxTfm]; 
  yTxTfm => yTxTfm_ptr !  double[nTxTfm]; 
  zTxTfm => zTxTfm_ptr !  double[nTxTfm]; 
  azimuthTxTfm => azimuthTxTfm_ptr !  double[nTxTfm]; 
  dipTxTfm => dipTxTfm_ptr !  double[nTxTfm];


  ! allocateDcopies(int nd)
  d_copy => d_copy_ptr !  float[nd];
  sd_copy => sd_copy_ptr !  float[nd];

  ! allocateDM(int nRxTfm)
  nRxTfm = nRxTfm_v
  sigsite => sigsite_ptr !  double[nRxTfm];
  RxInd => RxInd_ptr !  int[nRxTfm];
  rx1DtoRxIn => rx1DtoRxIn_ptr !  int[nRxTfm];

  !allocateL(int nlay1D)
  nlay1D = nlay1D_v
  layer_params => layer_params_ptr !  double[nlay1D][5];	  


  nRxRotSolve = nRxRotSolve_v 
  lRxRotSolve => lRxRotSolve_ptr

  !write(*,fmt='(A,(1000g15.4))') 'azimuthTxIn_ptr: ',(azimuthTxIn_ptr(i),i=1,n)
  !stop

end subroutine 

subroutine do_eb_field_allocs(n)
  use Occam
  use csem1d_mod
  use dipole1d

  implicit none

  integer, intent(in) :: n
  integer :: nAllocErr

  ! Complex allocations handled inside fortran....
  ! need to deallocate this later
  nRxTfm = n

  ! real vectros
  allocate( x1D(nRxTfm),y1D(nRxTfm),z1D(nRxTfm), stat=nAllocErr)
  if (nAllocErr .ne. 0) then
     write(*,*) 'Out of memory.  Too many receivers  (', nRxTfm, ')'
     stop 
  endif

  allocate( ex1D(nRxTfm),ey1D(nRxTfm),jz1D(nRxTfm), &
       & bx1D(nRxTfm),by1D(nRxTfm),bz1D(nRxTfm), &
       & stat=nAllocErr )  
			
  if (nAllocErr .ne. 0) then
     write(*,*) 'Out of memory.  Too many receivers  (', nRxTfm, ')'
     stop 
  endif
end subroutine do_EB_field_allocs


subroutine do_jacobian_allocs(nRxTfm, nlay)
  use dipole1d
  
  implicit none

  integer, intent(in) :: nRxTfm, nlay
  integer :: nAllocErr

  ! Complex allocations handled inside fortran....
  ! need to deallocate this later

  allocate( dexdsig(nRxTfm,nlay),deydsig(nRxTfm,nlay),djzdsig(nRxTfm,nlay), &
       & dbxdsig(nRxTfm,nlay),dbydsig(nRxTfm,nlay),dbzdsig(nRxTfm,nlay), &
       &  stat=nAllocErr )             
  
  if (nAllocErr .ne. 0) then
     write(*,*) 'Out of memory.  Too many receivers x 1D layers  (', nRxTfm,' x ',nlay, ')'
     stop 
  endif  
end 

subroutine free_jacobian_allocs()
  use dipole1d
  implicit none

  deallocate( dexdsig,deydsig,djzdsig,dbxdsig,dbydsig,dbzdsig)             

end 


subroutine set_occam_mod_pointers(np, n_data, cnDataParams, dp_ptr, d_ptr, sd_ptr,dm_ptr)
  use occam
  implicit none

  ! allocateD
  integer :: np,n_data,cnDataParams
  integer, target :: dp_ptr(n_data,cnDataParams)
  real, target :: d_ptr(n_data)
  real, target :: sd_ptr(n_data)
  real, target :: dm_ptr(n_data)

  nd = n_data
  nParams = np
  dp => dp_ptr
  d => d_ptr
  sd => d_ptr
  dm => dm_ptr

end subroutine

! Wrapper to make call to init_defaults_Dipole1D compiler-independent
! from c
subroutine set_defaults_dipole1D()
  use dipole1d

  call init_defaults_Dipole1D
end subroutine

subroutine set_dm_wj_pointer(dm_ptr,g_ptr)
  use Occam
  implicit none

  real, target :: dm_ptr(nd)
  real, target :: g_ptr(nd,nParams)

  ! point at passed in memory
  dm => dm_ptr
  wj => g_ptr

end subroutine
 
subroutine set_dm_pointer(dm_ptr)
  use Occam
  implicit none

  real, target :: dm_ptr(nd)

  ! point at passed in memory
  dm => dm_ptr

  ! zero thus for safety
  nullify(wj)

end subroutine

! Some extra needed functions removed from Occam.f90


!==============================================================================!  
!==================================================================! parseLine !  
!==============================================================================! 
subroutine parseLine( nLen, sLine, bComment )
  !   
  ! Subroutine to check if the sLine is blank or a comment line, and if it isn't
  ! then any comment at the end of the line is blanked. 
  !
  ! This is styled after D. Myer's parseCode function.
  !
  ! This is useful for reading in data tables that have user comment lines 
  ! or comments at the end of the line, denoted by the ! and % symbols.
  ! 
  ! If the entire line is a comment, bComment = .true.
  !
  ! Kerry Key
  ! Scripps Institution of Oceanography
  ! kkey@ucsd.edu
  !
  ! Version 1.0.   April, 2008.
  !
  implicit none

  ! Args
  integer, intent(in)     :: nLen
  character(nLen)         :: sLine
  logical, intent(out)    :: bComment

  ! Local vars
  integer :: iFrom, iTo

  ! Init returns
  bComment = .false.


  ! Convert all tab characters to spaces
  forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '

  ! Skip any beginning blanks
  do iFrom = 1,nLen
     if (sLine(iFrom:iFrom) .ne. ' ') exit
  enddo

  ! If the first char is a comment char, then the whole line is a comment.
  ! DGM April 2008 Also, if the line is blank, consider it a comment.
  if( iFrom >= nLen .or. sLine(iFrom:iFrom) == '%' &
       .or. sLine(iFrom:iFrom) == '!' ) then
     bComment = .true.
     return
  endif

  ! Now trim off any comments at the end of the line  
  iTo = len_trim(sLine)
  iFrom = index(sLine,'%')
  if (iFrom > 0 .and. iFrom < iTo) then
     sLine(iFrom:iTo) = ' '
  endif
  iFrom = index(sLine,'!')
  if (iFrom > 0 .and. iFrom < iTo) then
     sLine(iFrom:iTo) = ' '
  endif

end subroutine parseLine
