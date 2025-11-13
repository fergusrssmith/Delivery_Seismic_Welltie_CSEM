!    This file is part of LegacyScrippsDipole1D
!    LegacyScrippsDipole1D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    LegacyScrippsDipole1D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with LegacyScrippsDipole1D.  If not, see <http://www.gnu.org/licenses/>.

! Forward modeling package for the 1D earth 
! Modification of 10 June 2004 version by Steve Constable (Scripps)
! Fixes and augmentations by James Gunning (CSIRO) 
! and John Hedditch 
!
! Legacy SC notes:
! tangle not yet defined in new file format 
! need to fill derivative matrix not by layer but by free parameter
!
! Hacks by James Gunning: forward modelling & sensitivities of Re & Im
! components, to prevent phase wrapping wierdness.
! John Hedditch - fix of phase derivatives

!-----------------------------------------------------------------------
subroutine formod(np,nd,pm,dp,dm,d)
  !----------------------------------------------------------------------- 
  ! common replacements....
  use fmodel
  use mode

  !
  !   thickf: array of layer thicknesses, layer 0 = seawater
  !   rho: array of layer resistivities, layer 0 = seawater
  !   freqs: array of frequencies to compute
  !   thite: transmitter height, in meters
  !   tangle: transmitter angle , in degrees
  !   ipen1: array of penalties between layers, for inversion
  !   ipen2: array of penalties for prejudice structure, for inversion
  !   iff: array denoting which layers are free parameters
  !   nfreq: number of frequencies in data array
  !   nlayer: number of layers, including 1/2 space but not seawater

  !     arguments

  integer np, nd
  ! JGHACK  real*4 pm(npp),dp(ndd,6),dm(ndd)
  ! new...
  real*4 pm(np),dp(nd,6),dm(nd)
  ! JH HACK
  real*4 d(nd)
  real p0,pp360,pm360

  include 'occamdim.inc'

  ! dehhed stuff:
  integer nlyrsp,nrngsp
  parameter (nlyrsp=npp ,nrngsp=ndd )
  integer lunin,lunout,lvlout,nel,nrange,ifdiffa
  real*8 asigma(0:nlyrsp),athick(0:nlyrsp)
  real*8 aranges(nrngsp),afreq,dphi,zxmit,zrecv,reqerr 
  complex*16 er(nrngsp),ep(nrngsp)
  complex*16 derds(nlyrsp,nrngsp),depds(nlyrsp,nrngsp)
  complex*16 derdt(nlyrsp,nrngsp),depdt(nlyrsp,nrngsp)

  ! local variables
  integer i, ii, k, lut(ndd), ityp
  character*10 how
  real pmax, sqrt2, phia, rot, rangle, rrot, unwrap, freq
  real ta, arho, aphi, xreal, yreal, ximag, yimag, xamp, yamp
  real rmax, rmin
  parameter (rmax = 7., rmin = -5.)
  logical csemflg
  complex c

  integer imin, imax,i_g  
  real zxmitt,tstart,tfinish

  ! set water conductivity and depth      
  if (rho(0) .ne. 0.0) then
     asigma(0) = 1/rho(0)
     csemflg = .true.
  else
     csemflg = .false.
  end if 
  athick(0) = thickf(0)
  nel = nlayer

  ! update the free parameters in the model
  !  if (rmode .eq. 'occam     ') then
  if (.true.) then
     ! write(*,*) ' FIXME 1: OCCAM mode hard coded: fix rmode passing...'
     do i = 1,np
        rho(iff(i)) = 10.**min(rmax,max(rmin,pm(i)))
     end do
     do  i = 1,nlayer-1
        athick(i) = thickf(i)
     end do
     do  i = 1, nlayer
        asigma(i) = 1./rho(i)
     end do
  else if (rmode .eq. 'marquardt ') then
     if (np .ne. (2*nlayer - 1)) &
       write(*,*) ' layers don''t match parms'
     do i = 1,nlayer
        rho(i) = 10.**min(rmax,max(rmin,pm(i)))
     end do
     do  i = 1,nlayer-1
        athick(i) = 10.**min(rmax,max(rmin,pm(i+nlayer)))
        thickf(i) = athick(i)
     end do
     do  i = 1, nlayer
        asigma(i) = 1./rho(i)
     end do
  else
     write(*,*) ' unknown calling program type'
  end if

  !  set level of debug-output, logical units for any i/o (csem code)
  lvlout=0
  lunin=5
  lunout=6
  !  set computation method, requested accuracy (csem code)
  !  (the possibilities for how are 'filter', 'quad', 'compare'. 
  !  set the method for the default by using whatever other string
  how = 'filter'
  reqerr=0.1e-12
  !  ask for forward modelling only (no derivatives) (csem code)
  ifdiffa=0
  nrange = 1
  ! set azimuth to 45 degrees to get both components of the field
  ! (will be decomposed later for our specific needs) (csem code)
  dphi = 45.
  ! set the receiver height (csem code)
  zrecv = 0.0
  sqrt2 = sqrt(2.)

  !  put model parameters etc. in common:
  call hedfix(nel,asigma,athick,lunin,lunout,lvlout)

  !     JG timing tests
  !call cpu_time(tstart) 

  !     JG, replace this with
  !     do i = 1, nd
  i=1
100 if (i.le.nd) then

     ityp = int(dp(i,1)+0.01)

     if (ityp .lt. 20) then
        !         compute CSEM response        
        afreq = dp(i,3)
        aranges(1) = dp(i,2)
        zxmit = max(1.0, dp(i,6))
        imin=i

        if (.not.use_range_lumping) then
           imax=i  
           !           Compute fields one at a time...
           !           compute the field & desired derivatives at the given coordinates

           call hedsfi(afreq,nrange,aranges,dphi,zxmit,zrecv &
           ,nlyrsp,nrngsp,how,reqerr,ifdiffa &
           ,er,ep,derds,derdt,depds,depdt)

! write(*,*) 'formod time (uS) = ',1.0d6*(tfinish-tstart) 
        else  
           !           Compute fields in batches of common frequency and Tx height..
           imin = i
           imax = nd
           !           Find a cluster...
           do ii = imin+1, nd 
              !              transmitter height - FIXME! should be dp(ii,6), or suitable cluster average....
              zxmitt = dp(ii,6)
              if ((abs(dp(ii,3) - afreq) > 1.0d-6) .or. (abs(max(1.0, zxmitt) - zxmit) > 1.0d-6) &
                    .or.(.not.((int(dp(ii,1)+0.01)) < 20))) then
                 !                different freq or transmitter height or not CSEM data
                 imax = ii-1 
                 goto 500    
              else
                 aranges(ii-imin+1) = dp(ii,2)
                 !                dumb hack for ranges monotonicity....
                 if (.not.(aranges(ii-imin+1).gt.aranges(ii-imin)))then
                    aranges(ii-imin+1) = aranges(ii-imin)*(1.d0+1.0d-15)
                 end if
              end if
           end do
500        continue
           !           Groupable data spans imin..imax
           nrange = imax - imin+1
           !           compute the field & desired derivatives at the given coordinates
           !           (several ranges can be done at once)

           !call cpu_time(tstart) 
           ! Dominant time in models from this call:
           call hedsfi(afreq,nrange,aranges,dphi,zxmit,zrecv &
           ,nlyrsp,nrngsp,how,reqerr,ifdiffa &
           ,er,ep,derds,derdt,depds,depdt)

           !call cpu_time(tfinish)
           !tTot = tTot+(tfinish-tstart) !JGHACK_TIMING
        end if

        !         loop over lumped group
        i_g = 1
        do ii = imin,imax
           phia = dp(ii,4)
           rot = dp(ii,5)
           ityp = int(dp(ii,1)+0.01)

           !            if (int(dp(i,1)) .le. 6) then
           !            JGHACK for Re + Im
           if ((ityp .le. 6).or.((ityp.ge.14).and.(ityp.le.17))) then

              !              the azimuth has been fixed at 45 degrees.
              !              now, if the data type is 1-6 (or 14..17, JG) we want the erho and ephi terms,
              !              obtained simply be multiplying all the components by sqrt(2):
              er(i_g) = er(i_g)*sqrt2
              ep(i_g) = ep(i_g)*sqrt2
           else
              !              otherwise we want the erho and ephi components, obtained by
              !              further correcting the fields for azimuthal dependence of dipoles:
              er(i_g) = cos(phia/57.29578)*er(i_g)*sqrt2
              ep(i_g) = sin(phia/57.29578)*ep(i_g)*sqrt2
           end if

           !            now do the right thing for this data type:
           !            goto (1,2,3,4,5,6,7,8,9,10,11,12,13) ityp
           goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) ityp
           write(*,*) 'component flag invalid in dp()'
           stop

           !            erho magnitude
1          dm(ii) = abs(er(i_g))
           goto 310

           !            erho phase
2          dm(ii) = atan2(dimag(er(i_g)),dreal(er(i_g)))*57.29578
           ! JH HACK
           ! + JG corrections: check pp360, pm360 defs
           p0 = abs( dm(ii) - d(ii) )
           pm360 = abs( dm(ii) - 360.0 - d(ii) )
           pp360 = abs( dm(ii) + 360.0 - d(ii) )
           if( pm360.lt.p0 ) then
             dm(ii) = dm(ii) - 360.0
           elseif ( pp360.lt.p0) then
             dm(ii) = dm(ii) + 360.0
           endif
           goto 310

           !            ephi magnitude
3          dm(ii) = abs(ep(i_g))
           goto 310

           !            ephi phase
4          dm(ii) = atan2(dimag(ep(i_g)),dreal(ep(i_g)))*57.29578
           ! JH HACK
           !+ JG corrections: check pp360, pm360 defs
           p0 = abs( dm(ii) - d(ii) )
           pm360 = abs( dm(ii) - 360.0 - d(ii) )
           pp360 = abs( dm(ii) + 360.0 - d(ii) )
           if( pm360.lt.p0 ) then
             dm(ii) = dm(ii) - 360.0
           elseif ( pp360.lt.p0) then
             dm(ii) = dm(ii) + 360.0
           endif
           goto 310

           !            erho log10 magnitude
5          dm(ii) = log10(abs(er(i_g)))
           goto 310

           !            ephi log10 magnitude
6          dm(ii) = log10(abs(ep(i_g)))
           goto 310

           !            log of major axis of pe
7          call pemax(er(i_g),ep(i_g),pmax)
           dm(ii) = log10(pmax)
           goto 310

           !            linear major axis of pe
8          call pemax(er(i_g),ep(i_g),pmax)
           dm(ii) = pmax
           goto 310

           !            log of erho component
9          dm(ii) = log10(abs(er(i_g)))
           goto 310

           !            log of ephi component
10         dm(ii) = log10(abs(ep(i_g)))
           goto 310

           !            log of x-component
11         continue
           ta = tangle - 180.
           arho = (ta + phia + 90.)/57.29578
           aphi = (ta + phia + 180.)/57.29578
           xreal = cos(arho)*real(er(i_g)) + cos(aphi)*real(ep(i_g))
           ximag = cos(arho)*dimag(er(i_g))+ cos(aphi)*dimag(ep(i_g))
           xamp = sqrt(xreal**2 + ximag**2)
           !            xphz = atan2(ximag,xreal)*57.29578
           dm(ii) = log10(xamp)
           goto 310
           ! log of y-component
12         continue
           ta = tangle - 180.
           arho = (ta + phia + 90.)/57.29578
           aphi = (ta + phia + 180.)/57.29578
           yreal = sin(arho)*real(er(i_g)) + sin(aphi)*real(ep(i_g))
           yimag = sin(arho)*dimag(er(i_g))+ sin(aphi)*dimag(ep(i_g))
           yamp = sqrt(yreal**2 + yimag**2)
           ! yphz = atan2(yimag,yreal)*57.29578
           dm(ii) = log10(yamp)
           goto 310

13         continue
           !            log of single component data
           !            rot = receiver rotation - angle +ve anti-clocwise to receiver dipole 
           !            from erho
           rrot = rot / 57.29578
           yreal = cos(rrot)*dreal(er(i_g))+ sin(rrot)*dreal(ep(i_g))
           yimag = cos(rrot)*dimag(er(i_g))+ sin(rrot)*dimag(ep(i_g))
           yamp = sqrt(yreal**2 + yimag**2)
           dm(ii) = log10(yamp)
           !            calculate receiver orientaion from n as a check.
           rangle = tangle - 180. + phia + rot
           rangle = unwrap(3,rangle)
           !            write(*,*)'receiver orientaion:', rangle
           goto 310

           !            JGHACKS for Re + Im 
           !            erho Re
14         dm(ii) = dreal(er(i_g))
           goto 310

           !            erho Im
15         dm(ii) = dimag(er(i_g))
           goto 310
           !            ephi Re
16         dm(ii) = dreal(ep(i_g))
           goto 310

           !            ephi Im
17         dm(ii) = dimag(ep(i_g))
           goto 310


310        continue

           i_g = i_g+1
        end do                ! loop over group ii=imin...imax

     else if (ityp .le. 30) then
        !         MT data          
        freq = 1./dp(i,2)
        call mtlay(freq,c,nlayer,rho,thickf)
        goto (21, 22, 23, 24, 25, 26, 27) ityp-20      
        write(*,*) 'component flag invalid in dp()'
        stop
21      dm(i) = log10(7.8957e-06*(abs(c)**2)*freq)
        goto 1000
22      dm(i) = 57.29578*atan2(real(c),-aimag(c))
        goto 1000
23      dm(i) = real(c)
        goto 1000
24      dm(i) = -aimag(c)
        goto 1000
25      dm(i) = log10(real(c))
        goto 1000
26      dm(i) = log10(-aimag(c))
        goto 1000
27      dm(i) = 7.8957e-06*(abs(c)**2)*freq

1000    continue

     else
        write(*,*) 'component flag invalid in dp()'
        stop        
     end if

     if (imax.gt.imin) i=imax ! makes next increment good   
     i = i+1
     goto 100  
  endif
  !     end do

! write(*,*) 'formod time (uS) = ',1.0d6*(tfinish-tstart) 
  return
end subroutine

!-----------------------------------------------------------------------
subroutine fordiv(np,nd,pm,dp,dm,g,d)
  ! common replacements....
  use fmodel
  use mode

  implicit none

  integer np, nd
  ! JGHACK   real*4 pm(npp),dp(ndd,6),dm(ndd),g(ndd,npp)
  ! new...
  real*4 pm(np),dp(nd,6),dm(nd),g(nd,np)
  ! JH HACK
  real*4 d(nd)
  real p0,pp360,pm360
  real tstart,tfinish

  ! JGFIXME NEED TO REMOVE remaining npp dependencies etc..
  include 'occamdim.inc'



  ! dehhed stuff:
  integer nlyrsp,nrngsp
  parameter (nlyrsp=npp ,nrngsp=ndd )
  integer lunin,lunout,lvlout,nel,nrange,ifdiffa
  real*8 asigma(0:nlyrsp),athick(0:nlyrsp)
  real*8 aranges(nrngsp),afreq,dphi,zxmit,zrecv,reqerr 
  complex*16 er(nrngsp),ep(nrngsp)
  complex*16 derds(nlyrsp,nrngsp),depds(nlyrsp,nrngsp)
  complex*16 derdt(nlyrsp,nrngsp),depdt(nlyrsp,nrngsp)

  ! local variables
  integer i, j, k, lut(ndd), ityp
  character*10 how
  real pmax, dpmax, sqrt2, phia, rot, rrot
  real ta, arho, aphi, xreal, yreal, ximag, yimag, xamp, yamp
  real ddreal, ddimag, freq
  real rmax, rmin
  parameter (rmax = 7., rmin = -5.)
  logical csemflg
  complex c, dcdp(npp)
  real ac2, realc, aimgc, reald, aimgd

  ! JG additions for multi fwd models
  integer ii,i_g,imin,imax
  real*8 zxmitt

  ! JH additions for phase derivatives
  real*8 yoverx, dxds, dyds

  !     JG test comparison: accuracy tests of multi-range calls to hedsfi()
  !     call fordiv_test(np,nd,pm,dp,dm,g) 

  !write(*,fmt='(A,(1000f15.8))') 'thickf: ',(thickf(i),i=0,nlayer)
  !write(*,fmt='(A,(1000f15.8))') 'rho: ',(rho(i),i=0,nlayer)
  !write(*,fmt='(A,(1000f15.8))') 'ipen1: ',(ipen1(i),i=0,nlayer)
  !write(*,fmt='(A,(1000f15.8))') 'ipen2: ',(ipen2(i),i=0,nlayer)
  !write(*,fmt='(A,(1000i4))') 'iff: ',(iff(i),i=1,np)
  !write(*,fmt='(A,(1000f15.8))') 'pm: ',(pm(i),i=1,np)
  !stop 

  ! set water conductivity and depth      
  if (rho(0) .ne. 0.0) then
     asigma(0) = 1/rho(0)
     csemflg = .true.
  else
     csemflg = .false.
  end if
  athick(0) = thickf(0)
  nel = nlayer
  ! update the free parameters in the model
  do i = 1,np
     rho(iff(i)) = 10.**min(rmax,max(rmin,pm(i)))
  end do

  ! copy across the model into the csem code vectors
  do  i = 1,nlayer-1
     athick(i) = thickf(i)
  end do
  do  i = 1, nlayer
     asigma(i) = 1./rho(i)
  end do

  !  set level of debug-output, logical units for any i/o (csem code)
  lvlout=0
  lunin=5
  lunout=6
  !  set computation method, requested accuracy (csem code)
  !  (the possibilities for how are 'filter', 'quad', 'compare'. 
  !  set the method for the default by using whatever other string
  how = 'filter'
  reqerr = 0.1e-12
  !  ask for forward modelling only (no derivatives) (csem code)
  ifdiffa = 1
  nrange = 1
  ! set azimuth to 45 degrees to get both components of the field
  ! (will be decomposed later for our specific needs) (csem code)
  dphi = 45.
  ! set the receiver height (csem code)
  zrecv = 0.0
  sqrt2 = sqrt(2.)
  !  put model parameters etc. in common:

  call hedfix(nel,asigma,athick,lunin,lunout,lvlout)

    !     JG, replace this with
  !     do i = 1, nd
  i=1
100 if (i.le.nd) then

     ityp = int(dp(i,1)+0.01)

     if (ityp .lt. 20) then
        ! compute CSEM response        

        afreq = dp(i,3)
        aranges(1) = dp(i,2)
        zxmit = max(1.0, dp(i,6))
        imin = i

        if (.not.use_range_lumping) then  

           imax=i   
           !              Compute fields one at a time...
           !  compute the field & desired derivatives at the given coordinates
           !  (several ranges can be done at once)
           call hedsfi(afreq,nrange,aranges,dphi,zxmit,zrecv &
                   ,nlyrsp,nrngsp,how,reqerr,ifdiffa &
                   ,er,ep,derds,derdt,depds,depdt)
        else  
           !            Compute fields in batches of common frequency and Tx height..
           imin = i
           imax = nd
           !            Find a cluster...
           do ii = imin+1, nd 
              ! transmitter height - FIXME! should be dp(ii,6), or suitable cluster average....
              zxmitt = dp(ii,6)
              if ((abs(dp(ii,3) - afreq) > 1.0d-6).or.(abs(max(1.0, zxmitt) - zxmit) > 1.0d-6) &
                   .or.(.not.((int(dp(ii,1)+0.01)) < 20))) then
                 !                 different freq or transmitter height or not CSEM data
                 imax = ii-1
                 goto 500    
              else
                 aranges(ii-imin+1) = dp(ii,2)
                 ! dumb hack for ranges monotonicity....
                 if (.not.(aranges(ii-imin+1).gt.aranges(ii-imin)))then
                    aranges(ii-imin+1) = aranges(ii-imin)*(1.d0+1.0d-15)
                 end if
              end if
           end do
500        continue
           !            Groupable data spans imin..imax
           nrange = imax - imin+1
           !            compute the field & desired derivatives at the given coordinates
           !            (several ranges done at once)

           !call cpu_time(tstart) 
           call hedsfi(afreq,nrange,aranges,dphi,zxmit,zrecv &
                   ,nlyrsp,nrngsp,how,reqerr,ifdiffa &
                   ,er,ep,derds,derdt,depds,depdt)

           !call cpu_time(tfinish)
           !tTot = tTot+(tfinish-tstart) !JGHACK_TIMING
        end if

        !     loop over lumped group
        i_g = 1

        do ii = imin,imax
           phia = dp(ii,4)
           rot = dp(ii,5)
           ityp = int(dp(ii,1)+0.01)

           !     if (ityp .le. 6) then
           !     JGHACK for Re + Im
           if ((ityp .le. 6).or.((ityp.ge.14).and.(ityp.le.17))) then
              !     the azimuth has been fixed at 45 degrees.
              !     now, if the data type is 1-6 we want the erho and ephi terms,
              !     obtained simply be multiplying all the components by sqrt(2):
              er(i_g) = er(i_g)*sqrt2
              ep(i_g) = ep(i_g)*sqrt2
              do j = 1, nel
                 derds(j,i_g) = derds(j,i_g)*sqrt2
                 depds(j,i_g) = depds(j,i_g)*sqrt2
              end do
           else
              !     otherwise we want the erho and ephi components, obtained by
              !     further correcting the fields for azimuthal dependence of dipoles:
              er(i_g) = cos(phia/57.29578)*er(i_g)*sqrt2
              ep(i_g) = sin(phia/57.29578)*ep(i_g)*sqrt2
              do j = 1, nel
                 derds(j,i_g) = cos(phia/57.29578)*derds(j,i_g)*sqrt2
                 depds(j,i_g) = sin(phia/57.29578)*depds(j,i_g)*sqrt2
              end do
           end if

           !     now do the right thing for this data type:
           !     goto (1,2,3,4,5,6,7,8,9,10, 11, 12, 13) ityp
           !     JGHACKS for Re + Im
           goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) ityp
           write(*,*) 'component flag invalid in dp()'
           stop

           !     erho magnitude
1          dm(ii) = abs(er(i_g))
           do j = 1,nel
              g(ii,j) = (real(er(i_g))*real(derds(j,i_g)) + &
                dimag(er(i_g))*dimag(derds(j,i_g)))/dm(ii)
              g(ii,j) = -g(ii,j)/rho(j)/0.43429448190325183
           end do
           goto 310

           !     erho phase
2          dm(ii) = atan2(dimag(er(i_g)),dreal(er(i_g)))*57.29578
           ! JH additions + JG corrections: check pp360, pm360 defs
           p0 = abs( dm(ii) - d(ii) )
           pp360 = abs( dm(ii) + 360.0 - d(ii) )
           pm360 = abs( dm(ii) - 360.0 - d(ii) )
           if( pm360.lt.p0 ) then
             dm(ii) = dm(ii) - 360.0
           elseif ( pp360.lt.p0) then
             dm(ii) = dm(ii) + 360.0
           endif

           ! JH, JG: original SC derivatives here bogus. New versions..
           ! NB. 131.928433 = (180/pi)log10
           yoverx = dimag(er(i_g))/dreal(er(i_g))
           do j = 1, nel
              dyds = dimag(derds(j,i_g))
              dxds = dreal(derds(j,i_g))
              g(ii,j) = (dyds - (dxds*yoverx)) / (dreal(er(i_g))*(1+yoverx**2))
              ! WRONG!-> g(ii,j) = -131.928433 * g(j,i_g) / rho(j)
              g(ii,j) = -131.928433 * g(ii,j) / rho(j)
           end do

           goto 310

           !     ephi magnitude
3          dm(ii) = abs(ep(i_g))
           do j = 1,nel
              g(ii,j) = (real(ep(i_g))*real(depds(j,i_g)) + &
                dimag(ep(i_g))*dimag(depds(j,i_g)))/dm(ii)
              g(ii,j) = -g(ii,j)/rho(j)/0.43429448190325183
           end do
           goto 310

           !     ephi phase
4          dm(ii) = atan2(dimag(ep(i_g)),dreal(ep(i_g)))*57.29578
           p0 = abs( dm(ii) - d(ii) )
           ! JG corrections: check pp360, pm360 defs
           pm360 = abs( dm(ii) - 360.0 - d(ii) )
           pp360 = abs( dm(ii) + 360.0 - d(ii) )
           if( pm360.lt.p0 ) then
             dm(ii) = dm(ii) - 360.0
           elseif ( pp360.lt.p0) then
             dm(ii) = dm(ii) + 360.0
           endif


           ! JH, JG: original derivatives here bogus. New versions..
           yoverx = dimag(ep(i_g))/dreal(ep(i_g))
           do j = 1, nel
              dyds = dimag(depds(j,i_g))
              dxds = dreal(depds(j,i_g))
              g(ii,j) = (dyds - (dxds*yoverx)) / (dreal(ep(i_g))*(1+yoverx**2))
              ! WRONG!->g(ii,j) = -131.928433 * g(j,i_g) / rho(j)
              g(ii,j) = -131.928433 * g(ii,j)/ rho(j)
           end do
           goto 310

           !     erho log10 magnitude
5          dm(ii) = abs(er(i_g))
           do j = 1,nel
              g(ii,j) = (real(er(i_g))*real(derds(j,i_g)) + &
                dimag(er(i_g))*dimag(derds(j,i_g)))/dm(ii)
              g(ii,j) = -g(ii,j)/rho(j)/dm(ii)
           end do
           dm(ii) = log10(dm(ii))
           goto 310

           !     ephi log10 magnitude
6          dm(ii) = abs(ep(i_g))
           do j = 1,nel
              g(ii,j) = (real(ep(i_g))*real(depds(j,i_g)) + &
                dimag(ep(i_g))*dimag(depds(j,i_g)))/dm(ii)
              g(ii,j) = -g(ii,j)/rho(j)/dm(ii)
           end do
           dm(ii) = log10(dm(ii))
           goto 310

           !     log of major axis of pe
7          do j = 1, nel
              call pederiv(er(i_g),ep(i_g),derds(j,i_g), &
                depds(j,i_g),pmax,dpmax)
              g(ii,j) = -dpmax/rho(j)/pmax
           end do
           dm(ii) = log10(pmax)
           goto 310

           !     linear major axis of pe
8          do j = 1, nel
              call pederiv(er(i_g),ep(i_g),derds(j,i_g), &
                depds(j,i_g),pmax,dpmax)
              g(ii,j) = -dpmax/rho(j)/0.43429448190325183
           end do
           dm(ii) = pmax
           goto 310

           !     log of erho component
9          dm(ii) = abs(er(i_g))
           do  j = 1,nel
              g(ii,j) = (real(er(i_g))*real(derds(j,i_g)) + &
                dimag(er(i_g))*dimag(derds(j,i_g)))/dm(ii)
              g(ii,j) = -g(ii,j)/rho(j)/dm(ii)
           end do
           dm(ii) = log10(dm(ii))
           goto 310

           !     log of ephi component
10         dm(ii) = abs(ep(i_g))
           do j = 1,nel
              g(ii,j) = (real(ep(i_g))*real(depds(j,i_g)) + &
                dimag(ep(i_g))*dimag(depds(j,i_g)))/dm(ii)
              g(ii,j) = -g(ii,j)/rho(j)/dm(ii)
           end do
           dm(ii) = log10(dm(ii))
           goto 310

           !     log of x-component
11         continue
           ta = tangle - 180.
           arho = (ta + phia + 90.)/57.29578
           aphi = (ta + phia + 180.)/57.29578
           xreal = cos(arho)*real(er(i_g)) + cos(aphi)*real(ep(i_g))
           ximag = cos(arho)*dimag(er(i_g))+cos(aphi)*dimag(ep(i_g))
           xamp = sqrt(xreal**2 + ximag**2)
           !     xphz = atan2(ximag,xreal)*57.29578
           dm(ii) = log10(xamp)
           do  j = 1,nel
              ddreal=cos(arho)*real(derds(j,i_g)) &
                +cos(aphi)*real(depds(j,i_g))
              ddimag=cos(arho)*dimag(derds(j,i_g)) &
                +cos(aphi)*dimag(depds(j,i_g))
              g(ii,j)=(xreal*ddreal + ximag*ddimag)/xamp
              g(ii,j) = -g(ii,j)/rho(j)/xamp
           end do
           goto 310

           !     log of y-component
12         continue
           ta = tangle - 180.
           arho = (ta + phia + 90.)/57.29578
           aphi = (ta + phia + 180.)/57.29578
           yreal = sin(arho)*real(er(i_g)) + sin(aphi)*real(ep(i_g))
           yimag = sin(arho)*dimag(er(i_g))+ sin(aphi)*dimag(ep(i_g))
           yamp = sqrt(yreal**2 + yimag**2)
           !     yphz = atan2(yimag,yreal)*57.29578
           dm(ii) = log10(yamp)
           do j = 1,nel
              ddreal=sin(arho)*real(derds(j,i_g)) &
                +sin(aphi)*real(depds(j,i_g))
              ddimag=sin(arho)*dimag(derds(j,i_g)) &
                +sin(aphi)*dimag(depds(j,i_g))
              g(ii,j) = (yreal*ddreal + yimag*ddimag)/yamp
              g(ii,j) = -g(ii,j)/rho(j)/yamp
           end do
           goto 310

13         continue
           !     log of single component data
           !     rot = receiver rotation - angle +ve anti-clocwise to receiver dipole
           !     from erho
           rrot = rot / 57.29578
           yreal = cos(rrot)*dreal(er(i_g))+ sin(rrot)*dreal(ep(i_g))
           yimag = cos(rrot)*dimag(er(i_g))+ sin(rrot)*dimag(ep(i_g))
           yamp = sqrt(yreal**2 + yimag**2)
           dm(ii) = log10(yamp)
           do  j=1,nel
              ddreal=cos(rrot)*dreal(derds(j,i_g)) +  &
                sin(rrot)*dreal(depds(j,i_g))
              ddimag=cos(rrot)*dimag(derds(j,i_g)) +  &
                sin(rrot)*dimag(depds(j,i_g))
              g(ii,j) = (yreal*ddreal + yimag*ddimag)/yamp
              g(ii,j) = -g(ii,j)/rho(j)/yamp
           end do
           goto 310

           !     JGHACKS for Re + Im
           !     erho Re
14         continue
           dm(ii) = dreal(er(i_g))
           do j = 1,nel
              g(ii,j) = dreal(derds(j,i_g))
              !     Steve's original line has something like 
              !     g(i,j) = -g(i,j)/rho(j)/0.43429448190325183
              !     Here's my suggestion....(we want for some Q, 
              !     dQ/dlog_10(rho) = (dQ/dsigma).(-log(10)/rho)
              !     NB log(10) = 2.3025850929940457 = 1/0.43429448190325183
              g(ii,j) = -2.3025850929940457 * g(ii,j)/rho(j)
           end do
           goto 310

           !     erho Im   Same derivation as above
15         continue
           dm(ii) = dimag(er(i_g))
           do j = 1, nel
              g(ii,j) = dimag(derds(j,i_g))
              g(ii,j) = -2.3025850929940457 * g(ii,j)/rho(j)
           end do
           goto 310

           !     ephi Re.   Same derivation as above
16         continue
           dm(ii) = dreal(ep(i_g))
           do j = 1,nel
              g(ii,j) = dreal(depds(j,i_g))
              g(ii,j) = -2.3025850929940457 * g(ii,j)/rho(j)
           end do
           goto 310

           !     ephi Im.   Same derivation as above
17         continue
           dm(ii) = dimag(ep(i_g))
           do j = 1, nel
              g(ii,j) = dimag(depds(j,i_g))
              g(ii,j) = -2.3025850929940457 * g(ii,j)/rho(j)
           end do
           goto 310
           !     ........eo JGHACKS



310        continue
           i_g = i_g+1
        end do               ! loop over group ii=imin...imax

     else if (ityp .le. 30) then
        ! MT case
        freq = 1./dp(i,2)
        call mtdiv(freq,c,dcdp,nlayer,rho,thickf)
        ac2 = abs(c)**2
        realc = real(c)
        aimgc = aimag(c)
        goto (21, 22, 23, 24, 25, 26, 27) ityp-20
        write(*,*) 'component flag invalid in dp()'
        stop
21      dm(i) = log10(7.8957e-06*(abs(c)**2)*freq)
        do  j = 1, np
           reald = real(dcdp(j))
           aimgd = aimag(dcdp(j))
           g(i,j) = 2.*rho(j)*(realc*reald + aimgc*aimgd)/ac2
        end do
        goto 1000
22      dm(i) = 57.29578*atan2(real(c),-aimag(c))
        do  j = 1, np
           reald = real(dcdp(j))
           aimgd = aimag(dcdp(j))
           g(i,j) = 131.9284*rho(j)*(realc*aimgd - aimgc*reald)/ac2
        end do
        goto 1000
23      dm(i) = realc
        do  j = 1, np
           reald = real(dcdp(j))
           aimgd = aimag(dcdp(j))
           g(i,j) = rho(j)*reald/0.43429
        end do
        goto 1000
24      dm(i) = -aimgc
        do  j = 1, np
           reald = real(dcdp(j))
           aimgd = aimag(dcdp(j))
           g(i,j) = -rho(j)*aimgd/0.43429
        end do
        goto 1000
25      dm(i) = log10(realc)
        do  j = 1, np
           reald = real(dcdp(j))
           aimgd = aimag(dcdp(j))
           g(i,j) = rho(j)*reald/realc
        end do
        goto 1000
26      dm(i) = log10(-aimgc)
        do  j = 1, np
           reald = real(dcdp(j))
           aimgd = aimag(dcdp(j))
           g(i,j) = rho(j)*aimgd/aimgc
        end do
        goto 1000
27      dm(i) = 7.8957e-06*(abs(c)**2)*freq
        do  j = 1, np
           reald = real(dcdp(j))
           aimgd = aimag(dcdp(j))
           g(i,j) = 3.63610e-05*rho(j)*(realc*reald+aimgc*aimgd)*freq
        end do
1000    continue

     else
        write(*,*) 'component flag invalid in dp()'
        stop
     end if

     if (imax.gt.imin) i=imax ! makes next increment good   
     i = i+1
     goto 100  
  endif
  !     end do
  return
end subroutine fordiv


subroutine showTotalTiming()
  use fmodel
  implicit none

  ! write(*,*) 'formod time (S) = ',1.0*tTot 
end subroutine


!-----------------------------------------------------------------------
subroutine pederiv(e1, e2, de1, de2, pmax, dpmax)
  !-----------------------------------------------------------------------
  implicit none
  complex*16 e1, e2, de1, de2
  real pmax, dpmax
  real*8 pi2, x1, x2, y1, y2, dx1, dx2, dy1, dy2
  real*8 a, b, phi, s, c, p1, d1, p2, d2

  pi2 = 3.141592654/2.

  ! bust up the complex components of the electric field:
  x1 = real(e1)
  y1 = dimag(e1)
  x2 = real(e2)
  y2 = dimag(e2)
  dx1 = real(de1)
  dy1 = dimag(de1)
  dx2 = real(de2)
  dy2 = dimag(de2)

  ! find the critical angle
  a = -(x1*y1 + x2*y2)
  b = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
  phi = atan2(a,b)/2.

  ! find the two critical axes of the pe and derivatives
  s = sin(phi)
  c = cos(phi)
  p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
  p1 = sqrt(abs(p1))
  d1 = c*c*(x1*dx1+x2*dx2) &
   + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
  d1 = d1/p1

  s = sin(phi+pi2)
  c = cos(phi+pi2)
  p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
  p2 = sqrt(abs(p2))
  d2 = c*c*(x1*dx1+x2*dx2) &
   + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
  d2 = d2/p2

  if (p1 .ge. p2) then
     pmax = p1
     dpmax = d1
  else
     pmax = p2
     dpmax = d2
  end if

  return
end subroutine pederiv


!-----------------------------------------------------------------------
subroutine pemax(e1, e2, pmax)
  !-----------------------------------------------------------------------
  implicit none
  complex*16 e1, e2
  real pmax
  real*8 pi2, x1, x2, y1, y2
  real*8 a, b, phi, s, c, p1, p2

  pi2 = 3.141592654/2.

  ! bust up the complex components of the electric field:
  x1 = real(e1)
  y1 = dimag(e1)
  x2 = real(e2)
  y2 = dimag(e2)

  ! find the critical angle
  a = -(x1*y1 + x2*y2)
  b = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
  phi = atan2(a,b)/2.

  ! find the two critical axes of the pe 
  s = sin(phi)
  c = cos(phi)
  p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
  p1 = sqrt(abs(p1))

  s = sin(phi+pi2)
  c = cos(phi+pi2)
  p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
  p2 = sqrt(abs(p2))

  if (p1 .ge. p2) then
     pmax = p1
  else
     pmax = p2
  end if

  return
end subroutine pemax

!
!-----------------------------------------------------------------------
function unwrap(icom, angle)
  !-----------------------------------------------------------------------

  ! function to unwrap angles
  ! on input:
  ! icom = 1: unwrap to -180 to + 180
  !      = 2: unwrap to 0 to 360
  !      = 3: fold into first and second quadrant (0 to 180)
  ! angle = angle to unwrap
  !
  implicit none
  real unwrap, angle
  integer icom

  unwrap = angle

  if (icom .eq. 1) then
5    if (unwrap .lt. -180.) then
        unwrap = unwrap + 360.
        goto 5
     end if
10   if (unwrap .gt. 180.) then
        unwrap = unwrap - 360.
        goto 10
     end if
  else
20   if (unwrap .lt. 0.) then
        unwrap = unwrap + 360.
        goto 20
     end if
30   if (unwrap .gt. 360.) then
        unwrap = unwrap - 360.
        goto 30
     end if
  end if
  if (icom .eq. 3) then
     if (unwrap .gt. 180.) then
        unwrap = unwrap - 180.
     end if
  end if

  return
end function unwrap

!-----------------------------------------------------------------------
subroutine mtlay(f,carg,nlayer,rho,thick)
  !-----------------------------------------------------------------------
  ! occam1dmt: steven constable 10 april 1992
  !
  ! mtlay computes the complex c value over a layered earth
  ! s.c. constable, s.i.o. la jolla ca92093, august 1986
  !
  ! calls: no other routines
  ! on input:
  !   f = frequency in hertz
  !   nlayer = number of layers in model, including terminating half space
  !   rho = array of layer resistivities, ohm-m (linear domain)
  !   thick = array of layer thicknesses, m (linear domain)
  ! on output:
  !   c = complex c value, m (linear domain)
  !
  include 'occamdim.inc'
  dimension thick(0:nlayer-1),rho(0:nlayer)
  complex carg
  ! local variables are double precision.  this is necessary.
  complex*16 ck,omi,c
  dimension wrho(npp), depth(npp)

  omega = 6.283185308*f
  omi = dcmplx(omega)*1.256637062e-06*(0.,1.)

  ! depth(i) is the depth to the top of the i'th layer:
  depth(1) = 0.
  do i = 1,nlayer-1
     depth(i+1) = depth(i) + thick(i)
  end do
  ! map real (spherical earth) depths to flat earth equivalents using weidelt's
  ! algorithm
  re = 6371000.
  wrho(1) = rho(1)
  nl2 = nlayer
  do i = 2, nlayer
     if (depth(i) .ge. re) then
        ! don't consider layers deeper than the center of earth:
        nl2 = i-1
        goto 50
     else 
        r = (re - depth(i))/re
        depth(i) = re*(1./r - r**2 )/( 2./r + r**2)
        fm = (2./r + r**2)/3.
        wrho(i) = rho(i)/(fm**4)
     end if
  end do
50 continue

  c = 1./sqrt(omi/dcmplx(wrho(nl2)))
  do i = nl2-1,1,-1
     t = depth(i+1) - depth(i)
     ck = sqrt(omi/dcmplx(wrho(i)))
     c = exp(-2.*ck*dcmplx(t))*(ck*c-1.)/(ck*c+1.)
     c = (1.+c)/(1.-c)/ck
  end do
  carg = c
  return
end subroutine mtlay
!
!
!
!-----------------------------------------------------------------------
subroutine mtdiv(freq,carg,darg,nlayer,rho,thick)
  !-----------------------------------------------------------------------
  ! occam1dmt: steven constable 10 april 1992
  !
  ! mtdiv computes the complex c value over a layered earth and dc/drho for
  !   each layer (i.e. a row of the jacobian matrix)
  ! s.c. constable, s.i.o. la jolla ca92093, august 1986
  !
  ! calls: no other routines
  ! on input:
  !   freq = frequency in hertz
  !   nlayer = number of layers in model, including terminating half space
  !   rho = array of layer resistivities, ohm-m (linear domain)
  !   thick = array of layer thicknesses, m (linear domain)
  ! on output:
  !   c = complex c value, m (linear domain)
  !   d = array of complex derivatives of c w.r.t. layer resistivities, 
  !              m/ohm-m
  !
  include 'occamdim.inc'
  dimension thick(0:nlayer-1),rho(0:nlayer)
  complex darg(nlayer),carg
  ! local variables are double precision.  this is necessary.
  complex*16 d(npp),c
  complex*16 omi,ckm1,ckp1,opexp,omexp,ck,dciip1,theexp,ri
  dimension wrho(npp), depth(npp)
  !
  omega = 6.283185308*freq
  omi = dcmplx(omega)*1.256637062e-06*(0.,1.)

  ! depth(i) is the depth to the top of the i'th layer:
  depth(1) = 0.
  do i = 1,nlayer-1
     depth(i+1) = depth(i) + thick(i)
  end do
  ! map real (spherical earth) depths to flat earth equivalents using weidelt's
  ! algorithm
  re = 6371000.
  wrho(1) = rho(1)
  nl2 = nlayer
  do i = 2, nlayer
     if (depth(i) .ge. re) then
        ! don't consider layers deeper than the center of earth:
        nl2 = i-1
        goto 50
     else 
        r = (re - depth(i))/re
        depth(i) = re*(1./r - r**2 )/( 2./r + r**2)
        fm = (2./r + r**2)/3.
        wrho(i) = rho(i)/(fm**4)
     end if
  end do
50 continue

  ! first of all dc/drho for i=1,2,..nlayer will be stored in d()
  ! then propagated to the surface using dc(i)/dc(i-1) at each layer.
  !
  ! establish the values for the top of the terminating half-space
  ck = sqrt(omi/dcmplx(wrho(nl2)))
  c = 1./ck
  d(nl2) = 1./(ck*dcmplx(wrho(nl2))*2.)
  ! now sweep up through the layers
  do  i = nl2-1,1,-1
     t = depth(i+1) - depth(i)
     ck = sqrt(omi/dcmplx(wrho(i)))
     ckp1 = c*ck + 1.
     ckm1 = c*ck - 1.
     ri = ckm1/ckp1
     theexp = exp(-2.*ck*dcmplx(t))
     opexp = theexp*ri
     omexp = 1. - opexp
     opexp = 1. + opexp
     ! compute dc/drho for the (i)th layer
     d(i) = 2.*theexp*(c/ckp1**2 - dcmplx(t)*ri)/omexp**2
     d(i) = opexp/omexp/ck/2. - d(i)
     d(i) = d(i)/dcmplx(wrho(i))
     ! compute dc(i)/dc(i+1) for the (i)th layer
     dciip1 = 4.*theexp/(ckp1*omexp)**2
     ! now compute the value of c at the top of this layer
     c = opexp/omexp/ck
     ! multiply the product sums for all the layers below this one, which have
     !   already been partially accumulated.
     do j = i+1,nl2
        d(j) = d(j)*dciip1
     end do
  end do
  ! copy double complex results into single complex arguments.
  carg = c
  do i=1,nlayer
     if (i .le. nl2) then
        darg(i) = d(i)
     else
        darg(i) = 0
     end if
  end do
  ! and that's it
  return
end subroutine mtdiv


