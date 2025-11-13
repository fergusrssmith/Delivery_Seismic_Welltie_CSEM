c    This file is part of LegacyScrippsDipole1D
c    LegacyScrippsDipole1D is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    LegacyScrippsDipole1D is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    
c    You should have received a copy of the GNU General Public License
c    along with LegacyScrippsDipole1D.  If not, see <http://www.gnu.org/licenses/>.



c    ----------------------------------------------------------------------------------
c    Original code:
c    Agusta Flosadottir, 1993

c    These are the unmodified routines - will crash if you call fdrlpm
c    on the mac with derivatives set.
c                    hedfix
c
c  this routine must be called before hedsfi
c  (which computes seafloor horizontal electric field values).
c  hedfix feeds msc. "fixed" parameters into common blocks
c  in file hedcom.inc. 
c  these are later be used by hedsfi and its subroutines.
c  the parameters filled in here are: 
c    io logical units and debug-output level, 
c    source length and strength,
c    the layered seafloor structure (scaled for chave's kernels)
c  still to be set up for each field-evaluation (hedsfi) call are:
c    frequency,range,azimuth in dehhed's units and precision,
c    source and receiver heights,
c    level of parameter-derivatives to be computed
c    msc. control parameters for the numerical evaluation.
c  all arguments are in mks units.
c
       subroutine hedfix(nel,asigma,athick,lunin,lunout,lvlout)
c
c  agusta flosadottir, 1993
c  last modified june 93
c
c
c incoming arguments:
c      nel      no. earth layers
c                (counting deep halfspace, but not sea)
c      asigma(i), i = 0,nel
c        layer conductivities (starting with sea) in s/m
c        layers are numbered from the top (sea, i=0) down. 
c        final conductivity corresponds to the deep halfspace. 
c      athick(i), i = 0,nel-1
c        layer thickness in m
c        no thickness value for the deep halfspace. 
c      lunin,lunout, logical units for debug,error-message i/o. 
c      lvldbg debug output level. 
c           0  for no output
c        >  5  for output from this routine (hedfix)
c        > 10  for output from dehhed
c        >100  for output from hedmods
c
c  called from:            user supplied
c  subprograms required:   none
c  include files required: hedcom.inc
c                          heddat.inc
c
c omit double precision declarations and revise exponents in heddat.inc
c if using the cray version of dehhed
c
c
c--internal
c   common blocks for communication with & among subroutines
      include 'hedcom.inc'
c--arguments
c   double precision arguments
      integer nel 
      real*8 asigma(0:nel),athick(0:nel)
c
c   msc. numbers, math & phys constants
      include 'heddat.inc'
c
c----copy io logical units and debug output level into common, 
c    (used for error-messages and debugging output)
      ni=lunin
      no=lunout
      iw=lvlout
c----set local threshold debug output level:
      iw0=5
      namepu='hedfix'
       if (iw.gt.iw0) then
         write (no,*) ' ',namepu,' entering:'
         write (no,*) '    nel=',nel
         write (no,*) '    how=',how
       end if
c---make sure arrays are big enough:
      if (nel.gt.npl1) then
        write (no,*) nel,' earth layers requested,'
        write (no,*) ' max dimnsn in hedcom.inc, npl1=',npl1
        msg=' too many layers'
        go to 911
      end if
c
c---copy msc arguments into internal (double) precision variables
c   and scale as called for by dehhed:
c
c    argument no. earth layers => no. finite layers:
      nl=nel-1
      nl1=nel
c   scale earth structure
      do 30 i = 0,nl1
        sigma(i) = asigma(i)*tent9
30    continue
      h=athick(0)*tentm3
      do 35 i = 1,nl
        thick(i) = athick(i)*tentm3
35    continue
c
c---dipole transmitter selected, extended source not implemented
      sl = zero
c---transmitter dipole moment set to unity
      sc = one
c
       if (iw.gt.iw0) then
           write (no,*) ' ',namepu,': now stored in common hed'
           write (no,*) '   nppar=',nppar,',npdeh=',npdeh
           write (no,*) '   nl=',nl,',nl1=',nl1
           write (no,*) '   sigma(0)=',sigma(0)
           write (no,*) '   sigma: ',(sigma(n),n=0,nl1)
           write (no,*) '   thick: ',(thick(n),n=1,nl)
           write (no,*) '   h=',h,',zp=',zp,',zm=',zm
           write (no,*) '   sl=',sl,',sc=',sc
       end if
c
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c
c                    hedsfi
c
c  computes seafloor horizontal electric field
c  for a unit horizontal electric dipole over a layered seafloor.
c  also (optionally) computes derivatives with respect to
c  model conductivities and thicknesses.
c  takes a vector of ranges so advantage is taken of w.l.anderson's
c  scheme for lagged convolutions and spline interpolation.
c  the ranges must be in ascending order or an error condition results.
c
c  hedfix must have been called before this routine !
c  (in order to set up most of /deh/ and /inoutc/ in hedcom.inc).
c
c
       subroutine hedsfi(afreq,nrange,aranges,aphi,zxmit,zrecv
     $     ,nlyrsp,nrngsp
     $     ,how,reqerr,ifdiffa
     $     ,er,ep,derds,derdt,depds,depdt)
c
c  agusta flosadottir, 1993
c  last modified june 93
c
c     this subroutine 
c     (horizontal electric dipole seafloor/sommerfeld integral)
c
c     1. sets up msc. numerical parameters,
c     fills in some variables in hedcom.inc /hed/
c     (others have already been filled in by routine hedfix).
c     uses variables from /hed/ too (nl1 etc.)
c     all arguments are in mks units.
c
c     2. calls dehhed, my derivatives-added version
c     of a.d.chave's sommerfeld integral kernels (modified ehhed)
c     with w.l.anderson/chave hybrid digital filter or quadrature 
c     hankel transform (hy2fht).
c
c incoming arguments:
c      afreq    frequency in hz, assuming exp(-i omega t)
c      nrange   actual number of ranges
c      aranges  source <-> fieldpoint distance in m
c      aphi     azimuth in degrees
c      zxmit    transmitter height above seafloor in m
c      zrecv    receiver height above seafloor in m
c      nlyrsp   dimension for layers
c      nrngsp   dimension for ranges
c      how 
c        'filter'     to enforce digital filter integration
c        'quadrature' to enforce quadrature integration
c        'compare'    to compare the above
c        'ahed'  not-implmtd old version, no-derivatives filter integration
c        default      hybrid algorithm (quadrature if filter fails)
c      reqerr   requested error, try .1d-12 ?
c      ifdiffa  requested level of derivative computation
c        0:  no derivatives, field value only
c        1:  field value & conductivity derivatives
c        2:  field value, conductivity & thickness derivatives
c outgoing arguments:
c      er      radial electric field, various ranges.
c      ep      azimuthal electric field
c      derds   derho/dsigma, 1.column: successive layers, 2.: ranges
c      derdt   derho/dthickness,    --"--
c      depds   dephi/dsigma,        --"--
c      depdt   dehi/dthickness,     --"--
c
c
c  called from:            user supplied
c  subprograms required:   hy2fht-w.l. anderson's hybrid hankel-transform,
c                               using either a fast filter transform or
c                               a.d. chave's algorithm,
c                               double precision version.
c                        hedmods-collection of subprogram modules
c                        bessj-bessel functions (can use a dummy routine
c                        if only plannning to use filter integration)
c  include files required: hedcom.inc
c                          heddat.inc
c
c
c omit double precision declarations and revise exponents in heddat.inc
c if using the cray version of dehhed
c
c
c   common block for communication with & among subroutines
c   (and with hedfix, assumed to have been called already)
      include 'hedcom.inc'
c
c--arguments
c   double precision arguments
      integer nel 
      real*8 aranges(nrngsp),afreq,aphi,rerr,aerr 
      complex*16 er(nrngsp),ep(nrngsp)
      complex*16 derds(nlyrsp,nrngsp),depds(nlyrsp,nrngsp)
      complex*16 derdt(nlyrsp,nrngsp),depdt(nlyrsp,nrngsp)
      character*(*) how
c
c--internal
c   ranges for dehhed
      dimension rangkm(maxrng)
c   1d arrays receiving total field and derivatives from dehhed:
      complex*16 derho(npdeh,maxrng),dephi(npdeh,maxrng)
c   polarizations returned separately:
      complex*16 derpm(npdeh,maxrng),deppm(npdeh,maxrng)
      complex*16 dertm(npdeh,maxrng),deptm(npdeh,maxrng)
c   error flag from dehhed
      dimension ierrd(2)
c
c   msc. numbers, math & phys constants
      include 'heddat.inc'
c
c----copy debug output level into common, set local threshold debug output level:
      iw0=5
      namepu='hedsfi'
       if (iw.gt.iw0) then
         write (no,*) ' ',namepu,' entering:'
         write (no,*) '    arange,aphi,afreq=',arange,aphi,afreq
         write (no,*) '    zxmit,zrecv=',zxmit,zrecv
         write (no,*) '    nel,ifdiffa=',nel,ifdiffa
         write (no,*) '    how=',how
         write (no,*) '    reqerr=',reqerr
       end if
c---make sure arrays are big enough:
c   (nl1 and npl1 both from hedcom.inc /deh/)
      if (nl1.gt.npl1 .or. nl1.gt.nlyrsp) then
        write (no,*) nl1,' earth layers requested,'
        write (no,*) ' max dimnsn in hedcom.inc, npl1=',npl1
        write (no,*) ' max argumt dimnsn, nlyrsp=',nlyrsp
        msg='too many layers'
        go to 911
       else if (nl1.ge.1) then
        continue
       else
        write (no,*) ' nl1 =',nl1
        msg=' has hedfix been called ?'
        go to 911
      end if
      if (nrange.gt.nrngsp .or. nrange.gt.maxrng) then
        write (no,*) nrange,' ranges requested,'
        write (no,*) ' max argumt dimnsn, nrngsp=',nrngsp
        write (no,*) ' max local dimnsn, maxrng=',maxrng
        msg=' too many ranges at once'
        go to 911
      end if
c
c---copy azimuth and differentiation level into hedcom.inc /deh/: 
      phi=aphi
      ifdiff=ifdiffa
c---xmitter and receiver heights above the seafloor, zp,zm into /deh/:
      zp = abs(zxmit+zrecv)*tentm3
      zm = abs(zxmit-zrecv)*tentm3
      if (zp.lt.tentm3) then
        write (no,*) '  transmitter height over seafloor=',zxmit,' m.'
        write (no,*) ' and receiver height over seafloor=',zrecv,' m.'
        write (no,*) 
     $  ' non-zero zp=abs(zxmit+zrecv) better for num. stability'
        msg= 
     $' try receiver on bottom, transmitter 1 m off the bottom'
        go to 911
      end if
c---convert frequency, ranges to dehhed's precision, units:
      freq=afreq
      do nr=1,nrange
        rangkm(nr)=aranges(nr)*tentm3
      enddo
c
c---decide whether to analytically integrate part of tm filter xform:
c   (previously separate ifs in dehhed and integrand)
      if (zp.eq.0.0) then
        ianatm=1
        write (no,*) '**** warning ',namepu,':'
        write (no,*) 
     $ '    zp=0. partly analytical integration not tested in dehhed.'
        write (no,*) ' zp > 0 is better for numerical stability anyway'
        msg=' try transmitter 1 m off seafloor'
        go to 911
       else
        ianatm=-1
      end if
c
c---control parameters for hankel transform routine hy2fht:
      if ( (index(how,'filter').ne.0) 
     $.or. (index(how,'filter').ne.0)) then
        ipath=1
       else if ( (index(how,'quad').ne.0) 
     $      .or. (index(how,'quad').ne.0)
     $      .or. (index(how,'quad').ne.0)) then
        ipath=2
       else if ( (index(how,'compare').ne.0) 
     $      .or. (index(how,'compare').ne.0)
     $      .or. (index(how,'compare').ne.0)) then
        ipath=-2
       else 
        ipath=-1
      end if
c   digital filter transform parameters:
c    tol = approximately 0.01 * requested relative error
c    anderson mentions ntol= 2 or 3 for oscillating kernels.
        tol=0.01*reqerr
        ntol=2
c   quadrature parameters:
c    rerr,aerr are relative,absolute requested errors
c    npcs is no.pieces for partial integration, 
c    try > 1 for difficult kernels ?
        rerr=reqerr
        aerr=0.01*reqerr
        npcs=1
c
c
       if (iw.gt.iw0) then
           write (no,*) ' ',namepu,': before call dehhed, arguments'
           write (no,*) '    rangkm,freq=',rangkm,freq
           write (no,*) '    tol,ntol,rerr,aerr=',tol,ntol,rerr,aerr
           write (no,*) '    npcs,ipath=',npcs,ipath
           write (no,*) ' in common:'
           write (no,*) '   nppar=',nppar,',npdeh=',npdeh
           write (no,*) ' hed:'
           write (no,*) '   nl=',nl,',nl1=',nl1
           write (no,*) '   sigma(0)=',sigma(0)
           write (no,*) '   sigma: ',(sigma(n),n=0,nl1)
           write (no,*) '   thick: ',(thick(n),n=1,nl)
           write (no,*) '   h=',h,',zp=',zp,',zm=',zm
           write (no,*) '   sl=',sl,',sc=',sc,',phi=',phi
           write (no,*) '   ifdiff=',ifdiff
       end if
c
c---compute the field and its derivatives
         call dehhed(rangkm,nrange,nrngsp,freq,derho,dephi
     $            ,ierrd,derpm,deppm,dertm,deptm
     $            ,tol,ntol,rerr,aerr,npcs,ipath)
         if (iw.gt.iw0) then
           write (no,*) ' ',namepu,': back from dehhed'
         end if
c
c---copy results into multi-d arrays, argument precision. 
c   separate arrays for field, sigma-derivatives, thick-derivatives:
      do nr=1,nrange
        er(nr)=derho(1,nr)
        ep(nr)=dephi(1,nr)
        do 100 l=1,nl1
          j=l+1
          derds(l,nr)=derho(j,nr)
          depds(l,nr)=dephi(j,nr)
100     continue
        do 200 l=1,nl
          j=1+nl1+l
          derdt(l,nr)=derho(j,nr)
          depdt(l,nr)=dephi(j,nr)
200     continue
c  zero the derivative with respect to the deep halfspace thickness
c  to allow double loops over ist and l in calling routines
c  without crash due to undefined numbers:
        derdt(nl1,nr)=zero
        depdt(nl1,nr)=zero
c
        if (iw.gt.iw0) then
          write (no,*) '   ',namepu,' debugging:'
          write (no,*) '    er(',nr,')=',er(nr)
          write (no,*) '    ep(',nr,')=',ep(nr)
        end if
      enddo
c
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c
c                    dehhed.f77    
      subroutine dehhed(ranges,nrange,nrngsp,freq,derho,dephi
     $                 ,ierrd,derpm,deppm,dertm,deptm
     $                 ,tol,ntol,rerr,aerr,npcs,ipath)
c
c---- to be called by hedsfi, which fills in necessary common block /hed/. ----
c
c  computes horizontal electric field for an hed source.
c  does the computation for a vector of ranges.
c  see comments in hedsfi
c
c  input parameters:
c
c  ranges- horizontal source-receiver ranges
c  freq- source frequency
c  tol,ntol,rerr,aerr,ipath-controls for hankel transform, see routine hy2fht
c                           defaults supplied 
c  remaining parameters are assumed to have been placed in common /hed/
c  passing them to all subprograms -- 
c
c  output parameters:
c
c  derho-derho(1)=erho, radial electric field 
c       -derho(m)=d(erho )/d(sigma(m)), m=2,nl+2=jl
c       -derho(m)=d(erho )/d(thick(m)), m=jl+1,jl+nl=jmax
c  dephi-azimuthal electric field, otherwise as above
c  ierrd-vector containing an error flag, 0 for normal return
c  derpm-radial field and derivatives, pm contribution
c  dertm-radial field and derivatives, tm contribution
c  deppm-azimuthal field and derivatives, pm contribution
c  deptm-azimuthal field and derivatives, tm contribution
c
c
c  note: all lengths must be in km, time in seconds, mass in kg, but
c        subroutine returns electric field per unit source moment
c        in v/(amp-m**2)
c        and derivatives in mks (v/(amp-m**2))/(s/m) and (v/(amp-m**2))/m
c
c  called from:          hedsfi
c  subroutines required: hy2fht-w.l. anderson's hybrid hankel-transform,
c                               using either a fast filter transform or
c                               a.d. chave's algorithm,
c                               double precision version.
c                        hedmods-collection of subprogram modules
c                        bessj-bessel function (can use a dummy routine
c                        if only plannning to use filter integration)
c  include files required: hedcom.inc
c                          heddat.inc
c
c  common /hed/ parameter list
c
c  sigma-vector of layer conductivities, sigma(0) is source region
c        conductivity.  these are the actual conductivity in s/m
c        scaled (multiplied) by 1.e9
c  thick-vector of layer thicknesses in km
c  h-thickness of source region in km
c  zp-sum of source and receiver heights in km, where both the
c     source and receiver heights must be .ge. 0
c  zm-absolute value of difference of source and receiver heights in
c     km
c  sl-source length in km, 0 for point source
c  sc-source current in amps, used only when sl.ne.0
c  phi-azimuth in degrees
c  nl-number of layers, excluding sea and terminating halfspace
c  nl1-  nl1=nl+1, including the halfspace, but not the sea
c  ifdiff-0 to compute only electric field
c         1 to compute field and d(field)/d(sigma(m)),m=1,2,..,nl+1
c         2 to compute field and d(field)/d(thick(m)),m=1,2,..,nl
c  phi-azimuthal angle, positive ccw from the x-axis toward the y-axis
c  ---------------------------------------------------------
      include 'hedcom.inc'
c  arrays passed to the hankel-transform hy2fht
      integer nord(maxrel),nfun1(2)
      real*8 b(maxrng),bwork(maxrng,9)
      complex*16 zwork(801,maxrel),zans(maxrng,maxrel),fwork(maxrel)
c  input arrays
      dimension ranges(nrngsp)
c  arrays returning the answer, total and polarizations separately:
      complex*16 derho(npdeh,nrngsp),dephi(npdeh,nrngsp)
      dimension ierrd(2)
      complex*16 derpm(npdeh,nrngsp),deppm(npdeh,nrngsp)
      complex*16 dertm(npdeh,nrngsp),deptm(npdeh,nrngsp)
      complex*16 derww,depww
c  subprograms used as integrands, passed to hy2fht
      external fun1pm,fun1tm,rel1pm,rel1tm
      complex*16 fun1pm,fun1tm 
c  functions returning analytical integrals
      complex*16 ehsr1,ehsr2,ehsr3
c  local scalars
      complex*16 t0,t1,t2,dej1,dej0
c  for sample integrand calls
      complex*16 f1
c
      include 'heddat.inc'
c
      iw0=10
      namepu='dehhed'
      if (iw.gt.iw0) then
        write (no,*) ' entering dehhed. arguments:'
        write (no,*) '   ranges(1)=',ranges(1),' km, freq=',freq,' hz'
        write (no,*) '   ranges(',nrange,')=',ranges(nrange)
        write (no,*) '   tol,ntol,rerr,aerr=',tol,ntol,rerr,aerr
        write (no,*) '   npcs,ipath=',npcs,ipath
        write (no,*) ' here are the whole contents of commons:'
        write (no,*) '   nppar=',nppar,',npdeh=',npdeh
        write (no,*) ' hed:'
        write (no,*) '   nl=',nl,',nl1=',nl1
        write (no,*) '   sigma(0)=',sigma(0),' s/m multiplied by 1.e9'
        write (no,*) '   sigma: ',(sigma(n),n=0,nl1)
        write (no,*) '   thick: ',(thick(n),n=1,nl),' km'
        write (no,*) '   h=',h,',zp=',zp,',zm=',zm,' km'
        write (no,*) '   sl=',sl,',sc=',sc,',phi=',phi,' degrees'
        write (no,*) '   ifdiff=',ifdiff
      end if
      if (nrange.gt.maxrng) then
        write (no,*) ' requested nrange=',nrange,' ranges'
        write (no,*) ' but dimnsn in hedcom.inc is only maxrng=',maxrng
        msg=' too many ranges'
        go to 911
      end if
c   default control parameters for hankel transform routine hy2fht:
      if (.not.((tol.ge.0.0).and.(tol.lt.1))) then
        tol=0.1d-14
        write (no,*) ' tol < 0.0, set to default tol=',tol
      end if
      if ((ntol.lt.1).or.(ntol.gt.10)) then
        ntol=2
        write (no,*) 
     $  ' .not.(ntol.lt.1).or.(ntol.gt.10), set to default ntol=',ntol
      end if
      if (.not.((rerr.ge.0.0).and.(rerr.lt.1))) then
        rerr=0.1d-12
        write (no,*) 
     $  ' .not.((rerr.ge.0.0).and.(rerr.lt.1)), set to default rerr='
     $    ,rerr
      end if
      if (aerr.lt.0.0) then
        aerr=0.01*rerr
        write (no,*) ' aerr < 0.0, set to default aerr=0.01*rerr=',aerr
      end if
      if (.not.((npcs.gt.0).and.(npcs.lt.50))) then
        npcs=1
        write (no,*) 
     $ ' npcs .not.((npcs.gt.0).and.(npcs.lt.50)), set to default npcs='
     $   ,npcs
      end if
      if (.not.((ipath.gt.-3).and.(ipath.lt.3))) then
        ipath=-1
        write (no,*) 
     $ ' ipath .not.between -+3, set to default ipath=',ipath
      end if
c  extended  source  not implemented (af 02/90)
      if (sl.gt.zero) then
        write (no,*) ' **error dehhed: sl=',sl
        msg='  only dipole source (sl=0.0) is implemented'
        go to 911
      end if
      cosp=cos(rad*phi)
      sinp=sin(rad*phi)
      sigscl=(sigma(1)-sigma(0))/(sigma(1)+sigma(0))
c    d(sigscl)/d(sigma1)
      dsgscl=two*sigma(0)/(sigma(0)+sigma(1))**2
      if ((ifdiff.ne.0).and.(nl.eq.0)) then
        write (no,*) ' **error dehhed: no halfspace derivatives'
        msg=' use simpler model, more layers or no derivs'
        stop
      end if
      if (ifdiff.eq.0) then
c    how many field- and derivative-values to compute
c    1 undiff, nl+1 d/dsigmas, nl d/dthicks
        jl=1
        jmax=1
       else if (ifdiff.eq.1) then
        jl=1+nl1
        jmax=jl
       else if (ifdiff.eq.2) then
        jl=1+nl1
        jmax=jl+nl
       else
        write (no,*) ' **error dehhed: ifdiff=',ifdiff
        msg=' ifdiff should be 0,1, or 2'
        go to 911
      end if
c  check array physical dimensions:
      nvals=max(nl1,jmax)
      if (nvals.gt.nppar) then
        write (no,*) ' error dehhed: ',nvals,' field and '
     $  ,'derivative values requested'
     $  ,', or layers needed'
        write (no,*) '   but dimension is only nppar=',nppar
        msg=' check dimensions in hedcom.inc'
        go to 911
      end if
      if (nvals.gt.maxrel) then
        write (no,*) ' **error dehhed: ',nvals,' integrations '
     $  ,'may be required'
        write (no,*) '   but dimension is only maxrel =',maxrel
        msg=' check dimensions in hedcom.inc'
        go to 911
      end if
c  hankel transforms or order 1 and 0 to be computed by hy2fht:
      do 11 j=1,jmax
        nord(j)=1
        nord(jmax+j)=0
11    continue
c  ---------------------------------------------------------
c  compute fields for point source using digital filter integration
c    total number of related integrations
      nj1pm=jmax
      nj1tm=jmax
      if ((phi.eq.0.).or.(phi.eq.180.)) then
        nj0pm=0
       else
        nj0pm=jmax
      end if
      if ((phi.eq.90.).or.(phi.eq.270.)) then
        nj0tm=0
       else
        nj0tm=jmax
      end if
      nrelpm=nj1pm+nj0pm
      nreltm=nj1tm+nj0tm
      if (iw.gt.iw0) then
        write (no,*) ' nj1pm =',nj1pm,',nj1tm =',nj1tm
        write (no,*) ' nj0pm =',nj0pm,',nj0tm =',nj0tm
        write (no,*) ' nrelpm=',nrelpm,',nreltm=',nreltm
      end if
        ommu=tpi*freq*permu
c **double**
        ciommu=dcmplx(zero,ommu)
        gamma2=-ciommu*sigma(0)
        nb=nrange
        do nr=1,nrange
          b(nr)=ranges(nr)
        enddo
c
      if (iw.gt.iw0) then
         write (no,*) ' pass:'
         write (no,*) '   ommu=',ommu
         write (no,*) '   ciommu=',ciommu
         write (no,*) '   sigscl=',sigscl
         write (no,*) '   dsgscl=',dsgscl
         write (no,*) '   gamma2=',gamma2
         write (no,*) '   beta0=',beta0
         write (no,*) '   jl=',jl
         write (no,*) '   jmax=',jmax
         write (no,*) ' direct arguments:'
         write (no,*) '   b(1)=',b(1),', b(',nb,')=',b(nb)
         write (no,*) '   nrelpm=',nrelpm
         write (no,*) '   rerr=',rerr
         write (no,*) '   maxrel=',maxrel
         write (no,*) '   maxrng=',maxrng
         write (no,*) '  tol,ntol,rerr,aerr=',tol,ntol,rerr,aerr 
         write (no,*) '  npcs,ipath=',npcs,ipath 
         write (no,*) '  maxrng,maxrel=',maxrng,maxrel 
         write (no,*) '  array element sampler:'
         write (no,*) '   nord(',maxrel,')=',nord(maxrel)
c  	 pause 
c   calls to integrands so can test them before entering hy2fht.
          x=15.6
          f1=fun1pm(x)
          write (no,*) '  sampled  f1=fun1pm(',x,')=',f1
c		  pause
          call rel1pm(x,f1,fwork,jmax)
          write (no,*) 
     $ '  sampled rel1pm, (fwork(j),j=1,',jmax,')=',(fwork(j),j=1,jmax)
c	      pause
          f1=fun1tm(x)
          write (no,*) '  sampled  f1=fun1tm(',x,')=',f1
c		  pause
          call rel1tm(x,f1,fwork,jmax)
          write (no,*) 
     $ '  sampled rel1tm, (fwork(j),j=1,',jmax,')=',(fwork(j),j=1,jmax)
c
        end if
        if (iw.gt.100) then
          write (no,*) ' dehhed, iw=',iw,', debugging hedmods,' 
          write (no,*) ' stop before hy2fht call to limit output'
          stop
        end if
c
c
c  hankel-transform for the pm mode
          if (iw.gt.iw0) then
            write (no,*) namepu,', next call hy2fht for the pm mode:'
            write (no,*) '  ' 
          end if
          call hy2fht(b,nb,nrelpm,tol,ntol,rerr,aerr,nord,fun1pm,rel1pm
     *    ,zwork,fwork,bwork,npcs,ipath,maxrng,maxrel,zans,nfun1,ierr)
c
          if (iw.ge.iw0) then
            write (no,*) namepu,', back from pm hy2fht, ipath=',ipath
            write (no,*) ' no. hy2fht fun1 evals pm, filter  =',nfun1(1)
            write (no,*) ' no. hy2fht fun1 evals pm, quadrat =',nfun1(2)
                   write (no,*) ' after hy2fht, pm :'
            write (no,*) 
     $      '(zans(1,j),j=1,2*jmax)=',(zans(1,j),j=1,2*jmax)
            write (no,*) 
     $      '(zans(nb,j),j=1,2*jmax)=',(zans(nb,j),j=1,2*jmax)
          end if
c
          ierrd(1)=ierr
          if (ierr.ne.0) then
            write (no,*) ' **error dehhed,pm: ierr=',ierr
            write (no,*) '  ierr=0: no input parameter error detected.'
            write (no,*) '  ierr=1: nb<1 or nb>maxdim not allowed.'
            write (no,*) '  ierr=2: nrel<1 or nrel>maxrel not allowed.'
            write (no,*) '  ierr=3: b(i-1)<=0 or b(i)<=b(i-1), i=2,nb.'
            write (no,*) '    (b(i)=0.0 or not increasing i=1,nb.)'
            write (no,*) '  ierr=4 error ipath=+-1: '
            write (no,*) '    mb=dint(10.*dlog(b(nb)/b(1)))+5>maxdim.'
            write (no,*) '    (mb too large, or maxdim too small.)'
            write (no,*) '  ierr=5: ipath=+-1: dspln1 error return '
            write (no,*) '    after call dspln1.'
            write (no,*) '    (should never occur due to prior checks.)'
            write (no,*) '    ierr=6: only ipath= +- 1,2 allowed'
            write (no,*) '  ierr=11: error ipath=+-1 or -2: error '
            write (no,*) '    return after call hankl2. (see hankl2 '
            write (no,*) '    ierr=1 causes and corrective action)'
            write (no,*) '  ierr=21: ipath=2: error return after'
            write (no,*) '     call besax2. (see besax2 err=1 causes;  '
            write (no,*) '     possible correction is to try '
            write (no,*) '     using npcs>1 and increase rerr, aerr.)'
            write (no,*) '  ierr=31: ipath=-1 or -2: error return after'
            write (no,*) '     call besax2. (see besax2 err=1 causes;  '
            write (no,*) '     possible correction is to try '
            write (no,*) '     using npcs>1 and increase rerr, aerr.)'
            msg=' error-return from hankel-transform hy2fht, pm-mode.'
            go to 911
          end if
c
          do nr=1,nrange
            r=ranges(nr)
            do 25 j=1,jmax
              dej1=zans(nr,j)
              derpm(j,nr)=-gamma2*cosp*dej1/r
              if (nj0pm.eq.jmax) then
                dej0=zans(nr,jmax+j)
                deppm(j,nr)= gamma2*sinp*(dej0-dej1/r)
               else
               deppm(j,nr)=zero
              end if
25          continue
            if (iw.gt.iw0) then
               write (no,*) ' derpm=',(derpm(j,nr),j=1,jmax)
               write (no,*) ' deppm=',(deppm(j,nr),j=1,jmax)
               write (no,*) ' now call hy2fht for tm '
             end if
          enddo
c
c  hankel-transform for the tm mode
          call hy2fht(b,nb,nreltm,tol,ntol,rerr,aerr,nord,fun1tm,rel1tm
     *    ,zwork,fwork,bwork,npcs, ipath,maxrng,maxrel, zans,nfun1,ierr)
c
          if (iw.ge.iw0) then
            write (no,*) namepu,' back from tm hy2fht, ipath=',ipath
            write (no,*) ' no. hy2fht fun1 evals tm, filter  =',nfun1(1)
            write (no,*) ' no. hy2fht fun1 evals tm, quadrat =',nfun1(2)
          end if

c
          if (iw.gt.iw0) then
                    write (no,*) ' after hy2fht, tm :'
             write (no,*) 
     $       '(zans(1,j),j=1,2*jmax)=',(zans(1,j),j=1,2*jmax)
             write (no,*) 
     $       '(zans(nb,j),j=1,2*jmax)=',(zans(nb,j),j=1,2*jmax)
           end if
          ierrd(2)=ierr
          if (ierr.ne.0) then
            write (no,*) ' **dehhed,tm: ierr=',ierr
            write (no,*) '  ierr=0: no input parameter error detected.'
            write (no,*) '  ierr=1: nb<1 or nb>maxdim not allowed.'
            write (no,*) '  ierr=2: nrel<1 or nrel>maxrel not allowed.'
            write (no,*) '  ierr=3: b(i-1)<=0 or b(i)<=b(i-1), i=2,nb.'
            write (no,*) '    (b(i)=0.0 or not increasing i=1,nb.)'
            write (no,*) '  ierr=4 error ipath=+-1: '
            write (no,*) '    mb=dint(10.*dlog(b(nb)/b(1)))+5>maxdim.'
            write (no,*) '    (mb too large, or maxdim too small.)'
            write (no,*) '  ierr=5: ipath=+-1: dspln1 error return '
            write (no,*) '    after call dspln1.'
            write (no,*) '    (should never occur due to prior checks.)'
            write (no,*) '    ierr=6: only ipath= +- 1,2 allowed'
            write (no,*) '  ierr=11: error ipath=+-1 or -2: error '
            write (no,*) '    return after call hankl2. (see hankl2 '
            write (no,*) '    ierr=1 causes and corrective action)'
            write (no,*) '  ierr=21: ipath=2: error return after'
            write (no,*) '     call besax2. (see besax2 err=1 causes;  '
            write (no,*) '     possible correction is to try '
            write (no,*) '     using npcs>1 and increase rerr, aerr.)'
            write (no,*) '  ierr=31: ipath=-1 or -2: error return after'
            write (no,*) '     call besax2. (see besax2 err=1 causes;  '
            write (no,*) '     possible correction is to try '
            write (no,*) '     using npcs>1 and increase rerr, aerr.)'
	    msg=' error-return from hankel-transform hy2fht, tm mode'
            go to 911
          end if
c
          do nr=1,nrange
            r=ranges(nr)
            do 27 j=1,jmax
              dej1=zans(nr,j)
              deptm(j,nr)=-sinp*dej1/r
              if (nj0tm.eq.jmax) then
                dej0=zans(nr,jmax+j)
                dertm(j,nr)= cosp*(dej0-dej1/r)
               else
                dertm(j,nr)=zero
              end if
27          continue
          enddo
            if (iw.gt.iw0) then
              do nr=1,nrange
               write (no,*) ' dertm=',(dertm(j,nr),j=1,jmax)
               write (no,*) ' deptm=',(deptm(j,nr),j=1,jmax)
              enddo
            end if
          
c
c   add an analytically integrated tm contribution
          if (ianatm.eq.1) then
            do nr=1,nrange
              r=ranges(nr)
              r3=r**3
              dertm(1,nr)=dertm(1,nr)-two*sigscl*cosp/r3
              deptm(1,nr)=deptm(1,nr)-sigscl*sinp/r3
              dertm(2,nr)=dertm(2,nr)-two*dsgscl*cosp/r3
              deptm(2,nr)=deptm(2,nr)-dsgscl*sinp/r3
              if (iw.gt.iw0) then
                 xxr =-two*sigscl*cosp/r3
                 write (no,*) 
     $           ' as ianatm=1: add ',xxr,' to dertm(1,nr)'
                 dxxr=-two*dsgscl*cosp/r3
                 write (no,*) 
     $           '          and add ',dxxr,' to dertm(2,nr)'
                 write (no,*) 
     $           ' after modification involving sigscl etc:'
                 write (no,*) ' dertm=',(dertm(j,nr),j=1,jmax)
                 write (no,*) ' deptm=',(deptm(j,nr),j=1,jmax)
              end if
            enddo
          end if
c
c   add the pm and tm contributions, then add the "ww"
          do nr=1,nrange
            r=ranges(nr)
            do 29 j=1,jmax
              derho(j,nr)=derpm(j,nr)+dertm(j,nr)
              dephi(j,nr)=deppm(j,nr)+deptm(j,nr)
29          continue
c    "ww" does not depend on seabed structure (only r,freq), so only contrib to j=1
            t0=ehsr1(r)
            t1=ehsr2(r)/r
            t2=ehsr3(r)
            derww=cosp*(t1-t2)
            depww=sinp*(gamma2*t0+t1)
            derho(1,nr)=derho(1,nr)+derww
            dephi(1,nr)=dephi(1,nr)+depww
            if (iw.gt.iw0) then
               write (no,*) ' t0=',t0,',t1=',t1,',t2=',t2
               write (no,*) ' derww=',derww,',depww=',depww
               write (no,*) ' derho(1,nr)=',derho(1,nr)
               write (no,*) ' dephi(1,nr)=',dephi(1,nr)
             end if
          enddo
      goto 250
c  ---------------------------------------------------------
c  restore mks units:
250   scale=one/(two*tpi*sigma(0))
      if (iw.gt.iw0) then 
        write (no,*) 'scale = ',scale
        write (no,*) ' one=',one,',two=',two,',tpi=',tpi,',sigma(0)='
     $    ,sigma(0)
      end if
      if(sl.ne.0.) then
        scale=scale/(sl*sc)
        if (iw.gt.iw0) then 
          write (no,*) ' as sl.ne.0, modifield scale = ',scale
          write (no,*) ' sl=',sl,',sc=',sc
        end if
      end if
      sscale=scale*tent9
      if (iw.gt.iw0) write (no,*) 'sscale = ',sscale
      tscale=scale*tentm3
      if (iw.gt.iw0) write (no,*) 'tscale = ',tscale
      
      do nr=1,nrange
c                  undifferentiated fields:
            derho(1,nr)=scale*derho(1,nr)
            dephi(1,nr)=scale*dephi(1,nr)
            derpm(1,nr)=scale*derpm(1,nr)
            deppm(1,nr)=scale*deppm(1,nr)
            dertm(1,nr)=scale*dertm(1,nr)
            deptm(1,nr)=scale*deptm(1,nr)
c                sigma-derivatives, for sigma in s/m
        do 265 m=2,jl
            derho(m,nr)=sscale*derho(m,nr)
            dephi(m,nr)=sscale*dephi(m,nr)
            derpm(m,nr)=sscale*derpm(m,nr)
            deppm(m,nr)=sscale*deppm(m,nr)
            dertm(m,nr)=sscale*dertm(m,nr)
            deptm(m,nr)=sscale*deptm(m,nr)
265     continue
c                  thick-derivatives, for thick in m
        do 266 m=jl+1,jmax
            derho(m,nr)=tscale*derho(m,nr)
            dephi(m,nr)=tscale*dephi(m,nr)
            derpm(m,nr)=tscale*derpm(m,nr)
            deppm(m,nr)=tscale*deppm(m,nr)
            dertm(m,nr)=tscale*dertm(m,nr)
            deptm(m,nr)=tscale*deptm(m,nr)
266     continue
        if (iw.gt.iw0) then
          write (no,*) ' dehrho, nr=',nr,':'
           write (no,*) (derho(i,nr),i=1,2*nl+2)
          write (no,*) ' dehphi, nr=',nr,':'
           write (no,*) (dephi(i,nr),i=1,2*nl+2)
        end if
      enddo
 
      if (iw.gt.iw0) write (no,*) 'leaving ',namepu,'.'
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c
c################################################################################
      subroutine fdrlpm(x)
c  zfdpm(0) is the undifferentiated pm mode reflection coefficient
c  zfdpm(m), m=1,2,...,nl1 are its derivatives with respect to the
c  conductivity of layer m, nl1 being the final half-space.
c  the derivative with respect to the source-region conductivity is
c  not computed.
c  the vector of derivatives is passed in common /fdpm/.
c  adapted from reflpm, lines 567-585 in finhed.f77.
c  for tm mode: only definition of alpha,onema2,dsaim1,dsai changes.
c  onema2(i)=1-alpha(i)**2
c  dsai     =derivative of alpha(i)   w.resp.to sigma(i)
c  dsaim1   =  ---"---     alpha(i-1)    --"--   --"--
c  dar(i)   =  ---"---       rl  (i)     --"--  alpha(i)
c  dsri     =  ---"---       rl  (i)     --"--   --"--
c  dsrim1(i)=  ---"---       rl (i-1)    --"--   --"--
c  derr(i)  =  ---"---       rl (i-1)    --"--  e(i)*rl(i)
c  drr(i)   =  ---"---       rl (i-1)    --"--  rl(i)
c  dseri(i) =  ---"---      e(i)*rl(i)   --"--  sigma(i)
c  chain(i) =  ---"---       rl  (0)     --"--   --"-- =chain(i-1)*drr(i)
      include 'hedcom.inc'
      complex*16 dsai,dsaim1,dsrim1,deni2,dsb,dseri
      complex*16 beta(0:nppar),betsum(0:nppar),alpha(0:nppar) 
      complex*16 e(nppar),rl(0:nppar),erl(0:nppar),onema2(0:nppar)
      complex*16 derr(0:nppar),drr(0:nppar),dar(0:nppar) 
      complex*16 chain(0:nppar)
      complex*16 exparg
c
c
      include 'heddat.inc'
      namepu='fdrlpm'
      iw0=100
c  
        if (iw.gt.iw0) then
          write (no,*) ' entering fdrlpm, nl=',nl,',x=',x,':'
        end if
      x2=x*x 
      beta(0)=beta0
c   -----------------------------------------------------------------
      do 10 i=1,nl1
         beta(i)=sqrt(x2-ciommu*sigma(i))
         betsum(i)=beta(i-1)+beta(i)
         alpha(i-1)=(beta(i-1)-beta(i))/betsum(i)
         onema2(i-1)=four*beta(i-1)*beta(i)/betsum(i)**2
         if (iw.gt.iw0) then
            write (no,*) ' beta(',i,')=',beta(i)
            write (no,*) 'betsum(',i,')=',betsum(i)
            write (no,*) 'alpha(',i-1,')=',alpha(i-1)
            write (no,*) 'onema2(',i-1,')=',onema2(i-1)
         end if
 10   continue
      do 11 i=1,nl
         exparg=-two*beta(i)*thick(i)
c     -af explicitly avoid underflows and check for overflows
         if (real(exparg).lt.-exneg) then
            e(i)=zero
         else if (real(exparg).gt.+expos) then
            write (*,*) 
     $           namepu,'        beta,thick(',i,')=',beta(i),thick(i)
            write (*,*) 
     $           ' ***warning ',',exparg=-2*beta(i)*thick(i)=',exparg
            e(i)=zero
         else
            e(i)=exp(exparg)
         end if
         if (iw.gt.iw0) then
            write (no,*) 
     $           ' e(',i,')=exp(-two*beta(i)*thick(i))=exp(',exparg,')=',e(i)
         end if
 11   continue
c     -----------------------------------------------------------------
c           reflection coeff at 1/2-space boundary
      rl(nl)=alpha(nl)
      dar(nl)=one
      erl(nl)=e(nl)*rl(nl)
      if (iw.gt.iw0) then
         write (no,*) ' rl(',nl,')=',rl(nl)
         write (no,*) ' erl(',nl,')=',erl(nl)
         write (no,*) ' dar(',nl,')=',dar(nl)
      end if
c-incomplete          nl=0, just 1/2-space under source-region
c-incomplete           if (nl .lt. 1) then
c-incomplete         d(beta1)/d(sigma1)
c-incomplete             dsb=-half*ciommu/beta(1)
c-incomplete         d(alpha(0))/d(sigma1)
c-incomplete             zfdpm(1)= two*(beta(1)/betsum( 1 )**2)*dsb
c-incomplete             go to 9999
c-incomplete             end if
c     -----------------------------------------------------------------
c     recursion in reflection coeffs up from half-space boundary
      if (nl .lt. 1) go to 201
      do 20 i=nl-1,0,-1
         erl(i+1)=e(i+1)*rl(i+1)
         rl(i)=(alpha(i)+erl(i+1))/(1.+alpha(i)*erl(i+1))
         if (iw.gt.iw0) then
            write (no,*) ' erl(',i+1,')=',erl(i+1)
            write (no,*) ' rl(',i,')=',rl(i)
         end if
 20   continue
 201  continue
c     -----------------------------------------------------------------
c                    refl coeff as returned by reflpm(x,rl)
      zfdpm(0)=rl(0)
c                    only undifferentiated field required
      if (ifdiff.eq.0) go to 9999
c     -----------------------------------------------------------------
c          prepare for chain-rule differentiation of rl(0):
      chain(0)=one
      if (nl .lt. 1) go to 301
      do 30 i=1,nl
        deni2=(one+alpha(i-1)*erl(i))
        deni2=deni2*deni2
        derr(i)=onema2(i-1)/deni2
        drr(i) =e(i)*derr(i)
        dar(i-1)=(one-erl(i)*erl(i))/deni2
        chain(i)=chain(i-1)*drr(i)
          if (iw.gt.iw0) then
            write (no,*) ' derr(',i,')=',derr(i)
            write (no,*) ' drr(',i,')=',drr(i)
            write (no,*) ' dar(',i-1,')=',dar(i-1)
            write (no,*) ' chain(',i,')=',chain(i)
c			pause
           end if
 30     continue
 301    continue
c     -----------------------------------------------------------------
c     compute sigma-derivatives
c     last link in chain-rule,d(rl(i-1)/sigma(i))
c     ds means derivative w.resp.to sigma(i)
        do 40 i=1,nl
c     ds(beta(i))
           dsb=-half*ciommu/beta(i)
c     ds(alpha(i-1)
           dsaim1=-two*(beta(i-1)/betsum( i )**2)*dsb
c     ds(alpha(i)
           dsai  = two*(beta(i+1)/betsum(i+1)**2)*dsb
c     ds(e(i)*r(i))
           dseri=e(i)*(dsai*dar(i)-two*thick(i)*dsb*rl(i))
c     ds(r(i-1))
           dsrim1=dsaim1*dar(i-1)+dseri*derr(i)
c     ds(re(0))
           zfdpm(i)=chain(i-1)*dsrim1
           if (iw.gt.iw0) then
              write (no,*) ' zfdpm(',i,')=',zfdpm(i)
c             pause
           end if
 40     continue
c     --------------------------------------
c     for the last layer, the halfspace, the thickness 
c     is infinite. assuming the resistivity to be finite,
c     e(nl1) will be 0, no matter what sigma(nl1) is.
c     this makes dseri=0 for the halfspace, so skip
c     the second part of dsrim1:
        dsb=-half*ciommu/beta(nl1)
        dsaim1=-two*(beta(nl)/betsum(nl1)**2)*dsb
        dsrim1=dsaim1*dar(nl)
        zfdpm(nl1)=chain(nl)*dsrim1
        if (iw.gt.iw0) then
           write (no,*) ' dsb=',dsb
           write (no,*) ' dsaim1=',dsaim1
           write (no,*) ' zfdpm(',nl1,')=',zfdpm(nl1)
c           pause
        end if
c     -----------------------------------------------------------------
c     compute thickness-derivatives
        if (ifdiff.eq.2) then
           do 50 j=jl+1,jmax
              layer=j-jl
              zfdpm(j-1)=chain(layer-1)*derr(layer)*erl(layer)
     $             *(-two*beta(layer))
              if (iw.gt.iw0) then
                 write (no,*) ' in 50-loop, j=',j,', layer=',layer
                 write (no,*) '   zfdpm(',j-1,')=',zfdpm(j-1)
              end if
 50        continue
        end if
        if (iw.gt.iw0) then
           write (no,*) ' leaving fdrlpm, zfdpm='
           write (no,*) (zfdpm(j),j=0,jmax-1)
        end if
c     
 9999   return
 911    continue
        write (no,*) ' stop due to error in ',namepu
        write (no,*) msg
        stop
        end
c ###################################################################
      subroutine fdrltm(x)
c computes the tm mode reflection coefficient (rl(0)=zfdtm(0)) for layered
c structure. if jmax>1, it also computes its derivatives with respect to
c the conductivities of layers 1 through nl1 (zfdtm(1),zfdtm(2),...,zfdtm(nl1))
c and if jmax>jl, the derivatives with respect to layer thicknesses
c (zfdtm(nl+2),...,zfdtm(2*nl1).
c does not differentiate with respect to the source region conductivity.
c the vector of derivatives is passed in common /fdtm/.
c adapted from refltm, lines 5..-5.. in finhed.f77.
c  modified: alpha,onema2,dsaim1,dsai.  all else is same as in fdrlpm.
c  onema2(i)=1-alpha(i)**2
c  dsai     =derivative of alpha(i)   w.resp.to sigma(i)
c  dsaim1   =  ---"---     alpha(i-1)    --"--   --"--
c  dar(i)   =  ---"---       rl  (i)     --"--  alpha(i)
c  dsri     =  ---"---       rl  (i)     --"--   --"--
c  dsrim1(i)=  ---"---       rl (i-1)    --"--   --"--
c  derr(i)  =  ---"---       rl (i-1)    --"--  e(i)*rl(i)
c  drr(i)   =  ---"---       rl (i-1)    --"--  rl(i)
c  dseri(i) =  ---"---      e(i)*rl(i)   --"--  sigma(i)
c  chain(i) =  ---"---       rl  (0)     --"--   --"-- =chain(i-1)*drr(i)
      include 'hedcom.inc'
      complex*16 dsai,dsaim1,dsrim1
      complex*16 deni2,dsb,xxxx,dseri
      complex*16 beta(0:nppar),betsum(0:nppar),alpha(0:nppar) 
      complex*16 e(nppar),rl(0:nppar),erl(0:nppar),onema2(0:nppar)
      complex*16 derr(0:nppar),drr(0:nppar),dar(0:nppar)
      complex*16 chain(0:nppar)
      complex*16 exparg
c
      include 'heddat.inc'
      namepu='fdrltm'
      iw0=100
c
      x2=x*x
      beta(0)=beta0
c     -----------------------------------------------------------------
      do 10 i=1,nl1
        beta(i)=sqrt(x2-ciommu*sigma(i))
        betsum(i) = beta(i-1)*sigma(i)+beta(i)*sigma(i-1)
        alpha(i-1)=(beta(i-1)*sigma(i)-beta(i)*sigma(i-1))/betsum(i)
        onema2(i-1)=one - alpha(i-1)**2
 10   continue
      do 11 i=1,nl
          exparg=-two*beta(i)*thick(i)
c-af explicitly avoid underflows and check for overflows
          if (real(exparg).lt.-exneg) then
            e(i)=zero
           else if (real(exparg).gt.+expos) then
            write (*,*) 
     $      namepu,'        beta,thick(',i,')=',beta(i),thick(i)
            write (*,*) 
     $      ' ***warning ',',exparg=-2*beta(i)*thick(i)=',exparg
            e(i)=zero
           else
            e(i)=exp(exparg)
          end if
11    continue
c     -----------------------------------------------------------------
c              reflection coeff at 1/2-space boundary
      rl(nl)=alpha(nl)
      dar(nl)=one
      erl(nl)=e(nl)*rl(nl)
c-incomplete     nl=0, just 1/2-space under source-region
c-incomplete      if (nl .lt. 1) then
c-incomplete    d(beta1)/d(sigma1)
c-incomplete        dsb=-half*ciommu/beta(1)
c-incomplete    d(alpha(0))/d(sigma1)
c-incompletec     zfdtm(1)= two*(beta(1)/betsum( 1 )**2)*dsb
c-incomplete        go to 9999
c-incomplete        end if
c     -----------------------------------------------------------------
c             recursion in reflection coeffs up from half-space boundary
      if (nl .lt. 1) go to 201
      do 20 i=nl-1,0,-1
        erl(i+1)=e(i+1)*rl(i+1)
        rl(i)=(alpha(i)+erl(i+1))/(1.+alpha(i)*erl(i+1))
   20   continue
201   continue
c    refl coeff as returned by refltm(x,rl)
      zfdtm(0)=rl(0)
c    only undifferentiated field required
      if (ifdiff.eq.0) go to 9999
c     -----------------------------------------------------------------
c             prepare for chain-rule differentiation of rl(0):
      chain(0)=one
      if (nl .lt. 1) go to 301
      do 30 i=1,nl
        deni2=(one+alpha(i-1)*erl(i))
        deni2=deni2*deni2
c     d(r(i-1))/d(e(i)*r(i))
        derr(i)=onema2(i-1)/deni2
c    d(r(i-1))/dr(i)
        drr(i) =e(i)*derr(i)
c     d(r(i-1))/d(a(i-1))
        dar(i-1)=(one-erl(i)*erl(i))/deni2
c     d(r(0))/d(r(i))
        chain(i)=chain(i-1)*drr(i)
   30   continue
301   continue
c     -----------------------------------------------------------------
c                compute sigma-derivatives
c             last link in chain-rule,d(rl(i-1)/sigma(i))
c           ds means derivative w.resp.to sigma(i)
      do 40 i=1,nl
c    ds(beta(i))
        dsb=-half*ciommu/beta(i)
        xxxx=sigma(i)*dsb-beta(i)
c    ds(alpha(i-1)
        dsaim1=-two*sigma(i-1)*(beta(i-1)/betsum( i )**2)*xxxx
c    ds(alpha(i)
        dsai  = two*sigma(i+1)*(beta(i+1)/betsum(i+1)**2)*xxxx
c    ds(e(i)*r(i))
        dseri=e(i)*(dsai*dar(i)-two*thick(i)*dsb*rl(i))
c    ds(r(i-1))
        dsrim1=dsaim1*dar(i-1)+dseri*derr(i)
c    ds(re(0))
        zfdtm(i)=chain(i-1)*dsrim1
   40   continue
c      --------------------------------------
c             for the last layer, the halfspace, the thickness 
c             is infinite. assuming the resistivity to be finite,
c             e(nl1)=0, no matter what sigma(nl1) is.
c             this makes dseri=0 for the halfspace, so skip
c             the second part of dsrim1:
      dsb=-half*ciommu/beta(nl1)
      xxxx=sigma(nl1)*dsb-beta(nl1)
      dsaim1=-two*sigma(nl)*(beta(nl)/betsum(nl1)**2)*xxxx
      dsrim1=dsaim1*dar(nl)
      zfdtm(nl1)=chain(nl)*dsrim1
c     -----------------------------------------------------------------
c            compute thickness-derivatives
      if (ifdiff.eq.2) then
        do 50 j=jl+1,jmax
          layer=j-jl
          zfdtm(j-1)=chain(layer-1)*derr(layer)*erl(layer)
     $           *(-two*beta(layer))
50      continue
      end if
c      write (no,*) ' from fdrltm, zfdtm:'
c      do 991 i = 1,jmax 
c 991     write (no,*) zfdtm(i)
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c ###################################################################
      complex*16 function ehsr1(r)
c  analytic solution of integral from 0 to infinity
c  j0(k*r)*k/beta0*exp(-beta0*zm)*dk
      include 'hedcom.inc'
      complex*16 gamma
c
      include 'heddat.inc'
      namepu='ehsr1'
      iw0=100
c
c
      arg=sqrt(r*r+zm*zm)
      gamma=sqrt(gamma2)
      ehsr1=exp(-gamma*arg)/arg
chat           write (no,*) 'ehsr1',r,ehsr1
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c  ################################################################################
      complex*16 function ehsr2(r)
c  analytic solution of integral from 0 to infinity
c  j1(k*r)*k**2/beta0*exp(-beta0*zm)*dk
      include 'hedcom.inc'
      complex*16 gamma
c
      include 'heddat.inc'
      namepu='ehsr2'
      iw0=100
c
c
      arg=sqrt(r*r+zm*zm)
      gamma=sqrt(gamma2)
      ehsr2=r*exp(-gamma*arg)*(gamma*arg+one)/arg**3
chat            write (no,*) 'ehsr2',r,ehsr2
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c  ################################################################################
      complex*16 function ehsr3(r)
c  analytic solution of integral from 0 to infinity
c  j0(k*r)*k*beta0*exp(-beta0*zm)*dk
      include 'hedcom.inc'
      complex*16 gamma
c
      include 'heddat.inc'
      namepu='ehsr3'
      iw0=100
c
c
      arg=sqrt(r*r+zm*zm)
      gamma=sqrt(gamma2)
      ehsr3=exp(-gamma*arg)*((gamma*arg*zm)**2+(gamma*arg+one)*
     $      (two*zm**2-r**2))/arg**5
chat            write (no,*) 'ehsr3',r,ehsr3
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c ###################################################################
      complex*16 function fun1pm(x)
c  computes kernel for first order pm mode hankel transform
c  saves d(deh1pm)/d(rl)=d(kernel)/d(refl coeff) for use by rel1pm
c  adapted from ehj1pm, lines 391-431 in a.d.chave's finhed.f77, 29 sept 86  -af
c  last modified june 93
c  agusta flosadottir, sio & uw
c
c  called by hy2fht as fun1
c  calls     fdrlpm and uses common-blocks set up by it
c            dedrpm in common
c
c  nparam is the parameter being differentiated with respect to:
c   nparam=0:             undifferentiated
c   nparam=1,nl+1         d/d(sigma) from (rock) layer 1 to halfspace
c   nparam=nl+2,jmax      d/d(thick) from (rock) layer 1 to layer nl
c  jmax is the total number of field and derivative values
c  computed, jmax=1+number of parameters=1+(nl+1)+nl=jl+nl
c
      include 'hedcom.inc'
      complex*16 rl,e,ra,a,b,c,d,denom
c-faster      complex*16 emb0h2,epb0zp,epb0zm
c
      include 'heddat.inc'
      namepu='fun1pm'
      iw0=100
        if (iw.gt.iw0) then
          write (no,*) ' entering ',namepu,': x=',x
          write (no,*) ' common hed contains:'
          write (no,*) '   nl=',nl,',nl1=',nl1
          write (no,*) '   sigma: ',(sigma(n),n=0,nl1)
          write (no,*) '   thick: ',(thick(n),n=1,nl)
          write (no,*) '   h=',h,',zp=',zp,',zm=',zm
          write (no,*) '   sl=',sl,',sc=',sc,',phi=',phi
          write (no,*) '   ifdiff=',ifdiff
          write (no,*) '   rerr=',rerr,',aerr=',aerr,',ng=',ng
c		  pause
        end if
c
      beta0=sqrt(x*x+gamma2)
c    compute refl coeff & its sigma-derivatives
      call fdrlpm(x)
        if (iw.gt.iw0) then
          write (no,*) ' ',namepu,': back from fdrlpm-call at x=',x,':'
          write (no,*) ' common pass contains:'
          write (no,*) '   ommu=',ommu
          write (no,*) '   ciommu=',ciommu
          write (no,*) '   sigscl=',sigscl
          write (no,*) '   dsgscl=',dsgscl
          write (no,*) '   gamma2=',gamma2
          write (no,*) '   beta0=',beta0
          write (no,*) '   jl=',jl
          write (no,*) '   jmax=',jmax
          write (no,*) '  other results:'
          write (no,*) '   zfdpm:',(zfdpm(j),j=1,jmax-1)
c		  pause 
        end if
      rl=zfdpm(0)
      if(zp.ne.0.)then
        if(h.eq.0.)then
          if (iw.gt.iw0) write (no,*) ' zp.ne.0 and h.eq.0 :'
          dedrpm=exp(-beta0*zp)/beta0
          fun1pm=dedrpm*rl
        else
         if (iw.gt.iw0) write (no,*) ' not:  zp.ne.0 and h.eq.0 :'
         if (iw.gt.iw0) write (no,*) '   one=',one,',two=',two
          ra=(beta0-x)/(beta0+x)
          h2=h+h
c
c-faster          emb0h2=exp(-beta0*h2)
c-faster          epb0zp=exp(beta0*zp)
c-faster          epb0zm=exp(beta0*zm)
c-faster          a=ra*epb0zp*emb0h2
c-faster          b=one/epb0zp+ra*emb0h2*(epb0zm+one/epb0zm)
c-faster          c=beta0
c-faster          d=-beta0*ra*emb0h2
c-faster          denom=c+d*rl
c-faster          fun1pm=(a+b*rl)/denom
c-faster          dedrpm=(b-d*fun1pm)/denom
c
c  slower version to avoid trouble at large beta0
c  only decaying exponentials, no divisions
          deepm=real(beta0)*(h2-zm)
          deepp=real(beta0)*(h2-zp)
          if ((deepm.gt.250.0).and.(deepp.gt.250)) then
c    avoid underflows & save time at very large beta0
            a=zero 
            b=exp(-beta0*zp)
            c=beta0
            d=zero
           else
            a=ra*exp(beta0*(-h2+zp)) 
            b=exp(-beta0*zp)
     $       +ra*(exp(beta0*(-h2+zm))+exp(beta0*(-h2-zm)))
            c=beta0
            d=-beta0*ra*exp(beta0*(-h2))
          end if
          denom=c+d*rl
          fun1pm=(a+b*rl)/denom
          dedrpm=(b-d*fun1pm)/denom
c
        endif
      else
        if(h.eq.0.)then
          dedrpm=one/beta0
          fun1pm=dedrpm*rl
        else
          ra=(beta0-x)/(beta0+x)
          h2=h+h
          e=exp(-beta0*h2)
          a=ra*e
          b=one+ra*e*two
          c=beta0
          d=-beta0*ra*e
          denom=c+d*rl
          fun1pm=(a+b*rl)/denom
          dedrpm=(b-d*fun1pm)/denom
        endif
      endif
c
        if (iw.gt.iw0) then
          write (no,*) ' at end of  ',namepu,'.  '
          write (no,*) ' fun1pm(',x,')=',fun1pm,',dedrpm=',dedrpm
        end if
c
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end




      complex*16 function fun1tm(x)
c  computes kernel for first order tm mode hankel transform
c  saves d(fun1tm)/d(rl)=d(kernel)/d(refl coeff) for use by rel1tm
c  adapted from ehj1tm, lines 480-519
c  and from ehtmre, lines 520-527 in finhed.f77, 29sept86  -af
c  last modified june 93
c  agusta flosadottir, sio & uw

c sign errors in air wave fixed by lucy macgregor may 2000
c
c  nparam is the parameter being differentiated with respect to:
c   nparam=0:             undifferentiated
c   nparam=1,nl+1         d/d(sigma) from (rock) layer 1 to halfspace
c   nparam=nl+2,jmax      d/d(thick) from (rock) layer 1 to layer nl
c  jmax is the total number of field and derivative values
c  computed, jmax=1+number of parameters=1+(nl+1)+nl=jl+nl
c
      include 'hedcom.inc'
      complex*16 rl,e,a,b,c,d,denom
c-faster      complex*16 emb0h2,epb0zp,epb0zm
c
      include 'heddat.inc'
      namepu='fun1tm'
      iw0=100
c
      beta0=sqrt(x*x+gamma2)
c    compute refl coeff & its sigma-derivatives
      call fdrltm(x)
      rl=zfdtm(0)
      if(zp.ne.0.)then
        if(h.eq.0.)then
          dedrtm=beta0*exp(-beta0*zp)
          fun1tm=dedrtm*rl
        else
          h2=h+h

c-faster          emb0h2=exp(-beta0*h2)
c-faster          epb0zp=exp(beta0*zp)
c-faster          epb0zm=exp(beta0*zm)
c-fasterc    const part of numerator
c-faster          a=epb0zp*emb0h2
c-faster          b=one/epb0zp
c-fasterc    lin part of do
c-faster     $     -emb0h2*(epb0zm+one/epb0zm)
c-fasterc    const  part of denominator
c-faster          c=one
c-fasterc    linear part
c-faster          d=emb0h2
c-fasterc
c  slower version to avoid trouble at large beta0, 
c  only decaying exponentials, no divisions
          decexp=real(beta0)*(h2-zm-zp)
          if (decexp.gt.250.0) then
c    avoid underflows and save time at very large beta0
            a=zero
            b=exp(-beta0*zp) 
            c=one
            d=zero
           else
cxxxc            a=exp(beta0*(-h2+zp))
cxxxc            b=exp(-beta0*zp) 
cxxxc     $       - (exp(beta0*(-h2+zm))+exp(beta0*(-h2-zm)))
            a=exp(beta0*(-h2+zp))
            b=exp(-beta0*zp)
     $       + (exp(beta0*(-h2+zm))+exp(beta0*(-h2-zm)))
            c=one
            d=exp(-beta0*h2)
          end if
c
          denom=one+d*rl
cxxxc          fun1tm=beta0*(a+b*rl)/denom
	  fun1tm=beta0*(b*rl-a)/denom
          dedrtm=(beta0*b-d*fun1tm)/denom
        endif
      else
        if(h.eq.0.)then
          dedrtm=beta0
          fun1tm=beta0*rl
        else
          h2=h+h
          e=exp(-beta0*h2)
          a=-e
          b=one+two*e
          d=e
          denom=one+rl*d
          fun1tm=beta0*(a+b*rl)/denom
          dedrtm=(beta0*b-d*fun1tm)/denom
        endif
c   only for tm, does not contribute to d(eh)/d(rl), see fun1tm
c        if (method.eq.0) fun1tm=fun1tm-sigscl*x
        fun1tm=fun1tm-sigscl*x
      endif
c  ------------------------------------------
        if (iw.gt.iw0) then
          write (no,*) ' at end of  ',namepu,'.  '
          write (no,*) ' fun1tm(',x,')=',fun1tm,',dedrtm=',dedrtm
        end if
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c ###################################################################
      subroutine rel1pm(x,f1,fj,jarg)
c  adapted from a.d.chave's ehpmre, lines 520-527 in finhed.f77, 29sept86  -af
c  last modified june 93 
c
c  this function is used as "rel1" in hy2fht. 
c  it uses zfdpm, the vector
c  of partial derivatives of the pm reflection coeff rl with respect to
c  the conductivities of layers m=1,2,..,nl+1, ( zfdpm(0)=rl(0)=the
c  undifferentiated reflection coefficient),
c  and dedrpm, the derivative of the kernel deh1pm with respect to rl,
c  which were stored during an earlier call to "fun1"=fun1pm.
c  the undiff kernel is z=deh1pm(x) [ =zwork(ir,j=1)=fun1(x) in zhankl ]
c
c  jmax is the total number of field and derivative values computed, 
c  it is 1+number of parameters = 1+(nl+1)+nl=jl+nl
c  (from common)
c
c  nparam is the parameter being differentiated with respect to:
c   nparam=0:             undifferentiated
c   nparam=1,nl+1         d/d(sigma) from (rock) layer 1 to halfspace
c   nparam=nl+2,jmax      d/d(thick) from (rock) layer 1 to layer nl
c
      include 'hedcom.inc'
      complex*16 f1,fj(*)
      complex*16 deh1pm
c
      include 'heddat.inc'
      namepu='rel1pm'
      iw0=100
      if (iw.gt.iw0) then
        write (no,*) ' entering ',namepu,': x=',x
      end if
      deh1pm=f1
c
      jaabs=abs(jarg)
      if (jaabs.gt.2*jmax) then
       write (no,*) ' **error ',namepu,': jaabs=',jaabs
     $   ,', but jmax is only ',jmax
        msg=' jaabs > 2*jmax'
        go to 911
      end if
c  -----------------------------
c  branch for evaluating one at a time
c  do not use anderson's computed goto structure, 
c  this was easier to program 
      if (jarg.lt.0) then
        if (jaabs.eq.1) then
          fj(1)=deh1pm
         else if (jaabs.le.jmax) then
          nparam=jaabs-1
          fj(jaabs)=dedrpm*zfdpm(nparam)
         else if (jaabs.le.jmax+1) then
          fj(jaabs)=x*deh1pm
         else if (jaabs.le.2*jmax) then
          nparam=jaabs-jmax-1
          fj(jaabs)=x*dedrpm*zfdpm(nparam)
         else
          write (no,*) '**error: jaabs,jmax=',jaabs,jmax
          msg=' jaabs > 2*jmax'
          go to 911
        end if
        if (iw.gt.iw0) then
          write (no,*) '  at x=,x:'
          write (no,*) '  deh1pm,dedrpm=',deh1pm,dedrpm
          write (no,*) '  fj(',jaabs,'): ',fj(jaabs)
        end if
c  -----------------------------
c  branch for evaluating all related in one call
       else
c  the initial order:
        fj(1)=deh1pm
        n1=min(jmax,jaabs)
        do 160 j=2,n1
          nparam=j-1
          fj(j)=dedrpm*zfdpm(nparam)
160     continue
        if (iw.gt.iw0) write (no,*) '  fj after initial order:'
c  the next order is the previous times x:
        fj(n1+1)=x*deh1pm
        do 165 j=n1+2,jaabs
          nparam=j-jmax-1
          fj(j)=x*dedrpm*zfdpm(nparam)
165     continue
        if (iw.gt.iw0) then
          write (no,*) '  at x=,x:'
          write (no,*) '   deh1pm,dedrpm=',deh1pm,dedrpm
          write (no,*) '   fj(j),j=1,n1: ',(fj(j),j=1,n1)
          write (no,*) '   fj(j),j=n1+1,nrel: ',(fj(j),j=n1+1,nrel)
        end if
c  -----------------------------
      end if
c
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c ###################################################################
      subroutine rel1tm(x,f1,fj,jarg)
c  adapted from a.d.chave's ehtmre in finhed.f77, 29sept86  -af
c  last modified june 93 
c
c  this function is used as "rel1" in hy2fht. 
c  it uses zfdtm, the vector
c  of partial derivatives of the tm reflection coeff rl with respect to
c  the conductivities of layers m=1,2,..,nl+1, ( zfdtm(0)=rl(0)=the
c  undifferentiated reflection coefficient),
c  and dedrtm, the derivative of the kernel deh1tm with respect to rl,
c  which were stored during an earlier call to "fun1"=fun1tm.
c  the undiff kernel is z=deh1tm(x) [ =zwork(ir,j=1)=fun1(x) in zhankl ]
c
c  jmax is the total number of field and derivative values computed, 
c  it is 1+number of parameters = 1+(nl+1)+nl=jl+nl
c  (from common)
c
c  nparam is the parameter being differentiated with respect to:
c   nparam=0:             undifferentiated
c   nparam=1,nl+1         d/d(sigma) from (rock) layer 1 to halfspace
c   nparam=nl+2,jmax      d/d(thick) from (rock) layer 1 to layer nl
c
      include 'hedcom.inc'
      complex*16 f1,fj(*)
      complex*16 deh1tm
c
      include 'heddat.inc'
      namepu='rel1tm'
      iw0=100
      if (iw.gt.iw0) then
        write (no,*) ' entering ',namepu,': x=',x
      end if
      deh1tm=f1
c
      jaabs=abs(jarg)
      if (jaabs.gt.2*jmax) then
       write (no,*) ' **error ',namepu,': jaabs=',jaabs
     $   ,', but jmax is only ',jmax
        msg=' jaabs > 2*jmax'
        go to 911
      end if
c  -----------------------------
c  branch for evaluating one at a time
c  do not use anderson's computed goto structure, 
c  this was easier to program 
      if (jarg.lt.0) then
        if (jaabs.eq.1) then
          fj(1)=deh1tm
         else if (jaabs.le.jmax) then
          nparam=jaabs-1
          fj(jaabs)=dedrtm*zfdtm(nparam)
         else if (jaabs.le.jmax+1) then
          fj(jaabs)=x*deh1tm
         else if (jaabs.le.2*jmax) then
          nparam=jaabs-jmax-1
          fj(jaabs)=x*dedrtm*zfdtm(nparam)
         else
          write (no,*) '**error: jaabs,jmax=',jaabs,jmax
          msg=' jaabs > 2*jmax'
          go to 911
        end if
        if (iw.gt.iw0) then
          write (no,*) '  at x=,x:'
          write (no,*) '  deh1tm,dedrtm=',deh1tm,dedrtm
          write (no,*) '  fj(',jaabs,'): ',fj(jaabs)
        end if
c  -----------------------------
c  branch for evaluating all related in one call
       else
c  the initial order:
        fj(1)=deh1tm
        n1=min(jmax,jaabs)
        do 160 j=2,n1
          nparam=j-1
          fj(j)=dedrtm*zfdtm(nparam)
160     continue
        if (iw.gt.iw0) write (no,*) '  fj after initial order:'
c  the next order is the previous times x:
        fj(n1+1)=x*deh1tm
        do 165 j=n1+2,jaabs
          nparam=j-jmax-1
          fj(j)=x*dedrtm*zfdtm(nparam)
165     continue
        if (iw.gt.iw0) then
          write (no,*) '  at x=,x:'
          write (no,*) '   deh1tm,dedrtm=',deh1tm,dedrtm
          write (no,*) '   fj(j),j=1,n1: ',(fj(j),j=1,n1)
          write (no,*) '   fj(j),j=n1+1,nrel: ',(fj(j),j=n1+1,nrel)
        end if
c  -----------------------------
      end if
c
9999  return
911    continue
       write (no,*) ' stop due to error in ',namepu
       write (no,*) msg
       stop
      end
c ###################################################################
c ###################################################################
      double precision function bessj(x,order)
c
c     called by hy2fht (w.anderson's hybrid code), quadrature option only.
c     requires order 0,1 bessel function subroutines
c     if only using digital-filter method, dummy routines can be
c     loaded instead.
c     nb:  bessel functions are included with anderson's single
c          precision version, but he suggests slatec for use with
c          hy2fht.
c
      include 'hedcom.inc'
      double precision x,order
c-amos      double precision zr,zi,fnu,cyr(1),cyi(1)
      double precision besj0,besj1
c
      include 'heddat.inc'
      namepu='bessj'
      iw0=100
c
c-amos      iwhich=0
c-amos      if (iwhich.eq.1) then
c-amosc  this branch calls zbesj (from amos library, obtained 
c-amosc  through netlib@ornl.gov).
c-amosc  zbesj allows complex-argument and non-neg.order.
c-amosc  here the argument & the result will always be real.
c-amos        zr=x
c-amos        zi=0.0d+00
c-amos        fnu=order
c-amos        nordrs=1
c-amos        kode=1
c-amos        call zbesj(zr,zi,fnu,kode,nordrs,cyr,cyi,nz,ierbes)
c-amos        if (ierbes.ne.0) then
c-amos          write (no,*) 
c-amos     $ ' zbesj failed, ierbes=',ierbes,', at zr,zi,fnu=',zr,zi,fnu
c-amos          go to 911
c-amos        end if
c-amos        bessj=cyr(1)
c-amos       else
c  this branch calls r.parker's bessel functions (personal comm., 05/93)
c  they are based on slatec routines, improved large-argumt handling.
        if (abs(order).lt.0.001) then
          bessj=besj0(x)
         else if (abs(order-one).lt.0001) then
          bessj=besj1(x)
         else
          write (no,*) ' ###error bessj: order =',order,', not 0 or 1'
          stop
        end if
c-amos      end if
c
9999  return
911   write (no,*) ' **stop due to error in ',namepu
      write (no,*) msg
      stop
      end




c  besj0(x), besj1(x), besy0(x), besy1(x), besi0(x), besi0e(x), besi1(x),
c  besi1e(x), besk0(x), besk0e(x), besk1(x), besk1e(x),
c
c  double-precision routines for bessel functions of real argument:
c  besj0:  j0 = bessel funtion 1st kind, 0th order
c  besj1:  j1 = bessel funtion 1st kind, 1st order
c  besy0:  y0 = bessel funtion 2nd kind, 0th order
c  besy1:  y1 = bessel funtion 2nd kind, 1st order
c  besi0:  i0 = modified bessel funtion 1st kind, 0th order
c  besi1:  i1 = modified bessel funtion 1st kind, 1st order
c  besk0:  k0 = modified bessel funtion 2nd kind, 0th order
c  besk1:  k1 = modified bessel funtion 2nd kind, 1st order
c
c [useful because size remains reasonable for large aruments]
c  besi0e:  i0 * exp(-x) 
c  besi1e:  i1 * exp(-x)
c  besk0e:  k0 * exp(x)
c  besk1e:  k1 * exp(x)
c
c  accuracy of about 3 in 10**16; checks for overflow, underflow and
c  loss of precision.  adopted for sun sparc station from slatec 
c  library.
c
c  modification af:  
c  "function bes.."   => "double precision function bes.."
c  just to be on the safe side of unknown compilers
c
      double precision function besj0 (x)
c$$$$ calls csevl
      implicit double precision (a-h, o-z)
      dimension bj0cs(13), bm0cs(21), bth0cs(24)
      data bj0cs/                 .1002541619 68939137d+0,
     $ -.6652230077 64405132d+0,  .2489837034 98281314d+0,
     $ -.0332527231 700357697d+0, .0023114179 304694015d+0,
     $ -.0000991127 741995080d+0, .0000028916 708643998d+0,
     $ -.0000000612 108586630d+0, .0000000009 838650793d+0,
     $ -.0000000000 124235515d+0, .0000000000 001265433d+0,
     $ -.0000000000 000010619d+0, .0000000000 000000074d+0/
      data bm0cs/
     $ .0928496163 7381644d+0, -.0014298770 7403484d+0,
     $ .0000283057 9271257d+0, -.0000014330 0611424d+0,
     $ .0000001202 8628046d+0, -.0000000139 7113013d+0,
     $ .0000000020 4076188d+0, -.0000000003 5399669d+0,
     $ .0000000000 7024759d+0, -.0000000000 1554107d+0,
     $ .0000000000 0376226d+0, -.0000000000 0098282d+0,
     $ .0000000000 0027408d+0, -.0000000000 0008091d+0,
     $ .0000000000 0002511d+0, -.0000000000 0000814d+0,
     $ .0000000000 0000275d+0, -.0000000000 0000096d+0,
     $ .0000000000 0000034d+0, -.0000000000 0000012d+0,
     $ .0000000000 0000004d+0/
      data bth0cs /
     $   -.2463916377 4300119d+0,  .0017370983 07508963d+0,
     $   -.0000621836 33402968d+0, .0000043680 50165742d+0,
     $   -.0000004560 93019869d+0, .0000000621 97400101d+0,
     $   -.0000000103 00442889d+0, .0000000019 79526776d+0,
     $   -.0000000004 28198396d+0, .0000000001 02035840d+0,
     $   -.0000000000 26363898d+0, .0000000000 07297935d+0,
     $   -.0000000000 02144188d+0, .0000000000 00663693d+0,
     $   -.0000000000 00215126d+0, .0000000000 00072659d+0,
     $   -.0000000000 00025465d+0, .0000000000 00009229d+0,
     $   -.0000000000 00003448d+0, .0000000000 00001325d+0,
     $   -.0000000000 00000522d+0, .0000000000 00000210d+0,
     $   -.0000000000 00000087d+0, .0000000000 00000036d+0/
      data pi4 / 0.7853981633 9744831d+0 /
      data ntj0, ntm0, ntth0, xsml,      xmax / 
     $  12,        21,    24, 2.107d-08, 4.5036d+15/
ce
      y = abs(x)
      if (y.le.4.0) then
c  small x expansion
        besj0 = 1.0
        if (y.gt.xsml) besj0 = csevl (y*y/8.0-1.0, bj0cs, ntj0)
c  large x expansion
      elseif (y.lt.xmax) then
        z = 32.0/y**2 - 1.0
        ampl = (0.75d+0 + csevl (z, bm0cs, ntm0)) / sqrt(y)
        theta = y - pi4 + csevl (z, bth0cs, ntth0) / y
        besj0 = ampl * cos (theta)
c  x too large
      else
        print *,' in besj0 no precision because abs(x) is big: x=',x
        stop
      endif
c
      return
      end
c______________________________________________________________________
      double precision function besj1 (x)
c$$$$ calls csevl
      implicit double precision (a-h, o-z)
      dimension bj1cs(12), bm1cs(21), bth1cs(24)
      data bj1 cs /   
     $ -.1172614151 3332787d+0, -.2536152183 0790640d+0,
     $  .0501270809 84469569d+0, -.0046315148 09625081d+0,
     $  .0002479962 29415914d+0, -.0000086789 48686278d+0,
     $  .0000002142 93917143d+0, -.0000000039 36093079d+0,
     $  .0000000000 55911823d+0, -.0000000000 00632761d+0,
     $  .0000000000 00005840d+0, -.0000000000 00000044d+0/
      data bm1 cs /    .1047362510 931285d+0,
     $  .0044244389 3702345d+0, -.0000566163 9504035d+0,
     $  .0000023134 9417339d+0, -.0000001737 7182007d+0,
     $  .0000000189 3209930d+0, -.0000000026 5416023d+0,
     $  .0000000004 4740209d+0, -.0000000000 8691795d+0,
     $  .0000000000 1891492d+0, -.0000000000 0451884d+0,
     $  .0000000000 0116765d+0, -.0000000000 0032265d+0,
     $  .0000000000 0009450d+0, -.0000000000 0002913d+0,
     $  .0000000000 0000939d+0, -.0000000000 0000315d+0,
     $  .0000000000 0000109d+0, -.0000000000 0000039d+0,
     $  .0000000000 0000014d+0, -.0000000000 0000005d+0/
      data bth1cs /    
     $   .7406014102 6313850d+0, -.0045717556 59637690d+0,
     $  .0001198185 10964326d+0, -.0000069645 61891648d+0,
     $  .0000006554 95621447d+0, -.0000000840 66228945d+0,
     $  .0000000133 76886564d+0, -.0000000024 99565654d+0,
     $  .0000000005 29495100d+0, -.0000000001 24135944d+0,
     $  .0000000000 31656485d+0, -.0000000000 08668640d+0,
     $  .0000000000 02523758d+0, -.0000000000 00775085d+0,
     $  .0000000000 00249527d+0, -.0000000000 00083773d+0,
     $  .0000000000 00029205d+0, -.0000000000 00010534d+0,
     $  .0000000000 00003919d+0, -.0000000000 00001500d+0,
     $  .0000000000 00000589d+0, -.0000000000 00000237d+0,
     $  .0000000000 00000097d+0, -.0000000000 00000040d+0/
      data pi4 / 0.7853981633 9744831d+0 /
      data ntj1, ntm1, ntth1, xsml,          xmin, xmax / 
     $       12,  21,   24, 2.981d-08, 9.889d-324, 4.503d+15/
ce
c
      y = abs(x)
      if (y.le.4.0) then
        besj1 = 0.
        if (y.eq.0.0) return
        if (y.lt.xmin) 
     $  print *,' in besj1 abs(x) so small j1 underflows: x=',x
        if (y.gt.xmin) besj1 = x/2.0
        if (y.gt.xsml) besj1 = x * (.25d0+0 +
     $  csevl(y*y/8.0-1.0, bj1cs, ntj1))
c
      elseif (y.lt.xmax) then
        z = 32.0/y**2 - 1.0
        ampl = (0.75d+0 + csevl (z, bm1cs, ntm1)) / sqrt(y)
        theta = y - 3.0*pi4 + csevl (z, bth1cs, ntth1) / y
        besj1 = sign (ampl, x) * cos (theta)
      else
        print *,' in besj1 no precision because abs(x) is big: x=',x
      endif
c
      return
      end
c_______________________________________________________________________
      double precision function besy0 (x)
c$$$$ calls csevl, besj0
      implicit double precision (a-h, o-z)
      dimension by0cs(13), bm0cs(21), bth0cs(24)
      data by0 cs /   -.0112778393 92865573d+0,
     $ -.1283452375 6042035d+0,  -.1043788479 9794249d+0,
     $  .0236627491 83969695d+0, -.0020903916 47700486d+0,
     $  .0001039754 53939057d+0, -.0000033697 47162423d+0,
     $  .0000000772 93842676d+0, -.0000000013 24976772d+0,
     $  .0000000000 17648232d+0, -.0000000000 00188105d+0,
     $  .0000000000 00001641d+0, -.0000000000 00000011d+0/
      data bm0cs /    .0928496163 7381644d+0,
     $ -.0014298770 7403484d+0, .0000283057 9271257d+0,
     $ -.0000014330 0611424d+0, .0000001202 8628046d+0,
     $ -.0000000139 7113013d+0, .0000000020 4076188d+0,
     $ -.0000000003 5399669d+0, .0000000000 7024759d+0,
     $ -.0000000000 1554107d+0, .0000000000 0376226d+0,
     $ -.0000000000 0098282d+0, .0000000000 0027408d+0,
     $ -.0000000000 0008091d+0, .0000000000 0002511d+0,
     $ -.0000000000 0000814d+0, .0000000000 0000275d+0,
     $ -.0000000000 0000096d+0, .0000000000 0000034d+0,
     $ -.0000000000 0000012d+0, .0000000000 0000004d+0/
      data bth0cs /
     $ -.2463916377 4300119d+0,  .0017370983 07508963d+0,
     $ -.0000621836 33402968d+0, .0000043680 50165742d+0,
     $ -.0000004560 93019869d+0, .0000000621 97400101d+0,
     $ -.0000000103 00442889d+0, .0000000019 79526776d+0,
     $ -.0000000004 28198396d+0, .0000000001 02035840d+0,
     $ -.0000000000 26363898d+0, .0000000000 07297935d+0,
     $ -.0000000000 02144188d+0, .0000000000 00663693d+0,
     $ -.0000000000 00215126d+0, .0000000000 00072659d+0,
     $ -.0000000000 00025465d+0, .0000000000 00009229d+0,
     $ -.0000000000 00003448d+0, .0000000000 00001325d+0,
     $ -.0000000000 00000522d+0, .0000000000 00000210d+0,
     $ -.0000000000 00000087d+0, .0000000000 00000036d+0/
      data twodpi / 0.6366197723 6758134d+0 /
      data pi4 / 0.7853981633 9744831d+0 /
      data nty0, ntbm0, ntth0,    xsml,      xmax /
     $       13,    21,    24, 2.1074d-08, 4.5035d+15/
ce
c
      if (x.le.0.) then
        print *,' in besy0 x is zero or negative: x=',x
c
      elseif (x.le.4.0) then
        y = 0.
        if (x.gt.xsml) y = x*x
        besy0 = twodpi*log(x/2.0)*besj0(x) + .375d+0 +
     $    csevl (y/8.0-1.0, by0cs, nty0)
c
      elseif (x.lt.xmax) then
        z = 32.0/x**2 - 1.0
        ampl = (0.75d+0 + csevl (z, bm0cs, ntbm0)) / sqrt(x)
        theta = x - pi4 + csevl (z, bth0cs, ntth0) / x
        besy0 = ampl * sin (theta)
c
      else
        print *,' in besy0 no precision because x is big: x=',x
      endif
      return
      end
c_______________________________________________________________________
      double precision function besy1 (x)
      implicit double precision (a-h, o-z)
      dimension by1cs(14), bm1cs(21), bth1cs(24)
      data by1 cs /    
     $  .0320804710 0611908629d+0,   1.2627078974 33500450d+0,
     $  .0064999618 9992317500d+0, -.0893616452 8860504117d+0,
     $  .0132508812 2175709545d+0, -.0008979059 1196483523d+0,
     $  .0000364736 1487958306d+0, -.0000010013 7438166600d+0,
     $  .0000000199 4539657390d+0, -.0000000003 0230656018d+0,
     $  .0000000000 0360987815d+0, -.0000000000 0003487488d+0,
     $  .0000000000 0000027838d+0, -.0000000000 0000000186d+0/
      data bm1 cs /    .1047362510 931285d+0,
     $  .0044244389 3702345d+0, -.0000566163 9504035d+0,
     $  .0000023134 9417339d+0, -.0000001737 7182007d+0,
     $  .0000000189 3209930d+0, -.0000000026 5416023d+0,
     $  .0000000004 4740209d+0, -.0000000000 8691795d+0,
     $  .0000000000 1891492d+0, -.0000000000 0451884d+0,
     $  .0000000000 0116765d+0, -.0000000000 0032265d+0,
     $  .0000000000 0009450d+0, -.0000000000 0002913d+0,
     $  .0000000000 0000939d+0, -.0000000000 0000315d+0,
     $  .0000000000 0000109d+0, -.0000000000 0000039d+0,
     $  .0000000000 0000014d+0, -.0000000000 0000005d+0/
      data bth1cs /    
     $  .7406014102 6313850d+0,  -.0045717556 59637690d+0,
     $  .0001198185 10964326d+0, -.0000069645 61891648d+0,
     $  .0000006554 95621447d+0, -.0000000840 66228945d+0,
     $  .0000000133 76886564d+0, -.0000000024 99565654d+0,
     $  .0000000005 29495100d+0, -.0000000001 24135944d+0,
     $  .0000000000 31656485d+0, -.0000000000 08668640d+0,
     $  .0000000000 02523758d+0, -.0000000000 00775085d+0,
     $  .0000000000 00249527d+0, -.0000000000 00083773d+0,
     $  .0000000000 00029205d+0, -.0000000000 00010534d+0,
     $  .0000000000 00003919d+0, -.0000000000 00001500d+0,
     $  .0000000000 00000589d+0, -.0000000000 00000237d+0,
     $  .0000000000 00000097d+0, -.0000000000 00000040d+0/
      data twodpi / 0.6366197723 6758134d+0 /
      data pi4 / 0.7853981633 9744831d+0 /
      data nty1, ntm1, ntth1,        xmin,       xsml, xmax /
     $       14,   21,    24, 1.7654d-308, 2.1077d-08, 4.5035d+15/
ce
      if (x.le.0.) then
        print *,' in besy1 x is zero or negative: x=',x
c
      elseif (x.lt.xmin) then
        print *,' in besy1 x so small y1 overflows: x=',x
c
      elseif (x.le.4.0) then
c
        y = 0.
        if (x.gt.xsml) y = x*x
        besy1 = twodpi*log(x/2.0)*besj1(x) +
     $    (0.5d+0 + csevl (y/8.0-1.0, by1cs, nty1))/x
c
      elseif (x.lt.xmax) then
        z = 32.0/x**2 - 1.0
        ampl = (0.75d+0 + csevl (z, bm1cs, ntm1)) / sqrt(x)
        theta = x - 3.0*pi4 + csevl (z, bth1cs, ntth1) / x
        besy1 = ampl * sin (theta)
c
      else
        print *,' in besy1 no precision because x is big: x=',x
      endif
c
      return
      end
c_______________________________________________________________________
      double precision function besi0 (x)
c$$$$ calls csevl, besi0e
      implicit double precision (a-h, o-z)
      dimension bi0cs(12)
      data bi0 cs /
     $ -.0766054725 2839144951d+0, 1.9273379539 93808270d+0,
     $  .2282644586 920301339d+0, .0130489146 6707290428d+0,
     $ .0004344270 9008164874d+0, .0000094226 5768600193d+0,
     $ .0000001434 0062895106d+0, .0000000016 1384906966d+0,
     $ .0000000000 1396650044d+0, .0000000000 0009579451d+0,
     $ .0000000000 0000053339d+0, .0000000000 0000000245d+0/
      data nti0,        xsml,  xmax / 
     $       12, 2.1077d-08 ,  709.08/
ce
      y = abs(x)
c
      if (y.le.3.0) then
        besi0 = 1.0
        if (y.gt.xsml) besi0 = 2.75d+0 +
     $  csevl (y*y/4.5d+0-1.0, bi0cs, nti0)
      elseif (y.lt.xmax) then
        besi0 = exp(y) * besi0e(x)
c
      else
        print *,' in besi0 abs(x) so big i0 overflows: x=',x
      endif
      return
      end
      double precision function besi0e (x)
c$$$$ calls no other routines
      implicit double precision (a-h, o-z)
      dimension bi0cs(12), ai0cs(21), ai02cs(22)
      data bi0 cs /  
     $ -.0766054725 2839144951d+0, 1.9273379539 93808270d+0,
     $ .2282644586 920301339d+0,  .0130489146 6707290428d+0,
     $ .0004344270 9008164874d+0, .0000094226 5768600193d+0,
     $ .0000001434 0062895106d+0, .0000000016 1384906966d+0,
     $ .0000000000 1396650044d+0, .0000000000 0009579451d+0,
     $ .0000000000 0000053339d+0, .0000000000 0000000245d+0/
        data ai0 cs /           .0757599449 4023796d+0,
     $ .0075913808 1082334d+0,  .0004153131 3389237d+0,
     $ .0000107007 6463439d+0, -.0000079011 7997921d+0,
     $-.0000007826 1435014d+0,  .0000002783 8499429d+0,
     $ .0000000082 5247260d+0, -.0000000120 4463945d+0,
     $ .0000000015 5964859d+0,  .0000000002 2925563d+0,
     $-.0000000001 1916228d+0,  .0000000000 1757854d+0,
     $ .0000000000 0112822d+0, -.0000000000 0114684d+0,
     $ .0000000000 0027155d+0, -.0000000000 0002415d+0,
     $-.0000000000 0000608d+0,  .0000000000 0000314d+0,
     $-.0000000000 0000071d+0,  .0000000000 0000007d+0/
        data ai02cs / 
     $ .0544904110 1410882d+0,  .0033691164 7825569d+0,
     $ .0000688975 8346918d+0,  .0000028913 7052082d+0,
     $ .0000002048 9185893d+0,  .0000000226 6668991d+0,
     $ .0000000033 9623203d+0,  .0000000004 9406022d+0,
     $ .0000000000 1188914d+0, -.0000000000 3149915d+0,
     $-.0000000000 1321580d+0, -.0000000000 0179419d+0,
     $ .0000000000 0071801d+0,  .0000000000 0038529d+0,
     $ .0000000000 0001539d+0, -.0000000000 0004151d+0,
     $-.0000000000 0000954d+0,  .0000000000 0000382d+0,
     $ .0000000000 0000176d+0, -.0000000000 0000034d+0,
     $-.0000000000 0000027d+0,  .0000000000 0000003d+0/
      data nti0, ntai0, ntai02, xsml /
     $       12,    21,  22,   2.1077d-08/
ce
 10   y = abs(x)
      if (y.le.3.0) then
c
        besi0e = 1.0
        if (y.gt.xsml) besi0e = exp(-y) * ( 2.75d+0 +
     $  csevl (y*y/4.5d+0-1.0, bi0cs, nti0) )
c
      elseif (y.le.8.0) then
        besi0e = (.375d+0 + csevl ((48.0/y-11.)/5., ai0cs, ntai0)
     $  ) / sqrt(y)
c
      else
         besi0e = (.375d+0 + csevl (16./y-1., ai02cs, ntai02))
     $  / sqrt(y)
c
      endif
      return
      end
c_______________________________________________________________________
      double precision function besi1 (x)
c$$$$ calls csevl, besi1e
      implicit double precision (a-h, o-z)
      dimension bi1cs(11)
      data bi1 cs /             -.0019717132 61099859d+0,
     $  .4073488766 7546481d+0,  .0348389942 99959456d+0,
     $  .0015453945 56300123d+0, .0000418885 21098377d+0,
     $  .0000007649 02676483d+0, .0000000100 42493924d+0,
     $  .0000000000 99322077d+0, .0000000000 00766380d+0,
     $  .0000000000 00004741d+0, .0000000000 00000024d+0/
      data nti1, xmin,           xsml, xmax / 
     $       11, 9.89d-324, 2.981d-08, 709.0/
ce
      y = abs(x)
      besi1 = 0.0
      if (y.eq.0.0) return
      if (y.le.3.0) then
c
        if (y.lt.xmin) 
     $  print *,' in besi1 abs(x) so small i1 underflows: x=',x

        if (y.gt.xmin) besi1 = 0.5*x
        if (y.gt.xsml) besi1 = x * (.875d+0 +
     $     csevl(y*y/4.5d+0 -1.0, bi1cs, nti1))
c
      elseif (y.lt.xmax) then
        besi1 = exp(y) * besi1e(x)
      else
        print *,' in besi1 abs(x) so big i1 overflows: x=',x
      endif
c
      return
      end
c_______________________________________________________________________
      double precision function besi1e (x)
c$$$$ calls csevl
      implicit double precision (a-h, o-z)
      dimension bi1cs(11), ai1cs(21), ai12cs(22)
      data bi1 cs /   -.0019717132 61099859d+0,
     $  .4073488766 7546481d+0,  .0348389942 99959456d+0,
     $  .0015453945 56300123d+0,  .0000418885 21098377d+0,
     $  .0000007649 02676483d+0,  .0000000100 42493924d+0,
     $  .0000000000 99322077d+0,  .0000000000 00766380d+0,
     $  .0000000000 00004741d+0,  .0000000000 00000024d+0/
      data ai1 cs /   -.0284674418 1881479d+0,
     $ -.0192295323 1443221d+0, -.0006115185 8579437d+0,
     $ -.0000206997 1253350d+0,  .0000085856 1914581d+0,
     $  .0000010494 9824671d+0, -.0000002918 3389184d+0,
     $ -.0000000155 9378146d+0,  .0000000131 8012367d+0,
     $ -.0000000014 4842341d+0, -.0000000002 9085122d+0,
     $  .0000000001 2663889d+0, -.0000000000 1664947d+0,
     $ -.0000000000 0166665d+0,  .0000000000 0124260d+0,
     $ -.0000000000 0027315d+0,  .0000000000 0002023d+0,
     $  .0000000000 0000730d+0, -.0000000000 0000333d+0,
     $  .0000000000 0000071d+0, -.0000000000 0000006d+0/
      data ai12cs /    
     $  .0285762350 1828014d+0, -.0097610974 9136147d+0,
     $ -.0001105889 3876263d+0, -.0000038825 6480887d+0,
     $ -.0000002512 2362377d+0, -.0000000263 1468847d+0,
     $ -.0000000038 3538039d+0, -.0000000005 5897433d+0,
     $ -.0000000000 1897495d+0,  .0000000000 3252602d+0,
     $  .0000000000 1412580d+0,  .0000000000 0203564d+0,
     $ -.0000000000 0071985d+0, -.0000000000 0040836d+0,
     $ -.0000000000 0002101d+0,  .0000000000 0004273d+0,
     $  .0000000000 0001041d+0, -.0000000000 0000382d+0,
     $ -.0000000000 0000186d+0,  .0000000000 0000033d+0,
     $  .0000000000 0000028d+0, -.0000000000 0000003d+0/
      data nti1, ntai1, ntai12, xmin,         xsml / 
     $      11,     21,     22, 9.882d-324, 2.9803d-08/
ce
c
 10   y = abs(x)
      if (y.le.3.0) then
        besi1e = 0.0
        if (y.eq.0.0)  return
c
        if (y.lt.xmin) 
     $  print *,' in besi1e abs(x) so small i1 underflows: x=',x
        if (y.gt.xmin) besi1e = x/2.0
        if (y.gt.xsml) besi1e = x * (.875d+0 +
     $    csevl(y*y/4.5d+0-1., bi1cs, nti1))
        besi1e = exp(-y) * besi1e
c
      elseif (y.le.8.) then
        besi1e = (.375d+0 + csevl ((48./y-11.)/5., ai1cs, ntai1))/
     $   sqrt(y)
c
      else
        besi1e = (.375d+0 + csevl (16./y-1.0, ai12cs, ntai12))/
     $   sqrt(y)
      endif
      besi1e = sign (besi1e, x)
      return
      end
c_______________________________________________________________________
      double precision function besk0 (x)
c$$$$ calls csevl, besk0e, besi0
      implicit double precision (a-h, o-z)
      dimension bk0cs(11)
      data bk0cs /              - .0353273932 3390276872d+0,
     $ .3442898999 246284869d+0,  .0359799365 1536150163d+0,
     $ .0012646154 1144692592d+0, .0000228621 2103119451d+0,
     $ .0000002534 7910790261d+0, .0000000019 0451637722d+0,
     $ .0000000000 1034969525d+0, .0000000000 0004259816d+0,
     $ .0000000000 0000013744d+0, .0000000000 0000000035d+0/
c
      data ntk0,       xsml,      xmax / 
     $       11,     2.107d-08, 741.126/
ce
      if (x.le.0.) then
        print *,' in besk0 x is zero or negative: x=',x
      elseif (x .le. 2.0) then
c
        y = 0.
        if (x.gt.xsml) y = x*x
        besk0 = -log(x/2.0)*besi0(x) - .25d+0 +
     $    csevl (y/2.0-1.0, bk0cs, ntk0)
c
      elseif (x.le.xmax) then
        besk0 = exp(-x) * besk0e(x)
c
      else
        besk0 = 0.
        print *,' in besk0 x so big k0 underflows: x=',x
c
      endif
      return
      end
c_______________________________________________________________________
      double precision function besk0e (x)
c$$$$ calls csevl, besi0
      implicit double precision (a-h, o-z)
c
      dimension bk0cs(11), ak0cs(17), ak02cs(14)
      data bk0 cs /             -.0353273932 3390276872d+0,
     $ .3442898999 246284869d+0, .0359799365 1536150163d+0,
     $ .0012646154 1144692592d+0, .0000228621 2103119451d+0,
     $ .0000002534 7910790261d+0, .0000000019 0451637722d+0,
     $ .0000000000 1034969525d+0, .0000000000 0004259816d+0,
     $ .0000000000 0000013744d+0, .0000000000 0000000035d+0/
      data ak0 cs /            -.0764394790 3327941d+0,
     $ -.0223565260 5699819d+0, .0007734181 1546938d+0,
     $ -.0000428100 6688886d+0, .0000030817 0017386d+0,
     $ -.0000002639 3672220d+0, .0000000256 3713036d+0,
     $ -.0000000027 4270554d+0, .0000000003 1694296d+0,
     $ -.0000000000 3902353d+0, .0000000000 0506804d+0,
     $ -.0000000000 0068895d+0, .0000000000 0009744d+0,
     $ -.0000000000 0001427d+0, .0000000000 0000215d+0,
     $ -.0000000000 0000033d+0, .0000000000 0000005d+0/
      data ak02cs /
     $ -.0120186982 6307592d+0, -.0091748526 9102569d+0,
     $ .0001444550 9317750d+0, -.0000040136 1417543d+0,
     $ .0000001567 8318108d+0, -.0000000077 7011043d+0,
     $ .0000000004 6111825d+0, -.0000000000 3158592d+0,
     $ .0000000000 0243501d+0, -.0000000000 0020743d+0,
     $ .0000000000 0001925d+0, -.0000000000 0000192d+0,
     $ .0000000000 0000020d+0, -.0000000000 0000002d+0/
      data ntk0, ntak0, ntak02, xsml / 
     $       11,    17,     14, 2.107d-08/
ce
      if (x.le.0.) then
        print *,' in besk0e x is zero or negative: x=',x
c
      elseif (x.le.2.0) then
        y = 0.
        if (x.gt.xsml) y = x*x
        besk0e = exp(x) * (-log(x/2.0)*besi0(x)
     $    - .25d+0 + csevl (y/2.0-1.0, bk0cs, ntk0) )
c
      elseif (x.le.8.0) then
        besk0e = (1.25d+0 + csevl ((16./x-5.)/3., ak0cs, ntak0))/
     $  sqrt(x)
c
      elseif (x.gt.8.0) then
        besk0e = (1.25d+0 + csevl (16./x-1., ak02cs, ntak02))/
     $  sqrt(x)
c
      endif
      return
      end
c_______________________________________________________________________
      double precision function besk1 (x)
c$$$$ calls csevl, besk1e
      implicit double precision (a-h, o-z)
      dimension bk1cs(11)
      data bk1 cs /    .0253002273 389477705d+0,
     $ -.3531559607 76544876d+0,   -.1226111808 22657148d+0,
     $ -.0069757238 596398643d+0, -.0001730288 957513052d+0,
     $ -.0000024334 061415659d+0, -.0000000221 338763073d+0,
     $ -.0000000001 411488392d+0, -.0000000000 006666901d+0,
     $ -.0000000000 000024274d+0, -.0000000000 000000070d+0/
      data ntk1,       xmin,       xsml,   xmax /
     $       11, 1.124d-308, 2.1074d-08,  741.13/
ce
      if (x.le.0.) then
        print *,' in besk1 x is zero or negative: x=',x
      elseif (x.le.2.0) then
c
        if (x.lt.xmin) 
     $  print *,' in besk1 x so small k1 overflows: x=',x
        y = 0.
        if (x.gt.xsml) y = x*x
        besk1 = log(x/2.0)*besi1(x) +
     $  (0.75d+0 + csevl (y/2.0-1.0, bk1cs, ntk1))/x
c
      elseif (x.lt.xmax) then
        besk1 = exp(-x) * besk1e(x)
c
      else
        besk1 = 0.
        print *,' in besk1 x so big k1 underflows: x=',x
      endif
      return
      end
c_______________________________________________________________________
      double precision function besk1e (x)
c$$$$ calls csevl, besi1
      implicit double precision (a-h, o-z)
      dimension bk1cs(11), ak1cs(17), ak12cs(14)
      data bk1 cs /    .0253002273 389477705d+0,
     $ -.3531559607 76544876d+0, -.1226111808 22657148d+0,
     $ -.0069757238 596398643d+0, -.0001730288 957513052d+0,
     $ -.0000024334 061415659d+0, -.0000000221 338763073d+0,
     $ -.0000000001 411488392d+0, -.0000000000 006666901d+0,
     $ -.0000000000 000024274d+0, -.0000000000 000000070d+0/
      data ak1 cs /    .2744313406 973883d+0,
     $  .0757198995 3199368d+0, -.0014410515 5647540d+0,
     $  .0000665011 6955125d+0, -.0000043699 8470952d+0,
     $  .0000003540 2774997d+0, -.0000000331 1163779d+0,
     $  .0000000034 4597758d+0, -.0000000003 8989323d+0,
     $  .0000000000 4720819d+0, -.0000000000 0604783d+0,
     $  .0000000000 0081284d+0, -.0000000000 0011386d+0,
     $  .0000000000 0001654d+0, -.0000000000 0000248d+0,
     $  .0000000000 0000038d+0, -.0000000000 0000006d+0/
      data ak12cs /    
     $  .0637930834 3739001d+0, .0283288781 3049721d+0,
     $ -.0002475370 6739052d+0, .0000057719 7245160d+0,
     $ -.0000002068 9392195d+0, .0000000097 3998344d+0,
     $ -.0000000005 5853361d+0, .0000000000 3732996d+0,
     $ -.0000000000 0282505d+0, .0000000000 0023720d+0,
     $ -.0000000000 0002176d+0, .0000000000 0000215d+0,
     $ -.0000000000 0000022d+0, .0000000000 0000002d+0/
      data ntk1, ntak1, ntak12, xmin, xsml /
     $       11,    17,     14, 1.1238d-308, 2.1074d-08/
ce
      if (x.le.0.) then
        print *,' in besk1e x is zero or negative: x=',x
c
      elseif (x.le.2.0) then
        if (x.lt.xmin) 
     $  print *,' in besk1e x so small k1 overflows: x=',x
        y = 0.
        if (x.gt.xsml) y = x*x
        besk1e = exp(x) * (log(x/2.0)*besi1(x) +
     $  (0.75d+0 + csevl (y/2.0-1.0, bk1cs, ntk1))/x )
c
      elseif (x.le.8.)  then
        besk1e = (1.25d+0 + csevl ((16./x-5.)/3., ak1cs, ntak1))/
     $   sqrt(x)
      else 
         besk1e = (1.25d+0 + csevl (16./x-1., ak12cs, ntak12))/
     $   sqrt(x)
      endif
      return
      end
c_______________________________________________________________________
      function csevl (x, cs, n)
c$$$$ calls no other routines
      implicit double precision (a-h, o-z)
c   evaluate the n-term chebyshev series cs at x.
c
c evaluate the n-term chebyshev series cs at x. adapted from
c r broucke, algorithm 446, c.a.c.m., 16, 254 (1973). also see fox
c and parker, chebyshev polys in numerical analysis, oxford press, p.56.
c
c       input arguments --
c x    value at which the series is to be evaluated.
c cs   array of n terms of a chebyshev series. in evaluating cs,
c      only half the first coef is summed.
c n    number of terms in array cs.
c
c
       dimension cs(*)
c***first executable statement  csevl
c
       b1=0.
       b0=0.
       twox=2.0*x
       do 10 i=1,n
         b2=b1
         b1=b0
         b0=twox*b1-b2+cs(n+1-i)
 10    continue
c
       csevl = (b0-b2)/2.0
c
       return
       end
c______________________________________________________________________
