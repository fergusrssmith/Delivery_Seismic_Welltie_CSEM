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


      subroutine hy2fht(b,nb,nrel,tol,ntol,rerr,aerr,nord,fun1,rel1,
     *  zwork,fwork,bwork,npcs, ipath,maxdim,maxrel, zans,nfun1,ierr)
c=======================================================================
c     by w.l.anderson, u.s.geological survey, denver, colorado, usa.
c=======================================================================
c
c     hy2fht - a hybrid fast hankel transform subprogram for rapid and
c     accurate evaluation of hankel transforms of orders 0 and/or 1.
c     sub-methods within hybfht are related or lagged convolutions using
c     anderson's (1982) digital filtering subprogram (hankl2), and/or
c     chave's (1983) direct quadrature subprogram (besax2), where the
c     complex kernel functions can be highly oscillatory as well as
c     decreasing, increasing, or of bounded variation.
c
c     subprogram hy2fht has several options that the user may select
c     to control how the algorithm performs in speed and/or accuracy
c     (see definition of parameter ipath).  the idea was to use two
c     different algorithms in such a way as to take advantage of the
c     best features of each algorithm for total performance, and to
c     solve a wider class of problems involving hankel transforms.
c
c     note that the original sub-method subprograms have been renamed
c     and somewhat modified to coexist under control of hy2fht in
c     double-precision complex*16 arithmetic.  (for a sinlge-precision
c     complex version, see subroutine hybfht.)
c
c=======================================================================
c
c  references:
c
c  anderson, w.l., 1982, fast hankel transforms using related and lagged
c    convolutions: acm trans. on math. software, vol.8, no.4, p.344-368.
c
c  --------, 1984, discussion on: "numerical integration of related
c    hankel transforms by quadrature and continued fraction expansion"
c    (by a.d. chave), with reply by author: geophysics, vol.49, no.10,
c    p.1811-1813.
c
c  chave, a.d., 1983, numerical integration of related hankel transforms
c    by quadrature and continued fraction expansion: geophysics, vol.48,
c    no.12, p.1671-1686.
c
c=======================================================================
c
c  input parameters (use of "real" implies "real*8" in comments; also
c                    use of "complex" implies "complex*16" below):
c                    -------------------------------------------
c      b() - (real) array of hankel transform arguments b(i), i=1,..,nb.
c            for only a single argument b(1), use nb=1; nb>1 implies
c            subroutines hankl2 and/or besax2 to be invoked as required
c            using the specific b-arguments as listed in the input array
c            b(i), where b(i)>b(i-1)>0.0 for all i=2,...,nb. note that
c            nb>1 option is further controlled by parameter ipath below.
c            b() is assumed dimensioned b(maxdim), where maxdim is
c            defined in hy2fht as an adjustable dimension.  when hankl2
c            is used (ipath=+-1), array b() is used to interpolate
c            a cubic spline after lagged convolutions are performed;
c            when besax2 is used (ipath=+-2), each b() is used directly
c            to obtain results via a separate call to besax2.
c       nb - (integer) number of b's given in array b(i), i=1,nb>=1.
c            note that nb is the exact number of increasing arguments,
c            and is not the number of lagged convolutions (mb) as
c            defined in subroutine hankl.  in general, 1<=nb<maxdim
c            is required, with the further requirement that
c            mb = dint(10.*dlog(b(nb)/b(1)))+5<=maxdim if ipath=+-1.
c     nrel - (integer) number of related kernel functions supplied in
c            user written subroutine "rel1", which gives the nrel
c            relations as a function of "fun1".  for unrelated use,
c            set nrel=1 (see coding requirements of rel1 below.)
c            1<=nrel<=maxrel is required (see maxrel below).
c      tol - (real) tolerance factor passed to subroutine hankl2, as
c            required, depending on parameter ipath selected.  normally,
c            tol=.1d-14 is the suggested tolerance that should be used
c            to yield the approximate relative error .1d-12 from hankl2
c            using 32-bit floating-point words in double-precision.
c            tol=0.0 can be set for the maximum accuracy possible.
c            (see more information on "tol" in subroutine hankl2.)
c     ntol - (integer) number of consecutive times the tolerance is met
c            in subroutine hankl2, as required, depending on parameter
c            ipath selected.  normally, for smooth non-oscillatory
c            kernels, ntol=1 is sufficient, and ntol=2 or 3 possibly for
c            oscillating kernels to avoid a premature cut-off of the
c            adaptive convolution algorithm in hankl2.
c            (see more information on "ntol" in subroutine hankl2.)
c     rerr - (real) requested relative error passed to subroutine besax2
c            as required, depending on parameter ipath selected.  rerr
c            should not be chosen much smaller than about .1d-12 because
c            this is a 32bit floating-point complex*16 version of besax2
c            to coexist with hankl2 in double-precision complex*16.
c            (see more information on "rerr" in subroutine besax2.)
c     aerr - (real) requested absolute error passed to subroutine besax2
c            as required, depending on parameter ipath selected.  aerr
c            should not be chosen too close to zero, but this depends
c            on the particular "fun1" and/or "rel1" being used.  a good
c            rule of thumb is to generally choose aerr<=1.d-2*rerr.
c            (see more information on "aerr" in subroutine besax2, and
c            also in chave, 1983, p.1673.)
c   nord() - (integer) array of hankel transform orders nord(i), i=1,
c            nrel>=1.  note that nord(i)=0 or 1 is always required for
c            i=1,nrel. only orders 0 or 1 may be used in hy2fht, in any
c            combination; however, not all nrel transforms can be found
c            as related transforms if orders are unequal and when ipath
c            is selected such that subroutine besax2 is called.  this
c            is not the case if subroutine hankl2 is called, since the
c            digital filters can obtain related transforms of either
c            orders 0 or 1 in any nord() sequence desired.
c            nord() is assumed dimensioned nord(maxrel), where maxrel
c            is defined in hy2fht as an adjustable dimension.
c     fun1 - name of an external declared complex*16 function written by
c            the user to evaluate fun1(x) for x>0.0.  an external
c            statement must appear in the user's calling program as:
c                  external fun1
c            and the complex*16 function must have the general form:
c            (replace "<...>" below as required.)
c
c                  complex*16 function fun1(x)
c                  implicit real*8 (a-h,o-z)
c                  common/pass/<...include parameters as needed...>
c                  <...>
c                  fun1=dcmplx(<...>,<...>)
c                  return
c                  end
c
c            notes:  if parameters are needed in fun1 (for example a,c),
c            then use "common/pass/a,c" in calling program and fun1.
c            fun1 may be oscillating, increasing, decreasing, or of
c            bounded variation as the argument x>0 increases. (see more
c            information on "fun1" in subroutines hankl2 and besax2.)
c     rel1 - name of an external declared subroutine written by the user
c            to evaluate any (and all) nrel related kernels fj(nrel)
c            as a function of f1=fun1(x), assumed previously evaluated.
c            an external statement must appear in the user's program as:
c                  external rel1
c            and the subroutine must have the general form:
c            (replace "<...>" below as required.)
c
c                  subroutine rel1(x,f1,fj,j)
c                  implicit real*8 (a-h,o-z)
c                  complex*16 f1,fj(1)
c                  common/pass/<...include parameters as needed...>
c                  if(j.lt.0) go to (1,2,<...up to nrel...>), iabs(j)
c            1     fj(1)=f1
c                    if(j.lt.0) return
c            2     fj(2)=  <...expression in terms of x,f1, etc...>
c                    if(j.lt.0) return
c                  <...>
c                  return
c                  end
c
c        definition of rel1(x,f1,fj,j) parameters:
c        ----------------------------------------
c         x = real*8 input variable (controlled by hy2fht).
c        f1 = complex*16 input function fun1(x) value given by hy2fht.
c       fj()= complex*16 output array of related kernel functions, which
c       ****  the user must code;  note that each fj(j) may or may not
c             be a simple function of x, f1 and/or other parameters as
c             passed in common/pass/ from the user's calling program.
c         j = integer input variable (controlled by subroutine hy2fht).
c             option: when j<0, then only the iabs(j) related kernel
c             is returned in fj(iabs(j)); when j>0, all related kernels
c             fj(i), i=1,j are returned, where fj(1)=f1 always required.
c             (see rel1 example included with the test driver of hy2fht)
c  zwork() - (complex*16) array of dimension zwork(801,maxrel), where
c            maxrel is defined in hy2fht as an adjustable dimension.
c            zwork() is only used by subroutine hankl2, as required,
c            depending on parameter ipath selection.
c  fwork() - (complex*16) array of dimension fwork(maxrel), where
c            maxrel is defined in hy2fht as an adjustable dimension.
c            fwork() is used by either subroutines hankl2 or hy2fht,
c            as required, depending on parameter ipath selection.
c  bwork() - (real*8) array of dimension bwork(maxdim,9), where maxdim
c            is defined in hy2fht as an adjustable dimension.  bwork()
c            is only used whenever subroutine hankl2 is called, as
c            required, depending on parameter ipath selection.
c     npcs - (integer) number of pieces for partial integration as
c            used by subroutine besax2, as required, depending on
c            parameter ipath selection.  normally set npcs=1, but npcs
c            can be set greater than 1 for difficult kernels.
c            (see more information on "npcs" in subroutine besax2.)
c    ipath - (integer) option control path in subroutine hy2fht defined
c            as follows:      ---------------------------------
c        = 1 to force use of only subroutine hankl2, and as controlled
c            by input parameters above related to subroutine hankl2.
c        = 2 to force use of only subroutine besax2, and as controlled
c            by input parameters above related to subroutine besax2.
c       = -1 to begin with subroutine hankl2, but check during adaptive
c            convolution if a highly oscillatory kernel is being used
c            such that the oscillations exceed the filter's nyquist
c            frequency as transformed to log-space, i.e., 1/(2*0.1).
c            if this occurs, then hy2fht will automatically abort hankl2
c            and switch to subroutine besax2.  the number of function
c            fun1 evaluations will be noted in output array nfun1(1) for
c            hankl2 (until aborted) and nfun1(2) at exit.  if nfun1(2)=0
c            at exit, then no abort occurred, and subroutine hankl2 was
c            only used.
c       = -2 perform a relative error check using both subroutines
c            hankl2 & besax2 at the end points b(1) and b(nb) only when
c            nb>2 is selected.  the assumption is that besax2 is more
c            accurate, in general, and hence can be used to determine
c            if hankl2 is sufficiently accurate for lagged convolution
c            to be performed for the entire range of b(i),i=1,nb>2.
c            if the relative errors at b(1) or b(nb) are less than or
c            equal to parameter rerr*10, then subroutine hankl2 is used;
c            else subroutine besax2 is used.  the number of function
c            fun1 evaluations will be noted in output array nfun1(1)
c            and nfun1(2), where some extra evaluations are sometimes
c            necessary to perform the relative error tests, and when
c            hankl2 is followed to completion.  however, if nb>>2,
c            then many fun1 evaluations are saved with hankl2 than with
c            besax2 due to lagged convolutions.  thus ipath=-2 is safe
c            and usually not time-consuming in most cases when nb>>2.
c            (if ipath=-2 and nb<=2, then only besax2 is used.)
c   maxdim - (integer) maximum dimension used in program that calls
c            subroutine hy2fht for arrays: b(maxdim),bwork(maxdim,9),
c            zans(maxdim,maxrel).  nb<=maxdim is required if ipath=+-2
c            and mb=dint(10.*dlog(b(nb)/b(1)))+5<=maxdim if ipath=+-1.
c            it is recommended the latter relation always be used such
c            that maxdim>>nb. (this will allow for any ipath option.)
c   maxrel - (integer) maximum dimension used in program that calls
c            subroutine hy2fht for arrays: zwork(801,maxrel),
c            zans(maxdim,maxrel), and fwork(maxrel).  note that
c            nrel<=maxrel is required.
c
c-----------------------------------------------------------------------
c
c  output parameters:
c
c   zans() - (complex*16) array of hankel transform results zans(i,j),
c            i=1,nb and j=1,nrel.  the complex*16 array is assumed to be
c            dimensioned in the calling program as: zans(maxdim,maxrel).
c            note that 1<=nb<=maxdim and 1<=nrel<=maxrel are required.
c  nfun1() - (integer) array of number of function fun1 evaluations
c            taken by either subroutine hankl2 (stored in nfun1(1))
c            and/or subroutine besax2 (stored in nfun1(2)).  usually,
c            nfun1(1) is not much larger than the number of weights
c            needed for a single direct convolution (for any nb, nrel).
c            with regard to ipath selected, nfun1(2) is generally larger
c            than nfun1(1), but this depends on parameters tol, rerr,
c            aerr selected.  note that nfun1() reflects the outcome of
c            parameter ipath, and as a result of the fun1 and rel1
c            subprograms used in a call to subroutine hy2fht.
c     ierr - (integer) error return code defined as follows:
c        = 0 no error detected in input parameters.
c        = 1 error return: nb<1 or nb>maxdim not allowed.
c        = 2 error return: nrel<1 or nrel>maxrel not allowed.
c        = 3 error return: b(i-1)<=0 or b(i)<=b(i-1), i=2,nb.
c            (array b(i)=0.0 or not increasing for i=1,nb.)
c        = 4 error ipath=+-1: mb=dint(10.*dlog(b(nb)/b(1)))+5>maxdim.
c            (mb too large, or maxdim too small has been selected.)
c        = 5 error ipath=+-1: dspln1 error return after call dspln1.
c            (this should never occur in practice due to prior checks.)
c        = 6 error return: ipath<-2, ipath=0, or ipath>2 not allowed.
c       = 11 error ipath=+-1 or -2: error return after call hankl2.
c            (see subroutine hankl2 ierr=1 causes and corrective action)
c       = 21 error ipath=2: error return after call besax2.
c            (see subroutine besax2 ierr=1 causes; possible correction
c             is to try using npcs>1 and increase rerr, aerr.)
c       = 31 error ipath=-1 or -2: error return after call besax2.
c            (see subroutine besax2 ierr=1 causes; possible correction
c             is to try using npcs>1 and increase rerr, aerr.)
c
c=======================================================================
c
c  hy2fht usage:  (replace "<...>" below as required.)
c
c    the user's program should contain the following statements:
c    ----------------------------------------------------------
c     implicit real*8 (a-h,o-z)
c     parameter (maxdim= <...>,maxrel= <...>)
c     integer nord(maxrel),nfun1(2)
c     real*8 b(maxdim),bwork(maxdim,9)
c     complex*16 zwork(801,maxrel),zans(maxdim,maxrel),fwork(maxrel)
c     external fun1,rel1 <,...possibly include fun2,rel2, etc...>
c     common/pass/<...use to pass other parameters to fun1, etc...>
c     <...>
c     call hy2fht(<...see calling sequence above...>)
c     if(ierr.ne.0) stop <...>
c     <...>
c     end
c     complex*16 function fun1(x)
c     <...see form of fun1 above...>
c     end
c     subroutine rel1(x,f1,fj,j)
c     <...see form of rel1 above...>
c     end
c
c=======================================================================
c
      implicit real*8 (a-h,o-z)
      integer nord(maxrel),nfun1(2)
      real*8 b(maxdim),bwork(maxdim,9),der(2)
      complex*16 zwork(801,maxrel),zans(maxdim,maxrel),fwork(maxrel)
      double precision besr,besi
      external fun1,rel1
      common/contrl/jpath,j
      common/test/ng,nf2,ni
      data der/2*0.0d0/
c
c -af
c -af	    write (*,*) ' entering hy2fht, some argumts:'
c -af	    write (*,*) '  b,nb,nrel=',b,nb,nrel  
c -af	    write (*,*) '  tol,ntol,rerr,aerr=',tol,ntol,rerr,aerr 
c -af	    write (*,*) '  npcs,ipath=',npcs,ipath 
c -af	    write (*,*) '  maxdim,maxrel=',maxdim,maxrel 
c -af	    write (*,*) '  ' 
c -af
c  preset and test some parameters (others tested in subr. calls)
      ierr=0
      nfun1(1)=0
      nfun1(2)=0
      if(nb.lt.1.or.nb.gt.maxdim) then
        ierr=1
        return
      else if(nrel.lt.1.or.nrel.gt.maxrel) then
        ierr=2
        return
      else if (nb.gt.1) then
        do 5 i=2,nb
           if(b(i).le.b(i-1).or.b(i-1).le.0.0d0) then
c          if(b(i).lt.b(i-1).or.b(i-1).lt.0.0d0) then
            ierr=3
            return
          endif
5       continue
      endif
c  process via ipath selection (see ipath definitions)
      jpath=ipath
10    kpath=iabs(jpath)
c--beginning of "if(kpath.eq.1) then" loop:
      if(kpath.eq.1) then
        if(nb.gt.1) then
          bmax=b(nb)
          mb=dint(10.*dlog(bmax/b(1)))+5
          if(mb.gt.maxdim) then
            ierr=4
            return
          endif
        else
          bmax=b(1)
          mb=1
        endif
c  fht filter method (jpath=-1 or +1)
c  **********************************
        call hankl2(bmax,mb,nrel,tol,ntol,nord,fun1,rel1,zwork,fwork,
     1 maxdim,maxrel, zans,bwork(1,1),nf1,ier)
        nfun1(1)=nfun1(1)+nf1
        if(ier.ne.0) then
          ierr=ier+10
          return
        endif
c  jpath=-2 can be set automatically by hankl2 if ipath=-1 was set and
c  poor accuracy obtained by hankl2 due to highly oscillating kernels.
        if(jpath.eq.-2) then
          i1=1
          inb=nb
          go to 110
        endif
        if(nb.eq.1) return
c  nb>1, get spline interpolated zans(i,nrel) at b(i),i=1,nb.
        do 100 j=1,nrel
          k=8
          do 20 i=1,mb
            bwork(i,8)=0.0d0
            bwork(i,9)=0.0d0
            if(j.eq.1) bwork(i,1)=dlog(bwork(i,1))
            bwork(i,2)=dreal(zans(i,j))
            if(dimag(zans(i,j)).ne.0.0d0) k=-8
20        continue
30        kk=iabs(k)
          call dspln1(mb,0.0d0,bwork(1,1),bwork(1,2),bwork(1,3),
     1     bwork(1,4),bwork(1,5),0,der,bwork(1,6),bwork(1,7))
          if(mb.lt.0) then
            ierr=5
            return
          endif
          do 40 i=1,nb
            blog=dlog(b(i))
            if(blog.lt.bwork(1,1)) then
              bwork(i,kk)=bwork(i,2)
            else if(blog.gt.bwork(mb,1)) then
              bwork(i,kk)=bwork(mb,2)
            else
              call dsplnp(mb,bwork(1,1),bwork(1,2),bwork(1,3),
     1         bwork(1,4),bwork(1,5),blog,bwork(i,kk))
            endif
40        continue
          if(k.eq.-8) then
            k=9
            do 50 i=1,mb
              bwork(i,2)=dimag(zans(i,j))
50          continue
            go to 30
          else
            do 60 i=1,nb
              zans(i,j)=dcmplx(bwork(i,8),bwork(i,9))
60          continue
          endif
100     continue
c  done with spline interpolation when nb>1 and nrel>=1.
        return
c--rest of "if(kpath.eq.1) then" loop above:
      else if(kpath.eq.2) then
        i1=1
        inb=nb
        if(ipath.eq.2.or.nb.lt.3) go to 110
c  special preprocessing only when ipath=-2 and nb>2
        rerr10=10.*rerr
        do 105 i=nb,1,1-nb
          call hankl2(b(i),1,nrel,tol,ntol,nord,fun1,rel1,zwork,fwork,
     1   maxdim,maxrel, zans,bwork(1,1),nf1,ier)
          nfun1(1)=nfun1(1)+nf1
          if(ier.ne.0) then
            ierr=ier+10
            return
          endif
          do 105 j=1,nrel
            if(j.gt.1.and.nord(j).eq.nord(1)) then
              new=2
              call relat2(fwork,maxrel,rel1,j,0)
            else
              if(nrel.gt.1) then
                new=1
              else
                new=0
              endif
            endif
            testr=dreal(zans(1,j))
            testi=dimag(zans(1,j))
c  get assumed true values from besax2 at this (i,j),i=nb or 1, j=1,nrel
c  and compare with hankl2 test values at the end points b(nb), b(1).
c  switch to ipath=1 or 2 depending on relative error tests.
            call besax2(besr,besi,dfloat(nord(j)),1,7,b(i),fun1,
     1     rel1, fwork,maxrel,rerr,aerr,npcs,new,ier)
            nfun1(2)=nfun1(2)+nf2
            if(ier.ne.0) then
              ierr=ier+30
c             returns with ierr=31
            endif
            zans(i,j)=dcmplx(besr,besi)
            if(j.eq.1.and.nrel.gt.1) call relat2(fwork,maxrel,rel1,j,1)
            true=besr
            if(true.ne.0.0d0) then
              relerr=abs((testr-true)/true)
              if(relerr.gt.rerr10) then
                if(i.eq.nb) then
                  inb=nb-1
                else if(i.eq.1) then
                  i1=2
                  inb=nb-1
                endif
                jpath=2
                if(j.eq.nrel) go to 110
              endif
            endif
            true=besi
            if(true.ne.0.0d0) then
              relerr=abs((testi-true)/true)
              if(relerr.gt.rerr10) then
                if(i.eq.nb) then
                  inb=nb-1
                else if(i.eq.1) then
                  i1=2
                  inb=nb-1
                endif
                jpath=2
                if(j.eq.nrel) go to 110
              endif
            endif
105     continue
        jpath=1
        go to 10
c--rest of "if(kpath.eq.1) then" loop above:
      else
        ierr=6
        return
      endif
c  quadrature method (jpath=-2 or +2)
c  **********************************
110   do 200 i=i1,inb
      do 200 j=1,nrel
        if(j.gt.1.and.nord(j).eq.nord(1)) then
          new=2
          call relat2(fwork,maxrel,rel1,j,0)
        else
          if(nrel.gt.1) then
            new=1
          else
            new=0
          endif
        endif
        call besax2(besr,besi,dfloat(nord(j)),1,7,b(i),fun1,
     1 rel1, fwork,maxrel,rerr,aerr,npcs,new,ier)
        nfun1(2)=nfun1(2)+nf2
        if(ier.ne.0) then
          if(jpath.eq.-2) then
            ierr=ier+30
          else
            ierr=ier+20
          endif
c         returns with ierr>=21
        endif
        zans(i,j)=dcmplx(besr,besi)
        if(j.eq.1.and.nrel.gt.1) call relat2(fwork,maxrel,rel1,j,1)
200   continue
      return
      end
      subroutine relat2(fwork,maxrel,rel1,jrel,isave)
c----------------------------------------------------------------------
c  get related kernels as required for subroutine besaux via call rel1.
c  note all i/o parameters controlled by hybfht.  the user only needs
c  to code external subroutine rel1 to compute any (or all) related
c  kernels. (see comments about subroutine rel1 in hybfht.)
c----------------------------------------------------------------------
      implicit double precision (a-h,k,o-z)
c  nterm is maximum number of bessel function loops stored if new.ne.0
      parameter (nterm=50)
      save ksave,kasave,kbsave,nksave,npsave,npsav
      real*8 x,fre,fim
      dimension karg(255,nterm),kern(510,nterm),nk(nterm),
     1 kbes(510,nterm),ksave(510,nterm),kasave(255,nterm),
     2 kbsave(510,nterm),nksave(nterm)
      complex*16 fwork(maxrel)
      external rel1
      common/besint/nk,np,nps,karg,kern,kbes
      if(isave.eq.1) then
        npsav=np
        npsave=nps
c-af        do 1 i=1,min(nps,100)
        do 1 i=1,min(nps,nterm)
          nksave(i)=nk(i)
1       continue
      else
        np=npsav
        nps=npsave
c-af        do 2 i=1,min(nps,100)
        do 2 i=1,min(nps,nterm)
          nk(i)=nksave(i)
2       continue
      endif
c
c-af      do 10 i=1,min(nps,100)
      do 10 i=1,min(nps,nterm)
        lk=1
        do 10 j=1,2**(nk(i)+1)-1
          if(isave.eq.1) then
            tbes=kbes(j,i)
            if(tbes.ne.0.0) then
              rbes=1.0d0/tbes
            else
              rbes=0.0d0
            endif
            ksave(lk,i)=kern(lk,i)*rbes
            ksave(lk+1,i)=kern(lk+1,i)*rbes
            kasave(j,i)=karg(j,i)
            kbsave(j,i)=kbes(j,i)
            lk=lk+2
            go to 10
          endif
c--get saved argument, bessel, and kernel function value for f1.
          x=kasave(j,i)
          tbes=kbsave(j,i)
          fre=ksave(lk,i)
          fim=ksave(lk+1,i)
c--get related kernel from saved kernel as required via call rel1.
          call rel1(x,dcmplx(fre,fim),fwork,-jrel)
          kern(lk,i)=dreal(fwork(jrel))*tbes
          kern(lk+1,i)=dimag(fwork(jrel))*tbes
          lk=lk+2
10    continue
      return
      end
      subroutine besax2(besr,besi,order,nl,nu,r,fun1,rel1,fwork,maxrel,
     $  rerr,aerr, npcs,new,ierr)
c
c  performs automatic calculation of bessel transform to specified
c  relative and absolute error
c
c  argument list:
c
c  besr,besi-real and imaginary parts returned by besax2
c  order-order of the bessel function
c  nl-lower limit for gauss order to start computation
c  nu-upper limit for gauss order
c     nu,nl=1,...7 selects 3,7,15,31,63,127,and 255 point gauss
c     quadrature between the zero crossings of the bessel function
c  r-argument of the bessel function
c///////////////////////////////////////////////////////////////////////
c     the following parameters are included so that besax2 can be run
c     under the general control of subroutine hy2fht:
c///////////////////////////////////////////////////////////////////////
c     fun1 - name of an external declared complex*16 function written by
c            the user to evaluate fun1(g) for g>0.0.  an external
c            statement must appear in the user's calling program as:
c                  external fun1
c            and the complex*16 function must have the general form:
c            (replace "<...>" below as required.)
c
c                  complex*16 function fun1(g)
c                  implicit real*8 (a-h,o-z)
c                  common/pass/<...include parameters as needed...>
c                  <...>
c                  fun1=dcmplx(<...>,<...>)
c                  return
c                  end
c
c        notes: see comments in subroutines hy2fht and hankl2 on fun1.
c        the external funct call was changed from the form:
c        call fun1(x,yr,yi)  --  original method in besaut -- to
c        z=fun1(x),  where complex*16 fun1 used in besax2, bestx2, and
c        besqx2.  this method was used here to be compatible with
c        external complex*16 fun1 used in call hankl2 and hy2fht.
c     rel1 - name of an external declared subroutine written by the user
c            to evaluate any (and all) nrel related kernels fj(nrel)
c            as a function of f1=fun1(g), assumed previously evaluated.
c            an external statement must appear in the user's program as:
c                  external rel1
c            and the subroutine must have the general form:
c            (replace "<...>" below as required.)
c
c                  subroutine rel1(g,f1,fj,j)
c                  implicit real*8 (a-h,o-z)
c                  complex*16 f1,fj(1)
c                  common/pass/<...include parameters as needed...>
c                  if(j.lt.0) go to (1,2,<...up to nrel...>), iabs(j)
c            1     fj(1)=f1
c                    if(j.lt.0) return
c            2     fj(2)=  <...expression in terms of g,f1, etc...>
c                    if(j.lt.0) return
c                  <...>
c                  return
c                  end
c           notes:
c             see comments in subroutines hy2fht and hankl2 on rel1.
c             rel1 parameters g,f1,j are controlled by subroutine besax2
c             where the user must only code the fj(j),j=1,nrel related
c             kernels as a function of g,f1, and any other parameters
c             that can be passed in common/pass/, as required.
c     fwork - complex*16 work array fwork(maxrel), which is used in call
c             to rel1 and passed to besax2 by user's calling program.
c    maxrel - maximum dimension of fwork(maxrel) in calling program.
c///////////////////////////////////////////////////////////////////////
c  rerr,aerr-relative and absolute error for termination
c       besax2 terminates when increasing the gauss order does not
c       change the result by more than rerr or when the absolute error
c       is less than aerr or when a gauss order of nu is reached.
c  npcs-number of pieces into which each partial integrand is divided,
c       ordinarily set to one. for very small values of r where
c       the kernel function is appreciable only over the first few
c       loops of the bessel function, npcs may be increased to achieve
c       reasonable accuracy.
c  new  if new=1, the integrands are computed and saved at each gauss
c       order. if new=2, previously computed integrands are used. note
c       that order,r, and npcs must not be changed when setting new=2.
c  ierr-error parameter
c       ierr=0--normal return
c       ierr=1--result not accurate to rerr due to too low a gauss
c       order or convergence not achieved in bestx2
c
c  subroutines required:
c    bestx2,besqx2,bessj,padecf,cf,zeroj,dot
c
c  a.chave igpp/ucsd
c
      implicit double precision (a-h,o-z)
      dimension xsum(1)
      parameter (nsum=0)
c  common /test/ contains run statistics
c  ng=highest gauss order used
c  nf=total function evaluations  for all calls to bestx2
c  ni=total number of partial integrands on last call
c  to bestx2 only
      common /test/ng,nf,ni
c//
      external fun1,rel1
      complex*16 fwork(maxrel)
c//
      nf=0
      if(nl.gt.nu)then
        besr=0.
        besi=0.
        ierr=1
        return
      endif
      nw=max(new,1)
      call bestx2(besr,besi,order,nl,r,fun1,rel1,fwork,maxrel,
     $ .1*rerr,.1*aerr, npcs,xsum,nsum,nw,ierr)
      if((ierr.ne.0).and.(nl.eq.7))then
        ng=nl
        return
      else
        oldr=besr
        oldi=besi
        do 10 n=nl+1,nu
          call bestx2(besr,besi,order,n,r,fun1,rel1,fwork,maxrel,
     $     .1*rerr,.1*aerr, npcs,xsum,nsum,2,ierr)
          if((ierr.ne.0).and.(n.eq.7))then
            besr=oldr
            besi=oldi
            ng=n
            return
          elseif((abs(besr-oldr).le.rerr*abs(besr)+aerr).and.(abs(besi-
     $           oldi).le.rerr*abs(besi)+aerr))then
            ng=n
            return
          else
            oldr=besr
            oldi=besi
          endif
   10   continue
      endif
      ng=7
      ierr=1
      return
      end
      subroutine bestx2(besr,besi,order,ng,r,fun1,rel1,fwork,maxrel,
     $  rerr,aerr,npcs, xsum,nsum,new,ierr)
c
c  computes bessel transform of specified order defined as
c  integral(funct(x)*j-sub-order(x*r)*dx) from x=0 to infinity
c  computation is achieved by integration between the asymptotic
c  zero crossings of the bessel function using gauss quadrature.
c  the resulting series of partial integrands is summed by calculating
c  the pade approximants to speed up convergence.
c
c  argument list:
c
c  besr,besi  real and imaginary parts returned by bestx2
c  order  order of the bessel function
c  ng  number of gauss points to use in the quadrature routine.
c          ng=1 through 7 selects 3,7,15,31,63,126,and 255 terms.
c  r     argument of the bessel function
c///////////////////////////////////////////////////////////////////////
c   parameters "fun1,rel1,fwork,maxrel" are included so that besax2 can
c   be run under the general control of subroutine hy2fht.
c   (see subroutines besax2 and hy2fht for details on these parameters.)
c///////////////////////////////////////////////////////////////////////
c  rerr,aerr  specified relative and absolute error for the
c          calculation.  the integration
c          terminates when an additional term does not change the
c          result by more than rerr*result+aerr
c  npcs  number of pieces into which each partial integrand is divided,
c        ordinarily set to one. for very small values of range where
c        the kernel function is appreciable only over the first few
c        loops of the bessel function, npcs may be increased to achieve
c        reasonable accuracy. note that npcs affects only the pade
c        sum portion of the integration, over x(nsum) to infinity.
c  xsum  vector of values of the kernel argument of funct for which
c        explicit calculation of the integral is desired, so that the
c        integral over 0 to xsum(nsum) is added to the integral over
c        xsum(nsum) to infinity with the pade method invoked only for
c        the latter. this allows the pade summation method to be
c        overridden and some types of singularities to be handled.
c  nsum  number of values in xsum, may be zero.
c  new   determines method of kernel calculation
c        new=0 means  calculate but do not save integrands
c        new=1 means calculate kernel by calling funct-save kernel
c              times bessel function
c        new=2 means use saved kernels times bessel functions in
c              common /besint/. note that order,r,npcs,xsum, and
c              nsum may not be changed when setting new=2.
c  ierr  error parameter
c        0 normal return-integral converged
c        1 means no convergence after nstop terms in the pade sum
c
c  subroutines required:
c  besqx2,padecf,cf,zeroj,dot,bessj
c
c  a.chave igpp/ucsd
c
c  nterm is maximum number of bessel function loops stored if new.ne.0
c  nstop is maximum number of pade terms
c
      parameter (nterm=50, nstop=100)
      implicit double precision (a-h,o-z)
      double precision karg,kern,lastr,lasti,kbes
      dimension karg(255,nterm),kern(510,nterm),sr(nstop),si(nstop),
     $          nk(nterm),xsum(*),kbes(510,nterm)
c  common /test/ contains run statistics
c  ngauss=highest gauss order used
c  nf=total function evaluations for all calls to bestx2
c  ni=total number of partial integrands used in padecf on last call
c  to bestx2
      common /test/ngauss,nf,ni
c  if new.gt.0, common /besint/ contains the kernel function arguments
c  and values karg(i,n),kern(i,n) for n=1 to nps and i=1 to
c  2**(nk(n)+1)-1.
c  kern(i,n), kern(i+1,n) contain the real and imaginary parts of the
c  kernel function for the corresponding argument value karg(i,n).
c  each value of n corresponds to one entry in xsum or to one interval
c  between zero crossings of the bessel function divided by npcs.
c  note that up to nterm values are stored the integration proceeds
c  for n.gt.nterm as if new=0. np is an index passed internally by
c  bestx2.
      common /besint/nk,np,nps,karg,kern,kbes
c//
      external fun1,rel1
      complex*16 fwork(maxrel)
c//
c -af     
c -af     	    write (*,*) ' entering bestx2, some argumts:'
c -af	    write (*,*) '  order,ng,r=',order,ng,r  
c -af	    write (*,*) '  rerr,aerr=',rerr,aerr 
c -af	    write (*,*) '  npcs=',npcs 
c -af	    write (*,*) '  maxrel=',maxrel 
c -af	    write (*,*) '  ' 
c -af
c
      if(new.eq.2)then
        npo=nps
      else
        do 5 i=1,nterm
          nk(i)=0
    5   continue
        nps=0
        npo=nterm
      endif
c  check for trivial case
      if((order.ne.0.0).and.(r.eq.0.))then
        besr=0.
        besi=0.
        ierr=0
        return
      endif
      ni=0
      nw=new
      np=1
      npb=1
      l=1
      b=0.0
      sumr=0.0
      sumi=0.0
      xsumr=0.0
      xsumi=0.0
      if(nsum.gt.0)then
c  compute bessel transform explicitly on (0,xsum(nsum))
        lastr=0.0
        lasti=0.0
        do 10 n=1,nsum
          if((nw.eq.2).and.(np.gt.npo))nw=1
          if(np.gt.nterm)nw=0
          a=b
          b=xsum(n)
          call besqx2(a,b,termr,termi,ng,nw,order,r,fun1,rel1,
     1  fwork,maxrel)
          xsumr=xsumr+termr
          xsumi=xsumi+termi
          if((abs(xsumr-lastr).le.rerr*abs(xsumr)+aerr).and.
     $      (abs(xsumi-lasti).le.rerr*abs(xsumi)+aerr))then
            besr=xsumr
            besi=xsumi
            ierr=0
            nps=max(np,nps)
            return
          else
            np=np+1
            lastr=xsumr
            lasti=xsumi
          endif
   10   continue
c  find first zero crossing of bessel function beyond xsum(nsum)
   15   continue
        if(zeroj(npb,order).gt.xsum(nsum)*r)goto 20
        npb=npb+1
        goto 15
      endif
c  entry point for pade summation of partial integrands
   20 continue
      lastr=0.0
      lasti=0.0
      a=b
      b=zeroj(npb,order)/r
c  calculate terms and sum with padecf, quitting when convergence is
c  obtained
      do 30 n=1,nstop
        if(npcs.eq.1)then
          if((nw.eq.2).and.(np.gt.npo))nw=1
          if(np.gt.nterm)nw=0
          call besqx2(a,b,termr,termi,ng,nw,order,r,fun1,rel1,
     1  fwork,maxrel)
          np=np+1
        else
          termr=0.
          termi=0.
          xinc=(b-a)/npcs
          aa=a
          bb=a+xinc
          do 25 i=1,npcs
            if((nw.eq.2).and.(np.gt.npo))nw=1
            if(np.gt.nterm)nw=0
            call besqx2(aa,bb,tr,ti,ng,nw,order,r,fun1,rel1,
     1    fwork,maxrel)
            termr=termr+tr
            termi=termi+ti
            aa=bb
            bb=bb+xinc
            np=np+1
   25     continue
        endif
        ni=ni+1
        sr(l)=termr
        si(l)=termi
        call padecf(sumr,sumi,sr,si,l)
        if((abs(sumr-lastr).le.rerr*abs(sumr)+aerr).and.
     $     (abs(sumi-lasti).le.rerr*abs(sumi)+aerr))then
          besr=xsumr+sumr
          besi=xsumi+sumi
          ierr=0
          nps=max(np-1,nps)
          return
        else
          lastr=sumr
          lasti=sumi
          npb=npb+1
          a=b
          b=zeroj(npb,order)/r
          l=l+1
        endif
   30 continue
      besr=xsumr+sumr
      besi=xsumi+sumi
      ierr=1
      nps=max(np,nps)
      return
      end
      subroutine besqx2(a,b,besr,besi,ng,new,order,r,fun1,rel1,
     1 fwork,maxrel)
c
c  this subroutine calculates the integral of f(x)*j-sub-n(x*r)
c  over the interval a to b at a specified gauss order
c  the result is obtained using a sequence of 1,3,7,15,31,63,
c  127, and 255 point interlacing gauss formulae so that no integrand
c  evaluations are wasted. the kernel functions may be saved so that
c  bessel transforms of similar kernels are computed without new
c  evaluation of the kernel.
c  details on the formulae are given in 'the optimum addition of points
c  to quadrature formulae' by t.n.l. patterson, maths.comp.
c  22,847-856 (1968). gauss weights taken from
c  comm. a.c.m. 16,694-699 (1973)
c
c argument list:
c
c a      lower limit of integration
c b      upper limit of integration
c besr,besi returned integral value real and imaginary parts
c ng number of points in the gauss formula.  ng=1,...7
c        selects 3,7,15,31,63,127,and 255 point quadrature.
c new selects method of kernel evaluation
c     new=0 calculates kernels by calling f - nothing saved
c     new=1 calculates kernels by calling f and saves kernel times
c            bessel function in common /besint/
c     new=2 uses saved kernel times bessel functions in
c           common /besint/
c order order of the bessel function
c r    argument of the bessel function
c///////////////////////////////////////////////////////////////////////
c   parameters "fun1,rel1,fwork,maxrel" are included so that besax2 can
c   be run under the general control of subroutine hy2fht.
c   (see subroutines besax2 and hy2fht for details on these parameters.)
c///////////////////////////////////////////////////////////////////////
c//
c//  note:
c//  double-presision complex for besax2 at fun1 level only/
c//
c//  (see use of zfun1 and fun1 below as bounded by c// lines.)
c///////////////////////////////////////////////////////////////////////
c
c  a.chave igpp/ucsd
c
c  maximum number of bessel function loops that can be saved
      parameter (nterm=50)
      implicit double precision (a-h,o-z)
c//
      complex*16 zfun1,fun1,fwork(maxrel)
      common/contrl/jpath,jrel
      external fun1,rel1
c//
      double precision karg,kern,bessj,kbes
      dimension wt(254),wa(127),nwa(7),nwt(7),karg(255,nterm),
     $          kern(510,nterm),funct(254),fr1(64),fi1(64),fr2(64),
     $          fi2(64),bes1(64),bes2(64),nk(nterm),kbes(510,nterm)
c  common /test/ contains run statistics
c  ngauss=highest gauss order used
c  nf=total function evaluations  for all calls to bestx2
c  ni=total number of partial integrands on last call
c  to bestx2 only
      common /test/ngauss,nf,ni
      common /besint/nk,np,nps,karg,kern,kbes
c  starting addresses in arrays wt and wa at each order
      data nwt/1,3,7,15,31,63,127/,nwa/1,2,4,8,16,32,64/
c  array wt contains gauss quadrature weights stored as needed at
c  each order.  array wa contains the corresponding argument weights.
      data (wt(i),i=1,20)/
     $0.55555555555555555556d+00,0.88888888888888888889d+00,
     $0.26848808986833344073d+00,0.10465622602646726519d+00,
     $0.40139741477596222291d+00,0.45091653865847414235d+00,
     $0.13441525524378422036d+00,0.51603282997079739697d-01,
     $0.20062852937698902103d+00,0.17001719629940260339d-01,
     $0.92927195315124537686d-01,0.17151190913639138079d+00,
     $0.21915685840158749640d+00,0.22551049979820668739d+00,
     $0.67207754295990703540d-01,0.25807598096176653565d-01,
     $0.10031427861179557877d+00,0.84345657393211062463d-02,
     $0.46462893261757986541d-01,0.85755920049990351154d-01/
      data (wt(i),i=21,40)/
     $0.10957842105592463824d+00,0.25447807915618744154d-02,
     $0.16446049854387810934d-01,0.35957103307129322097d-01,
     $0.56979509494123357412d-01,0.76879620499003531043d-01,
     $0.93627109981264473617d-01,0.10566989358023480974d+00,
     $0.11195687302095345688d+00,0.11275525672076869161d+00,
     $0.33603877148207730542d-01,0.12903800100351265626d-01,
     $0.50157139305899537414d-01,0.42176304415588548391d-02,
     $0.23231446639910269443d-01,0.42877960025007734493d-01,
     $0.54789210527962865032d-01,0.12651565562300680114d-02,
     $0.82230079572359296693d-02,0.17978551568128270333d-01/
      data (wt(i),i=41,60)/
     $0.28489754745833548613d-01,0.38439810249455532039d-01,
     $0.46813554990628012403d-01,0.52834946790116519862d-01,
     $0.55978436510476319408d-01,0.36322148184553065969d-03,
     $0.25790497946856882724d-02,0.61155068221172463397d-02,
     $0.10498246909621321898d-01,0.15406750466559497802d-01,
     $0.20594233915912711149d-01,0.25869679327214746911d-01,
     $0.31073551111687964880d-01,0.36064432780782572640d-01,
     $0.40715510116944318934d-01,0.44914531653632197414d-01,
     $0.48564330406673198716d-01,0.51583253952048458777d-01,
     $0.53905499335266063927d-01,0.55481404356559363988d-01/
      data (wt(i),i=61,80)/
     $0.56277699831254301273d-01,0.56377628360384717388d-01,
     $0.16801938574103865271d-01,0.64519000501757369228d-02,
     $0.25078569652949768707d-01,0.21088152457266328793d-02,
     $0.11615723319955134727d-01,0.21438980012503867246d-01,
     $0.27394605263981432516d-01,0.63260731936263354422d-03,
     $0.41115039786546930472d-02,0.89892757840641357233d-02,
     $0.14244877372916774306d-01,0.19219905124727766019d-01,
     $0.23406777495314006201d-01,0.26417473395058259931d-01,
     $0.27989218255238159704d-01,0.18073956444538835782d-03,
     $0.12895240826104173921d-02,0.30577534101755311361d-02/
      data (wt(i),i=81,100)/
     $0.52491234548088591251d-02,0.77033752332797418482d-02,
     $0.10297116957956355524d-01,0.12934839663607373455d-01,
     $0.15536775555843982440d-01,0.18032216390391286320d-01,
     $0.20357755058472159467d-01,0.22457265826816098707d-01,
     $0.24282165203336599358d-01,0.25791626976024229388d-01,
     $0.26952749667633031963d-01,0.27740702178279681994d-01,
     $0.28138849915627150636d-01,0.50536095207862517625d-04,
     $0.37774664632698466027d-03,0.93836984854238150079d-03,
     $0.16811428654214699063d-02,0.25687649437940203731d-02,
     $0.35728927835172996494d-02,0.46710503721143217474d-02/
      data (wt(i),i=101,120)/
     $0.58434498758356395076d-02,0.70724899954335554680d-02,
     $0.83428387539681577056d-02,0.96411777297025366953d-02,
     $0.10955733387837901648d-01,0.12275830560082770087d-01,
     $0.13591571009765546790d-01,0.14893641664815182035d-01,
     $0.16173218729577719942d-01,0.17421930159464173747d-01,
     $0.18631848256138790186d-01,0.19795495048097499488d-01,
     $0.20905851445812023852d-01,0.21956366305317824939d-01,
     $0.22940964229387748761d-01,0.23854052106038540080d-01,
     $0.24690524744487676909d-01,0.25445769965464765813d-01,
     $0.26115673376706097680d-01,0.26696622927450359906d-01/
      data (wt(i),i=121,140)/
     $0.27185513229624791819d-01,0.27579749566481873035d-01,
     $0.27877251476613701609d-01,0.28076455793817246607d-01,
     $0.28176319033016602131d-01,0.28188814180192358694d-01,
     $0.84009692870519326354d-02,0.32259500250878684614d-02,
     $0.12539284826474884353d-01,0.10544076228633167722d-02,
     $0.58078616599775673635d-02,0.10719490006251933623d-01,
     $0.13697302631990716258d-01,0.31630366082226447689d-03,
     $0.20557519893273465236d-02,0.44946378920320678616d-02,
     $0.71224386864583871532d-02,0.96099525623638830097d-02,
     $0.11703388747657003101d-01,0.13208736697529129966d-01/
      data (wt(i),i=141,160)/
     $0.13994609127619079852d-01,0.90372734658751149261d-04,
     $0.64476204130572477933d-03,0.15288767050877655684d-02,
     $0.26245617274044295626d-02,0.38516876166398709241d-02,
     $0.51485584789781777618d-02,0.64674198318036867274d-02,
     $0.77683877779219912200d-02,0.90161081951956431600d-02,
     $0.10178877529236079733d-01,0.11228632913408049354d-01,
     $0.12141082601668299679d-01,0.12895813488012114694d-01,
     $0.13476374833816515982d-01,0.13870351089139840997d-01,
     $0.14069424957813575318d-01,0.25157870384280661489d-04,
     $0.18887326450650491366d-03,0.46918492424785040975d-03/
      data (wt(i),i=161,180)/
     $0.84057143271072246365d-03,0.12843824718970101768d-02,
     $0.17864463917586498247d-02,0.23355251860571608737d-02,
     $0.29217249379178197538d-02,0.35362449977167777340d-02,
     $0.41714193769840788528d-02,0.48205888648512683476d-02,
     $0.54778666939189508240d-02,0.61379152800413850435d-02,
     $0.67957855048827733948d-02,0.74468208324075910174d-02,
     $0.80866093647888599710d-02,0.87109650797320868736d-02,
     $0.93159241280693950932d-02,0.98977475240487497440d-02,
     $0.10452925722906011926d-01,0.10978183152658912470d-01,
     $0.11470482114693874380d-01,0.11927026053019270040d-01/
      data (wt(i),i=181,200)/
     $0.12345262372243838455d-01,0.12722884982732382906d-01,
     $0.13057836688353048840d-01,0.13348311463725179953d-01,
     $0.13592756614812395910d-01,0.13789874783240936517d-01,
     $0.13938625738306850804d-01,0.14038227896908623303d-01,
     $0.14088159516508301065d-01,0.69379364324108267170d-05,
     $0.53275293669780613125d-04,0.13575491094922871973d-03,
     $0.24921240048299729402d-03,0.38974528447328229322d-03,
     $0.55429531493037471492d-03,0.74028280424450333046d-03,
     $0.94536151685852538246d-03,0.11674841174299594077d-02,
     $0.14049079956551446427d-02,0.16561127281544526052d-02/
      data (wt(i),i=201,220)/
     $0.19197129710138724125d-02,0.21944069253638388388d-02,
     $0.24789582266575679307d-02,0.27721957645934509940d-02,
     $0.30730184347025783234d-02,0.33803979910869203823d-02,
     $0.36933779170256508183d-02,0.40110687240750233989d-02,
     $0.43326409680929828545d-02,0.46573172997568547773d-02,
     $0.49843645647655386012d-02,0.53130866051870565663d-02,
     $0.56428181013844441585d-02,0.59729195655081658049d-02,
     $0.63027734490857587172d-02,0.66317812429018878941d-02,
     $0.69593614093904229394d-02,0.72849479805538070639d-02,
     $0.76079896657190565832d-02,0.79279493342948491103d-02/
      data (wt(i),i=221,240)/
     $0.82443037630328680306d-02,0.85565435613076896192d-02,
     $0.88641732094824942641d-02,0.91667111635607884067d-02,
     $0.94636899938300652943d-02,0.97546565363174114611d-02,
     $0.10039172044056840798d-01,0.10316812330947621682d-01,
     $0.10587167904885197931d-01,0.10849844089337314099d-01,
     $0.11104461134006926537d-01,0.11350654315980596602d-01,
     $0.11588074033043952568d-01,0.11816385890830235763d-01,
     $0.12035270785279562630d-01,0.12244424981611985899d-01,
     $0.12443560190714035263d-01,0.12632403643542078765d-01,
     $0.12810698163877361967d-01,0.12978202239537399286d-01/
      data (wt(i),i=241,254)/
     $0.13134690091960152836d-01,0.13279951743930530650d-01,
     $0.13413793085110098513d-01,0.13536035934956213614d-01,
     $0.13646518102571291428d-01,0.13745093443001896632d-01,
     $0.13831631909506428676d-01,0.13906019601325461264d-01,
     $0.13968158806516938516d-01,0.14017968039456608810d-01,
     $0.14055382072649964277d-01,0.14080351962553661325d-01,
     $0.14092845069160408355d-01,0.14094407090096179347d-01/
      data (wa(i),i=1,20)/
     $0.77459666924148337704d+00,0.96049126870802028342d+00,
     $0.43424374934680255800d+00,0.99383196321275502221d+00,
     $0.88845923287225699889d+00,0.62110294673722640294d+00,
     $0.22338668642896688163d+00,0.99909812496766759766d+00,
     $0.98153114955374010687d+00,0.92965485742974005667d+00,
     $0.83672593816886873550d+00,0.70249620649152707861d+00,
     $0.53131974364437562397d+00,0.33113539325797683309d+00,
     $0.11248894313318662575d+00,0.99987288812035761194d+00,
     $0.99720625937222195908d+00,0.98868475754742947994d+00,
     $0.97218287474858179658d+00,0.94634285837340290515d+00/
      data (wa(i),i=21,40)/
     $0.91037115695700429250d+00,0.86390793819369047715d+00,
     $0.80694053195021761186d+00,0.73975604435269475868d+00,
     $0.66290966002478059546d+00,0.57719571005204581484d+00,
     $0.48361802694584102756d+00,0.38335932419873034692d+00,
     $0.27774982202182431507d+00,0.16823525155220746498d+00,
     $0.56344313046592789972d-01,0.99998243035489159858d+00,
     $0.99959879967191068325d+00,0.99831663531840739253d+00,
     $0.99572410469840718851d+00,0.99149572117810613240d+00,
     $0.98537149959852037111d+00,0.97714151463970571416d+00,
     $0.96663785155841656709d+00,0.95373000642576113641d+00/
      data (wa(i),i=41,60)/
     $0.93832039777959288365d+00,0.92034002547001242073d+00,
     $0.89974489977694003664d+00,0.87651341448470526974d+00,
     $0.85064449476835027976d+00,0.82215625436498040737d+00,
     $0.79108493379984836143d+00,0.75748396638051363793d+00,
     $0.72142308537009891548d+00,0.68298743109107922809d+00,
     $0.64227664250975951377d+00,0.59940393024224289297d+00,
     $0.55449513263193254887d+00,0.50768775753371660215d+00,
     $0.45913001198983233287d+00,0.40897982122988867241d+00,
     $0.35740383783153215238d+00,0.30457644155671404334d+00,
     $0.25067873030348317661d+00,0.19589750271110015392d+00/
      data (wa(i),i=61,80)/
     $0.14042423315256017459d+00,0.84454040083710883710d-01,
     $0.28184648949745694339d-01,0.99999759637974846462d+00,
     $0.99994399620705437576d+00,0.99976049092443204733d+00,
     $0.99938033802502358193d+00,0.99874561446809511470d+00,
     $0.99780535449595727456d+00,0.99651414591489027385d+00,
     $0.99483150280062100052d+00,0.99272134428278861533d+00,
     $0.99015137040077015918d+00,0.98709252795403406719d+00,
     $0.98351865757863272876d+00,0.97940628167086268381d+00,
     $0.97473445975240266776d+00,0.96948465950245923177d+00,
     $0.96364062156981213252d+00,0.95718821610986096274d+00/
      data (wa(i),i=81,100)/
     $0.95011529752129487656d+00,0.94241156519108305981d+00,
     $0.93406843615772578800d+00,0.92507893290707565236d+00,
     $0.91543758715576504064d+00,0.90514035881326159519d+00,
     $0.89418456833555902286d+00,0.88256884024734190684d+00,
     $0.87029305554811390585d+00,0.85735831088623215653d+00,
     $0.84376688267270860104d+00,0.82952219463740140018d+00,
     $0.81462878765513741344d+00,0.79909229096084140180d+00,
     $0.78291939411828301639d+00,0.76611781930376009072d+00,
     $0.74869629361693660282d+00,0.73066452124218126133d+00,
     $0.71203315536225203459d+00,0.69281376977911470289d+00/
      data (wa(i),i=101,120)/
     $0.67301883023041847920d+00,0.65266166541001749610d+00,
     $0.63175643771119423041d+00,0.61031811371518640016d+00,
     $0.58836243444766254143d+00,0.56590588542365442262d+00,
     $0.54296566649831149049d+00,0.51955966153745702199d+00,
     $0.49570640791876146017d+00,0.47142506587165887693d+00,
     $0.44673538766202847374d+00,0.42165768662616330006d+00,
     $0.39621280605761593918d+00,0.37042208795007823014d+00,
     $0.34430734159943802278d+00,0.31789081206847668318d+00,
     $0.29119514851824668196d+00,0.26424337241092676194d+00,
     $0.23705884558982972721d+00,0.20966523824318119477d+00/
      data (wa(i),i=121,127)/
     $0.18208649675925219825d+00,0.15434681148137810869d+00,
     $0.12647058437230196685d+00,0.98482396598119202090d-01,
     $0.70406976042855179063d-01,0.42269164765363603212d-01,
     $0.14093886410782462614d-01/
c  check for trivial case
      if(a.ge.b)then
        besr=0.
        besi=0.
        return
      endif
      if(new.eq.2)goto 200
c  scale factors
      sum=(b+a)/2.
      diff=(b-a)/2.
c  one point gauss
c//
c//      call f(sum,fzeror,fzeroi)
      zfun1=fun1(sum)
      if(new.ge.1) then
        call rel1(sum,zfun1,fwork,-jrel)
        zfun1=fwork(jrel)
      endif
      fzeror=dreal(zfun1)
      fzeroi=dimag(zfun1)
c//
      nf=nf+1
      besf=bessj(sum*r,order)
      fzeror=besf*fzeror
      fzeroi=besf*fzeroi
      if(new.eq.1)then
        karg(1,np)=sum
        kern(1,np)=fzeror
        kern(2,np)=fzeroi
c//
      kbes(1,np)=besf
c//
        la=2
        lk=3
      endif
      n=1
      kn=1
    5 continue
c  step through gauss orders
      do 40 k=kn,ng
c  compute new function values
        na=nwa(k)
        do 10 j=1,nwa(k)
          x=wa(na)*diff
          na=na+1
          sump=sum+x
          summ=sum-x
c//
c//          call f(sump,fr1(j),fi1(j))
        zfun1=fun1(sump)
      if(new.ge.1) then
        call rel1(sump,zfun1,fwork,-jrel)
        zfun1=fwork(jrel)
      endif
        fr1(j)=dreal(zfun1)
        fi1(j)=dimag(zfun1)
c//          call f(summ,fr2(j),fi2(j))
        zfun1=fun1(summ)
      if(new.ge.1) then
        call rel1(summ,zfun1,fwork,-jrel)
        zfun1=fwork(jrel)
      endif
        fr2(j)=dreal(zfun1)
        fi2(j)=dimag(zfun1)
c//
          nf=nf+2
          bes1(j)=bessj(sump*r,order)
          bes2(j)=bessj(summ*r,order)
          if(new.ge.1)then
            karg(la,np)=sump
            karg(la+1,np)=summ
c//
          kbes(la,np)=bes1(j)
          kbes(la+1,np)=bes2(j)
c//
            la=la+2
          endif
   10   continue
c  compute products of kernels and bessel functions
        do 20 j=1,nwa(k)
          fr1(j)=bes1(j)*fr1(j)
          fi1(j)=bes1(j)*fi1(j)
          fr2(j)=bes2(j)*fr2(j)
          fi2(j)=bes2(j)*fi2(j)
          funct(n)=fr1(j)+fr2(j)
          funct(n+1)=fi1(j)+fi2(j)
          n=n+2
   20   continue
        if(new.ge.1)then
          do 30 j=1,nwa(k)
            kern(lk,np)=fr1(j)
            kern(lk+1,np)=fi1(j)
            kern(lk+2,np)=fr2(j)
            kern(lk+3,np)=fi2(j)
            lk=lk+4
   30     continue
        endif
   40 continue
c  compute dot product of weights with integrand values
      nw=nwt(ng)
c  dot should be replaced with sdot from scilib on cray-1 if d.p. is
c  not used
      acumr=dot(nw,wt(nw),1,funct(1),2)
      acumi=dot(nw,wt(nw),1,funct(2),2)
      besr=(acumr+wt(2*nw)*fzeror)*diff
      besi=(acumi+wt(2*nw)*fzeroi)*diff
      if(np.le.nterm)nk(np)=ng
      return
  200 continue
c  construct funct from saved kernels
      diff=(b-a)/2.
      fzeror=kern(1,np)
      fzeroi=kern(2,np)
      k=min(nk(np),ng)
      nw=nwt(k)
      lk=3
      do 210 n=1,2*nw,2
        funct(n)=kern(lk,np)+kern(lk+2,np)
        funct(n+1)=kern(lk+1,np)+kern(lk+3,np)
        lk=lk+4
  210 continue
      if(nk(np).ge.ng)then
c  no additional orders required-compute dot product of weights with
c  integrand values and return
        acumr=dot(nw,wt(nw),1,funct(1),2)
        acumi=dot(nw,wt(nw),1,funct(2),2)
        besr=(acumr+wt(2*nw)*fzeror)*diff
        besi=(acumi+wt(2*nw)*fzeroi)*diff
        return
      else
c  compute additional orders before taking dot product
        sum=(b+a)/2.
        kn=k+1
        la=lk/2+1
        goto 5
      endif
      end
      subroutine padecf(sumr,sumi,sr,si,n)
c
c  computes sum(s(i)),i=1,...n by computation of pade approximant
c  using continued fraction expansion.  function is designed to be
c  called sequentially as n is incremented from 1 to its final value.
c  the nth continued fraction coefficient is calculated and
c  stored and the nth convergent returned.  it is up to the user to
c  stop the calculation when the desired accuracy is achieved.
c  algorithm from hanggi et al., z.naturforsch. 33a,402-417 (1977)
c  in their notation, vectors cfcor,cfcoi are lower case d, vectors dr,
c  di are upper case d, vectors xr,xi are x, and vectors sr,si are s
c
c  a.chave igpp/ucsd
c
      parameter (nc=100)
      implicit double precision (a-h,o-z)
      save xr,xi,dr,di,cfcor,cfcoi
      dimension sr(*),si(*),xr(nc),xi(nc),dr(0:nc),di(0:nc),cfcor(nc),
     $          cfcoi(nc)
      data dr(0)/-1./,di(0)/0./
      if(n.lt.3)then
c  initialize for recursive calculations
        if(n.eq.1)then
          do 10 i=1,nc
            xr(i)=0.
            xi(i)=0.
   10     continue
        endif
        dr(n)=sr(n)
        di(n)=si(n)
        denom=dr(n-1)**2+di(n-1)**2
        cfcor(n)=-(dr(n-1)*dr(n)+di(n-1)*di(n))/denom
        cfcoi(n)=-(dr(n-1)*di(n)-dr(n)*di(n-1))/denom
        call cf(sumr,sumi,cfcor,cfcoi,n)
        return
      else
        l=2*int((n-1)/2)
c  update x vectors for recursive calculation of cf coefficients
        do 20 k=l,4,-2
          xr(k)=xr(k-1)+cfcor(n-1)*xr(k-2)-cfcoi(n-1)*xi(k-2)
          xi(k)=xi(k-1)+cfcor(n-1)*xi(k-2)+cfcoi(n-1)*xr(k-2)
   20   continue
        xr(2)=xr(1)+cfcor(n-1)
        xi(2)=xi(1)+cfcoi(n-1)
c  interchange odd and even parts
        do 30 k=1,max(1,l-1),2
          t1=xr(k)
          t2=xi(k)
          xr(k)=xr(k+1)
          xi(k)=xi(k+1)
          xr(k+1)=t1
          xi(k+1)=t2
   30   continue
c  compute first coefficients
        dr(n)=sr(n)
        di(n)=si(n)
        do 40 k=1,max(1,l/2)
          dr(n)=dr(n)+sr(n-k)*xr(2*k-1)-si(n-k)*xi(2*k-1)
          di(n)=di(n)+si(n-k)*xr(2*k-1)+sr(n-k)*xi(2*k-1)
   40   continue
c  compute new cf coefficient
        denom=dr(n-1)**2+di(n-1)**2
        cfcor(n)=-(dr(n)*dr(n-1)+di(n)*di(n-1))/denom
        cfcoi(n)=-(dr(n-1)*di(n)-dr(n)*di(n-1))/denom
c  evaluate continued fraction
        call cf(sumr,sumi,cfcor,cfcoi,n)
        return
      endif
      end
      subroutine cf(resr,resi,cfcor,cfcoi,n)
c
c  evaluates a complex continued fraction by recursive division
c  starting at the bottom, as used by padecf
c  resr,resi are real and imaginary parts returned
c  cfcor,cfcoi are real and imaginary vectors of continued fraction
c  coefficients
c
      parameter (nc=100,ncm=nc-1)
      implicit double precision (a-h,o-z)
      dimension cfcor(*),cfcoi(*),one(nc)
      data one/0.,ncm*1./
      resr=one(n)+cfcor(n)
      resi=cfcoi(n)
      do 10 k=n-1,1,-1
        denom=resr**2+resi**2
        resro=resr
        resr=one(k)+(resr*cfcor(k)+resi*cfcoi(k))/denom
        resi=(resro*cfcoi(k)-resi*cfcor(k))/denom
   10 continue
      return
      end
      double precision function zeroj(nzero,order)
c
c  computes zero of bessel function of the first kind from
c  mcmahon's asymptotic expansion
c  nzero-number of the zero
c  order-order of the bessel function (0 or 1)
c
      implicit double precision (a-h,o-z)
      parameter (zt1=-1.d0/8.d0,zt2=124.d0/1536.d0,zt3=-120928.d0/
     $           491520.d0,zt4=401743168.d0/220200960.d0)
      parameter (ot1=3.d0/8.d0,ot2=-36.d0/1536.d0,ot3=113184.d0/
     $           491520.d0,ot4=-1951209.d0/220200960.d0)
      parameter (pi=3.141592653589793d0)
      if(order.eq.0.0)then
        beta=(nzero-.25)*pi
        zeroj=beta-zt1/beta-zt2/beta**3-zt3/beta**5-zt4/beta**7
      elseif(order.eq.1.0)then
        beta=(nzero+.25)*pi
        zeroj=beta-ot1/beta-ot2/beta**3-ot3/beta**5-ot4/beta**7
      endif
      return
      end
      double precision function dot(n,x1,inc1,x2,inc2)
c
c  computes dot product of two d.p. vectors with nonunit
c  incrementing allowed. replacement for blas subroutine sdot.
c
      implicit double precision (a-h,o-z)
      dimension x1(*),x2(*)
      if(inc2.gt.0)then
        k=1
      else
        k=n*abs(inc2)
      endif
      dot=0.0
      if(inc1.gt.0)then
        do 10 i=1,n,inc1
          dot=dot+x1(i)*x2(k)
          k=k+inc2
   10   continue
      else
        do 20 i=n,1,inc1
          dot=dot+x1(i)*x2(k)
          k=k+inc2
   20   continue
      endif
      return
      end
c-af
c-af  using freestanding routine bessj to call
c-af  toms bessel function jbess.
c-af  have changed name of this routine from jbess
c-af  to avoid conflict with toms library.
c-af
c-af      double precision function bessj(x,order)
c-afc
c-afc  computes bessel function of order "order" and argument x by calling
c-afc  imsl routines mmbsj0 and mmbsj1 (real*8 accuracy, in general).
c-afc
c-afc  if imsl library is not availble, then bessj could be changed to call
c-afc  nbs routines j0x and j1x (real*8 but approximately real*4 accuracy).
c-afc  << see c// below >>
c-afc
c-af      implicit double precision (a-h,o-z)
c-afc//   double precision j0x,j1x
c-af      double precision mmbsj0,mmbsj1
c-af      if(order.eq.0.0) then
c-af        bessj=mmbsj0(x,ier)
c-afc//       can replace with nbs routine:
c-afc//     bessj=j0x(x)
c-af      else if(order.eq.1.0) then
c-af        bessj=mmbsj1(x,ier)
c-afc//       can replace with nbs routine:
c-afc//     bessj=j1x(x)
c-af      endif
c-af      return
c-af      end
c
      subroutine dspln1(m,h,x,y,a,b,c,t,d,p,s)
c--one dimensional cubic spline interpolation
c
c  by w.l.anderson, u.s. geological survey, denver, colorado.
c
c  parms--- m= number of data points .gt. 2
c           h= equal interval option when h.gt.0.0d0 (use dummy x here),
c              unequal intervals if h=0.0d0 (x required storage)
c           x= indep.var when h=0.0d0 (dim .ge. m).
c           y= dependent variable  (dim .ge. m).
c           a,b,c=coeff.arrays (each dim .ge. m)
c                 results are returned in 1st(m-1) elements of a,b,&c.
c                 also used as work arrays during execution.
c           t= type of boundary condition supplied in d array. use
c              t=1 if 1st derivatives given at end points, or
c              t=0 if 2nd derivatives given at end points.
c              (t is integer*4)
c           d= boundary array (dim 2) at point 1 and m respectively.
c           p,s= work arrays (each dim=m).
c--error return with m=-(abs(m)) if any parm out of range.
c  the resulting cubic spline is of the form:
c     y=y(i)+a(i)*(x-x(i))+b(i)*(x-x(i))**2+c(i)*(x-x(i))**3
c       for i=1,2,...,m-1
c
      implicit real*8 (a-h,o-z)
      real*8  x(*),y(*),a(*),b(*),c(*),d(2),p(*),s(*),mul
      integer*4 t
      if(t.lt.0.or.t.gt.1.or.h.lt.0..or.m.lt.3) go to 999
      n=m-1
      if(t.eq.0) go to 20
c--1st derivative boundaries given
      ne=n-1
      if (h .lt. 0.0) then
        goto  999
      else if (h .gt. 0.0) then
        goto 1
      else
        goto 11
      endif

c--equal spacing h .gt. 0. and t=1
    1 hh=3.0d0/h
      do 2 i=1,ne
      b(i)=4.0d0
      c(i)=1.0d0
      a(i)=1.0d0
    2 p(i)=hh*(y(i+2)-y(i))
      p(1)=p(1)-d(1)
      p(ne)=p(ne)-d(2)
c--solution of tridiagonal matrix eq. of order ne
    3 c(1)=c(1)/b(1)
      p(1)=p(1)/b(1)
      do 4 i=2,ne
      mul=1.0d0/(b(i)-a(i)*c(i-1))
      c(i)=mul*c(i)
    4 p(i)=mul*(p(i)-a(i)*p(i-1))
c--obtain spline coefficients
      a(ne+t)=p(ne)
      i=ne-1
    5 a(i+t)=p(i)-c(i)*a(i+t+1)
      i=i-1
       if(i.ge.1) go to 5
      if(t.eq.0) go to 6
      a(1)=d(1)
      a(m)=d(2)
    6 if(h.eq.0.0d0) go to 14
      hh=1.0d0/h
      do 7 i=1,n
      mul=hh*(y(i+1)-y(i))
      b(i)=hh*(3.0d0*mul-(a(i+1)+2.0d0*a(i)))
    7 c(i)=hh*hh*(-2.0d0*mul+a(i+1)+a(i))
      return
c--unequal spacing h=0.. and t=1
   11 do 12 i=1,n
   12 s(i+1)=x(i+1)-x(i)
      do 13 i=1,ne
      b(i)=2.0d0*(s(i+1)+s(i+2))
      c(i)=s(i+1)
      a(i)=s(i+2)
   13 p(i)=3.0d0*(s(i+1)**2*(y(i+2)-y(i+1))+s(i+2)**2*(y(i+1)-y(i)))/
     *     (s(i+1)*s(i+2))
      p(1)=p(1)-s(3)*d(1)
      p(ne)=p(ne)-s(n)*d(2)
      go to 3
   14 do 15 i=1,n
      hh=1.0d0/s(i+1)
      mul=(y(i+1)-y(i))*hh**2
      b(i)=3.0d0*mul-(a(i+1)+2.0d0*a(i))*hh
   15 c(i)=-2.0d0*mul*hh+(a(i+1)+a(i))*hh**2
      return
c--2nd derivative boundaries given
   20 ne=n+1
      if (h .lt. 0.0) then
        goto  999
      else if (h .gt. 0.0) then
        goto 21
      else
        goto 31
      endif

c--equal spacing h .gt. 0 and t=0
   21 hh=3.0d0/h
      do 22 i=2,n
      b(i)=4.0d0
      c(i)=1.0d0
      a(i)=1.0d0
   22 p(i)=hh*(y(i+1)-y(i-1))
      b(1)=2.0d0
      b(ne)=2.0d0
      c(1)=1.0d0
      c(ne)=1.0d0
      a(ne)=1.0d0
      p(1)=hh*(y(2)-y(1))-0.5d0*h*d(1)
      p(ne)=hh*(y(m)-y(n))+0.5d0*h*d(2)
      go to 3
c--unequal spacing h=0 and t=0
   31 do 32 i=1,n
   32 s(i+1)=x(i+1)-x(i)
      n1=n-1
      do 33 i=1,n1
      b(i+1)=2.0d0*(s(i+1)+s(i+2))
      c(i+1)=s(i+1)
      a(i+1)=s(i+2)
   33 p(i+1)=3.0d0*(s(i+1)**2*(y(i+2)-y(i+1))+s(i+2)**2*(y(i+1)-y(i)))/
     *     (s(i+1)*s(i+2))
      b(1)=2.0d0
      b(ne)=2.0d0
      c(1)=1.0d0
      c(ne)=1.0d0
      a(ne)=1.0d0
      p(1)=3.0d0*(y(2)-y(1))/s(2)-0.5d0*s(2)*d(1)
      p(ne)=3.0d0*(y(m)-y(n))/s(m)+0.5d0*s(m)*d(2)
      go to 3
  999 m=-iabs(m)
      return
      end
      subroutine dsplnp(m,x,y,a,b,c,xx,yy)
c--given cubic spline coeff's a,b,c,and m obs.data arrays x,y
c  dsplnp evaluates the piecewise cubic spline ordinate yy at the
c  abscissa xx, where xx is in the closed interval (x(1),x(m)).
c  note: if computing over equal intervals, use the subr 'dcubic'
c  which requires only one call.
c
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),a(*),b(*),c(*)
      if(xx.lt.x(1).or.xx.gt.x(m)) go to 9
      m1=m-1
      do 1 i=1,m1
      j=i
      if(xx.le.x(i+1)) go to 2
1     continue
9     write(6,60) xx,x(1),x(m)
60    format('0error in dsplnp call--xx=',d16.8,'not in closed',
     $ ' interval (',d16.8,',',d16.8,')')
      return
2     z=xx-x(j)
      yy=y(j)+((c(j)*z+b(j))*z+a(j))*z
      return
      end
      subroutine hankl2(bmax,mb,nrel,tol,ntol,nord,fun1,rel1,
     *           zwork,fwork,maxdim,maxrel, zans,arg,nofun1,ierr)
c=======================================================================
c   a general algorithm for fast complex*16 hankel transforms of orders
c   0 and 1 using related and lagged convolutions.
c
c   by w.l.anderson, u.s.geological survey, denver, colorado, usa.
c=======================================================================
c
c   hankl2 is a double-precision (complex*16) version that uses the
c   double-precison (real*8) filters from subroutine dhankl as published
c   in the reference:
c
c   anderson, w.l., fast hankel transforms using related and lagged
c     convolutions: acm trans. on math. softw., v.8, n.4 (dec. 1982),
c     p.344-368.
c
c   hankl2 has been modified to exclude parameter ijrel & other changes
c   described in anderson (1982, p.365-366); also note the new variables
c   maxdim and maxrel as used in zans(maxdim,maxrel), maxdim>=mb and
c   maxrel>=nrel. (see definitions of new parms: rel1, maxdim, maxrel.)
c
c   hankl2 calls subprograms fun1 and rel1, where rel1 computes all
c   nrel>0 related input (kernel) functions (see anderson, 1982, p.365).
c   fun1 is the user's external input function of the hankel transform.
c   subroutine rel1 is only called whenever nrel>1, and has been defined
c   sightly different from subroutine relfun in anderson (1982, p.365).
c   (see definition of rel1 below.)
c=======================================================================
c
c      four general cases are possible using subprogram hankl2.
c
c    case 1.  single direct convolution at b=bmax=bmin (mb=1,nrel=1).
c    case 2.  related convolutions at a constant b=bmax (mb=1,nrel>1).
c    case 3.  lagged convolutions in (bmin,bmax) (mb>1,nrel=1).
c    case 4.  both related and lagged convolutions (mb>1,nrel>1).
c
c  machine dependent remarks
c
c      this subprogram was implemented and tested on a 64-bit double-
c    word machine with exp-range approximately 10**-38 to 10**+38 and
c    56-bit mantissa (about 16-decimal digits).  only double-precision
c    (real*8 and complex*16) and integer*4 operations are used.
c    for machines with other double-word sizes, changes in the number of
c    digits retained in some data statements may be required.
c
c  description of parameters
c
c    input (use of "real" implies "real*8" in comments below; also
c           use of "complex" implies "complex*16" in all comments):
c           -----------------------------------------------------
c      bmax - initial hankel transform argument b=bmax>0.0 (any case),
c             used in integral from 0 to infinity of
c             fun1(g)*jn(g*b)*dg, where jn=bessel function of order n,
c             n=0 or 1, and b>0.0. (see fun1 definition below).
c        mb - number of lagged convolutions desired (mb.ge.1).  use
c             mb=1 if b=bmin=bmax (i.e., case 1 or 2).  use
c             mb>1 if b is lagged in (bmin,bmax), where
c             bmin=bmax*exp(-.1*(mb-1)) does not underflow the exponent
c             range.  the b-lagged spacing is .1 in log-space. for
c             convenience in spline interpolation later, each b in
c             (bmin,bmax) is returned in array arg(i),i=1,mb, where
c             arg(i+1)/arg(i)=exp(.1) for all i.  if bmax>bmin>0 is
c             given, then an effective value of mb is determined as
c             mb=aint(10.*alog(bmax/bmin))+i, where i>1 is recommended,
c             particularly if using subsequent spline interpolation for
c             a different b-spacing than used in the sampled filters. if
c             spline interpolation is to be used later, it is generally
c             best to use dlog(arg(i)) instead of arg(i) -vs- zans(i,j),
c             for i=1,mb, and for any given j between 1 and nrel.  note
c             mb is used as an adjustable dimension in zans(mb,nrel).
c          ** (also see below new parameter maxdim>=mb now required.)
c      nrel - number of related convolutions desired (nrel.ge.1).  use
c             nrel=1 if only a single complex hankel transform is used.
c             nrel>1 requires user subroutine rel1 (see below).
c             note nrel is used as adjustable dimensions in arrays
c             zans(mb,nrel),zwork(801,nrel),nord(nrel).
c       tol - requested real truncation tolerance at both filter tails
c             for adaptive convolution for all mb*nrel transforms.  the
c             truncation criterion is established during convolution in
c             a fixed abscissa range (using weights 299-339) of either
c             order filter as the maximum absolute convolved product
c             times tol.  the convolution summation is terminated
c             on either side of the fixed range whenever the absolute
c             product .le. the truncation criterion.  both real and
c             imaginary parts of the complex summation must satisfy
c             the truncation criterion independently. in general, a
c             decreasing tolerance will produce higher accuracy since
c             more filter weights are used (unless exponent underflow
c             occurs in the transform input function evaluation).
c             one may set tol=0.0 to obtain maximum accuracy for all
c             mb*nrel complex hankel transforms in zans(mb,nrel).
c             however, the actual relative errors cannot be expected
c             to be smaller than about .1d-13 regardless of the tol
c             value used, since double-precision filter weights and
c             double-precision complex functions are used. in any event,
c             one should always choose tol<<desired relative error.
c          ** accuracy warning **  some highly oscillatory functions
c             fun1(g) and (or) limiting cases of b near machine-zero
c             (or infinity) should be avoided, otherwise unsatisfactory
c             results (e.g., relative & absolute errors>>tol) may occur.
c      ntol - number of consecutive times the truncation criterion (tol)
c             is to be met at either filter tail before filter
c             truncation occurs.  ntol=1 should be used for input
c             functions that do not have many zeros in (0,infinity). for
c             oscillatory functions with many zeros, ntol>1 may be used
c             to insure a premature cutoff does not occur for truncation
c             (see use of itol,ntol,tol in the code below).
c      nord - integer array nord(nrel) giving the nrel orders (0 or 1)
c             of each related hankel transform. if any nord(i),i=1,nrel,
c             is not 0 or 1, then order 1 will be assumed.
c     fun1 - name of an external declared complex*16 function written by
c            the user to evaluate fun1(g) for g>0.0.  an external
c            statement must appear in the user's calling program as:
c                  external fun1
c            and the complex*16 function must have the general form:
c            (replace "<...>" below as required.)
c
c                  complex*16 function fun1(g)
c                  implicit real*8 (a-h,o-z)
c                  common/pass/<...include parameters as needed...>
c                  <...>
c                  fun1=dcmplx(<...>,<...>)
c                  return
c                  end
c           notes:
c             fun1(g) must be a continuous decreasing complex function
c             for all real g>0.0.  the value of g must be unchanged
c             upon return from fun1. a multiple-pole of fun1(g) at g=0.0
c             can exist, provided the hankel transform converges (note
c             fun1(0.0) is not used). a single real function f1(g) may
c             be coded as fun1=dcmplx(f1(g),0.0). two independent real
c             functions f1(g),f2(g) may be integrated in parallel by
c             coding fun1=dcmplx(f1(g),f2(g)).  generally, fun1(g)
c             is defined analytically for all g>0.0.  however,
c             discretely defined functions may be used if fun1(g)
c             returns a smooth interpolation value (e.g., via cubic
c             splines) which satisfies the continuity condition for all
c             g>0, and provided the proper limiting value of fun1(g) is
c             given as g tends to infinity.  parameters other
c             than g needed in fun1(g) may be included by using labeled
c             common in fun1 and in the users calling program.  if
c             fun1(g) is an oscillating function, then the highest
c             frequency component (in log-space) should not exceed the
c             filter nyquist frequency, 1/(2*0.1).  in general,
c             subprogram hankl2 performs best when using smooth, well-
c             behaved functions fun1(g), that are characterized as
c             absolutely decreasing functions with relatively few
c             zeros for g>0. (see the accuracy warning under tol, and
c             error condition (4).)
c     rel1 - name of an external declared subroutine written by the user
c            to evaluate any (and all) nrel related kernels fj(nrel)
c            as a function of f1=fun1(g), assumed previously evaluated.
c            an external statement must appear in the user's program as:
c                  external rel1
c            and the subroutine must have the general form:
c            (replace "<...>" below as required.)
c
c                  subroutine rel1(g,f1,fj,j)
c                  implicit real*8 (a-h,o-z)
c                  complex*16 f1,fj(1)
c                  common/pass/<...include parameters as needed...>
c                  if(j.lt.0) go to (1,2,<...up to nrel...>), iabs(j)
c            1     fj(1)=f1
c                    if(j.lt.0) return
c            2     fj(2)=  <...expression in terms of g,f1, etc...>
c                    if(j.lt.0) return
c                  <...>
c                  return
c                  end
c           notes:
c             rel1 parameters g,f1,j are controlled by subroutine hankl2
c             where the user must only code the fj(j),j=1,nrel related
c             kernels as a function of g,f1, and any other parameters
c             that can be passed in common/pass/, as required.
c     zwork - complex work array zwork(801,nrel), which is used to hold
c             various computed functional values during related and
c             lagged convolutions.  a storage roll feature using
c             zwork(801,nrel) and internal array key(801) allows for
c             any b range (bmin,bmax) to be used during convolution.
c     fwork - complex work array fwork(maxrel), which is used in call
c             to rel1 when nrel>1.
c    maxdim - maximum dimension of zans(maxdim,nrel) in calling program,
c             where maxdim>=mb is required (else err=1 return,see below)
c    maxrel - maximum dimension of zans(maxdim,maxrel),zwork(801,maxrel)
c             and fwork(maxrel) in calling program, where maxrel>=nrel
c             is required.
c
c  output
c
c      zans - the complex*16 array zans(mb,nrel) is returned giving the
c             mb*nrel complex hankel transforms, with corresponding
c             b arguments given in real*8 array arg(mb).
c       arg - real*8 array arg(mb) is returned giving the resulting
c             b arguments in (bmin,bmax), where arg(i+1)/arg(i)=exp(.1),
c             i=1,mb-1 (this array could be eliminated to save storage
c             and regenerated after the call hankl2, if desired).
c    nofun1 - number of direct fun1 evaluations used for all mb*nrel
c             complex hankel transforms.  nofun1 is usually not more
c             than the number of weights needed for a single direct
c             convolution, for any mb and nrel.
c      ierr - error return code.  the following codes are possible --
c             = 0, no error in input parameters. (zans,arg computed.)
c             = 1, improper input parameters; i.e., mb<1,nrel<1,bmax<=0,
c               mb>maxdim, nrel>maxrel, or
c               bmax*exp(-.1*(mb-1))<=0.0.  (zans,arg not computed.)
c
c
c  error conditions
c
c    (1)  improper input parameters given (see ierr=1 above).
c    (2)  underflow conditions are possible during convolution, due to
c         the behavior of fun1, value of b in (bmin,bmax), tol, and
c         ntol. exponent and (or) arithmetic underflow traps must return
c         a value of 0.0 for the computer system being used.  note that
c         underflow may also occur in the users external function
c         fun1(g) for any value of g as set by subprogram hankl2.
c    (3)  an unrecoverable overflow condition can occur in executing
c         subr. rel1, depending on the value of b in (bmin,bmax),
c         tol, or the coding used in subroutine rel1.
c         in general, extremely small or large values of b should be
c         avoided (see accuracy warning under tol above).  also, in many
c         cases, exponent overflow can be avoided by proper choice of
c         fun1 and the related input function ordering defined in
c         subroutine rel1.
c    (4)  undetected errors are possible if fun1 or rel1 are
c         improperly coded, or if fun1(g)
c         does not yield doublee-precision complex accuracy, or is
c         not a continuous decreasing complex*16 function for all g>0.
c  usage
c
c    subprogram hankl2 is called as follows (use numerical values for
c    <expression>, excluding < and >, in declarations) --
c
c      implicit real*8 (a-h,o-z)
c      complex*16 zans(<maxdim>,<maxrel>),zwork(801,<maxrel>),
c     1        fwork(maxrel)
c      dimension arg(<mb>),nord(<nrel>)
c      external myfun1,myrel1
cc-----read/load input parameters for hankl2 as required
c      ...
c      call hankl2(bmax,mb,nrel,tol,ntol,nord,myfun1,myrel1,zwork,
c     * fwork,maxdim,maxrel, zans,arg,nofun1,ierr)
c      if(ierr.eq.1) stop
c      ...
c      end
c      complex*16 function myfun1(g)
c      implicit real*8 (a-h,o-z)
cc-----insert user supplied code for evaluation of myfun1(g),g>0.0
c      end
c      subroutine myrel1(g,f1,fj,j)
c      implicit real*8 (a-h,o-z)
c      complex*16 f1,fj(1)
cc-----insert user supplied code as described above for computing all
cc     nrel related kernels fj(j) as a function of g,f1=myfun1(g).
c      end
c
c=======================================================================
      implicit real*8 (a-h,o-z)
      integer mb,nrel,ntol,nofun1,ierr,maxdim,maxrel,nord(maxrel)
      real*8 arg(maxdim)
      complex*16 zwork(801,maxrel),zans(maxdim,maxrel),fwork(maxrel)
      complex*16 c,cmax,zsum,zero,fun1
      integer key(801)
      dimension t(2),tmax(2)
      dimension wt0(801),wt1(801)
      external fun1,rel1
c-----common/contrl/ preset by possible prior call to subroutine hy2fht.
c     to ignore use of jpath below, set jpath=0 (vax default).
      common/contrl/jpath,junk
c-----we define complex c,cmax to be equivalent to real*8 arrays t(2),
c     tmax(2), respectively, for use in the truncation criterion tests,
c     where c is any convolution product and cmax is the maximum
c     convolved product in the fixed abscissa range (see parameter tol).
      equivalence (c,t(1)),(cmax,tmax(1))
      data zero/(0.0d0,0.0d0)/
c-----abscis=base constant for filter abscissa generation
      data abscis/0.7059431685223780d0/
c-----e=dexp(.1d0), er=1.0d0/e (also used in abscissa generation)
      data e/1.10517091807564762d0/,er/.904837418035959573d0/
c-----wt0(i)=j0 hankel transform filter weights for i=1,801
      data
     *wt0(  1)/ 2.103562053838982d-29/,wt0(  2)/-1.264469361608894d-14/,
     *wt0(  3)/ 4.615731256788567d-14/,wt0(  4)/-2.798703374257668d-14/,
     *wt0(  5)/ 5.465764965410841d-14/,wt0(  6)/-2.652933109928729d-14/,
     *wt0(  7)/ 5.674913434067321d-14/,wt0(  8)/-2.157276828977208d-14/,
     *wt0(  9)/ 5.831846086773976d-14/,wt0( 10)/-1.546589284868783d-14/,
     *wt0( 11)/ 6.057302455652974d-14/,wt0( 12)/-8.502531259083065d-15/,
     *wt0( 13)/ 6.388018061147645d-14/,wt0( 14)/-5.659657635010288d-16/,
     *wt0( 15)/ 6.848500604791407d-14/,wt0( 16)/ 8.572897732168277d-15/,
     *wt0( 17)/ 7.465068154681814d-14/,wt0( 18)/ 1.920837293261338d-14/,
     *wt0( 19)/ 8.269345428975771d-14/,wt0( 20)/ 3.170116562922900d-14/,
     *wt0( 21)/ 9.300004039695208d-14/,wt0( 22)/ 4.649069639417992d-14/,
     *wt0( 23)/ 1.060441944490564d-13/,wt0( 24)/ 6.411216589597457d-14/,
     *wt0( 25)/ 1.224060834001701d-13/,wt0( 26)/ 8.521776751507022d-14/,
     *wt0( 27)/ 1.427957940487172d-13/,wt0( 28)/ 1.106026606968463d-13/,
     *wt0( 29)/ 1.680820203098405d-13/,wt0( 30)/ 1.412367028159546d-13/,
     *wt0( 31)/ 1.993271011776308d-13/,wt0( 32)/ 1.783032042964119d-13/,
     *wt0( 33)/ 2.378298196799422d-13/,wt0( 34)/ 2.232462602765097d-13/,
     *wt0( 35)/ 2.851776817107034d-13/,wt0( 36)/ 2.778285575785995d-13/,
     *wt0( 37)/ 3.433107735233557d-13/,wt0( 38)/ 3.442019764139949d-13/
      data
     *wt0( 39)/ 4.145997612327839d-13/,wt0( 40)/ 4.249938198224969d-13/,
     *wt0( 41)/ 5.019411631949951d-13/,wt0( 42)/ 5.234121310773422d-13/,
     *wt0( 43)/ 6.088737193247560d-13/,wt0( 44)/ 6.433743253810275d-13/,
     *wt0( 45)/ 7.397205280796933d-13/,wt0( 46)/ 7.896642978865920d-13/,
     *wt0( 47)/ 8.997626559623209d-13/,wt0( 48)/ 9.681243129571424d-13/,
     *wt0( 49)/ 1.095451187471558d-12/,wt0( 50)/ 1.185889375491455d-12/,
     *wt0( 51)/ 1.334666226164323d-12/,wt0( 52)/ 1.451673490118585d-12/,
     *wt0( 53)/ 1.627033241780647d-12/,wt0( 54)/ 1.776119296526221d-12/,
     *wt0( 55)/ 1.984309459865773d-12/,wt0( 56)/ 2.172225112712636d-12/,
     *wt0( 57)/ 2.420855801356264d-12/,wt0( 58)/ 2.655866524621122d-12/,
     *wt0( 59)/ 2.954213313000758d-12/,wt0( 60)/ 3.246433455110429d-12/,
     *wt0( 61)/ 3.605807223054322d-12/,wt0( 62)/ 3.967608279821300d-12/,
     *wt0( 63)/ 4.401806878720721d-12/,wt0( 64)/ 4.848316218220904d-12/,
     *wt0( 65)/ 5.374176077884770d-12/,wt0( 66)/ 5.923886142128417d-12/,
     *wt0( 67)/ 6.561955948855011d-12/,wt0( 68)/ 7.237468388830243d-12/,
     *wt0( 69)/ 8.012831864792648d-12/,wt0( 70)/ 8.841766480402105d-12/,
     *wt0( 71)/ 9.785047278800097d-12/,wt0( 72)/ 1.080115224902499d-11/,
     *wt0( 73)/ 1.194974128877562d-11/,wt0( 74)/ 1.319424925552153d-11/,
     *wt0( 75)/ 1.459380374689331d-11/,wt0( 76)/ 1.611708818260074d-11/
      data
     *wt0( 77)/ 1.782336249944076d-11/,wt0( 78)/ 1.968696083966218d-11/,
     *wt0( 79)/ 2.176804271234846d-11/,wt0( 80)/ 2.404712745375440d-11/,
     *wt0( 81)/ 2.658616922424555d-11/,wt0( 82)/ 2.937256616666064d-11/,
     *wt0( 83)/ 3.247112071587450d-11/,wt0( 84)/ 3.587699548548411d-11/,
     *wt0( 85)/ 3.965909071112399d-11/,wt0( 86)/ 4.382145152220635d-11/,
     *wt0( 87)/ 4.843856688602439d-11/,wt0( 88)/ 5.352476425684008d-11/,
     *wt0( 89)/ 5.916190912377549d-11/,wt0( 90)/ 6.537635327328955d-11/,
     *wt0( 91)/ 7.225949098391667d-11/,wt0( 92)/ 7.985185650562155d-11/,
     *wt0( 93)/ 8.825697213255412d-11/,wt0( 94)/ 9.753221923111176d-11/,
     *wt0( 95)/ 1.077963949370144d-10/,wt0( 96)/ 1.191270094182892d-10/,
     *wt0( 97)/ 1.316619519054340d-10/,wt0( 98)/ 1.455028951566729d-10/,
     *wt0( 99)/ 1.608114581091983d-10/,wt0(100)/ 1.777184270673645d-10/,
     *wt0(101)/ 1.964147916871289d-10/,wt0(102)/ 2.170665216346822d-10/,
     *wt0(103)/ 2.399008451839042d-10/,wt0(104)/ 2.651263504640323d-10/,
     *wt0(105)/ 2.930148720448538d-10/,wt0(106)/ 3.238267179640594d-10/,
     *wt0(107)/ 3.578885297833856d-10/,wt0(108)/ 3.955234710219298d-10/,
     *wt0(109)/ 4.371254308993604d-10/,wt0(110)/ 4.830940473937580d-10/,
     *wt0(111)/ 5.339056350072159d-10/,wt0(112)/ 5.900529573690070d-10/,
     *wt0(113)/ 6.521132758098962d-10/,wt0(114)/ 7.206928333934838d-10/
      data
     *wt0(115)/ 7.964924450372341d-10/,wt0(116)/ 8.802567084675067d-10/,
     *wt0(117)/ 9.728375895186283d-10/,wt0(118)/ 1.075148437456224d-09/,
     *wt0(119)/ 1.188226062693121d-09/,wt0(120)/ 1.313189706258061d-09/,
     *wt0(121)/ 1.451302163665563d-09/,wt0(122)/ 1.603933943511600d-09/,
     *wt0(123)/ 1.772624063293560d-09/,wt0(124)/ 1.959049733219870d-09/,
     *wt0(125)/ 2.165087540667258d-09/,wt0(126)/ 2.392789115986871d-09/,
     *wt0(127)/ 2.644443536014764d-09/,wt0(128)/ 2.922559573439235d-09/,
     *wt0(129)/ 3.229930291248597d-09/,wt0(130)/ 3.569622651576204d-09/,
     *wt0(131)/ 3.945045448172924d-09/,wt0(132)/ 4.359947261255963d-09/,
     *wt0(133)/ 4.818489091363688d-09/,wt0(134)/ 5.325251901762953d-09/,
     *wt0(135)/ 5.885315583343685d-09/,wt0(136)/ 6.504277635547376d-09/,
     *wt0(137)/ 7.188340419242570d-09/,wt0(138)/ 7.944342903082902d-09/,
     *wt0(139)/ 8.779858562959109d-09/,wt0(140)/ 9.703242578022093d-09/,
     *wt0(141)/ 1.072374322768940d-08/,wt0(142)/ 1.185156747840095d-08/,
     *wt0(143)/ 1.309800933225301d-08/,wt0(144)/ 1.447553742402158d-08/,
     *wt0(145)/ 1.599794451372024d-08/,wt0(146)/ 1.768046154055318d-08/,
     *wt0(147)/ 1.953993335487118d-08/,wt0(148)/ 2.159496468450500d-08/,
     *wt0(149)/ 2.386612830616206d-08/,wt0(150)/ 2.637614961034523d-08/,
     *wt0(151)/ 2.915015476269831d-08/,wt0(152)/ 3.221590205565798d-08/
      data
     *wt0(153)/ 3.560407926098509d-08/,wt0(154)/ 3.934859178954436d-08/,
     *wt0(155)/ 4.348692045365787d-08/,wt0(156)/ 4.806047869438007d-08/,
     *wt0(157)/ 5.311504443749314d-08/,wt0(158)/ 5.870120138001748d-08/,
     *wt0(159)/ 6.487486163571268d-08/,wt0(160)/ 7.169780940885934d-08/,
     *wt0(161)/ 7.923833480505012d-08/,wt0(162)/ 8.757190229426675d-08/,
     *wt0(163)/ 9.678192055835553d-08/,wt0(164)/ 1.069605631204864d-07/,
     *wt0(165)/ 1.182097045925483d-07/,wt0(166)/ 1.306419269237675d-07/,
     *wt0(167)/ 1.443816591198488d-07/,wt0(168)/ 1.595664099835792d-07/,
     *wt0(169)/ 1.763481565722233d-07/,wt0(170)/ 1.948948533650376d-07/,
     *wt0(171)/ 2.153921247351847d-07/,wt0(172)/ 2.380451115468297d-07/,
     *wt0(173)/ 2.630805351444854d-07/,wt0(174)/ 2.907489558998551d-07/,
     *wt0(175)/ 3.213272911558462d-07/,wt0(176)/ 3.551215767529809d-07/,
     *wt0(177)/ 3.924700396067688d-07/,wt0(178)/ 4.337464734078393d-07/,
     *wt0(179)/ 4.793639887921159d-07/,wt0(180)/ 5.297791390370199d-07/,
     *wt0(181)/ 5.854964979982190d-07/,wt0(182)/ 6.470737017046583d-07/,
     *wt0(183)/ 7.151270374758372d-07/,wt0(184)/ 7.903376040582035d-07/,
     *wt0(185)/ 8.734581359370587d-07/,wt0(186)/ 9.653205295304359d-07/,
     *wt0(187)/ 1.066844176299330d-06/,wt0(188)/ 1.179045157323447d-06/,
     *wt0(189)/ 1.303046419382764d-06/,wt0(190)/ 1.440089007182073d-06/
      data
     *wt0(191)/ 1.591544490556806d-06/,wt0(192)/ 1.758928685376697d-06/,
     *wt0(193)/ 1.943916830388476d-06/,wt0(194)/ 2.148360347694608d-06/,
     *wt0(195)/ 2.374305378111399d-06/,wt0(196)/ 2.624013254094558d-06/,
     *wt0(197)/ 2.899983137292810d-06/,wt0(198)/ 3.204977025773239d-06/,
     *wt0(199)/ 3.542047402097807d-06/,wt0(200)/ 3.914567778667102d-06/,
     *wt0(201)/ 4.326266465748067d-06/,wt0(202)/ 4.781263881007758d-06/,
     *wt0(203)/ 5.284113792545945d-06/,wt0(204)/ 5.839848890150422d-06/,
     *wt0(205)/ 6.454031158395685d-06/,wt0(206)/ 7.132807538712648d-06/,
     *wt0(207)/ 7.882971454044484d-06/,wt0(208)/ 8.712030795792633d-06/,
     *wt0(209)/ 9.628283069078871d-06/,wt0(210)/ 1.064089843325841d-05/,
     *wt0(211)/ 1.176001148348469d-05/,wt0(212)/ 1.299682267761921d-05/,
     *wt0(213)/ 1.436371043746991d-05/,wt0(214)/ 1.587435503282041d-05/,
     *wt0(215)/ 1.754387550120777d-05/,wt0(216)/ 1.938898096105121d-05/,
     *wt0(217)/ 2.142813784487508d-05/,wt0(218)/ 2.368175471430264d-05/,
     *wt0(219)/ 2.617238651818130d-05/,wt0(220)/ 2.892496032668667d-05/,
     *wt0(221)/ 3.196702481167932d-05/,wt0(222)/ 3.532902595927333d-05/,
     *wt0(223)/ 3.904461178455365d-05/,wt0(224)/ 4.315096908756330d-05/,
     *wt0(225)/ 4.768919563102351d-05/,wt0(226)/ 5.270471145419851d-05/,
     *wt0(227)/ 5.824771344935570d-05/,wt0(228)/ 6.437367773954276d-05/
      data
     *wt0(229)/ 7.114391489568718d-05/,wt0(230)/ 7.862618353777022d-05/,
     *wt0(231)/ 8.689536847209423d-05/,wt0(232)/ 9.603423013681158d-05/,
     *wt0(233)/ 1.061342328695059d-04/,wt0(234)/ 1.172964602657162d-04/,
     *wt0(235)/ 1.296326268075215d-04/,wt0(236)/ 1.432661958546299d-04/,
     *wt0(237)/ 1.583336152150400d-04/,wt0(238)/ 1.749856826064412d-04/,
     *wt0(239)/ 1.933890547221590d-04/,wt0(240)/ 2.137279149064782d-04/,
     *wt0(241)/ 2.362058162174614d-04/,wt0(242)/ 2.610477181410775d-04/,
     *wt0(243)/ 2.885022375064333d-04/,wt0(244)/ 3.188441357870562d-04/,
     *wt0(245)/ 3.523770680013592d-04/,wt0(246)/ 3.894366200728159d-04/,
     *wt0(247)/ 4.303936656695837d-04/,wt0(248)/ 4.756580748770003d-04/,
     *wt0(249)/ 5.256828130316798d-04/,wt0(250)/ 5.809684683619367d-04/,
     *wt0(251)/ 6.420682561143772d-04/,wt0(252)/ 7.095935446942389d-04/,
     *wt0(253)/ 7.842199637491256d-04/,wt0(254)/ 8.666941465924691d-04/,
     *wt0(255)/ 9.578411834759010d-04/,wt0(256)/ 1.058572843543458d-03/,
     *wt0(257)/ 1.169896665384001d-03/,wt0(258)/ 1.292925974983324d-03/,
     *wt0(259)/ 1.428890965739396d-03/,wt0(260)/ 1.579150889627630d-03/,
     *wt0(261)/ 1.745207548623177d-03/,wt0(262)/ 1.928720102656342d-03/,
     *wt0(263)/ 2.131521473093367d-03/,wt0(264)/ 2.355636277673631d-03/,
     *wt0(265)/ 2.603300731229735d-03/,wt0(266)/ 2.876984272542283d-03/
      data
     *wt0(267)/ 3.179413630050106d-03/,wt0(268)/ 3.513598723334281d-03/,
     *wt0(269)/ 3.882861625673649d-03/,wt0(270)/ 4.290867254031690d-03/,
     *wt0(271)/ 4.741657973486860d-03/,wt0(272)/ 5.239689338467119d-03/,
     *wt0(273)/ 5.789870985188969d-03/,wt0(274)/ 6.397607072000340d-03/,
     *wt0(275)/ 7.068843783102744d-03/,wt0(276)/ 7.810112796625723d-03/,
     *wt0(277)/ 8.628584975720666d-03/,wt0(278)/ 9.532112534855741d-03/,
     *wt0(279)/ 1.052928697105230d-02/,wt0(280)/ 1.162947042684567d-02/,
     *wt0(281)/ 1.284285301027502d-02/,wt0(282)/ 1.418045400408490d-02/,
     *wt0(283)/ 1.565416842684913d-02/,wt0(284)/ 1.727670025860997d-02/,
     *wt0(285)/ 1.906157878467033d-02/,wt0(286)/ 2.102295173845023d-02/,
     *wt0(287)/ 2.317553622087861d-02/,wt0(288)/ 2.553413681316817d-02/,
     *wt0(289)/ 2.811347056705691d-02/,wt0(290)/ 3.092716134472897d-02/,
     *wt0(291)/ 3.398734108220733d-02/,wt0(292)/ 3.730266904318316d-02/,
     *wt0(293)/ 4.087756584601726d-02/,wt0(294)/ 4.470845458227711d-02/,
     *wt0(295)/ 4.878245655564208d-02/,wt0(296)/ 5.307046388837388d-02/,
     *wt0(297)/ 5.752521479656082d-02/,wt0(298)/ 6.206888924914453d-02/,
     *wt0(299)/ 6.659098690582325d-02/,wt0(300)/ 7.092687097703928d-02/,
     *wt0(301)/ 7.485762155329997d-02/,wt0(302)/ 7.807466421172765d-02/,
     *wt0(303)/ 8.018887211338042d-02/,wt0(304)/ 8.067640670918657d-02/
      data
     *wt0(305)/ 7.891767306422777d-02/,wt0(306)/ 7.412406301630496d-02/,
     *wt0(307)/ 6.545864753141331d-02/,wt0(308)/ 5.195771725733346d-02/,
     *wt0(309)/ 3.284797274184859d-02/,wt0(310)/ 7.497076325831270d-03/,
     *wt0(311)/-2.386612869894549d-02/,wt0(312)/-6.017494378476118d-02/,
     *wt0(313)/-9.817899798850670d-02/,wt0(314)/-1.328147797272611d-01/,
     *wt0(315)/-1.554628569772562d-01/,wt0(316)/-1.563982157987450d-01/,
     *wt0(317)/-1.243049866529050d-01/,wt0(318)/-5.486815986343697d-02/,
     *wt0(319)/ 4.686255899170307d-02/,wt0(320)/ 1.511218295806206d-01/,
     *wt0(321)/ 2.119315534410599d-01/,wt0(322)/ 1.695134135887796d-01/,
     *wt0(323)/ 1.386197420331480d-02/,wt0(324)/-1.869350451838135d-01/,
     *wt0(325)/-2.455889606925336d-01/,wt0(326)/-5.309269401899814d-02/,
     *wt0(327)/ 2.519998415798595d-01/,wt0(328)/ 1.968224876057428d-01/,
     *wt0(329)/-2.014333618969160d-01/,wt0(330)/-2.458450827258603d-01/,
     *wt0(331)/ 3.433559314076636d-01/,wt0(332)/-4.770066510626292d-02/,
     *wt0(333)/-2.096628473507452d-01/,wt0(334)/ 2.509085123794248d-01/,
     *wt0(335)/-1.676412661335192d-01/,wt0(336)/ 7.495199786896762d-02/,
     *wt0(337)/-1.695135102751146d-02/,wt0(338)/-8.828524201374996d-03/,
     *wt0(339)/ 1.616787356183550d-02/,wt0(340)/-1.556411755992529d-02/,
     *wt0(341)/ 1.251311125216370d-02/,wt0(342)/-9.299315246978098d-03/
      data
     *wt0(343)/ 6.650525996779493d-03/,wt0(344)/-4.667707408671884d-03/,
     *wt0(345)/ 3.248885545556855d-03/,wt0(346)/-2.255529108544675d-03/,
     *wt0(347)/ 1.566875211100915d-03/,wt0(348)/-1.091132486368091d-03/,
     *wt0(349)/ 7.625383049040949d-04/,wt0(350)/-5.352540577820408d-04/,
     *wt0(351)/ 3.777104567214670d-04/,wt0(352)/-2.682528211535297d-04/,
     *wt0(353)/ 1.920265274996755d-04/,wt0(354)/-1.388213723208944d-04/,
     *wt0(355)/ 1.015998917810005d-04/,wt0(356)/-7.549725062631259d-05/,
     *wt0(357)/ 5.714172920052936d-05/,wt0(358)/-4.419160967680634d-05/,
     *wt0(359)/ 3.501776199578583d-05/,wt0(360)/-2.848496954139990d-05/,
     *wt0(361)/ 2.380104817324115d-05/,wt0(362)/-2.041274924032166d-05/,
     *wt0(363)/ 1.793363489791860d-05/,wt0(364)/-1.609367780777977d-05/,
     *wt0(365)/ 1.470394692017765d-05/,wt0(366)/-1.363208540087839d-05/,
     *wt0(367)/ 1.278540851480116d-05/,wt0(368)/-1.209910323276592d-05/,
     *wt0(369)/ 1.152780746372085d-05/,wt0(370)/-1.103964249967709d-05/,
     *wt0(371)/ 1.061214724478712d-05/,wt0(372)/-1.022954672255804d-05/,
     *wt0(373)/ 9.880801917269063d-06/,wt0(374)/-9.558136268373164d-06/,
     *wt0(375)/ 9.255989462874843d-06/,wt0(376)/-8.970368039834961d-06/,
     *wt0(377)/ 8.698439226217636d-06/,wt0(378)/-8.438199899378174d-06/,
     *wt0(379)/ 8.188185039708096d-06/,wt0(380)/-7.947275976752668d-06/
      data
     *wt0(381)/ 7.714620489277084d-06/,wt0(382)/-7.489590334670121d-06/,
     *wt0(383)/ 7.271712845563131d-06/,wt0(384)/-7.060595202720040d-06/,
     *wt0(385)/ 6.855890958634290d-06/,wt0(386)/-6.657307826365352d-06/,
     *wt0(387)/ 6.464610635947412d-06/,wt0(388)/-6.277596705248834d-06/,
     *wt0(389)/ 6.096069004564178d-06/,wt0(390)/-5.919835014912806d-06/,
     *wt0(391)/ 5.748721389173238d-06/,wt0(392)/-5.582576428368341d-06/,
     *wt0(393)/ 5.421255729073172d-06/,wt0(394)/-5.264611483233620d-06/,
     *wt0(395)/ 5.112497827757443d-06/,wt0(396)/-4.964780650195476d-06/,
     *wt0(397)/ 4.821336526159266d-06/,wt0(398)/-4.682043530422997d-06/,
     *wt0(399)/ 4.546777587820041d-06/,wt0(400)/-4.415418000076004d-06/,
     *wt0(401)/ 4.287852503112911d-06/,wt0(402)/-4.163974668109847d-06/,
     *wt0(403)/ 4.043678407829827d-06/,wt0(404)/-3.926857526517075d-06/,
     *wt0(405)/ 3.813409851360361d-06/,wt0(406)/-3.703239239834851d-06/,
     *wt0(407)/ 3.596252927182967d-06/,wt0(408)/-3.492358600906126d-06/,
     *wt0(409)/ 3.391465179914198d-06/,wt0(410)/-3.293485431241940d-06/,
     *wt0(411)/ 3.198336356446399d-06/,wt0(412)/-3.105937110348406d-06/,
     *wt0(413)/ 3.016207679069661d-06/,wt0(414)/-2.929069911562291d-06/,
     *wt0(415)/ 2.844448954094881d-06/,wt0(416)/-2.762272935820411d-06/,
     *wt0(417)/ 2.682471568415877d-06/,wt0(418)/-2.604975709140556d-06/
      data
     *wt0(419)/ 2.529718233664062d-06/,wt0(420)/-2.456634690985054d-06/,
     *wt0(421)/ 2.385662793855871d-06/,wt0(422)/-2.316741584624232d-06/,
     *wt0(423)/ 2.249811415880699d-06/,wt0(424)/-2.184814552602377d-06/,
     *wt0(425)/ 2.121695378987067d-06/,wt0(426)/-2.060399926949900d-06/,
     *wt0(427)/ 2.000875435893980d-06/,wt0(428)/-1.943070483347481d-06/,
     *wt0(429)/ 1.886935339579619d-06/,wt0(430)/-1.832421953831803d-06/,
     *wt0(431)/ 1.779483599068986d-06/,wt0(432)/-1.728074671861562d-06/,
     *wt0(433)/ 1.678150840818710d-06/,wt0(434)/-1.629669220867866d-06/,
     *wt0(435)/ 1.582588274939725d-06/,wt0(436)/-1.536867577717458d-06/,
     *wt0(437)/ 1.492467744182388d-06/,wt0(438)/-1.449350544031671d-06/,
     *wt0(439)/ 1.407478963215944d-06/,wt0(440)/-1.366817091697906d-06/,
     *wt0(441)/ 1.327329980403758d-06/,wt0(442)/-1.288983628354532d-06/,
     *wt0(443)/ 1.251745052743577d-06/,wt0(444)/-1.215582289623042d-06/,
     *wt0(445)/ 1.180464297731392d-06/,wt0(446)/-1.146360878232118d-06/,
     *wt0(447)/ 1.113242681755176d-06/,wt0(448)/-1.081081241360643d-06/,
     *wt0(449)/ 1.049848946307979d-06/,wt0(450)/-1.019518969919736d-06/,
     *wt0(451)/ 9.900652269648073d-07/,wt0(452)/-9.614623823829486d-07/,
     *wt0(453)/ 9.336858590396268d-07/,wt0(454)/-9.067118039187765d-07/,
     *wt0(455)/ 8.805170379415897d-07/,wt0(456)/-8.550790335506564d-07/
      data
     *wt0(457)/ 8.303759185665738d-07/,wt0(458)/-8.063864695073603d-07/,
     *wt0(459)/ 7.830900795968434d-07/,wt0(460)/-7.604667251660013d-07/,
     *wt0(461)/ 7.384969528074763d-07/,wt0(462)/-7.171618775580271d-07/,
     *wt0(463)/ 6.964431695613518d-07/,wt0(464)/-6.763230273695452d-07/,
     *wt0(465)/ 6.567841555251206d-07/,wt0(466)/-6.378097555683096d-07/,
     *wt0(467)/ 6.193835199757372d-07/,wt0(468)/-6.014896169720738d-07/,
     *wt0(469)/ 5.841126695566173d-07/,wt0(470)/-5.672377400509635d-07/,
     *wt0(471)/ 5.508503223216395d-07/,wt0(472)/-5.349363333399711d-07/,
     *wt0(473)/ 5.194820986824491d-07/,wt0(474)/-5.044743364413132d-07/,
     *wt0(475)/ 4.899001458890472d-07/,wt0(476)/-4.757469999761242d-07/,
     *wt0(477)/ 4.620027360551677d-07/,wt0(478)/-4.486555431391317d-07/,
     *wt0(479)/ 4.356939495368206d-07/,wt0(480)/-4.231068139067277d-07/,
     *wt0(481)/ 4.108833178901058d-07/,wt0(482)/-3.990129570237777d-07/,
     *wt0(483)/ 3.874855299952715d-07/,wt0(484)/-3.762911289578004d-07/,
     *wt0(485)/ 3.654201320223786d-07/,wt0(486)/-3.548631961765860d-07/,
     *wt0(487)/ 3.446112489455190d-07/,wt0(488)/-3.346554794873414d-07/,
     *wt0(489)/ 3.249873307934947d-07/,wt0(490)/-3.155984931421240d-07/,
     *wt0(491)/ 3.064808974928670d-07/,wt0(492)/-2.976267081233949d-07/,
     *wt0(493)/ 2.890283152696976d-07/,wt0(494)/-2.806783286653362d-07/
      data
     *wt0(495)/ 2.725695717362821d-07/,wt0(496)/-2.646950756050221d-07/,
     *wt0(497)/ 2.570480727284082d-07/,wt0(498)/-2.496219907786807d-07/,
     *wt0(499)/ 2.424104471687142d-07/,wt0(500)/-2.354072438896643d-07/,
     *wt0(501)/ 2.286063621819694d-07/,wt0(502)/-2.220019570969968d-07/,
     *wt0(503)/ 2.155883523600859d-07/,wt0(504)/-2.093600356601861d-07/,
     *wt0(505)/ 2.033116540781646d-07/,wt0(506)/-1.974380094205501d-07/,
     *wt0(507)/ 1.917340535876188d-07/,wt0(508)/-1.861948842172547d-07/,
     *wt0(509)/ 1.808157405985162d-07/,wt0(510)/-1.755919996494557d-07/,
     *wt0(511)/ 1.705191718700136d-07/,wt0(512)/-1.655928973946416d-07/,
     *wt0(513)/ 1.608089422670497d-07/,wt0(514)/-1.561631948835526d-07/,
     *wt0(515)/ 1.516516624771349d-07/,wt0(516)/-1.472704676262496d-07/,
     *wt0(517)/ 1.430158448819191d-07/,wt0(518)/-1.388841375627619d-07/,
     *wt0(519)/ 1.348717946589134d-07/,wt0(520)/-1.309753677747614d-07/,
     *wt0(521)/ 1.271915081246862d-07/,wt0(522)/-1.235169636417206d-07/,
     *wt0(523)/ 1.199485762099940d-07/,wt0(524)/-1.164832789732080d-07/,
     *wt0(525)/ 1.131180936862105d-07/,wt0(526)/-1.098501281312535d-07/,
     *wt0(527)/ 1.066765736321394d-07/,wt0(528)/-1.035947026599537d-07/,
     *wt0(529)/ 1.006018664975167d-07/,wt0(530)/-9.769549295058651d-08/,
     *wt0(531)/ 9.487308412474974d-08/,wt0(532)/-9.213221428331440d-08/
      data
     *wt0(533)/ 8.947052777472484d-08/,wt0(534)/-8.688573700943115d-08/,
     *wt0(535)/ 8.437562048443172d-08/,wt0(536)/-8.193802086874367d-08/,
     *wt0(537)/ 7.957084315467929d-08/,wt0(538)/-7.727205286392824d-08/,
     *wt0(539)/ 7.503967429747869d-08/,wt0(540)/-7.287178883148888d-08/,
     *wt0(541)/ 7.076653326692979d-08/,wt0(542)/-6.872209823264986d-08/,
     *wt0(543)/ 6.673672663335082d-08/,wt0(544)/-6.480871213716780d-08/,
     *wt0(545)/ 6.293639770567428d-08/,wt0(546)/-6.111817417008917d-08/,
     *wt0(547)/ 5.935247885127626d-08/,wt0(548)/-5.763779421770836d-08/,
     *wt0(549)/ 5.597264657919956d-08/,wt0(550)/-5.435560481863042d-08/,
     *wt0(551)/ 5.278527916288309d-08/,wt0(552)/-5.126031999018991d-08/,
     *wt0(553)/ 4.977941667022704d-08/,wt0(554)/-4.834129643623350d-08/,
     *wt0(555)/ 4.694472329047102d-08/,wt0(556)/-4.558849694284832d-08/,
     *wt0(557)/ 4.427145178025873d-08/,wt0(558)/-4.299245586445184d-08/,
     *wt0(559)/ 4.175040995827635d-08/,wt0(560)/-4.054424658081963d-08/,
     *wt0(561)/ 3.937292909065092d-08/,wt0(562)/-3.823545079527342d-08/,
     *wt0(563)/ 3.713083408552377d-08/,wt0(564)/-3.605812959486601d-08/,
     *wt0(565)/ 3.501641538356351d-08/,wt0(566)/-3.400479614676831d-08/,
     *wt0(567)/ 3.302240244515663d-08/,wt0(568)/-3.206838995736202d-08/,
     *wt0(569)/ 3.114193875409117d-08/,wt0(570)/-3.024225259359901d-08/
      data
     *wt0(571)/ 2.936855823761587d-08/,wt0(572)/-2.852010478676213d-08/,
     *wt0(573)/ 2.769616303496246d-08/,wt0(574)/-2.689602484264953d-08/,
     *wt0(575)/ 2.611900252829978d-08/,wt0(576)/-2.536442827753140d-08/,
     *wt0(577)/ 2.463165356908288d-08/,wt0(578)/-2.392004861730570d-08/,
     *wt0(579)/ 2.322900183088905d-08/,wt0(580)/-2.255791928733344d-08/,
     *wt0(581)/ 2.190622422255186d-08/,wt0(582)/-2.127335653510145d-08/,
     *wt0(583)/ 2.065877230473273d-08/,wt0(584)/-2.006194332493974d-08/,
     *wt0(585)/ 1.948235664905881d-08/,wt0(586)/-1.891951414942437d-08/,
     *wt0(587)/ 1.837293208920017d-08/,wt0(588)/-1.784214070660012d-08/,
     *wt0(589)/ 1.732668381117887d-08/,wt0(590)/-1.682611839179446d-08/,
     *wt0(591)/ 1.634001423585344d-08/,wt0(592)/-1.586795355952947d-08/,
     *wt0(593)/ 1.540953064868929d-08/,wt0(594)/-1.496435151022419d-08/,
     *wt0(595)/ 1.453203353344868d-08/,wt0(596)/-1.411220516125227d-08/,
     *wt0(597)/ 1.370450557074419d-08/,wt0(598)/-1.330858436314438d-08/,
     *wt0(599)/ 1.292410126264865d-08/,wt0(600)/-1.255072582398340d-08/,
     *wt0(601)/ 1.218813714839246d-08/,wt0(602)/-1.183602360782992d-08/,
     *wt0(603)/ 1.149408257713515d-08/,wt0(604)/-1.116202017395058d-08/,
     *wt0(605)/ 1.083955100614394d-08/,wt0(606)/-1.052639792651900d-08/,
     *wt0(607)/ 1.022229179461678d-08/,wt0(608)/-9.926971245406734d-09/
      data
     *wt0(609)/ 9.640182464660921d-09/,wt0(610)/-9.361678970810806d-09/,
     *wt0(611)/ 9.091221403102823d-09/,wt0(612)/-8.828577315878212d-09/,
     *wt0(613)/ 8.573520978800615d-09/,wt0(614)/-8.325833182853652d-09/,
     *wt0(615)/ 8.085301051938885d-09/,wt0(616)/-7.851717859915976d-09/,
     *wt0(617)/ 7.624882852931747d-09/,wt0(618)/-7.404601076884107d-09/,
     *wt0(619)/ 7.190683209868762d-09/,wt0(620)/-6.982945399464102d-09/,
     *wt0(621)/ 6.781209104717332d-09/,wt0(622)/-6.585300942697731d-09/,
     *wt0(623)/ 6.395052539483568d-09/,wt0(624)/-6.210300385452473d-09/,
     *wt0(625)/ 6.030885694751333d-09/,wt0(626)/-5.856654268826800d-09/,
     *wt0(627)/ 5.687456363899640d-09/,wt0(628)/-5.523146562267662d-09/,
     *wt0(629)/ 5.363583647325650d-09/,wt0(630)/-5.208630482195534d-09/,
     *wt0(631)/ 5.058153891863612d-09/,wt0(632)/-4.912024548723426d-09/,
     *wt0(633)/ 4.770116861425010d-09/,wt0(634)/-4.632308866934630d-09/,
     *wt0(635)/ 4.498482125712853d-09/,wt0(636)/-4.368521619921419d-09/,
     *wt0(637)/ 4.242315654571151d-09/,wt0(638)/-4.119755761525383d-09/,
     *wt0(639)/ 4.000736606276379d-09/,wt0(640)/-3.885155897415073d-09/,
     *wt0(641)/ 3.772914298716528d-09/,wt0(642)/-3.663915343765296d-09/,
     *wt0(643)/ 3.558065353046999d-09/,wt0(644)/-3.455273353434971d-09/,
     *wt0(645)/ 3.355451000003056d-09/,wt0(646)/-3.258512500097394d-09/
      data
     *wt0(647)/ 3.164374539601736d-09/,wt0(648)/-3.072956211332787d-09/,
     *wt0(649)/ 2.984178945504153d-09/,wt0(650)/-2.897966442199276d-09/,
     *wt0(651)/ 2.814244605795313d-09/,wt0(652)/-2.732941481281439d-09/,
     *wt0(653)/ 2.653987192416816d-09/,wt0(654)/-2.577313881675167d-09/,
     *wt0(655)/ 2.502855651924429d-09/,wt0(656)/-2.430548509791296d-09/,
     *wt0(657)/ 2.360330310661908d-09/,wt0(658)/-2.292140705271406d-09/,
     *wt0(659)/ 2.225921087836540d-09/,wt0(660)/-2.161614545686758d-09/,
     *wt0(661)/ 2.099165810350450d-09/,wt0(662)/-2.038521210054267d-09/,
     *wt0(663)/ 1.979628623594711d-09/,wt0(664)/-1.922437435542387d-09/,
     *wt0(665)/ 1.866898492740427d-09/,wt0(666)/-1.812964062059653d-09/,
     *wt0(667)/ 1.760587789374176d-09/,wt0(668)/-1.709724659722168d-09/,
     *wt0(669)/ 1.660330958617607d-09/,wt0(670)/-1.612364234479730d-09/,
     *wt0(671)/ 1.565783262147887d-09/,wt0(672)/-1.520548007450427d-09/,
     *wt0(673)/ 1.476619592797189d-09/,wt0(674)/-1.433960263766035d-09/,
     *wt0(675)/ 1.392533356654690d-09/,wt0(676)/-1.352303266970005d-09/,
     *wt0(677)/ 1.313235418827548d-09/,wt0(678)/-1.275296235235242d-09/,
     *wt0(679)/ 1.238453109235514d-09/,wt0(680)/-1.202674375881138d-09/,
     *wt0(681)/ 1.167929285020692d-09/,wt0(682)/-1.134187974870230d-09/,
     *wt0(683)/ 1.101421446348478d-09/,wt0(684)/-1.069601538153474d-09/
      data
     *wt0(685)/ 1.038700902559248d-09/,wt0(686)/-1.008692981911727d-09/,
     *wt0(687)/ 9.795519858036667d-10/,wt0(688)/-9.512528689090021d-10/,
     *wt0(689)/ 9.237713094575600d-10/,wt0(690)/-8.970836883316330d-10/,
     *wt0(691)/ 8.711670687664494d-10/,wt0(692)/-8.459991766370928d-10/,
     *wt0(693)/ 8.215583813149323d-10/,wt0(694)/-7.978236770771082d-10/,
     *wt0(695)/ 7.747746650530930d-10/,wt0(696)/-7.523915356928136d-10/,
     *wt0(697)/ 7.306550517412669d-10/,wt0(698)/-7.095465317049981d-10/,
     *wt0(699)/ 6.890478337962292d-10/,wt0(700)/-6.691413403408390d-10/,
     *wt0(701)/ 6.498099426367938d-10/,wt0(702)/-6.310370262500165d-10/,
     *wt0(703)/ 6.128064567350563d-10/,wt0(704)/-5.951025657682865d-10/,
     *wt0(705)/ 5.779101376817111d-10/,wt0(706)/-5.612143963858077d-10/,
     *wt0(707)/ 5.450009926701687d-10/,wt0(708)/-5.292559918710250d-10/,
     *wt0(709)/ 5.139658618950525d-10/,wt0(710)/-4.991174615891697d-10/,
     *wt0(711)/ 4.846980294463296d-10/,wt0(712)/-4.706951726376001d-10/,
     *wt0(713)/ 4.570968563611062d-10/,wt0(714)/-4.438913934986793d-10/,
     *wt0(715)/ 4.310674345713248d-10/,wt0(716)/-4.186139579848742d-10/,
     *wt0(717)/ 4.065202605574395d-10/,wt0(718)/-3.947759483205284d-10/,
     *wt0(719)/ 3.833709275859127d-10/,wt0(720)/-3.722953962705749d-10/,
     *wt0(721)/ 3.615398354722752d-10/,wt0(722)/-3.510950012884991d-10/
      data
     *wt0(723)/ 3.409519168717552d-10/,wt0(724)/-3.311018647143929d-10/,
     *wt0(725)/ 3.215363791563124d-10/,wt0(726)/-3.122472391091240d-10/,
     *wt0(727)/ 3.032264609905071d-10/,wt0(728)/-2.944662918626933d-10/,
     *wt0(729)/ 2.859592027691788d-10/,wt0(730)/-2.776978822639370d-10/,
     *wt0(731)/ 2.696752301275724d-10/,wt0(732)/-2.618843512650130d-10/,
     *wt0(733)/ 2.543185497794985d-10/,wt0(734)/-2.469713232177689d-10/,
     *wt0(735)/ 2.398363569815098d-10/,wt0(736)/-2.329075189002498d-10/,
     *wt0(737)/ 2.261788539610461d-10/,wt0(738)/-2.196445791904288d-10/,
     *wt0(739)/ 2.132990786842053d-10/,wt0(740)/-2.071368987808525d-10/,
     *wt0(741)/ 2.011527433743496d-10/,wt0(742)/-1.953414693624216d-10/,
     *wt0(743)/ 1.896980822262837d-10/,wt0(744)/-1.842177317380848d-10/,
     *wt0(745)/ 1.788957077923634d-10/,wt0(746)/-1.737274363579317d-10/,
     *wt0(747)/ 1.687084755467094d-10/,wt0(748)/-1.638345117961281d-10/,
     *wt0(749)/ 1.591013561618254d-10/,wt0(750)/-1.545049407174433d-10/,
     *wt0(751)/ 1.500413150584345d-10/,wt0(752)/-1.457066429068751d-10/,
     *wt0(753)/ 1.414971988143619d-10/,wt0(754)/-1.374093649601637d-10/,
     *wt0(755)/ 1.334396280418735d-10/,wt0(756)/-1.295845762558885d-10/,
     *wt0(757)/ 1.258408963651246d-10/,wt0(758)/-1.222053708514437d-10/,
     *wt0(759)/ 1.186748751503467d-10/,wt0(760)/-1.152463749655566d-10/
      data
     *wt0(761)/ 1.119169236611821d-10/,wt0(762)/-1.086836597292220d-10/,
     *wt0(763)/ 1.055438043302322d-10/,wt0(764)/-1.024946589050432d-10/,
     *wt0(765)/ 9.953360285547502d-11/,wt0(766)/-9.665809129205519d-11/,
     *wt0(767)/ 9.386565284680683d-11/,wt0(768)/-9.115388754922578d-11/,
     *wt0(769)/ 8.852046476362457d-11/,wt0(770)/-8.596312118607507d-11/,
     *wt0(771)/ 8.347965889923978d-11/,wt0(772)/-8.106794348344660d-11/,
     *wt0(773)/ 7.872590218244558d-11/,wt0(774)/-7.645152212241484d-11/,
     *wt0(775)/ 7.424284858301879d-11/,wt0(776)/-7.209798331980772d-11/,
     *wt0(777)/ 7.001508293831718d-11/,wt0(778)/-6.799235732265070d-11/,
     *wt0(779)/ 6.602806812690510d-11/,wt0(780)/-6.412052735069306d-11/,
     *wt0(781)/ 6.226809604991307d-11/,wt0(782)/-6.046918330329692d-11/,
     *wt0(783)/ 5.872224571634798d-11/,wt0(784)/-5.702578811834639d-11/,
     *wt0(785)/ 5.537836697680089d-11/,wt0(786)/-5.377860007122251d-11/,
     *wt0(787)/ 5.222519065363231d-11/,wt0(788)/-5.071698520520750d-11/,
     *wt0(789)/ 4.925310917157278d-11/,wt0(790)/-4.783328375532226d-11/,
     *wt0(791)/ 4.645856316442775d-11/,wt0(792)/-4.513304831495452d-11/,
     *wt0(793)/ 4.386786828115582d-11/,wt0(794)/-4.269042848860555d-11/,
     *wt0(795)/ 4.166589074066043d-11/,wt0(796)/-4.094706131972530d-11/,
     *wt0(797)/ 4.089025606285087d-11/,wt0(798)/-4.232439520086882d-11/
      data
     *wt0(799)/ 4.717597028661840d-11/,wt0(800)/-5.992051472297203d-11/,
     *wt0(801)/ 9.095360729014628d-11/
c-----wt1(i)=j1 hankel transform filter weights for i=1,801
      data
     *wt1(  1)/-2.377900110058238d-29/,wt1(  2)/ 7.972119293423644d-28/,
     *wt1(  3)/ 3.832566186338494d-28/,wt1(  4)/ 8.493403321392027d-28/,
     *wt1(  5)/ 8.232164371695149d-28/,wt1(  6)/ 1.100368717019849d-27/,
     *wt1(  7)/ 1.325428732821688d-27/,wt1(  8)/ 1.592246428628784d-27/,
     *wt1(  9)/ 1.998601965098251d-27/,wt1( 10)/ 2.368546581454427d-27/,
     *wt1( 11)/ 2.979879339217500d-27/,wt1( 12)/ 3.543101964370716d-27/,
     *wt1( 13)/ 4.425741130789444d-27/,wt1( 14)/ 5.318079447058835d-27/,
     *wt1( 15)/ 6.555619006665152d-27/,wt1( 16)/ 7.994776219905923d-27/,
     *wt1( 17)/ 9.707158322300619d-27/,wt1( 18)/ 1.200497993577129d-26/,
     *wt1( 19)/ 1.440665350963159d-26/,wt1( 20)/ 1.797046711269925d-26/,
     *wt1( 21)/ 2.145390397713542d-26/,wt1( 22)/ 2.681614737113612d-26/,
     *wt1( 23)/ 3.203783629081399d-26/,wt1( 24)/ 3.992219752668547d-26/,
     *wt1( 25)/ 4.793780377063212d-26/,wt1( 26)/ 5.935124745092138d-26/,
     *wt1( 27)/ 7.177334354418767d-26/,wt1( 28)/ 8.825226769587027d-26/,
     *wt1( 29)/ 1.073646008000832d-25/,wt1( 30)/ 1.313923424181352d-25/,
     *wt1( 31)/ 1.603838995307793d-25/,wt1( 32)/ 1.958633055912045d-25/,
     *wt1( 33)/ 2.393489595056268d-25/,wt1( 34)/ 2.921690713642289d-25/,
     *wt1( 35)/ 3.570455561097498d-25/,wt1( 36)/ 4.359157800465477d-25/,
     *wt1( 37)/ 5.325840709133780d-25/,wt1( 38)/ 6.503774912184682d-25/
      data
     *wt1( 39)/ 7.944611784834684d-25/,wt1( 40)/ 9.702980761516990d-25/,
     *wt1( 41)/ 1.185162186264197d-24/,wt1( 42)/ 1.447534903673577d-24/,
     *wt1( 43)/ 1.768047044621983d-24/,wt1( 44)/ 2.159464672432506d-24/,
     *wt1( 45)/ 2.637627922624120d-24/,wt1( 46)/ 3.221525590342443d-24/,
     *wt1( 47)/ 3.934898912364447d-24/,wt1( 48)/ 4.805930563105760d-24/,
     *wt1( 49)/ 5.870195744781909d-24/,wt1( 50)/ 7.169601400607034d-24/,
     *wt1( 51)/ 8.757287433857846d-24/,wt1( 52)/ 1.069583032926807d-23/,
     *wt1( 53)/ 1.306426797568440d-23/,wt1( 54)/ 1.595639689868334d-23/,
     *wt1( 55)/ 1.948949017485241d-23/,wt1( 56)/ 2.380425857645875d-23/,
     *wt1( 57)/ 2.907479573699533d-23/,wt1( 58)/ 3.551187751264232d-23/,
     *wt1( 59)/ 4.337441481530691d-23/,wt1( 60)/ 5.297756660121858d-23/,
     *wt1( 61)/ 6.470696377141974d-23/,wt1( 62)/ 7.903329008477021d-23/,
     *wt1( 63)/ 9.653141660115531d-23/,wt1( 64)/ 1.179038165380973d-22/,
     *wt1( 65)/ 1.440079827014168d-22/,wt1( 66)/ 1.758917653247009d-22/,
     *wt1( 67)/ 2.148347376559587d-22/,wt1( 68)/ 2.623996157283884d-22/,
     *wt1( 69)/ 3.204958075285148d-22/,wt1( 70)/ 3.914542160210567d-22/,
     *wt1( 71)/ 4.781235486108842d-22/,wt1( 72)/ 5.839810952393910d-22/,
     *wt1( 73)/ 7.132764830587759d-22/,wt1( 74)/ 8.711974589175415d-22/,
     *wt1( 75)/ 1.064083434862305d-21/,wt1( 76)/ 1.299673924385260d-21/
      data
     *wt1( 77)/ 1.587425911126261d-21/,wt1( 78)/ 1.938885693339251d-21/,
     *wt1( 79)/ 2.368161151141687d-21/,wt1( 80)/ 2.892477586360832d-21/,
     *wt1( 81)/ 3.532881273577120d-21/,wt1( 82)/ 4.315069519949315d-21/,
     *wt1( 83)/ 5.270439518682773d-21/,wt1( 84)/ 6.437327295395713d-21/,
     *wt1( 85)/ 7.862571844418139d-21/,wt1( 86)/ 9.603363814287760d-21/,
     *wt1( 87)/ 1.172957893208765d-20/,wt1( 88)/ 1.432653521781635d-20/,
     *wt1( 89)/ 1.749847566646075d-20/,wt1( 90)/ 2.137267893871469d-20/,
     *wt1( 91)/ 2.610465830070486d-20/,wt1( 92)/ 3.188429016442547d-20/,
     *wt1( 93)/ 3.894357412134842d-20/,wt1( 94)/ 4.756577134959383d-20/,
     *wt1( 95)/ 5.809698595129943d-20/,wt1( 96)/ 7.095979207053098d-20/,
     *wt1( 97)/ 8.667051908059619d-20/,wt1( 98)/ 1.058595695124327d-19/,
     *wt1( 99)/ 1.292972220966051d-19/,wt1(100)/ 1.579239186764568d-19/,
     *wt1(101)/ 1.928887914036691d-19/,wt1(102)/ 2.355947994119052d-19/,
     *wt1(103)/ 2.877562665708650d-19/,wt1(104)/ 3.514661357021616d-19/,
     *wt1(105)/ 4.292819113549590d-19/,wt1(106)/ 5.243258540458620d-19/,
     *wt1(107)/ 6.404133671732880d-19/,wt1(108)/ 7.822022466379286d-19/,
     *wt1(109)/ 9.553844929891425d-19/,wt1(110)/ 1.166908611051965d-18/,
     *wt1(111)/ 1.425266206309626d-18/,wt1(112)/ 1.740823055783105d-18/,
     *wt1(113)/ 2.126247365091890d-18/,wt1(114)/ 2.597002781312310d-18/
      data
     *wt1(115)/ 3.171988392298073d-18/,wt1(116)/ 3.874272813855638d-18/,
     *wt1(117)/ 4.732050718776532d-18/,wt1(118)/ 5.779735750225090d-18/,
     *wt1(119)/ 7.059390282316958d-18/,wt1(120)/ 8.622352349869446d-18/,
     *wt1(121)/ 1.053137301019041d-17/,wt1(122)/ 1.286303788941302d-17/,
     *wt1(123)/ 1.571096273151045d-17/,wt1(124)/ 1.918939713846172d-17/,
     *wt1(125)/ 2.343800282007509d-17/,wt1(126)/ 2.862721583723883d-17/,
     *wt1(127)/ 3.496539240998871d-17/,wt1(128)/ 4.270678642844850d-17/,
     *wt1(129)/ 5.216223744747131d-17/,wt1(130)/ 6.371103687921334d-17/,
     *wt1(131)/ 7.781691646330707d-17/,wt1(132)/ 9.504569536565720d-17/,
     *wt1(133)/ 1.160892016025935d-16/,wt1(134)/ 1.417915110596655d-16/,
     *wt1(135)/ 1.731847439832610d-16/,wt1(136)/ 2.115280706868521d-16/,
     *wt1(137)/ 2.583612876753947d-16/,wt1(138)/ 3.155627883343085d-16/,
     *wt1(139)/ 3.854297646671823d-16/,wt1(140)/ 4.707643426724462d-16/,
     *wt1(141)/ 5.749936655677889d-16/,wt1(142)/ 7.022978436839158d-16/,
     *wt1(143)/ 8.577897883808816d-16/,wt1(144)/ 1.047705221625211d-15/,
     *wt1(145)/ 1.279672050432180d-15/,wt1(146)/ 1.562992451553929d-15/,
     *wt1(147)/ 1.909046462722709d-15/,wt1(148)/ 2.331710624424154d-15/,
     *wt1(149)/ 2.847962809294981d-15/,wt1(150)/ 3.478503311969800d-15/,
     *wt1(151)/ 4.248661490032314d-15/,wt1(152)/ 5.189316858234294d-15/
      data
     *wt1(153)/ 6.338258511920505d-15/,wt1(154)/ 7.741550588468052d-15/,
     *wt1(155)/ 9.455571172640662d-15/,wt1(156)/ 1.154903563040828d-14/,
     *wt1(157)/ 1.410605553114248d-14/,wt1(158)/ 1.722913542283047d-14/,
     *wt1(159)/ 2.104376349301754d-14/,wt1(160)/ 2.570284789891445d-14/,
     *wt1(161)/ 3.139360843035402d-14/,wt1(162)/ 3.834414037586157d-14/,
     *wt1(163)/ 4.683376407804295d-14/,wt1(164)/ 5.720273100006960d-14/,
     *wt1(165)/ 6.986777175159004d-14/,wt1(166)/ 8.533643955954289d-14/,
     *wt1(167)/ 1.042304766775628d-13/,wt1(168)/ 1.273069965554705d-13/,
     *wt1(169)/ 1.554936139359936d-13/,wt1(170)/ 1.899197032959808d-13/,
     *wt1(171)/ 2.319692366825060d-13/,wt1(172)/ 2.833268748946759d-13/,
     *wt1(173)/ 3.460574729307441d-13/,wt1(174)/ 4.226739834760657d-13/,
     *wt1(175)/ 5.162571427996572d-13/,wt1(176)/ 6.305554147709428d-13/,
     *wt1(177)/ 7.701652476024983d-13/,wt1(178)/ 9.406780256740293d-13/,
     *wt1(179)/ 1.148951682734905d-12/,wt1(180)/ 1.403326528650366d-12/,
     *wt1(181)/ 1.714034726435434d-12/,wt1(182)/ 2.093516885212530d-12/,
     *wt1(183)/ 2.557039701256063d-12/,wt1(184)/ 3.123159736572392d-12/,
     *wt1(185)/ 3.814655555107956d-12/,wt1(186)/ 4.659206104975550d-12/,
     *wt1(187)/ 5.690798281976283d-12/,wt1(188)/ 6.950717591317659d-12/,
     *wt1(189)/ 8.489674870123959d-12/,wt1(190)/ 1.036925035218241d-11/
      data
     *wt1(191)/ 1.266510893213373d-11/,wt1(192)/ 1.546910089464080d-11/,
     *wt1(193)/ 1.889412592258518d-11/,wt1(194)/ 2.307718220992490d-11/,
     *wt1(195)/ 2.818672942203327d-11/,wt1(196)/ 3.442710316037366d-11/,
     *wt1(197)/ 4.204966817018048d-11/,wt1(198)/ 5.135919134271227d-11/,
     *wt1(199)/ 6.273074786797763d-11/,wt1(200)/ 7.661889201263122d-11/,
     *wt1(201)/ 9.358330171196327d-11/,wt1(202)/ 1.143019267753882d-10/,
     *wt1(203)/ 1.396099167826153d-10/,wt1(204)/ 1.705183920059801d-10/,
     *wt1(205)/ 2.082735788845278d-10/,wt1(206)/ 2.543834767842491d-10/,
     *wt1(207)/ 3.107077590612055d-10/,wt1(208)/ 3.794954396016648d-10/,
     *wt1(209)/ 4.635216515055463d-10/,wt1(210)/ 5.661404892806104d-10/,
     *wt1(211)/ 6.914932735763221d-10/,wt1(212)/ 8.445820788220271d-10/,
     *wt1(213)/ 1.031587101304291d-09/,wt1(214)/ 1.259967952077770d-09/,
     *wt1(215)/ 1.538947681019300d-09/,wt1(216)/ 1.879650592174711d-09/,
     *wt1(217)/ 2.295841053082368d-09/,wt1(218)/ 2.804108039269473d-09/,
     *wt1(219)/ 3.424993797807946d-09/,wt1(220)/ 4.183235823312564d-09/,
     *wt1(221)/ 5.109492567813164d-09/,wt1(222)/ 6.240651650854630d-09/,
     *wt1(223)/ 7.622470723630660d-09/,wt1(224)/ 9.309953703014609d-09/,
     *wt1(225)/ 1.137139562283767d-08/,wt1(226)/ 1.388881160360346d-08/,
     *wt1(227)/ 1.696413754044688d-08/,wt1(228)/ 2.071966056701105d-08/
      data
     *wt1(229)/ 2.530753298545001d-08/,wt1(230)/ 3.091008276019684d-08/,
     *wt1(231)/ 3.775442398336976d-08/,wt1(232)/ 4.611239491146397d-08/,
     *wt1(233)/ 5.632301501798812d-08/,wt1(234)/ 6.879156104001002d-08/,
     *wt1(235)/ 8.402411522137500d-08/,wt1(236)/ 1.026248703632157d-07/,
     *wt1(237)/ 1.253493263032310d-07/,wt1(238)/ 1.530981849317344d-07/,
     *wt1(239)/ 1.869993328577540d-07/,wt1(240)/ 2.283954331971998d-07/,
     *wt1(241)/ 2.789703820965002d-07/,wt1(242)/ 3.407255717347445d-07/,
     *wt1(243)/ 4.161751163160901d-07/,wt1(244)/ 5.083021657140526d-07/,
     *wt1(245)/ 6.208605585127383d-07/,wt1(246)/ 7.582965471628845d-07/,
     *wt1(247)/ 9.262152932479593d-07/,wt1(248)/ 1.131243348429578d-06/,
     *wt1(249)/ 1.381750677989317d-06/,wt1(250)/ 1.687612654971312d-06/,
     *wt1(251)/ 2.061328504830195d-06/,wt1(252)/ 2.517614221705925d-06/,
     *wt1(253)/ 3.075136493715169d-06/,wt1(254)/ 3.755823026355034d-06/,
     *wt1(255)/ 4.587552788474290d-06/,wt1(256)/ 5.602996673366901d-06/,
     *wt1(257)/ 6.843794758360242d-06/,wt1(258)/ 8.358620144706847d-06/,
     *wt1(259)/ 1.020967008292656d-05/,wt1(260)/ 1.246945025422475d-05/,
     *wt1(261)/ 1.523086887102219d-05/,wt1(262)/ 1.860192083272764d-05/,
     *wt1(263)/ 2.272139580071132d-05/,wt1(264)/ 2.775012613215501d-05/,
     *wt1(265)/ 3.389544779882663d-05/,wt1(266)/ 4.139684113647092d-05/
      data
     *wt1(267)/ 5.056404631634212d-05/,wt1(268)/ 6.175357572143610d-05/,
     *wt1(269)/ 7.542814171915086d-05/,wt1(270)/ 9.211833800758343d-05/,
     *wt1(271)/ 1.125153468149895d-04/,wt1(272)/ 1.374085645976000d-04/,
     *wt1(273)/ 1.678302935465055d-04/,wt1(274)/ 2.049544027080777d-04/,
     *wt1(275)/ 2.503222397064943d-04/,wt1(276)/ 3.056782429102864d-04/,
     *wt1(277)/ 3.733227575027256d-04/,wt1(278)/ 4.558456931739724d-04/,
     *wt1(279)/ 5.566775872841220d-04/,wt1(280)/ 6.796584276296773d-04/,
     *wt1(281)/ 8.298982533418059d-04/,wt1(282)/ 1.013079388381259d-03/,
     *wt1(283)/ 1.236799879071809d-03/,wt1(284)/ 1.509444128879166d-03/,
     *wt1(285)/ 1.842278291863989d-03/,wt1(286)/ 2.247620855799141d-03/,
     *wt1(287)/ 2.742106207595588d-03/,wt1(288)/ 3.343716405309918d-03/,
     *wt1(289)/ 4.076852784788600d-03/,wt1(290)/ 4.967507347193811d-03/,
     *wt1(291)/ 6.051117627655160d-03/,wt1(292)/ 7.364681573940744d-03/,
     *wt1(293)/ 8.958863753805689d-03/,wt1(294)/ 1.088505212285743d-02/,
     *wt1(295)/ 1.321388763870185d-02/,wt1(296)/ 1.601387913655473d-02/,
     *wt1(297)/ 1.937949548174205d-02/,wt1(298)/ 2.339556612555552d-02/,
     *wt1(299)/ 2.817934173583488d-02/,wt1(300)/ 3.382079060619734d-02/,
     *wt1(301)/ 4.044464120111045d-02/,wt1(302)/ 4.810974644754662d-02/,
     *wt1(303)/ 5.689934182811928d-02/,wt1(304)/ 6.675168423575902d-02/
      data
     *wt1(305)/ 7.759114407311523d-02/,wt1(306)/ 8.904761297446566d-02/,
     *wt1(307)/ 1.006540591855666d-01/,wt1(308)/ 1.114010542073269d-01/,
     *wt1(309)/ 1.200693293287907d-01/,wt1(310)/ 1.245859746888521d-01/,
     *wt1(311)/ 1.226876554954843d-01/,wt1(312)/ 1.111594205067191d-01/,
     *wt1(313)/ 8.731995478664586d-02/,wt1(314)/ 4.846500030224664d-02/,
     *wt1(315)/-4.956156637375708d-03/,wt1(316)/-6.947108182097863d-02/,
     *wt1(317)/-1.335126953345830d-01/,wt1(318)/-1.795750728565332d-01/,
     *wt1(319)/-1.819110629504884d-01/,wt1(320)/-1.209847646201776d-01/,
     *wt1(321)/ 5.791427226916448d-03/,wt1(322)/ 1.534812881363312d-01/,
     *wt1(323)/ 2.352209731069856d-01/,wt1(324)/ 1.515040026640386d-01/,
     *wt1(325)/-8.279085782145131d-02/,wt1(326)/-2.711543104015952d-01/,
     *wt1(327)/-1.331279806376289d-01/,wt1(328)/ 2.242058279665979d-01/,
     *wt1(329)/ 2.304292248782090d-01/,wt1(330)/-2.572813803130152d-01/,
     *wt1(331)/-1.475549438192476d-01/,wt1(332)/ 3.611935526731244d-01/,
     *wt1(333)/-2.435371552288723d-01/,wt1(334)/ 3.811454725614711d-02/,
     *wt1(335)/ 8.273650321459338d-02/,wt1(336)/-1.091467781119748d-01/,
     *wt1(337)/ 8.964328522861212d-02/,wt1(338)/-6.097788815932206d-02/,
     *wt1(339)/ 3.780123181837066d-02/,wt1(340)/-2.243862121105374d-02/,
     *wt1(341)/ 1.315627818991557d-02/,wt1(342)/-7.789693311677149d-03/
      data
     *wt1(343)/ 4.737032461705975d-03/,wt1(344)/-2.996529369596568d-03/,
     *wt1(345)/ 1.987854486046977d-03/,wt1(346)/-1.386982713556120d-03/,
     *wt1(347)/ 1.015740041573320d-03/,wt1(348)/-7.763183089745113d-04/,
     *wt1(349)/ 6.146040263274792d-04/,wt1(350)/-5.002069052247894d-04/,
     *wt1(351)/ 4.157061680153099d-04/,wt1(352)/-3.508644775971589d-04/,
     *wt1(353)/ 2.994936403367197d-04/,wt1(354)/-2.577357341490762d-04/,
     *wt1(355)/ 2.231048844827660d-04/,wt1(356)/-1.939428674265101d-04/,
     *wt1(357)/ 1.691035669098717d-04/,wt1(358)/-1.477663053511028d-04/,
     *wt1(359)/ 1.293230643516167d-04/,wt1(360)/-1.133088634015300d-04/,
     *wt1(361)/ 9.935787087369453d-05/,wt1(362)/-8.717517169548774d-05/,
     *wt1(363)/ 7.651820346514748d-05/,wt1(364)/-6.718422692060818d-05/,
     *wt1(365)/ 5.900162548547053d-05/,wt1(366)/-5.182370043717230d-05/,
     *wt1(367)/ 4.552413144789780d-05/,wt1(368)/-3.999355821930881d-05/,
     *wt1(369)/ 3.513692405823457d-05/,wt1(370)/-3.087135626461035d-05/,
     *wt1(371)/ 2.712444449271282d-05/,wt1(372)/-2.383282162264958d-05/,
     *wt1(373)/ 2.094097402272335d-05/,wt1(374)/-1.840022834334601d-05/,
     *wt1(375)/ 1.616788131114735d-05/,wt1(376)/-1.420645034584182d-05/,
     *wt1(377)/ 1.248302592309649d-05/,wt1(378)/-1.096870811434504d-05/,
     *wt1(379)/ 9.638113747637144d-06/,wt1(380)/-8.468945120923146d-06/
      data
     *wt1(381)/ 7.441613142631071d-06/,wt1(382)/-6.538907689594796d-06/,
     *wt1(383)/ 5.745708304710312d-06/,wt1(384)/-5.048729864716120d-06/,
     *wt1(385)/ 4.436299286415519d-06/,wt1(386)/-3.898159774165185d-06/,
     *wt1(387)/ 3.425299055300640d-06/,wt1(388)/-3.009798361246221d-06/,
     *wt1(389)/ 2.644699605727544d-06/,wt1(390)/-2.323888735130234d-06/,
     *wt1(391)/ 2.041993364758909d-06/,wt1(392)/-1.794292864274530d-06/,
     *wt1(393)/ 1.576639283761010d-06/,wt1(394)/-1.385387831873578d-06/,
     *wt1(395)/ 1.217335820899756d-06/,wt1(396)/-1.069669061011715d-06/,
     *wt1(397)/ 9.399147533431300d-07/,wt1(398)/-8.259000682217990d-07/,
     *wt1(399)/ 7.257157398665060d-07/,wt1(400)/-6.376840939530841d-07/,
     *wt1(401)/ 5.603309658957991d-07/,wt1(402)/-4.923610193285461d-07/,
     *wt1(403)/ 4.326360479415582d-07/,wt1(404)/-3.801559094539231d-07/,
     *wt1(405)/ 3.340417795380239d-07/,wt1(406)/-2.935214398955969d-07/,
     *wt1(407)/ 2.579163472426386d-07/,wt1(408)/-2.266302675538730d-07/,
     *wt1(409)/ 1.991392901933410d-07/,wt1(410)/-1.749830560542944d-07/,
     *wt1(411)/ 1.537570500204901d-07/,wt1(412)/-1.351058266793548d-07/,
     *wt1(413)/ 1.187170569339987d-07/,wt1(414)/-1.043162975874935d-07/,
     *wt1(415)/ 9.166239627647569d-08/,wt1(416)/-8.054345360702865d-08/,
     *wt1(417)/ 7.077327442848185d-08/,wt1(418)/-6.218824944698919d-08/
      data
     *wt1(419)/ 5.464461550533503d-08/,wt1(420)/-4.801604844364623d-08/,
     *wt1(421)/ 4.219154778426137d-08/,wt1(422)/-3.707357778095240d-08/,
     *wt1(423)/ 3.257643398813585d-08/,wt1(424)/-2.862480813390940d-08/,
     *wt1(425)/ 2.515252712714758d-08/,wt1(426)/-2.210144493659511d-08/,
     *wt1(427)/ 1.942046882050334d-08/,wt1(428)/-1.706470371165871d-08/,
     *wt1(429)/ 1.499470045020363d-08/,wt1(430)/-1.317579520863365d-08/,
     *wt1(431)/ 1.157752900730261d-08/,wt1(432)/-1.017313762883149d-08/,
     *wt1(433)/ 8.939103427741988d-09/,wt1(434)/-7.854761525494339d-09/,
     *wt1(435)/ 6.901953769038680d-09/,wt1(436)/-6.064724650448546d-09/,
     *wt1(437)/ 5.329054110861765d-09/,wt1(438)/-4.682622765482983d-09/,
     *wt1(439)/ 4.114605612289784d-09/,wt1(440)/-3.615490760147216d-09/,
     *wt1(441)/ 3.176920140880805d-09/,wt1(442)/-2.791549544501460d-09/,
     *wt1(443)/ 2.452925636291220d-09/,wt1(444)/-2.155377892876758d-09/,
     *wt1(445)/ 1.893923644634260d-09/,wt1(446)/-1.664184635389759d-09/,
     *wt1(447)/ 1.462313704613378d-09/,wt1(448)/-1.284930364574083d-09/,
     *wt1(449)/ 1.129064192077534d-09/,wt1(450)/-9.921050860675488d-10/,
     *wt1(451)/ 8.717595588633507d-10/,wt1(452)/-7.660123299045589d-10/,
     *wt1(453)/ 6.730925786777026d-10/,wt1(454)/-5.914442910883912d-10/,
     *wt1(455)/ 5.197002025874266d-10/,wt1(456)/-4.566589020735369d-10/
      data
     *wt1(457)/ 4.012647133825055d-10/,wt1(458)/-3.525900173049568d-10/,
     *wt1(459)/ 3.098197178898629d-10/,wt1(460)/-2.722375929324441d-10/,
     *wt1(461)/ 2.392143002251476d-10/,wt1(462)/-2.101968387708856d-10/,
     *wt1(463)/ 1.846992883946521d-10/,wt1(464)/-1.622946726187645d-10/,
     *wt1(465)/ 1.426078085707891d-10/,wt1(466)/-1.253090242443995d-10/,
     *wt1(467)/ 1.101086379038488d-10/,wt1(468)/-9.675210714907983d-11/,
     *wt1(469)/ 8.501576639446501d-11/,wt1(470)/-7.470308140061063d-11/,
     *wt1(471)/ 6.564135815534322d-11/,wt1(472)/-5.767885098561786d-11/,
     *wt1(473)/ 5.068222145975428d-11/,wt1(474)/-4.453430552590912d-11/,
     *wt1(475)/ 3.913215150436401d-11/,wt1(476)/-3.438529608259476d-11/,
     *wt1(477)/ 3.021424943923932d-11/,wt1(478)/-2.654916412468056d-11/,
     *wt1(479)/ 2.332866540856460d-11/,wt1(480)/-2.049882351089608d-11/,
     *wt1(481)/ 1.801225050660065d-11/,wt1(482)/-1.582730677850963d-11/,
     *wt1(483)/ 1.390740372895185d-11/,wt1(484)/-1.222039107416370d-11/,
     *wt1(485)/ 1.073801846245205d-11/,wt1(486)/-9.435462400371624d-12/,
     *wt1(487)/ 8.290910564095319d-12/,wt1(488)/-7.285196534694753d-12/,
     *wt1(489)/ 6.401478841143664d-12/,wt1(490)/-5.624958936714699d-12/,
     *wt1(491)/ 4.942633385892820d-12/,wt1(492)/-4.343076111680976d-12/,
     *wt1(493)/ 3.816247056854583d-12/,wt1(494)/-3.353324055246731d-12/
      data
     *wt1(495)/ 2.946555097697522d-12/,wt1(496)/-2.589128518642891d-12/,
     *wt1(497)/ 2.275058929421262d-12/,wt1(498)/-1.999086988193079d-12/,
     *wt1(499)/ 1.756591328117270d-12/,wt1(500)/-1.543511168961402d-12/,
     *wt1(501)/ 1.356278316168970d-12/,wt1(502)/-1.191757408630524d-12/,
     *wt1(503)/ 1.047193414582761d-12/,wt1(504)/-9.201654964381461d-13/,
     *wt1(505)/ 8.085464719703853d-13/,wt1(506)/-7.104671929827940d-13/,
     *wt1(507)/ 6.242852449483966d-13/,wt1(508)/-5.485574434836370d-13/,
     *wt1(509)/ 4.820156670945356d-13/,wt1(510)/-4.235456214910329d-13/,
     *wt1(511)/ 3.721681798548054d-13/,wt1(512)/-3.270229865892332d-13/,
     *wt1(513)/ 2.873540499880031d-13/,wt1(514)/-2.524970825616983d-13/,
     *wt1(515)/ 2.218683770237039d-13/,wt1(516)/-1.949550316531413d-13/,
     *wt1(517)/ 1.713063613517743d-13/,wt1(518)/-1.505263505679493d-13/,
     *wt1(519)/ 1.322670217058939d-13/,wt1(520)/-1.162226079680830d-13/,
     *wt1(521)/ 1.021244330498886d-13/,wt1(522)/-8.973641194353509d-14/,
     *wt1(523)/ 7.885109750931469d-14/,wt1(524)/-6.928620661062564d-14/,
     *wt1(525)/ 6.088156763982168d-14/,wt1(526)/-5.349643831874401d-14/,
     *wt1(527)/ 4.700714885865928d-14/,wt1(528)/-4.130503101260246d-14/,
     *wt1(529)/ 3.629459834036019d-14/,wt1(530)/-3.189194721306437d-14/,
     *wt1(531)/ 2.802335178098858d-14/,wt1(532)/-2.462402937628328d-14/
      data
     *wt1(533)/ 2.163705567637205d-14/,wt1(534)/-1.901241146150821d-14/,
     *wt1(535)/ 1.670614500365970d-14/,wt1(536)/-1.467963606029139d-14/,
     *wt1(537)/ 1.289894914807883d-14/,wt1(538)/-1.133426526661138d-14/,
     *wt1(539)/ 9.959382555826870d-15/,wt1(540)/-8.751277525286314d-15/,
     *wt1(541)/ 7.689719507742115d-15/,wt1(542)/-6.756931880732829d-15/,
     *wt1(543)/ 5.937294383089198d-15/,wt1(544)/-5.217081541398457d-15/,
     *wt1(545)/ 4.584232826172386d-15/,wt1(546)/-4.028150688809716d-15/,
     *wt1(547)/ 3.539523097325829d-15/,wt1(548)/-3.110167599068186d-15/,
     *wt1(549)/ 2.732894299122300d-15/,wt1(550)/-2.401385459874979d-15/,
     *wt1(551)/ 2.110089705537303d-15/,wt1(552)/-1.854129059999304d-15/,
     *wt1(553)/ 1.629217261290832d-15/,wt1(554)/-1.431587984759346d-15/,
     *wt1(555)/ 1.257931773006046d-15/,wt1(556)/-1.105340616423372d-15/,
     *wt1(557)/ 9.712592562913251d-16/,wt1(558)/-8.534423949641643d-16/,
     *wt1(559)/ 7.499170965983152d-16/,wt1(560)/-6.589497487925016d-16/,
     *wt1(561)/ 5.790170318870357d-16/,wt1(562)/-5.087804097802586d-16/,
     *wt1(563)/ 4.470637150906226d-16/,wt1(564)/-3.928334533103688d-16/,
     *wt1(565)/ 3.451814961285280d-16/,wt1(566)/-3.033098741093958d-16/,
     *wt1(567)/ 2.665174140678718d-16/,wt1(568)/-2.341879973739602d-16/,
     *wt1(569)/ 2.057802425625493d-16/,wt1(570)/-1.808184394757190d-16/
      data
     *wt1(571)/ 1.588845831226758d-16/,wt1(572)/-1.396113738580090d-16/,
     *wt1(573)/ 1.226760666607451d-16/,wt1(574)/-1.077950665155628d-16/,
     *wt1(575)/ 9.471917939160380d-17/,wt1(576)/-8.322943929278698d-17/,
     *wt1(577)/ 7.313344150029316d-17/,wt1(578)/-6.426212060448496d-17/,
     *wt1(579)/ 5.646691937188211d-17/,wt1(580)/-4.961730103765172d-17/,
     *wt1(581)/ 4.359856336499476d-17/,wt1(582)/-3.830991786613391d-17/,
     *wt1(583)/ 3.366280201995165d-17/,wt1(584)/-2.957939622304163d-17/,
     *wt1(585)/ 2.599132063935437d-17/,wt1(586)/-2.283849012615860d-17/,
     *wt1(587)/ 2.006810806115271d-17/,wt1(588)/-1.763378222157133d-17/,
     *wt1(589)/ 1.549474790998944d-17/,wt1(590)/-1.361518531745643d-17/,
     *wt1(591)/ 1.196361969267083d-17/,wt1(592)/-1.051239427255871d-17/,
     *wt1(593)/ 9.237207148064648d-18/,wt1(594)/-8.116704309607398d-18/,
     *wt1(595)/ 7.132122057424711d-18/,wt1(596)/-6.266972788687375d-18/,
     *wt1(597)/ 5.506768899643017d-18/,wt1(598)/-4.838780179316216d-18/,
     *wt1(599)/ 4.251820632131413d-18/,wt1(600)/-3.736061159634172d-18/,
     *wt1(601)/ 3.282864964518770d-18/,wt1(602)/-2.884642920675208d-18/,
     *wt1(603)/ 2.534726487315339d-18/,wt1(604)/-2.227256039024718d-18/,
     *wt1(605)/ 1.957082741736910d-18/,wt1(606)/-1.719682331465223d-18/,
     *wt1(607)/ 1.511079351949921d-18/,wt1(608)/-1.327780582522198d-18/
      data
     *wt1(609)/ 1.166716541382233d-18/,wt1(610)/-1.025190084716866d-18/,
     *wt1(611)/ 9.008312409083032d-19/,wt1(612)/-7.915575234972099d-19/,
     *wt1(613)/ 6.955390583093758d-19/,wt1(614)/-6.111679407722684d-19/,
     *wt1(615)/ 5.370313102870695d-19/,wt1(616)/-4.718876907676534d-19/,
     *wt1(617)/ 4.146462011511975d-19/,wt1(618)/-3.643482877427443d-19/,
     *wt1(619)/ 3.201516724850009d-19/,wt1(620)/-2.813162483527929d-19/,
     *wt1(621)/ 2.471916856791182d-19/,wt1(622)/-2.172065418717432d-19/,
     *wt1(623)/ 1.908586921555019d-19/,wt1(624)/-1.677069210972190d-19/,
     *wt1(625)/ 1.473635341074691d-19/,wt1(626)/-1.294878651921690d-19/,
     *wt1(627)/ 1.137805722345515d-19/,wt1(628)/-9.997862427642855d-20/,
     *wt1(629)/ 8.785089685693581d-20/,wt1(630)/-7.719430165048124d-20/,
     *wt1(631)/ 6.783038559324186d-20/,wt1(632)/-5.960234254861157d-20/,
     *wt1(633)/ 5.237238746867620d-20/,wt1(634)/-4.601944907799879d-20/,
     *wt1(635)/ 4.043714244041481d-20/,wt1(636)/-3.553198745789200d-20/,
     *wt1(637)/ 3.122184346997662d-20/,wt1(638)/-2.743453374142427d-20/,
     *wt1(639)/ 2.410663680492846d-20/,wt1(640)/-2.118242441910240d-20/,
     *wt1(641)/ 1.861292835641397d-20/,wt1(642)/-1.635512039302220d-20/,
     *wt1(643)/ 1.437119176831859d-20/,wt1(644)/-1.262792004764438d-20/,
     *wt1(645)/ 1.109611278531706d-20/,wt1(646)/-9.750118671903616d-21/
      data
     *wt1(647)/ 8.567397980805881d-21/,wt1(648)/-7.528145122423124d-21/,
     *wt1(649)/ 6.614956985587211d-21/,wt1(650)/-5.812541511365789d-21/,
     *wt1(651)/ 5.107461617652140d-21/,wt1(652)/-4.487910185482492d-21/,
     *wt1(653)/ 3.943512339081440d-21/,wt1(654)/-3.465151709308229d-21/,
     *wt1(655)/ 3.044817771736880d-21/,wt1(656)/-2.675471703023044d-21/,
     *wt1(657)/ 2.350928509369901d-21/,wt1(658)/-2.065753454081455d-21/,
     *wt1(659)/ 1.815171050414263d-21/,wt1(660)/-1.594985094847466d-21/,
     *wt1(661)/ 1.401508399670877d-21/,wt1(662)/-1.231501046949432d-21/,
     *wt1(663)/ 1.082116130537259d-21/,wt1(664)/-9.508520796066248d-22/,
     *wt1(665)/ 8.355107675059663d-22/,wt1(666)/-7.341607053512435d-22/,
     *wt1(667)/ 6.451047024956878d-22/,wt1(668)/-5.668514489455267d-22/,
     *wt1(669)/ 4.980905412312518d-22/,wt1(670)/-4.376705343564336d-22/,
     *wt1(671)/ 3.845796559678223d-22/,wt1(672)/-3.379288625754182d-22/,
     *wt1(673)/ 2.969369534739024d-22/,wt1(674)/-2.609174906141295d-22/,
     *wt1(675)/ 2.292673041856001d-22/,wt1(676)/-2.014563919875952d-22/,
     *wt1(677)/ 1.770190438740656d-22/,wt1(678)/-1.555460422525817d-22/,
     *wt1(679)/ 1.366778079672955d-22/,wt1(680)/-1.200983777985825d-22/,
     *wt1(681)/ 1.055301137906078d-22/,wt1(682)/-9.272905558580399d-23/,
     *wt1(683)/ 8.148083664639517d-23/,wt1(684)/-7.159709489187570d-23/
      data
     *wt1(685)/ 6.291231754727859d-23/,wt1(686)/-5.528106802589832d-23/,
     *wt1(687)/ 4.857554932085610d-23/,wt1(688)/-4.268346388064600d-23/,
     *wt1(689)/ 3.750613425528731d-23/,wt1(690)/-3.295685226039235d-23/,
     *wt1(691)/ 2.895942787876105d-23/,wt1(692)/-2.544691315369192d-23/,
     *wt1(693)/ 2.236048039475426d-23/,wt1(694)/-1.964843688844607d-23/,
     *wt1(695)/ 1.726535979915737d-23/,wt1(696)/-1.517133602902566d-23/,
     *wt1(697)/ 1.333129356639782d-23/,wt1(698)/-1.171441336528463d-23/,
     *wt1(699)/ 1.029361299133160d-23/,wt1(700)/-9.045094266181886d-24/,
     *wt1(701)/ 7.947947207266389d-24/,wt1(702)/-6.983802597398395d-24/,
     *wt1(703)/ 6.136526075969871d-24/,wt1(704)/-5.391947793289710d-24/,
     *wt1(705)/ 4.737623457434066d-24/,wt1(706)/-4.162624052945447d-24/,
     *wt1(707)/ 3.657351341563861d-24/,wt1(708)/-3.213375054547230d-24/,
     *wt1(709)/ 2.823287629843971d-24/,wt1(710)/-2.480574676640252d-24/,
     *wt1(711)/ 2.179502108849238d-24/,wt1(712)/-1.915020796639139d-24/,
     *wt1(713)/ 1.682686753114865d-24/,wt1(714)/-1.478592609598710d-24/,
     *wt1(715)/ 1.299306570857684d-24/,wt1(716)/-1.141816901768704d-24/,
     *wt1(717)/ 1.003481245022785d-24/,wt1(718)/-8.819804660957208d-25/,
     *wt1(719)/ 7.752768477813106d-25/,wt1(720)/-6.815767858002522d-25/,
     *wt1(721)/ 5.992985900665248d-25/,wt1(722)/-5.270459063734392d-25/
      data
     *wt1(723)/ 4.635862411179793d-25/,wt1(724)/-4.078328923510551d-25/,
     *wt1(725)/ 3.588287019442058d-25/,wt1(726)/-3.157312905756895d-25/,
     *wt1(727)/ 2.778003181153750d-25/,wt1(728)/-2.443872105993099d-25/,
     *wt1(729)/ 2.149270912330287d-25/,wt1(730)/-1.889318548953116d-25/,
     *wt1(731)/ 1.659828344765530d-25/,wt1(732)/-1.457219798327836d-25/,
     *wt1(733)/ 1.278418979568063d-25/,wt1(734)/-1.120761165538443d-25/,
     *wt1(735)/ 9.819063583398694d-26/,wt1(736)/-8.597706957547351d-26/,
     *wt1(737)/ 7.524761802539771d-26/,wt1(738)/-6.583229692872437d-26/,
     *wt1(739)/ 5.757815308418353d-26/,wt1(740)/-5.034908817517119d-26/,
     *wt1(741)/ 4.402475606598968d-26/,wt1(742)/-3.849819197101101d-26/,
     *wt1(743)/ 3.367329111558304d-26/,wt1(744)/-2.946341179122434d-26/,
     *wt1(745)/ 2.579121784723996d-26/,wt1(746)/-2.258873118181162d-26/,
     *wt1(747)/ 1.979667076757990d-26/,wt1(748)/-1.736329567227616d-26/, 
     *wt1(749)/ 1.524356928842631d-26/,wt1(750)/-1.339878795333426d-26/,
     *wt1(751)/ 1.179604582048012d-26/,wt1(752)/-1.040714749090827d-26/,
     *wt1(753)/ 9.207194488904074d-27/,wt1(754)/-8.173124929559144d-27/,
     *wt1(755)/ 7.282515970348745d-27/,wt1(756)/-6.513100644712918d-27/,
     *wt1(757)/ 5.843030126248539d-27/,wt1(758)/-5.251393008101120d-27/,
     *wt1(759)/ 4.718791967930222d-27/,wt1(760)/-4.228406914229241d-27/
      data
     *wt1(761)/ 3.767695474750273d-27/,wt1(762)/-3.329913400113583d-27/,
     *wt1(763)/ 2.914300433863011d-27/,wt1(764)/-2.524366823406565d-27/,
     *wt1(765)/ 2.164690315762398d-27/,wt1(766)/-1.837673402805216d-27/,
     *wt1(767)/ 1.542301736687901d-27/,wt1(768)/-1.275894419058812d-27/,
     *wt1(769)/ 1.037166283039543d-27/,wt1(770)/-8.274817350070647d-28/,
     *wt1(771)/ 6.490942784063843d-28/,wt1(772)/-5.023466637656253d-28/,
     *wt1(773)/ 3.845211097313078d-28/,wt1(774)/-2.909163620505845d-28/,
     *wt1(775)/ 2.167800946716682d-28/,wt1(776)/-1.587457861911164d-28/,
     *wt1(777)/ 1.153739423502031d-28/,wt1(778)/-8.690036418978894d-29/,
     *wt1(779)/ 7.435572366263697d-29/,wt1(780)/-7.847142831542739d-29/,
     *wt1(781)/ 9.894985212873060d-29/,wt1(782)/-1.342720280592225d-28/,
     *wt1(783)/ 1.817424273402145d-28/,wt1(784)/-2.376361960950870d-28/,
     *wt1(785)/ 2.977819810711053d-28/,wt1(786)/-3.586836968065507d-28/,
     *wt1(787)/ 4.184476666502118d-28/,wt1(788)/-4.766085106688037d-28/,
     *wt1(789)/ 5.328777445258146d-28/,wt1(790)/-5.858824985834212d-28/,
     *wt1(791)/ 6.329409859706737d-28/,wt1(792)/-6.709268971860476d-28/,
     *wt1(793)/ 6.975416172179417d-28/,wt1(794)/-7.122929261229363d-28/,
     *wt1(795)/ 7.167821355397955d-28/,wt1(796)/-7.142209321177599d-28/,
     *wt1(797)/ 7.084785743816030d-28/,wt1(798)/-7.032486471839363d-28/
      data
     *wt1(799)/ 7.017650479151096d-28/,wt1(800)/-7.069608472164614d-28/,
     *wt1(801)/ 7.214920505613761d-28/
c
      nofun1=0
c-----error checks
      if(mb.lt.1.or.nrel.lt.1.or.bmax.le.0.0.or.
     *   mb.gt.maxdim.or.nrel.gt.maxrel) go to 9999
      y=bmax*er**(mb-1)
      if(y.le.0.0d0) go to 9999
      ierr=0
c-----initialize related convolution within lagged convolution loops
      do 10 i=1,801
   10   key(i)=0
      mb1=mb+1
      lag=-1
c-----preset initial filter abscissa for starting bmax, the argument
c     used in the complex*16 function fun1(g).  note the abscissas
c     are equally spaced (e=dexp(.1d0), er=1.0d0/e) in log-space.
      y1=abscis/bmax
c-----lagged convolution, outermost loop 1010
      do 1010 ilag=1,mb
        lag=lag+1
        istore=mb1-ilag
        if(lag.gt.0) y1=y1*e
        arg(istore)=abscis/y1
c-------related convolution, innermost loop 1000
        do 1000 jrel=1,nrel
c---------special case flag none=1 is set if fun1(g)=0 for all g in
c         filter fixed range (using weights 299-339).
          none=0
          itol=ntol
          zsum=zero
          cmax=zero
          y=y1
c---------begin right side convolution at weight 299 (m=return label)
          assign 20 to m
          i=299
          y=y*e
c---------call pseudo subroutine at 100 (return to 20 via m assigned)
          go to 100
   20     tmax(1)=dmax1(dabs(t(1)),tmax(1))
          tmax(2)=dmax1(dabs(t(2)),tmax(2))
          i=i+1
          y=y*e
c---------call pseudo subroutine at 100 (return to 20 via m assigned)
          if(i.le.339) go to 100
          if(tmax(1).eq.0.0d0.and.tmax(2).eq.0.0d0) none=1
c---------establish truncation criterion (cmax=dcmplx(tmax(1),tmax(2))
          cmax=tol*cmax
          assign 30 to m
c---------call pseudo subroutine at 100 (return to 30 via m assigned)
          go to 100
c---------check for filter truncation at right end
   30     if(dabs(t(1)).le.tmax(1).and.dabs(t(2)).le.tmax(2)) go to 50
          itol=ntol
   40     i=i+1
          y=y*e
c---------call pseudo subroutine at 100 (return to 30 via m assigned)
          if(i.le.801) go to 100
   50     itol=itol-1
          if(itol.gt.0.and.i.lt.801) go to 40
          itol=ntol
          y=y1
c---------continue with left side convolution at weight 298
          assign 60 to m
          i=298
c---------call pseudo subroutine at 100 (return to 60 via m assigned)
          go to 100
c---------check for filter truncation at left end
   60     if(dabs(t(1)).le.tmax(1).and.dabs(t(2)).le.tmax(2).and.
     *     none.eq.0) go to 80
          itol=ntol
   70     i=i-1
          y=y*er
c---------call pseudo subroutine at 100 (return to 60 via m assigned)
          if(i.gt.0) go to 100
   80     itol=itol-1
          if(itol.gt.0.and.i.gt.1) go to 70
c---------normalize zsum by arg(istore) to account for integration
c         range change, and store in zans(istore,jrel)
          zans(istore,jrel)=zsum/arg(istore)
c---------skip over pseudo subroutine to end of do 1000 innermost loop
          go to 1000
c
c=======================================================================
c=====store/retrieve pseudo subroutine for related/lagged convolution.
c     the internal (pseudo) subroutine entry is label 100, and returns
c     to the label assigned to m.  this calling mechanism could occur
c     a maximum of 801*mb*nrel times, where parameters mb>0 and nrel>0
c     can be arbitrarily large.  if a more-structured standard fortran
c     subroutine call was used, then the usual compiler linkage
c     convention could generate a maximum of 801*mb*nrel machine-
c     language instructions for register saves/restores and other
c     memory references.  for most compilers, timing tests reveal that
c     the pseudo-call method used here generated faster machine code
c     than with using external subroutine calls (e.g., call linkage
c     versus pseudo-call ratio was 2.6:1 on a vax-11/780 using
c     mb=50,nrel=61, and nofun1=199).
c
c=========pseudo-call entry point at 100 (returns via go to m below)
  100     look=i+lag
          iq=look/802
          ir=mod(look,802)
          if(ir.eq.0) ir=1
          iroll=iq*801
          if(key(ir).le.iroll) go to 150
c=========use existing saved functional values in zwork(ir,jrel)
 110      if (nord(jrel) .lt. 0.0) then
            goto  130
          elseif (nord(jrel) .gt. 0.0) then
            goto 130
          else
            goto 120
          endif

  120     c=zwork(ir,jrel)*wt0(i)
          go to 140
  130     c=zwork(ir,jrel)*wt1(i)
  140     zsum=zsum+c
c=========return convolution control via assigned m value, and with
c         the last convolution product c=dcmplx(t(1),t(2)).
          go to m,(20,30,60)
c=========compute external fun1 directly only when necessary
  150     key(ir)=iroll+ir
          g=y
c=========get direct evaluation of 1st kernel fun1(g), store in zwork.
c -af      if (ir.eq.299 .or. ir.eq.406 .or. ir.eq.407) then
c -af	    iwsave=iwaf
c -af	    iwaf=1.0e+06
c -af	    write (*,*) ' special ir, get debug output from kernels'
c -af	  end if
        zwork(ir,1)=fun1(g)
          nofun1=nofun1+1
c -af	  write (*,*) 
c -af     $    '    hy2fht:  zwork(',ir,',1)=fun1(',g,')=',zwork(ir,1) 
c>>>>>>>>>check for possible abort only if jpath=-1 (see subr. check2)
          if(jpath.eq.-1) then
            call check2(zwork,ir)
            if(jpath.eq.-2) return
          endif
c=========get all related zwork(ir,j),j=2,nrel for this ir if nrel>1
        if(nrel.gt.1) then
          call rel1(g,zwork(ir,1),fwork,nrel)
          do 160 j=2,nrel
            zwork(ir,j)=fwork(j)
160       continue
c -af	  write (*,*) '      zwork(',ir,',2)=',zwork(ir,2) 
c -af	  write (*,*) '      zwork(',ir,',',nrel,')=',zwork(ir,nrel) 
c -af      if (ir.eq.299 .or. ir.eq.406 .or. ir.eq.407) iwaf=iwsave
        endif
          go to 110
c=========end of pseudo subroutine (entry 100, return go to m above)
c=======================================================================
c
c-------end loop 1000 (get remaining related convolutions for this arg)
 1000   continue
c-----end loop 1010 (get remaining lagged convolutions for next arg)
 1010 continue
c-----exit with zans(mb,nrel),arg(mb) completed with minimal fun1 calls
      return
 9999 ierr=1
      return
      end
      subroutine check2(zwork,ir)
c-----------------------------------------------------------------------
c check2 is called by subroutine hankl2 only if jpath=-1 in common
c and to check if a highly oscillatory kernel function is being used.
c if the kernel is oscillating greater than the filter nyquist frequency
c then jpath=-2 is set and subroutine check2 immediately returns.  this
c in turn will force subroutine hy2fht to use jpath=-2 and thus the
c quadrature subroutine besax2 will then be invoked (i.e., subroutine
c hankl2 will be aborted dynamically for this type of condition, but
c only when ipath=-1 in call hybfht).
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 zwork(801,1)
      common/contrl/jpath,junk
      if(ir.eq.299.or.ir.eq.300) then
        return
      else if(ir.gt.300) then
        ir2=ir-2
        ir1=ir-1
      else if(ir.lt.299) then
        ir2=ir+2
        ir1=ir+1
      endif
        r2=dreal(zwork(ir2,1))
        if(dabs(r2).lt.1.d-10) go to 10
        r1=dreal(zwork(ir1,1))
        if(dabs(r1).lt.1.d-10) go to 10
        r=dreal(zwork(ir,1))
        if(dabs(r).lt.1.d-10) go to 10
        if(dsign(1.d0,r2).ne.dsign(1.d0,r1).and.
     1    dsign(1.d0,r1).ne.dsign(1.d0,r)) then
          jpath=-2
          return
        endif
10      r2=dimag(zwork(ir2,1))
        if(dabs(r2).lt.1.d-10) return
        r1=dimag(zwork(ir1,1))
        if(dabs(r1).lt.1.d-10) return
        r=dimag(zwork(ir,1))
        if(dabs(r).lt.1.d-10) return
        if(dsign(1.d0,r2).ne.dsign(1.d0,r1).and.
     1    dsign(1.d0,r1).ne.dsign(1.d0,r)) then
          jpath=-2
          return
        endif
      return
      end
