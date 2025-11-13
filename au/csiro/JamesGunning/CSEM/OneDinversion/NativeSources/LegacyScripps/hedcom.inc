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


c                   hedcom.inc
      implicit double precision (a-h,o-z)
c   parameters and common blocks for programs associated with hedsfi
c   data stmts in file heddat.inc to satisfy standard statement ordering.
c
c   ===================================================================
c   =====   edit npl to change max no. layers     =====================
c   =====   edit maxrng to change max no. ranges  =====================
c   ==  edit expos,exneg to change under/overflow checks in kernels ===
c   ===================================================================
c
c   npl=no layers (finite layers only, not including sea, earth halfspaces)
c   maxrng=no ranges to do at once (fewer lagged hankel-xfms may be performed)
c JGHACK npl=100 HERE!!!!!!!!!
      parameter (npl=1100,maxrng=500)
      parameter (expos=30.0,exneg=300.0)
c
c   ===================================================================
c   =====     no need to edit below here !          ===================
c   ===================================================================
c
c  logical units for i/o, write-out level
      common /inoutc/ni,no,iw
      save /inoutc/
c  for debug-output; program unit name and quit-message
      character namepu*8,msg*72
c
c   npl1 is no. layers including the bottom halfspace (but not ocean)
c   nppar=max no of model parameters (sigma's and thick's)
      parameter (npl1=npl+1,nppar=npl1+npl)
c   npdeh=max no of field & derivative values potentially returned by hedsfi
      parameter (npdeh=2*npl1)
c
c  hankel-xfm workspace dimnsn: maxrel=no of related transforms
c  (maxrng is used for no of lagged transforms maxlag)
      parameter (maxrel=2*npdeh)
c
c  the variables in common-blocks hed and passc 
c    are passed in common so that workspace arrays and
c    variables that remain fixed need not be passed as
c    arguments into the kernels
c
c--common /hed/ contains model parameters
c   some of these variables (sigma, thick, etc.)
c   are filled in by routine hedfix,
c   which must be called once before hedsfi is called.
c   others (phi,ifdiff,..) are filled in by hedsfi or dehhed.
      common /hed/sigma(0:npl1),thick(npl)
     $           ,h,zp,zm,sl,sc,phi,nl,nl1
     $           ,ifdiff,ianatm 
      save /hed/
c
c-- common /passc/ is used internally by dehhed and its subprograms
      complex*16 ciommu,gamma2,beta0
      common /passc/ ciommu,sigscl,dsgscl
     $     ,gamma2,beta0,ommu,jl,jmax 
      save /passc/
c
c--common /fdc/ conpains refl coeffs with parameter-derivs,from fdrl and fun1
      complex*16 zfdpm,zfdtm,dedrpm,dedrtm
      common /fdc/zfdpm(0:nppar),zfdtm(0:nppar),dedrpm,dedrtm 
      save /fdc/
c
c--weights from block data zwtdat for use by digital filter hankel transform routine zhanks
c  zhanks does not include this file
c-ahed      common /zwtcom/ wt0(801),wt1(801) 
c-ahed      save /zwtcom/
