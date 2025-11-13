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
!    
!    You should have received a copy of the GNU General Public License
!    along with LegacyScrippsDipole1D.  If not, see <http://www.gnu.org/licenses/>.


module fmodel
  !     include 'occamdim.inc'
  !     commons for seafloor controlled source (with ipen1,2 now reals, JGHACK) (and MT)
  !      real thickf, rho, freqs, thite, tangle
  !      integer nlayer, iff, nfreq
  !      real ipen1, ipen2
  !      common /fmodel/ thickf(0:npp),rho(0:npp),freqs(ndd),thite,tangle,ipen1(0:npp),ipen2(0:npp),iff(0:npp),nfreq,nlayer
  !      save /fmodel/
  !    
  !     Allocation sizes used in module:
  !       thickf(0:nlayer), rho(0:nlayer), freqs(nd)
  !       ipen1(0:nlayer), ipen2(0:nlayer), iff(1:np)

  real, pointer ::  thickf(:), rho(:), freqs(:), ipen1(:), ipen2(:)
  integer, pointer :: iff(:)

  real :: thite, tangle
  integer::  nlayer, nfreq

!$omp threadprivate(thickf, rho,freqs,ipen1,ipen2,thite,tangle,nlayer,nfreq)

end module fmodel

module mode
  !     c new commons to allow both marquardt and occam to use same routines
  character*10 rmode
!$omp threadprivate(rmode)
end module mode

subroutine set_module_pointers(np, nd, nlayer_in, nfreq_in, thickf_ptr, rho_ptr, freqs_ptr, ipen1_ptr,ipen2_ptr, iff_ptr)
  use fmodel
  integer :: np,nd, nlayer_in, nfreq_in

  
  real, target :: thickf_ptr(nlayer_in+1)
  real, target :: rho_ptr(nlayer_in+1), freqs_ptr(nd) 
  real, target :: ipen1_ptr(nlayer_in+1), ipen2_ptr(nlayer_in+1)
  integer, target :: iff_ptr(1:np)

  nlayer = nlayer_in
  nfreq = nfreq_in

  thickf => thickf_ptr(2:nlayer_in+1)
  rho => rho_ptr(2:nlayer_in+1)      
  freqs => freqs_ptr     
  ipen1 => ipen1_ptr(2:nlayer_in+1)     
  ipen2 => ipen2_ptr(2:nlayer_in+1)   
  iff  => iff_ptr          

end subroutine set_module_pointers


subroutine set_simple_module_pointers(np,nlayer_in, thickf_ptr, rho_ptr, iff_ptr)
  use fmodel
  integer :: np,nlayer_in
  
  real, target :: thickf_ptr(nlayer_in+1)
  real, target :: rho_ptr(nlayer_in+1)
  integer, target :: iff_ptr(1:np)

  nlayer = nlayer_in

  thickf => thickf_ptr(2:nlayer_in+1)
  rho => rho_ptr(2:nlayer_in+1)      
  iff  => iff_ptr          

end subroutine set_simple_module_pointers
