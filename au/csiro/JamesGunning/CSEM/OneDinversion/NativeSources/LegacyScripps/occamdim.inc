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

! OCCAM 2.0: Steven Constable IGPP/SIO La Jolla CA 92093-0225
! Commons Revision 2.00, 13 May 1992
!-----------------------------------------------------------------------
	integer npp, ndd, nnp
!   old (SC value) npp=100 
	parameter (npp = 1000, ndd = 500, nnp = 2*npp)
!   CSIRO (JG) Additions
	logical use_range_lumping
	parameter (use_range_lumping=.true.)
!
!   npp = maximum number of parameters used
!   ndd = maximum number of data used
!   nnp = number of rows in penalty matrix.  
!   For 1D, nnp = npp or greater
!   For 2D, nnp = 2*npp or greater
