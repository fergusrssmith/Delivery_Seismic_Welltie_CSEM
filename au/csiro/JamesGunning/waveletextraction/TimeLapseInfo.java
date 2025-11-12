/**********************************************************************
DELIVERY - a program for stochastic seismic inversion
This program module Copyright (C) 2006  BHP Billiton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA or visit the link http://www.gnu.org/licenses/gpl.txt.

PLEASE NOTE: One or more BHP Billiton patents or patent applications may
cover this software or its use. BHP Billiton hereby grants a limited,
personal, irrevocable, perpetual, royalty-free license under such
patents but only for things you are allowed to do under the terms of
the General Public License.  All other rights are expressly reserved.

To contact BHP Billiton about this software you can e-mail
James.Gunning@csiro.au or visit http://www.petroleum.csiro.au to learn
more.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 



package au.csiro.JamesGunning.waveletextraction;
import au.csiro.JamesGunning.lib.*;

/**
   Basic entries, constants etc for TimeLapse signature estimation
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/11/2017
*/

public class TimeLapseInfo
{
  String base_contacts_filename;
  String monitor_contacts_filename;
  boolean base_timeunitssecs = true;
  boolean monitor_timeunitssecs = true;
  // sand props
  double vp_s,vs_s,phi_s;
  // shale props
  double vp_m,vs_m,rho_m;
  // moving fluid (hydrocarbon..) props
  double vp_h,rho_h;
  // displaced fluid props (brine)
  double vp_b,rho_b;
  // sand grain props
  double rho_g,K_g;
}
