/*************************************************************************
DELIVERY - a program for stochastic seismic inversion 
This program module Copyright (C) 2006, BHP Billiton and is licensed
under the following BSD License.

Copyright (C) 2006, BHP Billiton All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of BHP Billiton Petroleum nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************************/
//------------------------- EOCOPYRIGHT --------------------------------- 


package au.csiro.JamesGunning.lib;

/**
   Delivery modules' licence text
*/
public class DeliveryLicense
{ 
  private static String trailer = "Copyright (C) 2006 BHP Billiton%n"
	+"DELIVERY comes with ABSOLUTELY NO WARRANTY.%n"
	+"This is free software, and you are welcome to redistribute it under%n"
	+"certain conditions. Portions of DELIVERY are licensed under the GNU%n"
	+"General Public License Version 2 (see http://www.gnu.org/licenses/gpl.txt).%n"
	+"Other portions are licensed under the BSD License (see%n"
	+"http://www.opensource.org/licenses/bsd-license.php).  Please refer to%n"
	+"the headers of each source code file for more information.%n";// Licence text 

  // Licence text 
  public static String license_text = "DELIVERY - a program for stochastic seismic inversion%n"+trailer;
  // For CSEM  
  public static String csem_license_text = "DELIVERYCSEM - a program for Bayesian CSEM inversion%n"+trailer;
}
