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
package script;

import java.util.Vector;

/**
  * Abstract class which represents a scatter parameter.  It's main 
  * functionality is to handle the warnings and enforce the resolve and 
  * resolveLater methods be defined in any non-abstract subclass. 
  */
public abstract class AScatterParam {
    /**
      * The name of this scatter parameter
      */
    public String name;
    protected Vector warnings = new Vector(); // a String vector of warning
    // messages

    /**
      * turns any AScatterParam into a ResolvedScatterParam -- the type 
      * that Scatters hold and the type that can be used to get the ith 
      * value for the scatter.
      * @param currThreads the current number of threads running
      * @param maxThreads the maximum threads generateable without a warning
      * (also parameters which can be ranged will compress to fit under
      * maxThreads)
      * @return a ResolvedScatterParam which can then be used to get the 
      * needed values for each thread in the scatter
      */
    public abstract ResolvedScatterParam resolve(int currThreads,
						 int maxThreads) 
	throws MalformedScriptException;

    /**
      * Adds a warning about this scatter parameter to the list of warnings.
      * @param warning the warning to add
      */
    public final void addWarning(String warning) {
	warnings.add(warning);
    }

    /**
      * Returns as a vector the list of warnings accrued for this scatter 
      * parameter so far.
      * @return the list of warnings for this scatter parameter
      */
    public final Vector getWarnings() {
	return warnings;
    }

    /**
      * Returns true if this scatter parameter should be resolved after the 
      * initial round of resolving.  This basically means it can be compressed 
      * so it is wise to do the uncompressable ones first in an attempt to 
      * stay under max_threads threads.
      * @return true if this scatter parameter should be resolved after the 
      * initial round of resolving
      */
    public abstract boolean resolveLater();

}
