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

import model.DEditDebug;

/**
  * Class which represents a scatter parameter that covers a range 
  * (e.g. 2-26[4]).
  */
public class RangeScatterParam extends AScatterParam {
    /**
      * minimum value in the range to scatter on
      */
    public int min;
    /**
      * maximum value in the range to scatter on
      */
    public int max;
    /**
      * increment at which to proceed through the range to scatter on
      */
    public int increment;
    /**
      * minimum number of values covered per thread.  -1 indicates that 
      * the range can not be compressed and one and only one value must be 
      * used per thread.
      */
    public int min_number_per_thread;
    
    /**
      * Constructs a ranged scatter parameter with given min and max values 
      * and increment, plus given name and min number of values per thread
      * @param name the name of the scatter parameter
      * @param min minimum value in the range to scatter on
      * @param max maximum value in the range to scatter on
      * @param increment increment at which to proceed through the range 
      * to scatter on
      * @param min_number_per_thread minimum number of values covered per 
      * thread.  -1 indicates that the range can not be compressed and one 
      * and only one value must be used per thread.
      */
    public RangeScatterParam(String name, int min, int max, int increment, 
			     int min_number_per_thread) {
	this.name = name;
	this.min = min;
	int temp = max - min;
	this.max = max - (temp % increment);
	this.increment = increment;
	this.min_number_per_thread = min_number_per_thread;
	DEditDebug.println(this,6,this.name + " " + this.min + " " 
			   + this.max + " " + this.increment + " " 
			   + this.min_number_per_thread);
    }    

    /**
      * Returns whether or not this scatter parameter should be resolved after 
      * the first round of resolving or with the first round.  If the parameter
      * is uncompressable (minimum number of values per thread == -1) then 
      * it should be resolved right away
      */
    public boolean resolveLater() {
	return (min_number_per_thread != -1);
    }

    /**
      * Resolves this RangeScatterParam into a ResolvedScatterParam which can 
      * be used to divide up values for a scatter.  If it is not compressable, 
      * it just calculates the values and puts them in an array.
      * If it is compressable, it expands it as much as possible, but it can 
      * stay as compressed as even not splitting at all
      * @param currThreads the current number of threads running
      * @param maxThreads the number of threads not to exceed (only 
      * guaranteed not to exceed if parameter is compressable)
      * @return a ResolvedScatterParam generated from this RangeScatterParam
      */
    public ResolvedScatterParam resolve(int currThreads, int maxThreads) 
	throws MalformedScriptException {
	if (increment == 0) {
	    throw new MalformedScriptException("Increment cannot equal 0 for "
					       +"ranged scatter parameter.");
	}
	// every thread gets one and only one value
	if (min_number_per_thread == -1) {
	    int rangeSize = (max+1-min);
	    int newMulFact = rangeSize / increment;
	    if (rangeSize % increment != 0) {
		newMulFact++;
	    }
	    String[] values = new String[newMulFact];
	    for (int j = 0; j < newMulFact; j++) {
		values[j] = String.valueOf(min + increment*j);
	    }
	    return new ResolvedScatterParam(name,values,newMulFact);
	} else if (min_number_per_thread <= 0) {
	    throw new MalformedScriptException("Illegal minimum number per "
					       + "thread" 
					       + min_number_per_thread 
					       + "for scatter parameter "
					       + name + ".");
	} else {
	    // each value is a range
	    int rangeSize = (max+1 - min);
	    int denom = (increment * min_number_per_thread);
	    int numScripts = rangeSize / denom;
	    if (rangeSize % denom != 0) {
		numScripts++;
	    }
	    DEditDebug.println(this,6,"> " + numScripts + " " + rangeSize
			       + "/" + denom);
	    int orig_num = numScripts;
	    if (maxThreads != -1) {
		while (numScripts*currThreads > maxThreads) {    
		    DEditDebug.println(this,6,"> " + numScripts);
		    min_number_per_thread++;
		    denom = (increment * min_number_per_thread);
		    numScripts = rangeSize / denom;
		    if (rangeSize % denom != 0) {
			numScripts++;
		    }
		    if (numScripts == 1) {
			break;
		    }
		}
	    }
	    String[] values = new String[numScripts];
	    for (int j = 0; j < numScripts; j++) {
		int total = (max+1 - min)/(increment);
		int thisOne = j*min_number_per_thread*increment;
		int start = min + thisOne;
		DEditDebug.println(this,3,name + " " + j + "_start " + min 
				   + " " + thisOne + " ");
		int finish = start + min_number_per_thread*increment - 1;
		finish = Math.min(finish,max);
		int temp = finish - start;
		finish -= (temp % increment);
		values[j] = start + "-" + finish 
		    + "[" + increment + "]";
	    }
	    ResolvedScatterParam ret =  new ResolvedScatterParam(name,values,
								 numScripts);
	    if (numScripts != orig_num) {
		ret.addWarning("Split on " + name + " reduced from " + orig_num
			       + " threads to " + numScripts + " thread(s).");
		DEditDebug.println(this,2,"add warning");
	    }
	    return ret;
	}
    }
}
