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
  * Represents a scatter parameter which has been resolved such that 
  * it has x number of values and the scripts just have to get the proper
  * index into the array of values and ask the ResolvedScatterParam for the
  * value.
  */
public class ResolvedScatterParam extends AScatterParam {
    
    private String[] values; // the array of values to place in scripts
    private int mulFact; // the expansion factor on the number of threads

    /**
      * Creates a "ready-to-go" AScatterParam with the right name,
      * values, and mulFact (== values.length I believe).  Unless you 
      * somehow know all values already, create a different type of scatter
      * param and call .resolve() to get a ResolvedScatterParam instead of 
      * using this constructor.
      * @param name the name of the scatter parameter
      * @param values the different values of the scatter parameter
      * @param mulFact the multiplication factor to the total number of 
      * "threads" for this scatter parameter
      */
    public ResolvedScatterParam(String name,String[] values, int mulFact) {
	this.name = name;
	this.values = values;
	this.mulFact = mulFact;
	String valuesArr = "";
	for (int i = 0; i < values.length; i++) {
	    valuesArr += values[i]+",";
	}
	DEditDebug.println(this,6,"name: " + name + ", mulFact: " + mulFact
			   + ", values: " + valuesArr);
    }
    
    /**
      * Gets the factor by which this scatter parameter increases the number 
      * of threads, must be >= 1.
      * @return the factor by which this scatter paraemter increases the 
      * number of threads
      */
    public int getMulFact() {
	return mulFact;
    }

    /**
      * Gets the value at a certain index into the resolved values of this
      * ResolvedScatterParameter.
      * @param i the index into the value array at which to retrieve the value
      * @param commandOpt the source of the request for a value, used to check 
      * if it is a legal place for a range (if the value is a range)
      */
    public String getValue(int i, String commandOpt) 
	throws MalformedScriptException {
	if (values == null || i < 0 || i >= values.length) {
	    throw new MalformedScriptException("Internal indexing problem: "
					       + "null values or outOfRange");
	}
	String ret = values[i];
	// sometimes you may want "[" in the string, in which case add 
	// && !commandOpt.equals(appropriate_string_for_command)
	// to if clause
	if (ret.indexOf("[") != -1
	    && !commandOpt.equals("bhpread=key->range")
	    && !commandOpt.equals("unix_command=text")) {
	    throw new MalformedScriptException("Ranged parameter not allowed "
					       + "for " + commandOpt + ".");
	}
	return ret;
    }
    
    /**
      * Returns true if this scatter parameter should be resolved after the 
      * initial round of resolving instead of with the first round.
      * returns false in this case since it is already resolved.
      * @return false
      */
    public boolean resolveLater() {
	return false;
    }

    /**
      * Resolves the scatter parameter, but since a ResolvedScatterParam is 
      * already resolved, just returns "this".
      * @param currThreads not used
      * @param maxThreads not used
      * @return "this"
      */
    public ResolvedScatterParam resolve(int currThreads, int maxThreads) {
	return this;
    }
}
