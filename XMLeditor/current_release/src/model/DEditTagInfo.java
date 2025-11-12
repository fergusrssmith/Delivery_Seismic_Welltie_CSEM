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
package model;
import com.sun.msv.grammar.Expression;
import com.sun.msv.grammar.NameClassAndExpression;

/**
  * Class which represents information about a tag or element (but not
  * its content).
  */
public class DEditTagInfo {
    
    /**
      * Special value for maximum number of occurences which means unlimited
      * occurances are allowed.
      */
    public static final int UNLIMITED = -1;
    /**
      * Special value used to say that this is an epsilon or 
      * non-existant tag.  It means that if it is the result of a 
      * query of TagInfo's on the children of a choice node, then 
      * null or no tag is a valid option.
      */
    public static final int EPSILON = -2;
    
    private Expression exp; // the Expression corresponding to this tag
    private int minOccurs; // minimum number of times this tag can occur
    private int maxOccurs; // maximum number of times this tag can occur

    /**
      * Creates a DEditTagInfo with given information
      * @param exp the Expression this tag has information about
      * @param min the minimum number of times this tag's Expression must occur
      * @param max the maximum number of times this tag's Expression can occur
      */
    public DEditTagInfo(Expression exp, int min, int max) {
	this.exp = exp;
	this.minOccurs = min;
	this.maxOccurs = max;
    }

    /**
      * Retruns the minimum number of occurences of this tag's Expression.
      * @return the minimum number of occurences of this tag's Expression
      */
    public int getMinOccurs() {
	return minOccurs;
    }

    /**
      * Retruns the maximum number of occurences of this tag's Expression.
      * @return the maximum number of occurences of this tag's Expression
      */
    public int getMaxOccurs() {
	return maxOccurs;
    }

    /**
      * Retruns the Expression associated with this DEditTagInfo.
      * @return the Expression associated with this DEditTagInfo
      */
    public Expression getExpression() {
	return exp;
    }

    /**
      * Overrides Object.toString() to give some useful information about 
      * the DEditTagInfo.
      * @return a String representation of the DEditTagInfo
      */
    public String toString() {
	return "min: " + minOccurs + " max: " + maxOccurs + " name: "
	    + ((NameClassAndExpression)exp).getNameClass().toString();
    }

}
