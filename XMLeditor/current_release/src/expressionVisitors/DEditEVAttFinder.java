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
package expressionVisitors;

import com.sun.msv.grammar.relax.*;
import com.sun.msv.grammar.*;
import java.util.Vector;
import model.DEditTagInfo;
import model.DEditDebug;

/**
  * Class which finds and returns as a Vector of DEditTagInfos, all
  * attributes for the passed in Expression.
  */
public final class DEditEVAttFinder implements RELAXExpressionVisitor {
  
    private String indent = "";
    
    public Object onAttPool(AttPoolClause exp) {
	return null;
    }
	
    public Object onElementRules(ElementRules exp) {
	return null;
    }
	
    public Object onHedgeRules(HedgeRules exp) {
	return null;
    }
	
    public Object onTag(TagClause exp) {
	return null;
    }

    public Object onAnyString() {
	return null;
    }
                   
    public Object onAttribute(AttributeExp exp) {
	String myIndent = indent;
	indent += " ";
	Vector ret = new Vector();
	ret.add(0,new DEditTagInfo(exp,1,1));
	return ret;
    } 
	
    public Object onChoice(ChoiceExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	Vector ret = new Vector();
	boolean hasEpsilon = false;
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    Vector subRet = (Vector)(children[i].visit(this));	    
	    if (subRet == null || subRet.size() == 0) {

	    } else if (subRet.size() > 1) {
		ret.addAll(subRet);
	    } else {
		DEditTagInfo ti = (DEditTagInfo)(subRet.get(0));
		if (ti.getMinOccurs() == DEditTagInfo.EPSILON) {
		    hasEpsilon = true;
		} else {
		    ret.addAll(subRet);
		}
	    }
	}
	if (hasEpsilon) {
	    for (int i = 0; i < ret.size(); i++) {
		DEditTagInfo ti = (DEditTagInfo)(ret.get(i));
		ret.setElementAt(new DEditTagInfo(ti.getExpression(),
						  0,ti.getMaxOccurs()),i);
	    }
	}
	return ret;
    } 
	
    public Object onConcur(ConcurExp p) {
	return null;
    } 
                   
    public Object onData(DataExp exp) {
	return null;
    }   
                   
    public Object onElement(ElementExp exp) {
	return null;
    }
                   
    public Object onEpsilon() {
	Vector vec = new Vector();
	vec.add(new DEditTagInfo(null,DEditTagInfo.EPSILON,0));
	return vec;
    } 
                   
    public Object onInterleave(InterleaveExp p) {
	return null;
    }   
                   
    public Object onList(ListExp exp) {
	return null;
    }  
                   
    public Object onMixed(MixedExp exp) {
	return null;
    }   
                   
    public Object onNullSet() {
	return null;
    }  
                   
    public Object onOneOrMore(OneOrMoreExp exp) {
	String myIndent = indent;
	indent += " ";
	Vector ret = (Vector)(exp.exp.visit(this));
	if (ret.size() == 1) {
	    DEditTagInfo ti = (DEditTagInfo)(ret.get(0));
	    ret.setElementAt(new DEditTagInfo(ti.getExpression(),1,
					      DEditTagInfo.UNLIMITED),0);
	} else if (ret.size() == 0) {
	    return ret;
	} else {
	    DEditDebug.println(this,4,"Houston we have a problem.");
	}
	return ret;
    }   
                   
    public Object onOther(OtherExp exp) {
	return null;
    }   
                   
    public Object onRef(ReferenceExp exp) {
	String myIndent = indent;
	indent += " ";
	return exp.exp.visit(this);
    } 
                   
    public Object onSequence(SequenceExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	Vector ret = new Vector();
	boolean hasEpsilon = false;
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    Vector subRet = (Vector)(children[i].visit(this));	    
	    if (subRet == null || subRet.size() == 0) {
		
	    } else if (subRet.size() > 1) {
		ret.addAll(subRet);
	    } else {
		DEditTagInfo ti = (DEditTagInfo)(subRet.get(0));
		if (ti.getMinOccurs() == DEditTagInfo.EPSILON) {
		    hasEpsilon = true;
		} else {
		    ret.addAll(subRet);
		}
	    }
	}
	if (hasEpsilon) {
	    for (int i = 0; i < ret.size(); i++) {
		DEditTagInfo ti = (DEditTagInfo)(ret.get(i));
		ret.setElementAt(new DEditTagInfo(ti.getExpression(),0,
						  ti.getMaxOccurs()),i);
	    }
	}
	return ret;	
    }   
                   
    public Object onValue(ValueExp exp)  {
	return null;
    }  
}
