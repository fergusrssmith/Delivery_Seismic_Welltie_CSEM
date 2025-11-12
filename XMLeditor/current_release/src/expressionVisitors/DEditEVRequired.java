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
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import java.util.Vector;
import model.DEditTagInfo;
import model.DEditDebug;

/**
  * Retrieves, as a Vector of DEditTagInfo's, all required children of the 
  * visited Expression.  The reason this cannot be done with DEditEVOptions is 
  * that DEditEVRequired makes a choice on onChoice() (epsilon -- no child 
  * -- if possible, otherwise the first choice).  But I believe it may mess
  * up sometimes?  If (A | B), and they have B, it might return A.  However,
  * given that everything is checked when they are adding nodes, it won't 
  * break then, and it will only pick A or B when there is nothing there 
  * (during the initial default XML generation).
  */
public final class DEditEVRequired implements RELAXExpressionVisitor {
  
    private NodeList nl = null;
   
    public void setNodeList(NodeList nl) {
	this.nl = nl;
    }

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
	return null;
    } 
	
    public Object onChoice(ChoiceExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	boolean hasEpsilon = false;
	Vector subRet = null;
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    subRet = (Vector)(children[i].visit(this));	    
	    if (subRet == null || subRet.size() == 0) {

	    } else {
		int oldSize = subRet.size();
		if (oldSize == 1) {
		    DEditTagInfo ti = (DEditTagInfo)(subRet.get(0));
		    if (ti.getMinOccurs() == DEditTagInfo.EPSILON) {
			return null;
		    } 
		}
		if (nl != null) {
		    int childIndex = 0;
		    for (int j = 0; j < subRet.size(); j++) {
			DEditTagInfo ti = (DEditTagInfo)(subRet.get(j));
			int minOccurs = ti.getMinOccurs();
			if (minOccurs != DEditTagInfo.EPSILON) {
			    for (int k = childIndex; k < nl.getLength(); k++) {
				ElementExp eexp=(ElementExp)ti.getExpression();
				String name = eexp.getNameClass().toString();
				Node n = nl.item(k);
				if (n.getNodeType() == Node.ELEMENT_NODE 
				    && name.equals(n.getNodeName())) {
				    minOccurs--;
				    childIndex = k+1;
				    if (minOccurs <= 0) {
					subRet.remove(j);
					j--;
				    }
				}
			    }
			}
		    }
		}
		if (subRet.size() == 0) {
		    return null;
		} else if (subRet.size() < oldSize) {
		    return subRet; // has part of this one satisfied
		}
	    }
	}
	return subRet; // nobody has nothing, return last choice.
    } 
	
    public Object onConcur(ConcurExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	Vector ret = new Vector();
	boolean hasEpsilon = false;
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    Vector subRet = (Vector)(children[i].visit(this));	    
	    if (subRet == null) {

	    } else {
		if (subRet.size() == 1) {
		    DEditTagInfo ti = (DEditTagInfo)(subRet.get(0));
		    if (ti.getMinOccurs() == DEditTagInfo.EPSILON) {
			hasEpsilon = true;
		    } else {
			ret.addAll(subRet);
		    }
		} else {
		    ret.addAll(subRet);
		}
	    }
	    
	}
	return ret;
    } 
                   
    public Object onData(DataExp exp) {
	return null;
    }   
                   
    public Object onElement(ElementExp exp) {
	String myIndent = indent;
	indent += " ";
	Vector ret = new Vector();
	ret.add(0,new DEditTagInfo(exp,1,1));
	return ret;
    }
                   
    public Object onEpsilon() {
	Vector vec = new Vector();
	vec.add(new DEditTagInfo(null,DEditTagInfo.EPSILON,0));
	return vec;
    } 
                   
    public Object onInterleave(InterleaveExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	Vector ret = new Vector();
	boolean hasEpsilon = false;
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    Vector subRet = (Vector)(children[i].visit(this));	    
	    if (subRet == null) {

	    } else {
		if (subRet.size() == 1) {
		    DEditTagInfo ti = (DEditTagInfo)(subRet.get(0));
		    if (ti.getMinOccurs() == DEditTagInfo.EPSILON) {
			hasEpsilon = true;
		    } else {
			ret.addAll(subRet);
		    }
		} else {
		    ret.addAll(subRet);
		}
	    }
	    
	}
	return ret;
    }   
                   
    public Object onList(ListExp exp) {
	return exp.exp.visit(this);
    }  
                   
    public Object onMixed(MixedExp exp) {
	return exp.exp.visit(this);
    }   
                   
    public Object onNullSet() {
	return null;
    }  
                   
    public Object onOneOrMore(OneOrMoreExp exp) {
	String myIndent = indent;
	indent += " ";
	Vector ret = (Vector)(exp.exp.visit(this));
	if (ret == null || ret.size() == 0) {
	    return null;
	} 
	DEditDebug.println(this,2,"onOneOrMore");
	for (int i = 0; i < ret.size(); i++) {
	    DEditTagInfo ti = (DEditTagInfo)(ret.get(i));
	    DEditTagInfo retTi = new DEditTagInfo(ti.getExpression(),
						  ti.getMinOccurs(),
						  DEditTagInfo.UNLIMITED);
	    DEditDebug.println(this,2,retTi.toString());
	    ret.setElementAt(retTi,i);
	}
	DEditDebug.println(this,2,"end onOneOrMore");
	return ret;
    }   
                   
    public Object onOther(OtherExp exp) {
	return exp.exp.visit(this);
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
	    if (subRet == null) {

	    } else {
		if (subRet.size() == 1) {
		    DEditTagInfo ti = (DEditTagInfo)(subRet.get(0));
		    if (ti.getMinOccurs() == DEditTagInfo.EPSILON) {
			hasEpsilon = true;
		    } else {
			for (int j = 0; j < subRet.size(); j++) {
			    DEditTagInfo newTi = (DEditTagInfo)(subRet.get(j));
			    boolean hasNewTi = false;
			    for (int k = 0; k < ret.size(); k++) {
				DEditTagInfo oldTi=(DEditTagInfo)(ret.get(k));
				Expression oldExp = oldTi.getExpression();
				if (oldExp == newTi.getExpression()) {
				    hasNewTi = true;
				    int newMin = oldTi.getMinOccurs()
					+ newTi.getMinOccurs();
				    int newMax = DEditTagInfo.UNLIMITED;
				    if ((oldTi.getMaxOccurs() 
					 != DEditTagInfo.UNLIMITED)
					&& (newTi.getMaxOccurs() 
					    != DEditTagInfo.UNLIMITED)) {
					newMax = oldTi.getMaxOccurs()
					    + newTi.getMaxOccurs();
				    }
				    ret.setElementAt(new DEditTagInfo(oldExp,
								      newMin,
								      newMax),
						     k);
				    break;
				}
			    }
			    if (!hasNewTi) {
				ret.add(newTi);
			    }
			}
		    }
		} else {
		    for (int j = 0; j < subRet.size(); j++) {
			DEditTagInfo newTi = (DEditTagInfo)(subRet.get(j));
			boolean hasNewTi = false;
			for (int k = 0; k < ret.size(); k++) {
			    DEditTagInfo oldTi=(DEditTagInfo)(ret.get(k));
			    Expression oldExp = oldTi.getExpression();
			    if (oldExp == newTi.getExpression()) {
				hasNewTi = true;
				int newMin = oldTi.getMinOccurs()
				    + newTi.getMinOccurs();
				int newMax = DEditTagInfo.UNLIMITED;
				if ((oldTi.getMaxOccurs() 
				     != DEditTagInfo.UNLIMITED)
				    && (newTi.getMaxOccurs() 
					!= DEditTagInfo.UNLIMITED)) {
				    newMax = oldTi.getMaxOccurs()
					+ newTi.getMaxOccurs();
				}
				ret.setElementAt(new DEditTagInfo(oldExp,
								  newMin,
								  newMax),
						 k);
				break;
			    }
			}
			if (!hasNewTi) {
			    ret.add(newTi);
			}
		    }
		}
	    }
	
	}
	return ret;
    }   
                   
    public Object onValue(ValueExp exp)  {
	return null;
    }  
}
