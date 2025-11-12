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
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
  * Returns, as a Vector of DEditTagInfos, all possible children of
  * the visited Expression.
  */
public final class DEditEVOptions implements RELAXExpressionVisitor {
  
    private String indent = "";
    private NodeList nl = null;
   
    public void setNodeList(NodeList nl) {
	this.nl = nl;
    }

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
	Vector ret = new Vector();
	boolean satisfied = false;
	Vector spare = null;
	boolean hasEpsilon = false;
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    Vector subRet = (Vector)(children[i].visit(this));	    
	    if (subRet == null || subRet.size() == 0) {

	    } else {
		if (subRet.size() == 1 
		    && ((DEditTagInfo)subRet.get(0)).getMinOccurs() 
		    == DEditTagInfo.EPSILON) {
		    hasEpsilon = true;
		} else {
		    if (nl != null && !satisfied && isSatisfied(subRet)) {
			spare = subRet;
			satisfied = true;
		    }
		    for (int j = 0; j < subRet.size(); j++) {
			DEditTagInfo ti = (DEditTagInfo)(subRet.get(j));
			boolean addOn = true;
			for (int k = 0; k < ret.size(); k++) {
			    DEditTagInfo other = (DEditTagInfo)(ret.get(k));
			    Expression current = ti.getExpression();
			    if (current == other.getExpression()) {
				addOn = false;
				int newMax = DEditTagInfo.UNLIMITED;
				if ((ti.getMaxOccurs() 
				     != DEditTagInfo.UNLIMITED)
				    && (other.getMaxOccurs() 
					!= DEditTagInfo.UNLIMITED)) {
				    newMax = ti.getMaxOccurs() 
					+ other.getMaxOccurs();
				}
				int newMin = ti.getMinOccurs() 
				    + other.getMinOccurs();
				ret.setElementAt(new DEditTagInfo(current,
								  newMin,
								  newMax),k);
			    }
			}
			if (addOn) {
			    ret.add(ti);
			}
		    }
		}
	    }
	}
	if (satisfied) {
	    for (int i = 0; i < ret.size(); i++) {
		DEditTagInfo old = (DEditTagInfo)ret.get(i);
		boolean spareThisTi = false;
		for (int j = 0; j < spare.size(); j++) {
		    DEditTagInfo spareTi = (DEditTagInfo)(spare.get(j));
		    if (old.getExpression() == spareTi.getExpression()) {
			spareThisTi = true;
			break;
		    }
		}
		if (!spareThisTi) {
		    ret.setElementAt(new DEditTagInfo(old.getExpression(),
						      old.getMinOccurs(),0),i);
		}
	    }
	}
	if (hasEpsilon) {
	    for (int i = 0; i < ret.size(); i++) {
		DEditTagInfo old = (DEditTagInfo)ret.get(i);
		ret.setElementAt(new DEditTagInfo(old.getExpression(),
						  0,old.getMaxOccurs()),i);
	    }
	}
	return ret;
    } 

    private boolean isSatisfied(Vector tagInfos) {
	int childIndex = 0;
	int debugLevel = 0;
	DEditTagInfo first = (DEditTagInfo)(tagInfos.get(0));
	ElementExp exp = (ElementExp)first.getExpression();
	if (exp.getNameClass().toString().equals("bhpread")) {
	    debugLevel = 4;
	}
	for (int i = 0; i < tagInfos.size(); i++) {
	    DEditTagInfo ti = (DEditTagInfo)(tagInfos.get(i));
	    DEditDebug.println(this,debugLevel,ti.toString());
	    int minOccurs = ti.getMinOccurs();
	    ElementExp eexp=(ElementExp)ti.getExpression();
	    String name = eexp.getNameClass().toString();
	    for (int j = childIndex; j < nl.getLength(); j++) {
		Node n = nl.item(j);
		if (n.getNodeType() == Node.ELEMENT_NODE) {
		    String nodeName = n.getNodeName();
		    DEditDebug.println(this,debugLevel,nodeName);
		    if (name.equals(nodeName)) {
			minOccurs = Math.max(0,minOccurs-1);
			childIndex = j+1;
		    }
		}
	    }
	    if (minOccurs > 0) {
		DEditDebug.println(this,debugLevel,"we want more>>>");
		return false;
	    }
	}
	DEditDebug.println(this,debugLevel,"satisfied");
	return true;
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
	if (hasEpsilon) {
	    for (int i = 0; i < ret.size(); i++) {
		DEditTagInfo ti = (DEditTagInfo)(ret.get(i));
		ret.setElementAt(new DEditTagInfo(ti.getExpression(),0,
						  ti.getMaxOccurs()),i);
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
	DEditDebug.println(this,4,"launching Epsilon");
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
	if (hasEpsilon) {
	    for (int i = 0; i < ret.size(); i++) {
		DEditTagInfo ti = (DEditTagInfo)(ret.get(i));
		ret.setElementAt(new DEditTagInfo(ti.getExpression(),0,
						  ti.getMaxOccurs()),i);
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
	DEditDebug.println(this,4,"onOneOrMore");
	for (int i = 0; i < ret.size(); i++) {
	    DEditTagInfo ti = (DEditTagInfo)(ret.get(i));
	    DEditTagInfo retTi = new DEditTagInfo(ti.getExpression(),
						  ti.getMinOccurs(),
						  DEditTagInfo.UNLIMITED);
	    DEditDebug.println(this,4,retTi.toString());
	    ret.setElementAt(retTi,i);
	}
	DEditDebug.println(this,4,"end onOneOrMore");
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
