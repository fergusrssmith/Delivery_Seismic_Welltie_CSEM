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
import model.DEditDebug;

/**
  * Class which finds the Expression corresponding to a passed in chain 
  * of Nodes and returns it
  */
public final class DEditEVExpression 
    implements RELAXExpressionVisitorExpression {
  
    private String indent = "";
    private Node[] nodes = null;
    private int nodeIndex = 1;  // want to start with 1 instead of zero,
    // but there might be a cleaner, more obvious way of doing it

    /**
      * Sets the node chain to use by this DEditEVExpression.  Must be 
      * done before every new traversal.
      * @param nodes the node chain to use to find an Expression
      */
    public void setNodes(Node[] nodes) {
	this.nodes = nodes;
	this.nodeIndex = 1;
    }

    public Expression onAttPool(AttPoolClause exp) { return null; }	
    public Expression onElementRules(ElementRules exp) { return null; }
    public Expression onHedgeRules(HedgeRules exp) { return null; }	
    public Expression onTag(TagClause exp) { return null; }
    public Expression onAnyString() { return null; }                   
    public Expression onAttribute(AttributeExp exp) { return null; } 
	
    public Expression onChoice(ChoiceExp exp) {
	int curIndex = nodeIndex;
	Expression ret = null;
	Vector possibilities = new Vector();
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    nodeIndex = curIndex;
	    ret = children[i].visit(this);
	    if (ret != null) {
		possibilities.add(ret);
	    }
	}
	if (possibilities.size() == 0) {
	    return null;
	} else if (possibilities.size() == 1) {
	    return (Expression)possibilities.get(0);
	} else {
	    return resolveMatch(false,possibilities,curIndex);
	}
    } 
	
    public Expression onConcur(ConcurExp exp) {
	int curIndex = nodeIndex;
	Expression ret = null;
	Vector possibilities = new Vector();
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    nodeIndex = curIndex;
	    ret = children[i].visit(this);
	    if (ret != null) {
		possibilities.add(ret);
	    }
	}
	if (possibilities.size() == 0) {
	    return null;
	} else if (possibilities.size() == 1) {
	    return (Expression)possibilities.get(0);
	} else {
	    return resolveMatch(true,possibilities,curIndex);
	}   
    } 
                   
    public Expression onData(DataExp exp) { return null; }   
                   
    public Expression onElement(ElementExp exp) {
	NameClass nmc = exp.getNameClass();
	String ncString = nmc.toString();
	String nodeString = nodes[nodeIndex].getNodeName();
	if (ncString == null || nodeString == null) {
	    return null;
	} else if (ncString.equals(nodeString)) {
	    nodeIndex++;
	    if (nodeIndex == nodes.length) {
		return exp; // found it
	    } else {
		return exp.getContentModel().visit(this);
	    }
	} else {
	    return null;
	} 
    }
                   
    public Expression onEpsilon() { return null; } 
                   
    public Expression onInterleave(InterleaveExp exp) { 
	int curIndex = nodeIndex;
	Expression ret = null;
	Vector possibilities = new Vector();
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    nodeIndex = curIndex;
	    ret = children[i].visit(this);
	    if (ret != null) {
		possibilities.add(ret);
	    }
	}
	if (possibilities.size() == 0) {
	    return null;
	} else if (possibilities.size() == 1) {
	    return (Expression)possibilities.get(0);
	} else {
	    return resolveMatch(true,possibilities,curIndex);
	}
    }   
                   
    public Expression onList(ListExp exp) {
	return exp.exp.visit(this);
    }  
                   
    public Expression onMixed(MixedExp exp) { 
	return exp.exp.visit(this);
    }   
                   
    public Expression onNullSet() { return null; }  
                   
    public Expression onOneOrMore(OneOrMoreExp exp) {
	return exp.exp.visit(this);
    }   
                   
    public Expression onOther(OtherExp exp) { 
	return exp.exp.visit(this); 
    }   
                   
    public Expression onRef(ReferenceExp exp) {
	return exp.exp.visit(this);
    } 
                   
    public Expression onSequence(SequenceExp exp) {
	int curIndex = nodeIndex;
	Expression ret = null;
	Vector possibilities = new Vector();
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    nodeIndex = curIndex;
	    ret = children[i].visit(this);
	    if (ret != null) {
		possibilities.add(ret);
	    }
	}
	if (possibilities.size() == 0) {
	    return null;
	} else if (possibilities.size() == 1) {
	    return (Expression)possibilities.get(0);
	} else {
	    return resolveMatch(true,possibilities,curIndex);
	}   
    }

    public Expression onValue(ValueExp exp)  { return null; }
                   
    /**
      * Decides which of the possible Expressions matches the searched 
      * for one and returns it.
      * @param sequence if this resolve should be treated as a sequence
      * (as far as I can tell this doesn't matter, but maybe later)
      * @param possibilities an array of Expressions, one should match
      * nodes[index]
      * @param index the index into the node chain to find a match for
      * @return the Expression which matches nodes[index]
      */
    private Expression resolveMatch(boolean sequence,
				    Vector possibilities,
				    int index) {
	Node parent = nodes[index-1];
	Node relevant = nodes[index];
	NodeList sibs = parent.getChildNodes(); // item(i) getLength()
	Node temp = sibs.item(0);
	int j = 0;
	while (temp.getNodeType() != Node.ELEMENT_NODE 
	       || !temp.getNodeName().equals(relevant.getNodeName())) {
	    j++;
	    temp = sibs.item(j);
	}
	for (int i = 0; i < possibilities.size(); i++) {
	    if (sibs.item(j) == relevant) {
		return (Expression)possibilities.get(i);
	    }
	    j++;
	}
	DEditDebug.println(this,1,"no good");
	return (Expression)possibilities.get(possibilities.size()-1);
    }
  
}
