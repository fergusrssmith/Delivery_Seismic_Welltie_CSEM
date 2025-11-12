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
import com.sun.msv.grammar.xmlschema.*;
import com.sun.msv.datatype.xsd.*;
import java.util.Iterator;
import org.relaxng.datatype.Datatype;

/**
  * Class which prints out a representation of the passed in Expression and 
  * its children to standard out.  Call resetIndent() before traversing
  * a new tree.
  */
public final class DEditEVPrintout implements RELAXExpressionVisitorVoid {
  
    private String indent = "";
   
    /**
      * Resets the indent.  Call before each new .visit() after the first.
      */
    public void resetIndent() {
	indent = "";
    }
    
    public void onAttPool(AttPoolClause exp) {
	System.out.println(indent + "<attPool>");
	onRef(exp);
	System.out.println(indent + "</attPool>");
    }
    
    public void onElementRules(ElementRules exp) {
	System.out.println(indent + "<elementRules>");
	onRef(exp);
	System.out.println(indent + "</elementRules>");
    }	
    
    public void onHedgeRules(HedgeRules exp) {
	System.out.println(indent + "<hedgeRules>");
	onRef(exp);
	System.out.println(indent + "</hedgeRules>");
    }	
    
    public void onTag(TagClause exp) {
	System.out.print(indent + "<onTag> name = ");
	System.out.println(exp.nameClass.toString());
	onRef(exp);
	System.out.println(indent + "</onTag>");
    }	
    
    public void onAnyString() {
	System.out.println(indent + "<anyString/>");
    }
    
    public void onAttribute(AttributeExp exp) {
	String myIndent = indent;
	indent += " ";
	System.out.print(myIndent + "<attribute>");
	NameClass nmc = exp.getNameClass();
	System.out.println(nmc.toString());
	exp.getContentModel().visit(this);
	System.out.println(myIndent + "</attribute>");
    }
    
    public void onChoice(ChoiceExp exp) {
	String myIndent = indent;
	System.out.print(myIndent + "<choice> ");
	Expression[] children = exp.getChildren();
	System.out.println(children.length);
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    children[i].visit(this);
	    //System.out.println(children[i]);
	}
	System.out.println(myIndent + "</choice>");
    }
    
    public void onConcur(ConcurExp exp) {
	String myIndent = indent;
	System.out.println(myIndent + "<concur>");
	Expression[] children = exp.getChildren();
	System.out.println(children.length);
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    children[i].visit(this);
	    //System.out.println(children[i]);
	}
	System.out.println(myIndent + "</concur>");
    }
    
    public void onData(DataExp exp) {
	String myIndent = indent;
	indent += " ";
	System.out.println(myIndent + "<data> "+ exp);
	System.out.println(indent + exp.name.localName);
	System.out.println(indent + exp.name.namespaceURI);
	System.out.println(indent + exp.dt);
	Datatype dt = exp.dt;
	if (dt instanceof XSDatatype) { // should always be true for us
	    viewXSDT((XSDatatype)dt);
	}
	exp.except.visit(this);
	System.out.println(myIndent + "</data>");
    }
    
    private void viewXSDT(XSDatatype dt) {
	EnumerationFacet ef = (EnumerationFacet)dt.getFacetObject("enumeration");
	if (ef != null && ef.values != null) {
	    Iterator it = ef.values.iterator();
	    System.out.println(indent + "enumeration");
	    while (it.hasNext()) {
		System.out.println(indent + it.next());
	    }
	}
	MaxInclusiveFacet mxif = (MaxInclusiveFacet)dt.getFacetObject("maxInclusive");
	if (mxif != null) {
	    System.out.println(indent + "maxInclusive " + mxif.limitValue);
	}
	MinInclusiveFacet mnif = (MinInclusiveFacet)dt.getFacetObject("minInclusive");
	if (mnif != null) {
	    System.out.println(indent + "minInclusize " + mnif.limitValue);
	}
    }

    public void onElement(ElementExp exp) {
	String myIndent = indent;
	indent += " ";
	System.out.print(myIndent + "<element>");
	NameClass nmc = exp.getNameClass();
	System.out.println("name = " + nmc.toString());
	if (exp instanceof ElementDeclExp.XSElementExp) {
	    ElementDeclExp.XSElementExp xsexp 
		= (ElementDeclExp.XSElementExp)exp;
	    ElementDeclExp ede = xsexp.parent;
	    ReferenceExp body = ede.body;
	    //System.out.println(indent + body);
	    for (int i = 0; i < xsexp.identityConstraints.size(); i++) {
		IdentityConstraint ic = (IdentityConstraint)(xsexp.identityConstraints.get(i));
		System.out.println("ic: " + ic.localName);
	    }
	} else {
	    //System.out.println(exp);
	}
	exp.getContentModel().visit(this);
	System.out.println(myIndent + "</element>");
    } 
    
    public void onEpsilon() {
	System.out.println(indent + "<epsilon/>");
    } 

    public void onInterleave(InterleaveExp exp) {
	String myIndent = indent;
	System.out.println(myIndent + "<interleave>");
	Expression[] children = exp.getChildren();
	System.out.println(children.length);
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    children[i].visit(this);
	    //System.out.println(children[i]);
	}
	System.out.println(myIndent + "</interleave>");
    }   
                   
    public void onList(ListExp exp) {
	String myIndent = indent;
	System.out.println(myIndent + "<list>");
	indent += " ";
	exp.exp.visit(this);
	System.out.println(myIndent + "</list>");
    }

    public void onMixed(MixedExp exp) {
	String myIndent = indent;
	System.out.println(myIndent + "<mixed>");
	indent += " ";
	exp.exp.visit(this);
	System.out.println(myIndent + "</mixed>");
    }   
    
    public void onNullSet() {
	System.out.println(indent + "<nullSet/>");
    }  

    public void onOneOrMore(OneOrMoreExp exp) {
	String myIndent = indent;
	indent += " ";
	System.out.println(myIndent + "<oneOrMore>");
	exp.exp.visit(this);
	System.out.println(myIndent + "</oneOrMore>");
    }   
    
    public void onOther(OtherExp exp) {
	String myIndent = indent;
	indent += " ";
	System.out.println(myIndent + "<other>");
	exp.exp.visit(this);
	System.out.println(myIndent + "</other>");
    }   
                   
    public void onRef(ReferenceExp exp) {
	//String myIndent = indent;
	//indent += " ";
	//System.out.print(myIndent + "<refexp>");
	//System.out.println("  refexp name " + exp.name);
	//System.out.println(myIndent + "ref to " + exp.exp);
	exp.exp.visit(this);
	//System.out.println(myIndent + "</refexp>");
    } 
    
    public void onSequence(SequenceExp exp) {
	String myIndent = indent;
	System.out.print(myIndent + "<sequence> ");
	Expression[] children = exp.getChildren();
	System.out.println(children.length);
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    children[i].visit(this);
	    //System.out.println(children[i]);
	}
	System.out.println(myIndent + "</sequence>");
    }   
    
    public void onValue(ValueExp exp)  {
	String myIndent = indent;
	indent += " ";
	System.out.println(myIndent + "<value>");
	System.out.println(indent + "name= " + exp.name.localName);
	System.out.println(indent + "value= " + exp.value);
	System.out.println(myIndent + "</value>");
    }
}
