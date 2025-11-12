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
import model.DEditDebug;

/**
  * Determines if visited Expression is a leaf in the XML tree.
  */
public final class DEditEVIsLeaf implements RELAXExpressionVisitorBoolean {
  
    public boolean onAttPool(AttPoolClause exp) {
	return true;
    }
	
    public boolean onElementRules(ElementRules exp) {
	return true;
    }
	
    public boolean onHedgeRules(HedgeRules exp) {
	return true;
    }
	
    public boolean onTag(TagClause exp) {
	return true;
    }

    public boolean onAnyString() {
	return true;
    }
                   
    public boolean onAttribute(AttributeExp exp) {
	return exp.getContentModel().visit(this);
    } 
	
    public boolean onChoice(ChoiceExp exp) {
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    if (!children[i].visit(this)) {
		DEditDebug.println(this,0,"it's me, onChoice");
		return false;
	    }
	}
	return true;
    } 
	
    public boolean onConcur(ConcurExp exp) {
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    if (!children[i].visit(this)) {
		DEditDebug.println(this,0,"it's me, onConcur");
		return false;
	    }
	}
	return true;
    } 
                   
    public boolean onData(DataExp exp) {
	return exp.except.visit(this); // ???
    }   
                   
    public boolean onElement(ElementExp exp) {
	DEditDebug.println(this,0,"it's me, onElement");
	return false;
    } 
                   
    public boolean onEpsilon() {
	return true;
    } 
                   
    public boolean onInterleave(InterleaveExp exp) {
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    if (!children[i].visit(this)) {
		DEditDebug.println(this,0,"it's me, onInterleave");
		return false;
	    }
	}
	return true;
    }   
                   
    public boolean onList(ListExp exp) {
	return exp.exp.visit(this);
    }  
                   
    public boolean onMixed(MixedExp exp) {
	return exp.exp.visit(this);
	    
    }   
                   
    public boolean onNullSet() {
	return true;
    }  
                   
    public boolean onOneOrMore(OneOrMoreExp exp) {
	return exp.exp.visit(this);
    }   
                   
    public boolean onOther(OtherExp exp) {
	return true;
    }   
                   
    public boolean onRef(ReferenceExp exp) {
	return exp.exp.visit(this);
    } 
                   
    public boolean onSequence(SequenceExp exp) {
	Expression[] children = exp.getChildren();
	for (int i = 0; i < children.length; i++) {
	    if (!children[i].visit(this)) {
		DEditDebug.println(this,0,"it's me, onSequence");
		return false;
	    }
	}
	return true;
    }   
                   
    public boolean onValue(ValueExp exp)  {
	return true;
    }  
}
