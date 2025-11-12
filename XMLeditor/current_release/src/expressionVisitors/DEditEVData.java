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
import com.sun.msv.datatype.xsd.*;
import com.sun.msv.grammar.*;
import java.util.Vector;
import java.util.Iterator;
import model.DEditDataType;
import model.DEditDebug;

/**
  * Class which finds the DEditDataType for a given Expression and returns it.
  */
public final class DEditEVData implements RELAXExpressionVisitor {
  
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
	
    /**
      * A choice for data should come up primarily to interpret the hacked
      * method of including the default value in the msv representation of a
      * Schema.
      * @param exp a ChoiceExp encountered searching for the kind of data 
      * that belongs in a node
      * @return an appropriate DEditDataType
      */
    public Object onChoice(ChoiceExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	Vector holder = new Vector();
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    DEditDataType subRet = (DEditDataType)(children[i].visit(this));
	    if (subRet != null) {
		holder.add(subRet);
	    }
	}
	if (holder.size() > 1) {
	    DEditDataType[] dtArr = new DEditDataType[holder.size()];
	    holder.toArray(dtArr);
	    return synthesize(dtArr);
	} else if (holder.size() == 1) {
	    DEditDataType temp = (DEditDataType)holder.get(0);
	    if (temp.getType() == DEditDataType.EPSILON) {
		return null;
	    } else {
		return temp;
	    }
	} else {
	    return null;
	}
    }

    /**
      * Creates a DEditDataType out of other data types.  Mostly for getting
      * default values right.
      * @param dtArr an array of DEditDataType's to combine into one
      * @return a DEditDataType which is the combination of input ones
      */
    private DEditDataType synthesize(DEditDataType[] dtArr) {
	boolean hasFixed = false;
	boolean hasEpsilon = false;
	boolean hasConstraint = false;
	String defVal = null;
	int dataType = -1;
	Vector values = null;
	boolean allStrings = false;
	Vector choices = new Vector();
	DEditDebug.println(this,0,"start " + dtArr.length);
	for (int i = 0; i < dtArr.length; i++) {
	    DEditDebug.println(this,0,"loop " + i);
	    DEditDataType dt = dtArr[i];
	    int type = dt.getType();
	    if (type == DEditDataType.FIXED) {
		String temp = dt.getValueAt(0);
		hasFixed = true;
		defVal = temp;
		if (!choices.contains(temp)) {
		    choices.add(temp);
		}
		DEditDebug.println(this,0,"fixed val " + temp);
	    } else if (type == DEditDataType.EPSILON) {
		hasEpsilon = true;
	    } else if (type == DEditDataType.CONSTRAINT) { // no more than 1
		hasConstraint = true;
		dataType = dt.getDataType();
		values = new Vector(); 
		for (int j = 0; j < dt.numValues(); j++) {
		    values.add(dt.getValueAt(j));
		}
		DEditDebug.println(this,0,"constraint at " + i);
	    } else if (type == DEditDataType.CHOICE) {
		allStrings = (dt.getDataType() == DEditDataType.STRING);
		for (int j = 0; j < dt.numValues(); j++) {
		    String temp = dt.getValueAt(j);
		    if (!choices.contains(temp)) {
			choices.add(temp);
		    }
		}
		DEditDebug.println(this,0,"choice at " + i);
	    } else {
		DEditDebug.println(this,0,"what the flivem flavem?");
	    }
	}
	if (hasConstraint) {
		return new DEditDataType(DEditDataType.CONSTRAINT,
					 dataType,defVal,values);
	} 
	if (hasFixed) {
	    if (allStrings) {
		return new DEditDataType(DEditDataType.CHOICE,
					 DEditDataType.STRING,defVal,
					 choices);
	    } else {
		return new DEditDataType(DEditDataType.CHOICE,
					 DEditDataType.BOOLEAN,
					 // marked as temporary for 
					 // unknown reason
					 defVal,choices);
	    }   
	}
	DEditDebug.println(this,4,"Shouldn't ever have this!  ERROR!");
	return null; // 
    }

    /**
      * I have never encountered this sort of expression so it might be
      * handled incorrectly here.
      */
    public Object onConcur(ConcurExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	Vector holder = new Vector();
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    DEditDataType subRet = (DEditDataType)(children[i].visit(this));   
	    if (subRet != null) {
		holder.add(subRet);
	    }
	}
	// how to handle fusing???
	if (holder.size() > 0) {
	    return holder.get(0);
	} else {
	    return null;
	}
    } 
            
    // look at datatype always!!!
    /**
      * The function in which the actual data specification is found.
      * @param exp the DataExp telling the program what kind of data should be
      * found in the node
      * @return a DEditDataType corresponding to the correct kind of data
      */
    public Object onData(DataExp exp) {
	Vector vec = null;
	String typeOrVal = exp.name.localName;
	int type = DEditDataType.CONSTRAINT;
	int datatype = DEditDataType.STRING;
	if (typeOrVal.startsWith("boolean")) {
	    vec = new Vector();
	    vec.add("false");
	    vec.add("true");
	    type = DEditDataType.CHOICE;
	    datatype = DEditDataType.BOOLEAN;
	} else if (typeOrVal.startsWith("float")) {
	    type = DEditDataType.CONSTRAINT;
	    datatype = DEditDataType.FLOAT;
	} else if (typeOrVal.startsWith("integer")) {
	    type = DEditDataType.CONSTRAINT;
	    datatype = DEditDataType.INTEGER;
	} else if (typeOrVal.startsWith("anyURI")) {
	    type = DEditDataType.CONSTRAINT;
	    datatype = DEditDataType.URI;
	} else { // add on here
	    type = DEditDataType.CONSTRAINT;
	    datatype = DEditDataType.STRING;
	}
	
	if (typeOrVal.endsWith("-derived")) {
	    XSDatatype dt = (XSDatatype)exp.dt;
	    return resolveDT(dt,datatype);
	} else {
	    return new DEditDataType(type,datatype,null,vec);
	}
    }   
                   
    /**
      * A somewhat ugly helper function to resolve more complicated data types
      * (integers with a range, only certain strings, etc.).
      * @param dt an XSDataType which specifies how the data type is restricted
      * @param datatype the general type of data (int, float, string, etc.)
      * @return a DEditDataType which the appropriate restrictions
      */
    private DEditDataType resolveDT(XSDatatype dt,int datatype) {
	Vector vals = new Vector();
	int type = DEditDataType.CONSTRAINT;
	EnumerationFacet ef = (EnumerationFacet)dt.getFacetObject("enumeration");
	if (ef != null && ef.values != null) {
	    Iterator it = ef.values.iterator();
	    while (it.hasNext()) {
		vals.add(it.next().toString());
	    }
	    type = DEditDataType.CHOICE;
	    return new DEditDataType(type,datatype,null,vals);
	}
	MinInclusiveFacet mnif = (MinInclusiveFacet)dt.getFacetObject("minInclusive");
	if (mnif != null) {
	    vals.add(mnif.limitValue.toString());
	    type = DEditDataType.CONSTRAINT;
	}
	MaxInclusiveFacet mxif = (MaxInclusiveFacet)dt.getFacetObject("maxInclusive");
	if (mxif != null) {
	    vals.add(mxif.limitValue.toString());
	    type = DEditDataType.CONSTRAINT;
	}

	if (vals.size() == 0) {
	    return null;
	} else return new DEditDataType(type,datatype,null,vals);
    }

    public Object onElement(ElementExp exp) {
	return null;
    }
                   
    public Object onEpsilon() {
	return new DEditDataType(DEditDataType.EPSILON,
				 DEditDataType.EMPTY,null,null);
    } 
              
    /**
      * I have not encountered this type of node looking for data so this 
      * function might not handle it correctly all the time.
      */
    public Object onInterleave(InterleaveExp exp) {
	String myIndent = indent;
	Expression[] children = exp.getChildren();
	Vector holder = new Vector();
	for (int i = 0; i < children.length; i++) {
	    indent = myIndent + " ";
	    DEditDataType subRet = (DEditDataType)(children[i].visit(this));   
	    if (subRet != null) {
		holder.add(subRet);
	    }
	}
	// how to handle fusing???
	if (holder.size() > 0) {
	    return holder.get(0);
	} else {
	    return null;
	}
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
	DEditDebug.println(this,4,"shouldn't hit this -- one or more");
	return exp.exp.visit(this);
    }   
    
    public Object onOther(OtherExp exp) {
	DEditDebug.println(this,4,"shouldn't hit this -- other");
	return exp.exp.visit(this);
    }   
    
    public Object onRef(ReferenceExp exp) {
	return exp.exp.visit(this);
    } 
    
    /**
      * Generally, this sort of node should not be encountered looking for
      * data, but if it were this function only pays attention to the first
      * thing in the sequence.
      */
    public Object onSequence(SequenceExp exp) {
	Expression[] children = exp.getChildren();
	Vector holder = new Vector();
	for (int i = 0; i < children.length; i++) {
	    DEditDataType subRet = (DEditDataType)(children[i].visit(this));   
	    if (subRet != null) {
		holder.add(subRet);
	    }
	}
	// how to handle fusing???
	if (holder.size() > 0) {
	    return holder.get(0);
	} else {
	    return null;
	}
    }   

    /**
      * Returns a DEditDataType for "fixed" data -- usually meaning that it 
      * is in a choice so that it functions as default.
      * @param exp the ValueExp from whch the string value is retrieved
      * @return a DEditDataType which has a fixed value
      */
    public Object onValue(ValueExp exp)  {
	String name = exp.name.localName;
	String value = exp.value.toString();
	if (name.equals("token")) {
	    Vector vec = new Vector();
	    vec.add(value);
	    return new DEditDataType(DEditDataType.FIXED,DEditDataType.STRING,
				     null,vec);
	} else {
	    // name should always be "token" when looking for data
	    return null;
	}
    }  
}
