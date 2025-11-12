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

import java.util.Vector;

/**
  * XML Editor version of a schema simple data type w/ or w/o restriction
  */
public class DEditDataType {

    /**
      * Type to be used for a fixed value (fixed="...").
      */
    public static final int FIXED = 1;
    /**
      * Type to be used for a constrained value (type="xsd:integer" etc.) 
      * which may have limits (minimum value, maximum value, etc.) 
      * specified in the constructor.
      */
    public static final int CONSTRAINT = 2;
    /**
      * Type to be used when value must be one of a list of values specified 
      * in the constructor.
      */
    public static final int CHOICE = 3;
    /**
      * Special type which represents no value and is used only to 
      * facilitate recognizing default values via my dirty hack 
      * (see EditorSchemaREADME)
      */
    public static final int EPSILON = 4;
    private final int type; // state variable for which of the above types 
    // this DEditDataType is

    /**
      * An integer signifying that this DEditDataType's actual data is of type
      * String.
      */
    public static final int STRING = 100;
    /**
      * An integer signifying that this DEditDataType's actual data is of type
      * boolean.
      */
    public static final int BOOLEAN = 101;
    /**
      * An integer signifying that this DEditDataType's actual data is of type
      * integer.
      */
    public static final int INTEGER = 102;
    /**
      * An integer signifying that this DEditDataType's actual data is of type
      * float.
      */
    public static final int FLOAT = 103;
    /**
      * An integer signifying that this DEditDataType's actual data is of type
      * URI (although this is treated == file by the editor).
      */
    public static final int URI = 104;
    /**
      * An integer signifying that this DEditDataType represents a null value, 
      * or no data.
      */
    public static final int EMPTY = 105;
    private final int datatype;// state variable for which of the above 
    // datatypes this DEditDataType is

    private final Vector values; // either constraint values or choice values
    private final String defaultValue; // the default value when this 
    // DEditDataType occurs
    private String constrMessage = null; // a message explaining constraints 
    // on the value

    /**
      * Creates a DEditDataType given the pertinent information
      * @param type the general type of the value to be edited with respect 
      * to how it is edited, one of DEditDataType. FIXED, CONSTRAINT, 
      * CHOICE, or EPSILON.
      * @param datatype the basic type of data in terms of int, boolean, 
      * float, string, etc.  Currently supported are:  
      *    DEditDataType.STRING (xsd:string)
      *    DEditDataType.BOOLEAN (xsd:boolean)
      *    DEditDataType.INTEGER (xsd:integer)
      *    DEditDataType.FLOAT (xsd:float)
      *    DEditDataType.URI (xsd:anyURI)
      * For now, anything else is treated as an xsd:string.
      */
    public DEditDataType(int type, int datatype, String defaultValue,
			 Vector values) {
	this.type = type;
	this.datatype = datatype;
	this.defaultValue = defaultValue;
	this.values = values;
    }
    
    /**
      * Gets the editing type of this DEditDataType - FIXED, CONSTRAINT,
      * CHOICE, or EPSILON
      * @return the editing type of this DEditDataType
      */
    public int getType() {
	return type;
    }

    /**
      * Gets the datatype of this DEditDataType - STRING, INTEGER, FLOAT, 
      * URI, BOOLEAN, or EMPTY
      * @return the datatype of this DEditDataType
      */
    public int getDataType() {
	return datatype;
    }

    /**
      * Gets the constraint or choice value (depends on this.type) at the 
      * index i in the Vector used for choices or constraints, or null if
      * there is none (no exceptions are thrown).
      * @param i the index of the value to get
      * @return a String representing a value constraint or choice
      */
    public String getValueAt(int i) {
	if (values == null || i < 0 || i > values.size()) {
	    return null;
	}
	return (String)(values.get(i));
    }

    /**
      * Gets the default value for this DEditDataType as a String.
      * @return the default value for this DEditDataType
      */
    public String getDefaultValue() {
	return defaultValue;
    }

    /**
      * Converts and returns the constraint or choice values into an array 
      * which of type String.
      * @param arr the array of Strings to put the values into
      * @return ignore it unless you like casting individual array elements,
      * it is the result of values.toArray(arr) where values is a Vector
      */
    public Object[] valsToArray(String[] arr) {
	if (values != null) {
	    return values.toArray(arr);
	} else {
	    return null;
	}
    }
    
    /**
      * Gets the number of constraint or choice values available.
      * @return the number of constraint or choice values available.
      */
    public int numValues() {
	if (values != null) {
	    return values.size();
	} else {
	    return 0;
	}
    }

    /**
      * Checks to see if the value provided would be valid in for 
      * this DEditDataType.
      * @param value the value to check
      * @return true if value would be valid for this DEditDataType.
      */
    public boolean isValid(String value) {
	if (type == FIXED) {
	    return value.equals(getValueAt(0));
	} else if (type == CONSTRAINT) {
	    try {
		if (datatype == FLOAT) {
		    float floatVal = Float.parseFloat(value);
		    if (numValues() < 2) {
			return true;
		    }
		    float lowerBound = Float.parseFloat(getValueAt(0));
		    float upperBound = Float.parseFloat(getValueAt(1));
		    return (lowerBound < floatVal && floatVal < upperBound);
		} else if (datatype == INTEGER) {
		    int intVal = Integer.parseInt(value);
		    if (numValues() < 2) {
			return true;
		    }
		    int lowerBound = Integer.parseInt(getValueAt(0));
		    int upperBound = Integer.parseInt(getValueAt(1));
		    return (lowerBound < intVal && intVal < upperBound);
		} else if (datatype == BOOLEAN) {
		    return (value.equals("false") || value.equals("true"));
		} else if (datatype == URI) {
		    return true; 
		    //for now, can't mess up since FileChooser is used
		} else if (datatype == STRING) {
		    return true; // can't set constraint on a string
		} else if (datatype == EMPTY) {
		    return (value == null || value.equals(""));
		}
	    } catch (NumberFormatException nfe) {
		return false;
	    }
	    return false;
	} else if (type == CHOICE) {
	    boolean ret = false;
	    for (int i = 0; i < numValues(); i++) {
		if (value.equals(getValueAt(i))) {
		    ret = true;
		}
	    }
	    return ret; 
	} else {
	    return value.equals(""); // type EPSILON
	}
    }

    /**
      * Gets the constraint message for this DeditDataType.
      * @return the constraint message for this DEditDataType.
      */
    public String getConstraintMessage() {
	if(constrMessage == null) {
	    constrMessage = makeConstraintMessage();
	}
	return constrMessage;
    }

    /**
      * Lazy construction of constraint message.  Returns a constraint message 
      * created for this DEditDataType.
      * @return a constraint message made for this DEditDataType.
      */
    private String makeConstraintMessage() {
	String ret = "";
	if (datatype == DEditDataType.STRING) {
	    ret += "value must be a string";
	} else if (datatype == DEditDataType.FLOAT) {
	    ret += "value must be a float";
	} else if (datatype == DEditDataType.INTEGER) {
	    ret += "value must be an integer";
	} else if (datatype == DEditDataType.BOOLEAN) {
	    ret += "value must be a boolean";
	} else if (datatype == DEditDataType.URI) {
	    ret += "value must be a URI";
	} else if (datatype == DEditDataType.EMPTY) {
	    ret += "node can have no value.";
	}
	
	if (type == DEditDataType.FIXED) {
	    ret = "Value must equal " + getValueAt(0);
	} else if (type == DEditDataType.CONSTRAINT) {
	    if (numValues() >= 2) {
		ret += ", between " + getValueAt(0) + " and " 
		    + getValueAt(1);
	    } 
	} else if (type == DEditDataType.CHOICE) {
	    ret += ", one of < ";
	    for (int i = 0; i < numValues(); i++) {
		ret += getValueAt(i) + " ";
	    }
	    ret  += ">";
	}
	ret += ".";
	return ret;
    }
	
}
