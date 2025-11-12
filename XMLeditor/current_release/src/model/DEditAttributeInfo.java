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

/**
  * Class which represents information about attributes.  This class might 
  * be extended to include DEditTagInfo and / or other relevant information 
  * as it comes up.
  */
public class DEditAttributeInfo {

    /**
      * The DEditDatatype associated with this DEditAttributeInfo
      */
    public DEditDataType datatype;
    /**
      * The name of the attribute
      */
    public String name;
    /**
      * The value of the attribute
      */
    public String value;
    /**
      * Whether or not the corresponding node has this attribute, as in 
      * actually appears with the attribute, not can have
      */
    public boolean exists; 

    /**
      * Makes a new DEditAttributeInfo with given paramet.  In the XML it 
      * would appear <tag name="value">...</tag> where name is String name 
      * and value is String value.  The type of value should be one of 
      * DEditAttributeInfo.STRING, .BOOLEAN, .INTEGER, or .FLOAT.
      * @param type the DEditDataType of the attribute
      * @param name the name of the attribute
      * @param value the value of the attribute
      * @param exists true if the attribute appears in the tag
      */
    public DEditAttributeInfo(DEditDataType type, String name, String value,
			      boolean exists) {
	this.name = name;
	this.value = value;
	this.datatype = type;
	this.exists = exists;
    }

    /**
      * Returns a deep copy of this DEditAttributeInfo.
      * @return a deep copy of this DEditAttributeInfo
      */
    public DEditAttributeInfo copy() {
	return new DEditAttributeInfo(datatype,name,value,exists);
    }
}
