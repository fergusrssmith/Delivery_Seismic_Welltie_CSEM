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

import expressionVisitors.*;
import java.io.File;
import java.util.Vector;
import java.util.Hashtable;
import java.net.URL;
import java.lang.Math;
import javax.xml.parsers.*;    
import javax.swing.tree.DefaultTreeModel;
import org.iso_relax.verifier.*;
import org.w3c.dom.*;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;  
import org.xml.sax.SAXParseException;
import com.sun.msv.reader.util.GrammarLoader;
import com.sun.msv.grammar.relax.*;
import com.sun.msv.grammar.*;
import com.sun.msv.verifier.jarv.TheFactoryImpl;

/**
  * Class which wraps a Schema and verifier and all the various things 
  * DEditor uses it for.
  * @author William Ryan
  */
public final class DEditSchema {

    private VerifierFactory factory = new TheFactoryImpl();
    private Schema schema = null;
    private Verifier verifier = null;
    private Grammar grammar = null;
    private String schemaFilename = null;
    private String verifyString = null;
    private final DEditor editor;
    private boolean relaxed = false;
    // key => com.sun.msv.grammar.Expression
    // value => java.util.Vector containing DEditTagInfo
    private Hashtable attTable = new Hashtable();
    // key => com.sun.msv.grammar.Expression
    // value => java.util.Vector containing java.lang.Boolean
    private Hashtable leafTable = new Hashtable();
    // key => com.sun.msv.grammar.NameClassAndExpression
    // value => DEditDataType
    private Hashtable dataTable = new Hashtable();
    /**
      * Static string which is used to say that XML is valid
      */
    public static final String VALID_XML = "XML document is valid.";
    /**
      * Static string which is used to say that XML is invalid
      */
    public static final String INVALID_XML = "XML document is invalid.";
    /**
      * Static string which is used to say that no Schema is loaded
      */
    public static final String NOT_LOADED = "No Schema loaded.";
    // singleton
    private static DEditEVExpression expVisitor = new DEditEVExpression();
    private static DEditEVOptions optVisitor = new DEditEVOptions();
    private static DEditEVIsLeaf leafVisitor = new DEditEVIsLeaf();
    private static DEditEVAttFinder attVisitor = new DEditEVAttFinder();
    private static DEditEVRequired reqVisitor = new DEditEVRequired();
    private static DEditEVData dataVisitor = new DEditEVData();
    // booleans used by helper function to disambiguate generic functionality
    // when necessary
    private static final boolean MIN = true;
    private static final boolean MAX = false;
    
    /**
      * Constructor for DEditSchema which keeps track of what DEditor created 
      * it, and register whether or not it should operate under "relaxed" mode
      * which may allow some illegal choices.
      * @param deditor the DEditor "owner" of this DEditSchema
      * @param relaxed true if the Schema should operate under "relaxed" mode
      */
    public DEditSchema(DEditor deditor, boolean relaxed) {
	this.editor = deditor;
	this.relaxed = relaxed;
    }

    /**
      * Clears all the Hashtable's that this DEditSchema uses.
      */
    private void clearHashtables() {
	attTable.clear();
	leafTable.clear();
	dataTable.clear();
    }

    /**
      * Gets the String name associated with the top level element in an 
      * XML following this DEditSchema's XML Schema.
      * It Assumes top level is first required <xsd:element> tag inside 
      * of <xsd:schema> tag, or last tag if all are optional
      * @return the String name for the top level element in an XML given 
      * this DEditSchema's associated XML Schema
      */
    public String getTopName() {
	Expression top = grammar.getTopLevel();
	reqVisitor.setNodeList(null);
	Vector potentials = (Vector)(top.visit(reqVisitor));
	if (potentials == null || potentials.size() == 0) {
	    return null;
	}
	int i = 1;
	DEditTagInfo ti = (DEditTagInfo)(potentials.get(0));
	while (ti.getMinOccurs() == 0 && i < potentials.size()) {
	    ti = (DEditTagInfo)(potentials.get(i));
	}
	ElementExp eexp = (ElementExp)(ti.getExpression());
	return eexp.getNameClass().toString();
    }

    /**
      * Returns whether or not this DEditSchema is ready to do any validating
      * or generating work.
      * @return true if this DEditSchema can validate and generate, 
      * should be true after successful calls to DEditSchema.loadSchema().
      */
    public boolean isReady() {
	return (verifier != null && grammar != null);
    }

    // This could be changed or a function added to allow loading of remote
    // files for Schemas !!!
    /**
      * Loads a schema and returns true for success although it shouldn't 
      * need to; a DEditSchemaException should be thrown on any error
      * @param schemaFilename the filename of the XML Schema file to be used
      * by this DEditSchema.
      * @return true if this call was a success
      */
    public boolean loadSchema(String schemaFilename) 
    throws DEditSchemaException { 
	try {
	    InputSource schemaIS;
	    if (editor.isApplet()) {
		ClassLoader cl = getClass().getClassLoader();
		URL schemaURL = cl.getResource(schemaFilename);
		schemaIS = new InputSource(schemaURL.toString());
		//String contents = editor.readRemoteFile(schemaFilename);
		//byte[] bytes = contents.getBytes();
		//schemaIS = new InputSource(new java.io.ByteArrayInputStream(bytes));
	    } else {
		ClassLoader cl = getClass().getClassLoader();
		URL schemaURL = cl.getResource(schemaFilename);
		if (schemaURL == null) {
		    System.out.println(schemaFilename);
		    File f = new File(schemaFilename);
		    schemaURL = f.toURL();
		}
		//String schemaIS = schemaURL.toString();
		schemaIS = new InputSource(schemaURL.toString());
	    }
	    grammar = GrammarLoader.loadSchema(schemaIS);
	    schema = factory.compileSchema(schemaIS);
	    verifier = schema.newVerifier();
	    this.schemaFilename = schemaFilename;
	} catch (Exception e) { // restore???
	    schema = null;
	    verifier = null;
	    schemaFilename = null;
	    clearHashtables();
	    e.printStackTrace();
	    throw new DEditSchemaException("Error loading schema: "  
					   + e.toString());
	}
	clearHashtables();
	return (grammar != null);
    }
    
    /**
      * Returns the filename of the XML Schema document used by this
      * DEditSchema.
      */
    public String getFilename() {
	return schemaFilename;
    }

    /**
      * Returns a String updated by any verification attempt which can be
      * used to inform the user of what went wrong.
      * @return the string result of the last verify attempt.
      */
    public String getVerifyString() {
	return verifyString;
    }

    /**
      * Checks to see if passed in DOM XML document is valid and sets 
      * verifyString that may be used to get more information
      * @param d the Document checked for validity
      * @return true if Document d is valid given schema
      */
    public boolean isValid(Document d) {
	if (verifier == null) {
	    verifyString = DEditSchema.NOT_LOADED;
	    return false;
	}
	try {
	    if (verifier.verify(d)) {
		verifyString = DEditSchema.VALID_XML;
		return true;
	    } else {
		verifyString = DEditSchema.INVALID_XML;
		return false;
	    }
	} catch (SAXParseException spe) {
	    verifyString = spe.getMessage();
	    return false;
	} catch (Exception e) {
	    String err = e.toString();
	    verifyString = err;
	    return false;
	}
    }

    /**
      * Queries the Schema to see if the passed-in Node represents a leaf
      * in the XML tree of a proper document.
      * @param node the node being checked to see if it is a leaf
      * @return true if node should be a leaf, false otherwise
      */
    public boolean isLeaf(Node node) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	return isLeaf(exp);
    }

    /**
      * Queries the Schema to see if the passed-in ElementExp is a 
      * schema representation of a leaf in the XML tree.
      * @param exp the ElementExp being checked to see if it is a leaf
      * @return true if node is be a leaf, false otherwise
      */
    public boolean isLeaf(ElementExp exp) {
	if (exp == null) {
	    return true; //???
	}
	Boolean temp = (Boolean)leafTable.get(exp);
	if (temp != null) {
	    return temp.booleanValue();
	}
	boolean ret = exp.getContentModel().visit(leafVisitor);
	leafTable.put(exp,new Boolean(ret));
	return ret;
    }

    /**
      * Returns, as an array of Strings, all potential names for elements 
      * which are children of the element represented by the node passed in.
      * @param node the node to find potential sub-element names for
      * @return an array of potential sub-element names, or null if none are 
      * found
      */
    public String[] getPotentialSubElementNames(Node node) {
	DEditTagInfo[] tagArr = getSubCategoryElems(node,DEditSchema.MAX);
	if (tagArr == null) {
	    return null;
	}
	DEditDebug.println(this,4,"here come the potentials");
	for (int i = 0; i < tagArr.length; i++) {
	    DEditDebug.println(this,4,tagArr[i].toString());
	}
	DEditDebug.println(this,4,"end potentials");
	String[] ret = new String[tagArr.length];
	for (int i = 0; i < ret.length; i++) {
	    ElementExp elexp = (ElementExp)(tagArr[i].getExpression());
	    ret[i] = elexp.getNameClass().toString();
	}
	return ret;
    }

    /**
      * Returns, as an array of Strings, all required names for elements 
      * which are children of the element represented by the node passed in.
      * @param node the node to find required sub-element names for
      * @return an array of required sub-element names, or null if none are 
      * found
      */
    public String[] getRequiredSubElementNames(Node node) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	return getRequiredSubElementNames(node);
    }
    
    /**
      * Returns, as an array of Strings, all required names for elements 
      * which are children of the element represented by the node passed in, 
      * which presumably lines up with the ElementExp passed in.  This function
      * is a helper that prevents retraversal of the tree to get the 
      * relevant expression (in this case ElementExp exp) if you are 
      * recursively running this function.
      * @param node the node to find required sub-element names for
      * @param exp the ElementExp corresponding to node, should == 
      * getRelevantExpression(node)
      * @return an array of required sub-element names, or null if none are 
      * found
      */
    public String[] getRequiredSubElementNames(Node node, ElementExp exp) {
	DEditTagInfo[] tagArr = getSubCategoryElems(node, DEditSchema.MIN,exp);
	if (tagArr == null) {
	    return null;
	}
	String[] ret = new String[tagArr.length];
	for (int i = 0; i < ret.length; i++) {
	    ElementExp elexp = (ElementExp)(tagArr[i].getExpression());
	    ret[i] = elexp.getNameClass().toString();
	}
	return ret;
    }

    /**
      * Gets the actual Expressions required as children of the input Node 
      * and assuming the input ElementExp is the same as  you would get from 
      * calling getRelevantExpression(node).
      * @param node the node to find required sub-elements for
      * @param exp the ElementExp corresponding to node, should == 
      * getRelevantExpression(node)
      * @return an array of required sub-elements, or null if none are 
      * found
      */
    private DEditTagInfo[] getRequiredSubElements(Node node, ElementExp exp) {
	DEditTagInfo[] ret = getSubCategoryElems(node,DEditSchema.MIN,exp);
	if (ret != null) {
	    DEditDebug.println(this,4,"here come the required");
	    for (int i = 0; i < ret.length; i++) {
		DEditDebug.println(this,4,ret[i].toString());
	    }
	    DEditDebug.println(this,4,"end required");
	}
	return ret;
    }
    
    /**
      * Helper which in turn calls the major helper by adding the last piece 
      * of the puzzle -- the relevant Expression given the Node node.
      * See getSubCategoryElems(Node node, boolean elements, boolean min,
      * ElementExp exp) for details.
      */
    private DEditTagInfo[] getSubCategoryElems(Node node, boolean min) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	return getSubCategoryElems(node,min,exp);
    }

    /**
      * Gets those Expressions (which must implement NameClassAndExpression
      * and currently should only be ElementExp or AttributeExp) which are 
      * sub elements of Node node and ElementExp exp (which should == 
      * getRelevantExpression(node)), then pares the array of expressions to 
      * include only those required (minOccurs - occurences in node) or 
      * potential (maxOccurs - occurences in node) depending on the value of 
      * min (true for required).
      * @param node the node to get sub elements for
      * @param min true for querying required, false for potential
      * @param exp == getRelevantExpression(node)
      * @return an array of Expressions which are either potential or 
      * required children of Node node.
      */
    private DEditTagInfo[] getSubCategoryElems(Node node, boolean min,
					       ElementExp exp) {
	DEditDebug.println(this,0,"start getSub...");
	if (exp == null) {
	    return null;
	}
	DEditDebug.println(this,0,"exp non null");
	Vector potentials;
	if (min) {
	    reqVisitor.setNodeList(node.getChildNodes());
	    potentials = (Vector)(exp.getContentModel().visit(reqVisitor));  
	} else {
	    optVisitor.setNodeList(node.getChildNodes());
	    potentials = (Vector)(exp.getContentModel().visit(optVisitor));
	}
	if (potentials != null && potentials.size() != 0) {
	    DEditTagInfo[] taginfos = new DEditTagInfo[potentials.size()];
	    DEditDebug.println(this,4,"we started with....");
	    for (int i = 0; i < taginfos.length; i++) {
		taginfos[i] = (DEditTagInfo)(potentials.get(i));
		DEditDebug.println(this,4,taginfos[i].toString());
	    }
	    NodeList children = node.getChildNodes();
	    for (int i = 0; i < children.getLength(); i++) {
		Node child = children.item(i);
		String name = child.getNodeName();
		if (name != null) {
		    for (int j = 0; j < taginfos.length; j++) {
			DEditTagInfo inspected = taginfos[j];
			NameClassAndExpression nmne 
			    =(NameClassAndExpression)inspected.getExpression();
			String tagname = nmne.getNameClass().toString();
			int currentMin = inspected.getMinOccurs();
			int currentMax = inspected.getMaxOccurs();
			if (tagname.equals(name)) {
			    if (min) {
				taginfos[j]=new DEditTagInfo((Expression)nmne,
							     currentMin-1,
							     currentMax);
			    } else {
				if (currentMax != DEditTagInfo.UNLIMITED) {
				    int newMax = Math.max(0,currentMax-1);
				    taginfos[j] 
					= new DEditTagInfo((Expression)nmne,
							   currentMin,
							   newMax);
				}
			    }
			}
		    }
		}
	    }
	    Vector retVec = new Vector();
	    for (int i = 0; i < taginfos.length; i++) {
		DEditTagInfo ti = taginfos[i];
		DEditDebug.println(this,0,"min " + ti.getMinOccurs() 
				   + " max " + ti.getMaxOccurs());
		if ((min && ti.getMinOccurs() > 0)
		    || (!min && ti.getMaxOccurs() != 0)) {
		    retVec.add(ti);
		}
	    }
	    DEditTagInfo[] retAr = new DEditTagInfo[retVec.size()];
	    retVec.toArray(retAr);
	    for (int i = 0; i < retAr.length; i++) {
		DEditDebug.println(this,2,"get opt adding " 
				   + ((ElementExp)(retAr[i].getExpression())).getNameClass().toString());
	    }
	    DEditDebug.println(this,2,"ret length: " + retAr.length);
	    return retAr;
	} else {
	    DEditDebug.println(this,3,"returning null: no possibilities");
	    return null;
	}
    }

    /**
      * Get the pertinent constraint message for the input Node.
      * @param node the Node to find the constraint message for.
      * @return the constraint message for node's values
      */
    public String getConstraintMessage(Node node) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	return getConstraintMessage(exp);
    }

    /**
      * Gets the pertinent constraint message for the input 
      * NameClassAndExpression.
      * @param exp the NameClassAndExpression (i.e., ElementExp or 
      * AttributeExp) to find the the constraint message for.
      * @return the constraint message for exp's values
      */
    private String getConstraintMessage(NameClassAndExpression exp) {
	if (exp == null) {
	    return null;
	}
	DEditDataType dt = getDataType(exp);
	if (dt == null) {
	    return null;
	}
	return dt.getConstraintMessage();
    }
    
    /**
      * Gets the DEditDataType associated with a Node.
      * @param node the Node whose DEditDataType will be returned.
      * @return the DEditDataType corresponding to node
      */
    public DEditDataType getDataType(Node node) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	return getDataType(exp);
    }

    /**
      * Gets the DEditDataType associated with a NameClassAndExpression.
      * @param nce the NameClassAndExpression whose DEditDataType will be 
      * returned.
      * @return the DEditDataType corresponding to nce
      */
    private DEditDataType getDataType(NameClassAndExpression nce) {
	if (nce == null) {
	    return null;
	}
	DEditDataType dt = (DEditDataType)dataTable.get(nce);
	if (dt == null) {
	    dt = (DEditDataType)(nce.getContentModel().visit(dataVisitor));
	    if (dt == null) {
		dataTable.put(nce,new DEditDataType(DEditDataType.EPSILON,
						    DEditDataType.EMPTY,
						    null,null));
	    } else {
		dataTable.put(nce,dt);
	    }
	}
	if (dt == null || dt.getType() == DEditDataType.EPSILON) {
	    DEditDebug.println(this,3,"returning null data type");
	    return null;
	}
	return dt;
    }

    /**
      * Gets the default value in String form of the leaf Node passed in.
      * This is currently unimplemented and currently unimplementable until
      * I complete the current Grammar code base or find another one.
      * @param node the leaf Node to find the default text value for
      * @return the default value for leaf Node node in String form
      */
    public String getDefaultTextValue(Node node) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	return getDefaultTextValue(exp);
    }

    /**
      * Gets the default value in String form of the leaf ElementExp passed in.
      * This is currently unimplemented and currently unimplementable until
      * I complete the current Grammar code base or find another one.
      * @param exp the leaf ElementExp to find the default text value for
      * @return the default value for leaf ElementExp exp in String form
      */
    public String getDefaultTextValue(ElementExp exp) {
	DEditDataType dt = getDataType(exp);
	if (dt == null) {
	    return null;
	}
	return dt.getDefaultValue();
    }

    /**
      * Gets, as an array, all attributes corresponding to node in the form 
      * of DEditAttributeInfo's.  From this they can be compared to what is 
      * in Node to be turned into the current node status plus information
      * about what else is allowable (though that step is not done here).
      * @param node the Node to get attributes for
      * @return an array of DEditAttributeInfo's containing information 
      * about various attributes for Node node
      */
    public DEditAttributeInfo[] getAttributes(Node node) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	Vector potentials = (Vector)attTable.get(exp);
	if (potentials == null) {
	    potentials = (Vector)(exp.getContentModel().visit(attVisitor));
	    if (potentials != null) {
		attTable.put(exp,potentials);
	    } else {
		attTable.put(exp,new Vector(0)); // so we don't keep looking
	    }
	}
	if (potentials != null && potentials.size() != 0) {
	    DEditAttributeInfo[] ret
		= new DEditAttributeInfo[potentials.size()];
	    for (int i = 0; i < ret.length; i++) {
		DEditTagInfo ti = (DEditTagInfo)(potentials.get(i));
		AttributeExp ae = (AttributeExp)ti.getExpression();
		DEditDataType dt = getDataType(ae);
		String name = ae.getNameClass().toString();
		ret[i] = new DEditAttributeInfo(dt,name,"",false);
	    }
	    return ret;
	} else {
	    return null;
	}
    }

    /**
      * Gets DEditTagInfo's for all possible sub-elements of Node node.
      * @param node the Node whose possible sub-elements' DEditTagInfo's 
      * will be returned.
      * @return an array of DEditTagInfo's corresponding to possible children 
      * of node
      */
    public DEditTagInfo[] getSubElementTagInfos(Node node) {
	ElementExp exp = (ElementExp)getRelevantExpression(node);
	return getSubElementTagInfos(node, exp);
    }
    
    /**
      * Gets DEditTagInfo's for all possible sub-elements of ElementExp.
      * @param exp the ElementExp whose possible sub-elements' DEditTagInfo's 
      * will be returned.
      * @return an array of DEditTagInfo's corresponding to possible children 
      * of exp
      */
    private DEditTagInfo[] getSubElementTagInfos(Node node, ElementExp exp) {
	if (exp == null) {
	    return null;
	}
	optVisitor.setNodeList(node.getChildNodes());
	Vector potentials = (Vector)(exp.getContentModel().visit(optVisitor));
	if (potentials != null && potentials.size() != 0) {
	    DEditTagInfo[] ret = new DEditTagInfo[potentials.size()];
	    potentials.toArray(ret);
	    DEditDebug.println(this,2,"" + ret.length);
	    return ret;
	}
	return null;
    }

    // does optVisitor remove some it shouldn't????
    /**
      * Gets an array of Strings corresponding to all the subelement names 
      * of node, whether they are currently filled or not.
      * @param node the Node to find sub-element names for
      * @return an array of Strings which are the names of all possible 
      * sub-elements of node
      */
    public String[] getAllSubElementNames(Node node) {
	DEditTagInfo[] tiArr = getSubElementTagInfos(node);
	if (tiArr == null) {
	    return null;
	}
	DEditDebug.println(this,3,"getAllSubElementNames");
	Vector holder = new Vector();
	for (int i = 0; i < tiArr.length; i++) {
	    String temp = ((ElementExp)tiArr[i].getExpression()).getNameClass().toString();
	    DEditDebug.println(this,4,"has" + temp);
	    if (!holder.contains(temp)) {
		holder.add(temp);
	    }
	}
	String[] ret = new String[holder.size()];
	holder.toArray(ret);
	return ret;
    }

    /**
      * A widely used helper function that one should be careful not to use 
      * when unneccessary as it has to go all the way up to the top of node,
      * then down it's tree again as well as the Grammar tree.
      * The relevant expression is the Expression in the Grammar tree that 
      * corresponds to node.
      * @param node the org.w3c.dom.Node to find the relevant Expression for
      * @return the relevant Expression for node
      */
    private Expression getRelevantExpression(Node node) {
	Vector nodes = new Vector();
	Node temp = node;
	while (temp != null) {
	    nodes.add(0,temp);
	    temp = temp.getParentNode();
	}
	Node[] nodesArray = new Node[nodes.size()];
	for (int i = 0; i < nodesArray.length; i++) {
	    nodesArray[i] = (Node)(nodes.get(i));
	}
	expVisitor.setNodes(nodesArray);
	Expression exp = grammar.getTopLevel().visit(expVisitor);
	return exp;
    }

    /**
      * DEditor.java uses this to generate a default JTree and DOM tree 
      * given this DEditSchema.  After calling this function the DOM tree 
      * should be valid under the input DEditNode given the schema.
      * @param parent the DEditNode to generate children for.
      * @param model the DefaultTreeModel used to modify the JTree that 
      * DEditNode parent is part of
      */
    public void generateRequiredChildren(DEditNode parent,
					 DefaultTreeModel model) {
	ElementExp exp = (ElementExp)getRelevantExpression(parent.getNode());
	generateRequiredChildren(parent,model,exp);
    }

    /**
      * Helper for generateRequiredChildren(DEditNode, DefaultTreeModel) which 
      * allows recursion without repeated calculation of the relevant 
      * Expression of the current DEditNode.
      * @param parent the DEditNode to generate children for.
      * @param model the DefaultTreeModel used to modify the JTree that 
      * DEditNode parent is part of
      * @param exp the relevant Expression for parent
      */
    private void generateRequiredChildren(DEditNode parent, 
					  DefaultTreeModel model,
					  ElementExp exp) {
	if (exp == null) {
	    return;
	}
	Node node = parent.getNode();
	Document owner = node.getOwnerDocument();
	// default attributes??? pass owner document???
	if (isLeaf(exp)) {
	    String defVal = getDefaultTextValue(exp);
	    if (defVal == null) {
		defVal = "";
	    }
	    DEditDebug.println(this,2,node.getNodeName() + " " + defVal);
	    Text textNode = owner.createTextNode(defVal);
	    node.appendChild(textNode);
	    return;
	}
	DEditTagInfo[] kids = getRequiredSubElements(node,exp);
	if (kids == null) {
	    return;
	}
	/*
	int numSpaces = (parent.getLevel()) * 2;
	char[] spaces = new char[numSpaces];
	for (int i = 0; i < spaces.length; i++) {
	    spaces[i] = ' ';
	}
	String indentSpaces = new String(spaces);		
	String indent = indentSpaces + "  ";
	Node parNode = parent.getNode();
	parNode.appendChild(owner.createTextNode("\n"));
	*/
	for (int i = 0; i < kids.length; i++) {
	    ElementExp kid = (ElementExp)kids[i].getExpression();
	    String kidName = kid.getNameClass().toString();
	    for (int j = 0; j < kids[i].getMinOccurs(); j++) {
		Element newNode = owner.createElement(kidName);
		DEditNode newDENode = new DEditNode(newNode,model);
		//Node textNode = owner.createTextNode(indent);
		//parNode.appendChild(textNode);
		parent.insertAtEnd(newDENode);
		//textNode = owner.createTextNode("\n");
		//parNode.appendChild(textNode);
		generateRequiredChildren(newDENode,model,kid);
	    }
	}
	//parent.getNode().appendChild(owner.createTextNode(indentSpaces));
    }

    /**
      * Gets DEditNode's for creating new tags before / above the input node.
      * @param denode returned DEditNode's will be created immediately before 
      * / above this DEditNode
      * @return DEditNode's for creating new tags before / above the input node
      */
    public String[] getNewTagBeforeNames(DEditNode denode) {
	if (!isReady()) {
	    return null; // free for all???
	}
	final DEditNode parent = (DEditNode)denode.getParent();
	if (parent == null) {
	    return null;
	}
	DefaultTreeModel model = denode.getModel();
	final int curIndex = parent.getIndex(denode);
	String[] opts = getPotentialSubElementNames(parent.getNode());
	if (!relaxed && isValid(denode.getNode().getOwnerDocument())) {
	    DEditNode[] temp = getNewTagItemsVerifyCheck(parent,curIndex,
							 opts,model);
	    String[] ret = new String[temp.length];
	    for (int i = 0; i < ret.length; i++) {
		ret[i] = temp[i].getNode().getNodeName();
	    } 
	    return ret;
	} else {
	    DEditDebug.println(this,6,"relaxed or not valid " 
			       + getVerifyString());
	}
	return opts;
    }
    
    // Streamline with getNewTagBeforeNodes() ???
    /**
      * Gets DEditNode's for creating new tags after / below the input node.
      * @param denode returned DEditNode's will be created immediately after /
      * below this DEditNode
      * @return DEditNode's for creating new tags after / below the input node
      */
    public String[] getNewTagAfterNames(DEditNode denode) {
	if (!isReady()) {
	    return null; // free for all???
	}
	final DEditNode parent = (DEditNode)denode.getParent();
	if (parent == null) {
	    return null;
	}
	DefaultTreeModel model = denode.getModel();
	final int curIndex = parent.getIndex(denode);
	String curName = denode.getNode().getNodeName();
	String[] opts = getPotentialSubElementNames(parent.getNode());
	if (!relaxed && isValid(denode.getNode().getOwnerDocument())) {
	    DEditNode[] temp = getNewTagItemsVerifyCheck(parent,curIndex+1,
							 opts,model);
	    String[] ret = new String[temp.length];
	    for (int i = 0; i < ret.length; i++) {
		ret[i] = temp[i].getNode().getNodeName();
	    } 
	    return ret;
	}
	DEditDebug.println(this,6,"relaxed or not valid");
	return opts;
    }

    /**
      * The newer, preffered way of determining what can be put where.
      * It verifies that the items created will be legal when inserted.  
      * However this is not without hefty calculation costs.
      * @param parent the DEditNode to look for new child DEditNode's for
      * @param index the index of parent's children at which to try inserting 
      * to see what's legal
      * @param opts the list of names of potential new DEditNode's
      * @param model the DefaultTreeModel used to create DEditNode's
      * @return an array of DEditNode's to put in a menu for inserting new 
      * nodes in the tree
      */
    private DEditNode[] getNewTagItemsVerifyCheck(DEditNode parent, int index,
						  String[] opts,
						  DefaultTreeModel model) {
	DEditDebug.println(this,3,"running getNew...VerifyCheck!");
	long allstart = System.currentTimeMillis();
	Vector showable = new Vector();
	Node next = null;
	long genTime = 0;
	long checkTime = 0;
	if (index < parent.getChildCount()) {
	    DEditNode nextDEN = (DEditNode)parent.getChildAt(index);
	    next = nextDEN.getNode();
	}
	Document d = parent.getNode().getOwnerDocument();
	for (int i = 0; i < opts.length; i++) {
	    DEditDebug.println(this,3,"request create >" + opts[i] + "<");
	    Node node = d.createElement(opts[i]);
	    DEditNode newDENode = new DEditNode(node,model);
	    parent.insertAtIndex(newDENode,index,next);
	    long start = System.currentTimeMillis();
	    generateRequiredChildren(newDENode,model);
	    long checkstart = System.currentTimeMillis();
	    genTime += (checkstart - start);
	    if (isValid(newDENode.getNode().getOwnerDocument())) {
		showable.add(newDENode);
	    }
	    checkTime += (System.currentTimeMillis() - checkstart);
	    newDENode.removeSelf();
	}
	DEditDebug.println(this,5,"cost of generation = " + genTime);
	DEditDebug.println(this,5,"cost of check = " + checkTime);
	DEditNode[] ret = new DEditNode[showable.size()];
	showable.toArray(ret);
	long stop = System.currentTimeMillis();
	DEditDebug.println(this,6,"cost of verify: " + (stop - allstart));
	return ret;
    }  
    
    /**
      * Gets the names as Strings of all new child tags which can be created 
      * as children of the input DEditNode
      * @param parent the returned tag names are potential children of this 
      * DEditNode
      * @return the names as Strings of all new child tags which can be created
      */
    public String[] getNewChildTagNames(DEditNode parent) {
	if (parent == null) {
	    return null;
	} else if (!isReady()) {
	    // free for all???
	    return null;
	}
	return getPotentialSubElementNames(parent.getNode());
    }

    // remove ???
    /**
      * Creates and returns an array of DEditNode's with the provided names 
      * to be inserted as children of DEditNode parent at the index provided
      * @param names the names of the DEditNodes to make
      * @param parent the DEditNode that the returned DEditNodes can be 
      * inserted as children of
      * @param index the index into parent's children at which to create 
      * the nodes
      * @param model the DefaultTreeModel used to make and insert DEditNode's
      * @return an array of DEditNode's which can be inserted as children of 
      * DEditNode parent and the provided index
      */
    private DEditNode[] createNodes(String[] names, DEditNode parent,
				    int index, DefaultTreeModel model) {
	long start = System.currentTimeMillis();
	DEditNode[] ret = new DEditNode[names.length];
	Document d = parent.getNode().getOwnerDocument();
	for (int i = 0; i < ret.length; i++) {
	    Node node = d.createElement(names[i]);
	    ret[i] = new DEditNode(node,model);
	    Node nextNode = null;
	    if (index < parent.getChildCount()) {
		DEditNode denode = (DEditNode)parent.getChildAt(index);
		if (denode != null) {
		    nextNode = denode.getNode();
		}
	    }
	    parent.insertAtIndex(ret[i],index,nextNode);
	    generateRequiredChildren(ret[i],model);
	    ret[i].removeSelf();
	} 
	long stop = System.currentTimeMillis();
	DEditDebug.println(this,6,"cost of creating: " + (stop - start));
	return ret;
    }

    /**
      * Gets DEditNode's that DEditNode current can switch to.
      * @param current the node examined to find alternative nodes
      * @param model the DefaultTreeModel used to create DEditNode's
      * @return DEditNode's that DEditNode current can switch to
      */
    public String[] getSwitchTagToNames(final DEditNode current,
					DefaultTreeModel model) {
	final DEditNode parent = (DEditNode)current.getParent();
	if (parent == null) {
	    return null;
	}
	final int index = parent.getIndex(current);
	Document d = current.getNode().getOwnerDocument();
	if (!isValid(d)) {
	    return null;
	}
	current.removeNode();
	String[] opts = getPotentialSubElementNames(parent.getNode());
	Node next = null;
	if (index+1 < parent.getChildCount()) {
	    DEditNode nextDEN = (DEditNode)parent.getChildAt(index+1);
	    next = nextDEN.getNode();
	}
	if (relaxed) {
	    parent.reinsertNodeBefore(current,next);
	    return opts;
	}
	long start = System.currentTimeMillis();
	Vector switchable = new Vector();
	for (int i = 0; i < opts.length; i++) {
	    Node n = d.createElement(opts[i]);
	    DEditNode newDENode = new DEditNode(n,model);
	    parent.insertAtIndex(newDENode,index,next);
	    generateRequiredChildren(newDENode,model);
	    if (isValid(d)) {
		switchable.add(newDENode);
	    }
	    newDENode.removeSelf();
	}
	parent.reinsertNodeBefore(current,next);
	String[] ret = new String[switchable.size()];
	for (int i = 0; i < ret.length; i++) {
	    ret[i] = ((DEditNode)switchable.get(i)).getNode().getNodeName();
	}
	long stop = System.currentTimeMillis();
	DEditDebug.println(this,6,"cost of verify: " + (stop - start));
	return ret;
    }

    /**
      *
      * Creates a new child element tag in the XML.
      * @param parent the parent to create a new child node for
      * @param newTagName the name of the new node to be created
      * @param model the DefaultTreeModel used to create and insert DEditNode's
      */
    public void createNewChildTag(DEditNode parent, String newTagName,
				   DefaultTreeModel model) {
	Document d = parent.getNode().getOwnerDocument();
	if (isValid(d)) {
	    createNewChildTagVerifyCheck(parent, newTagName, model);
	    return;
	}
	String[] full = getAllSubElementNames(parent.getNode());
	if (newTagName != null) {
	    Node node = d.createElement(newTagName);
	    DEditNode newDENode = new DEditNode(node,model);
	    boolean doneBuild = false;
	    String[] existing = new String[parent.getChildCount()];
	    for (int i = 0; i < existing.length; i++) {
		DEditNode child = (DEditNode)parent.getChildAt(i);
		existing[i] = child.getNode().getNodeName();
	    }
	    int j = 0;
	    if (existing.length == 0) {
		insertEndGenerate(parent,newDENode,model);
		return;
	    }
	    for (int i = 0; i < full.length; i++) {
		if (full[i].equals(newTagName)) {
		    Node n = ((DEditNode)parent.getChildAt(j)).getNode();
		    insertIndexGenerate(parent,newDENode,j,n,model);
		    return;
		} else {
		    while (full[i].equals(existing[j])) {
			j++;
			if (j == existing.length) {
			    insertEndGenerate(parent,newDENode,model);
			    return;
			}
		    }
		} 
	    }
	    DEditDebug.println(this,7,"Failed to create node w/bad logic.");
	}
    }

    /**
      * Helper which inserts newDENode as a child of parent at the end of 
      * its children and generates required sub elements.
      * @param parent the node in which to insert at the end of its children
      * @param newDENode the node to insert
      * @param model the DefaultTreeModel used for DEditNode creation and 
      * insertion
      */
    private void insertEndGenerate(DEditNode parent, DEditNode newDENode,
				   DefaultTreeModel model) {
	parent.insertAtEnd(newDENode);
	generateRequiredChildren(newDENode,model);
    }

    /**
      * Helper which inserts newDENode as a child of parent at a certain index
      * of its children and generates required sub elements.
      * @param parent the node in which to insert at a certain index
      * @param newDENode the node to insert
      * @param index the index at which to insert newDENode
      * @param ref the Node == parent.getChildAt(index).getNode() before 
      * insertion of newDENode
      * @param model the DefaultTreeModel used for DEditNode creation and 
      * insertion
      */
    private void insertIndexGenerate(DEditNode parent, DEditNode newDENode,
				     int index, Node ref,
				     DefaultTreeModel model) {
	parent.insertAtIndex(newDENode,index,ref);
	generateRequiredChildren(newDENode,model);
    }

    /**
      * Creates a new child for DEditNode parent with the provided name 
      * newTagName and chooses location by making sure it is valid there.
      * @param parent the DEditNode who will have a child inserted
      * @param newTagName the name of the DEditNode that will be created
      * @param model the DefaultTreeModel used to create and insert DEditNode's
      */
    private void createNewChildTagVerifyCheck(DEditNode parent, 
					      String newTagName,
					      DefaultTreeModel model) {
	Document d = parent.getNode().getOwnerDocument();
	for (int i = 0; i < parent.getChildCount(); i++) {
	    Node node = d.createElement(newTagName);
	    DEditNode newDENode = new DEditNode(node,model);
	    DEditNode child = (DEditNode)parent.getChildAt(i);
	    parent.insertAtIndex(newDENode,i,child.getNode());
	    generateRequiredChildren(newDENode,model);
	    if (isValid(d)) {
		return;
	    }
	    newDENode.removeSelf();
	}
	Node node = d.createElement(newTagName);
	DEditNode newDENode = new DEditNode(node,model);
	parent.insertAtEnd(newDENode);
	generateRequiredChildren(newDENode,model);
	if (isValid(d)) {
	    return;
	}
	newDENode.removeSelf();
	DEditDebug.println(this,7,"Failed to create node w/expensive logic.");
    }

    /**
      * Prints out a tree representation of the XML Schema for this DEditSchema
      * to standard out.  Warning:  this can be quite long.  This function 
      * is called by DEditSchema.main() which functions as a way to inspect 
      * a Schema you have created.
      */
    public void printoutSchema() {
	Expression exp = grammar.getTopLevel();
	ExpressionPool ep = grammar.getPool();
	System.out.println(exp);
	System.out.println(ep);
	DEditEVPrintout ev = new DEditEVPrintout();
	exp.visit(ev);
    }

    /**
      * Prints out the input schema as an XML tree to standard out.
      * @param args arguments -- the first should be the schema to use
      */
    public static void main(String[] args) {
	DEditSchema es = new DEditSchema(null,false);
	String schemaFile = null;
	if (args.length > 0) {
	    schemaFile = args[0];
	}
	if (schemaFile == null) {
	    System.out.println("No schema file provided.  Quiting.");
	}
	try {
	    boolean huh = es.loadSchema(schemaFile);
	    System.out.println(huh);
	} catch (DEditSchemaException e) {
	    System.out.println(e.toString());
	    e.printStackTrace();
	    System.out.println("aborting");
	    System.exit(0);
	}
	es.printoutSchema();
    }
}
