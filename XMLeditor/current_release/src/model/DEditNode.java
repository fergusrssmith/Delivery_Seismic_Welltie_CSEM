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

import java.util.regex.*;
import org.w3c.dom.*;
import javax.swing.tree.*;

/**
  * Class which represents a node in the JTree that displays the Nodes in the 
  * DOM Document.  Only Nodes of type ELEMENT should be paired with a 
  * DEditNode.
  */
public class DEditNode extends DefaultMutableTreeNode {

    // the tree model to use
    private DefaultTreeModel model = null;
    // the node with actual XML data
    private Node node = null;
    /**
      * Singleton used to represent an empty node which is a holder until 
      * model.setRoot() is called with something useful.
      */
    public static final DEditNode EMPTY_NODE = new DEditNode();

    // private singleton constructor
    private DEditNode() {
	super("");
    }

    /**
      * Constructor for DEditNode that assigns a model and XML Node 
      * to this JTree usable node.  It auto expands and constructs children
      * elements given a Node.
      * @param node the Node to be paired with this JTree node.  Should be
      * of type ELEMENT
      * @param model the TreeModel to use in insertion operations.
      */
    public DEditNode(Node node, DefaultTreeModel model) {
	super(node.getNodeName());
	this.node = node;
	this.model = model;
	// recursively create
	NodeList nlist = node.getChildNodes();
	for (int i = 0; i < nlist.getLength(); i++) {
	    Node child = nlist.item(i);
	    if (child.getNodeType() == Node.ELEMENT_NODE) {
		model.insertNodeInto(new DEditNode(child,model),this,
				     this.getChildCount());
	    }
	}
    }

    /**
      * Returns the Node associated with this DEditNode.
      * @return the Node associated with this DEditNode.
      */
    public Node getNode() {
	return node;
    }

    /**
      * Returns the DefaultTreeModel associated with this DEditNode.
      * @return the DefaultTreeModel associated with this DEditNode.
      */
    public DefaultTreeModel getModel() {
	return model; // keep this function ???
    }

    /**
      * Inserts a DEditNode as a child of this node (their org.w3c.dom.Node's
      * must have the same owner Document) at the end of it's child list.
      * @param denode the DEditNode to add to the end of this Node's list
      * of children.
      */
    public void insertAtEnd(DEditNode denode) {
	model.insertNodeInto(denode,this,this.getChildCount());
	int numSpaces = getLevel() * 2;
	char[] spaces = new char[numSpaces];
	for (int i = 0; i < spaces.length; i++) {
	    spaces[i] = ' ';
	}
	String indentSpaces = new String(spaces);		
	NodeList nl = node.getChildNodes();
	Document d = node.getOwnerDocument();
	if (nl == null || nl.getLength() == 0) {
	    node.appendChild(d.createTextNode("\n"+indentSpaces));
	}
	Node n = node.getLastChild();
	node.insertBefore(d.createTextNode("\n  "+indentSpaces),n);
	node.insertBefore(denode.getNode(),n);
    }

    /**
      * Inserts a DEditNode as a child of this node (their org.w3c.dom.Node's 
      * must have the same owner Document) at a certain index in the 
      * DEditNode's children and before a certain one of it's children based 
      * on it's org.w3c.dom.Node.  These should of course match.  TODO: 
      * make function which takes only one or the other of these parameters.
      * @param denode the DEditNode to insert
      * @param index the index at which to insert it (as a TreeNode)
      * @param refChild the org.w3c.dom.Node to insert it's own 
      * org.w3c.dom.Node before
      */
    public void insertAtIndex(DEditNode denode,int index,Node refChild) {
	model.insertNodeInto(denode,this,index);
	int numSpaces = getLevel() * 2;
	char[] spaces = new char[numSpaces];
	for (int i = 0; i < spaces.length; i++) {
	    spaces[i] = ' ';
	}
	String indentSpaces = new String(spaces);		
	NodeList nl = node.getChildNodes();
	Document d = node.getOwnerDocument();
	if (nl == null || nl.getLength() == 0) {
	    node.appendChild(d.createTextNode("\n"+indentSpaces));
	}
	if (refChild == null) {
	    refChild = node.getLastChild();
	    if (refChild.getNodeType() != Node.TEXT_NODE) {
		DEditDebug.println(this,6,"tree spacing error");
		node.appendChild(denode.getNode());
		node.appendChild(d.createTextNode("\n" + indentSpaces));
		return;
	    }
	}
	Node temp = refChild.getPreviousSibling();
	while (temp != null && temp.getNodeType() == Node.TEXT_NODE) {
	    refChild = temp;
	    temp = temp.getPreviousSibling();
	}
	node.insertBefore(d.createTextNode("\n  "+indentSpaces),refChild);
	node.insertBefore(denode.getNode(),refChild);
    }

    /**
      * Removes this DEditNode from its parents, both its parent DEditNode 
      * and it's org.w3c.dom.Node's parent.
      * @return the removed node
      */
    public DEditNode removeSelf() {
	model.removeNodeFromParent(this);
	Node parent = node.getParentNode();
	if (parent != null) {
	    Node temp = node.getPreviousSibling(); 
	    while (temp != null && temp.getNodeType() == Node.TEXT_NODE) {
		parent.removeChild(temp);
		temp = node.getPreviousSibling();
	    }
	    parent.removeChild(node);
	}
	return this;
    }

    /**
      * Just detaches the org.w3c.dom.Node from its parents, for validation 
      * checking purposes.
      */
    public void removeNode() {
	Node parent = node.getParentNode();
	if (parent != null) {
	    parent.removeChild(node);
	}
    }

    /**
      * Sews up the hole from removeNode() by reinserting the passed in 
      * DEditNode's org.w3c.dom.Node as a child of this DEditNode's at the 
      * end of it's list.
      * @param denode the DEditNode whose org.w3c.dom.Node is to be 
      * reinserted
      */
    public void reinsertNodeAtEnd(DEditNode denode) {
	node.appendChild(denode.getNode());
    }

    /**
      * Sews up the hole from removeNode() by reinserting the passed in 
      * DEditNode's org.w3c.dom.Node as a child of this DEditNode's before a 
      * given Node in the list.
      * @param denode the DEditNode whose org.w3c.dom.Node is to be 
      * reinserted
      * @param next the org.w3c.dom.Node to insert denode's Node before
      */
    public void reinsertNodeBefore(DEditNode denode, Node next) {
	node.insertBefore(denode.getNode(),next);
    }

    /**
      * Sets the model for this DEditNode for insertion purposes.
      * @param model the new model
      */
    public void setModel(DefaultTreeModel model) {
	this.model = model;
    }

    /**
      * Creates a deep copy of the current DEditNode.
      * @return a deep copy of the current DEditNode
      */
    public DEditNode copy() {
	return new DEditNode(node.cloneNode(true),model);
    }
   
    /**
      * Gets the "display name" for a DEditNode, that is, that name 
      * by which it should appear in the JTree.  Usually this is 
      * equal to this.node.getNodeName(), but if this DEditNode has a child 
      * called "name", then that name is appended inside quotation marks.
      * @return the display name for this DEditNode
      */
    public String getDisplayName() {
	String ret = node.getNodeName();
	NodeList nlist = node.getChildNodes();
	String nameAppend = "";
	String valueAppend = "";
	for (int i = 0; i < nlist.getLength(); i++) {
	    Node n = nlist.item(i);
	    if (n != null && n.getNodeType() == Node.ELEMENT_NODE) { 
		if (n.getNodeName().equals("name")) {
		    NodeList childlist = n.getChildNodes();
		    for (int j = 0; j < childlist.getLength(); j++) {
			Node child = childlist.item(j);
			if (child.getNodeType() == Node.TEXT_NODE) {
			    nameAppend = "  \"" + child.getNodeValue() + "\"";
			    break;
			}
		    }
		}
		/*
		if (n.getNodeName().equals("value")) {
		    NodeList childlist = n.getChildNodes();
		    for (int j = 0; j < childlist.getLength(); j++) {
			Node child = childlist.item(j);
			if (child.getNodeType() == Node.TEXT_NODE) {
			    valueAppend = "=\"" + child.getNodeValue() + "\"";
			    break;
			}
		    }
		}
		*/
	    }
	}
	return ret + nameAppend + valueAppend;
    }

}
