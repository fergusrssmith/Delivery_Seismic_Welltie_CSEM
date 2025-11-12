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

import javax.xml.parsers.*;    
import org.xml.sax.SAXException;  
import org.xml.sax.SAXParseException;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.io.*;
import java.util.Vector;
import java.util.Enumeration;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.event.*;
import java.awt.event.*;
import java.awt.Dimension;
import java.lang.reflect.Method;
import org.w3c.dom.*;
import javax.help.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.net.URL;
import java.net.MalformedURLException;

/**
  * Class which represents a window for editing the XML tree but also contains 
  * it's text frame and a number of other features.  It is essentially the 
  * graphical and logical components of an opened XML file from the
  * perspective of the editor.
  */
public final class DEditXMLWindow {
    
    final private JTree tree; // the tree displaying the XML
    final private JTextField docStatus;  // status bar at the bottom of 
    // XML window
    final private JTextArea text; // text of XML to be displayed when asked
    final private JInternalFrame treeFrame; // the window in which the xml 
    // tree is displayed
    final private JInternalFrame textFrame; // the window in which the xml 
    // text is displayed
    private Document document = null;  // the DOM document for the XML
    private String openFilename = null; // name of the file that is open
    private boolean unsaved = false; // true if the file has ever been saved
    private boolean modifiedSinceSaved = false; // true if the file has been 
    // modified since it was last saved
    private DEditor owner; // the DEditor that this DEditXMLWindow is inside of
    final private DefaultTreeModel model; //TreeModel for TreeNode manipulation
    final private DEditSchema schema; // the Schema being used
    private static int unsavedIndex = 1; // used to give names to unsaved files
    private static final int UNDO_BUFFER_SIZE = 4; // number of trees in the 
    // undo buffer == max number of undo's that can be performed at any given 
    // time
    private Vector undoBuffer = new Vector(); // the buffer of root DEditNodes
    // Map IDs for JavaHelp
    final private static String TREE_ID = "tree"; 
    final private static String STATUS_ID = "status";
    final private static String TEXT_ID = "text";

    private DEditPlugin plugin = null; // the plugin being used

    /**
      * Creates a DEditXMLWindow for an existing file by opening it as for 
      * a new file then changing the relevant information.
      * @param owner the DEditor this DEditXMLWindow belongs to
      * @param jdp the JDesktopPane that this DEditXMLWindow's JInternalFrames
      * will appear in
      * @param schema the DEditSchema used for this DEditXMLWindow (and all 
      * others as well)
      * @param filename the name of the file that produced document
      * @param renderer the DEditTreeCellRenderer that will be used to 
      * render the JTree in this DEditXMLWindow
      * @param document the DOM Document currently open
      * @param mls MouseListeners to install into treeFrame
      */
    public DEditXMLWindow(final DEditor owner, JDesktopPane jdp,
			  DEditSchema schema, String filename,
			  DEditTreeCellRenderer renderer, Document document,
			  MouseListener[] mls, DEditPlugin plugin) {
	this(owner,jdp,schema,renderer,document,mls,plugin,false);	
	treeFrame.setTitle(filename);
	openFilename = filename; // reset from generic to actual file name
	unsavedIndex--; // don't waste a generic
	unsaved = false; // it has been saved
	update(false); // update, but nothing's changed
    }

    /**
      * Creates a DEditXMLWindow by setting up it's appropriate JInternalFrames
      * and installing MouseListeners and a (DEdit)TreeCellRenderer.
      * @param owner the DEditor this DEditXMLWindow belongs to
      * @param jdp the JDesktopPane that this DEditXMLWindow's JInternalFrames
      * will appear in
      * @param schema the DEditSchema used for this DEditXMLWindow (and all 
      * others as well)
      * @param renderer the DEditTreeCellRenderer that will be used to 
      * render the JTree in this DEditXMLWindow
      * @param document the DOM Document currently open
      * @param mls MouseListeners to install into treeFrame
      * @param generate true if the schema should generate all required 
      * children for the root of the tree
      */
    public DEditXMLWindow(final DEditor owner, JDesktopPane jdp,
			  DEditSchema schema, DEditTreeCellRenderer renderer,
			  Document document, MouseListener[] mls,
			  DEditPlugin plugin, boolean generate) {
	// set fields
	this.owner = owner;
	this.schema = schema;
	this.document = document;
	this.unsaved = true;
	this.tree = new JTree();
	String treeID = null;
	if (plugin != null) {
	    treeID = plugin.getTreeHelpID();
	}
	if (treeID == null) {
	    treeID = TREE_ID;
	}
	CSH.setHelpIDString(tree,treeID);
	TreeSelectionModel tsm = tree.getSelectionModel();
	// allow contiguous selection for copy and cut
	tsm.setSelectionMode(TreeSelectionModel.CONTIGUOUS_TREE_SELECTION);
	model = new DefaultTreeModel(DEditNode.EMPTY_NODE);
	tree.setEditable(false);
	JScrollPane treeView = new JScrollPane(); // holds the JTree
	text = new JTextArea(33,40);
	CSH.setHelpIDString(text,TEXT_ID);
	text.setEditable(false);
	docStatus = new JTextField();
	CSH.setHelpIDString(docStatus,STATUS_ID);
	docStatus.setEditable(false);
	docStatus.addMouseListener(new MouseAdapter() {
		public void mouseClicked(MouseEvent me) {
		    if (me.getClickCount() > 1) {
			resetDocStatus();
			JOptionPane.showMessageDialog(owner,
						      docStatus.getText());
			
		    }
		}
	    });
	openFilename = "Untitled<" + unsavedIndex + ">";
	unsavedIndex++; // increment anonymous file index
	JPanel leftPanel = new JPanel(new BorderLayout());
	leftPanel.add(treeView,BorderLayout.CENTER);
	leftPanel.add(docStatus,BorderLayout.SOUTH);
	JScrollPane textView = new JScrollPane(text);
	treeView.setPreferredSize(new Dimension(440,511));
	treeFrame = new JInternalFrame(openFilename,true,true,true,true);
	treeFrame.setDefaultCloseOperation(JInternalFrame.DO_NOTHING_ON_CLOSE);
	treeFrame.addInternalFrameListener(new InternalFrameAdapter() {
		public void internalFrameActivated(InternalFrameEvent e) {
		    owner.changeEditOptions(DEditXMLWindow.this);
		    update(false);
		}
		public void internalFrameClosing(InternalFrameEvent e) {
		    close();
		}

		public void internalFrameClosed(InternalFrameEvent e) {
		    owner.removeXMLWindow(DEditXMLWindow.this);
		    textFrame.dispose();
		    //treeFrame.dispose();
		}
	    });
	treeFrame.getContentPane().add(leftPanel);
	textFrame = new JInternalFrame("XML Text",true,true,true,true);
	textFrame.getContentPane().add(textView);
	textFrame.setLocation(450,0);
	textFrame.setDefaultCloseOperation(JInternalFrame.HIDE_ON_CLOSE);
	tree.addTreeSelectionListener(new TreeSelectionListener() {
		public void valueChanged(TreeSelectionEvent e) {
		    owner.changeEditOptions(DEditXMLWindow.this);
		}
	    });
	treeView.setViewportView(tree);
	jdp.add(textFrame);
	jdp.add(treeFrame);
	textFrame.pack();
	//textFrame.show();
	treeFrame.pack();
	treeFrame.show();
	for (int i = 0; i < mls.length; i++) {
	    tree.addMouseListener(mls[i]);
	}
	DEditNode rootNode = new DEditNode(document.getDocumentElement(),
					   model);
	model.setRoot(rootNode);
	tree.setModel(model);
	if (generate) {
	    schema.generateRequiredChildren(rootNode,model);
	}
	model.reload();
	tree.setCellRenderer(renderer);	
	update(false);
    }

    /**
      * Resets the text in the document status text field.
      */
    private void resetDocStatus() {
	if (!schema.isReady()) {
	    docStatus.setText("No schema loaded.");
	} else if (openFilename == null) { // this shouldn't ever happen
	    docStatus.setText("No XML document opened.");
	} else {
	    isValid();
	    docStatus.setText(schema.getVerifyString());
	    docStatus.setCaretPosition(0);
	}
    }

    /**
      * Cuts the passed in subtrees and returns a copy of the cut trees.
      * @param toCut the subtrees to cut from the tree
      * @return the passed in array of subtrees which have now been removed
      * from the JTree
      */
    public DEditNode[] cutSubTree(DEditNode[] toCut) {
	addToUndoBuffer();
	for (int i = 0; i < toCut.length; i++) {
	    toCut[i].removeSelf();
	}
	update(true);
	return toCut;
    }

    /**
      * Undo back to the most recent tree stored in the undo buffer.
      */
    public void undo() {
	if (undoBuffer.size() == 0) {
	    return;
	}
	DEditNode oldRoot = (DEditNode)undoBuffer.remove(undoBuffer.size()-1);
	Node domRoot = document.getDocumentElement();
	document.replaceChild(oldRoot.getNode(),domRoot);
	model.setRoot(oldRoot);
	update(true);
    }

    /**
      * Add the current tree to the undo buffer (do this right before
      * cutting, pasting, inserting, etc.).
      */
    private void addToUndoBuffer() {
	if (undoBuffer.size() == UNDO_BUFFER_SIZE) {
	    undoBuffer.remove(0);
	}
	undoBuffer.add(((DEditNode)model.getRoot()).copy());
    }

    /**
      * Used to determine if undo button or menu item should be enabled.
      * @return true if this DEditXMLWindow can "undo."
      */
    public boolean canUndo() {
	return (undoBuffer.size() != 0);
    }

    /**
      * Pastes an altered version of the copied DEditNodes into the tree 
      * before the currently selected node.
      * @param copied the copied DEditNodes to create DEditNodes suitable 
      * for this JTree from and then paste in
      */
    public void pasteSibTreeBefore(DEditNode[] copied) {
	DEditNode denode = getSelected();
	if (denode == null) {
	    return;
	}
	addToUndoBuffer();
	pasteSibTreeBefore(denode,copied);
    }

    /**
      * Pastes an altered version of the copied DEditNodes into the tree 
      * before the input DEditNode denode.
      * @param denode the DEditNode to paste before
      * @param copied the copied DEditNodes to create DEditNodes suitable for 
      * this JTree from and then paste in
      */
    private DEditNode[] pasteSibTreeBefore(DEditNode denode,
					  DEditNode[] copied) {
	DEditNode parent = (DEditNode)(denode.getParent());
	if (parent == null) {
	    return null; // shouldn't happen
	}
	DEditNode[] inserted = new DEditNode[copied.length];
	for (int i = 0; i < copied.length; i++) {
	    Node imported = document.importNode(copied[i].getNode(),true);
	    DEditNode toInsert = new DEditNode(imported,model);	
	    int index = model.getIndexOfChild(parent,denode);
	    parent.insertAtIndex(toInsert,index,denode.getNode());
	    model.nodeChanged(toInsert);
	    inserted[i] = toInsert;
	}
	update(true);
	return inserted;
    }

    /**
      * Pastes an altered version of the copied DEditNode into the tree after 
      * the currently selected node.
      * @param copied the copied DEditNode to create a DEditNode suitable for 
      * this JTree from and then paste in
      */
    public void pasteSibTreeAfter(DEditNode[] copied) {
	DEditNode denode = getSelected();
	if (denode == null) {
	    return;
	}
	addToUndoBuffer();
	pasteSibTreeAfter(denode,copied);
    }
    
    /**
      * Pastes an altered version of the copied DEditNode into the tree after 
      * the input DEditNode denode.
      * @param denode the DEditNode to paste after
      * @param copied the copied DEditNode to create a DEditNode suitable for 
      * this JTree from and then paste in
      */
    private DEditNode[] pasteSibTreeAfter(DEditNode denode,
					 DEditNode[] copied) {
	DEditNode parent = (DEditNode)(denode.getParent());
	if (parent == null) {
	    return null; // shouldn't happen
	}
	DEditNode[] inserted = new DEditNode[copied.length];
	int index = model.getIndexOfChild(parent,denode);
	Node refChild = denode.getNode().getNextSibling();
	for (int i = 0; i < copied.length; i++) {
	    Node imported = document.importNode(copied[i].getNode(),true);
	    DEditNode toInsert = new DEditNode(imported,model);	
	    parent.insertAtIndex(toInsert,index+1+i,
				 refChild);
	    model.nodeChanged(toInsert);
	    inserted[i] = toInsert;
	}
	update(true);
	return inserted;
    }

    /**
      * Gets Action's which represent creating and inserting new sub trees 
      * before / above the currently selected subtree.
      * @return an array of Action's which insert DEditNodes into the JTree
      */
    public Action[] getNewTagBeforeActions() {
	DEditNode selected = getSelected();
	if (selected == null) {
	    return null;
	}
	return getNewTagBeforeActions(selected);
    }
    
    /**
      * Gets Action's which represent creating and inserting new sub trees 
      * before / above the input subtree.
      * @param selected the DEditNode to insert before / above
      * @return an array of Action's which insert DEditNodes into the JTree
      */
    public Action[] getNewTagBeforeActions(DEditNode selected) {
	String[] names = schema.getNewTagBeforeNames(selected);
	if (names != null) {
	    DEditNode parent = (DEditNode)selected.getParent();
	    int index = parent.getIndex(selected);
	    return createInsertActions(names,parent,index);
	}
	return null;
    }

    /**
      * Gets Action's which represent creating and inserting new sub trees 
      * after / below the currently selected subtree.
      * @return an array of Action's which insert DEditNodes into the JTree
      */
    public Action[] getNewTagAfterActions() {
	DEditNode selected = getSelected();
	if (selected == null) {
	    return null;
	}
	return getNewTagAfterActions(selected);
    }
    
    /**
      * Gets Action's which represent creating and inserting new sub trees 
      * after / below the input subtree.
      * @param selected the DEditNode to insert after / below
      * @return an array of Action's which insert DEditNodes into the JTree
      */
    public Action[] getNewTagAfterActions(DEditNode selected) {
	String[] names = schema.getNewTagAfterNames(selected);
	if (names != null) {
	    DEditNode parent = (DEditNode)selected.getParent();
	    int index = parent.getIndex(selected);
	    return createInsertActions(names,parent,index+1);
	}
	return null;
    }

    /**
      * Gets Action's which will create and insert a chosen subtree into 
      * the JTree as a child of the currently selected node on the JTree.
      * @return an array of Action's which create and insert subtrees into the
      * currently selected node on the JTree.
      */
    public Action[] getNewChildTagActions() {
	DEditNode parent = getSelected();
	if (parent == null) {
	    return null;
	}
	return getNewChildTagActions(parent);
    }

    /**
      * Gets Action's which will create and insert a chosen subtree into 
      * the JTree as a child of the input node on the JTree.
      * @param parent the DEditNode which the returned Action's will insert 
      * subtrees into as children
      * @return an array of Action's which create and insert subtrees into the
      * input node on the JTree.
      */
    public Action[] getNewChildTagActions(final DEditNode parent) {
	String[] names = schema.getNewChildTagNames(parent);
	if (names == null || names.length == 0) {
	    return null;
	}
	Action[] ret = new Action[names.length];
	for (int i = 0; i < ret.length; i++) {
	    final String name = names[i];
	    ret[i] = new AbstractAction(name) {
		    public void actionPerformed(ActionEvent ae) {
			createNewChildTag(parent,name);
			update(true);
		    } 

		    public String toString() {
			return (String)getValue(Action.NAME);
		    }
		};
	}
	return ret;
    }    
    
    /**
      * Creates a new child element tag in the XML by delegating the task to
      * the DEditSchema.
      * @param parent the parent to create a new child node for
      * @param newTagName the name of the new node to be created
      */
    private void createNewChildTag(DEditNode parent, String newTagName) {
	schema.createNewChildTag(parent,newTagName,model);
    }

    /**
      * Gets Action's which will switch the currently selected DEditNode
      * to a different DEditNode.
      * @return Action's which will switch the currently selected DEditNode
      * to a different DEditNode
      */
    public Action[] getSwitchTagToActions() {
	DEditNode current = getSelected();
	if (current == null) {
	    return null;
	}
	return getSwitchTagToActions(current);
    }
    
     /**
      * Gets Action's which will switch the input DEditNode
      * to a different DEditNode.
      * @param selected the returned Action's will change this DEditNode 
      * into a different DEditNode
      * @return Action's which will switch the input DEditNode
      * to a different DEditNode
      */
    public Action[] getSwitchTagToActions(final DEditNode selected) {
	String[] nodeNames = schema.getSwitchTagToNames(selected,model);
	if (nodeNames == null || nodeNames.length == 0) {
	    return null;
	}
	Node next = null;
	final DEditNode parent = (DEditNode)selected.getParent();
	if (parent == null) {
	    return null;
	}
	final int index = parent.getIndex(selected);
	if (index+1 < parent.getChildCount()) {
	    DEditNode nextDEN = (DEditNode)parent.getChildAt(index+1);
	    next = nextDEN.getNode();
	}
	Action[] ret = new Action[nodeNames.length];
	final Node insertBefore = next;
	for (int i = 0; i < ret.length; i++) {
	    final String name = nodeNames[i];
	    Action temp = new AbstractAction(name) {
		    public void actionPerformed(ActionEvent ae) {
			DEditNode addable = makeDEditNode(name);
			addToUndoBuffer();
			selected.removeSelf();
			parent.insertAtIndex(addable,index,insertBefore);
			schema.generateRequiredChildren(addable,model);
		    }
		    public String toString() {
			return (String)getValue(Action.NAME);
		    }
		};
	    ret[i] = temp;
	}
	return ret;
    }

    /**
      * Gets Action's which paste the currently copied DEditNodes in various 
      * locations as children of the currently selected node
      * @param copied the currently copied DEditNodes
      * @return Action's which paste the currently copied DEditNodes in 
      * various locations as children of the currently selected node
      */
    public Action[] getPasteChildActions(DEditNode[] copied) {
	DEditNode selected = getSelected();
	if (selected == null) {
	    return null;
	}
	return getPasteChildActions(selected,copied);
    }

    /**
      * Returns an array of Actions to be used to create a "where you 
      * can paste" submenu.
      * @param denode the node being checked to see where the copied
      * DEditNode can be pasted as one of it's children
      * @param copied the currently copied DEditNode's
      * @return an array of Actions detailing paste-able locations.
      */
    public Action[] getPasteChildActions(DEditNode denode,
					 final DEditNode[] copied) {
	if (denode == null || copied == null) {
	    // shouldn't happen, paste should be 
	    //disallowed in this case
	    return null;
	}
	if (denode.isLeaf()) {
	    return null;
	}
	boolean validBefore = schema.isReady() && isValid();
	Vector toBeAdded = new Vector();
	Enumeration kids = denode.children();
	while (kids.hasMoreElements()) {
	    final DEditNode child = (DEditNode)kids.nextElement();
	    int index = getIndexOfChild(denode,child);
	    if (!validBefore || pasteBeforeCheck(child,copied)) {
		Action act = new AbstractAction("Paste before " 
						   + child.getDisplayName()) {
			public void actionPerformed(ActionEvent ae) {
			    addToUndoBuffer();
			    pasteSibTreeBefore(child,copied);
			}
			public String toString() {
			    return (String)getValue(Action.NAME);
			}
		    };
		toBeAdded.add(act);
	    }	
	}
	final DEditNode lastChild = (DEditNode)denode.getLastChild();
	if (lastChild != null) {
	    int index = getIndexOfChild(denode,lastChild);
	    
	    if (!validBefore 
		|| pasteAfterCheck(lastChild,copied)) {
		Action act = new AbstractAction("Paste at end") {
			public void actionPerformed(ActionEvent ae) {
			    addToUndoBuffer();
			    pasteSibTreeAfter(lastChild,copied);
			}
			public String toString() {
			    return (String)getValue(Action.NAME);
			}
		    };
		toBeAdded.add(act);
	    }	
	}
	Action[] ret = new Action[toBeAdded.size()];
	toBeAdded.toArray(ret);
	return ret;
    }

    /**
      * Creates Action's which will be for inserting the DEditNodes in nodes 
      * into DEditNode parent at the index provided.  
      * (createInsertActions(...)[i] inserts nodes[i])
      * @param nodes the DEditNode's the returned Action's will insert
      * @param parent the DEditNode into which the returned Action's will 
      * insert DEditNode's as children
      * @param index the index of parent's children at which to insert
      */
    private Action[] createInsertActions(final String[] nodeNames, 
					 final DEditNode parent,
					 final int index) {
	Action[] ret = new Action[nodeNames.length];
	int numSpaces = (parent.getLevel()) * 2;
	char[] spaces = new char[numSpaces];
	for (int i = 0; i < spaces.length; i++) {
	    spaces[i] = ' ';
	}
	final String indentSpaces = new String(spaces);		
	for (int i = 0; i < ret.length; i++) {
	    final String name = nodeNames[i];
	    ret[i] = new AbstractAction(name) {
		    public void actionPerformed(ActionEvent ae) {
			addToUndoBuffer();
			DEditNode nextSibling = null;
			Node parNode = parent.getNode();
			Document d = parNode.getOwnerDocument();
			Node nextNode = null;
			DEditDebug.println(this,6,"i " + index + " cc "
					   + parent.getChildCount());
			
			if (index < parent.getChildCount()) {
			    nextSibling = (DEditNode)parent.getChildAt(index);
			}
			if (nextSibling != null) {
			    nextNode = nextSibling.getNode();
			}
			DEditNode node = makeDEditNode(name);
			parent.insertAtIndex(node,index,nextNode);
			schema.generateRequiredChildren(node,model);
			update(true);
			owner.changeEditOptions(DEditXMLWindow.this);
		    }
		    
		    public String toString() {
			return (String)getValue(Action.NAME);
		    }
		};
	}
	return ret;
    }

    /**
      * Creates a new DEditNode (useful because DEditNode always wants a
      * DefaultTreeModel and the correct one is stored in this class, and 
      * it must have a Node created from the proper Document, also stored 
      * in this class).
      * @param nodeName the name of the Node to create and put in the 
      * returned DEditNode
      * @return a DEditNode created with the provided name
      */
    public DEditNode makeDEditNode(String nodeName) {
	return new DEditNode(document.createElement(nodeName),model);
    }

    /**
      * Creates a new DEditNode (useful because DEditNode always wants a
      * DefaultTreeModel and the correct one is stored in this class).
      * @param node the Node to put in the new DEditNode
      * @return a DEditNode created with the provided Node
      */
    public DEditNode makeDEditNode(Node node) {
	return new DEditNode(node,model);
    }

    /**
      * A convenience method that asks the schema to generate the children 
      * with the right model, and then tells the model that the node who just 
      * had children created for it changed.
      * @param node the DEditNode to generate children for
      */
    public void generateRequiredChildren(DEditNode node) {
	schema.generateRequiredChildren(node,model);
	model.nodeChanged(node);
    }
    
    /**
      * Limited accessor to TreeModel knowledge which returns the root 
      * DEditNode in this tree.
      * @return the root DEditNode in this tree
      */
    public DEditNode getRoot() {
	return (DEditNode)(model.getRoot());
    }

    /**
      * Tells the TreeModel that a node in the tree has changed
      * @param node the node that has changed
      */
    public void nodeChanged(DEditNode node) {
	model.nodeChanged(node);
    }

    /**
      * Tells the TreeModel that a node in the tree has drastically changed
      * @param node the node that has changed
      */
    public void nodeStructureChanged(DEditNode node) {
	model.nodeStructureChanged(node);
    }

    /**
      * Returns the index of a child w/respect to it's parent in the tree. 
      * Again this is just the DEditXMLWindow asking the model for information 
      * to shield the model from exposure to other classes.
      * @param parent the parent whose child's index we are looking for
      * @param child the node whose index as a child of parent we are 
      * searching for
      * @return the index of child in the children of parent
      */
    public int getIndexOfChild(DEditNode parent, DEditNode child) {
	return model.getIndexOfChild(parent,child);
    }

    /**
      * Asks the DOM Document of this DEditXMLWindow to create an attribute 
      * typed Node with given name.
      * @param attName the name to give to the attribute Node
      * @return an attribute typed Node with the given name
      */
    public Attr createAttribute(String attName) {
	return document.createAttribute(attName);
    }

    /**
      * Returns the first DEditNode currently selected in the tree.
      * @return the first DEditNode currently selected in the tree
      */
    public DEditNode getSelected() {
	DEditNode[] selected = getMultipleSelected();
	if (selected == null || selected.length != 1) {
	    return null;
	}
	return selected[0];
    }

    /**
      * Returns the DEditNodes currently selected in the tree.
      * @return the DEditNodes currently selected in the tree
      */
    public DEditNode[] getMultipleSelected() {
	TreePath[] tp = tree.getSelectionPaths();
	if (tp == null || tp.length == 0) {
	    return null;
	}
	DEditNode[] ret = new DEditNode[tp.length];
	for (int i = 0; i < ret.length; i++) {
	    ret[i] = (DEditNode)tp[i].getLastPathComponent();
	}
	return ret;
    }

    /**
      * Checks to see if it would be legal to paste DEditNode copied before 
      * DeditNode selected.
      * @param selected the DEditNode to see if pasting before it with 
      * the copied DEditNode is legal
      * @param copied the DEditNode being pasted into the tree
      * @return true if it's legal to paste copied before selected
      */
    public boolean pasteBeforeCheck(DEditNode selected, DEditNode[] copied) {
	if (selected.isRoot()) {
	    return false;
	}
	boolean modified = modifiedSinceSaved;
	DEditNode[] newNodes = pasteSibTreeBefore(selected,copied);
	boolean ret = isValid();
	for (int i = 0; i < newNodes.length; i++) {
	    newNodes[i].removeSelf();
	}
	modifiedSinceSaved = modified;
	return ret;
    }
    
    /**
      * Checks to see if it would be legal to paste DEditNode copied after 
      * DeditNode selected.
      * @param selected the DEditNode to see if pasting after it with 
      * the copied DEditNode is legal
      * @param copied the DEditNode being pasted into the tree
      * @return true if it's legal to paste copied after selected
      */
    public boolean pasteAfterCheck(DEditNode selected, DEditNode[] copied) {
	if (selected.isRoot()) {
	    return false;
	}
	boolean modified = modifiedSinceSaved;
	DEditNode[] newNodes = pasteSibTreeAfter(selected,copied);
	boolean ret = isValid();
	for (int i = 0; i < newNodes.length; i++) {
	    newNodes[i].removeSelf();
	}
	modifiedSinceSaved = modified;
	return ret;
    }

    /**
      * Checks to see if XML document is valid which sets String in the 
      * DEditSchema that may be used to get more information
      * @return true if XML is valid given schema
      */
    public boolean isValid() {
	return schema.isValid(document);
    }

    /**
      * Updates some features of the DEditXMLWindow.
      */
    public void update(boolean changed) {
	modifiedSinceSaved = (modifiedSinceSaved || changed);
	resetDocStatus(); // what else????
	text.setText(toXMLString());
    }

    /**
      * Returns true if the file represented by the DEditXMLWindow is unsaved,
      * meaning "Save" functions as "Save As...".
      * @return true if the file represented by the DEditXMLWindow is unsaved
      */
    public boolean unsaved() {
	return unsaved;
    }

    /**
      * Returns true if the XML Document has been modified since it was last 
      * saved, which means before closing it should prompt for saving. 
      * @return true if the XML Document has been modified since it was last 
      * saved
      */
    public boolean modifiedSinceSaved() {
	return modifiedSinceSaved;
    }

    /**
      * Saves the file currently being edited.
      * @param filename the name of the file to save.
      */
    public void saveFile(String filename) {
	if (filename == null) {
	    return;
	}
	if (owner.isApplet()) { 
	    owner.writeRemoteFile(filename,toXMLString());
	    System.out.println("wrote file " + filename);
	} else {
	    try {
		FileWriter fw = new FileWriter(filename);
		fw.write(toXMLString());
		fw.close();
	    } catch (IOException ioe) {
		JOptionPane.showMessageDialog(owner,"Couldn't save file: " 
					      + ioe.toString());
		return;
	    }
	}
	modifiedSinceSaved = false;
	openFilename = filename;
	treeFrame.setTitle(openFilename);
    }

    /**
      * Converts the DOM document into an String for display or 
      * printing.  Uses a workaround hand-made function for the applet version.
      * @return the DOM document as a String
      */
    public String toXMLString() {
	if (owner.isApplet()) {
	    return domWalk(document,0);
	} else {
	    return toXMLString2();
	}
    }
    
    /**
      * This method walks the DOM tree and converts it into a string for use
      * in display or printing to file.  This is a modified form of a simple 
      * function which printed to the terminal that I found (w/no copyright)
      * on a DOM tutorial online.
      * @param node the current node in our recursive descent
      * @param indent the indent to use on this recursion
      * @return a string representation of everything under node plus what 
      * came before */
    private String domWalk(Node node,int length) {
        String ret = "";
	DEditNode selTreeNode = getSelected();
	Node selected = null;
	if (selTreeNode != null) {
	    selected = selTreeNode.getNode();
	} else {
	    //caretPos = 0;
	}
	if (node == selected) {
	    //caretPos = length;
	    //System.out.println("caretPos = " + length);
	}
	int type = node.getNodeType();
        
	if (type == Node.DOCUMENT_NODE) { // this first string maybe should be
	    // formatted differently !!!
	    ret += "<?xml version=\"1.0\" encoding=\"" +
		"UTF-8" + "\"?>";                   
	}//end of document
	else if (type == Node.ELEMENT_NODE) {
	    ret += '<' + encode(node.getNodeName());
	    NamedNodeMap nnm = node.getAttributes();
	    if(nnm != null ) {
		int len = nnm.getLength() ;
		Attr attr;
		for (int i = 0; i < len; i++) {
		    attr = (Attr)nnm.item(i);
		    ret += ' ' + encode(attr.getNodeName()) + "=\""
			+ encode(attr.getNodeValue()) +  '"';
		}
	    }
	    ret += '>';
	}//end of element
	else if (type == Node.ENTITY_REFERENCE_NODE) {		    
	    ret += '&' + node.getNodeName() + ';';		    
	}//end of entity
        else if (type == Node.CDATA_SECTION_NODE) {
	    ret += "<![CDATA[" + node.getNodeValue() + "]]>";
	} else if (type == Node.TEXT_NODE) {
	    ret += encode(node.getNodeValue());
	} else if (type == Node.PROCESSING_INSTRUCTION_NODE) {
	    ret += "<?" + node.getNodeName();
	    String data = node.getNodeValue();
	    if (data != null && data.length() > 0) {
		ret += ' ' + data;
	    }
	    ret += "?>";		    
	} else {
	    System.err.println("Invalid type of node encountered while "
			       +"traversing DOM."); // bail out of function ???
	}
	
	//recurse
	for(Node child = node.getFirstChild(); child != null;
	    child = child.getNextSibling()) {
		ret += domWalk(child,length + ret.length());
	}
	
	//without this the ending tags will miss
	if (type == Node.ELEMENT_NODE) {
	    ret +="</" + encode(node.getNodeName()) + ">";
	}
	return ret;
    }

    /**
      * Helper for applet-induced hand-made transition from DOM to String.
      * Turns all standard escape characters for XML into the appropriate 
      * sequences.
      * @param input the string to change and return
      * @return input with special characters escaped to &__; sequences
      */
    private static String encode(String input) {
	input = input.replaceAll("[&]","&amp;");
	input = input.replaceAll("[<]","&lt;");
	input = input.replaceAll("[>]","&gt;");
	input = input.replaceAll("[\"]","&quot;");
	input = input.replaceAll("[']","&apos;");
	return input;
    }

    /**
      * Turns standard XML escapes back into proper characters.
      * @param input the input String to transform.
      * @return the transformed String
      */
    public static String decode(String input) {
	input = input.replaceAll("[&]lt;","<");
	input = input.replaceAll("[&]gt;",">");
	input = input.replaceAll("[&]quot;","\"");
	input = input.replaceAll("[&]apos;","'");
	input = input.replaceAll("[&]amp;","&");
	return input;
    }

    /**
      * Converts the DOM document into an String for display or 
      * printing.
      * @return the DOM document as a String
      */
    public String toXMLString2() {
	Node match = null;
	DEditNode selTreeNode = getSelected();
	if (selTreeNode != null) {
	    match = selTreeNode.getNode();
	}
	try {
	    TransformerFactory tFactory = TransformerFactory.newInstance();
	    Transformer tranformer = tFactory.newTransformer();
	    DOMSource source = new DOMSource(document);
	    ByteArrayOutputStream baos = new ByteArrayOutputStream();
	    StreamResult result = new StreamResult(baos);
	    tranformer.transform(source,result);
	    return baos.toString();
	} catch (TransformerConfigurationException tce) {
	    JOptionPane.showMessageDialog(owner,"Couldn't convert to XML: "
					  + tce.toString());
	} catch (TransformerException te) {
	    JOptionPane.showMessageDialog(owner,"Couldn't convert to XML: "
					  + te.toString());
	}
	return null;
    }
    
    /**
      * Instructs this DEditXMLWindow to display it's JInternalFrame which 
      * contains the text of the XML that would be generated by the tree.
      */
    public void showTextFrame() {
	textFrame.show();
    }

    /**
      * Gets the name of the file that is currently opened, or null for an 
      * unsaved file.
      * @return the name of the file that is currently opened, null for unsaved
      */
    public String getOpenFilename() {
	return openFilename;
    }
    
    /**
      * Returns true if this DEditXMLWindow's JInternalFrame is displayed. 
      * Might change to include textFrame and/or ask Plugin???
      * @return true if this DEditXMLWindow's JInternalFrame is displayed
      */
    public boolean isSelected() {
	return treeFrame.isSelected();
    }

    /**
      * Perform closing operations on this DEditXMLWindow and then close it.
      * @return true if the window was actually closed
      */
    public boolean close() {
	if (modifiedSinceSaved) {
	    String currentName = openFilename;
	    String msg = "File " + openFilename + " changed but not saved.  "
		+"Save before closing?";
	    int choice = JOptionPane.showConfirmDialog(owner,msg,
						       "File Not Saved",
						       JOptionPane.YES_NO_CANCEL_OPTION);
	    if (choice == JOptionPane.YES_OPTION) {
		if (unsaved) {
		    openFilename = owner.getSaveFile();
		}
		if (openFilename == null) {
		    openFilename = currentName;
		    return false;
		}
		saveFile(openFilename);
	    }
	    if (choice == JOptionPane.CANCEL_OPTION) {
		return false;
	    }
	} 
	treeFrame.dispose();
	return true;
    }

    /**
      * Disposes the DEditXMLWindow.  Needs only dispose the treeFrame 
      * because on close it disposes the other frames.
      */
    public void dispose() {
	treeFrame.dispose();
    }

}
