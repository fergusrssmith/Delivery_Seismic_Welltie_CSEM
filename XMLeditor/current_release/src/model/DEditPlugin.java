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
import javax.swing.JMenuBar;
import javax.swing.JTree;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.JTextArea;
import javax.swing.JDesktopPane;
import java.util.Vector;
import java.awt.Component;

/**
  * Interface for adding Plugins to DEditor.
  */
public interface DEditPlugin {

    /**
      * Loads a plugin into the DEditor and its important final components.
      * @param editor the DEditor that will receive this plugin
      * @param jdp the JDesktopPane all JInternalFrame's / DEditXMLWindow's
      * will appear in
      * @param schema the DEditSchema for editor
      * @param menuBar the JMenuBar used by editor
      */
    public void load(final DEditor editor, JDesktopPane jdp, 
		     DEditSchema schema, JMenuBar menuBar);

    /**
      * Plugin potential override / add-on for editing nodes.  Returns 
      * true iff no more editing should occur after this function.  If it 
      * returns false the standard editing dialog will be displayed
      * @param selected the DEditNode that might be edited
      * @return true if editing should stop after this function
      */
    public boolean edit(DEditNode selected);

    /**
      * Returns whether or not the input node should be editable.
      * @param selected the DEditNode in question for editable / not editable
      * @param defaultEnabled true if Schema says it should be editable (i.e., 
      * has parameters or is a leaf element)
      * @return true if selected DEditNode should be editable
      */
    public boolean isEditable(DEditNode selected, boolean defaultEnabled);
    
    /**
      * Puts restrictions on what is cuttable above and beyond Schema 
      * restrictions.  Returns true when it should be cuttable if the Schema 
      * allows it.
      * @param selected the node in question for cuttable / uncuttable
      * @return true if it should be cuttable if allowed by schema
      */
    public boolean isCuttable(DEditNode selected);
    
    /**
      * Puts restrictions on what is copyable above and beyond Schema 
      * restrictions.  Returns true when it should be copyable if the Schema 
      * allows it.
      * @param selected the node in question for copyable / uncopyable
      * @return true if it should be copyable if allowed by schema
      */
    public boolean isCopyable(DEditNode selected);
    
    /**
      * Puts restrictions on where new tags may be be inserted above and 
      * beyond Schema restrictions.  The input DEditNode is the node that 
      * the new node will be inserted before.
      * @param selected the node queried to see if insert before it is 
      * ever allowed
      * @return true if given that it is valid under the Schema, a node 
      * may be inserted before this node
      */
    public boolean newTagBeforeAllowed(DEditNode selected);
    
    /**
      * Puts restrictions on where new tags may be be inserted above and 
      * beyond Schema restrictions.  The input DEditNode is the node that 
      * the new node will be inserted as a child of.
      * @param selected the node queried to see if insert as its child is 
      * ever allowed
      * @return true if given that it is valid under the Schema, a node 
      * may be inserted as a child of this node
      */
    public boolean newChildTagAllowed(DEditNode selected);
    
    /**
      * Puts restrictions on where new tags may be be inserted above and 
      * beyond Schema restrictions.  The input DEditNode is the node that 
      * the new node will be inserted after.
      * @param selected the node queried to see if insert after it is 
      * ever allowed
      * @return true if given that it is valid under the Schema, a node 
      * may be inserted after this node
      */
    public boolean newTagAfterAllowed(DEditNode selected);
    
    /**
      * Puts restrictions on what tags may be switched to other types of 
      * tags valid in the same place above.
      * @param selected the node in question for switchable / unswitchable
      * @return true if node is ever switchable
      */
    public boolean switchTagToAllowed(DEditNode selected);
    
    /**
      * Puts restrictions on whether or not a particular node should 
      * allow something to be pasted in as a child of the DEditNode passed in
      * if the Schema allows it
      * @param selected the node for which to decide if pasting a node as 
      * a child should be allowed if the Schema allows it
      * @return true if node can ever accept children to be pasted
      */
    public boolean pasteChildAllowed(DEditNode selected);
    
    /**
      * Puts restrictions on whether or not a particular node should 
      * allow something to be pasted in before the DEditNode passed in
      * if the Schema allows it
      * @param selected the node for which to decide if pasting a node before 
      * it should be allowed if the Schema allows it
      * @return true if node can ever accept a node to be pasted before it
      */
    public boolean pasteBeforeAllowed(DEditNode selected);
    
    /**
      * Puts restrictions on whether or not a particular node should 
      * allow something to be pasted in after the DEditNode passed in
      * if the Schema allows it
      * @param selected the node for which to decide if pasting a node after 
      * it should be allowed if the Schema allows it
      * @return true if node can ever accept a node to be pasted after it
      */
    public boolean pasteAfterAllowed(DEditNode selected);

    /**
      * Gets the String ID for JavaHelp that it should use to display help 
      * for the main tree.  All DEditXMLWindows query the plugin for this.
      * If you want to use default behavior, return null.
      * @return the String ID for JavaHelp to use to display help for the JTree
      */
    public String getTreeHelpID();

    /**
      * Essentially the same as the TreeCellRenderer function, except that 
      * Object value is replaced with DEditNode node (it's the same object) 
      * to save casting on the plugin side.  It also receives the renderer
      * from the DEditTreeCellRenderer's rendering, which you can return to 
      * have default behavior.
      * @param rndr the result of 
      * DEditTreeCellRenderer.getTreeCellRendererComponenent()
      * @param tree same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param node the Object called "value" from 
      * TreeCellRenderer.getTreeCellRendererComponent()
      * @param selected same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param expanded same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param leaf same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param row same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param hasFocus same as TreeCellRenderer.getTreeCellRendererComponent()
      * @return a Component with which to display DEditNode node
      */
    public Component getTreeCellRendererComponent(DEditTreeCellRenderer rndr,
						  JTree tree, 
						  DEditNode node, 
						  boolean selected, 
						  boolean expanded,
						  boolean leaf,
						  int row,
						  boolean hasFocus);
    
}
