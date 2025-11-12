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
import java.awt.Image;
import java.awt.Insets;
import java.io.*;
import java.util.Vector;
import java.util.Enumeration;
import java.util.Hashtable;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.event.*;
import java.awt.event.*;
import java.awt.Dimension;
import java.lang.reflect.Method;
import org.w3c.dom.*;
import javax.help.*;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import org.xml.sax.InputSource;

/**
  * The major class for the editor.  Sets up all GUI components and opens 
  * documents.
  */
public final class DEditor extends JFrame {
   
    final private static String GPLlicenceText=
	  "This program is free software; you can redistribute it and/or modify\n"+
	  "it under the terms of the GNU General Public License Version 2 as\n"+
	  "published by the Free Software Foundation.\n\n"+
	  "This program is distributed in the hope that it will be useful, but\n"+
	  "WITHOUT ANY WARRANTY; without even the implied warranty of\n"+
	  "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n"+
	  "General Public License for more details.\n\n"+
	  "You should have received a copy of the GNU General Public License\n"+
	  "along with this program; if not, write to the Free Software\n"+
	  "Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA\n"+
	  "02110-1301, USA or visit the link http://www.gnu.org/licenses/gpl.txt.\n\n"+
	  "PLEASE NOTE: One or more BHP Billiton patents or patent applications may\n"+
	  "cover this software or its use. BHP Billiton hereby grants a limited,\n"+
	  "personal, irrevocable, perpetual, royalty-free license under such\n"+
	  "patents but only for things you are allowed to do under the terms of\n"+
	  "the General Public License.  All other rights are expressly reserved.";

    final private static String releaseDate = "09/05/03";
    // this is the cutoff for how many items go in a sub menu
    final private static int DEFAULT_CUTOFF = 20;
    final private JDesktopPane jdp;
    final private JMenuBar menuBar;
    private File curXMLDir; // the directory to start looking for xml files in
    private File curXSDDir; // the directory to start looking for xsd files in
    private DEditNode[] copied = null; // the currently copied nodes
    private DEditSchema schema;
    private final JPopupMenu popup; // the right click popup menu
    private final JToolBar toolbar; // the toolbar with all the buttons
    private final static int BUTTON_ICON_SIZE = 33; // size of icon buttons
    private final Dimension buttonDim = new Dimension(BUTTON_ICON_SIZE,
						      BUTTON_ICON_SIZE);
    private final static Insets NO_MARGIN_INSETS = new Insets(0,0,0,0);
    private final JMenu newChildTagMenu; // all the sub menus ... 
    private final JMenu newChildTagPopupMenu; // rearrange for space?
    private final JMenu newTagBeforeMenu;
    private final JMenu newTagBeforePopupMenu;
    private final JMenu newTagAfterMenu;
    private final JMenu newTagAfterPopupMenu;
    private final JMenu switchTagToMenu;
    private final JMenu switchTagToPopupMenu;
    private final JMenu pasteChildMenu;
    private final JMenu pasteChildPopupMenu;
    private final EditAction editValAction; // all the actions
    private final EditAction copyAction;
    private final EditAction cutAction;
    private final EditAction pasteBeforeAction;
    private final EditAction pasteAfterAction;
    private final EditAction undoAction;
    private final Action nodeHelpAction;
    private final Action saveAction;
    private final JMenuItem verifyItem; // it's listener verifies the xml 
    private static String openErrString = null; // stores error messages on
    // open file commands
    private JMenuItem saveAsItem; // it's listener does "Save As..." on the xml
    private boolean allowNewFile; // true if making a new file is allowed
    private boolean allowOpenFile; // true if opening files is allowed
    private boolean allowLoadSchema; // true if opening schemas is allowed
    private DEditLeafDialog dialog; // dialog for editing values
    private DEditTreeCellRenderer renderer; // renders JTrees, gets handed off 
    // to every DEditXMLWindow that's created
    private DEditPlugin plugin = null; // the currently loaded plugin
    private MouseListener ml; // double click -> edit value dialog listener
    private MouseListener popupListener = new PopupListener(); // right click 
    // -> popup edit menu listener

    final private static int MAX_WINDOWS = 10; // max number of open xml files
    private Vector windows = new Vector(); // holds open DEditXMLWindow's
    private JApplet applet = null; // holds the applet if there is one
    private boolean isApplet = false; // state for if DEditor is running as
    // an applet
    private Hashtable properties = new Hashtable();
    
    // HelpSet, HelpBroker, and Map ID's for JavaHelp
    private HelpSet deditHelpSet;
    private HelpBroker deditHelpBroker;
    private static final String NEW_XML_ID = "new_xml";
    private static final String OPEN_XML_ID = "open_xml";
    private static final String SAVE_XML_ID = "save_xml";
    private static final String INSERT_ABOVE_ID = "insert_above";
    private static final String INSERT_CHILD_ID = "insert_child";
    private static final String INSERT_BELOW_ID = "insert_below";
    private static final String SWITCH_ID = "switch";
    private static final String EDIT_VALUE_ID = "edit_value";
    private static final String COPY_ID = "copy";
    private static final String CUT_ID = "cut";
    private static final String PASTE_ABOVE_ID = "paste_above";
    private static final String PASTE_CHILD_ID = "paste_child";
    private static final String PASTE_BELOW_ID = "paste_below";

    /**
      * Constructor for the top-level class. Sets up all GUI components.
      * @param curDirString the current directory the user is running in
      * @param configFile name of an XML file to use for configuration of 
      * DEditor
      * @param filename name of an XML file to open if non-null
      * @param schemaFilename name of an XML Schema file to open
      * @param allowNewFile true if user will be allowed to makenew files
      * @param allowOpenFile true if user will be allowed to open files
      * @param allowLoadSchema true if user will be allowed to load schemas
      * @param relaxed true if the Editor should run in relaxed verification
      * mode (see ???README)
      * @param applet the applet (if any) that this DEditor is running from
      */ 
    public DEditor(String curDirString, String configFile,
		   String filename, String schemaFilename,
		   boolean allowNewFile, boolean allowOpenFile,
		   boolean allowLoadSchema, boolean relaxed,
		   DEditPlugin plugin, JApplet applet) {
	super("XML EDITOR");
	this.applet = applet;
	this.isApplet = (applet != null);
	// initialize Schema
	this.schema = new DEditSchema(this,relaxed);
	// initialize working directories
	this.curXMLDir = new File(curDirString);
	this.curXSDDir = this.curXMLDir;
	// initialize availabilities
	this.allowNewFile = allowNewFile;
	this.allowOpenFile = allowOpenFile;
	this.allowLoadSchema = allowLoadSchema;
	// set up help
	String helpSetName = "help/DEditHelp.hs";
	try {
	    ClassLoader cl = this.getClass().getClassLoader();
	    URL hsURL = HelpSet.findHelpSet(cl, helpSetName);
	    deditHelpSet = new HelpSet(cl, hsURL);
	} catch (Exception e) {
	    System.err.println("Helpset " + helpSetName + " not found.");
	    // now what ???
	}
	// load up configuration file
	Document configDoc = null;
	if (configFile != null) {
	    configDoc = createDOMforXML(configFile);
	}
	String[][] leafPairs = null;
	String[][] rendererTextPairs = null; 
	String[][] rendererIconPairs = null;
	if (configDoc != null) {
	    leafPairs = getLeafPairs(configDoc);
	    rendererTextPairs = getRendererTextPairs(configDoc);
	    rendererIconPairs = getRendererIconPairs(configDoc);
	} else {
	    System.err.println(openErrString);
	}
	// initialize desktop pane, edit value dialog, TreeCellRenderer,
	// and the main frame
	jdp = new JDesktopPane();
	this.dialog = new DEditLeafDialog(this.curXMLDir,leafPairs,isApplet);
	renderer = new DEditTreeCellRenderer(rendererTextPairs,
					     rendererIconPairs);
	this.getContentPane().setLayout(new BorderLayout());
	jdp.setPreferredSize(new Dimension(900,600));
	this.getContentPane().add(jdp,BorderLayout.CENTER);
	/**** Menu Set-up ****/
	menuBar = new JMenuBar();
	// ----- File Menu -----
	JMenu fileMenu = new JMenu("File");
	// New
	JMenuItem newItem = new JMenuItem("New");
	newItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    makeNewFile();
		}
	    });
	newItem.setEnabled(allowNewFile);
	fileMenu.add(newItem);
	// Open
	JMenuItem openItem = new JMenuItem("Open XML File");
	openItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    openXML();
		}
	    });
	openItem.setEnabled(allowOpenFile);
	fileMenu.add(openItem);
	// Load Schema
	JMenuItem loadSchemaItem = new JMenuItem("Load Schema File");
	loadSchemaItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    // have user select .xsd Schema file
		    String filename;
		    if (isApplet()) {
			String msg = "Specify remote URL for Schema";
			RemoteFileChooser rfc 
			    = new RemoteFileChooser(DEditor.this,
						    "Open Configuration File",
						    "loadSchemaURL","Open");
			int reFchooser = rfc.showDialog();
			if (reFchooser == JFileChooser.APPROVE_OPTION) {
			    String pathString = rfc.getPath();
			    String name = rfc.getName();
			    String absoluteName = pathString + "/" + name;
			    System.out.println("(" + name + ")OOpenSession " + absoluteName);
			    filename = absoluteName;
			} else {
			    return;
			}
		    } else {
			JFileChooser chooser = new JFileChooser(curXSDDir);
			chooser.setDialogTitle("Load Schema File");
			chooser.setAcceptAllFileFilterUsed(true);
			DEditSimpleFileFilter filter 
			    = new DEditSimpleFileFilter(".xsd",
							"XML Schema Files");
			chooser.setFileFilter(filter);
			int returnVal = chooser.showOpenDialog(DEditor.this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
			    File file = chooser.getSelectedFile();
			    curXSDDir = file.getParentFile();
			    filename = file.getAbsolutePath();
			} else {
			    return;
			}
		    }
		    boolean enable;
		    try {
			enable = schema.loadSchema(filename);
			for (int i = 0; i < windows.size(); i++) {
			    DEditXMLWindow window 
				= (DEditXMLWindow)windows.get(i);
			    window.update(false); // validate etc.
			}
		    } catch (DEditSchemaException e) {
			JOptionPane.showMessageDialog(DEditor.this,
						      "Couldn't open "
						      + "Schema: " 
						      + e.toString());
			enable = false;
		    }
		    // reset affected components
		    resetFrameString();
		    changeEditOptions();
		    verifyItem.setEnabled(enable);
		}
	    });
	loadSchemaItem.setEnabled(allowLoadSchema && !isApplet);
	fileMenu.add(loadSchemaItem);
	// Save (== Save As if file is not yet saved)
	saveAction = new AbstractAction("Save") {
		public void actionPerformed(ActionEvent ae) {
		    saveXML();
		}
	    };
	JMenuItem saveItem = new JMenuItem(saveAction);
	saveAction.setEnabled(false); // not enabled until file is opened
	fileMenu.add(saveItem);
	// Save As
	saveAsItem = new JMenuItem("Save As...");
	saveAsItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow window = getSelectedWindow();
		    if (window == null) {
			return;
		    }
		    String openFilename = getSaveFile();
		    if (openFilename == null) {
			return;
		    }
		    window.saveFile(openFilename);
		}
	    });
	saveAsItem.setEnabled(false);
	fileMenu.add(saveAsItem);
	// Exit
	JMenuItem exitItem = new JMenuItem("Exit");
	exitItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    close();
		}
	    });
	fileMenu.add(exitItem);
	menuBar.add(fileMenu);
	// ----- Edit Menu -----
	JMenu editMenu = new JMenu("Edit");
	// clear submenus every time invoked
	editMenu.addMenuListener(new MenuListener() {
		public void menuCanceled(MenuEvent e) {}
		public void menuDeselected(MenuEvent e) {}
		public void menuSelected(MenuEvent e) { 
		    newTagBeforeMenu.removeAll();
		    newChildTagMenu.removeAll();
		    newTagAfterMenu.removeAll();
		    switchTagToMenu.removeAll();
		    pasteChildMenu.removeAll();
		}
	    });
	nodeHelpAction = new AbstractAction("Get Help For Item") {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow selWin = getSelectedWindow();
		    if (selWin == null) {
			return;
		    }
		    DEditNode selected = selWin.getSelected();
		    if (selected == null) {
			return;
		    }
		    String name = selected.getNode().getNodeName();
		    try {
			deditHelpBroker.setCurrentID(name);
		    } catch (BadIDException bie) {
			String message = "No help available.";
			JOptionPane.showMessageDialog(DEditor.this,
						      message);
			return;
		    }
		    try {
			deditHelpBroker.setDisplayed(true);
		    } catch (javax.help.UnsupportedOperationException e) {
			DEditDebug.println(this,6,"couldn't show help");
		    }
		}
	    };
	    
	// create sub menus and add listeners to determine contents 
	// on selection
	newTagBeforeMenu = new JMenu("Insert Item Above...");
	newChildTagMenu = new JMenu("Insert Sub-Item...");
	newTagAfterMenu = new JMenu("Insert Item Below...");
	switchTagToMenu = new JMenu("Switch Item To...");
	pasteChildMenu = new JMenu("Paste As Sub-Item...");
	EditMenuListener eml 
	    = new EditMenuListener(newTagBeforeMenu,
				   EditMenuListener.INSERT_ABOVE_TYPE);
	newTagBeforeMenu.addMenuListener(eml);
	eml = new EditMenuListener(newChildTagMenu,
				   EditMenuListener.INSERT_CHILD_TYPE);
	newChildTagMenu.addMenuListener(eml);
	eml = new EditMenuListener(newTagAfterMenu,
				   EditMenuListener.INSERT_BELOW_TYPE);
	newTagAfterMenu.addMenuListener(eml);
	eml = new EditMenuListener(switchTagToMenu,
				   EditMenuListener.SWITCH_TYPE);
	switchTagToMenu.addMenuListener(eml);
	eml = new EditMenuListener(pasteChildMenu,
				   EditMenuListener.PASTE_CHILD_TYPE);
	pasteChildMenu.addMenuListener(eml);
	// create Actions
	editValAction = new EditAction("Edit Value...",
				       EditAction.EDIT_TAG_TYPE);
	copyAction = new EditAction("Copy", EditAction.COPY_TYPE);
	cutAction = new EditAction("Cut", EditAction.CUT_TYPE);
	undoAction = new EditAction("Undo", EditAction.UNDO_TYPE);
	pasteBeforeAction = new EditAction("Paste Before",
					   EditAction.PASTE_BEFORE_TYPE);
	pasteAfterAction = new EditAction("Paste After",
					   EditAction.PASTE_AFTER_TYPE);
	// disable nearly everything until at least a file is opened
	newChildTagMenu.setEnabled(false);
	newTagBeforeMenu.setEnabled(false);
	newTagAfterMenu.setEnabled(false);
	switchTagToMenu.setEnabled(false);
	editValAction.setEnabled(false);
	copyAction.setEnabled(false);
	cutAction.setEnabled(false);
	undoAction.setEnabled(false);
	pasteChildMenu.setEnabled(false);
	pasteBeforeAction.setEnabled(false);
	pasteAfterAction.setEnabled(false);
	// add them all to the Edit menu
	editMenu.add(nodeHelpAction);
	editMenu.add(undoAction);
	editMenu.add(editValAction);
	editMenu.add(copyAction);
	editMenu.add(cutAction);
	editMenu.add(pasteBeforeAction);
	editMenu.add(pasteAfterAction);
	editMenu.add(pasteChildMenu);
	editMenu.add(newChildTagMenu);
	editMenu.add(newTagBeforeMenu);
	editMenu.add(newTagAfterMenu);
	editMenu.add(switchTagToMenu);
	menuBar.add(editMenu);
	// ----- Tools Menu -----
	JMenu toolsMenu = new JMenu("Tools");
	// Verify
	verifyItem = new JMenuItem("Verify");
	verifyItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    verify();
		}
	    });
	verifyItem.setEnabled(false); // disabled until Schema is loaded and 
	// file is open
	toolsMenu.add(verifyItem);
	// Show XML
	JMenuItem showXMLItem = new JMenuItem("Show XML Text");
	showXMLItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow window = getSelectedWindow();
		    if (window == null) {
			DEditDebug.println(this,3,"no active XML window");
			return;
		    }
		    window.showTextFrame();
		}
	    });
	toolsMenu.add(showXMLItem);
	menuBar.add(toolsMenu);
	// help menu inserted after plugin is loaded
	// ----- End Menus ------
	this.setJMenuBar(menuBar);
	
	/**** Popup Edit Menu ***/
	popup = new JPopupMenu();
	// clear all sub menus before invocation
	popup.addPopupMenuListener(new PopupMenuListener() {
		public void popupMenuCanceled(PopupMenuEvent e) {}
		public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {}
		public void popupMenuWillBecomeVisible(PopupMenuEvent e) { 
		    newTagBeforePopupMenu.removeAll();
		    newChildTagPopupMenu.removeAll();
		    newTagAfterPopupMenu.removeAll();
		    switchTagToPopupMenu.removeAll();
		    pasteChildPopupMenu.removeAll();
		}
	    });
	// create and disable sub menus
	newTagBeforePopupMenu = new JMenu("Insert Item Above...");
	newTagBeforePopupMenu.setEnabled(false);
	newChildTagPopupMenu = new JMenu("Insert Sub-Item...");
	newChildTagPopupMenu.setEnabled(false);
	newTagAfterPopupMenu = new JMenu("Insert Item Below...");
	newTagAfterPopupMenu.setEnabled(false);
	pasteChildPopupMenu = new JMenu("Paste Sub-Item...");
	pasteChildPopupMenu.setEnabled(false);
	switchTagToPopupMenu = new JMenu("Switch Item To...");
	switchTagToPopupMenu.setEnabled(false);
	// add listeners to fill menus on selection
	eml = new EditMenuListener(newTagBeforePopupMenu,
				   EditMenuListener.INSERT_ABOVE_TYPE);
	newTagBeforePopupMenu.addMenuListener(eml);
	eml = new EditMenuListener(newChildTagPopupMenu,
				   EditMenuListener.INSERT_CHILD_TYPE);
	newChildTagPopupMenu.addMenuListener(eml);
	eml = new EditMenuListener(newTagAfterPopupMenu,
				   EditMenuListener.INSERT_BELOW_TYPE);
	newTagAfterPopupMenu.addMenuListener(eml);
	eml = new EditMenuListener(switchTagToPopupMenu,
				   EditMenuListener.SWITCH_TYPE);
	switchTagToPopupMenu.addMenuListener(eml);
	eml = new EditMenuListener(pasteChildPopupMenu,
				   EditMenuListener.PASTE_CHILD_TYPE);
	pasteChildPopupMenu.addMenuListener(eml);
	// add them all to the popup Edit menu
	popup.add(nodeHelpAction);
	popup.add(undoAction);
	popup.add(editValAction);
	popup.add(copyAction);
	popup.add(cutAction);
	popup.add(pasteBeforeAction);
	popup.add(pasteAfterAction); 
	popup.add(pasteChildPopupMenu);
	popup.add(newChildTagPopupMenu);
	popup.add(newTagBeforePopupMenu);
	popup.add(newTagAfterPopupMenu);
	popup.add(switchTagToPopupMenu);
	// create double click invokes edit dialog listener.  It is added 
	// to every DEditXMLWindow's JTree.
	ml = new MouseAdapter() {
		public void mousePressed(MouseEvent e) {
		    // should this do more thorough click checking / change 
		    // selection / etc.???
		    if (e.getClickCount() == 2 && editValAction.isEnabled()) {
			editTag(); // change for edit in place ???
		    }
		}
	    };

	/*** The panel of icon buttons at the top of the window ***/
	JPanel toolbarWrapper = new JPanel(new BorderLayout());
	toolbar = new JToolBar();
	ClassLoader cl = DEditor.class.getClassLoader();
	// New
	ImageIcon newXMLIcon
	    = new ImageIcon(cl.getResource("icons/new_xml.jpg"));
	JButton newXMLButton = new JButton(newXMLIcon);
	newXMLButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    makeNewFile(); 
		}
	    });
	newXMLButton.setEnabled(allowNewFile);
	newXMLButton.setToolTipText("Make new XML file");
	CSH.setHelpIDString(newXMLButton,NEW_XML_ID);
	addIconButton(newXMLButton);
	// Open
	ImageIcon openXMLIcon
	    = new ImageIcon(cl.getResource("icons/open_xml.jpg"));
	JButton openXMLButton = new JButton(openXMLIcon);
	openXMLButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    openXML(); 
		}
	    });
	openXMLButton.setEnabled(allowOpenFile);
	openXMLButton.setToolTipText("Open XML file");
	CSH.setHelpIDString(openXMLButton,OPEN_XML_ID);
	addIconButton(openXMLButton);
	// Save
	ImageIcon saveXMLIcon
	    = new ImageIcon(cl.getResource("icons/save_xml.jpg"));
	JButton saveXMLButton = new JButton(saveAction);
	saveXMLButton.setIcon(saveXMLIcon);
	saveXMLButton.setToolTipText("Save XML file");
	CSH.setHelpIDString(saveXMLButton,SAVE_XML_ID);
	addIconButton(saveXMLButton);
	// Insert Above
	ImageIcon insBeforeIcon
	    = new ImageIcon(cl.getResource("icons/insert_before.jpg"));
	JButton insBeforeButton = new JButton(insBeforeIcon);
	insBeforeButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow window = getSelectedWindow();
		    if (window == null) {
			JOptionPane.showMessageDialog(DEditor.this,
					      "Select a window first.");
			return;
		    }
		    Action[] acts = window.getNewTagBeforeActions();
		    setUpActionList(acts,"Insert Item Above", "Insert"); 
		}
	    });
	insBeforeButton.setToolTipText("Insert item above");
	CSH.setHelpIDString(insBeforeButton,INSERT_ABOVE_ID);
	addIconButton(insBeforeButton);
	// Insert As Sub-Item
	ImageIcon insChildIcon
	    = new ImageIcon(cl.getResource("icons/insert_child.jpg"));
	JButton insChildButton = new JButton(insChildIcon);
	insChildButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow window = getSelectedWindow();
		    if (window == null) {
			JOptionPane.showMessageDialog(DEditor.this,
					      "Select a window first.");
			return;
		    }
		    Action[] acts = window.getNewChildTagActions();
		    setUpActionList(acts,"Insert As Sub-Item", "Insert"); 
		}
	    });
	insChildButton.setToolTipText("Insert sub-item");
	CSH.setHelpIDString(insChildButton,INSERT_CHILD_ID);
	addIconButton(insChildButton);
	// Insert Below
	ImageIcon insAfterIcon
	    = new ImageIcon(cl.getResource("icons/insert_after.jpg"));
	JButton insAfterButton = new JButton(insAfterIcon);
	insAfterButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow window = getSelectedWindow();
		    if (window == null) {
			JOptionPane.showMessageDialog(DEditor.this,
					      "Select a window first.");
			return;
		    }
		    Action[] acts = window.getNewTagAfterActions();
		    setUpActionList(acts,"Insert Item Below", "Insert"); 
		}
	    });
	insAfterButton.setToolTipText("Insert item below");
	CSH.setHelpIDString(insAfterButton,INSERT_BELOW_ID);
	addIconButton(insAfterButton);
	// Switch
	ImageIcon switchIcon
	    = new ImageIcon(cl.getResource("icons/switch.jpg"));
	JButton switchButton = new JButton(switchIcon);
	switchButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow window = getSelectedWindow();
		    if (window == null) {
			JOptionPane.showMessageDialog(DEditor.this,
					      "Select a window first.");
			return;
		    }
		    Action[] acts = window.getSwitchTagToActions();
		    setUpActionList(acts,"Switch Item To", "Switch"); 
		}
	    });
	switchButton.setToolTipText("Switch item to another item");
	CSH.setHelpIDString(switchButton,SWITCH_ID);
	addIconButton(switchButton);
	// Edit
	JButton editValButton = new JButton(editValAction);
	ImageIcon editIcon
	    = new ImageIcon(cl.getResource("icons/edit.jpg"));
	editValButton.setIcon(editIcon);
	editValButton.setToolTipText("Edit item");
	CSH.setHelpIDString(editValButton,EDIT_VALUE_ID);
	addIconButton(editValButton);
	// Copy
	JButton copyButton = new JButton(copyAction);
	ImageIcon copyIcon
	    = new ImageIcon(cl.getResource("icons/copy.jpg"));
	copyButton.setIcon(copyIcon);
	copyButton.setToolTipText("Copy");
	CSH.setHelpIDString(copyButton,COPY_ID);
	addIconButton(copyButton);
	// Cut
	JButton cutButton = new JButton(cutAction);
	ImageIcon cutIcon
	    = new ImageIcon(cl.getResource("icons/cut.jpg"));
	cutButton.setIcon(cutIcon);
	cutButton.setToolTipText("Cut");
	CSH.setHelpIDString(cutButton,CUT_ID);
	addIconButton(cutButton);
	// Paste Above
	JButton pasteBeforeButton = new JButton(pasteBeforeAction);
	ImageIcon pbIcon
	    = new ImageIcon(cl.getResource("icons/paste_before.jpg"));
	pasteBeforeButton.setIcon(pbIcon);
	pasteBeforeButton.setToolTipText("Paste copied item(s) above");
	CSH.setHelpIDString(pasteBeforeButton,PASTE_ABOVE_ID);
	addIconButton(pasteBeforeButton);
	// Paste As Sub-Item
	ImageIcon pasteChildIcon
	    = new ImageIcon(cl.getResource("icons/paste_child.jpg"));
	JButton pasteChildButton = new JButton(pasteChildIcon);
	pasteChildButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow window = getSelectedWindow();
		    if (window == null) {
			JOptionPane.showMessageDialog(DEditor.this,
					      "Select a window first.");
			return;
		    }
		    Action[] acts = getPasteChildActions(window);
		    setUpActionList(acts,"Paste Sub-Item", "Paste"); 
		}
	    });
	pasteChildButton.setToolTipText("Paste copied item(s) as sub-items");
	CSH.setHelpIDString(pasteChildButton,PASTE_CHILD_ID);
	addIconButton(pasteChildButton);
	// Paste After
	ImageIcon paIcon
	    = new ImageIcon(cl.getResource("icons/paste_after.jpg"));
	JButton pasteAfterButton = new JButton(pasteAfterAction);
	pasteAfterButton.setIcon(paIcon);
	pasteAfterButton.setToolTipText("Paste copied item(s) below");
	CSH.setHelpIDString(pasteAfterButton,PASTE_BELOW_ID);
	addIconButton(pasteAfterButton);
	//toolbarWrapper.add(toolbar,BorderLayout.WEST);
	this.getContentPane().add(toolbar,BorderLayout.NORTH);
	// ------ End button panel -------
	// set main window close behavior
	this.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
	this.addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e) {
		    close();
		}
	    });
	// load schema and xml file if passed in via command line
	if (schemaFilename != null) {
	    boolean enable = false;
	    try {
		enable = schema.loadSchema(schemaFilename);
	    } catch (DEditSchemaException e) {
		JOptionPane.showMessageDialog(this,
					      "Couldn't open Schema: " 
					      + e.toString());
		enable = false;
	    }
	    resetFrameString();
	    changeEditOptions();
	    verifyItem.setEnabled(enable);
	}
	if (plugin != null) {
	    loadPlugin(plugin);
	}
	if (filename != null) {
	    openFile(filename);
	}
    }
    
    /**
      * A menu listener designed for sub menus in the Edit and popup 
      * Edit menus.  Fired on menu invocation, it fills the chosen menu 
      * with all valid options available to it, or an unselectable "none 
      * available" option if there are none.
      */
    private class EditMenuListener implements MenuListener {
	private static final int INSERT_ABOVE_TYPE = 0;
	private static final int INSERT_CHILD_TYPE = 1;
	private static final int INSERT_BELOW_TYPE = 2;
	private static final int SWITCH_TYPE = 3;
	private static final int PASTE_CHILD_TYPE = 4;

	private int type; // which type of items to fill menu with
	private JMenu menu; // the menu being filled
	private EditMenuListener(JMenu menu,int type) {
	    this.type = type;
	    this.menu = menu;
	}
	public void menuCanceled(MenuEvent e) {}
	public void menuDeselected(MenuEvent e) {}
	public void menuSelected(MenuEvent e) {
	    if (menu.getItemCount() == 0) {
		setUpMenu(menu,type);
	    }
	}
    }

    /**
      * Class which implements action to provide uniform behavior for
      * standard edit menu and popup edit menu.
      */
    private class EditAction extends AbstractAction {
	
	private static final int NEW_TAG_BEFORE_TYPE = 0;
	private static final int NEW_CHILD_TAG_TYPE = 1;
	private static final int NEW_TAG_AFTER_TYPE = 2;
	private static final int EDIT_TAG_TYPE = 3;
	private static final int COPY_TYPE = 4;
	private static final int CUT_TYPE = 5;
	private static final int PASTE_BEFORE_TYPE = 6;
	private static final int PASTE_AFTER_TYPE = 7;
	private static final int UNDO_TYPE = 8;
	private int type;
	
	/**
	  * Constructor for edit action which in addition to performing
	  * AbstractAction's constructor assigns the type for this
	  * EditAction (which refers to an editing action).
	  * @param name the String name for the action
	  * @param thetype the type of this EditAction (i.e., cut, copy, etc.)
	  */
	private EditAction(String name, int thetype) {
	    super(name);
	    type = thetype;
	}

	public void actionPerformed(ActionEvent ae) {
	    DEditXMLWindow window = getSelectedWindow();
	    if (window == null) {
		JOptionPane.showMessageDialog(DEditor.this,
					      "Select a window first.");
		return;
	    }
	    switch(type) {
	    case EDIT_TAG_TYPE:editTag(window);break;
	    case COPY_TYPE:copySubTree(window);break;
	    case CUT_TYPE:cutSubTree(window);break;
	    case PASTE_BEFORE_TYPE:pasteSibTreeBefore(window);break;
	    case PASTE_AFTER_TYPE:pasteSibTreeAfter(window);break;
	    case UNDO_TYPE:undo(window);break;
	    }
	}
	
    }

    /**
      * Class used by DEditor to implement the display of the popup menu.
      */
    class PopupListener extends MouseAdapter {
	public void mousePressed(MouseEvent e) {
	    maybeShowPopup(e);
	}
	
	public void mouseReleased(MouseEvent e) {
	    maybeShowPopup(e);
	}
	
	private void maybeShowPopup(MouseEvent e) {
	    if (e.isPopupTrigger()) {
		popup.show(e.getComponent(),e.getX(), e.getY());
	    }
	}
    }

    /**
      * Adds a JButton (which should have an associated icon) to the icon 
      * button panel at the top of the window.  Any text will be removed and
      * the icon will be resized to fit existing buttons.  As this is done 
      * only once, before JFrame.show(), it's acceptable to use smooth scaling.
      * @param button the JButton to add to the button panel
      */
    public void addIconButton(JButton button) {
	Icon icon = button.getIcon();
	if (icon != null && icon instanceof ImageIcon) {
	    ImageIcon imic = (ImageIcon)icon;
	    Image image = imic.getImage();
	    image = image.getScaledInstance(BUTTON_ICON_SIZE,
					    BUTTON_ICON_SIZE,
					    Image.SCALE_SMOOTH);
	    imic.setImage(image);    
	}
	button.setText("");
	button.setMargin(NO_MARGIN_INSETS);
	button.setSize(buttonDim);
	toolbar.add(button);
    }
    
    /**
      * Convenient helper which constructs a menu with items such that 
      * when there are cutoff items in a menu, it adds a menu to that 
      * menu and starts putting more items there until cutoff is reached, 
      * and then puts another menu in that menu, etc.
      * @param menu the starting, main JMenu to fill
      * @param subItems the JMenuItem's to add to menu
      * @param cutoff maximum number of items in a menu - 1 (for the 
      * link to the next menu)
      */
    private static void fillMenuRecursively(JMenu menu, JMenuItem[] subItems,
					    int cutoff) {
	int n = 0;
	while (n < subItems.length) {
	    for (int j = 0; j < cutoff; j++) {
		if (n+j >= subItems.length) {
		    return;
		}
		menu.add(subItems[n+j]);		
	    }
	    n += cutoff;
	    JMenu temp = new JMenu("More..."); // change?  configureable???
	    menu.add(temp);
	    menu = temp;
	}
    }

    /**
      * Writes a remote file using the servlet from BhpViewer.
      * @param filename the name of the file to write
      * @param contents the contents of the file to write
      */
    public void writeRemoteFile(String filename, String contents) {
	try {
	    URL url = new URL(getProperty("servletSaveFilecfgURL"));
	    URLConnection conn = url.openConnection();
            conn.setDoOutput(true);
            PrintWriter out = new PrintWriter(new OutputStreamWriter(conn.getOutputStream()));
            out.print("&fileName=");
            out.print(filename);
            out.print("&content=");
	    //String filtered = contents.replaceAll("[&]",":::{052}");
	    //filtered = filtered.replaceAll("[+]",":::{053}");
	    //filtered = filtered.replaceAll("[-]",":::{054}");
	    //System.out.println("more replacesments: &,+,-");
	    //out.print(filtered);
	    out.print(URLEncoder.encode(contents,"UTF-8"));
            out.close();
            BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream()));
            System.out.println("writeRemoteFile " + reader.readLine());
	    //if (!filtered.equals(contents)) {
	    //	String loc = "$PRGROOT/xmledit/beta/";
		//executeRemote(loc + "convert.sh " + filename);
		//executeRemote("perl convert.pl " + filename);
	    //	executeRemote("perl " + loc +  "convert.pl " + filename);
		//executeRemote("date");
	    // }
	} catch (Exception e) {
	    System.err.println("writing exception for " + filename 
			       + ": " + e.toString());
	    e.printStackTrace();
	}
    }

    /**
      * Gets the InputStream containing data from the remote file specified 
      * by "filename" using the servlet from BhpViewer
      * @param filename the name of the remote file whose data will be 
      * retrievable through the returned InputStream
      * @return an InputStream containing the data from file "filename"
      */
    public InputStream getRemoteFileInputStream(String filename) {
	try {
	    String urlString = getProperty("servletReadFilecfgURL");
	    URL url = new URL(urlString + "?FileName=" + filename);
	    URLConnection conn = url.openConnection();
	    InputStream is = conn.getInputStream();
	    return is;
	} catch (Exception e) {
	    System.err.println("error getting input stream for " + filename
			       + ": " + e.toString());
	}
	return null;
    }

    /**
      * Reads the contents of a remote file and returns them as a String
      * using the servlet from BhpViewer to get the InputStream containing the
      * data (with getRemoteFileInputStream(String)).
      * @param filename the name of the file whose contents will be retrieved
      * @return the contents of "filename" as a String
      */
    public String readRemoteFile(String filename) {
	try {
	    InputStream is = getRemoteFileInputStream(filename);
	    if (is == null) {
		return null;
	    }
	    StringWriter out = new StringWriter();
	    int c;
	    
	    while ((c = is.read()) != -1)
		out.write(c);
	    
	    is.close();
	    return out.toString();
	} catch (Exception e) {
	    System.err.println("reading exception for " + filename 
			       + ": " + e.toString());
	}
	return null;
    }

    /**
      * Executes a command remotely using the servlet from BhpViewer. 
      * @param command the command to run remotely (note: it will be put in a 
      * script and run)
      * @param return 0 for success, -1 for failure -- eventually should be 
      * command's return value if possible.
      */
    public int executeRemote(String command) {
	System.out.println("Command: " + command);
	String urlString = getProperty("servletExecuteScriptURL");
        String outDir = getProperty("temporaryDirectory");
        try {
            URL url = new URL(urlString);
            URLConnection conn = url.openConnection();
            conn.setDoInput(true);
            conn.setDoOutput(true);
            conn.setRequestProperty("Content-type", "application/x-www-form-urlencoded");
            PrintWriter out = new PrintWriter(new OutputStreamWriter(conn.getOutputStream()));
            out.print("&scriptType=new");
            out.print("&outputDir=");
            out.print(java.net.URLEncoder.encode(outDir,"UTF-8"));
            out.print("&clientMessage=");
            out.write(java.net.URLEncoder.encode(command,"UTF-8"));
            out.flush();
            out.close();

            BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream()));
            String responseLine = reader.readLine();
            while(responseLine != null) {
                System.out.println("ExecuteRemote : " + responseLine);
                responseLine = reader.readLine();
            }
	    return 0;
        } catch (Exception ex) {
            System.out.println("DEditor.executeRemote Exception " 
			       + ex.toString());
            ex.printStackTrace();
        }
	return -1;
    }

    /**
      * Creates a DOM from an XML, locating the file locally only if it is 
      * not an applet.
      * @param xmlFileName the name of the file to open
      * @param return a DOM Document for the XML in "xmlFileName"
      */
    public Document createDOMforXML(String xmlFileName) {
	if (!isApplet) {
	    return createDOMforLocalXML(xmlFileName);
	}
	String contents = readRemoteFile(xmlFileName);
	//contents = DEditXMLWindow.decode(contents);
	//System.out.println("contents of file:\n" + contents);
	InputSource is = new InputSource(new StringReader(contents));
	DocumentBuilderFactory dbfactory=DocumentBuilderFactory.newInstance();
	try {
	    DocumentBuilder builder = dbfactory.newDocumentBuilder();
	    return builder.parse(is);
	} catch (SAXException sxe) {
	    // Error generated during parsing
	    Exception  x = sxe;
	    if (sxe.getException() != null)
		x = sxe.getException();
	    openErrString = x.toString();
	    return null;
	} catch (ParserConfigurationException pce) {
	    // Parser with specified options can't be built
	    openErrString = pce.toString();
	    return null;
	} catch (IOException ioe) {
	    // I/O error
	    openErrString = ioe.toString();
	    return null;
	} catch (Exception e) {
	    openErrString = e.toString();
	    return null;
	}
    }

    /**
      * Creates a DOM Document given the pathname of a local XML file.
      * @param xmlFileName the filename of the local XML to be converted to DOM
      * @return a Document derived from the input local XML file
      */
    public static Document createDOMforLocalXML(String xmlFileName) {
	DocumentBuilderFactory dbfactory=DocumentBuilderFactory.newInstance();
	try {
	    DocumentBuilder builder = dbfactory.newDocumentBuilder();
	    URL xmlURL = ClassLoader.getSystemResource(xmlFileName);
	    if (xmlURL != null) {
		xmlFileName = xmlURL.toString();
	    }
	    return builder.parse(xmlFileName);
	} catch (SAXException sxe) {
	    // Error generated during parsing
	    Exception  x = sxe;
	    if (sxe.getException() != null)
		x = sxe.getException();
	    openErrString = x.toString();
	    return null;
	} catch (ParserConfigurationException pce) {
	    // Parser with specified options can't be built
	    openErrString = pce.toString();
	    return null;
	} catch (IOException ioe) {
	    // I/O error
	    openErrString = ioe.toString();
	    return null;
	} catch (Exception e) {
	    openErrString = e.toString();
	    return null;
	}
    }

    /**
      * Changes the available options under the edit menu based on the 
      * currently selected DEditXMLWindow.
      */
    public void changeEditOptions() {
	DEditXMLWindow current = getSelectedWindow();
	if (current != null) {
	    changeEditOptions(current);
	} else {
	    changeEditOptions(null,null);
	}
    }

    /**
      * Helper or regularly called method which changes the edit options for
      * the passed in DEditXMLWindow which should in all cases ==
      * getSelectedWindow().
      * @param window the currently selected DEditXMLWindow, 
      * == getSelectedWindow()
      */
    public void changeEditOptions(DEditXMLWindow window) {
	DEditNode[] selectedNodes = window.getMultipleSelected();
	changeEditOptions(window,selectedNodes);
    }

    /**
      * Changes the available options under the edit menu based on the 
      * input DEditXMLWindow (which should be the currently selected 
      * DEditXMLWindow) and the given selected nodes.  If more than one 
      * node is selected, only copy and cut can be chooseable.
      * @param window the DEditXMLWindow to use to determine what's 
      * allowable -- should be the currently selected DEditXMLWindow
      * @param selectedNodes the currently selected DEditNodes
      */
    public void changeEditOptions(DEditXMLWindow window,
				  DEditNode[] selectedNodes) {
	long startTime = System.currentTimeMillis();
	undoAction.setEnabled(window != null && window.canUndo());
	if (selectedNodes == null || selectedNodes.length != 1) {
	    // # of selected items != 1, most stuff is disabled
	    copyAction.setEnabled(false);
	    cutAction.setEnabled(false);
	    pasteChildMenu.setEnabled(false);
	    pasteChildPopupMenu.setEnabled(false);
	    pasteBeforeAction.setEnabled(false);
	    pasteAfterAction.setEnabled(false);
	    editValAction.setEnabled(false);
	    newTagBeforeMenu.setEnabled(false);
	    newChildTagMenu.setEnabled(false);
	    newTagAfterMenu.setEnabled(false);
	    switchTagToMenu.setEnabled(false);
	    newTagBeforePopupMenu.setEnabled(false);
	    newChildTagPopupMenu.setEnabled(false);
	    newTagAfterPopupMenu.setEnabled(false);
	    switchTagToPopupMenu.setEnabled(false);
	    // re-enable copy and cut for multiple selection
	    if (selectedNodes != null && selectedNodes.length > 1) {
		copyAction.setEnabled(true);
		cutAction.setEnabled(true);
	    }
	} else {
	    // enable major menus
	    pasteChildMenu.setEnabled(true);
	    pasteChildPopupMenu.setEnabled(true);
	    newTagBeforeMenu.setEnabled(true);
	    newChildTagMenu.setEnabled(true);
	    newTagAfterMenu.setEnabled(true);
	    switchTagToMenu.setEnabled(true);
	    newTagBeforePopupMenu.setEnabled(true);
	    newChildTagPopupMenu.setEnabled(true);
	    newTagAfterPopupMenu.setEnabled(true);
	    switchTagToPopupMenu.setEnabled(true);
	    // only one node is selected, the first and only node in the array
	    DEditNode selected = selectedNodes[0];
	    // copy should be enabled unless plugin says otherwise
	    copyAction.setEnabled(true);
	    if (plugin != null) {
		copyAction.setEnabled(plugin.isCopyable(selected));
	    }
	    // check current state (to compare to later state)
	    boolean schemaReady = schema.isReady();
	    boolean validBefore = (schemaReady && window.isValid());
	    // determine if it is editable -- ask plugin after normal check
	    if (schemaReady) {
		Node n = selected.getNode();
		boolean enabled = ((schema.isLeaf(n) 
				    && schema.getDataType(n) != null)
				   || schema.getAttributes(n) != null);
		if (plugin != null) {
		    enabled = plugin.isEditable(selected,enabled);
		}
		editValAction.setEnabled(enabled);
	    } else {
		// schema's not ready, but we'll let them try and edit
		editValAction.setEnabled(true);
	    }
	    // can't cut root node
	    if (selected.isRoot()) {
		cutAction.setEnabled(false);
	    } else {
		// check if cutting is allowed (legal && plugin allows)
		DEditNode parent = (DEditNode)selected.getParent();
		int index = parent.getIndex(selected);
		Node next = selected.getNode().getNextSibling();
		selected.removeNode();
		boolean enable = !validBefore || window.isValid();
		if (plugin != null) {
		    enable = enable && plugin.isCuttable(selected);
		}
		cutAction.setEnabled(enable);
		if (next == null) {
		    parent.reinsertNodeAtEnd(selected);
		} else {
		    parent.reinsertNodeBefore(selected,next);
		}
	    }
	    // if nothing is copied then nothing is pasteable
	    if (copied == null) {
		pasteChildPopupMenu.setEnabled(false);
		pasteChildMenu.setEnabled(false);
		pasteAfterAction.setEnabled(false);
		pasteBeforeAction.setEnabled(false);
	    } else {
		// check paste before and after validity, leave paste child 
		// for when its menu is expanded
		boolean enablePasteBefore = true;
		boolean enablePasteAfter = true;
		if (plugin != null) {
		    enablePasteBefore = plugin.pasteBeforeAllowed(selected);
		    enablePasteAfter = plugin.pasteAfterAllowed(selected);
		}
		if (enablePasteBefore) {
		    enablePasteBefore = !validBefore 
			|| window.pasteBeforeCheck(selected,copied);
		}
		if (enablePasteAfter) {
		    enablePasteAfter = !validBefore 
			|| window.pasteAfterCheck(selected,copied);
		}
		pasteBeforeAction.setEnabled(enablePasteBefore);
		pasteAfterAction.setEnabled(enablePasteAfter);		
	    }
	    long totalTime = System.currentTimeMillis() - startTime;
	    if (selected != null) {
		DEditDebug.println(this,5,"total time = " + totalTime
				   + ", " + selected.getNode().getNodeName());
	    }
	}
	// update selected window
	if (window != null) {
	    window.update(false);
	}
    }

    /**
      * Sets up and displays a dialog containing a JList which contains 
      * all the Action's passed each.  All actions are assumed to have 
      * something useful provided by toString() in order to display 
      * nicely in the JList.  When the "Accept"-type button is chosen, 
      * the selected action is performed.
      * @param actions the Action's that will be displayed in the JList
      * @param title the title for the dialog with the JList in it
      * @param acceptMessage the message to put on the "Accept" button
      */
    private void setUpActionList(Action[] actions,String title,
				 String acceptMessage) {
	if (actions == null || actions.length == 0) {
	    JOptionPane.showMessageDialog(this,"No options available for "
					  + "\"" + title + "\".");
	    return;
	}
	//DEditDebug.println(this,6,"C'mon");
	final JDialog selectDialog = new JDialog(this,title,true);
	JPanel mainPanel = new JPanel(new BorderLayout());
	final JList items = new JList(actions);
	//items.setEditable(false);
	items.setSelectedIndex(0);
	JScrollPane itemsPane = new JScrollPane(items);
	mainPanel.add(itemsPane,BorderLayout.NORTH);
	JButton acceptButton = new JButton(acceptMessage);
	acceptButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    Action act = (Action)items.getSelectedValue();
		    if (act == null) {
			return;
		    }
		    act.actionPerformed(ae);
		    selectDialog.dispose();
		}
	    });
	mainPanel.add(acceptButton,BorderLayout.WEST);
	JButton cancelButton = new JButton("Cancel");
	cancelButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    selectDialog.dispose();
		}
	    });
	mainPanel.add(cancelButton,BorderLayout.EAST);
	selectDialog.getContentPane().add(mainPanel);
	selectDialog.pack();
	selectDialog.setLocationRelativeTo(this);
	selectDialog.show();
    }

    /**
      * Sets up a particular menu which is passed in and with items determined 
      * by the integer passed in.
      * @param menu the menu which will be filled
      * @param which determines what type of items to get (see code if 
      * you want to know which is which)
      */
    private void setUpMenu(JMenu menu, int which) {
	long start = System.currentTimeMillis();
	DEditXMLWindow window = getSelectedWindow();
	undoAction.setEnabled(window != null && window.canUndo());
	if (window == null) {
	    return;
	}
	DEditNode selected = window.getSelected();
	if (selected == null) {
	    menu.setEnabled(false);
	    return;
	}
	Action[] actions = null;
	boolean allowable = true;
	if (which == 0) {
	    if (plugin == null || plugin.newTagBeforeAllowed(selected)) {
		actions = window.getNewTagBeforeActions(selected);
	    }
	} else if (which == 1) {
	    if (plugin == null || plugin.newChildTagAllowed(selected)) {
		actions = window.getNewChildTagActions(selected);
	    } 
	} else if (which == 2) {
	    if (plugin == null || plugin.newTagAfterAllowed(selected)) {
		actions = window.getNewTagAfterActions(selected);
	    }
	} else if (which == 3) {
	    if (plugin == null || plugin.switchTagToAllowed(selected)) {
		actions = window.getSwitchTagToActions(selected);
	    }
	} else {
	    if (plugin == null || plugin.pasteChildAllowed(selected)) {
		actions = window.getPasteChildActions(selected,copied);
	    }
	}
	if (actions == null || actions.length == 0) { // nothing available
	    JMenuItem emptyItem = new JMenuItem("None available");
	    emptyItem.setEnabled(false);
	    menu.add(emptyItem);
	    menu.setEnabled(true);
	} else {
	    // create menu items from the actions...could we just 
	    // add actions to menu instead???
	    JMenuItem[] subItems = new JMenuItem[actions.length];
	    for (int j = 0; j < subItems.length; j++) {
		subItems[j] = new JMenuItem(actions[j]);
	    }
	    menu.setEnabled(true);
	    fillMenuRecursively(menu,subItems,DEFAULT_CUTOFF);
	}
	long stop = System.currentTimeMillis();
	DEditDebug.println(this,6,"set up menu " + which + ": "
			   + (stop - start));
    }

    /**
      * Creates a new XML file based on the currently loaded XML Schema.
      */
    private void makeNewFile() {
	if (!schema.isReady()) {
	    JOptionPane.showMessageDialog(this,"Please load an XML Schema template first.");
	    return;
	}
	// create new document
	DocumentBuilderFactory dbfactory=DocumentBuilderFactory.newInstance();
	DocumentBuilder builder;
	try {
	    builder = dbfactory.newDocumentBuilder();
	} catch (ParserConfigurationException pce) {
	    JOptionPane.showMessageDialog(this,"Can't make new file :"
					  + pce.toString());
	    return;
	}
	Document document = builder.newDocument();
	// create top of document, put it in, then create a DEditXMLWindow
	// from that and other important components
	Node n = document.createElement(schema.getTopName());
	document.appendChild(n); 
	MouseListener[] mls = new MouseListener[]{ml,popupListener};
	DEditXMLWindow newWindow = new DEditXMLWindow(this,jdp,schema,
						      renderer,document,
						      mls,plugin,true);
	// if there are too many windows open, dispose the new one
	if (windows.size() < MAX_WINDOWS) {
	    windows.add(newWindow);
	} else {
	    JOptionPane.showMessageDialog(this,"Too many windows open.  Close "
					  + "some before opening any new "
					  + "windows.");
	    newWindow.dispose();
	}
	// update the editing options to those of the new window and reset
	// affected compoenents
	changeEditOptions(newWindow);
	//copied = null;
	resetFrameString();
	saveAsItem.setEnabled(true);
	saveAction.setEnabled(true);
    }

    /**
      * Opens an XML file by displaying a file chooser, getting a file, and 
      * then calling openFile() on the filename chosen.  Uses 
      * RemoteFileChooser when running as an applet.
      */
    private void openXML() {
	if (isApplet) {
	    String msg = "Specify remote URL";
	    RemoteFileChooser rfc 
		= new RemoteFileChooser(this,"Open XML","loadXmlURL","Open");
	    int reFchooser = rfc.showDialog();
	    String filename = null;
	    if (reFchooser == JFileChooser.APPROVE_OPTION) {
		String pathString = rfc.getPath();
		String name = rfc.getName();
		filename = pathString + "/" + name;
		System.out.println("(" + name + ")OOpenSession " + filename);
	    } else {
		return;
	    }
	    openFile(filename);
	    return;
	}
	JFileChooser chooser = new JFileChooser(curXMLDir);
	chooser.setDialogTitle("Open XML File");
	chooser.setAcceptAllFileFilterUsed(true);
	DEditSimpleFileFilter filter 
	    = new DEditSimpleFileFilter(".xml","XML Files");
	chooser.setFileFilter(filter);
	int returnVal = chooser.showOpenDialog(DEditor.this);
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	    File f = chooser.getSelectedFile();
	    curXMLDir = f.getParentFile();
	    openFile(f.getAbsolutePath());
	}
    }

    /**
      * Saves an XML file. First it check to see if it has been saved yet.
      * If so, saves it to its current name.  Otherwise it prompts the 
      * user with a file chooser dialog.
      */
    private void saveXML() {
	 DEditXMLWindow window = getSelectedWindow();
	 if (window == null) {
	     JOptionPane.showMessageDialog(DEditor.this,
					   "Select a window first.");
	     return;
	 }
	 String openFilename = null;
	 if (window.unsaved()) {
	     openFilename = getSaveFile();
	     if (openFilename == null) {
		 return;
	     }
	 } else {
	     openFilename = window.getOpenFilename();
	 }
	 window.saveFile(openFilename);
    }

    /**
      * Opens an XML file given a filename.
      * @param filename the name of the XML file to open.
      */
    private void openFile(String filename) {
	Document document = createDOMforXML(filename);
	if (document == null) {
	    JOptionPane.showMessageDialog(this,"Couldn't open file: " 
					  + openErrString);
	    return;
	}
	// just like makeNewFile from here on
	MouseListener[] mls = new MouseListener[]{ml,popupListener};
	DEditXMLWindow newWindow = new DEditXMLWindow(this,jdp,schema,
						      filename,renderer,
						      document,mls,plugin);
	if (newWindow == null) {
	    return;
	}
	if (windows.size() < MAX_WINDOWS) {
	    windows.add(newWindow);
	} else {
	    JOptionPane.showMessageDialog(this,"Too many windows open.  Close "
					  + "some before opening any new "
					  + "windows.");
	    newWindow.dispose();
	}
	changeEditOptions(newWindow);
	resetFrameString();
	saveAction.setEnabled(true);
	saveAsItem.setEnabled(true);
	if (schema.isReady()) {
	    verifyItem.setEnabled(true);
	} else {
	    verifyItem.setEnabled(false);
	}
    }
    
    /**
      * Gets the name of a file to use to save the currently open XML.  Uses 
      * RemoteFileChooser when running as an applet.
      * @return a name of a file to use to save the currently open XML
      */
    public String getSaveFile() {
	if (isApplet) {
	    String msg = "Specify remote URL";
	    RemoteFileChooser rfc 
		= new RemoteFileChooser(this,"Save XML","SaveXmlURL","Save");
	    int reFchooser = rfc.showDialog();
	    String filename = null;
	    if (reFchooser == JFileChooser.APPROVE_OPTION) {
		String pathString = rfc.getPath();
		String name = rfc.getName();
		filename = pathString + "/" + name;
		System.out.println("(" + name + ")OOpenSession " + filename);
	    } else {
		return null;
	    }
	    return filename;
	}
	JFileChooser chooser = new JFileChooser(curXMLDir);
	chooser.setDialogTitle("Save XML File");
	chooser.setAcceptAllFileFilterUsed(true);
	DEditSimpleFileFilter filter 
	    = new DEditSimpleFileFilter(".xml","XML Files");
	chooser.setFileFilter(filter);
	int returnVal = chooser.showSaveDialog(this);
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	    File f = chooser.getSelectedFile();
	    curXMLDir = f.getParentFile();
	    return f.getAbsolutePath();
	} else {
	    return null;
	}
    }

    /**
      * Resets the string on top of the JFrame main window.
      */
    private void resetFrameString() {
	String frameString = "XML Editor";
	String schemaFilename = schema.getFilename();
	if (schemaFilename != null) {
	    frameString += ", schema " + schemaFilename;
	}
	setTitle(frameString);
    }
    
    /**
      * Calls undo on the last action performed in DEditXMLWindow window.
      * @param window the window to call undo on
      */
    private void undo(DEditXMLWindow window) {
	window.undo();
	undoAction.setEnabled(window.canUndo());
    }

    /**
      * Edits the currently selected node in the currently selected window.
      */
    private void editTag() {
	DEditXMLWindow window = getSelectedWindow();
	if (window != null) {
	    editTag(window);
	}
    }

    /**
      * Edits information about a tag and its element.
      * @param window the window to get the selected tag from and edit
      */
    private void editTag(DEditXMLWindow window) {
	/* if it is a leaf, show change value and/or set attributes dialog. */
	/* if it is a branch, show add element dialog if they can add plus
	 *      attributes. */
	//window.addToUndoBuffer(); ???
	DEditNode toEdit = window.getSelected();
	if (toEdit == null) {
	    DEditDebug.println(this,0,"improper selection");
	    return;
	}
	
	// give plugin a chance to intercept
	if (plugin != null) {
	    boolean stop = plugin.edit(toEdit);
	    if (stop) {
		DEditDebug.println(this,3,"editing stopped by plugin.");
		return;
	    }
	}

	// determine if it is a leaf 
	boolean isLeaf = toEdit.isLeaf();
	Node node = toEdit.getNode();
	String oldVal = null;
	NodeList nlist = node.getChildNodes();
	Node child = null;
	for (int i = 0; i < nlist.getLength(); i++) {
	    Node temp = nlist.item(i);
	    if (temp.getNodeType() == Node.TEXT_NODE) {
		child = temp;
		break;
	    }
	}
	// determine its data type
	DEditDataType dt = null;
	if (isLeaf) {
	    if (child == null) {
		child = node.getOwnerDocument().createTextNode("");
		node.appendChild(child);
	    }
	    oldVal = child.getNodeValue();
	    dt = schema.getDataType(node);
	}
	
	// get state information
	boolean validBefore = false;
	if (schema.isReady()) {
	    validBefore = window.isValid();
	}
	// no value, no attributes -- no editing
	DEditAttributeInfo[] attInfo = schema.getAttributes(node);
	if (!isLeaf && attInfo == null) {
	    DEditDebug.println(this,3,"nothing to edit.");
	    return;
	}
	
	DEditAttributeInfo[] backup = null; // holds old values should 
	// a reset be necessary
	// fill up attInfo with current attribute information and values
	NamedNodeMap nnm = null;
	if (attInfo != null) {
	    nnm = node.getAttributes();
	    if (nnm == null) {
		DEditDebug.println(this,3,"NOT an Element!!");
		return;
	    }
	    if (nnm.getLength() > 0) {
		for (int i = 0; i < attInfo.length; i++) {
		    Attr att = (Attr)nnm.getNamedItem(attInfo[i].name);
		    if (att != null) {
			attInfo[i].exists = true;
			attInfo[i].value = att.getValue();
		    } 
		}
	    }
	    backup = new DEditAttributeInfo[attInfo.length];
	    for (int i = 0; i < backup.length; i++) {
		backup[i] = attInfo[i].copy();
	    }
	}
	// idetermine which dialog to show and show it
	String newVal = null;
	String titleString = "Edit " + toEdit.getDisplayName();
	if (isLeaf && oldVal != null && dt != null) {
	    if (attInfo != null) {
		newVal = dialog.showAttValDialog(this,titleString,
						 dt,oldVal,attInfo);
	    } else {
		newVal = dialog.showValueDialog(this,titleString,dt,oldVal);
	    }
	    if (newVal == null) { // return if canceled
		return;
	    }
	} else if (attInfo != null) {
	    attInfo = dialog.showAttributeDialog(this,titleString,attInfo);
	} else { // not a leaf, no attributes -- no editing
	    return;
	}
	// update attributes if not canceled
	if (attInfo != null) {
	    for (int i = 0; i < attInfo.length; i++) {
		if (attInfo[i].exists) {
		    Attr att = window.createAttribute(attInfo[i].name);
		    att.setValue(attInfo[i].value);
		    nnm.setNamedItem(att);
		} else {
		    if (nnm.getNamedItem(attInfo[i].name) != null) {
			nnm.removeNamedItem(attInfo[i].name);
		    }
		}
	    }
	}
	// set node value if applicable
	if (isLeaf && child != null) {
	    child.setNodeValue(newVal);
	}
	// check to make sure xml is still valid.  If not, move back to 
	// old values and inform user
	if (validBefore && !window.isValid()) {
	    if (isLeaf) {
		child.setNodeValue(oldVal);
	    }
	    if (attInfo != null) {
		attInfo = backup;
		for (int i = 0; i < attInfo.length; i++) {
		    if (attInfo[i].exists) {
			Attr att=window.createAttribute(attInfo[i].name);
			DEditDebug.println(this,0,attInfo[i].value);
			att.setValue(attInfo[i].value);
			nnm.setNamedItem(att);
		    } else {
			if (nnm.getNamedItem(attInfo[i].name) != null) {
			    nnm.removeNamedItem(attInfo[i].name);
			}
		    }
		}
		
	    }
	    JOptionPane.showMessageDialog(this,
					  "Can't set tag: " 
					  + schema.getVerifyString());
	}
	window.update(true);
	window.nodeChanged(toEdit);
    }

    /**
      * Checks validity of XML document given schema and informs user.
      */
    private void verify() {
	DEditXMLWindow window = getSelectedWindow();
	if (window == null) {
	    JOptionPane.showMessageDialog(this,"Select a tree to verify.");
	} else if (!schema.isReady()) {
	    JOptionPane.showMessageDialog(this,"Please load schema first.");
	} else {
	    window.update(false);
	    JOptionPane.showMessageDialog(this,schema.getVerifyString());
	}
    }

    /**
      * Gets the currently selected DEditXMLWindow.
      * @return the currently selected DEditXMLWindow
      */
    public DEditXMLWindow getSelectedWindow() {
	for (int i = 0; i < windows.size(); i++) {
	    DEditXMLWindow temp = (DEditXMLWindow)windows.get(i);
	    if (temp.isSelected()) {
		return temp;
	    } 
	}
	return null;
    }

    /**
      * I'd put it in DEditXMLWindow but copied can't be final.  Returns 
      * an array of Actions which represent pasting the copied DEditNode's 
      * (in order) before various DEditNodes that are children of the currently
      * selected one, or after all of the currently selected one's children.
      * @return an array of Actions that will paste the copied DEditNode's as 
      * a child of the currently selected node.
      */
    private Action[] getPasteChildActions(DEditXMLWindow window) {
	return window.getPasteChildActions(copied);
    }

    /**
      * Copies the currently selected subtree.
      */
    private void copySubTree(DEditXMLWindow window) {
	DEditNode[] toCopy = window.getMultipleSelected();
	if (toCopy == null) {
	    return;
	}
	copied = new DEditNode[toCopy.length];
	for (int i = 0; i < copied.length; i++) {
	    copied[i] = toCopy[i].copy();
	}
	changeEditOptions(window,toCopy);	
    }

    /**
      * Cuts the currently selected subtree.
      */
    private void cutSubTree(DEditXMLWindow window) {
	DEditNode[] toCopy = window.cutSubTree(window.getMultipleSelected());
	if (toCopy == null) {
	    return;
	}
	copied = new DEditNode[toCopy.length];
	for (int i = 0; i < copied.length; i++) {
	    copied[i] = toCopy[i].copy();
	}
	changeEditOptions(window);
    }
    
    /**
      * Pastes the currently copied subtree as the sibling before the
      * currently selected one (vertically above).
      */
    private void pasteSibTreeBefore(DEditXMLWindow window) {
	window.pasteSibTreeBefore(copied);
	changeEditOptions(window);
    }
    
    /**
      * Pastes the currently copied subtree as the sibling after the
      * currently selected one (vertically below).
      */
    private void pasteSibTreeAfter(DEditXMLWindow window) {
	window.pasteSibTreeAfter(copied);
	changeEditOptions(window);
    }

    /**
      * Gets an array of 2-sized arrays from the Document d which correspond
      * to the key and value pairs for customization of the edit dialog box
      * @param d the Document to derived the pairs from
      * @return an array of pairs of strings which are key, value pairs
      */
    private String[][] getLeafPairs(Document d) {
	return getPairs("edit_dialog",d);
    }

    /**
      * Derives an array of pairs of strings which are key, value pairs 
      * for customizing some aspect of the editor.  It gets the pairs for 
      * the feature specified by tagName (i.e., <tagName>) found in Document
      * d.
      * @param tagName the tag name to search for customizations for
      * @param d the Document in which to search for pairs
      * @return an array of pairs of strings which are key, value pairs
      */
    private String[][] getPairs(String tagName, Document d) {
	NodeList nl = d.getElementsByTagName(tagName);
	if (nl != null && nl.getLength() > 0) {
	    Node n = nl.item(0);
	    NodeList customs = n.getChildNodes();
	    Vector forRet = new Vector();
	    for (int i = 0; i < customs.getLength(); i++) {
		Node custom = customs.item(i);
		if (custom.getNodeType() == Node.ELEMENT_NODE) {
		    String name = null, value = null, effect = null;
		    NodeList vals = custom.getChildNodes();
		    for (int j = 0; j < vals.getLength(); j++) {
			Node temp = vals.item(j);
			String nodeName = temp.getNodeName();
			if (nodeName.equals("name")) {
			    name = temp.getFirstChild().getNodeValue();
			} else if (nodeName.equals("value")) {
			    value = temp.getFirstChild().getNodeValue();
			} else if (nodeName.equals("effect")) {
			    effect = temp.getFirstChild().getNodeValue();
			} else if (nodeName.equals("any_value")) {
			    value = "";
			}
		    }
		    if (name != null && value != null && effect != null) {
			String[] customization = new String[] {name + "=" 
							       + value,
							       effect};
			forRet.add(customization);
		    } 
		}
	    }
	    String[][] ret = new String[forRet.size()][2];
	    forRet.toArray(ret);
	    return ret;
	} else {
	    return null;
	}
    }

    /**
      * Gets an array of 2-sized arrays from the Document d which correspond
      * to the key and value pairs for customization of the display of text
      * for elements in the JTree
      * @param d the Document to derived the pairs from
      * @return an array of pairs of strings which are key, value pairs
      */
    private String[][] getRendererTextPairs(Document d) {
	return getPairs("tree_text_display",d);
    }

    /**
      * Gets an array of 2-sized arrays from the Document d which correspond
      * to the key and value pairs for customization of the display of icons
      * for elements in the JTree
      * @param d the Document to derived the pairs from
      * @return an array of pairs of strings which are key, value pairs
      */
    private String[][] getRendererIconPairs(Document d) {
	NodeList nl = d.getElementsByTagName("tree_icon_display");
	if (nl != null && nl.getLength() > 0) {
	    Node n = nl.item(0);
	    NodeList customs = n.getChildNodes();
	    Vector forRet = new Vector();
	    for (int i = 0; i < customs.getLength(); i++) {
		Node custom = customs.item(i);
		String icon_name = null, node_names = null;
		NodeList vals = custom.getChildNodes();
		for (int j = 0; j < vals.getLength(); j++) {
		    Node temp = vals.item(j);
		    String nodeName = temp.getNodeName();
		    if (nodeName.equals("icon_name")) {
			icon_name = temp.getFirstChild().getNodeValue();
		    } else if (nodeName.equals("nodes_for_icon")) {
			node_names = temp.getFirstChild().getNodeValue();
		    } else {
			
		    }
		}
		if (icon_name != null && node_names != null) {
		    String[] customization = new String[]{icon_name,
							  node_names};
		    forRet.add(customization);
		}
	    }
	    String[][] ret = new String[forRet.size()][2];
	    forRet.toArray(ret);
	    return ret;
	} else {
	    return null;
	}
    }

    /**
      * Removes a DEditXMLWindow from the DEditor.
      * @param gone the DEditXMLWindow to remove.
      */
    public void removeXMLWindow(DEditXMLWindow gone) {
	windows.remove(gone);
	if (windows.size() == 0) {
	    saveAction.setEnabled(false);
	}
	changeEditOptions();
    }

    /**
      * Starts the application by adding the help menu, 
      * and packing and showing the editor.
      */
    public void start() {
	// ----- Help Menu ------
	JMenu helpMenu = new JMenu("Help");
	JMenuItem aboutItem = new JMenuItem("About...");
	aboutItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    JOptionPane.showMessageDialog(DEditor.this,
						  "Delivery XML Editor\n" +
						  "release: "
						  + releaseDate+"\nLicence:\n"+
										  GPLlicenceText);
		}
	    });
	helpMenu.add(aboutItem);
	deditHelpBroker = deditHelpSet.createHelpBroker();
	JMenuItem getHelpItem = new JMenuItem("Open Help...");
	getHelpItem.addActionListener(new CSH.DisplayHelpFromSource(deditHelpBroker));
	helpMenu.add(getHelpItem);
	menuBar.add(helpMenu);
	// ----- Help Button -----
	ClassLoader cl = DEditor.class.getClassLoader();
	ImageIcon helpIcon 
	    = new ImageIcon(cl.getResource("icons/help.jpg"));
	JButton helpButton = new JButton(helpIcon);
	helpButton.addActionListener(new CSH.DisplayHelpAfterTracking(deditHelpBroker));
	addIconButton(helpButton);
	// -----------------------
	this.pack();
	this.show();
    }

    /**
      * loads a plugin for this DEditor which modifies behavior.  In general, 
      * after making a new DEditor(), load any special plugin right before 
      * calling DEditor.start().  
      * @param plugin the plugin to load into DEditor
      */
    private void loadPlugin(DEditPlugin plugin) {
	if (this.plugin == null) {
	    plugin.load(this,jdp,schema,menuBar);
	    this.plugin = plugin;
	    renderer.setPlugin(plugin);
	} else {
	    JOptionPane.showMessageDialog(this,"Plugin already loaded.");
	}
    }

    /**
      * Returns the HelpSet displayed by the editor so that plugins may add
      * to it.
      * @return the HelpSet displayed by the editor so that plugins may add
      * to it
      */
    public HelpSet getHelpSet() {
	return deditHelpSet;
    }

    /**
      * Returns the HelpBroker displayed by the editor so that plugins may 
      * instruct it to display.
      * @return the HelpBroker displayed by the editor so that plugins may 
      * instruct it to display
      */
    public HelpBroker getHelpBroker() {
	return deditHelpBroker;
    }

    /**
      * Returns whether or not the editor is running as an applet.
      * @return whether or not the editor is running as an applet
      */
    public boolean isApplet() {
	return isApplet;
    }

    /**
      * Called when the user wants to close the application.  It checks to 
      * make sure all open XML documents have been saved, cancelling at 
      * any point if the user decides to do so.
      */
    private void close() {
	while (windows.size() > 0) {
	    DEditXMLWindow window = (DEditXMLWindow)windows.get(0);
	    boolean closed = window.close();
	    if (!closed) {
		return;
	    }
	}
	DEditDebug.println(this,2,"#windows " + windows.size());
	dispose();
	if (!isApplet) {
	    System.exit(0);
	}
    }

    /**
      * Returns a property value, which only has meaning when running as an 
      * applet or when that property has earlier been set with a call to 
      * setProperty().  First this function tries to get a value from a 
      * previous setProperty().  If that fails and it is an applet, it 
      * tries applet.getParameters(propName).
      * @param propName the name of the property to retrieve
      * @return the value for the input property
      */
    public String getProperty(String propName) {
	String value = (String)properties.get(propName);
	if (value != null) {
	    return value;
	}
	if (applet == null || !isApplet) {
	    return null;
	}
	return applet.getParameter(propName);
    }

    /**
      * Sets the value for a property so that the next call to getProperty() 
      * for that property name will return the value input here.  For an 
      * applet this overrides the property of the parameter by name "propName "
      * for that applet
      * @param propName the name of the property to set
      * @param value the value to set property "propName" to
      */
    public void setProperty(String propName, String value) {
	if (propName != null && value != null) {
	    properties.put(propName,value);
	}
    }

    /**
      * Command line invoked function which calls makeDEditor(args,null,null).
      * This means that the DEditPlugin and JApplet are null, but a command 
      * line option may specify the plugin.
      * @param args arguments to be used to configure the DEditor.
      */
    public static void main(String[] args) {
	DEditor dedit = makeDEditor(args,null,null);
	dedit.start();
    }

    /**
      * Makes a DEditor much which can be start()ed and places 
      * extra arguments not used by it into any plugins or
      * into the input vector if there are no plugins
      * @param args command line arguments
      * @param plugin a pre-existing plugin to be used if none is specified 
      * by args
      * @param applet the applet this DEditor should run in (runs as an 
      * application if null)
      * @return a DEditor ready to be used
      */
    public static DEditor makeDEditor(String[] args, DEditPlugin plugin,
				      JApplet applet) {
	//java.util.Properties p = System.getProperties();
	//p.list(System.out);
	String configFile = "configFiles/simpleconfig.xml";
	String schemaFilename = null;
	String xmlFilename = null;
	boolean allowLoadSchema = true;
	boolean allowOpenFile = true;
	boolean allowNewFile = true;
	boolean relaxed = false;
	String pluginName = null;
	Vector extraArgs = new Vector();
	for (int i = 0; i < args.length; i++) {
	    if (args[i].startsWith("xsd=")) {
		schemaFilename = args[i].substring(4);
	    } else if (args[i].startsWith("xml=")) {
		xmlFilename = args[i].substring(4);
	    } else if (args[i].startsWith("cfg=")) {
		configFile = args[i].substring(4);
	    } else if (args[i].startsWith("plugin=")) {
		pluginName = args[i].substring(7);
	    } else if (args[i].equals("-S")) {
		allowLoadSchema = false;
	    } else if (args[i].equals("-N")) {
		allowNewFile = false;
	    } else if (args[i].equals("-X")) {
		allowOpenFile = false;
	    } else if (args[i].equals("-V")) {
		relaxed = true;
	    } else if (args[i].equals("-pl")) {
		int priority = 0;
		i++;
		if (i < args.length) {
		    priority = Integer.parseInt(args[i]);
		}
		DEditDebug.setPriority(priority);
	    } else {
		extraArgs.add(args[i]);
	    }
	}
	//String dirString = System.getProperty("user.dir");
	String dirString = "";
	if (pluginName != null) {
	    try {
		String[] leftOvers = new String[extraArgs.size()];
		extraArgs.toArray(leftOvers);
		extraArgs.removeAllElements();
		Class c = Class.forName(pluginName);
		if (c != null) {
		    String[] temp = new String[0];
		    Method m=c.getDeclaredMethod("buildPlugin",
						 new Class[]{temp.getClass()});
		    if (m != null) {
			plugin =(DEditPlugin)m.invoke(null, 
						      new Object[]{leftOvers});
		    }
		}
	    } catch (Exception e) {
		System.err.println("load exception error: " + e.toString());
	    }
	}
	DEditor dedit = new DEditor(dirString,configFile,xmlFilename,
				    schemaFilename,allowNewFile,
				    allowOpenFile,allowLoadSchema,relaxed,
				    plugin,applet);
	return dedit;
    }
}
