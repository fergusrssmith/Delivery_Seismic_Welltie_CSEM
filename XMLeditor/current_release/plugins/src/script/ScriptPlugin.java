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
package script;

import javax.swing.tree.DefaultTreeModel;
import javax.swing.event.CellEditorListener;
import javax.swing.event.ChangeEvent;
import javax.swing.table.*;
import javax.swing.border.EmptyBorder;
import javax.swing.plaf.basic.BasicComboBoxEditor;
import javax.swing.*;
import org.w3c.dom.*;
import java.awt.Component;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;
import javax.help.*;
import java.net.URL;
import java.io.*;
import model.*;

/**
  * Class which handles GUI related aspects of the Script Plugin, and 
  * passes off Script generation/validation to ScriptGenerator
  */
public class ScriptPlugin implements DEditPlugin {

    private static final int FREE = 100;
    private static final int RESTRICTED = 101;
    private String dirString = "";
    private String fileString = "";
    private String optString = "";
    private File curSHDir = null;
    private int mode = FREE;
    private DEditor editor;
    private JDesktopPane jdp;
    private ScriptGenerator sg = null;
    private HelpSet scriptHelpSet = null;
    private static final String GENERATE_ID = "generate_script";
    private static final String EXECUTE_ID = "execute_script";
    private static final String PARAMETERS_ID = "editparameters";
    private String TREE_ID = null;
    private boolean isApplet = false;

    /**
      * Private constructor used by buildPlugin (a method users should 
      * use to construct any plugin).  
      * @param curDirString a String representing the current directory, 
      * as per System.getProperty("user.dir").
      * @param extraArgs an array of arguments, generally those left over by
      * DEditor.main();
      */
    private ScriptPlugin(String curDirString, String[] extraArgs) {
	this.curSHDir = new File(curDirString);
	this.mode = FREE;
	String helpSetName = "plugins/help/script/ScriptHelp.hs";
	try {
	    ClassLoader cl = this.getClass().getClassLoader();
	    URL hsURL = HelpSet.findHelpSet(cl,
					    helpSetName);
	    scriptHelpSet = new HelpSet(cl, hsURL);
	    TREE_ID = "scriptintro";
	} catch (Exception e) {
	    System.err.println("Helpset " + helpSetName 
			       + " not found.");
	}
	for (int i = 0; i < extraArgs.length; i++) {
	    if (extraArgs[i].equals("-r")) {
		this.mode = RESTRICTED;
	    } else if (extraArgs[i].startsWith("script-hs=")) {
		String subHSName = extraArgs[i].substring(10);
		DEditDebug.println(this,6,"HS:"+subHSName);
		try {
		    ClassLoader cl = this.getClass().getClassLoader();
		    URL shsURL = HelpSet.findHelpSet(cl, subHSName);
		    if (shsURL == null) {
			System.out.println("applet=>screwed");
			File f = new File(subHSName);
			shsURL = f.toURL();
			if (shsURL == null) {
			    DEditDebug.println(this,6,"dammit");
			} else {
			    HelpSet subHS = new HelpSet(cl, shsURL);
			    scriptHelpSet.add(subHS);
			    TREE_ID = "special_script_intro";
			    DEditDebug.println(this,6,"...added" + subHS);
			}
		    } else {
			HelpSet subHS = new HelpSet(null, shsURL);
			scriptHelpSet.add(subHS);
			TREE_ID = "special_script_intro";
			DEditDebug.println(this,6,"...added" + subHS);
		    }
		} catch (Exception e) {
		    System.err.println("Helpset " +  subHSName 
				       + " not found.");
		    e.printStackTrace();
		}
	    } else {
		DEditDebug.println(this,5,"useless arg " + extraArgs[i]);
	    }
	}
    }

    /**
      * Loads the sbub plugin.
      * @param deditor the DEditor that will receive this plugin
      * @param jdp the JDesktop pane that the windows appear in, should the 
      * user want to put in any more JInternalFrames, like the script text 
      * frame in this case
      * @param schema the DEditSchema for editor
      * @param menuBar the JMenuBar used by editor
      */
    public void load(final DEditor deditor, final JDesktopPane jdp, 
		     DEditSchema schema, JMenuBar menuBar) {
	this.editor = deditor;
	this.jdp = jdp;
	this.sg = new ScriptGenerator(deditor);
	ClassLoader cl = ScriptPlugin.class.getClassLoader();
	if (scriptHelpSet != null) {
	    HelpSet deditHelpSet = deditor.getHelpSet();
	    java.util.Enumeration enum = scriptHelpSet.getHelpSets();
	    boolean setSpecial = false;
	    while (enum.hasMoreElements()) {
		HelpSet child = (HelpSet)enum.nextElement();
		boolean check = scriptHelpSet.remove(child);
		DEditDebug.println(this,3,"" + check);
		deditHelpSet.add(child);
		setSpecial = true;
	    }
	    deditHelpSet.add(scriptHelpSet);
	    if (setSpecial) {
		deditHelpSet.setHomeID("special_script_intro");
	    } else {
		deditHelpSet.setHomeID("scriptintro");
	    }
	}
	JMenu scriptMenu = new JMenu("Scripts");
	JMenuItem showScriptItem = new JMenuItem("View Script...");
	JMenuItem generateScriptItem = new JMenuItem("Generate Script...");
	JMenuItem executeScriptItem = new JMenuItem("Execute Script...");
	showScriptItem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    String filename = chooseFile(false); 
		    if (filename == null) {
			return;
		    }
		    String file = "";
		    if (deditor.isApplet()) {
			file = deditor.readRemoteFile(filename);
		    } else {
			try {
			    file = getFileAsString(filename);
			} catch (MalformedScriptException mse) {
			    JOptionPane.showMessageDialog(deditor,
							  mse.toString());
			}
		    }
		    JTextArea scriptArea = new JTextArea(33,40);
		    scriptArea.setEditable(false);
		    JScrollPane scriptView = new JScrollPane(scriptArea);
		    JInternalFrame scriptFrame 
			= new JInternalFrame(filename,true,true,true,true);
		    scriptFrame.getContentPane().add(scriptView);
		    scriptFrame.setLocation(450,0);
		    scriptArea.setText(file);
		    jdp.add(scriptFrame);
		    scriptFrame.pack();
		    scriptFrame.show();
		}
	    });
	Action generateAction = new AbstractAction() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow curr = deditor.getSelectedWindow();
		    if (curr == null) {
			JOptionPane.showMessageDialog(deditor,
						      "Select a window "
						      +"first.");
			return;
		    }
		    final DEditNode root = curr.getRoot();
		    String openFilename = curr.getOpenFilename();
		    int lastDir = Math.max(0,openFilename.lastIndexOf("/"));
		    String xmlDir = openFilename.substring(0,lastDir+1);
		    generateAndExecuteScriptFile(deditor,xmlDir,root.getNode(),
						 false);
		}
	    };
	generateScriptItem.addActionListener(generateAction);
	ImageIcon generateIcon
	    = new ImageIcon(cl.getResource("icons/generate.jpg"));
	JButton generateButton = new JButton(generateAction);
	CSH.setHelpIDString(generateButton,GENERATE_ID);
	generateButton.setIcon(generateIcon);
	generateButton.setToolTipText("Generate script(s)");
	editor.addIconButton(generateButton);
	Action executeAction = new AbstractAction() {
		public void actionPerformed(ActionEvent ae) {
		    DEditXMLWindow curr = editor.getSelectedWindow();
		    if (curr == null) {
			JOptionPane.showMessageDialog(editor,
						      "Select a window "
						      +"first.");
			return;
		    }
		    final DEditNode root = curr.getRoot();
		    // get savefile
		    String openFilename = curr.getOpenFilename();
		    int lastDir = Math.max(0,openFilename.lastIndexOf("/"));
		    String xmlDir = openFilename.substring(0,lastDir+1);
		    String filename 
			= generateAndExecuteScriptFile(deditor,xmlDir,
						       root.getNode(),true);
		    if (filename == null) {
			return;
		    }
		    if (editor.isApplet()) {
			String msg = "Script "+filename+".sh finished,\n "
			    +"output going to "+filename+"_out.txt and "
			    +filename+"_err.txt.";
			JOptionPane.showMessageDialog(deditor,msg);
		    }
		}
	    };
	executeScriptItem.addActionListener(executeAction);
	ImageIcon executeIcon
	    = new ImageIcon(cl.getResource("icons/execute.jpg"));
	JButton executeButton = new JButton(executeAction);
	CSH.setHelpIDString(executeButton,EXECUTE_ID);
	executeButton.setIcon(executeIcon);
	executeButton.setToolTipText("Generate and execute script(s)");
	editor.addIconButton(executeButton);
	//executeScriptItem.setEnabled(false);
	scriptMenu.add(showScriptItem);
	scriptMenu.add(generateScriptItem);
	scriptMenu.add(executeScriptItem);
	menuBar.add(scriptMenu);
    }

    /**
      * Returns the JavaHelp Map ID which should be used for the XML tree.
      * @return the JavaHelp Map ID which should be used for the XML tree
      */
    public String getTreeHelpID() {
	return TREE_ID;
    }

    /**
      * Returns the contents of filename (which should be a script file but it 
      * doesn't really mattter) as a String.
      * @param filename the file whose contents will be returned
      * @return filename's contents as a String
      */
    public static String getFileAsString(String filename) 
	throws MalformedScriptException {
	String scriptString = null;
	try {
	    FileReader fr = new FileReader(filename);
	    StringWriter out = new StringWriter();
	    int c;
	    
	    while ((c = fr.read()) != -1)
		out.write(c);
	    
	    fr.close();
	    scriptString = out.toString();
	} catch (IOException ioe) {
	    throw new MalformedScriptException("Locate script " + filename
					       + " error: " + ioe.toString());
	}
	return scriptString;
    }

    /**
      * Convenience method which creates a Save File Dialog or an Open 
      * File Dialog for selecting a file, then returns the file name 
      * as a String or null if the "APPROVE_OPTION" was not chosen.
      * @param save true if it should show a Save File Dialog, false for 
      * Open File Dialog
      * @return the name of the file to save to or open
      */
    private String chooseFile(boolean save) {
	String prefix = "Open";
	if (save) {
	    prefix = "Save";
	}
	if (editor.isApplet()) {
	    String msg = "Specify remote URL";
	    RemoteFileChooser rfc 
		= new RemoteFileChooser(editor,prefix + " Script File",
					"ScriptFileURL", prefix);
	    int reFchooser = rfc.showDialog();
	    String filename = null;
	    if (reFchooser == JFileChooser.APPROVE_OPTION) {
		String pathString = rfc.getPath();
		String name = rfc.getName();
		filename = pathString + "/" + name;
		System.out.println("(" + name + ")OOpenSession " + filename);
		return filename;
	    } else {
		return null;
	    }
	}
	JFileChooser chooser = new JFileChooser(curSHDir);
	chooser.setDialogTitle(prefix + " Script File");
	chooser.setAcceptAllFileFilterUsed(true);
	DEditSimpleFileFilter filter 
	    = new DEditSimpleFileFilter(".sh","Script Files");
	chooser.setFileFilter(filter);
	String filename = null;
	int returnVal;
	if (save) {
	    returnVal = chooser.showSaveDialog(editor);
	} else {
	    returnVal = chooser.showOpenDialog(editor);
	}
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	    File f = chooser.getSelectedFile();
	    curSHDir = f.getParentFile();
	    filename = f.getAbsolutePath();
	    return filename;
	} else {
	    return null;
	}
    }

    /**
      * Generates and executes a script file as per directions specified on a
      * JDialog it provides.  Returns the name of the top-level script 
      * generated.
      * @param deditor the DEditor to display the dialog in
      * @param root the root Node (under the Document node) for the XML 
      * tree to be turned into a script
      * @return the name of the top-level script generated
      */
    private String generateAndExecuteScriptFile(final DEditor deditor,
						final String xmlDir,
						final Node root,
						final boolean execute) {
	String dialogTitle = "Generate ";
	if (execute) {
	    dialogTitle += "And Exectute ";
	}
	dialogTitle += "Script File(s)";
	final JDialog dialog = new JDialog(deditor,dialogTitle,true);
	JPanel runOptsPanel = new JPanel(new GridLayout(3,1));
	final JCheckBox forceLinear = new JCheckBox("Force linear run");
	final JTextField dirName = new JTextField(dirString,30);
	final JTextField scriptBase = new JTextField(fileString,15);
	final JTextField options = new JTextField(optString,40);
	forceLinear.setSelected(true);
	forceLinear.setEnabled(true);
	runOptsPanel.add(forceLinear);
	runOptsPanel.add(new JLabel("Command-line Options for bsub"));
	runOptsPanel.add(options);
	String buttonText = "Generate";
	if (execute) {
	    buttonText = "Execute";
	}
	JButton executeButton = new JButton(buttonText);
	executeButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    String dirBaseName = dirName.getText();
		    String fileName = scriptBase.getText();
		    setScriptLocation(dirBaseName,fileName);
		    try {
			generateScript(dirBaseName + "/",
				       fileName,xmlDir,root);
		    } catch (MalformedScriptException mse) {
			dialog.dispose();
			JOptionPane.showMessageDialog(deditor,mse.toString());
			return;
		    }
		    String baseName = dirBaseName + "/" + fileName;
		    File bsubVersion = new File(baseName + ".bsub.sh");
		    String opts = options.getText();
		    setOptString(opts);
		    String bsubScriptName = baseName + ".bsub.sh";
		    String bsubOptName = baseName + ".bsub.opt.sh";
		    if (!opts.equals("")) {
			try {
			    if (editor.isApplet()) {
				String contents = "#!/bin/sh\nsh "+baseName
				    +".bsub.sh "+opts+"\n";
				editor.writeRemoteFile(bsubOptName,contents);
				sg.chmod(bsubOptName,"775");
			    } else {
				FileWriter fw = new FileWriter(bsubOptName);
				fw.write("#!/bin/sh\n");
				fw.write("sh " + bsubVersion.getAbsolutePath()
					 + " " + opts + "\n");
				fw.close();
				sg.chmod(bsubOptName,"775");
			    }
			    bsubScriptName = bsubOptName;
			} catch (Exception e) {
			    String err = "Error writing options included"+
				" bsub script.\n";
			    JOptionPane.showMessageDialog(deditor,
							  err+e.toString());
			}
		    }
		    if (execute) {
			if (forceLinear.isSelected() 
			    || (!editor.isApplet() && !bsubVersion.exists())) {
			    if (editor.isApplet()) {
				String msg = "Script "+baseName+".sh started,"
				    +"\n output going to "+baseName+"_out.txt"
				    +" and "+baseName+"_err.txt.";
				JOptionPane.showMessageDialog(editor,msg);
			    }
			    runScript(baseName+".sh > "+baseName+"_out.txt 2>"
				      +baseName+"_err.txt",
				      false);
			} else {
			    if (editor.isApplet()) {
				String msg = "Script "+bsubScriptName
				    +" started,\n output going to "+baseName
				    +"_out.txt and "+baseName+"_err.txt.";
				JOptionPane.showMessageDialog(editor,msg);
			    }
			    runScript(bsubScriptName+" > "+baseName
				      +"_out.txt 2>"+baseName+"_err.txt",
				      true);
			}
		    }
		    dialog.dispose();
		    displayWarnings(getWarnings());
		}
	    });
	JPanel mainPanel = buildMainGeneratePanel(dialog,executeButton,dirName,
						  scriptBase);
	mainPanel.add(runOptsPanel,BorderLayout.CENTER);
	dialog.getContentPane().add(mainPanel);
	dialog.pack();
	dialog.setLocationRelativeTo(editor);
	dialog.show();
	String scriptName = dirName.getText() + "/" + scriptBase.getText();
	if (scriptName.equals("/")) {
	    return null;
	}
	return scriptName;
    }

    /**
      * Helper whose only utility is letting ScriptGenerator sg not be final 
      * so I can construct it in load() instead of the constructor.
      * @param dirBase the directory base for the scripts to be generated in
      * @param filename the base name for all scripts and job names
      * @param xmlDir the directory the xml that is being used to make the 
      * scripts resides in (used for includes)
      * @param root the room Node of the XML tree
      */
    private void generateScript(String dirBase, String filename,
				String xmlDir, Node root) 
	throws MalformedScriptException {
	if (sg != null) {
	    sg.generateScript(dirBase,filename,xmlDir,root);
	}
    }

    /**
      * Gets the list of warnings generated by the ScriptGenerator on the 
      * most recent generateScript() call.
      * @return a Vector of String warnings generated by the ScriptGenerator
      * on the most recent generateScript() call
      */
    private Vector getWarnings() {
	if (sg != null) {
	    return sg.getWarnings();
	} else {
	    return new Vector();
	}
    }

    /**
      * Another helper to save inner class "final" problems.  Sets the 
      * value stored to use for starting text in bsub execute options field.
      * @param newOpts the value to set the optString to.  optString is
      * the String that holds the last input options for bsub arguments.
      */
    private void setOptString(String newOpts) {
	optString = newOpts;
    }

    /**
      * Displays a message dialog containing all the Strings in warnings
      * separated by newlines.
      * @param warnings a Vector of Strings which are warning messages.
      */
    private void displayWarnings(Vector warnings) {
	if (warnings == null || warnings.size() == 0) {
	    return;
	}
	String message = "The following warnings were reported during script "
	    + "generation:\n";
	for (int i = 0; i < warnings.size(); i++) {
	    message += "  - " + (String)(warnings.get(i)) + "\n";
	}
	JOptionPane.showMessageDialog(editor,message);
    }

    /**
      * Helper function which build and manages the components common to the 
      * generate and execute panels.
      * @param dialog the JDialog to add components to, but which we need 
      * here in order to tell it to dispose in certain ActionListener's
      * @param acceptButton the button which will launch the "ACCEPT" type 
      * action 
      * @param dirname the JTextField holding the directory name
      * @param scriptBase the JTextField holding the base name for the script
      * (also the base name for submitted jobs)
      * @return a JPanel to be added to the main dialog
      */
    private JPanel buildMainGeneratePanel(final JDialog dialog,
					  JButton acceptButton, 
					  final JTextField dirName,
					  final JTextField scriptBase) {
	JPanel mainPanel = new JPanel(new BorderLayout());
	JButton browseButton = new JButton("Browse...");
	browseButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    String filename = chooseFile(true);
		    if (filename == null) {
			return;
		    }
		    String scriptBaseStr = null;
		    int index = filename.lastIndexOf("/"); // "\" ???
		    if (index != -1) {
			dirName.setText(filename.substring(0,index));
			scriptBaseStr = filename.substring(index+1);
		    } else {
			dirName.setText("");
			scriptBaseStr = filename;
		    }
		    index = scriptBaseStr.lastIndexOf(".sh");
		    if (index != -1 && index == scriptBaseStr.length() - 3) {
			scriptBaseStr = scriptBaseStr.substring(0,index);
		    }
		    scriptBase.setText(scriptBaseStr);
		}
	    });
	JPanel choosePanel = new JPanel(new BorderLayout());
	JPanel dirPanel = new JPanel();
	dirPanel.add(new JLabel("Directory Name"));
	dirPanel.add(dirName);
	choosePanel.add(dirPanel,BorderLayout.NORTH);
	JPanel baseBrowsePanel = new JPanel();
	baseBrowsePanel.add(new JLabel("Script Base Name"));
	baseBrowsePanel.add(scriptBase);
	baseBrowsePanel.add(browseButton);
	choosePanel.add(baseBrowsePanel,BorderLayout.SOUTH);
	JButton cancelButton = new JButton("Cancel");
	cancelButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    setScriptLocation(dirName.getText(),scriptBase.getText());
		    dialog.dispose();
		}
	    });
	JPanel buttonPanel = new JPanel();
	buttonPanel.add(acceptButton);
	buttonPanel.add(cancelButton);
	mainPanel.add(choosePanel,BorderLayout.NORTH);
	mainPanel.add(buttonPanel,BorderLayout.SOUTH);
	return mainPanel;
    }

    /**
      * Another helper to avoid inner class "final" issues.  Sets the starting 
      * location of scripts and base name.
      * @param dir the String to use for script location
      * @param file the file base name for scripts and job names
      */
    private void setScriptLocation(String dir, String file) {
	dirString = dir;
	fileString = file;
    }

    /**
      * Runs a script called command. 
      * @param command the name of the script to run
      * @param bsub true if there are bsub calls in the script.  This 
      * will then launch a "bkill 0" (kill all jobs for this user) command 
      * when "Kill Jobs" is selected.
      */
    private void runScript(String command, final boolean bsub) {
	if (editor != null && editor.isApplet()) {
	    editor.executeRemote(command);
	    return;
	}
	try {
	    Runtime rt = Runtime.getRuntime();
	    // chmod
	    // run
	    JTextArea outputArea = new JTextArea(33,40);
	    JScrollPane outputView = new JScrollPane(outputArea);
	    final JInternalFrame outputFrame 
		= new JInternalFrame("Standard Error for " + command,true,
				     true,true,true);
	    JPanel jp = new JPanel(new BorderLayout());
	    jp.add(outputView,BorderLayout.CENTER);
	    JButton dismissButton = new JButton("Kill Jobs");
	    final Process p = rt.exec("sh " + command);
	    dismissButton.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent ae) {
			if (bsub) {
			    try {
				Runtime rt = Runtime.getRuntime();
				Process kill = rt.exec("bkill 0");
			    } catch (IOException ioe) {
				String err = "Couldn't execute script: " 
				    + ioe.toString();
				JOptionPane.showMessageDialog(editor,err);
			    }
			}
			p.destroy();
		    }
		});
	    jp.add(dismissButton,BorderLayout.SOUTH);
	    outputFrame.getContentPane().add(jp);
	    outputFrame.setLocation(450,0);
	    jdp.add(outputFrame);
	    outputFrame.pack();
	    outputFrame.show();
	    if (p != null) {
		StreamGobbler errorGobbler 
		    = new StreamGobbler(p.getErrorStream(),"Error",
					null/*System.err*/,"",outputArea);
		StreamGobbler outputGobbler 
		    = new StreamGobbler(p.getInputStream(),"Output", 
					null/*System.out*/,"",outputArea);
		errorGobbler.start();
		outputGobbler.start();
	    }
	} catch (IOException ioe) {
	    JOptionPane.showMessageDialog(editor,"Couldn't execute script: " 
					  + ioe.toString());
	}
    }
    
    /**
      * Plugin potential override / add-on for editing nodes.  Returns 
      * true iff no more editing should occur after this function.  If it 
      * returns false the standard editing dialog will be displayed
      * @param selected the DEditNode that might be edited
      * @return true if editing should stop after this function
      */
    public boolean edit(DEditNode selected) {
	String name = selected.getNode().getNodeName();
	// !!! technically this only works flawlessly if there are no other
	// nodes called parameters or parameter -- could check parents here
	// to fix problem should the need arise
	DEditXMLWindow window = editor.getSelectedWindow();
	if (window == null) {
	    return true; // nothing to do
	}
	if (name.equals("parameters")) {
	    NodeList paramList = selected.getNode().getChildNodes();
	    Vector paramVec = new Vector();
	    for (int i = 0; i < paramList.getLength(); i++) {
		Node n = paramList.item(i);
		if (n.getNodeType() == Node.ELEMENT_NODE) {
		    paramVec.add(n);
		}
	    }
	    if (paramVec.size() == 0) {
		return true;
	    }
	    Node[] params = new Node[paramVec.size()];
	    paramVec.toArray(params);
	    editParameters(params);
	    window.nodeStructureChanged(selected);
	    window.update(true);
	    return true;
	} else if (name.equals("parameter")) {
	    editParameters(new Node[]{selected.getNode()});
	    window.nodeStructureChanged(selected);
	    window.update(true);
	    return true;
	} else if (mode == RESTRICTED) { // shouldn't ever be true
	    return true; // uneditable
	} else {
	    return false; // go about your business
	}
    }

    /**
      * Returns whether or not the input node should be editable.
      * @param selected the DEditNode in question for editable / not editable
      * @param defaultEnabled true if Schema says it should be editable (i.e., 
      * has parameters or is a leaf element)
      * @return true if selected DEditNode should be editable
      */
    public boolean isEditable(DEditNode selected, boolean defaultEnabled) {
	String name = selected.getNode().getNodeName();
	return (name.equals("parameters") || name.equals("parameter")
		|| (mode != RESTRICTED && defaultEnabled));
    }

    /**
      * Class which manages the table model for the JTable in the "Edit 
      * Parameter(s)" dialog.  Only really useful for isCellEditable so
      * maybe switch back to extending DefaultTableModel.  Switched during 
      * ongoing efforts to fix JTable editing problem (cancels on change 
      * focus).
      */
    private class ParameterTableModel extends AbstractTableModel {
	private Object[][] rowData;
	private String[] columnNames;
	private Node[][] nodes;

	public ParameterTableModel(Object[][] data, String[] names,
				   Node[][] nodes) {
	    this.rowData = data;
	    this.columnNames = names;
	    this.nodes = nodes;
	}
	public String getColumnName(int col) { 
	    return columnNames[col].toString(); 
	}
	public int getRowCount() { return rowData.length; }
	public int getColumnCount() { return columnNames.length; }
	public Object getValueAt(int row, int col) { 
	    return rowData[row][col]; 
	}
	public boolean isCellEditable(int row, int column) {
	    return (nodes[row][column] != null 
		    && (mode != RESTRICTED || column > 0));
	}
	public void setValueAt(Object value, int row, int col) {
		    rowData[row][col] = value;
		    fireTableCellUpdated(row, col);
	}
    }

    /**
      * Helper for edit which displays a JDialog for assigning an array of 
      * parameters.
      * @param paramNodes the Nodes corresponding to parameters
      */
    private void editParameters(final Node[] paramNodes) {
	
	final String[][] params = new String[paramNodes.length][2];
	final Node[][] temp = new Node[paramNodes.length][2];
	// fill from JTable
	for (int i = 0; i < paramNodes.length; i++) {
	    Node nameNode = null;// = paramNodes[i].getFirstChild();
	    NodeList paramChilds = paramNodes[i].getChildNodes();
	    int j = 0;
	    for (; j < paramChilds.getLength(); j++) {
		Node n = paramChilds.item(j);
		if (n.getNodeType() == Node.ELEMENT_NODE 
		    && n.getNodeName().equals("name")) {
		    nameNode = n;
		    j++;
		    break;
		}
	    }
	    if (nameNode == null) {
		continue;
	    }
	    Node nameTextNode = nameNode.getFirstChild();
	    temp[i][0] = nameTextNode;
	    params[i][0] = nameTextNode.getNodeValue();
	    Node valueNode = null;
	    for (; j < paramChilds.getLength(); j++) {
		Node n = paramChilds.item(j);
		if (n.getNodeType() == Node.ELEMENT_NODE 
		    && n.getNodeName().equals("value")) {
		    valueNode = n;
		    break;
		}
	    }
	    if (valueNode == null) {
		continue;
	    }
	    NamedNodeMap nnm = valueNode.getAttributes();
	    String subVal = null;
	    if (nnm != null) {
		Node valueTextNode = nnm.getNamedItem("substitutionString");
		if (valueTextNode != null) {
		    subVal = valueNode.getNodeValue();
		}
	    }
	    if (subVal == null || subVal.equals("")) {
		Node valueTextNode = valueNode.getFirstChild();
		if (valueTextNode == null) {
		    Document d = valueNode.getOwnerDocument();
		    valueTextNode = d.createTextNode("");
		    valueNode.appendChild(valueTextNode);
		}
		temp[i][1] = valueTextNode;
		params[i][1] = valueTextNode.getNodeValue();
	    } else {
		params[i][1] = subVal;
		temp[i][1] = null;
	    }
	}
	final JDialog dialog = new JDialog(editor,"Edit Parameter(s)",false);
	JPanel paramPanel = new JPanel(new BorderLayout());
	paramPanel.setBorder(new EmptyBorder(5,5,5,5));
	JPanel tablePanel = new JPanel(new BorderLayout());
	JLabel tableLabel = new JLabel("Parameters");
	tablePanel.add(tableLabel,BorderLayout.NORTH);
	String[] names = new String[]{"Name","Value"};
	ParameterTableModel ptm = new ParameterTableModel(params,names,temp);
	final JTable paramTable = new JTable(ptm);
	//paramTable.setIntercellSpacing(new Dimension(5,5));
	CellEditorListener cel = new CellEditorListener() {
		public void editingCanceled(ChangeEvent e) {
		    DEditDebug.println(this,4,"cancelled");
		}
		public void editingStopped(ChangeEvent e) {
		    DEditDebug.println(this,4,"stopped");
		    updateParameterInfo(temp,paramTable);
		} 
	    };
	JComboBox jcb = new JComboBox();
	jcb.setEditable(true);
	//ComboBoxEditor tempC = jcb.getEditor();
	//System.out.println(tempC);
	jcb.setEditor(new NoBoldComboBoxEditor());
	JTextField jtf = new JTextField();
	NoBoldTableCellEditor dce = new NoBoldTableCellEditor(jcb);
	dce.addCellEditorListener(cel);
	TableColumn nameCells = paramTable.getColumnModel().getColumn(0);
	TableColumn valueCells = paramTable.getColumnModel().getColumn(1);
	nameCells.setCellEditor(dce);
	valueCells.setCellEditor(dce);
	final String[][] backup = new String[params.length][params[0].length];
	for (int i = 0; i < params.length; i++) {
	    for (int j = 0; j < params[i].length; j++) {
		backup[i][j] = params[i][j];
	    }
	}
	JScrollPane jsp = new JScrollPane(paramTable);	
	jsp.setPreferredSize(new Dimension(350,300));
	tablePanel.add(jsp,BorderLayout.CENTER);
	paramPanel.add(tablePanel,BorderLayout.CENTER);
	JPanel buttonPanel = new JPanel();
	JButton setButton = new JButton("Accept");
	setButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditDebug.println(this,4,"accept");
		    updateParameterInfo(temp,paramTable);
		    dialog.dispose();
		}
	    });
	buttonPanel.add(setButton);
	JButton cancelButton = new JButton("Cancel");
	cancelButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    for (int i = 0; i < backup.length; i++) {
			for (int j = 0; j < backup[i].length; j++) {
			    paramTable.setValueAt(backup[i][j],i,j);
			}
		    }
		    updateParameterInfo(temp,paramTable);
		    dialog.dispose();
		}
	    });
	buttonPanel.add(cancelButton);
	JButton clearValuesButton = new JButton("Clear Values");
	clearValuesButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    for (int i = 0; i < paramNodes.length; i++) {
			paramTable.setValueAt("",i,1);
		    }
		}
	    });
	buttonPanel.add(clearValuesButton);
	JButton helpButton = new JButton("Help");
	helpButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    HelpBroker hb = editor.getHelpBroker();
		    dialog.setModal(false);
		    if (hb != null) {
			try {
			    hb.setCurrentID(PARAMETERS_ID);
			    hb.setDisplayed(true);
			} catch (Exception e) {
			    DEditDebug.println(this,6,"couldn't display help");
			}
		    }
		}
	    });
	buttonPanel.add(helpButton);
	paramPanel.add(buttonPanel,BorderLayout.SOUTH);
	dialog.getContentPane().add(paramPanel);
	dialog.pack();
	dialog.setLocationRelativeTo(editor);
	dialog.show();
    }

    /**
      * Updates the information for the parameters being edited.  As this 
      * is called on focus change, it is also called with backup values
      * reset on clicking "Cancel".
      * @param temp the nodes whose values must be updated (an array of 
      * arbitrary length containing arrays of length two with element 0 
      * being the node for the name of the parameter and element 1 being the 
      * node with the value of the parameter)
      * @param paramTable the JTable to get the values from
      */
    private void updateParameterInfo(Node[][] temp, JTable paramTable) {
	for (int i = 0; i < temp.length; i++) {
	    if (temp[i][1] == null) {
		continue;
	    }
	    String name = (String)paramTable.getValueAt(i,0);
	    temp[i][0].setNodeValue(name);
	    String value = (String)paramTable.getValueAt(i,1);
	    temp[i][1].setNodeValue(value);
	    DEditDebug.println(this,3,name+"="+value);
	}
    }

    /**
      * Puts restrictions on what is cuttable above and beyond Schema 
      * restrictions.  Returns true when it should be cuttable if the Schema 
      * allows it.
      * @param selected the node in question for cuttable / uncuttable
      * @return true if it should be cuttable if allowed by schema
      */
    public boolean isCuttable(DEditNode selected) {
	return (mode != RESTRICTED); 
		//|| selected.getNode().getNodeName().equals("parameter"));
    }

    /**
      * Puts restrictions on what is copyable above and beyond Schema 
      * restrictions.  Returns true when it should be copyable if the Schema 
      * allows it.
      * @param selected the node in question for copyable / uncopyable
      * @return true if it should be copyable if allowed by schema
      */
    public boolean isCopyable(DEditNode selected) {
	return (mode != RESTRICTED); 
	//|| selected.getNode().getNodeName().equals("parameter"));
    }

    /**
      * Puts restrictions on where new tags may be be inserted above and 
      * beyond Schema restrictions.  The input DEditNode is the node that 
      * the new node will be inserted before.
      * @param selected the node queried to see if insert before it is 
      * ever allowed
      * @return true if given that it is valid under the Schema, a node 
      * may be inserted before this node
      */
    public boolean newTagBeforeAllowed(DEditNode selected) {
	return (mode != RESTRICTED); 
	//	|| selected.getNode().getNodeName().equals("parameter"));
    }

    /**
      * Puts restrictions on where new tags may be be inserted above and 
      * beyond Schema restrictions.  The input DEditNode is the node that 
      * the new node will be inserted as a child of.
      * @param selected the node queried to see if insert as its child is 
      * ever allowed
      * @return true if given that it is valid under the Schema, a node 
      * may be inserted as a child of this node
      */
    public boolean newChildTagAllowed(DEditNode selected) {
	return (mode != RESTRICTED); 
	//	|| selected.getNode().getNodeName().equals("parameters"));
    }

    /**
      * Puts restrictions on where new tags may be be inserted above and 
      * beyond Schema restrictions.  The input DEditNode is the node that 
      * the new node will be inserted after.
      * @param selected the node queried to see if insert after it is 
      * ever allowed
      * @return true if given that it is valid under the Schema, a node 
      * may be inserted after this node
      */
    public boolean newTagAfterAllowed(DEditNode selected) {
	return (mode != RESTRICTED); 
	//|| selected.getNode().getNodeName().equals("parameter"));
    }

    /**
      * Puts restrictions on what tags may be switched to other types of 
      * tags valid in the same place above.
      * @param selected the node in question for switchable / unswitchable
      * @return true if node is ever switchable
      */
    public boolean switchTagToAllowed(DEditNode selected) {
	return (mode != RESTRICTED); 
	//	|| selected.getNode().getNodeName().equals("parameter"));
    }

    /**
      * Puts restrictions on whether or not a particular node should 
      * allow something to be pasted in as a child of the DEditNode passed in
      * if the Schema allows it
      * @param selected the node for which to decide if pasting a node as 
      * a child should be allowed if the Schema allows it
      * @return true if node can ever accept children to be pasted
      */
    public boolean pasteChildAllowed(DEditNode selected) {
	return (mode != RESTRICTED); 
	//	|| selected.getNode().getNodeName().equals("parameters"));
    }

    /**
      * Puts restrictions on whether or not a particular node should 
      * allow something to be pasted in before the DEditNode passed in
      * if the Schema allows it
      * @param selected the node for which to decide if pasting a node before 
      * it should be allowed if the Schema allows it
      * @return true if node can ever accept a node to be pasted before it
      */
    public boolean pasteBeforeAllowed(DEditNode selected) {
	return (mode != RESTRICTED); 
	//|| selected.getNode().getNodeName().equals("parameter"));
    }

    /**
      * Puts restrictions on whether or not a particular node should 
      * allow something to be pasted in after the DEditNode passed in
      * if the Schema allows it
      * @param selected the node for which to decide if pasting a node after 
      * it should be allowed if the Schema allows it
      * @return true if node can ever accept a node to be pasted after it
      */
    public boolean pasteAfterAllowed(DEditNode selected) {
	return (mode != RESTRICTED); 
	//|| selected.getNode().getNodeName().equals("parameter"));
    }
    
    // the class below was largely copied from the INT viewer
    private class StreamGobbler extends Thread {
        private InputStream _sgInput;
        private String _sgType;
        private OutputStream _sgOutput;
        private String _initMessage;
	private JTextArea _jta = null;

        public StreamGobbler(InputStream ins, String type, String initMsg) {
            this(ins, type, null, initMsg);
        }

        public StreamGobbler(InputStream ins, String type, OutputStream out,
			     String initMsg) {
	    this(ins,type,out,initMsg,null);
	}
	public StreamGobbler(InputStream ins, String type, OutputStream out,
			     String initMsg, JTextArea jta) {
            _sgInput = ins;
            _sgType = type;
            _sgOutput = out;
            _initMessage = initMsg + "\n\n";
	    _jta = jta;
        }


        public void run() {
            try {
                BufferedOutputStream output = null;
                if (_sgOutput != null) {
                    output = new BufferedOutputStream(_sgOutput);
                    output.write(_initMessage.getBytes());
                }
                BufferedInputStream input = new BufferedInputStream(_sgInput);
                byte[] readResult = new byte[1024];
                int currentReadLen = 0;
                currentReadLen = input.read(readResult, 0, 1024);
                while(currentReadLen > -1) {
                    if (output != null) output.write(readResult, 0, currentReadLen);
                    if (_jta != null) {
			byte[] readBytes = new byte[currentReadLen];
			System.arraycopy(readResult,0,readBytes,0,
					 currentReadLen);
			String readStr = new String(readBytes);
			_jta.append(readStr);
			DEditDebug.println(this,3,">>"+readStr+"<<");
		    }
		    currentReadLen = input.read(readResult, 0, 1024);
                }

                if (_sgOutput != null) _sgOutput.close();
                if (_sgInput != null) _sgInput.close();
            }
            catch (Exception ex) {
                if (_sgOutput != null) {
                    PrintWriter erWriter = new PrintWriter(_sgOutput);
                    erWriter.println("ServletExecuteScript.StreamGobbler exception :" + ex.toString());
                    erWriter.flush();
                    erWriter.close();
                }
            }
        }
    }

    /**
      * Makes things that aren't editable appear gray.
      * @param rndr the DEditTreeCellRenderer that may be modified and then 
      * will be returned
      * @param tree same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param denode same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param selected same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param expanded same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param leaf same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param row same as TreeCellRenderer.getTreeCellRendererComponent()
      * @param hasFocus same as TreeCellRenderer.getTreeCellRendererComponent()
      * @return same idea as TreeCellRenderer.getTreeCellRendererComponent()
      */
    public Component getTreeCellRendererComponent(DEditTreeCellRenderer rndr,
						  JTree tree, 
						  DEditNode denode, 
						  boolean selected, 
						  boolean expanded,
						  boolean leaf,
						  int row,
						  boolean hasFocus) {
	if (!isEditable(denode,true)) {
	    rndr.setForeground(Color.GRAY);
	}
	String name = denode.getNode().getNodeName();
	return rndr;
    }

    /**
      * Although unspecified in the interface (can't have static methods),
      * this method must be found in every plugin to construct it through
      * the command line argument "plugin=pluginName".
      * @param args the arguments used to build the returned DEditPlugin
      * @return a DEditPlugin that can be loaded
      */
    public static DEditPlugin buildPlugin(String[] args) {
	//String dirString = System.getProperty("user.dir");
	String dirString = "";
	return new ScriptPlugin(dirString,args);
    }

    /**
      * An alternative way of running the Plugin.
      * @param args arguments for DEditor and this Plugin
      */
    public static void main(String[] args) {
	String dirString = System.getProperty("user.dir");
	DEditPlugin pl = ScriptPlugin.buildPlugin(args);
	DEditor editor = DEditor.makeDEditor(args,pl,null);
	editor.start();
    }
    
}
