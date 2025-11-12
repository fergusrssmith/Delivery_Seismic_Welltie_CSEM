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

import java.io.File;
import java.util.Vector;
import java.util.Hashtable;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.border.EmptyBorder;
import javax.swing.*;


/**
  * Dialog for editing information in a leaf node (value, attributes, etc.).
  */
public class DEditLeafDialog {

    private static final int VALUE = 0;
    private static final int ATTR = 1;
    private static final int ALL = 2;
    private JDialog dialog;
    private DEditAttributeInfo[] attInfo; // information about attributes
    private String value; // value for leaf node
    private JComboBox valueBox; // text field that holds value
    private JTextField tf;
    private boolean set; //true if dialog was closed by choosing 
    //"Set" not "Quit"
    private DEditDataType dt;
    private boolean valueIsConstraint = false;
    private String[] constrMessages;
    private JLabel[] attrMessages;
    private JFrame frame;
    private String defValConstraint;
    private JTextArea valConstraintText;
    private Hashtable ht;
    private File suDir;
    private boolean supress;
    private DEditTable attTable;
    private boolean isApplet = false;

    /**
      * Constructor for DEditLeafDialog which sets up the current directory
      * and the special display properties.
      * @param curDir the File that will be the starting directory
      * @param configPairs key - value pairs to set up display properties.  
      * See README for detailed description.
      */
    public DEditLeafDialog(File curDir, String[][] configPairs,
			   boolean isApplet) {
	ht = new Hashtable();
	this.isApplet = isApplet;
	reset();
	suDir = curDir;
	if (configPairs == null) {
	    return;
	}
	for (int i = 0; i < configPairs.length; i++) {
	    if (configPairs[i] != null && configPairs[i].length == 2) {
		ht.put(configPairs[i][0],configPairs[i][1]);
	    }
	}
    }
    
    /**
      * Resets all fields for this class that need to be reset before
      * a new dialog can be created.
      */
    private void reset() {
	dialog = null;
	attInfo = null; // information about attributes
	value = null; // value for leaf node
	valueBox = null; // text field that holds value
	tf = null;
	set = false; //true if dialog was closed by choosing 
	dt = null;
	valueIsConstraint = false;
	constrMessages = null;
	attrMessages = null;
	frame = null;
	defValConstraint = null;
	valConstraintText = null;
	supress = false;
	attTable = null;
    }

    /**
      * Makes a JDialog with given parent JFrame, title, and initial value.
      * @param frame the parent JFrame for this DEditLeafDialog
      * @param title the String to go in the title of this dialog
      * @param dt the DEditDataType for the node this dialog will edit
      * @param initialValue the initial value (as a string) of the value 
      * in the leaf being edited
      * @return the new value if "Accept" is chosen, null otherwise
      */
    public String showValueDialog(JFrame frame, String title, 
				  DEditDataType dt,String initialValue) {
	reset();
	this.dt = dt;
	if (this.dt == null) {
	    //this.dt = new DEditDataType(DEditDataType.EPSILON,
	    //			DEditDataType.EMPTY,"",null);
	    return null;
	}
	this.value = initialValue;
	this.frame = frame;
	JPanel choosePanel = setUpChoosePanel(VALUE);
	JPanel valPanel = setUpValPanel(VALUE);
	JPanel jp = setUpMainPanel();
	jp.add(choosePanel,BorderLayout.SOUTH);
	jp.add(valPanel,BorderLayout.CENTER);
	dialog = new JDialog(frame,title,true);
	dialog.getContentPane().add(jp);
	dialog.pack();
	dialog.setLocationRelativeTo(frame);
	boolean setVal = showDialog();
	if (setVal) {
	    return value;
	} else {
	    return null;
	}
    }
    
    /**
      * Makes a JDialog with given parent JFrame, title, and information about
      * the attributes.
      * @param frame the parent JFrame for this DEditLeafDialog
      * @param title the String to go in the title of this dialog
      * @param atts a DEditAttributeInfo array used to store the information 
      * about attributes
      * @return the attribute array result (actually the same as atts) if 
      * "Accept" is chosen, null otherwise
      */
    public DEditAttributeInfo[] showAttributeDialog(JFrame frame, 
						    String title,
						    DEditAttributeInfo[] atts)
    {
	reset();
	this.attInfo = atts;
	this.frame = frame;
	JPanel choosePanel = setUpChoosePanel(ATTR);
	JPanel attPanel = setUpAttPanel();
	JPanel jp = setUpMainPanel();
	jp.add(choosePanel,BorderLayout.SOUTH);
	jp.add(attPanel,BorderLayout.NORTH);
	dialog = new JDialog(frame,title,true);
	dialog.getContentPane().add(jp);
	dialog.pack();
	dialog.setLocationRelativeTo(frame);
	boolean setVal = showDialog();
	if (setVal) {
	    return atts;
	} else {
	    return null;
	}
    }

    
    /**
      * Makes a JDialog with given parent JFrame, title, initial value,
      * and information about the attributes.
      * @param frame the parent JFrame for this DEditLeafDialog 
      * @param title the String to go in the title of this dialog
      * @param dt the DEditDataType for the node this dialog will edit
      * @param initialValue the initial value (as a string) of the value 
      * in the leaf being edited
      * @param atts a DEditAttributeInfo array used to store the information 
      * about attributes
      * @return the new value for the edited node if "Accept" is chosen, 
      * null otherwise.  atts has also changed to new values which should of 
      * course also not be set.
      */
    public String showAttValDialog(JFrame frame, String title,
				   DEditDataType dt, String initialValue,
				   DEditAttributeInfo[] atts) {
	reset();
	this.attInfo = atts;
	this.dt = dt;
	if (this.dt == null) {
	    //this.dt = new DEditDataType(DEditDataType.EPSILON,
	    //			DEditDataType.EMPTY,"",null);
	    return null;
	}
	this.value = initialValue;
	this.frame = frame;
	JPanel choosePanel = setUpChoosePanel(ALL);
	JPanel attPanel = setUpAttPanel();
	JPanel valPanel = setUpValPanel(ALL);
	JPanel jp = setUpMainPanel();
	jp.add(choosePanel,BorderLayout.SOUTH);
	jp.add(attPanel,BorderLayout.NORTH);
	jp.add(valPanel,BorderLayout.CENTER);
	dialog = new JDialog(frame,title,true);
	dialog.getContentPane().add(jp);
	dialog.pack();
	dialog.setLocationRelativeTo(frame);
	updateValueBox();
	boolean setVal = showDialog();
	if (setVal) {
	    return value;
	} else {
	    return null;
	}
    }

    /**
      * Creates and returns a JPanel with border to be used as the main 
      * panel in a DEditLeafDialog.
      * @return a JPanel to be used as the main panel in a DEditLeafDialog
      */
    private JPanel setUpMainPanel() {
	JPanel ret = new JPanel(new BorderLayout());
	ret.setBorder(new EmptyBorder(5,5,5,5));
	return ret;
    }

    /**
      * Creates, sets up, and returns a JPanel to display and edit the 
      * attributes of a node being edited by this DEditLeafDialog.
      * @return a JPanel which will be used to edit attributes
      */
    private JPanel setUpAttPanel() {
	JPanel attrPanel = new JPanel(new BorderLayout());
	attrPanel.add(new JLabel("Attributes"),BorderLayout.NORTH);
	for (int i = 0; i < attInfo.length; i++) {
	    DEditAttributeInfo inf = attInfo[i];
	    if (inf.exists == false) {
		inf.value = inf.datatype.getDefaultValue();
		inf.exists = true;
	    }
	}
	attTable = new DEditTable(this,attInfo);
	//attTable.setPreferredSize(new Dimension(400,100));
	JScrollPane jsp = new JScrollPane(attTable);
	jsp.setPreferredSize(new Dimension(400,100));
	attrPanel.add(jsp,BorderLayout.CENTER);
	return attrPanel;
    }

    /**
      * Creates, sets up, and returns a JPanel which will be used to
      * view and edit the value of the DEditNode being edited by this 
      * DEditLeafDialog.
      * @param whoToSet specifies which components to if the accept-type
      * action is launced (in this case via hitting enter in the JTextField)
      * @return a JPanel which will be used to edit a node's value
      */
    private JPanel setUpValPanel(final int whoToSet) {
	if (dt.getDataType() == DEditDataType.URI) {
	    return setUpFileChoosePanel(whoToSet);
	}
	tf = new JTextField(value);
	valueBox = new JComboBox();
	int type = dt.getType();
	JPanel textPanel;
	valueIsConstraint = (type == DEditDataType.CONSTRAINT);
	if (valueIsConstraint) {
	    textPanel = new JPanel(new BorderLayout());
	    defValConstraint = dt.getConstraintMessage();
	    valConstraintText = new JTextArea(defValConstraint,2,20);
	    valConstraintText.setEditable(false);
	    //valConstraintText.setLineWrap(true);
	    valConstraintText.setBackground(textPanel.getBackground());
	    textPanel.add(new JLabel("Value"),BorderLayout.NORTH);
	    textPanel.add(valConstraintText,BorderLayout.CENTER);
	    textPanel.add(tf,BorderLayout.SOUTH);
	    tf.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent ae) {
			boolean dispose = false;
			switch(whoToSet) {
			case VALUE:dispose = setValue();break;
			case ATTR:dispose = setAttributes();break;
			case ALL:dispose = setAll();break;
			default:dispose = false;
			}
			if (dispose) {
			    dialog.dispose();
			}
		    }
		});
	} else {
	    valueBox.setEditable(false);
	    for (int i = 0; i < dt.numValues(); i++) {
		valueBox.addItem(dt.getValueAt(i));
	    }
	    valueBox.setSelectedItem(value);
	    textPanel = new JPanel(new BorderLayout());
	    textPanel.add(new JLabel("Value"),BorderLayout.NORTH);
	    textPanel.add(valueBox,BorderLayout.SOUTH);
	}
	return textPanel;
    }

    /**
      * Creates and returns a JPanel used to choose a file and hand edit the 
      * path if need be
      * @param whoToSet says what to set (attributes, values, or both) when 
      * accept is chosen (in this case by hitting enter in the JTextField)
      * @return a JPanel used to choose a filename
      */
    private JPanel setUpFileChoosePanel(final int whoToSet) {
	tf = new JTextField(value,30);
	valueBox = new JComboBox();
	valueBox.setEditable(true);
	valueBox.addItem(value);
	JPanel textPanel = new JPanel(new BorderLayout());
	defValConstraint = dt.getConstraintMessage();
	valConstraintText = new JTextArea(defValConstraint,2,20);
	valConstraintText.setEditable(false);
	//valConstraintText.setLineWrap(true);
	valConstraintText.setBackground(textPanel.getBackground());
	textPanel.add(new JLabel("Value"),BorderLayout.NORTH);
	textPanel.add(valConstraintText,BorderLayout.CENTER);
	tf.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    boolean dispose = false;
		    switch(whoToSet) {
		    case VALUE:dispose = setValue();break;
		    case ATTR:dispose = setAttributes();break;
		    case ALL:dispose = setAll();break;
		    default:dispose = false;
		    }
		    if (dispose) {
			    dialog.dispose();
		    }
		}
	    });
	valueIsConstraint = true;
	JPanel browsePanel = new JPanel(new BorderLayout());
	browsePanel.add(tf,BorderLayout.WEST);
	JButton browseButton = new JButton("Browse...");
	browseButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    chooseFile();
		}
	    });
	browsePanel.add(browseButton,BorderLayout.EAST);
	textPanel.add(browsePanel,BorderLayout.SOUTH);
	return textPanel;
    }

    /**
      * Creates and returns a JPanel containing buttons for the ultimate user 
      * choices: accept and cancel
      * @param whoToSet if accept is chosen, this determines what to set 
      * (attributes, value, or both)
      * @return a JPanel used to dismiss the DEditLeafDialog either by 
      * accepting or canceling the changes made
      */
    private JPanel setUpChoosePanel(final int whoToSet) {
	JPanel buttonPanel = new JPanel(new GridLayout(1,2));
	JButton setButton = new JButton("Accept");
	setButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    DEditDebug.println(this,6,"accept");
		    boolean dispose = false;
		    switch(whoToSet) {
		    case VALUE:dispose = setValue();break;
		    case ATTR:dispose = setAttributes();break;
		    case ALL:dispose = setAll();break;
		    default:dispose = false;
		    }
		    if (dispose) {
			dialog.dispose();
		    }
		}
	    });
	buttonPanel.add(setButton);
	JButton cancelButton = new JButton("Cancel");
	cancelButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ae) {
		    cancel();
		}
	    });
	buttonPanel.add(cancelButton);
	return buttonPanel;
    }

    /**
      * Helper which launches a generic choose file dialog and sets the 
      * field tf (a JTextField) to have the filename as its text
      */
    private void chooseFile() {
	if (isApplet) {
	    String msg = "Specify remote URL";
	    RemoteFileChooser rfc 
		= new RemoteFileChooser(frame,"Choose File","browseFileURL");
	    int reFchooser = rfc.showDialog();
	    String filename = null;
	    if (reFchooser == JFileChooser.APPROVE_OPTION) {
		String pathString = rfc.getPath();
		String name = rfc.getName();
		filename = pathString + "/" + name;
		tf.setText(filename);
	    }
	    return;
	}
	JFileChooser chooser = new JFileChooser(suDir);
	chooser.setDialogTitle("Open File");
	chooser.setAcceptAllFileFilterUsed(true);
	int returnVal = chooser.showOpenDialog(dialog);
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	    File f = chooser.getSelectedFile();
	    suDir = f.getParentFile();
	    tf.setText(f.getAbsolutePath());
	}
    }
    
    /**
      * Updates the display of the value box.  This may mean changing it's 
      * information message, making it ineditable, or many other things.
      * Typically called as a result of changes to the attributes.
      */
    void updateValueBox() {
	if (valueIsConstraint) {
	    if (!tf.isEditable()) {
		tf.setText(value);
	    }
	    tf.setEditable(true);
	} else {
	    // what here???
	}
	boolean constraintTextSet = false;
	for (int i = 0; i < attInfo.length; i++) {
	    String val = attInfo[i].value;
	    if (!val.equals("")) {
		for (int j = 0; j < 2; j++) {
		    String key;
		    if (j == 0) {
			key = attInfo[i].name + "=";
		    } else {
			key = attInfo[i].name + "=" + val;
			
		    }
		    String restriction = (String)ht.get(key);
		    if (restriction != null) {
			int length = restriction.length();
			if (restriction.startsWith("ne") 
			    && (length == 2 || restriction.charAt(2) == '+')) {
			    if (valueIsConstraint) {
				value = tf.getText();
				tf.setEditable(false);
			    } else {
				// remove ???
				valueBox.setEditable(false);
			    }
			    if (length > 2) {
				restriction = restriction.substring(3);
			    }
			}
			if (restriction.startsWith("rcm:") 
			    && valueIsConstraint) {
			   valConstraintText.setText(restriction.substring(4));
			   constraintTextSet = true;
			} else if (restriction.startsWith("pcm:") 
				   && valueIsConstraint) {
			    valConstraintText.setText(restriction.substring(4)
						      + defValConstraint);
			    constraintTextSet = true;
			} else if (restriction.startsWith("acm:")
				   && valueIsConstraint) {
			    valConstraintText.setText(defValConstraint 
						   + restriction.substring(4));
			    constraintTextSet = true;
			} else if (restriction.startsWith("rvm:")
				   && valueIsConstraint) {
			    tf.setText(restriction.substring(4));
			}
		    }
		}
	    }
	}
	
	if (valConstraintText != null && !constraintTextSet 
	    && valueIsConstraint) {
	    valConstraintText.setText("\nConstant " + defValConstraint);
	}
    }

    /**
      * Sets the internal value to be used for the node's value with the 
      * value currently shown in the JTextField or JComboBox (as appropriate).
      * @return true if setting it to the new value is valid, false otherwise
      */
    private boolean setValue() {
	set = true;
	if (valueIsConstraint) {
	    if (tf.isEditable()) {
		value = tf.getText();
	    }
	} else {
	    value = (String)valueBox.getSelectedItem();
	}
	if (!dt.isValid(value)) {
	    JOptionPane.showMessageDialog(dialog,"Value for node is invalid.");
	    return false;
	}
	return true;
    }
    
    /**
      * Sets the internal value to be used for the values of the attributes of 
      * the node being edited with the current values in the table.
      * @return true if setting it to the new value is valid, false otherwise
      */
    private boolean setAttributes() {
	attTable.updateAttributeInfo();
	set = true;
	for (int i = 0; i < attInfo.length; i++) {
	    if (!attInfo[i].datatype.isValid(attInfo[i].value)) {
		JOptionPane.showMessageDialog(dialog,
					      "Value for attribute "
					      + attInfo[i].name
					      + " is invalid.");
		return false;
	    }
	}
	return true;
    }

    /**
      * Sets both attributes and values.
      * @return true if setting to these new attributes and values is valid,
      * false otherwise
      */
    private boolean setAll() {
	boolean check1 = setValue();
	boolean check2 = setAttributes();
	return (check1 && check2);
    }

    /**
      * Launched when the cancel button is hit.  Shouldn't really have to 
      * set set == false.
      */
    private void cancel() {
	set = false;
	dialog.dispose();
    }
    
    /**
      * Shows the dialog and returns the value of set which in turn will 
      * determine whether or not edited values should be set in the xml tree
      * or discarded.
      * @return true if the accept action was chosen
      */
    private boolean showDialog() {
	dialog.show();
	return set;
    }

}



