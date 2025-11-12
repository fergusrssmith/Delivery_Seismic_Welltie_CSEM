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

import java.io.*;
import java.awt.*;
import java.net.*;
import java.awt.event.*;
import java.util.StringTokenizer;
import java.util.Enumeration;
import javax.swing.*;
import javax.swing.event.*;

// this is just a slightly modified version of Synthia's original.  - WTR
/**
 * Title:        BHP Viewer <br><br>
 * Description:  This class provides general behavior and gui representation
 *               for choosing a remote file on the server. It uses servlet
 *               servletReadDirectory to read a directory and displays all its
 *               content in a list. When using this class, a subclass of this one
 *               needs to be created implementing the abstrace method execute.
 *               This method is where the selected file will be handled.
 *               This abstrace dialog can be used for both file openning and
 *               file saving.
 *               <br><br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class RemoteFileChooser extends JDialog {
    private static final int _GENERAL_FILE = 2;
    private static final int _DIRECTORY = 3;

    private JList _list;
    private DefaultListModel _listModel;
    private JTextField _pathField;
    private JTextField _nameField;
    private String _property;
    private String _homeDirectory;
    private int _returnValue;

    private ImageIcon _fileIcon = null;
    private ImageIcon _dirIcon = null;
    private ImageIcon _homeIcon = null;
    private ImageIcon _upIcon = null;
    private String _acceptButtonText = "  Open ";

    public RemoteFileChooser(Frame f, String title) {
        this(f, title, "xmleditFilecfgPath");
    }

    public RemoteFileChooser(Frame f, String title, String property,
			     String acceptButtonText) {
        super(f, title, true);
        _returnValue = JFileChooser.CANCEL_OPTION;
        _property = property;
        _homeDirectory = "";
	_acceptButtonText = acceptButtonText;
	//System.out.println("Accept set to " + _acceptButtonText);
        try {
	    ClassLoader cl = RemoteFileChooser.class.getClassLoader();
	    _fileIcon = new ImageIcon(cl.getResource("icons/GeneralFile.gif"));
	    _dirIcon = new ImageIcon(cl.getResource("icons/Directory.gif"));
	    _homeIcon = new ImageIcon(cl.getResource("icons/Home.gif"));
	    _upIcon = new ImageIcon(cl.getResource("icons/Up.gif"));
        }
        catch (Exception ex) {
            System.out.println("RemoteFileChooser.constructor failed to load images");
            System.out.println("    " + ex.toString());
        }
        setupGUI();
    }

    public RemoteFileChooser(Frame f, String title, String property) {
	this(f,title,property,"  Open ");
    }

    public String getPath() { return _pathField.getText(); }
    public String getName() { return _nameField.getText(); }

    public void updateList() {
        DEditor theViewer = null;
        try {
            theViewer = (DEditor) (SwingUtilities.windowForComponent(this));
        }
        catch (Exception ex) {
        }
        if (theViewer == null) {
            System.out.println("RemoteFileChooser.updateList Error : Cannot find DEditor instance");
            return;
        }

        String urlString = theViewer.getProperty("servletReadDirectoryURL");
        //String pathString = _pathField.getText() + "/" + _nameField.getText();
        String pathString = _pathField.getText();
        _listModel.removeAllElements();

        try {
            URL url = new URL(urlString + "?FileName=" + pathString);
            URLConnection conn = url.openConnection();
            BufferedReader in = new BufferedReader(new InputStreamReader(conn.getInputStream()));
            String line = in.readLine();
            if (line != null && line.length() > 0) {
                StringTokenizer stk = new StringTokenizer(line, ";");
                String token;
                while(stk.hasMoreTokens()) {
                    token = stk.nextToken().trim();
                    if (token.startsWith("3")) {
                        _listModel.addElement(new ListObject(_DIRECTORY, token.substring(1)));
                    }
                    else {
                        _listModel.addElement(new ListObject(_GENERAL_FILE, token.substring(1)));
                    }
                }
            }
            _pathField.setText(pathString);
            _nameField.setText("");
        }
        catch (Exception e) {
            System.out.println("RemoteFileChooser update Exception " 
			       + e.toString());
            e.printStackTrace();
        }
    }

    public boolean isExistingFile(String name) {
        Enumeration enum = _listModel.elements();
        ListObject obj = null;
        while(enum.hasMoreElements()) {
            obj = (ListObject) enum.nextElement();
            if (obj.getType() == _GENERAL_FILE && obj.toString().equals(name))
                return true;
        }
        return false;
    }

    public int showDialog() {
        this.updateList();
        this.show();
	if (_returnValue == JFileChooser.APPROVE_OPTION) {
	    DEditor theViewer = null;
	    try {
		theViewer=(DEditor)(SwingUtilities.windowForComponent(this));
	    } catch (Exception ex) {
	    }
	    if (theViewer == null) {
		System.err.println("RemoteFileChooser.updateList Error : Cannot find DEditor instance");
	    } else {
		theViewer.setProperty(_property,_pathField.getText());
	    }
	}
        return _returnValue;
    }

    public void hideDialog() {
	if (_returnValue == JFileChooser.APPROVE_OPTION) {
	    DEditor theViewer = null;
	    try {
		theViewer=(DEditor)(SwingUtilities.windowForComponent(this));
	    } catch (Exception ex) {
	    }
	    if (theViewer == null) {
		System.err.println("RemoteFileChooser.updateList Error : Cannot find DEditor instance");
	    } else {
		theViewer.setProperty(_property,_pathField.getText());
	    }
	}
        //this.setVisible(false);
        super.hide();
    }

    private void setupGUI() {
        this.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                _returnValue = JFileChooser.CANCEL_OPTION;
            }
        });

        JPanel extraButtonPanel = new JPanel(new FlowLayout());
        JButton upButton = new JButton(_upIcon);
        upButton.setPreferredSize(new Dimension(25, 25));
        JButton hmButton = new JButton(_homeIcon);
        hmButton.setPreferredSize(new Dimension(25, 25));
        upButton.addActionListener(new UpListener());
        hmButton.addActionListener(new HomeListener());
        extraButtonPanel.add(upButton);
        extraButtonPanel.add(hmButton);

        JPanel p1 = new JPanel(new BorderLayout());
        p1.add(new JLabel("  Look in:    "), BorderLayout.WEST);
        _pathField = new JTextField("");
        p1.add(_pathField, BorderLayout.CENTER);
        p1.add(extraButtonPanel, BorderLayout.EAST);

        //p1.add(new JLabel("File Name"));
        //_nameField = new JTextField("");
        //p1.add(_nameField);

        _pathField.addKeyListener(new MyPathKeyListener());

        DEditor theViewer = null;
        try {
            theViewer = (DEditor) (SwingUtilities.windowForComponent(this));
        }
        catch (Exception ex) {
        }
        if (theViewer != null) {
            _homeDirectory = theViewer.getProperty(_property);
	    if (_homeDirectory == null) {
		_homeDirectory = theViewer.getProperty("xmleditFilecfgPath");
	    }
            _pathField.setText(_homeDirectory);
        }

        JPanel p2 = new JPanel(new BorderLayout());
        _listModel = new DefaultListModel();
        _list = new JList(_listModel);
        _list.setCellRenderer(new MyCellRenderer());
        _list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        _list.addMouseListener(new MyListMouseAdapter());
        _list.addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (_list.getSelectedValue() != null)
                    _nameField.setText(_list.getSelectedValue().toString());
            }
        });
        p2.add(new JScrollPane(_list), BorderLayout.CENTER);

        //JPanel p3 = new JPanel(new FlowLayout());
        JButton b1 = new JButton(_acceptButtonText);
	//System.out.println("using text " + _acceptButtonText);
        b1.addActionListener(new ExecuteListener());
        JButton b3 = new JButton("Cancel");
        b3.addActionListener(new CancelListener());
        //p3.add(b1);
        //p3.add(b3);

        JPanel p4 = new JPanel(new BorderLayout(10, 10));
        p4.add(new JLabel("  File name:      "), BorderLayout.WEST);
        _nameField = new JTextField("");
        p4.add(_nameField, BorderLayout.CENTER);
        p4.add(b1, BorderLayout.EAST);

        JPanel p5 = new JPanel(new BorderLayout(10, 10));
        p5.add(new JLabel("  Files of type:  "), BorderLayout.WEST);
        String[] types = {"All Files(*.*)"};
        p5.add(new JComboBox(types));
        p5.add(b3, BorderLayout.EAST);

        JPanel panel = new JPanel(new BorderLayout(10, 10));
        panel.add(p2, BorderLayout.CENTER);
        panel.add(p1, BorderLayout.NORTH);

        Box panel2 = new Box(BoxLayout.Y_AXIS);
        panel2.add(panel2.createVerticalStrut(10));
        panel2.add(p4);
        panel2.add(panel2.createVerticalStrut(10));
        panel2.add(p5);
        panel2.add(panel2.createVerticalStrut(10));

        this.getContentPane().add(panel, BorderLayout.CENTER);
        this.getContentPane().add(panel2, BorderLayout.SOUTH);

        this.pack();
        this.setSize(500, 350);
    }

    private class MyListMouseAdapter extends MouseAdapter {
        public void mouseClicked(MouseEvent e) {
            if (SwingUtilities.isLeftMouseButton(e) == false) return;
            if (e.getClickCount() < 2) return;

            JList source = (JList) e.getSource();
            int index = source.locationToIndex(e.getPoint());
            source.setSelectedIndex(index);

            Object selectedObj = source.getSelectedValue();
            if (selectedObj instanceof ListObject) {
                ListObject listObj = (ListObject) selectedObj;
                int type = listObj.getType();
                if ( type == _DIRECTORY) {
                    String oldPath = _pathField.getText();
                    if (!oldPath.endsWith("/")) oldPath = oldPath + "/";
                    _pathField.setText(oldPath + listObj.toString());
                    updateList();
                }
                else {
                    _returnValue = JFileChooser.APPROVE_OPTION;
                    hideDialog();
                }
            }
        }
    }

    private class ExecuteListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            Object selectedObj = _list.getSelectedValue();
            if (selectedObj instanceof ListObject) {
                ListObject listObj = (ListObject) selectedObj;
                int type = listObj.getType();
                if ( type == _DIRECTORY) {
                    String oldPath = _pathField.getText();
                    if (!oldPath.endsWith("/")) oldPath = oldPath + "/";
                    _pathField.setText(oldPath + listObj.toString());
                    updateList();
                }
                else {
                    _returnValue = JFileChooser.APPROVE_OPTION;
                    hideDialog();
                }
            } else if (!_nameField.getText().equals("")) {
		_returnValue = JFileChooser.APPROVE_OPTION;
		hideDialog();
	    }
        }
    }

    private class CancelListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            _returnValue = JFileChooser.CANCEL_OPTION;
            hideDialog();
        }
    }

    private class UpListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            String curDir = _pathField.getText();
            String upDir = curDir;
            if (curDir != null && curDir.lastIndexOf("/") >= 0)
                upDir = curDir.substring(0, curDir.lastIndexOf("/"));
            if (upDir.length() == 0) upDir = "/";
            _pathField.setText(upDir);
            updateList();
        }
    }

    private class HomeListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            _pathField.setText(_homeDirectory);
            updateList();
        }
    }

    private class ListObject {
        private int _type;
        private String _name;

        public ListObject(int type, String name) {
            _type = type;
            _name = name;
        }

        public int getType() { return _type; }

        public String toString() {
            return _name;
        }

        public ImageIcon getSymbolIcon() {
            if (_type == _DIRECTORY) return _dirIcon;
            else return _fileIcon;
        }
    }

    private class MyCellRenderer extends JLabel implements ListCellRenderer {
        public Component getListCellRendererComponent(JList list, Object value,
                            int index, boolean isSelected, boolean cellHasFocus) {
            if (value == null) return this;
            String s = value.toString();
            setOpaque(true);
            setText(s);
            if (value instanceof ListObject)
                setIcon(((ListObject)value).getSymbolIcon());
            if (isSelected) setBackground(list.getSelectionBackground());
            else setBackground(list.getBackground());
            return this;
        }
    }

    private class MyPathKeyListener extends KeyAdapter {
        public void keyPressed(KeyEvent e) {
            if (e.getKeyCode() == KeyEvent.VK_ENTER)
                updateList();
        }
    }
}
