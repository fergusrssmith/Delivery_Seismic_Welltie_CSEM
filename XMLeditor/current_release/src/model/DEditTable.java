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

import javax.swing.plaf.basic.BasicComboBoxEditor;
import javax.swing.table.*;
import javax.swing.JTable;
import javax.swing.JComboBox;
import java.awt.Component;
import java.util.EventObject;
import javax.swing.event.*;

/**
  * Class for displaying the attributes of a node in the JTree as a JTable 
  * with a few special features.  This class could be contained inside 
  * DEditLeafDialog but it's big enough to warrant it's own space.
  */
public class DEditTable extends JTable {
    
    final private DEditAttributeInfo[] attInfo; //attribute names, values, etc.
    final private DEditLeafDialog owner; // the DEditLeafDialog that created 
    // this DEditTable
    final private String[] columnNames = new String[]{"name","value"}; 
    final private Object[][] rowData;

    /**
      * Creates a new DEditTable.
      * @param ownerDialog the DEditLeafDialog this DEditTable will be
      * displayed inside
      * @param atts the info about the attributes this DEditTable will 
      * represent
      */
    public DEditTable(DEditLeafDialog ownerDialog,
		      DEditAttributeInfo[] atts) {
	super();
	this.attInfo = atts;
	this.owner = ownerDialog;
	rowData = new Object[attInfo.length][2];
	setModel(new AbstractTableModel() {
		public String getColumnName(int col) { 
		    return columnNames[col].toString(); 
		}
		public int getRowCount() { return rowData.length; }
		public int getColumnCount() { return columnNames.length; }
		public Object getValueAt(int row, int col) { 
		    return rowData[row][col]; 
		}
		public boolean isCellEditable(int row, int col)
		{
		    if (col == 0) {
			return false;
		    } else {
			return true;
		    }
		}
		public void setValueAt(Object value, int row, int col) {
		    rowData[row][col] = value;
		    fireTableCellUpdated(row, col);
		}
	    });
       
	TableColumn values = getColumnModel().getColumn(1);
	DEditTableCellEditor cellEditor = new DEditTableCellEditor(attInfo);
	cellEditor.addCellEditorListenerToAll(new CellEditorListener() {
		public void editingCanceled(ChangeEvent e) {
		    DEditDebug.println(this,6,"canceled");
		}
		public void editingStopped(ChangeEvent e) {
		    DEditDebug.println(this,6,"stopped");
		    updateAttributeInfo();
		    owner.updateValueBox();
		} 
	    });
	values.setCellEditor(cellEditor);
	for (int i = 0; i < attInfo.length; i++) {
	    setValueAt(attInfo[i].name,i,0);
	    setValueAt(attInfo[i].value,i,1);
	}
    }

    /**
      * Called to update the information stored about attributes to reflect
      * what is shown on the screen
      */
    public void updateAttributeInfo() {
	for (int i = 0; i < attInfo.length; i++) {
	    attInfo[i].value = (String)getValueAt(i,1);
	    DEditDebug.println(this,6,attInfo[i].value);
	}
    }

    /**
      * Inner class to handling cell editing for the "value" column of the 
      * table.  Design is heavily influenced by swing examples at 
      * http://www2.gol.com/users/tame/swing/examples/SwingExamples.html.
      */
    private class DEditTableCellEditor implements TableCellEditor {
	
	private NoBoldTableCellEditor[] cellEditors;
	private int which = 0;

	public DEditTableCellEditor(DEditAttributeInfo[] atts) {
	    this.cellEditors = new NoBoldTableCellEditor[atts.length];
	    for (int i = 0; i < cellEditors.length; i++) {
		DEditAttributeInfo temp = atts[i];
		DEditDataType dt = temp.datatype;
		int type = dt.getType();
		JComboBox newBox = null;
		if (type == DEditDataType.CHOICE) {
		    newBox = new JComboBox();
		    for (int j = 0; j < dt.numValues(); j++) {
			newBox.addItem(dt.getValueAt(j));
		    }
		    newBox.setEditable(false);
		} else if (type == DEditDataType.CONSTRAINT) {
		    newBox = new JComboBox();
		    newBox.setEditable(true);
		} else if (type == DEditDataType.FIXED) {
		    newBox = new JComboBox();
		    newBox.setEditable(false);
		} else {
		    newBox = null;
		}
		newBox.setEditor(new BasicComboBoxEditor());
		cellEditors[i] = new NoBoldTableCellEditor(newBox);
	    }
	}
	
	public Component getTableCellEditorComponent(JTable table, 
						     Object value,
						     boolean isSelected,
						     int row, int column) {
	    which = row;
	    return cellEditors[which].getTableCellEditorComponent(table,value,
								  isSelected,
								  row,column);
	    
	}

	public Object getCellEditorValue() {
	    return cellEditors[which].getCellEditorValue();
	}
	public Component getComponent() {
	    return cellEditors[which].getComponent();
	}
	public boolean stopCellEditing() {
	    return cellEditors[which].stopCellEditing();
	}
	public void cancelCellEditing() {
	    cellEditors[which].cancelCellEditing();
	}
	public boolean isCellEditable(EventObject anEvent) {
	    return cellEditors[which].isCellEditable(anEvent);
	}
	public boolean shouldSelectCell(EventObject anEvent) {
	    return cellEditors[which].shouldSelectCell(anEvent);
	}
	public void addCellEditorListener(CellEditorListener l) {
	    cellEditors[which].addCellEditorListener(l);
	}
	public void addCellEditorListenerToAll(CellEditorListener l) {
	    for (int i = 0; i < cellEditors.length; i++) {
		cellEditors[i].addCellEditorListener(l);
	    }
	}
	public void removeCellEditorListener(CellEditorListener l) {
	    cellEditors[which].removeCellEditorListener(l);
	}
	public void removeCellEditorListenerFromAll(CellEditorListener l) {
	    for (int i = 0; i < cellEditors.length; i++) {
		cellEditors[i].removeCellEditorListener(l);
	    }
	}
	public void setClickCountToStart(int n) {
	    cellEditors[which].setClickCountToStart(n);
	}
	public int getClickCountToStart() {
	    return cellEditors[which].getClickCountToStart();
	}
    }
}
