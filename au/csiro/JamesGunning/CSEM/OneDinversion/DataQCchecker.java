/**********************************************************************
DELIVERYCSEM - a program for Bayesian CSEM inversion
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
package au.csiro.JamesGunning.CSEM.OneDinversion;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

import cern.colt.list.*;

import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.lib.graphutil.*;

public class DataQCchecker 
// implements  WindowListener
{
 
  // WindowListener
  /*
	public void windowClosed(WindowEvent e) {}
	public void windowClosing(WindowEvent e) {}
	public void windowOpened(WindowEvent e) {}    
	public void windowIconified(WindowEvent e) {}    
	public void windowDeiconified(WindowEvent e) {}    
	public void windowActivated(WindowEvent e) {}
	public void windowDeactivated(WindowEvent e) {}
  */

  /*
	Multiple plots of data types, with possible outliers hilighted
  */
  public void plotData(ArrayList<DataArray> dArrs)
  {
	int i,j,Nplots;
	DataArray dArr;
	CanvasPlotter[] CParr;
	CanvasPlotter CP;
	CanvasPlotter.TwoDPlotGroup pg;
	CanvasPlotter.TwoDPlot pN,pO;   // Normal data, O = outliers

	Nplots = dArrs.size();

	CParr = new CanvasPlotter[Nplots];
	ArrayList<CanvasPlotter.TwoDPlotGroup> Pgroups;

	for (i=0; i<Nplots; i++)
	  {
		dArr = dArrs.get(i);
		CParr[i] = new CanvasPlotter(false);
		CP =  CParr[i];
		Pgroups = new ArrayList<CanvasPlotter.TwoDPlotGroup>(1);
		pg = new CanvasPlotter.TwoDPlotGroup(CP);
		Pgroups.add(pg);



		// "normal data"
		pN = new CanvasPlotter.TwoDPlot();
		pN.c = Color.blue;
		pN.x = dArr.X.elements(); 
		pN.y = dArr.Y.elements();
		pN.e = dArr.E.elements();
		pN.symbol = 'r';
		pN.drawline = false;

		// "outliers"
		pO = new CanvasPlotter.TwoDPlot();
		pO.c = Color.red;
		pO.x = dArr.Xo.elements(); 
		pO.y = dArr.Yo.elements();
		pO.e = dArr.Eo.elements();
		pO.symbol = 'r';
		pO.drawline = false;

		if (dArr.logView) 
		  {
			double[] tmp;
			pg.GG.logPlotY = true;
			tmp = (double[])(pN.y.clone());
			for(j=0;j<tmp.length;j++) tmp[j]=Math.abs(pN.y[j]);
			pN.y = tmp;

			tmp = (double[])(pO.y.clone());
			for(j=0;j<tmp.length;j++) tmp[j]=Math.abs(pO.y[j]);
			pO.y = tmp;
		  }
		pg.plots.add(pN);
		pg.GG.computeRanges(pN.x,pN.y);
		pg.plots.add(pO);
		pg.GG.computeRanges(pO.x,pO.y);

		CP.windowTitle = String.format("Median filter outlier detection f=%.2f data-code=%d",dArr.f,dArr.code);
		// title
		pg.GG.title = CP.windowTitle;
		CP.run(Pgroups); 
		// CP.frame.addWindowListener(this); // to listen for all graphs being shut off		
	  }
  }

  public static class DataArray
  {
	float f;
	int code;
	boolean logView;
	DoubleArrayList X,Y,E;
	DoubleArrayList Yo,Xo,Eo;

	public DataArray()
	{
	  X = new DoubleArrayList();
	  Y = new DoubleArrayList();
	  E = new DoubleArrayList();
	  Xo = new DoubleArrayList();
	  Yo = new DoubleArrayList();
	  Eo = new DoubleArrayList();
	}

	void trim()
	{
	  X.trimToSize();
	  Y.trimToSize();
	  E.trimToSize();
	  Xo.trimToSize();
	  Yo.trimToSize();
	  Eo.trimToSize();
	}
  }
}