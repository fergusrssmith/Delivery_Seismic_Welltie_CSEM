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


package au.csiro.JamesGunning.lib.graphutil;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.font.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JFrame;
import javax.swing.KeyStroke;
import java.awt.image.*;
import java.awt.print.PrinterJob;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.ErrorHandler;



/**
   Basic canvas extension with draw zoom.
   @author James.Gunning@csiro.au
   @version 1.0, 30/4/2008
*/

public class ZoomableCanvas extends Canvas implements MouseListener, MouseMotionListener,
													  ComponentListener

{
  // graphics stuff
  public static int DEFAULT_SIZE = 900;	

  public int canvasHeight = DEFAULT_SIZE;  // component height
  public int canvasWidth = canvasHeight;      // component width

 // I/O & error stuff
  ErrorHandler EH;

  public JFrame frame;

  // Stuff for Zoom managing
  private boolean mouseDown=false;
  private Rectangle currentRect;
  private Rectangle rectToDraw = null;
  private Rectangle previousRectDrawn = new Rectangle();
  public AffineTransform at = new AffineTransform();  // current snapshot
  private AffineTransform xAT = new AffineTransform(); // for preconcats...
  private BasicStroke stroke;
  private double currentScale=1.0;

  public ZoomableCanvas()
  {
	// handle errors: can change how these are handled if desired
	EH = new ErrorHandler();

	// Add mouse handlers
	addMouseMotionListener(this);
	addMouseListener(this);
	// etc
	addComponentListener(this);
	// setup Zoom transforms
	at.setToIdentity();
  }


  /** 
	  Example  test driver to show how to view extraction info from typical
	  MLgraphing-file.*.txt files produced by the WabeletExtractor.
	  The display could also be driven from some other class that has a
	  fully populated MLdata object, and then calling the display() method.
  */
  public static void main(String argv[]) 
  {
	 ZoomableCanvas ZC = new ZoomableCanvas();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("ZoomableCanvas",ZC.EH);

	// view everything. This can be called as soon as data is fully populated
	ZC.display("ZoomableCanvas");
  }

  public void display(String title)
  {
	display(title,null);
  }

  public void display(String title,String epsFilename)
  {
	JMenuBar menuBar; // FM 
	JMenuItem menuItem;

	frame = new JFrame(title);
	//	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);


	setSize(new Dimension (canvasWidth,canvasHeight));
	setBackground(Color.lightGray);
	setForeground(Color.black);
	frame.getContentPane().add(this);



	/*
	  // Create the menu bar.
	  menuBar = new JMenuBar();
	  menuItem = new JMenuItem("Print",KeyEvent.VK_P);
	  menuItem.setAccelerator(KeyStroke.getKeyStroke('p'));
	  menuItem.addActionListener(this);
	  menu.add(menuItem);
	*/


	// frame.setLocation(0,0) ;
	frame.pack();
	frame.setVisible((epsFilename == null) ?  true : false); 

  } 

  private void setGraphTitle(String title)
  {
	frame.setTitle(title);
  }

  public void paint(Graphics g) 
  {
	Graphics2D g2 = (Graphics2D) g;


	// save context Xform
	AffineTransform saveXform = g2.getTransform();
    // apply new one
	g2.setTransform(at);
	// then draw away....

	
	// here , the main graphs we want to see
	paintContent(g2);

	// paint Zoom & mouse 3 in old coords
	g2.setTransform(saveXform);

	// zoom box
	paintZoomBox(g2);
  }

  private void paintContent(Graphics2D g2)
  {
	double ws,hs,w,h;
	
	// Here, scale everything to canvasHeight, canvasWidth
	ws = 0.1*canvasHeight;
	hs = 0.1*canvasWidth;
	h = 0.8*canvasHeight;
	w = 0.8*canvasWidth;

	g2.draw(new Line2D.Double(ws,hs,ws,hs+h));
	g2.draw(new Line2D.Double(ws,hs+h,ws+w,hs+h));
	g2.draw(new Line2D.Double(ws+w,hs+h,ws+w,hs));
	g2.draw(new Line2D.Double(ws+w,hs,ws,hs));
  }

  public void mousePressed(MouseEvent e) 
  {
	int modifiers = e.getModifiers();
	if ((modifiers & InputEvent.BUTTON1_MASK) == InputEvent.BUTTON1_MASK)
	  {
		// mouse 1 zoom
		currentRect = new Rectangle(e.getX(), e.getY(), 0, 0);
		updateDrawableRect(getWidth(), getHeight());
	  }
  }
  
  public void mouseDragged(MouseEvent e) 
  {
	if ((e.getModifiers() & InputEvent.BUTTON1_MASK) == InputEvent.BUTTON1_MASK)
	  { 
		// mouse 1 zoom
		// update canvas scale
		updateSize(e);
		mouseDown=true;
	  }
  }
    
  public void mouseReleased(MouseEvent e) 
  {
	if ((e.getModifiers() & InputEvent.BUTTON1_MASK) == InputEvent.BUTTON1_MASK)
	  {
		// mouse 1 zoom
		updateSize(e);
		if (mouseDown)
		  {
			// end of drag. Do A ZoomBox
			// translation defined by left edge of box
			xAT.setToTranslation(-rectToDraw.x,-rectToDraw.y);
			at.preConcatenate(xAT);
			
			// and scale by box ratios
			double Sx = (double)getWidth()/Math.max((double)rectToDraw.width,1);
			double Sy = (double)getHeight()/Math.max((double)rectToDraw.height,1);
			
			xAT.setToScale(Sx,Sy);
			at.preConcatenate(xAT);
			currentScale *= Math.sqrt(Sx*Sx+Sy*Sy);
		  }
		else
		  {
			// restore original
			at.setToIdentity();
			currentScale=1.0;
		  }
		mouseDown=false;
		
		// to clear box...
		currentRect=null;
		repaint();
	  }
  }
  
  /* 
   * Update the size of the current rectangle
   * and call repaint.  Because currentRect
   * always has the same origin, translate it
   * if the width or height is negative.
   * 
   * For efficiency (though
   * that isn't an issue for this program),
   * specify the painting region using arguments
   * to the repaint() call.
   * 
   */
  private void updateSize(MouseEvent e)
  {
	if (currentRect!=null)
	  {
		currentRect.setSize(e.getX() - currentRect.x, e.getY() - currentRect.y);
		updateDrawableRect(getWidth(), getHeight());
		Rectangle totalRepaint = rectToDraw.union(previousRectDrawn);
		repaint(totalRepaint.x, totalRepaint.y,	totalRepaint.width, totalRepaint.height);
	  }
  }

  public void paintZoomBox(Graphics g) 
  {
	//If currentRect exists, paint a box on top.
	if (currentRect != null) 
	  {
		Graphics2D g2 = (Graphics2D) g;
		// Draw a rectangle on top of the image.
		Paint oldp = g2.getPaint();
		g2.setPaint(Color.black);
		g2.setXORMode(Color.white); //Color of line varies
		// depending on image colors
		g2.drawRect(rectToDraw.x, rectToDraw.y, 
				   rectToDraw.width - 1, rectToDraw.height - 1);
		g2.setPaint(oldp);		
		updateLabel(rectToDraw);
	  }
  }

  private void updateLabel(Rectangle rect) 
  {
	int width = rect.width;
	int height = rect.height;
	
	//Make the coordinates look OK if a dimension is 0.
	if (width == 0) {
	  width = 1;
	}
	if (height == 0) {
	  height = 1;
	}
  }

  private void updateDrawableRect(int compWidth, int compHeight) 
  {
	int x = currentRect.x;
	int y = currentRect.y;
	int width = currentRect.width;
	int height = currentRect.height;
    
	//Make the width and height positive, if necessary.
	if (width < 0) {
	  width = 0 - width;
	  x = x - width + 1; 
	  if (x < 0) {
		width += x; 
		x = 0;
	  }
	}
	if (height < 0) {
	  height = 0 - height;
	  y = y - height + 1; 
	  if (y < 0) {
		height += y; 
		y = 0;
	  }
	}
    
	//The rectangle shouldn't extend past the drawing area.
	if ((x + width) > compWidth) {
	  width = compWidth - x;
	}
	if ((y + height) > compHeight) {
	  height = compHeight - y;
	}
          
	//Update rectToDraw after saving old value.
	if (rectToDraw != null) {
	  previousRectDrawn.setBounds(
								  rectToDraw.x, rectToDraw.y, 
								  rectToDraw.width, rectToDraw.height);
	  rectToDraw.setBounds(x, y, width, height);
	} else {
	  rectToDraw = new Rectangle(x, y, width, height);
	}
  }




  // This method required by MouseListener.
  public void mouseMoved(MouseEvent e){}
  
  // These methods are required by MouseMotionListener.
  public void mouseClicked(MouseEvent e){}
  public void mouseExited(MouseEvent e){}
  public void mouseEntered(MouseEvent e){}

  // Component Listener
  public void componentResized(ComponentEvent e) 
  {
	Rectangle R = e.getComponent().getBounds();
	// resize to smallest dim.
	// canvasHeight = (int)Math.min(R.getWidth(),R.getHeight());
	canvasHeight = (int)R.getHeight();
	canvasWidth = (int)R.getWidth();

	repaint(); 
  }


  public void componentHidden(ComponentEvent e) {}
  public void componentMoved(ComponentEvent e) {}
  public void componentShown(ComponentEvent e) {}

} 
