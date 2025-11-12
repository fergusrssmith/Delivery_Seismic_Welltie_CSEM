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
package bhpServlet.ospServlet;

import java.io.*;
import javax.servlet.*;
import javax.servlet.http.*;
import com.interactive.jgeo.seismic.data.cgSeismicConstants;

import external.ospDataAdapter.*;

/**
 * Title:        ServletOspDataQuaryer <br><br>
 * Description:  This servlet is responsible for interacting with openspirit,
 *               and quary seismic and horizon data. <br>
 *               <br> <br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ServletOspDataSource extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
	    doGet(request, response);
    }

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
              
    String commend = request.getParameter("ospkey");
    
    if(commend.equals("getHeaders"))
    	getHeaders(request, response);
    else if(commend.equals("volumeReader"))	
    	createSeismicReader(request, response);
    else if(commend.equals("horizonReader"))
    	createHorizonReader(request, response);
	}//doGet
	
	/** createSeismicReader infomation */
	private void createHorizonReader(HttpServletRequest request, HttpServletResponse response){
		String session = request.getParameter("session");
		String dataName = request.getParameter("dataName");
		String range1 = request.getParameter("range1");
		String range2 = request.getParameter("range2");
		String property = request.getParameter("property");
			
		//String sourceName = "SOURCE" + "=" + session + "=" + "Horizon" + "=" + dataName;
		HttpSession httpss = request.getSession();				
		OspDataSource dataSource = (OspDataSource)httpss.getAttribute("SOURCE");
		OspHorizonGrid hg = (OspHorizonGrid)dataSource.getInfo();
		
		double[] iRange = hg.passInputRange(true, range1);
		double[] xRange = hg.passInputRange(false, range2);
		hg.setInlineRange(iRange[0], iRange[1]);
		hg.setCrosslineRange(xRange[0], xRange[1]);
				
		float[] data = hg.getSubGrid(property);

		try{
 		  DataOutputStream out = new DataOutputStream(response.getOutputStream());
      out.writeInt(hg.getNumInline());
      out.writeInt(hg.getNumCrossline());
      out.writeInt(getHorizonSampleUnit(hg));
      float[] rg = hg.getValueRange(property);
			out.writeFloat(rg[0]);
			out.writeFloat(rg[1]);
			out.writeDouble(iRange[0]);
			//out.writeDouble(iRange[1]);
			out.writeDouble(xRange[0]);
			//out.writeDouble(xRange[1]);
			
			for(int i=0; i<data.length; i++)
				out.writeFloat(data[i]);
     	out.close();
     	
			}catch (Exception ex) { ex.printStackTrace();	}
			
	}//createHorizonReader

	
	
	/** createSeismicReader infomation */
	private void createSeismicReader(HttpServletRequest request, HttpServletResponse response){
		String session = request.getParameter("session");
		String dataName = request.getParameter("dataName");
		String range1 = request.getParameter("range1");
		String range2 = request.getParameter("range2");
					
		//String sourceName = "SOURCE" + "=" + session + "=" + "Seismic" + "=" + dataName;
		HttpSession httpss = request.getSession();
		OspDataSource dataSource = (OspDataSource)httpss.getAttribute("SOURCE");
		OspSeismicVolume svl = (OspSeismicVolume)dataSource.getInfo();
		
		double[] iRange = svl.passInputRange(true, range1);
		double[] xRange = svl.passInputRange(false, range2);
		org.openspirit.EpiSeismic.BinGridSubVolumeAccessor sass = 
				svl.getSubVolumeAccessor(iRange, xRange);	
				
		float[] data = svl.getSubVolumeData(sass);
		try{
 		  DataOutputStream out = new DataOutputStream(response.getOutputStream());
      out.writeInt(sass.numU());
      out.writeInt(sass.numV());
			out.writeInt(sass.numW());
			out.writeInt(svl.getSampleUnit());
			double[] rs = getSampleRateAndStart(svl);
			out.writeDouble(rs[0]);
			out.writeDouble(rs[1]);
			out.writeDouble(iRange[0]);
			out.writeDouble(xRange[0]);
			out.writeDouble(rs[2]); //dt
			out.writeDouble(rs[3]); //delrt

			
			for(int i=0; i<data.length; i++)
				out.writeFloat(data[i]);
			
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}
			
	}//createSeismicReader

	/** getHorizonSampleUnit */
	private int getHorizonSampleUnit(OspHorizonGrid hg){
		org.openspirit.EpUnitOfMeasureService.Unit unit = hg.getCoordinateSys().coordUnit3();

			int un = cgSeismicConstants.CG_SU_TIME;
        try{
        if (unit.isIdentical(com.openspirit.osp.util.UnitOfMeasure.feetUS()) ||
                unit.isIdentical(com.openspirit.osp.util.UnitOfMeasure.feet())) {
                un = cgSeismicConstants.CG_SU_DEPTH_FEET;
            }
            else if (unit.isIdentical(com.openspirit.osp.util.UnitOfMeasure.metre())) {
                un = cgSeismicConstants.CG_SU_DEPTH_METERS;
            }
            else if (unit.isEquivalent(com.openspirit.osp.util.UnitOfMeasure.second())) {
                un = cgSeismicConstants.CG_SU_TIME;
            }
            else if (unit.isEquivalent(com.openspirit.osp.util.UnitOfMeasure.metre())) {
                un = cgSeismicConstants.CG_SU_DEPTH_METERS;
            }
            }catch (Exception e) { //e.printStackTrace(); }
                System.out.println("Cannot handle the unit: " + unit.name());
            }
			return un;
	}
	/** getSampleRateAndStart */
	private double[] getSampleRateAndStart(OspSeismicVolume svl){
		
			// Determine sample start in user SBG space
            double startSample = svl.getMinSampleNumber();

            // Get sample rate
            double sr = svl.getSampleRate();
						
            try{
            // Set sample rate and start value based on sample unit
            switch(svl.getSampleUnit()) {
    	    	case cgSeismicConstants.CG_SU_DEPTH_FEET: // feet
                    sr = com.openspirit.osp.util.UnitOfMeasure.convert(
                        sr, svl.getVerticalUnit(),
                        com.openspirit.osp.util.UnitOfMeasure.feetUS());
                    startSample = com.openspirit.osp.util.UnitOfMeasure.convert(
                        startSample, svl.getVerticalUnit(),
                        com.openspirit.osp.util.UnitOfMeasure.feetUS());
	    	    	break;
    	    	case cgSeismicConstants.CG_SU_DEPTH_METERS: // meters
                    sr = com.openspirit.osp.util.UnitOfMeasure.convert(
                        sr, svl.getVerticalUnit(),
                        com.openspirit.osp.util.UnitOfMeasure.metre());
                    startSample = com.openspirit.osp.util.UnitOfMeasure.convert(
                        startSample, svl.getVerticalUnit(),
                        com.openspirit.osp.util.UnitOfMeasure.metre());
	    	    	break;
                case cgSeismicConstants.CG_SU_TIME:
                    sr = com.openspirit.osp.util.UnitOfMeasure.convert(
                        sr, svl.getVerticalUnit(),
                        com.openspirit.osp.util.UnitOfMeasure.second());
                    startSample = com.openspirit.osp.util.UnitOfMeasure.convert(
                        startSample, svl.getVerticalUnit(),
                        com.openspirit.osp.util.UnitOfMeasure.second());
                    break;
	        	default:
                    sr = 0.0;
                    startSample = 0.0;
                    System.err.println("Don't know the unit for the volume.  It's not foot, meter,");
                    System.err.println("millisecond, or second.");
            } // switch
				}catch (Exception ex) { ex.printStackTrace();	}

				double scale = checkVerticalUnit(svl);
				double dt = sr*scale*1000;
		    double delrt = startSample*scale*1000;

        double[] rt = {sr, startSample, dt, delrt};
				return	rt;
	}
	 
	/** Check to see if the Vertical Unit is time */
	private  double checkVerticalUnit( OspSeismicVolume smv){
		double scale = 0;
		org.openspirit.EpUnitOfMeasureService.UnitPackage.Conversion convsion = null;

		// get the unit of measure finder using the OpenSpirit utility class.
     try {
         org.openspirit.EpUnitOfMeasureService.Unit fUnit = 
         			smv.getVerticalUnit();
         						
				 org.openspirit.EpUnitOfMeasureService.Unit tUnit = 
				 	com.openspirit.osp.util.UnitOfMeasure.millisecond();
				
				if(fUnit.isEquivalent(tUnit)){
						convsion = fUnit.conversionTo (tUnit);
						scale = (convsion.a + convsion.b * 1)/(convsion.c + convsion.d * 1);
				}else{
					  tUnit = com.openspirit.osp.util.UnitOfMeasure.feet();
						convsion = fUnit.conversionTo (tUnit);
						scale = (convsion.a + convsion.b * 1)/(convsion.c + convsion.d * 1);
						scale = scale/5d;
				} 
      }catch (Exception e) {System.err.println(e); }
		return scale;
	}

	/** get the headers infomation */
	private void getHeaders(HttpServletRequest request, HttpServletResponse response){
		String type = request.getParameter("type");
		String session = request.getParameter("session");
		String dataName = request.getParameter("dataName");
		

		HttpSession httpss = request.getSession();
		OspObjectSelector oos = (OspObjectSelector)httpss.getAttribute(type);
		OspSessionController ospSc = (OspSessionController)httpss.getAttribute("ospSesCont");

		oos.setSelectedItem(dataName);
		Object obj = oos.getSelectedItem();
		
		OspDataSource dataSource = null;
		if(type.equals("Seismic")){
			org.openspirit.EpiSeismic.PostStack3d seismic = (org.openspirit.EpiSeismic.PostStack3d) obj;          
      dataSource = new OspDataSource(session, new OspSeismicVolume(seismic));

		}else if(type.equals("Horizon")){
			org.openspirit.EpiInterpretation.Representation horizon = 
			(org.openspirit.EpiInterpretation.Representation) obj;
			dataSource =
            new OspDataSource(session,new OspHorizonGrid(horizon, ospSc.getSession()));
		}
		
		String[] header = dataSource.getHeaders();
		String returnInfo = header[0] + "=" + header[1] + "=" + header[2];
		
		try{
 		  PrintWriter out = response.getWriter();
      out.println(returnInfo);
      out.println(dataSource.getProperties());
      returnInfo = dataSource.getHeaderSetting(header[0])
      	 + "=" + dataSource.getHeaderSetting(header[1])
      	 + "=" + 	dataSource.getHeaderSetting(header[2]);
			out.println(returnInfo);
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}
			
		//cleanDataSource(httpss);
		//String sourceName = "SOURCE" + "=" + session + "=" + type + "=" + dataName;
		httpss.setAttribute("SOURCE", dataSource);

	}//loadObject
	
	/** Clean the data source from the web session
	private void cleanDataSource(HttpSession httpss){
			java.util.Enumeration e = httpss.getAttributeNames();
			String name = null;
			
			for (; e.hasMoreElements() ;) {
           name = (String)e.nextElement();
					 if(name.startsWith("SOURCE")){
					 		httpss.removeAttribute(name);
					 }	
       }//for
	 }//cleanDataSource*/
}