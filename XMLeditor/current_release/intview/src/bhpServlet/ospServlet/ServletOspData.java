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

public class ServletOspData extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
	    doGet(request, response);
    }

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
              
    String commend = request.getParameter("ospkey");
    //System.out.println("commend = " + commend);

    if(commend.equals("loadObject"))
    	loadObject(request, response);
    else if(commend.equals("create"))
    	createSelector(request, response);
    else if(commend.equals("clearObject"))
    	clearObject(request, response);
    	
	}//doGet
	
	private void clearObject(HttpServletRequest request, HttpServletResponse response){
		String commend = request.getParameter("type");
		HttpSession httpss = request.getSession();
		OspObjectSelector oos = (OspObjectSelector)httpss.getAttribute(commend);
		oos.clearObject();
		
	}//createSelector

	
	private void createSelector(HttpServletRequest request, HttpServletResponse response){
		String commend = request.getParameter("type");
		HttpSession httpss = request.getSession();
			
		if(commend.equals("Seismic")){
			httpss.setAttribute(commend, new OspSeismicSelector());
		}else if(commend.equals("Horizon")){
			httpss.setAttribute(commend, new OspHorizonSelector());
		}
	}//createSelector

	
	private void loadObject(HttpServletRequest request, HttpServletResponse response){
		
		String commend = request.getParameter("type");
		HttpSession httpss = request.getSession();
		
		OspObjectSelector oos = (OspObjectSelector)httpss.getAttribute(commend);
		
		OspSessionController ospSc = (OspSessionController)httpss.getAttribute("ospSesCont");

		oos.loadObject(ospSc);
		String[] nms = oos.getAllItemNames();
		String returnInfo = "";
		
		if(nms.length >0){
			for(int i=0; i<nms.length -1; i++)
				nms[0] = nms[0] + "=" + nms[i+1];
			returnInfo = nms[0];
		}
		
		try{
 		  PrintWriter out = response.getWriter();
      out.println(returnInfo);
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}
		
	}//loadObject
}