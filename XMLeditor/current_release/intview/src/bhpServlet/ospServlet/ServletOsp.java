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
public class ServletOsp extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
	    doGet(request, response);
    }

  /**
  *	KeysIn: ospkey, host, uid, pw,
  * 
  * 
  * KeysOut: loginTrue, loginFalse, 
  */
	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {


		//response.setContentType("text/html"); //tmp
 		//PrintWriter out = response.getWriter();
    String commend = request.getParameter("ospkey");
    
    if(commend.equals("isInLogin"))
    	isInLogin(request, response);
    else if(commend.equals("login"))
    	handleLogin(request, response);
 		else if(commend.equals("getSessionList")) 
 			getSessionList(request, response);
 		else if(commend.equals("connect")) 
 			handleConnect(request, response);
		else if(commend.equals("isConnected"))
			isConnected(request, response);
		else if(commend.equals("disconnect"))
			disconnect(request, response);
		else if(commend.equals("logout"))
    	handleLogout(request, response);
  
	}//do get
	
	/** handleLogout */
	private void handleLogout(HttpServletRequest request, HttpServletResponse response){
			HttpSession httpss = request.getSession();						
			try{
			Object att = httpss.getAttribute("ospSesCont");
			OspSessionController ospSc = (OspSessionController)att;
			//ospSc.getSession().user().checkInLicense("BHP_BHPViewer");
			ospSc.disconnect();
 		  PrintWriter out = response.getWriter();
      out.println("True");
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}
			httpss.invalidate();
	}

	/** check to see if it is in login state */
	private void isInLogin(HttpServletRequest request, HttpServletResponse response){
			HttpSession httpss = request.getSession();
			String returnInfo = "False";
			
			if(httpss.getAttribute("ospSesCont") != null)	
				returnInfo = "True";
						
			try{
 		  PrintWriter out = response.getWriter();
      out.println(returnInfo);
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}
	}

	/** Get the session list and send to the applet */
	private void disconnect(HttpServletRequest request, HttpServletResponse response){
			HttpSession httpss = request.getSession();

			String returnInfo = "False";
							
			Object att = httpss.getAttribute("ospSesCont");
			OspSessionController ospSc = (OspSessionController)att;
			ospSc.disconnect();	
	}
	
	/** Get the session list and send to the applet */
	private void isConnected(HttpServletRequest request, HttpServletResponse response){
			HttpSession httpss = request.getSession();

			String returnInfo = "False";
							
			Object att = httpss.getAttribute("ospSesCont");
			OspSessionController ospSc = (OspSessionController)att;
			if(ospSc != null && ospSc.isConnected())
				returnInfo = "True";
						
			try{
 		  PrintWriter out = response.getWriter();
      out.println(returnInfo);
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}

	}
	
	/** Get the session list and send to the applet */
	private void handleConnect(HttpServletRequest request, HttpServletResponse response){
			HttpSession httpss = request.getSession();
			String ospSession = request.getParameter("ospSession");
							
			Object att = httpss.getAttribute("ospSesCont");
			OspSessionController ospSc = (OspSessionController)att;
			ospSc.connect(ospSession);
			
			String returnInfo = "True";
			if(!ospSc.isConnected())
				returnInfo = "False";

			try{
 		  PrintWriter out = response.getWriter();
      out.println(returnInfo);
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}
			
	}
	
	/** Get the session list and send to the applet */
	private void getSessionList(HttpServletRequest request, HttpServletResponse response){
			HttpSession httpss = request.getSession();
			String returnInfo = "True";
			
			if(httpss.isNew()){
				returnInfo = "False";
			}	
			Object att = httpss.getAttribute("ospSesCont");

			OspSessionController ospSc = (OspSessionController)att;
			String[] ssList = ospSc.getSessionList();
			if(ssList == null) System.out.println("getSessionList3");

			if(ssList == null) returnInfo = "False";
			
			for(int i=0; i<ssList.length; i++){
					returnInfo = returnInfo + "=" + ssList[i];
			}
			
			try{
 		  PrintWriter out = response.getWriter();
      out.println(returnInfo);
     	out.close();
			}catch (Exception ex) { ex.printStackTrace();	}

	}
	
	/** Get the host, uid, password, login the user account handleLogin */
	private void handleLogin(HttpServletRequest request, HttpServletResponse response){
		
			String returnInfo = "True";
			String host = request.getParameter("host");
			String uid  = request.getParameter("uid");
			String pw = request.getParameter("pw");
			
			//System.getProperty("java.class.path");
			//System.out.println(System.getProperty("java.class.path"));
			
			OspSessionController ospSc = OspSessionController.instance();
			
			org.openspirit.EpUser.UserProfile upf = null;
			try{
				upf = com.openspirit.osp.util.UserProfile.getUserProfile(host, uid, pw);
			  //upf.checkOutLicense("BHP_BHPViewer");
			}catch (Exception ex) {
				ex.printStackTrace();
				returnInfo = "False";
			}
			
			org.openspirit.EpSession.SessionManager smg = 
			org.openspirit.EpSession.SessionManagerHelper.narrow(upf.sessionMgr());		
			ospSc.setSessionManager(smg);
						
			HttpSession httpss = request.getSession();
			httpss.setAttribute("ospSesCont", ospSc);

			try{
 		  PrintWriter out = response.getWriter();
      out.println(returnInfo);
     	out.close();
			}catch (Exception ex) { System.out.println(ex); }

	}//handleLogin
	
	public void init(ServletConfig config){
		  
		  //set environment variables
			String ospenv = config.getInitParameter("OSPenv");
						
			String[] envs = split(ospenv, "#");//ospenv.split("#");
			String[] env = null;
			for(int i=0; i<envs.length; i++){		
				env = split(envs[i], "=");//envs[i].split("=");
				System.setProperty(env[0], env[1]);
			}
					  
		  try{
			super.init(config);
			}catch (Exception ex) { System.out.println("init false"); }
			
	}//init
	
	/** Use this one because BHP have not use jdk1.4 yet */
	private String[] split(String s, String d){
		java.util.StringTokenizer st = new java.util.StringTokenizer(s,d);
		int cont = st.countTokens();
		String[] rs = new String[cont];
		for(int i=0; i<cont; i++)
			rs[i] = st.nextToken();
		return rs;
	}
}