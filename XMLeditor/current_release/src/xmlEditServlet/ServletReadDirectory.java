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
package xmlEditServlet;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * Title:        BHP Viewer <br><br>
 * Description:  This servlet is responsible for listing a directory.
 *               The only parameter in the request, FileName, specifies the name
 *               of the directory.
 *               The content of the directory listing is stored in a semicolon
 *               separated string and written back to the client in the response
 *               as plain text.
 *               If threr is any problem, an error string is written in the response.
 *               <br> <br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ServletReadDirectory extends HttpServlet {
	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
		response.setContentType("text/plain");
		PrintWriter out = response.getWriter();

 		String dname = request.getParameter("FileName");
		File directory = new File(dname);

 		if (directory.exists() && directory.isDirectory()) {
 			StringBuffer msg = new StringBuffer();
/*
            String[] files = directory.list();
            java.util.Arrays.sort(files);
            for (int i=0; i<files.length; i++) {
                msg.append(files[i] + ";");
            }
*/
 			File[] files = directory.listFiles();
            String[] fileNames = new String[files.length];
 			for (int i=0; i<files.length; i++) {
                if (files[i].isDirectory())
                    fileNames[i] = new String("3" + files[i].getName() + ";");
                else fileNames[i] = new String("2" + files[i].getName() + ";");
 			}
            java.util.Arrays.sort(fileNames);
            for (int j=0; j<fileNames.length; j++)
                msg.append(fileNames[j]);

 			//out.println(msg);
            out.println(msg);
			out.close();
		}
		else {
			out.println("Failed reading " + dname + ";");
			out.close();
		}
	}
}
