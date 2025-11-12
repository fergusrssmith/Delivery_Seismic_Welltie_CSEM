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
package bhpServlet;

import java.io.*;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * Title:        BHP Viewer <br><br>
 * Description:  This servlet is responsible for saving the current working
 *               session in a configuration file.
 *               It has two parameter in the request, fileName and content.
 *               "fileName" is the name of the would-be configuration file.
 *               "content" is an XML string representing the working session.
 *               The response of the servlet consists of the information
 *               whether the saving is successful. <br><br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ServletSaveFilecfg extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        doGet(request, response);
	}

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        PrintWriter out = response.getWriter();
        String fname = null;
        try {
/*
            fname = request.getParameter("FileName");
            File file = new File(fname);
            if (file.exists()) file.delete();
            BufferedInputStream binput = new BufferedInputStream(request.getInputStream());
            BufferedOutputStream boutput = new BufferedOutputStream(new FileOutputStream(fname));
            if (binput != null) {
                out.println("Okay input " + binput.available());
                out.flush();
            }
            byte[] buffer = new byte[binput.available()];
            binput.read(buffer);
            boutput.write(buffer);
            binput.close();
            boutput.close();
            out.println("ServletSaveFilecfg successfully save sesseion to " + fname);
            out.close();
*/

            fname = request.getParameter("fileName");
            String content = request.getParameter("content");

            File file = new File(fname);
            if (file.exists()) file.delete();

            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(fname)));
            writer.write(content);
            writer.close();

            out.println("ServletSaveFilecfg successfully save sesseion to " + fname);
            out.close();

        }
        catch (Exception ex) {
            out.println("ServletSaveFilecfg exception saving " + fname + " : " + ex.toString());
            ex.printStackTrace(out);
            out.println("FinishStackTrace");
            out.close();
        }
	}
}