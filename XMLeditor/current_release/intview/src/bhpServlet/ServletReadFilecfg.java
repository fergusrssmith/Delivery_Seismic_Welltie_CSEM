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
 * Description:  This servlet is responsible for reading a configuration file.
 *               The only parameter in the request, fileName, specifies the name
 *               of the configuration file.
 *               The content of the file is written back to the client in the response.
 *               If threr is any problem, an error string is written in the response.
 *               <br> <br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ServletReadFilecfg extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        doGet(request, response);
	}

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        try {
            response.setContentType("text/html");
            //BufferedOutputStream out = new BufferedOutputStream(response.getOutputStream());
            PrintWriter out = response.getWriter();

            String fname = request.getParameter("FileName");
            File file = new File(fname);
            if (file.exists()==false || file.isDirectory()) {
                throw (new Exception(fname + " does not exist or is not a file"));
            }
            long len = file.length();
            int ilen = (int) len;
            response.setContentLength(ilen);
            char[] buf = new char[ilen];
            BufferedReader breader = new BufferedReader(new FileReader(file));
            breader.read(buf, 0, ilen);
            breader.close();
            out.print(buf);
            out.close();
        }
        catch (Exception ex) {
            response.getWriter().print("ServletReadFilecfgException " + ex.toString());
        }
	}
}