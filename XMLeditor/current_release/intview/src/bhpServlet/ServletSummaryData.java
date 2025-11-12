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
 * Description:  This servlet is responsible for running the bhpio command to
 *               summarize the data.
 *               It has only one parameter in the request named clientMessage,
 *               which is the whole bhpio command that it needs to run.
 *               The servlet uses Runtime.exec() to start the real bhpio process.
 *               The output of the process is written back to the client as
 *               plain text. <br><br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ServletSummaryData extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        doGet(request, response);
	}

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        response.setContentType("text/plain");
        PrintWriter owriter = response.getWriter();
        try {
            //BufferedOutputStream out = new BufferedOutputStream(response.getOutputStream());

            String command = request.getParameter("clientMessage");
            Process process = Runtime.getRuntime().exec(command);

            StreamGobbler errGobbler = new StreamGobbler(process.getErrorStream());
            errGobbler.start();
            InputStream summaryIps = process.getInputStream();
            BufferedReader reader = new BufferedReader(new InputStreamReader(summaryIps));
            String summaryLine = reader.readLine();
            while(summaryLine != null) {
                owriter.println(summaryLine);
                summaryLine = reader.readLine();
            }
            summaryIps.close();

            int exitValue = process.waitFor();
            if (exitValue != 0) {
                owriter.println("ServletSummayData process exit " + exitValue);
            }
            owriter.close();
/*
            process.waitFor();
            InputStream ins = process.getInputStream();
            int ilen = ins.available();
            byte[] buf = new byte[ilen];
            BufferedInputStream input = new BufferedInputStream(ins);
            input.read(buf, 0, ilen-1);
            input.close();

            response.setContentLength(ilen);
            out.write(buf, 0, ilen);
            out.close();
*/
        }
        catch (Exception ex) {
            owriter.println("ServletSummaryDataError");
            owriter.println("   " + ex.toString());
            owriter.close();
        }
	}
}