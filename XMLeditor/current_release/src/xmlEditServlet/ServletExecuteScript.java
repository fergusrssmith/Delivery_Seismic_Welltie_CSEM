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

import java.io.*;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * Title:        BHP Viewer <br><br>
 * Description:  This servlet is responsible for executing script.
 *               The request is firstly checked for its scriptType parameter. <br>
 *               If it is "existing", the request then has two more parameters.
 *               scriptName is the name of a script file that user wants to execute.
 *               clientMessage is the parameter for executing the script.
 *               If it is not "existing", the request then has two more parameters.
 *               outputDir is the directory where the new script file will be generated.
 *               clientMessage is the whole script command. The servlet will generat
 *               a new file in the specified directory.
 *               In both cases, the script is executed with Runtime.exec().
 *               The stdout and stderr of the process are written to log file
 *               $script$.out.log and $script$.err.log respectively.
 *               Information about running the script will be sent back to the client.
 *               <br> <br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ServletExecuteScript extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
	    doGet(request, response);
    }

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        Process process = null;
 		PrintWriter owriter = response.getWriter();
        String type = request.getParameter("scriptType");
        String logOutFname = "scriptOutput.log";
        String logErrFname = "scriptError.log";
        String initMessage = "Execute Script: ";
        try {
            if (type.equals("existing")) {
                String script = request.getParameter("scriptName");
                String parameter = request.getParameter("clientMessage");
                File scriptFile = new File(script);

                logOutFname = script + ".out.log";
                logErrFname = script + ".err.log";
                initMessage = initMessage + script + " " + parameter;

                if (scriptFile.exists() && scriptFile.isFile()) {
                    String[] tokens = new String[2];
                    tokens[0] = script;
                    tokens[1] = parameter;
                    process = Runtime.getRuntime().exec(tokens);
                }
            }
            else {  // new
                String outDir = request.getParameter("outputDir");
                String command = request.getParameter("clientMessage");
                if (outDir.endsWith("/") == false) outDir = outDir + "/";
                int fileIndex = 0;
                File tmp;
                while(true) {
                    tmp = new File(outDir + "userscript" + fileIndex);
                    if (tmp.exists() == false) break;
                    fileIndex++;
                }
                String commandFile = outDir + "userscript" + fileIndex;

                FileWriter writer = new FileWriter(new File(commandFile));
                writer.write(command);
                writer.close();

                logOutFname = commandFile + ".out.log";
                logErrFname = commandFile + ".err.log";
                initMessage = initMessage + command;

                process = Runtime.getRuntime().exec("sh " + commandFile);
            }

            if (process != null) {
                File of = new File(logOutFname);
                File ef = new File(logErrFname);
                FileOutputStream logFile = new FileOutputStream(of);
                FileOutputStream elogFile = new FileOutputStream(ef);
                StreamGobbler errorGobbler = new StreamGobbler(process.getErrorStream(),
                                            "Error", elogFile, initMessage);
                StreamGobbler outputGobbler = new StreamGobbler(process.getInputStream(),
                                            "Output", logFile, initMessage);
                errorGobbler.start();
                outputGobbler.start();

                int exitValue = process.waitFor();
                if (exitValue != 0) {
                    owriter.println("ServletExecuteScript process exit abnormal : " + exitValue);
                }
                else {
                    owriter.println("ServletExecuteScript process exit normal");
                }
                owriter.close();
            }

        }
        catch (Exception ex) {
            owriter.println("ServletExecuteScript Exception :" + ex.toString());
            owriter.close();
        }
	}
}
