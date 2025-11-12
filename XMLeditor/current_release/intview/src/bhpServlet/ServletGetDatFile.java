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
import java.util.StringTokenizer;

import javax.servlet.*;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * Title:        BHP Viewer <br><br>
 * Description:  This servlet is responsible for reading seismic data.
 *               It handles four different types of request: Begin, ReadTrace,
 *               Finish, and SaveEvnet. <br>
 *               In the case of Begin, the request has two more parameters,
 *               parameter1 is the bhpread command, parameter2 is a suggested name
 *               of the output file. The servlet will use Runtime.exec() to run
 *               the bhpread process and the output of that process is redirected
 *               the the output file. If this process finishs successfully, the
 *               servlet will continue to run the bhpread summary using the same
 *               technique. If both succeed, the actual name of the output file
 *               and the cotent of the read summary will be sent back to the client
 *               as plain text. Otherwise, an error message is sent. The output file
 *               generated here is a temporary file, it will be deleted eventually
 *               by the application. <br>
 *               In the case of ReadTrace, the request has three more parameters.
 *               parameter1 is the name of file for reading.
 *               parameter2 is the offset indicates where to start the reading.
 *               parameter3 indicates how much to read.
 *               The servlet will open the file as random access and read the
 *               specified segment of the file. Whatever is read will be sent
 *               to the client in the binary format. <br>
 *               In the case of Finish, the servlet will delete the temporary file
 *               it generated with bhpread. The request has parameter1 which is the
 *               name of the file that needs to be deleted. <br>
 *               SaveEvent will update the temporary file and use bhpwrite to update
 *               the real data. It is used to save the event after it is edit by the
 *               user. In this case, the request has five more parameters.
 *               parameter1 is the bhpwrite command needed to update the disk file.
 *               parameter2 is the name of the temporary file.
 *               parameter3 represents new values. it's a space separated string.
 *               parameter4 is the size of each trace in the temporary file.
 *               parameter5 is the offset in the trace where writing should begin.
 *               The servlet will open the temporary file for writing, update each
 *               trace with the new value. Then bhpwrite command will be run with
 *               Runtime.exec() to update the real disk file. The stderr of the process
 *               will be sent back to the client. If any problem happens, an error
 *               message will be sent back. <br> <br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ServletGetDatFile extends HttpServlet {
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
	    doGet(request, response);
    }

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
 		String reqType = request.getParameter("type");
        if (reqType.equals("Begin")) {
            beginRead(request, response);
        }
        else if (reqType.equals("ReadTrace")) {
            readTrace(request, response);
        }
        else if (reqType.equals("Finish")) {
            finishRead(request, response);
        }
        else if (reqType.equals("SaveEvent")) {
            saveEvent(request, response);
        }
	}

    private void saveEvent(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        response.setContentType("text/plain");
        PrintWriter owriter = response.getWriter();

        String command = request.getParameter("parameter1");
        String ofname = request.getParameter("parameter2");
        String valueString = request.getParameter("parameter3");
        String traceSizeString = request.getParameter("parameter4");
        String traceOffsetString = request.getParameter("parameter5");
        int traceSize = 0;
        int traceOffset = 0;
        StringTokenizer stk = new StringTokenizer(valueString);
        int valueLength = stk.countTokens();
        double[] values = new double[valueLength];
        try {
            int index = 0;
            while(stk.hasMoreTokens()) {
                values[index] = Double.parseDouble(stk.nextToken());
                index++;
            }
            traceSize = Integer.parseInt(traceSizeString);
            traceOffset = Integer.parseInt(traceOffsetString);
        }
        catch (Exception ex1) {
            owriter.println("ServletGetDatFile.saveEvent exception1 : " + ex1.toString());
            owriter.close();
            return;
        }

        RandomAccessFile wfile = null;
        try {
            wfile = new RandomAccessFile(ofname, "rw");
        }
        catch (Exception ex2) {
            owriter.println("ServletGetDatFile.saveData exception2 " + ex2.toString());
            owriter.close();
            return;
        }
        if (wfile == null) {
            owriter.println("ServletGetDatFile.saveData null file : " + ofname);
            owriter.close();
            return;
        }

        byte[] bvalue = new byte[4];
        try {
            for (int i=0; i<values.length; i++) {
                wfile.seek(i*traceSize + traceOffset);
                convertToBytes(values[i], bvalue);
                wfile.write(bvalue);
            }
        }
        catch (Exception ex3) {
            owriter.println("ServletGetDatFile.saveData exception3 :" + ex3.toString());
            owriter.close();
            return;
        }

        String[] cmdArray = new String[3];
        cmdArray[0] = "sh";
        cmdArray[1] = "-c";
        cmdArray[2] = command;
        try {
            Process wprocess = Runtime.getRuntime().exec(cmdArray);
            StreamGobbler errorGobbler = new StreamGobbler(wprocess.getInputStream());
            errorGobbler.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(wprocess.getErrorStream()));
            String output = reader.readLine();
            while(output != null) {
                owriter.println("    WW " + output);
                output = reader.readLine();
            }
            int exitValue = wprocess.waitFor();
            owriter.println("ServletGetDatFile.saveEvent exit " + exitValue);
            owriter.close();
        }
        catch (Exception wex) {
            owriter.println("ServletGetDatFile.saveEvent wexception : " + wex.toString());
            owriter.close();
            return;
        }
    }

    private void beginRead(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        response.setContentType("text/plain");
        PrintWriter owriter = response.getWriter();
        //ObjectOutputStream outputStream = new ObjectOutputStream(response.getOutputStream());

        String command = request.getParameter("parameter1");
        String ofname = request.getParameter("parameter2");
/*
        int suIndex = 0;
        File tmp;
        while(true) {
            tmp = new File(ofname + suIndex);
            if (tmp.exists() == false) break;
            suIndex++;
        }
*/
        try {
            //ofname = new String(ofname + suIndex);
            String tmpDir = ofname.substring(0, ofname.lastIndexOf("/"));
            ofname = File.createTempFile("BHP_", ".sutmp", new File(tmpDir)).getAbsolutePath();
            String fcommand = command + " > " + ofname;
            String[] cmdArray = new String[3];
            cmdArray[0] = "sh";
            cmdArray[1] = "-c";
            cmdArray[2] = fcommand;
            Process readProcess = Runtime.getRuntime().exec(cmdArray);

            StringBuffer errorOutput = new StringBuffer();
            BufferedReader reader = new BufferedReader(new InputStreamReader(readProcess.getErrorStream()));
            String line = reader.readLine();
            while(line != null) {
                errorOutput.append(line);
                line = reader.readLine();
            }

            int exitValue = readProcess.waitFor();
            //cmdFile.delete();

            if (exitValue != 0) {     // exit abnormal
                owriter.println("BHPMODELVIEWERERRORERROR");
                owriter.println("read abnormal exit");
                owriter.println(fcommand);
                owriter.println(errorOutput.toString());
                owriter.close();
                File file = new File(ofname);
                file.delete();
                return;
            }
        }
        catch (Exception ex) {
            owriter.println("BHPMODELVIEWERERRORERROR");
            owriter.println("read exception " + ex.toString());
            owriter.close();
            return;
        }

        try {
            int insertIndex = command.indexOf("keys=");
            String queryParameter = command.substring(0, insertIndex) +
                                    " request=summary " +
                                    command.substring(insertIndex);
            Process summaryProcess = Runtime.getRuntime().exec(queryParameter);

            StreamGobbler errorGobbler = new StreamGobbler(summaryProcess.getErrorStream());
            errorGobbler.start();

            InputStream summaryIps = summaryProcess.getInputStream();
            BufferedReader reader = new BufferedReader(new InputStreamReader(summaryIps));
            String summaryContent = null;
            String summaryLine = reader.readLine();
            while(summaryLine != null) {
                summaryContent = summaryLine;
                summaryLine = reader.readLine();
            }
            summaryIps.close();

            int exitValue = summaryProcess.waitFor();
            if (exitValue != 0) {  // abnormal exit
                owriter.println("BHPMODELVIEWERERRORERROR");
                owriter.println("summary abnormal exit");
                owriter.close();
                File file = new File(ofname);
                file.delete();
                return;
            }
            else {
                // write to the response
                //owriter.println(ofname);
                owriter.println("\n\n");
                owriter.println("BHPMODELVIEWERNAME " + ofname);
                owriter.println(summaryContent);
                owriter.close();
            }
        }
        catch (Exception ex) {
            owriter.println("BHPMODELVIEWERERRORERROR");
            owriter.println("summary exception " + ex.toString());
            owriter.close();
            return;
        }
    }

    private void readTrace(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        response.setContentType("binary");
        String fname = request.getParameter("parameter1");
        int offset = Integer.parseInt(request.getParameter("parameter2"));
        int length = Integer.parseInt(request.getParameter("parameter3"));

        RandomAccessFile file = new RandomAccessFile(fname, "r");
        byte[] data = new byte[length];
        file.seek(offset);
        file.readFully(data);

        response.setContentLength(length);
        BufferedOutputStream out = new BufferedOutputStream(response.getOutputStream());
        out.write(data, 0, length);

        out.close();
        file.close();
    }

    private void finishRead(HttpServletRequest request, HttpServletResponse response) throws IOException, ServletException {
        String fname = request.getParameter("parameter1");
        File file = new File(fname);
        PrintWriter out = response.getWriter();
        if (file.exists() && file.isFile()) {
            file.delete();
            out.println("File deleted " + fname);
            out.close();
            return;
        }
        out.println("Error deleting file " + fname);
        out.close();
    }

    private void convertToBytes(double value, byte[] bv) {
        int bits = Float.floatToIntBits((float)value);
        bv[0] = (byte)((bits & 0xff000000) >> 24);
        bv[1] = (byte)((bits & 0xff0000) >> 16);
        bv[2] = (byte)((bits & 0xff00) >> 8);
        bv[3] = (byte)(bits & 0xff);
    }

}