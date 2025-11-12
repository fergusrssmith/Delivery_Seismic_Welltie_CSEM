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

/**
 * Title:        BHP Viewer <br><br>
 * Description:  To comsume the stram. We need this class to consume the stdout and
 *               stderr of a process if they are not otherwise consumed.
 *               This prevents deadlock.
 *               User can provide an output stream in the constructor, so that
 *               the content of the input stream will be written to the output stream.
 *               Otherwise, the content is thrown away. <br><br>
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class StreamGobbler extends Thread {
    private InputStream _sgInput;
    private String _sgType;
    private OutputStream _sgOutput;
    private String _initMessage;

    public StreamGobbler(InputStream ins) {
        this(ins, "Out", null, "");
    }

    public StreamGobbler(InputStream ins, String type, String initMsg) {
        this(ins, type, null, initMsg);
    }

    public StreamGobbler(InputStream ins, String type, OutputStream out, String initMsg) {
        _sgInput = ins;
        _sgType = type;
        _sgOutput = out;
        _initMessage = "";
        if (initMsg != null && initMsg.length() > 0)
            _initMessage = initMsg + "\n\n";
    }

    public void run() {
        try {
/*
// for some reason, this does not work
            PrintWriter writer = null;
            if (_sgOutput != null)
                writer = new PrintWriter(_sgOutput);
            BufferedReader reader = new BufferedReader(new InputStreamReader(_sgInput));
            String line = reader.readLine();
            int index = 10;
            while(line != null || index < 10) {
                if (writer != null && line != null) writer.println(line);
                line = reader.readLine();
            }
*/

            BufferedOutputStream output = null;
            if (_sgOutput != null) {
                output = new BufferedOutputStream(_sgOutput);
                output.write(_initMessage.getBytes());
            }
            BufferedInputStream input = new BufferedInputStream(_sgInput);
            byte[] readResult = new byte[1024];
            int currentReadLen = 0;
            currentReadLen = input.read(readResult, 0, 1024);
            while(currentReadLen > -1) {
                if (output != null) output.write(readResult, 0, currentReadLen);
                currentReadLen = input.read(readResult, 0, 1024);
            }

            if (_sgOutput != null) _sgOutput.close();
            if (_sgInput != null) _sgInput.close();
        }
        catch (Exception ex) {
            if (_sgOutput != null) {
                PrintWriter erWriter = new PrintWriter(_sgOutput);
                erWriter.println("ServletExecuteScript.StreamGobbler exception :" + ex.toString());
                erWriter.flush();
                erWriter.close();
            }
        }
    }
}
