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
package script;

import java.io.*;
import java.util.Vector;
import model.DEditDebug;
import model.DEditor;

/**
  * Scatter parameter which represents a list of values found in a file.
  */
public class ListScatterParam extends AScatterParam {
    /**
      * the name of the file from which the values will be gathered
      */
    public String listfile;
    private DEditor editor;

    /**
      * Constructor which takes the name of the scatter parameter and the 
      * name of the file in which to find the list of values for the parameter.
      * @param name the name of the scatter parameter
      * @param listfile the name of the file which has the list of values for
      * this scatter parameter.
      * @param editor regrettable addition required so that in the case of 
      * an applet it can read the remote file
      */
    public ListScatterParam(String name, String listfile, DEditor editor) {
	this.name = name;
	this.listfile = listfile;
	this.editor = editor;
    }

    /**
      * Returns true if this scatter parameter should be resolved after the 
      * initial round of resolving.  For ListScatterParam it always returns 
      * false since the list values are predetermined and there is no way to 
      * compress them.
      * @return false
      */
    public boolean resolveLater() {
	return false;
    }

    /**
      * Resolves this ListScatterParam into a ResolvedScatterParam and returns
      * it.  In this case that means reading the values off the list from the
      * file and ignoring the maximum number of threads.
      * @param currThreads not used
      * @param maxThreads not used
      * @return a ResolvedScatterParam for this ListScatterParam
      */
    public ResolvedScatterParam resolve(int currThreads, int maxThreads) 
	throws MalformedScriptException {
	String valuesString = null;
	if (editor != null && editor.isApplet()) {
	    valuesString = editor.readRemoteFile(listfile);
	} else {
	    try {
		FileReader fr = new FileReader(listfile);
		StringWriter out = new StringWriter();
		int c;
		
		while ((c = fr.read()) != -1)
		    out.write(c);
		
		fr.close();
		valuesString = out.toString();
	    } catch (IOException ioe) {
		throw new MalformedScriptException("Locate listfile " 
						   + listfile + " error: " 
						   + ioe.toString());
	    }
	}
	Vector valVec = new Vector();
	int index = valuesString.indexOf('\n');
	String toAdd;
	while (index != -1) {
	    toAdd = valuesString.substring(0,index).trim();
	    if (!toAdd.equals("")) {
		valVec.add(valuesString.substring(0,index));
		valuesString = valuesString.substring(index+1);
		index = valuesString.indexOf('\n');
	    }
	}
	toAdd = valuesString.trim();
	if (!toAdd.equals("")) {
	    valVec.add(valuesString);
	}
	int mulFact = valVec.size();
	if (mulFact == 0) {
	    throw new MalformedScriptException("No values found in file "
					       + listfile + ".");
	}
	String[] values = new String[mulFact];
	for (int j = 0; j < mulFact; j++) {
	    values[j] = (String)valVec.get(j);
	}
	return new ResolvedScatterParam(name,values,mulFact);
    }
}
