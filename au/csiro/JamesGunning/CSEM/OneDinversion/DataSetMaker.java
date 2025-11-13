/**********************************************************************
DELIVERYCSEM - a program for Bayesian CSEM inversion
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

package au.csiro.JamesGunning.CSEM.OneDinversion;

import java.io.*;
import java.nio.*;
import java.util.*;
import static java.lang.Math.*;

import cern.colt.map.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;  
import cern.colt.matrix.linalg.*;
import cern.colt.list.*; 

 

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 



/** 
	Simple code to make CSEM synthetic data sets and models
	
	@author James.Gunning@csiro.au
	@version 1.0, 10/2/2008
*/


public class DataSetMaker
{

  double dumE = 1e-10;
  double dumErr = 1e-11;
  double[] offset=null;
  double[] freq=null;
  double transHeight = 30.0;

  boolean doAmp=true; 
  boolean doPhase=false,doReIm=false;

  boolean makeMod=false,makeDat=false,makeStart=false;

  int Nl;
  ModData[] mod=null;
  double start_val=0.0;
  String key=null;

  ErrorHandler EH = new ErrorHandler();

  public static void main(String argv[]) 
  {
    DataSetMaker DSM = new DataSetMaker();
	ErrorHandler _EH = DSM.EH;
	int i,n;
	String s;
	String[] opts;


	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("DataSetMaker",DSM.EH);
	options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber);
	options.addSpielText("Makes data sets and models");

	options.addOption("--model",null,1,"String","thick,rho:thick,rho...",
					  "input model-data as thick,rho pairs");
	options.addOption("--key",null,1,"String","model, startup file and data file naming string");

	options.addOption("-s",null,0,null,null,"with --model, make start file too");

	options.addOption("-f",null,1,"String","f","comma-separated frequencies");
	options.addOption("--off",null,1,"String","Xmin,dx,Xmax","offsets");
	options.addOption("--phase",null,0,null,null,"add phase");
	options.addOption("--ReIm",null,0,null,null,"do Re & Im");

	// leave this to last!
	options.addExtraHelpOption("-X","Show help with all options (there are lots more!)");


	// parse options
	options.parseOptions(argv,true);

	if (options.isOptionSet("--model"))
	  {
		DSM.makeMod=true;
		s = ((String)options.getOptionValues("--model")[0]);
		DSM.setModel(s);
	  	if (options.isOptionSet("-s")) DSM.makeStart=true;
	  }

	if (options.isOptionSet("--key"))
	  {
		DSM.key = ((String)options.getOptionValues("--key")[0]);
	  }

	if (options.isOptionSet("-f") || options.isOptionSet("--off")) DSM.makeDat = true;
 
	if (options.isOptionSet("--phase")) DSM.doPhase=true;
	if (options.isOptionSet("--ReIm")) 
	  {
		DSM.doReIm=true;
		DSM.doPhase=false;
		DSM.doAmp=false;
	  }

	if (options.isOptionSet("-f")) 
	  {	
		s = ((String)options.getOptionValues("-f")[0]);
		opts = s.split(",");
		DSM.freq = new double[opts.length];
		for(i=0;i<opts.length;i++)
		  {
			try {DSM.freq[i] = Double.parseDouble(opts[i]);}
			catch (NumberFormatException e) {_EH.die_err("-f has not a number: "+s);}
		  }
	  }
	else if (DSM.makeDat)
	  {
		_EH.die_err("-f not specified");
	  }
	

	if (options.isOptionSet("--off")) 
	  {	
		s = ((String)options.getOptionValues("--off")[0]);
		opts = s.split(",");
		n = opts.length;
		if (n!=3) _EH.die_err("--off option needs Xmin,deltaX,Xmax");
		double Xmax=0,dX=0,Xmin=0;
		try 
		  {
			Xmin = Double.parseDouble(opts[0]);
			dX = Double.parseDouble(opts[1]);
			Xmax = Double.parseDouble(opts[2]);
		  }
		catch (NumberFormatException e) {_EH.die_err("--off has not a number: "+s);}
		if (dX<=0) _EH.die_err("--off option needs deltaX>0");
		if (Xmin<0 || Xmax<0) _EH.die_err("--off option needs Xmin,Xmax>0");
		if (Xmax<Xmin) _EH.die_err("--off option needs Xmax>=Xmin");
		n = (int)rint((Xmax-Xmin)/dX)+1;
		if (n<=0) _EH.die_err("--off yields -ve offsets: n="+n);

		DSM.offset = new double[n];
		for (i=0;i<n;i++) DSM.offset[i] = Xmin+i*dX;
	  }
	else if (DSM.makeDat) _EH.die_err("--off not specified");

	if (!DSM.makeDat && !DSM.makeMod) options.Usage();

	DSM.run();
  }
  
  void run()
  {
	if (makeMod) 
	  {
		makeModelFile();
		if (makeStart) 
		  {
			makeStartupFile();
			// makeProfileFile();
		  }
	  }
	if (makeDat) makeDataFile();

  }

  void setModel(String s)
  {
	int i;
	String[] toks1,toks2;
	ModData m;

	toks1 = s.split(":");
	Nl = toks1.length;
	EH.warn("Using "+Nl+" layers");
	mod = new ModData[Nl];
	for(i=0;i<Nl;i++)
	  {
		toks2 = toks1[i].split(",");
		m = mod[i] = new ModData();
		if (toks2.length!=2) EH.die_err("In "+s+", model must be : separated thick,rho pairs");
		try 
		  {
			m.th = Double.parseDouble(toks2[0]);
			m.rho = Double.parseDouble(toks2[1]);
		  }
		catch (NumberFormatException e) {EH.die_err("-f has not a number: "+s);}
	  }
  }



  void makeModelFile()
  {
	int i,j,Nd,No,Nf,Npos;
	double rho;
	PrintWriter out;
	String filename;

	filename = modelFilename();

	ModData _m;
	// write header
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename + " for output");
	  }
	out.printf( "FORMAT:           OCCAM1DMODEL_1.0%n"
			   +"MODEL NAME:       Model%n"
			   +"NUM LAYERS:       %d%n"
			   +"THICKNESS  RESISTIVITY  PENALTY    PREJUDICE   PENALTY%n",Nl-1);
	for(i=0;i<Nl;i++)
	  {
		_m = mod[i];
		// write file for inversion, assume sea fixed
		rho = (i==0) ? _m.rho : -1;
		out.printf("%-10.5g %-10.5g	%-10.5g %-10.5g	%-10.5g%n",
				   _m.th,rho,_m.rho_pen,_m.prej,_m.pen);
	  }
	out.close();
	EH.println("Wrote file "+filename);
  }

  void makeStartupFile()
  {
	PrintWriter out;
	String filename;
	int i;

	filename = startupFilename();

	// write header
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename + " for output");
	  }
	out.printf("FORMAT:            OCCAMITER_1.0%n"+
			   "DESCRIPTION:       test%n"+
			   "MODEL FILE:        %s%n"+
			   "DATA FILE:         %s%n"+                                                                       
			   "DATE/TIME:         %s%n"+                     
			   "MAX ITER:          1000%n"+
			   "REQ TOL:           0.9%n"+
			   "IRUF:              1%n"+
			   "DEBUG LEVEL:       0%n"+
			   "ITERATION:         0%n"+
			   "PMU:               -10%n"+
			   "RLAST:             1%n"+
			   "TLAST:             1%n"+   
			   "IFFTOL:            0%n"+
			   "NO. PARMS:         %d%n",
			   modelFilename(),
			   dataFilename(),
			   (new Date()).toString(),
			   Nl-1);
	for(i=1;i<Nl;i++) out.printf("%-5.2g ",mod[i].rho);
	out.println();
	out.close();
	EH.println("Wrote file "+filename);

  }


  void makeProfileFile()
  {
	PrintWriter out;
	String filename;
	int i;

	filename = profileFilename();

	// write header
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename + " for output");
	  }
	ModData _m;
	float z=0f;
	for(i=0;i<Nl;i++)
	  {
		_m = mod[i];
		out.printf("%-10.5g %-10.5g%n",z,_m.rho);
		z += _m.th;
		out.printf("%-10.5g %-10.5g%n",z,_m.rho);
	  }
	out.close();
	EH.println("Wrote file "+filename);

  }

  void makeDataFile()
  {
	PrintWriter out;
	String filename;
	int i,j,Nd,No,Nf,Npos;

	filename = dataFilename();

	// write header
	try 
	  {
		out = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	  } 
	catch (IOException ioe) 
	  {
		throw new RuntimeException("couldn't open " + filename + " for output");
	  }

	No = offset.length;
	Nf = freq.length;
	Npos= Nf * No;
	Nd=0;
	if (doAmp) 
	  {
		Nd+=Npos;
		if (doPhase) Nd+=Npos;
	  }
	else if (doReIm) 
	  {
		Nd = 2*Npos;
	  }

	// Now write
	out.printf( "FORMAT:           OCCAM1DDATA_1.0%n"
			   +"DESCRIPTION:      CSEM data%n"
			   +"NUM DATA:         %d%n",Nd);
	out.printf("Type Datum Error range/period frq az rot theight%n");
	// frequencies
	for(i=0;i<Nf;i++)
	  {
		// offsets
		for (j=0;j<No;j++)
		  {
			if (doAmp) 
			  {
				out.printf("%d %-10.5g %-10.5g  %-10.5g   %-10.5g   %-10.5g   %-10.5g   %-10.5g%n",
						   1,dumE,dumErr,offset[j],freq[i],0.0,0.0,transHeight);
			   if (doPhase)
				 {
				   out.printf("%d %-10.5g %-10.5g  %-10.5g   %-10.5g   %-10.5g   %-10.5g   %-10.5g%n",
						   2,0.,5.,offset[j],freq[i],0.,0.,transHeight);
				 }
			  }
			else if (doReIm) 
			  {
				out.printf("%d %-10.5g %-10.5g  %-10.5g   %-10.5g   %-10.5g   %-10.5g   %-10.5g%n",
						   14,dumE,dumErr,offset[j],freq[i],0.,0.,transHeight);
				out.printf("%d %-10.5g %-10.5g  %-10.5g   %-10.5g   %-10.5g   %-10.5g   %-10.5g%n",
						   15,dumE,dumErr,offset[j],freq[i],0.,0.,transHeight);
			  }
		  }
	  }
	out.close();
	EH.println("Wrote file "+filename);
  }

  String modelFilename()
  {
	return (key==null) ? "model.mod": ("model."+key+".mod");
  }

  String startupFilename()
  {
	return (key==null) ? "startup": ("startup."+key);
  }

  String dataFilename()
  {
	return  (key==null) ? "data.dat": ("data."+key+".dat");
  }

  String profileFilename()
  {
	return  (key==null) ? "profile_final.txt": ("profile_final."+key+".txt");
  }

  static class ModData
  {
	double th,rho,rho_pen=0,prej=0,pen=1;
  }
}


