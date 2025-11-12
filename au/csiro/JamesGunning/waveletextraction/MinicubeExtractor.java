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


package au.csiro.JamesGunning.waveletextraction;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;

import au.csiro.JamesGunning.lib.OptionsParser;
import au.csiro.JamesGunning.lib.*; 
import au.csiro.JamesGunning.sulib.*; 

import edu.wlu.cs.levy.CG.KDTree;

/**
   Utility code to extract SU minicubes from 3D cubes by header word.
   Usage:

   MinicubeExtractor filename.su wellXY1.txt[,wellXY2,txt,...]  minicube_suffix.su

   Here, wellXYs.txt is a 2-column ascii file of x y values (from e.g. a survey)
   Assuming x,y,inline,crossline are in header words in filename.su, this code
   will dump to minicube_outputfile.su a minicube su file "covering" the xy set of points.
   We assume here that wellXYs are "interior" to the 3D survey.
   Default header words are gx,gy, cdp,ep, but -XY and -IC can override these. E.g.

   MinicubeExtractor filename.su wellXYs.txt minicube_outputfile.su -XY sx,sy -IC tracr,tracl

   Ouput filenames are strings of form 
*/
public class MinicubeExtractor
{
  ErrorHandler EH = new ErrorHandler();
  boolean fromLE= false,showG=false;
  String filename,filenameWell,filenameO;
  String Xword="gx",Yword="gy";
  String Iword="ep",Cword="cdp";
  int Xcol=0,Ycol=1;

  // buffer size - number of traces of "padding" in inline/crossline coords about the well
  public int buf=2;

  String[] filenameOutElimSubstrings=null;

  public static void main(String argv[]) 
  {
	String s;
	String[] sa;
	MinicubeExtractor ME = new MinicubeExtractor();
	OptionsParser options = new OptionsParser("MinicubeExtractor",ME.EH);
	options.addMandatoryArg(1,"filename.su","Big endian su file");
	options.addMandatoryArg(2,"wellXY1.txt[,wellXY2,txt,...]","2-column ascii file with x y well-trajectory points");
	options.addMandatoryArg(3,"suffix.su","destination .su file is wellXY1.txt_suffix.su");

	options.addOption("-XY",null,1,"String","gx,gy","use named XY header words");
	options.addOption("-IC",null,1,"String","ep,cdp","use named inline/crossline IC header words");
	options.addOption("-c",null,1,"String","col1,col2","column numbers for x,y in line.txt file");
	options.addOption("--LE",null,0,null,null,"Flip seismic from little endian when reading");
	options.addOption("-b",null,1,"int","buf_size","Buffer padding size in trace-spacing");
	options.addOption("-N",null,1,"String","str1,str2,...","substrings to remove from output filenames");

	// parse options
	options.parseOptions(argv,true);

	ME.filename = options.getMandatoryArg(1);
	ME.filenameWell = options.getMandatoryArg(2);
	ME.filenameO= options.getMandatoryArg(3);
	if (options.isOptionSet("-XY")) 
	  {
		s = options.getOptionString("-XY",0);
		sa = s.split(",");
		if(sa.length!=2) ME.EH.die_err("-XY needs wordX,wordY");
		ME.Xword=sa[0]; 
		ME.Yword=sa[1];
	  }
	if (options.isOptionSet("-N")) 
	  {
		ME.filenameOutElimSubstrings = options.getOptionString("-N",0).split(",");
	  }
	if (options.isOptionSet("-IC")) 
	  {
		s = options.getOptionString("-IC",0);
		sa = s.split(",");
		if(sa.length!=2) ME.EH.die_err("-IC needs inlineword,xlineword");
		ME.Iword=sa[0]; 
		ME.Cword=sa[1];
	  }
	if (options.isOptionSet("-c")) 
	  {
		s = options.getOptionString("-c",0);
		sa = s.split(",");
		if(sa.length!=2) ME.EH.die_err("-c needs col1,col2");
		ME.Xcol=Integer.parseInt(sa[0])-1; 
		ME.Ycol=Integer.parseInt(sa[1])-1;
	  }
	if (options.isOptionSet("--LE")) ME.fromLE = true;
	if (options.isOptionSet("-b")) 
	  {
		ME.buf = options.getOptionInt("-b",0);
		if(ME.buf<1) ME.EH.die_err("buffer must be >=1");
	  }

	ME.run();
  }
  
  void run()
  {
	int i,j,k,X,Y,I,C,dI=0,dC=0,Iold=0,Cold=0;
	double scalco=1;
	SuReader sureader = new SuReader();
	SuSegyFormat tr = new SuSegyFormat();
	FileInputStream f=null;
	BufferedReader in=null;
	try
	  {
		f = (FileInputStream) new FileInputStream(filename);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open su file "+filename);
	  }

	// kdtree from seismic {x,y}s
	KDTree<Integer> kd = new KDTree<Integer>(2);
	ArrayList<double[]> trxys = new ArrayList<double[]>();
	ArrayList<int[]> trICs = new ArrayList<int[]>();
	
	EH.println("Reading seismic "+filename);
	try
	  {
		for(i=0,j=0;;i++)
		  {
			if (sureader.fgettr(f,tr,fromLE)==0) break;
		  
			X=tr.get_int(Xword);
			Y=tr.get_int(Yword);
			I=tr.get_int(Iword);
			C=tr.get_int(Cword);
			if (i!=0) 
			  {
				if (dI==0) if ((k=abs(I-Iold))!=0) dI=k;
				if (dC==0) if ((k=abs(C-Cold))!=0) dC=k;
				if (I==Iold && C==Cold) continue; // e.g. angle gathers
			  }
			else
			  {
				scalco = (double)tr.get_short("scalco");
				if (scalco<0) scalco = 1.0/abs(scalco);
				if (scalco==0.0) scalco = 1.0;
			  }
		
			double[] keys = new double[]{X*scalco,Y*scalco};
			try
			  {
				kd.insert(keys,j);
				trxys.add(keys); 
				trICs.add(new int[]{I,C});
				Iold=I; Cold=C;
				j++;
			  }
			catch (edu.wlu.cs.levy.CG.KeyDuplicateException e)
			  {
				continue;
				// ignore these (eg angle gathers)
			  }
			if(i%100000==0) EH.printf("Read unique XY trace %d%n",i);
		  }
		f.close();
	  }
	catch (Exception e) 
	  {
		EH.println("Error: "+e);
	  }
	EH.printf("Inline/xline increments %d %d%n",dI,dC);


 	EH.println("Reading well XYs");
	// read files of well xys 
	String s;
	String[] sa,sw;
	int[] IC;
	double[] pt;
	sw = filenameWell.split(",");
	int Nwells = sw.length;
	int[] Imin=new int[Nwells];
	int[] Imax=new int[Nwells];
	int[] Cmin=new int[Nwells];
	int[] Cmax=new int[Nwells];
	
	for(i=0;i<Nwells;i++)
	  {
		try
		  {
			in = new BufferedReader(new FileReader(sw[i]));
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't open file "+sw[i]); 
		  }
		double[] xy,xyt;
		ArrayList<double[]> xys = new ArrayList<double[]>();
		try
		  {
			while ((s = in.readLine()) != null)
			  {
				s = s.trim();
				if(s.startsWith("#")) continue;
				if (s.equals("")) continue;
				sa = s.split("\\s+");     // split just at white space
				xy = new double[2];
				xy[0] = Double.parseDouble(sa[Xcol]);
				xy[1] = Double.parseDouble(sa[Ycol]);
				xys.add(xy);
			  }
		  }
		catch (Exception e) 
		  {
			EH.die_err("Can't read file  "+sw[i]); 
		  }
		xys.trimToSize();
	  
		// form inline/xline bounding box
		EH.println("Computing bounding box "+sw[i]);

		Cmin[i]=Imin[i]=Integer.MAX_VALUE;
		Cmax[i]=Imax[i]=-Integer.MAX_VALUE;
		try
		  {
			for(j=0;j<xys.size();j++)
			  {
				pt = xys.get(j);		   
				k = kd.nearest(pt);
				//EH.printf("%d\n",k);
				IC = trICs.get(k);	
				Imin[i] = min(Imin[i],IC[0]);
				Imax[i] = max(Imax[i],IC[0]);
				Cmin[i] = min(Cmin[i],IC[1]);
				Cmax[i] = max(Cmax[i],IC[1]);
			  }
		  }
		catch (edu.wlu.cs.levy.CG.KeySizeException e)
		  {
			EH.die_err("KeySizeException "+e);
		  }

		Imin[i]-=buf*dI;
		Imax[i]+=buf*dI;
		Cmin[i]-=buf*dC;
		Cmax[i]+=buf*dC;

		// Minicube in/cross box
		EH.printf("BBbox %d %d %d %d%n",Imin[i],Imax[i],Cmin[i],Cmax[i]);
	  }

 	EH.println("Assembling minicube (second pass read)");
	try
	  {
		f = (FileInputStream) new FileInputStream(filename);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open su file "+filename);
	  }
	try
	  {
		// write select traces (nearest to pts on line) out
		FileOutputStream[] out = new FileOutputStream[Nwells];
		String[] fout_names = new String[Nwells];
		for(j=0;j<Nwells;j++) 
		  {
			// eg. wellname.su
			fout_names[j] = sw[j]+"."+filenameO;
			if (filenameOutElimSubstrings!=null && filenameOutElimSubstrings.length>0)
			  {
				for(k=0;k<filenameOutElimSubstrings.length;k++)
				  {
					fout_names[j] = fout_names[j].replaceAll(filenameOutElimSubstrings[k],"");
				  }
			  }
			out[j] = new FileOutputStream(fout_names[j],false);
		  }
		SuWriter suwriter = new SuWriter();
		i=0;
		while (sureader.fgettr(f,tr,fromLE) != 0)
		  {
			I=tr.get_int(Iword);
			C=tr.get_int(Cword);
			for(j=0;j<Nwells;j++)
			  {
				if (I>=Imin[j] && I<=Imax[j]  && C>=Cmin[j]  && C<=Cmax[j])
				  {
					suwriter.fputtr(out[j],tr);
				  }
			  }
			if(i%100000==0) EH.printf("Read trace %d%n",i);
			i++;
		  }
		for(j=0;j<Nwells;j++)
		  {
			suwriter.flush(out[j]); // must do this!
			out[j].close();
			EH.printf("Wrote %s%n",fout_names[j]);
		  }
		f.close();
	  }
	catch (Exception e) 
	  {
		System.err.println(e);	 
	  }
  }
}