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

import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.MathUtil;

import cern.colt.matrix.impl.*;
import cern.colt.map.OpenDoubleIntHashMap;

/**
   Data class for storing seismic stack data for minicubes and allied information
   
   @author James.Gunning@csiro.au
   @version 1.0, 12/2/2003
*/


public class SeismicStackInfo
{
  /**
	 If we use header words like gx, sx etc, sometimes the spacing is not perfect...
	 Let's tolerate this much slop.
  */
  public static final double TOLERABLE_SPACING_IRREGULARITY_EPSILON = 0.05;

  /**
	 If the registraion time is uncertain, this is the allowable maximum of this
	 when specified as a standard deviation (seconds)
  */
  public static final double MAXIMUM_SYSTEMATIC_TIME_SCALE_SECS = 0.5;
  /**
	 Number of standard deviations safety margin at edges of cube
	 if registration errors are modelled
  */
  public static final double SYSTEMATIC_T_SHIFT_SAFETY_MARGIN = 10;

  /**
	Factors that control auto-detected of powers of ten multipliers in
	seismic x,y header words Auto detection is disabled
	if AUTO_DETECT_SEISMIC_XY_WORD_UNITS = false.  We disable this, as
	scaling is more robust if we can assume the scalco header word is
	supplied.
  */
  public static final boolean AUTO_DETECT_SEISMIC_XY_WORD_UNITS = false;
  /**
	 Ratio of well midpoint x to header-word average |x| is likely to be several orders of magnitude
	 if units have been changed. Sets the size of this detection threshold
  */
  public static final double AUTO_DETECT_SEISMIC_XY_WORD_THRESHOLD = 5.0;
  /** 
	  And if there is a units issue, we assume it's a power of ten. This epsilon
	  is a tolerance for establishing the nearest power of ten
  */
  public static final double AUTO_DETECT_SEISMIC_XY_WORD_POWER_EPS = 0.1;


  /**
	 Filename of data
  */
  String filename;

  // Stack parameters
  // specified stack angle if constant stack. Converted to radians as soon as parsed in degrees
  double theta_specified = -Double.MAX_VALUE;       
  double sinsq_theta_specified = -Double.MAX_VALUE;  // specified sin^2(theta) of stack angle theta if constant stack
  double offsetXmin;
  double offsetXmax;
  double T_event;
  double V_stack;
  double thetaApprox=0;  // approximate stack angle
  boolean isPSdata = false;
  private double V_denom,sqrt_V_denom;

  /**
	 Property for modelling the systematic time shift;
  */
  Property systematicTProp;
  

  // number of samples in SU minicube traces
  private int ns;
  // first time sample in SU minicube traces
  private double f1;
  // first sample time of data in stored cube data[][][]
  double tStartData=0;
  
  /**
	 amplitude data in seismic minicube
  */
  au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D data;

  /**
	 Backup'ed amplitude data in seismic minicube used in bootstrapping 
  */
  au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D dataBk;

  /**
	 Array of x,y's as they are read in
  */
  ArrayList XYValues;
  /**
	 Minicube parameters. cube occupies coordinates
	 Xmin, Xmin+inlineSpacing,.....Xmax etc
  */
  int inline_max;
  int inline_min;
  int xline_max;
  int xline_min;
  int inlineSpacing;
  int xlineSpacing;

  char fastDirection = 'i';
  /*
	Set to true if we use 1D snapping
  */
  boolean is1DdataAndUseSnapping=false;
  /*
	Set to true if we use 2D snapping
  */
  boolean is2DdataAndUseSnapping=false;
  /*
	Set to true if we use 2D snapping
  */
  boolean is2DdataAndUseXYwords=true;

  /**
	 Local 2D utm info if needed
  */
  UTMConversionInfo utmInfo=null;

  // Parameters of coordinate rotation between inline,xline and x,y
  private double X0,Y0;
  double mean_dInline, mean_dXline;
  private double inLine,xLine;

  private double A11,A12,A21,A22;

  /**
	 Constructor. Precompute V_denom and other useful things at construction.
	 Xmin/max, T_event and V_stack must be in consistent units
  */
  public SeismicStackInfo(String filename, double offsetXmax, double offsetXmin,double T_event,double V_stack) 
  {

   this.filename = filename;

   if (offsetXmin>offsetXmax)
	 throw new RuntimeException("Stack must have offsetXmin<=offsetXmax");
   else
	 {
	   this.offsetXmax=offsetXmax;
	   this.offsetXmin=offsetXmin;
	 }

   this.T_event = T_event;
   this.V_stack = V_stack;

   // time shift property 
   systematicTProp = new Property(0.0, false, false);

   computeStackParams();
  }

  void computeStackParams()
  {
	double Xms;
	// mean square offset <X^2>
	if (offsetXmax > offsetXmin)
	  {
		Xms = (MathUtil.cube(offsetXmax) - MathUtil.cube(offsetXmin))/(3.0*(offsetXmax - offsetXmin));
	  }
	else Xms = (offsetXmax*offsetXmax);

	// denominator of coefficient of theta^2
	// DIX equation   <theta^2> = <X^2> Vint^2/(Vrms^4 Tevent^2)
	// Vdenom = (Vrms^4 Tevent^2)/<X^2>
	if (Xms!=0)
	  {
		V_denom = MathUtil.sqr(MathUtil.sqr(V_stack)) * MathUtil.sqr(T_event)/Xms;
		sqrt_V_denom = Math.sqrt(V_denom);
	  }
	else 
	  {
		V_denom = sqrt_V_denom = Double.MAX_VALUE;
	  }

	if (theta_specified < 0)
	  {
		// approximate stack angle
		thetaApprox =  Math.sqrt(Xms/(MathUtil.sqr(V_stack) * MathUtil.sqr(T_event)));
	  }
	else
	  {
		// compute this useful guy
		sinsq_theta_specified = MathUtil.sqr(Math.sin(theta_specified));		
	  }
  }
  
  public double getVdenom()
  {
	return V_denom;
  }

  public double getsqrtVdenom()
  {
	return sqrt_V_denom;
  }

  public boolean isStackAnglePrescribed()
  {
	return (theta_specified >= 0);
  }

  public boolean isPSdata()
  {
	return isPSdata;
  }

  void loadSeismicData(WaveletExtractor WE, int wellNumber)
  {
	int i,j,k,kk;
	int trace_number;
	int ns_this;
	int inline=0,xline=0,inlineOld=0,xlineOld=0;
	int Ninline,Nxline,NzFirst=0,NzLast=0,Nz;
	int inline_first=0,xline_first=0,inline_last=0,xline_last=0;
	int dInline=0,dXline=0,Dinline=0,Dxline=0,deltaInline=0,deltaXline=0;
	double scalco=1.0;
	double X,Y,tmp,tStartCube,tEndCube;
	int count=0;
	boolean gotTrace;
	double eps;
	double OEx,OEy,ECx,ECy,xProduct,cubeXYangle,Dx_c,Dy_c,diagonalSpan;
	double inLineSpan,xLineSpan;
	RandomAccessFile in_strm;
	UTMConversionInfo.UTMxyValues seisCoordsFromUTM=null;

	// Open input file
	try 
	  {
		in_strm = new RandomAccessFile(filename,"r");
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " +filename);
	  }
	// check inline, xline words
	if (WE.tr.headerMap.get(WE.seismicInlineWord) == null)
	  {
		WE.EH.die_err(WE.seismicInlineWord +" is not a valid SU header word");
	  }
	if (WE.tr.headerMap.get(WE.seismicXlineWord) == null)
	  {
		WE.EH.die_err(WE.seismicInlineWord +" is not a valid SU header word");
	  }
	// check x y words
	if (WE.useUTMconversion)
	  {
		WE.seismicXword = WE.seismicInlineWord;
		WE.seismicYword = WE.seismicXlineWord;
	  }
	else
	  {
		if (WE.tr.headerMap.get(WE.seismicXword) == null)
		  {
			WE.EH.die_err(WE.seismicXword +" is not a valid SU header word");
		  }
		if (WE.tr.headerMap.get(WE.seismicYword) == null)
		  {
			WE.EH.die_err(WE.seismicYword +" is not a valid SU header word");
		  }
	  }



	// read the first 2 traces to determine fastest direction
	for (trace_number=0;trace_number<=1;trace_number++)
	  {
		gotTrace = (WE.sureader.fgettr(in_strm,WE.tr)!=0);
		if (!gotTrace) break;
		if (trace_number==0) inline = inline_first = get_value(WE.tr,WE.seismicInlineWord);
		if (trace_number==1) inline = get_value(WE.tr,WE.seismicInlineWord);
	  }
	if (inline != inline_first)
	  {
		fastDirection = 'i';
	  }
	else
	  {
		WE.EH.println("        Taking "+ WE.seismicXlineWord+" as fastest direction in minicube: "+filename);
		fastDirection = 'x';
		if (false)
		  {
			// we need to swap the meanings of inline/xline. Code expects inline fastest
			String tmpWord = new String(WE.seismicInlineWord);
			WE.seismicInlineWord = WE.seismicXlineWord;
			WE.seismicXlineWord = tmpWord;
			if (WE.useUTMconversion)
			  {
				WE.seismicXword = WE.seismicInlineWord;
				WE.seismicYword = WE.seismicXlineWord;
			  }
		  }
	  }

	// Rewind. Read the traces. We have to get system sizes, so this must be two pass
	// data must exist on a regular spaced lattice.
	// minicube must exist on lattice *aligned* with inline,xline
	try 
	  {
		in_strm.seek(0);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't rewind file ptr to 0 in " +filename);
	  }
	// OK find biggest cube;
	inline_min = xline_min = Integer.MAX_VALUE;
	inline_max = xline_max = -Integer.MAX_VALUE;

	// XY values stored here, inline the fastest direction
	XYValues = new ArrayList();
	for (trace_number=0;;trace_number++)
	  {
		gotTrace = (WE.sureader.fgettr(in_strm,WE.tr)!=0);
		if (!gotTrace) break;

		// assume f1 is constant. This could also be specified by delrt
		// assume dt is constant
		if (trace_number==0) 
		  {
			f1 = (double) WE.tr.get_float("f1");
			if (WE.tr.get_short("delrt") != 0)
			  {
				// looks like delrt field defines window start (ms)
				f1 = WE.tr.get_short("delrt")*1e-3; // to secs
			  }
			// x,y scale factor
			scalco = (double)WE.tr.get_short("scalco");
			if (scalco<0) scalco = 1.0/Math.abs(scalco);
			if (scalco==0.0) scalco = 1.0;
			// dt
			WE.dt = WE.tr.get_short("dt")/1.0e6;
			if (WE.dt<=0) WE.EH.die_err(filename+": dt not set correctly");

			// BBANDTEST_ISSUES
		    // if (WE.assumeNoReflectionsOutsideLogs) //OLD. Not needed
			if (true)
			  {
				// may need safe data margin
				// tStartCube = WE.tStart - WE.tCodaMax;
				// tEndCube = WE.tEnd - WE.tPrecursorMin;
				// 4pt Lagrange interpolation means we need samples -1,+2 above window also
				tStartCube = WE.tStart - WE.dt;
				tEndCube = WE.tEnd + 2*WE.dt;
			  }
			else
			  {
				// old
				tStartCube = WE.tStart;
				tEndCube = WE.tEnd;
			  }
			
			// additional margin for registration errors
			if(systematicTProp.mean != 0.0)
			  {
				tStartCube += systematicTProp.mean;
				tEndCube += systematicTProp.mean;
			  }
			if (systematicTProp.sigma>0)
			  {
				tStartCube -= SYSTEMATIC_T_SHIFT_SAFETY_MARGIN*systematicTProp.sigma;
				tEndCube += SYSTEMATIC_T_SHIFT_SAFETY_MARGIN*systematicTProp.sigma;
			  }
			NzLast = (int)Math.ceil((tEndCube - f1)/WE.dt);
			if (NzLast<0) WE.EH.die_err(filename+": time range of data does not extend to time "+tEndCube);
			NzFirst = (int)Math.floor((tStartCube - f1)/WE.dt);
			if (NzFirst<0) WE.EH.die_err(filename+": time range of data does not extend to time "+tStartCube);
			tStartData = f1+NzFirst*WE.dt;
			inline_first = get_value(WE.tr,WE.seismicInlineWord);
			xline_first = get_value(WE.tr,WE.seismicXlineWord);
			if (WE.useUTMconversion)
			  {
				seisCoordsFromUTM = WE.utmInfo.getUTMSfromInlineXline(inline_first,xline_first);
				X0 = seisCoordsFromUTM.x;
				Y0 = seisCoordsFromUTM.y;
			  }
			else
			  {
				X0 = get_value(WE.tr,WE.seismicXword)*scalco;
				Y0 = get_value(WE.tr,WE.seismicYword)*scalco;
			  } 
		  }
		// check ns is sane
		ns_this = (int)WE.tr.get_short("ns");
		if (ns==0) ns=ns_this;
		if (ns_this!=ns) WE.EH.die_err(filename+": ns is not regularly sampled");
 
		// get coords
		inline = get_value(WE.tr,WE.seismicInlineWord);
		xline = get_value(WE.tr,WE.seismicXlineWord);
		inline_min = Math.min(inline_min,inline);
		xline_min = Math.min(xline_min,xline);
		inline_max = Math.max(inline_max,inline);
		xline_max = Math.max(xline_max,xline);
		if (WE.useUTMconversion)
		  {
			seisCoordsFromUTM = WE.utmInfo.getUTMSfromInlineXline(inline,xline);
			X = seisCoordsFromUTM.x;
			Y = seisCoordsFromUTM.y;
		  }
		else
		  {
			X = get_value(WE.tr,WE.seismicXword)*scalco;
			Y = get_value(WE.tr,WE.seismicYword)*scalco;
		  } 

		// Store X,Y's
		XYValues.add(new SeismicXYValues(X,Y));

		// System.out.println("POS: "+ inline+" "+xline+" "+X+" "+Y);

		// successive inline,xline increments must be "regular". Only permissible increments
		// are 0, dInline, Dinline
		if (trace_number>0)
		  {
			if ((deltaInline = inline-inlineOld) !=0)
			  {
				if (trace_number==1) fastDirection = 'i';
				if (deltaInline!=0)
				  {
					if (dInline==0) dInline = deltaInline;                  // first increment
					if (deltaInline != dInline && Dinline==0) Dinline = deltaInline;  // first end-of line back sweep
					// henceforth, must have  deltaInline=dInline || deltaInline=Dinline
					if (deltaInline != dInline && deltaInline != Dinline)
					  {
						eps = Math.abs((deltaInline-dInline)/dInline);
						if (eps > TOLERABLE_SPACING_IRREGULARITY_EPSILON)
						  {
							WE.EH.die_err("Seismic minicube " + filename +" is not regularly spaced in inline direction. "+
										  "To see successive differences, try:%n sugethw verbose=0 key="+
										  WE.seismicInlineWord+" < "+filename + 
										  " | grep -v '^$' | awk '{gsub(/"+WE.seismicInlineWord+"=/,\"\");print $1-p;p=$1}'");
						  }
					  }
				  }
			  }
			if ((deltaXline = xline-xlineOld) !=0)
			  {
				if (trace_number==1) fastDirection = 'x';
				if (deltaXline!=0)
				  {
					if (dXline==0) dXline = deltaXline;                  // first increment
					if (deltaXline != dXline && Dxline==0) Dxline = deltaXline;  // first end-of line back sweep
					// henceforth, must have deltaXline=dXline || deltaXline=Dxline
					if (deltaXline != dXline && deltaXline != Dxline)
					  {
						WE.EH.die_err("Seismic minicube " + filename +" is not regularly spaced in xline direction. "+
									  "To see successive differences, try:%n sugethw verbose=0 key="+
									  WE.seismicXlineWord+" < "+filename + 
									  " | grep -v '^$' | awk '{gsub(/"+WE.seismicXlineWord+"=/,\"\");print $1-p;p=$1}'");
					  }
				  }
			  }
		  }
		inlineOld=inline;
		xlineOld=xline;
	  }
	inline_last= inline;
	xline_last= xline;
	XYValues.trimToSize();

	if (fastDirection != 'i' && (XYValues.size() > 1))
	  {
		// this should never happen on account of pre-emptive determination of fastest direction above
		WE.EH.warn("        N.B: Seismic minicube " + filename +" does not have fast direction corresponding to inline word "
				   + WE.seismicInlineWord);
	  }

	// Checks
	// 1) Does the user specify the dimension correctly?
	if ((inline_max==inline_min) && (xline_max==xline_min))
	  {
		// Data set is single trace
		if (!is1DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is single trace (according to specified keywords). Switch on \"is-single-trace-only\"" +
									   " XML attribute if trace is OK");	
		  }
		if (is2DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a trace (according to specified keywords). Remove two-D-utm-conversions entry for 2D line");	
		  }
	  }
	else if ((inline_max==inline_min) || (xline_max==xline_min))
	  {
		// 2D
		if (is1DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 2D line (according to specified keywords). Switch off \"is-single-trace-only\" XML attribute entry");
		  }
		if (!is2DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 2D line (according to specified keywords). Add in \"two-D-utm-conversions\" entry for 2D line");	
		  }
	  }
	else
	  {
		if (is1DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 3D cube (according to specified keywords). Switch off \"is-single-trace-only\" XML entry");
		  }
		if (is2DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 3D cube (according to specified keywords). Remove two-D-utm-conversions entry for 2D line");	
		  }
	  }
	
	// Is 2D utm spec OK?
	if (is2DdataAndUseSnapping && is2DdataAndUseXYwords)
	  {
		utmInfo = new UTMConversionInfo();

		utmInfo.x[0] = ((SeismicXYValues)XYValues.get(0)).x;
		utmInfo.y[0] = ((SeismicXYValues)XYValues.get(0)).y;
		utmInfo.inLine[0] = inline_first;
		utmInfo.xLine[0] = xline_first;
		
		utmInfo.x[1] = ((SeismicXYValues)XYValues.get(XYValues.size()-1)).x;
		utmInfo.y[1] = ((SeismicXYValues)XYValues.get(XYValues.size()-1)).y;
		utmInfo.inLine[1] = inline_last;
		utmInfo.xLine[1] = xline_last;

		if (fastDirection == 'i')
		  {
			if (utmInfo.xLine[0] != utmInfo.xLine[1])
			  {
				throw new RuntimeException("2D data set doesn't have a common xLine in: " + filename);
			  }
		  }
		else
		  {
			if (utmInfo.inLine[0] != utmInfo.inLine[1])
			  {
				throw new RuntimeException("2D data set doesn't have a common inLine in: " + filename);
			  }
		  }

		// set the third point by trivial 90deg rotation, preserving squareness of geometry
		// Points then end up like this:          2    1
		//                                          0
		utmInfo.x[2] = utmInfo.x[0]+utmInfo.y[0]-utmInfo.y[1];
		utmInfo.y[2] = utmInfo.y[0]+utmInfo.x[1]-utmInfo.x[0];
		utmInfo.inLine[2] = utmInfo.inLine[0]+utmInfo.xLine[0]-utmInfo.xLine[1];
		utmInfo.xLine[2]  = utmInfo.xLine[0]+utmInfo.inLine[1]-utmInfo.inLine[0];
		// setup UTMS and return false if they're bogus
		if (!utmInfo.setUpUTMmapping(WE.EH))
		  {
			throw new RuntimeException("2D data set cannot be used to set geometry property from header words: " + filename);
		  }
	  }


	// Allocate required storage for SU minicube. Data will be stored in 
	// matrix with sense [inline][xline][t]
	Ninline = 1 + (inline_max-inline_min)/(dInline!=0 ? Math.abs(dInline) : 1);
	Nxline = 1 + (xline_max-xline_min)/(dXline!=0 ? Math.abs(dXline) : 1);
	Nz = NzLast - NzFirst+1;

	//	data = new DenseDoubleMatrix3D(Ninline,Nxline,Nz);  JGDENSEHACK
	data = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D(Ninline,Nxline,Nz);

	// second pass over data. Read it and stuff into data cube
	// set file ptr to start
	try 
	  {
		in_strm.seek(0);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't rewind file ptr to 0 in " +filename);
	  }
	X0 = X = ((SeismicXYValues)XYValues.get(0)).x;
	Y0 = Y = ((SeismicXYValues)XYValues.get(0)).y;
	SeismicXYValues	XY;
	for (trace_number=0;;trace_number++)
	  {
		gotTrace = (WE.sureader.fgettr(in_strm,WE.tr)!=0);
		if (!gotTrace) break;
		// get coords
		inline = get_value(WE.tr,WE.seismicInlineWord);
		xline = get_value(WE.tr,WE.seismicXlineWord);
		// matrix indices
		i = (inline-inline_min)/((dInline!=0) ? Math.abs(dInline) : 1);
		j = (xline-xline_min)/(dXline!=0 ? Math.abs(dXline) : 1);
		if (i<0 || i>=Ninline)  WE.EH.die_err("inline indx out of range");
		if (j<0 || j>=Nxline)  WE.EH.die_err("xline indx out of range");

		// Store X,Y's in corresponding unrolled indexed location of [i][j]
		// but with inline the fastest direction. On completion of the loop
		// XYValues is stored *as if* 
		// inline were the fastest direction, in both cases.
		if (WE.useUTMconversion)
		  {
			seisCoordsFromUTM = WE.utmInfo.getUTMSfromInlineXline(inline,xline);
			X = seisCoordsFromUTM.x;
			Y = seisCoordsFromUTM.y;
		  }
		else
		  {
			X = get_value(WE.tr,WE.seismicXword)*scalco;
			Y = get_value(WE.tr,WE.seismicYword)*scalco;
		  } 
		// unrolled [j][i] index
		XY = (SeismicXYValues)XYValues.get(j*Ninline+i);
		XY.x = X;
		XY.y = Y;


		// and samples
		k = (int)Math.rint((tStartData-f1)/WE.dt);
		for (kk=0;kk<Nz;kk++,k++) 
		  {
			data.setQuick(i,j,kk,(double)WE.tr.get_data_value(k));
		  }

	  }

	try 
	  {
		in_strm.close();
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " +filename);
	  }

	// set object fields
	inlineSpacing = Math.abs(dInline);
	xlineSpacing = Math.abs(dXline);

	// Compute average dX, dY values (to prevent flt. pt error)

	SeismicXYValues originPt = (SeismicXYValues)XYValues.get(0);
	SeismicXYValues oppositeCornerPt = (SeismicXYValues)XYValues.get(XYValues.size()-1);
	SeismicXYValues endOfFirstSweepPt=null;
	endOfFirstSweepPt	= (SeismicXYValues)XYValues.get(Ninline-1);

	X0 = X = originPt.x;
	Y0 = Y = originPt.y;

	// length of minicube side along inline direction from Pythagoras...
	OEx = endOfFirstSweepPt.x - X0;
	OEy = endOfFirstSweepPt.y - Y0;

	// inline spans directly from inline side of minicube. Fast but
	// sensitive to flt pt error of endOfFirstSweepPt point.
	inLineSpan = Math.sqrt(OEx*OEx+OEy*OEy);
	mean_dInline = (Ninline>1) ? inLineSpan/(Ninline-1) : 0.0;


	ECx = oppositeCornerPt.x - endOfFirstSweepPt.x;
	ECy = oppositeCornerPt.y - endOfFirstSweepPt.y;

	// Compute mean_dXline from Pythagoras over full diagonal span of minicube
	Dx_c = oppositeCornerPt.x - X0;
	Dy_c = oppositeCornerPt.y - Y0;
	diagonalSpan = Math.sqrt(Dx_c*Dx_c+Dy_c*Dy_c);
	xLineSpan = Math.sqrt(ECx*ECx+ECy*ECy);
	cubeXYangle = Math.toDegrees(Math.acos((inLineSpan*xLineSpan>0) ? 
										 (inLineSpan*inLineSpan+xLineSpan*xLineSpan-diagonalSpan*diagonalSpan)/
										 (2.0*inLineSpan*xLineSpan) : 0));
	mean_dXline = (Nxline>1) ? xLineSpan/(Nxline-1) : 0.0;
	
	if (WE.verbosity>2)
	  {
		WE.EH.println("        Minicube: "+filename +", Nxline="+Nxline+" Ninline="+Ninline);
		WE.EH.println("                  trace spacing inline: %-2g crossline: %-2g",
					  mean_dInline, mean_dXline);
		WE.EH.println("                  minicube x/y axis angle (degrees): %-3g", cubeXYangle);
		if (is1DdataAndUseSnapping) WE.EH.println("                  Using 1D minicube (single-trace) snapping");
		if (is2DdataAndUseSnapping) WE.EH.println("                  Using 2D minicube (seismic line) snapping");
	  }

	// Now set up rotation geometry (which may be trivial)

	// use three points at (inline,xline) = 
	// (inline_min,xline_min),(inline_max,xline_min),(inline_max,xline_max)
	// so form linear system for
	// (inLine-inline_min, xline-xline_min) = A.(X-X0,Y-Y0)
	
	ECx = endOfFirstSweepPt.x - X0;
	ECy = endOfFirstSweepPt.y - Y0;
	
	// protection not really required, but be safe.
	double tmp1,tmp2;
	tmp1 = Dy_c*ECx - Dx_c*ECy;
	tmp2 = -(Dy_c*ECx) + Dx_c*ECy;
	// entries in A matrix
	A11 = (tmp1!=0) ? (((Dy_c - ECy)*(inline_max - inline_min))/tmp1) : 0;
	A12 = (tmp2!=0) ? ((Dx_c - ECx)*(inline_max - inline_min))/tmp2: 0;
	A21 = (tmp2!=0) ? ((ECy*(xline_max - xline_min))/tmp2) : 0;
	A22 = (tmp1!=0) ? (ECx*(xline_max - xline_min))/tmp1 : 0;

	// Silly checks
	if (!true)
	  {
		System.out.println("Seismic transformation checks....");
		computeIndices(oppositeCornerPt.x,oppositeCornerPt.y);
		System.out.println("Corner pt inline xline = "+inLine+" "+xLine);
		computeIndices(X0,Y0);
		System.out.println("First trace inline xline = "+inLine+" "+xLine);

		SeismicXYValues pt = (SeismicXYValues)XYValues.get(1);
		computeIndices(pt.x,pt.y);
		System.out.println("Second trace inline xline = "+inLine+" "+xLine);

		pt = (SeismicXYValues)XYValues.get(Ninline-1);
		computeIndices(pt.x,pt.y);
		System.out.println(Ninline +"-th trace inline xline = "+inLine+" "+xLine);

		pt = (SeismicXYValues)XYValues.get(XYValues.size()-Ninline);
		computeIndices(pt.x,pt.y);
		System.out.println((XYValues.size()-Ninline+1) +"-th trace inline xline = "+inLine+" "+xLine);

		System.exit(0);
	  }
  }

  /*
	This is not used (its the old version), but we'll leave it here until any bugs in the new
	one get sorted out..
  */
  void OLD_loadSeismicData(WaveletExtractor WE, int wellNumber)
  {
	int i,j,k,kk;
	int trace_number;
	int ns_this;
	int inline=0,xline=0,inlineOld=0,xlineOld=0;
	int Ninline,Nxline,NzFirst=0,NzLast=0,Nz;
	int inline_first=0,xline_first=0,inline_last=0,xline_last=0;
	int dInline=0,dXline=0,Dinline=0,Dxline=0,deltaInline=0,deltaXline=0;
	double scalco=1.0;
	double X,Y,tmp,tStartCube,tEndCube;
	int count=0;
	boolean gotTrace;
	double eps;
	double OEx,OEy,ECx,ECy,xProduct,cubeXYangle,Dx_c,Dy_c,diagonalSpan;
	double inLineSpan,xLineSpan;
	char fastDirection = 'i';
	RandomAccessFile in_strm;

	// Open input file
	try 
	  {
		in_strm = new RandomAccessFile(filename,"r");
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open " +filename);
	  }
	// check inline, xline words
	if (WE.tr.headerMap.get(WE.seismicInlineWord) == null)
	  {
		WE.EH.die_err(WE.seismicInlineWord +" is not a valid SU header word");
	  }
	if (WE.tr.headerMap.get(WE.seismicXlineWord) == null)
	  {
		WE.EH.die_err(WE.seismicInlineWord +" is not a valid SU header word");
	  }
	// check x y words
	if (WE.useUTMconversion)
	  {
		WE.seismicXword = WE.seismicInlineWord;
		WE.seismicYword = WE.seismicXlineWord;
	  }
	else
	  {
		if (WE.tr.headerMap.get(WE.seismicXword) == null)
		  {
			WE.EH.die_err(WE.seismicXword +" is not a valid SU header word");
		  }
		if (WE.tr.headerMap.get(WE.seismicYword) == null)
		  {
			WE.EH.die_err(WE.seismicYword +" is not a valid SU header word");
		  }
	  }



	// read the first 2 traces to determine fastest direction
	for (trace_number=0;trace_number<=1;trace_number++)
	  {
		gotTrace = (WE.sureader.fgettr(in_strm,WE.tr)!=0);
		if (!gotTrace) break;
		if (trace_number==0) inline = inline_first = get_value(WE.tr,WE.seismicInlineWord);
		if (trace_number==1) inline = get_value(WE.tr,WE.seismicInlineWord);
	  }
	if (inline != inline_first)
	  {
		fastDirection = 'i';
	  }
	else
	  {
		WE.EH.println("        Taking "+ WE.seismicXlineWord+" as fastest direction in minicube: "+filename);
		fastDirection = 'x';
		// we need to swap the meanings of inline/xline. Code expects inline fastest
		String tmpWord = new String(WE.seismicInlineWord);
		WE.seismicInlineWord = WE.seismicXlineWord;
		WE.seismicXlineWord = tmpWord;
		if (WE.useUTMconversion)
		  {
			WE.seismicXword = WE.seismicInlineWord;
			WE.seismicYword = WE.seismicXlineWord;
		  }
	  }

	// Rewind. Read the traces. We have to get system sizes, so this must be two pass
	// data must exist on a regular spaced lattice.
	// minicube must exist on lattice *aligned* with inline,xline
	try 
	  {
		in_strm.seek(0);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't rewind file ptr to 0 in " +filename);
	  }
	// OK find biggest cube;
	inline_min = xline_min = Integer.MAX_VALUE;
	inline_max = xline_max = -Integer.MAX_VALUE;

	// XY values stored here, inline the fastest direction
	XYValues = new ArrayList();
	for (trace_number=0;;trace_number++)
	  {
		gotTrace = (WE.sureader.fgettr(in_strm,WE.tr)!=0);
		if (!gotTrace) break;

		// assume f1 is constant. This could also be specified by delrt
		// assume dt is constant
		if (trace_number==0) 
		  {
			f1 = (double) WE.tr.get_float("f1");
			if (WE.tr.get_short("delrt") != 0)
			  {
				// looks like delrt field defines window start (ms)
				f1 = WE.tr.get_short("delrt")*1e-3; // to secs
			  }
			// x,y scale factor
			scalco = (double)(WE.tr.get_short("scalco"));
			if (scalco<0) scalco = 1.0/Math.abs(scalco);
			if (scalco==0.0) scalco = 1.0;
			// dt
			WE.dt = WE.tr.get_short("dt")/1.0e6;
			if (WE.dt<=0) WE.EH.die_err(filename+": dt not set correctly");
			if (WE.assumeNoReflectionsOutsideLogs)
			  {
				// may need safe data margin
				tStartCube = WE.tStart - WE.tCodaMax;
				tEndCube = WE.tEnd - WE.tPrecursorMin;
			  }
			else
			  {
				tStartCube = WE.tStart;
				tEndCube = WE.tEnd;
			  }
			// additional margin for registration errors
			if(systematicTProp.mean != 0.0)
			  {
				tStartCube += systematicTProp.mean;
				tEndCube += systematicTProp.mean;
			  }
			if (systematicTProp.sigma>0)
			  {
				tStartCube -= SYSTEMATIC_T_SHIFT_SAFETY_MARGIN*systematicTProp.sigma;
				tEndCube += SYSTEMATIC_T_SHIFT_SAFETY_MARGIN*systematicTProp.sigma;
			  }
			NzLast = (int)Math.ceil((tEndCube - f1)/WE.dt);
			if (NzLast<0) WE.EH.die_err(filename+": time range of data does not extend to time "+tEndCube);
			NzFirst = (int)Math.floor((tStartCube - f1)/WE.dt);
			if (NzFirst<0) WE.EH.die_err(filename+": time range of data does not extend to time "+tStartCube);
			tStartData = f1+NzFirst*WE.dt;
			X0 = get_value(WE.tr,WE.seismicXword)*scalco;
			Y0 = get_value(WE.tr,WE.seismicYword)*scalco; 
			inline_first = get_value(WE.tr,WE.seismicInlineWord);
			xline_first = get_value(WE.tr,WE.seismicXlineWord);
		  }
		// check ns is sane
		ns_this = (int)WE.tr.get_short("ns");
		if (ns==0) ns=ns_this;
		if (ns_this!=ns) WE.EH.die_err(filename+": ns is not regularly sampled");
 
		// get coords
		inline = get_value(WE.tr,WE.seismicInlineWord);
		xline = get_value(WE.tr,WE.seismicXlineWord);
		inline_min = Math.min(inline_min,inline);
		xline_min = Math.min(xline_min,xline);
		inline_max = Math.max(inline_max,inline);
		xline_max = Math.max(xline_max,xline);
		X = get_value(WE.tr,WE.seismicXword)*scalco;
		Y = get_value(WE.tr,WE.seismicYword)*scalco; 
		// Store X,Y's
		XYValues.add(new SeismicXYValues(X,Y));

		// System.out.println("POS: "+ inline+" "+xline+" "+X+" "+Y);

		// successive inline,xline increments must be "regular". Only permissible increments
		// are 0, dInline, Dinline
		if (trace_number>0)
		  {
			if ((deltaInline = inline-inlineOld) !=0)
			  {
				if (trace_number==1) fastDirection = 'i';
				if (deltaInline!=0)
				  {
					if (dInline==0) dInline = deltaInline;                  // first increment
					if (deltaInline != dInline && Dinline==0) Dinline = deltaInline;  // first end-of line back sweep
					// henceforth, must have  deltaInline=dInline || deltaInline=Dinline
					if (deltaInline != dInline && deltaInline != Dinline)
					  {
						eps = Math.abs((deltaInline-dInline)/dInline);
						if (eps > TOLERABLE_SPACING_IRREGULARITY_EPSILON)
						  {
							WE.EH.die_err("Seismic minicube " + filename +" is not regularly spaced in inline direction. "+
										  "To see successive differences, try:%n sugethw verbose=0 key="+
										  WE.seismicInlineWord+" < "+filename + 
										  " | grep -v '^$' | awk '{gsub(/"+WE.seismicInlineWord+"=/,\"\");print $1-p;p=$1}'");
						  }
					  }
				  }
			  }
			if ((deltaXline = xline-xlineOld) !=0)
			  {
				if (trace_number==1) fastDirection = 'x';
				if (deltaXline!=0)
				  {
					if (dXline==0) dXline = deltaXline;                  // first increment
					if (deltaXline != dXline && Dxline==0) Dxline = deltaXline;  // first end-of line back sweep
					// henceforth, must have deltaXline=dXline || deltaXline=Dxline
					if (deltaXline != dXline && deltaXline != Dxline)
					  {
						WE.EH.die_err("Seismic minicube " + filename +" is not regularly spaced in xline direction. "+
									  "To see successive differences, try:%n sugethw verbose=0 key="+
									  WE.seismicXlineWord+" < "+filename + 
									  " | grep -v '^$' | awk '{gsub(/"+WE.seismicXlineWord+"=/,\"\");print $1-p;p=$1}'");
					  }
				  }
			  }
		  }
		inlineOld=inline;
		xlineOld=xline;
	  }
	inline_last= inline;
	xline_last= xline;
	XYValues.trimToSize();

	if (fastDirection != 'i' && (XYValues.size() > 1))
	  {
		// this should never happen on account of pre-emptive determination of fastest direction above
		WE.EH.die_err("Seismic minicube " + filename +" does not have fast direction corresponding to inline word "
					  + WE.seismicInlineWord +".%nSuggest you swap inline/xline words");
	  }



	// Auto detect funky x/y header words with extra powers of ten if required
	// Not required: scalco should contain all we need
	if (AUTO_DETECT_SEISMIC_XY_WORD_UNITS)
	  {
		double xMid,yMid,test,corr,xMean=0,yMean=0;
		SeismicXYValues pt;
		// Compare typical x'y's to midpoint of longest logged interval (seems a reasonable choice)
		LogData lD = WE.logData[wellNumber];
		// x,y's of midpoint of survey file
		xMid = (lD.survey.X[0]+lD.survey.X[lD.survey.X.length-1])/2.0;
		yMid = (lD.survey.Y[0]+lD.survey.Y[lD.survey.Y.length-1])/2.0;
		for (i=0;i<XYValues.size();i++)
		  {
			pt = (SeismicXYValues)XYValues.get(i);
			xMean += Math.abs(pt.x);
			yMean += Math.abs(pt.y);
		  }
		xMean/=XYValues.size();
		yMean/=XYValues.size();
		if (xMid!=0 && xMean!=0 && yMid!=0 && yMean!=0)
		  {
			// test X
			test = xMean/Math.abs(xMid);
			if (test>AUTO_DETECT_SEISMIC_XY_WORD_THRESHOLD || test<(1.0/AUTO_DETECT_SEISMIC_XY_WORD_THRESHOLD))
			  {
				// and Y
				test = yMean/Math.abs(yMid);
				if (test>AUTO_DETECT_SEISMIC_XY_WORD_THRESHOLD || test<(1.0/AUTO_DETECT_SEISMIC_XY_WORD_THRESHOLD))
				  {
					// looks credible. Check for power of 10 factor
					test = (Math.abs(xMean)/Math.abs(xMid)) / (Math.abs(yMean)/Math.abs(yMid));
					if ((test>(1-AUTO_DETECT_SEISMIC_XY_WORD_POWER_EPS)) 
						&& (test<(1+AUTO_DETECT_SEISMIC_XY_WORD_POWER_EPS)))
					  {
						// same ratios of x and y
						// find this ratio and reset the x/y words
						test = Math.abs(xMean)/Math.abs(xMid);
						int pow = (int)Math.rint(cern.jet.math.Arithmetic.log10(test));
						corr = Math.pow(10,pow);
						WE.EH.println("     WARNING! Correcting seismic x,y header words with factor 10^"+pow+". Check this!");
						for (i=0;i<XYValues.size();i++)
						  {
							pt = (SeismicXYValues)XYValues.get(i);
							pt.x /= corr;
							pt.y /= corr;
						  }
					  }
				  }
			  }
		  }
	  }

	// Checks
	// 1) Does the user specify the dimension correctly?
	if ((inline_max==inline_min) && (xline_max==xline_min))
	  {
		// Data set is single trace
		if (!is1DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is single trace. Switch on \"is-single-trace-only\"" +
									   " XML attribute if trace is OK");	
		  }
		if (is2DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a trace. Remove two-D-utm-conversions entry for 2D line");	
		  }
	  }
	else if ((inline_max==inline_min) || (xline_max==xline_min))
	  {
		// 2D
		if (is1DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 2D line . Switch off \"is-single-trace-only\" XML attribute entry");
		  }
		if (!is2DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 2D line. Add in \"two-D-utm-conversions\" entry for 2D line");	
		  }
	  }
	else
	  {
		if (is1DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 3D cube. Switch off \"is-single-trace-only\" XML entry");
		  }
		if (is2DdataAndUseSnapping)
		  {
			throw new RuntimeException(filename + 
									   " data set is a 3D cube. Remove two-D-utm-conversions entry for 2D line");	
		  }
	  }
	
	// Is 2D utm spec OK?
	if (is2DdataAndUseSnapping && is2DdataAndUseXYwords)
	  {
		utmInfo = new UTMConversionInfo();

		utmInfo.x[0] = ((SeismicXYValues)XYValues.get(0)).x;
		utmInfo.y[0] = ((SeismicXYValues)XYValues.get(0)).y;
		utmInfo.inLine[0] = inline_first;
		utmInfo.xLine[0] = xline_first;
		
		utmInfo.x[1] = ((SeismicXYValues)XYValues.get(XYValues.size()-1)).x;
		utmInfo.y[1] = ((SeismicXYValues)XYValues.get(XYValues.size()-1)).y;
		utmInfo.inLine[1] = inline_last;
		utmInfo.xLine[1] = xline_last;

		if (utmInfo.xLine[0] != utmInfo.xLine[1])
		  {
			throw new RuntimeException("2D data set doesn't have a common xLine in: " + filename);
		  }
		// set the third point by trivial 90deg rotation, preserving squareness of geometry
		// Points then end up like this:          2    1
		//                                          0
		utmInfo.x[2] = utmInfo.x[0]+utmInfo.y[0]-utmInfo.y[1];
		utmInfo.y[2] = utmInfo.y[0]+utmInfo.x[1]-utmInfo.x[0];
		utmInfo.inLine[2] = utmInfo.inLine[0]+utmInfo.xLine[0]-utmInfo.xLine[1];
		utmInfo.xLine[2]  = utmInfo.xLine[0]+utmInfo.inLine[1]-utmInfo.inLine[0];
		
		// setup UTMS and return false if they're bogus
		if (!utmInfo.setUpUTMmapping(WE.EH))
		  {
			throw new RuntimeException("2D data set cannot be used to set geometry property from header words: " + filename);
		  }
	  }



	// allocate required storage for SU minicube
	Ninline = 1 + (inline_max-inline_min)/(dInline!=0 ? Math.abs(dInline) : 1);
	Nxline = 1 + (xline_max-xline_min)/(dXline!=0 ? Math.abs(dXline) : 1);
	Nz = NzLast - NzFirst+1;

	//	data = new DenseDoubleMatrix3D(Ninline,Nxline,Nz);  JGDENSEHACK
	data = new au.csiro.JamesGunning.lib.Numerical.DenseDoubleMatrix3D(Ninline,Nxline,Nz);

	// second pass over data. Read it and stuff into data cube
	// set file ptr to start
	try 
	  {
		in_strm.seek(0);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't rewind file ptr to 0 in " +filename);
	  }
	X0 = X = ((SeismicXYValues)XYValues.get(0)).x; // FIXME REMOVE	
	Y0 = Y = ((SeismicXYValues)XYValues.get(0)).y; // FIXME REMOVE	
	for (trace_number=0;;trace_number++)
	  {
		gotTrace = (WE.sureader.fgettr(in_strm,WE.tr)!=0);
		if (!gotTrace) break;
		// get coords
		inline = get_value(WE.tr,WE.seismicInlineWord);
		xline = get_value(WE.tr,WE.seismicXlineWord);
		// matrix indices
		i = (inline-inline_min)/((dInline!=0) ? Math.abs(dInline) : 1);
		j = (xline-xline_min)/(dXline!=0 ? Math.abs(dXline) : 1);
		if (i<0 || i>=Ninline)  WE.EH.die_err("inline indx out of range");
		if (j<0 || j>=Nxline)  WE.EH.die_err("xline indx out of range");

		// Store X,Y's in corresponding unrolled indexed location of [i][j]
		// but with inline the fastest direction
		X = get_value(WE.tr,WE.seismicXword)*scalco;
		Y = get_value(WE.tr,WE.seismicYword)*scalco; 
		SeismicXYValues	XY = (SeismicXYValues)XYValues.get(j*Ninline+i);
		XY.x = X;
		XY.y = Y;
		// System.out.println("SIXY: "+(X-X0)+" "+(Y-Y0)); // FIXME REMOVE		

		// and samples
		k = (int)Math.rint((tStartData-f1)/WE.dt);
		for (kk=0;kk<Nz;kk++,k++) 
		  {
			data.setQuick(i,j,kk,(double)WE.tr.get_data_value(k));
		  }
	  }
	try 
	  {
		in_strm.close();
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " +filename);
	  }

	// set object fields
	inlineSpacing = Math.abs(dInline);
	xlineSpacing = Math.abs(dXline);

	// Compute average dX, dY values (to prevent flt. pt error)

	SeismicXYValues originPt = (SeismicXYValues)XYValues.get(0);
	SeismicXYValues oppositeCornerPt = (SeismicXYValues)XYValues.get(XYValues.size()-1);
	SeismicXYValues endOfFirstSweepPt=null;
	endOfFirstSweepPt	= (SeismicXYValues)XYValues.get(Ninline-1);

	X0 = X = originPt.x;
	Y0 = Y = originPt.y;

	// length of minicube side along inline direction from Pythagoras...
	OEx = endOfFirstSweepPt.x - X0;
	OEy = endOfFirstSweepPt.y - Y0;

	if (true)
	  {
		// inline spans directly from inline side of minicube. Fast but
		// sensitive to flt pt error of endOfFirstSweepPt point.
		inLineSpan = Math.sqrt(OEx*OEx+OEy*OEy);
		mean_dInline = (Ninline>1) ? inLineSpan/(Ninline-1) : 0.0;
	  }
	else
	  {
		// statistical average method
		mean_dInline=0.0;
		count=0;
		for (i=0,k=0;i<Nxline;i++)
		  {
			for (j=0;j<Ninline;j++,k++)
			  {
				SeismicXYValues nextPt = (SeismicXYValues)XYValues.get(k);
				// System.out.println(i+" "+j+" "+nextPt.x+" "+ nextPt.y);
				if (j>0)
				  {
					// average  distance hops in the inline direction
					tmp = Math.sqrt((nextPt.x-X)*(nextPt.x-X) + (nextPt.y-Y)*(nextPt.y-Y));
					mean_dInline += tmp;
					count++;
					// System.out.println("POS: "+ nextPt.x+" "+ nextPt.y+" " + tmp);
				  }
				X=nextPt.x;
				Y=nextPt.y;
			  }
		  }
		if (count>0)
		  {
			mean_dInline /=(double)count;
		  }
		else
		  {
			mean_dInline = 0.0;
		  }
		inLineSpan = mean_dInline*(Ninline-1);
	  }


	
	ECx = oppositeCornerPt.x - endOfFirstSweepPt.x;
	ECy = oppositeCornerPt.y - endOfFirstSweepPt.y;

	// Compute mean_dXline from Pythagoras over full diagonal span of minicube
	Dx_c = oppositeCornerPt.x - X0;
	Dy_c = oppositeCornerPt.y - Y0;
	diagonalSpan = Math.sqrt(Dx_c*Dx_c+Dy_c*Dy_c);
	xLineSpan = Math.sqrt(ECx*ECx+ECy*ECy);
	cubeXYangle = Math.toDegrees(Math.acos((inLineSpan*xLineSpan>0) ? 
										 (inLineSpan*inLineSpan+xLineSpan*xLineSpan-diagonalSpan*diagonalSpan)/
										 (2.0*inLineSpan*xLineSpan) : 0));
	mean_dXline = (Nxline>1) ? xLineSpan/(Nxline-1) : 0.0;
	
	if (WE.verbosity>2)
	  {
		WE.EH.println("        Minicube: "+filename +", Nxline="+Nxline+" Ninline="+Ninline);
		WE.EH.println("                  trace spacing inline: %-2g crossline: %-2g",
					  mean_dInline, mean_dXline);
		WE.EH.println("                  minicube x/y axis angle (degrees): %-3g",cubeXYangle);
	  }

	// Now set up rotation geometry (which may be trivial)
	// use three points at (inline,xline) = 
	// (inline_min,xline_min),(inline_max,xline_min),(inline_max,xline_max)
	// so form linear system for
	// (inLine-inline_min, xline-xline_min) = A.(X-X0,Y-Y0)

	ECx = endOfFirstSweepPt.x - X0;
	ECy = endOfFirstSweepPt.y - Y0;

	// protection not really required, but be safe.
	double tmp1,tmp2;
	tmp1 = Dy_c*ECx - Dx_c*ECy;
	tmp2 = -(Dy_c*ECx) + Dx_c*ECy;

	// entries in A matrix
	A11 = (tmp1!=0) ? (((Dy_c - ECy)*(inline_max - inline_min))/tmp1) : 0;
	A12 = (tmp2!=0) ? ((Dx_c - ECx)*(inline_max - inline_min))/tmp2: 0;
	A21 = (tmp2!=0) ? ((ECy*(xline_max - xline_min))/tmp2) : 0;
	A22 = (tmp1!=0) ? (ECx*(xline_max - xline_min))/tmp1 : 0;


	// Silly checks
	if (!true)
	  {
		System.out.println("Seismic transformation checks....");
		computeIndices(oppositeCornerPt.x,oppositeCornerPt.y);
		System.out.println("Corner pt inline xline = "+inLine+" "+xLine);
		computeIndices(X0,Y0);
		System.out.println("First trace inline xline = "+inLine+" "+xLine);

		SeismicXYValues pt = (SeismicXYValues)XYValues.get(1);
		computeIndices(pt.x,pt.y);
		System.out.println("Second trace inline xline = "+inLine+" "+xLine);

		pt = (SeismicXYValues)XYValues.get(Ninline-1);
		computeIndices(pt.x,pt.y);
		System.out.println(Ninline +"-th trace inline xline = "+inLine+" "+xLine);

		pt = (SeismicXYValues)XYValues.get(XYValues.size()-Ninline);
		computeIndices(pt.x,pt.y);
		System.out.println((XYValues.size()-Ninline+1) +"-th trace inline xline = "+inLine+" "+xLine);

		System.exit(0);
	  }
  }

  /**
	 After loadSeismicData() is called, this method computes the effective
	 (double) indices of the point x,y, which can then be retrieved with
	 getInLine() getXline()
  */
  void computeIndices(double x,double y)
  {
	double dX,dY;

	dX = x - X0;
	dY = y - Y0;
	
	inLine = inline_min + A11*dX + A12*dY;
	xLine = xline_min + A21*dX + A22*dY;
	// System.out.println("inline/xline = " +inLine+ " "+xLine);
  }

  /**
	 Returns effective inLine index of a point x,y. Sensible only after
	 computeIndices(x,y) is called.
  */
  double getInLine()
  {
	//System.out.println("inline = " +inLine);
	return inLine;
  }
  /**
	 Returns effective xLine index of a point x,y. Sensible only after
	 computeIndices(x,y) is called.
  */
  double getXLine()
  {
	//System.out.println("xline = " +xLine);
	return xLine;
  }

  private int get_value(SuSegyFormat tr, String keyword)
  {
	int value=0;

	SuSegyFormat.SuHeaderEntry entry = (SuSegyFormat.SuHeaderEntry) tr.headerMap.get(keyword);

	switch (entry.type) 
	  {
	  case SuSegyFormat.SuHeadTypes.short_type:
		value = (int)tr.get_short(keyword);
		break;
	  case SuSegyFormat.SuHeadTypes.int_type:
		value = tr.get_int(keyword);
		break;
	  }
	return value;
  }
  
  public class SeismicXYValues
  {
	double x;
	double y;

	public SeismicXYValues(double X,double Y)
	{
	  this.x = X;
	  this.y = Y;
	}
  }
}
