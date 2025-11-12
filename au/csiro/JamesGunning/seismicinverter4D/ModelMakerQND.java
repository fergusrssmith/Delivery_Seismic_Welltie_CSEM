package au.csiro.JamesGunning.seismicinverter4D;

import java.lang.*;
import java.text.*;
import java.io.*;
import java.util.*;
import static java.lang.Math.*;


import cern.colt.list.*; 

import au.csiro.JamesGunning.lib.*;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;


public final class ModelMakerQND
{
  public static final double PHI_LWR_BND=0.25; // Lower bound on porosity

  boolean noPad=false;
  boolean useInsituDepth=false;
  boolean mergeMode=false;
  String mergeFiles=null;
  Model m;
  ErrorHandler EH;
  String filename;
  String oFile;
  int dumpDZlayer=-1;
  double T_TOP_DEFAULT = 1500;
  float SgasThresh=0;
  int chunk=0;

  ModelMakerQND()
  {
	EH = new ErrorHandler();
  }


  /** ModelMakerQND
	  Sample driver from input streams etc
  */
  public static void main(String argv[])
  {
	ModelMakerQND MM = new ModelMakerQND();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("ModelMakerQND",MM.EH);
	//options.addLicense(DeliveryLicense.license_text);
	options.addVersion(VersionNumber.versionNumber); 

	options.addSpielText(
"Program to assemble model.su file suitable for input to delivery, using ascii file(s) to populate the model.su file. %n%n"+
"Supplied model-file.txt has ascii columns for model properties. It can have "+ 
"# marked header & comments. It is constructed to refer to Nl layers. The model.su file "+
"written out by this program is padded above and below by 1 layer, unless --NP is supplied. "+
"Each ascii line corresponds to a trace location. The columns must be supplied in 3 groups: (1) "+
"the first group of columns is fixed \"info\" for the trace & model, using some selection from: %n%n"+
    "X_grid or X            : typical grid or cmp X location %n"+
    "Y_grid or Y            : typical grid or cmp Y location %n"+
    "X_map                  : typical grid or cmp utm alternative X location %n"+
    "Y_map                  : typical grid or cmp utm alternative Y location %n"+
    "Ecl_i                  : inline or eclipse i grid index %n"+
    "Ecl_j                  : xline or eclipse j grid index %n"+
    "TWT                    : Two way time of top layer%n%n"+
"Then (2) there are a sequence of Nl column-sets of block properties, formed by some selection from: %n"+
	"Z     : depth of layer %n"+
	"t     : time of layer %n"+
	"phi   : porosity %n"+
	"pr_i  : pressure %n"+
	"NG    : Net to gross %n"+
	"S     : saturation %n"+
	"vpg   : special(gas) velocity %n"+
	"rho   : density  %n"+
	"rock  : rock-type %n"+
"(Typically the entries are like NG_1 Z_1 t_1 NG_2 Z_2 t_2 .....) %n%n"+
"Finally (3) the last column must be one of %n"+
"zbot                   : bottom depth of lowest layer %n"+
"t_base                 : bottom time of lowest layer %n%n"+
"A mini-example file might look like %n%n"+
"#My minimodel %n"+
"Ecl_i %n"+
"Ecl_j %n"+
"NG %n"+
"t %n"+
"t_base %n"+ 
"1 1 .12 900 .23 920 .31 925 .42 930 950 %n"+
"1 2 .15 900 .25 920 .35 925 .41 931 955");



   	options.addMandatoryArg(1,"model-file.txt","ascii file with properties in columns required to form model.su file");
	options.addOption("-b",null,1,"int","delivery_layer_num","dump block thicknesses to stdout, named delivery_layer");
	options.addOption("--TWT-top",null,1,"float","TWT","specify (constant) top-layer two-way time (ms) (default "+MM.T_TOP_DEFAULT+")");
	options.addOption("--min-Sgas",null,1,"float","TWT","specify S_gas threshold below which to zero S_gas values (default 0)");
	options.addOption("-o",null,1,"String","filename.su","output .su filename (default model.su)");
	options.addOption("--NP",null,0,null,null,"Don't add padding layer above & below");
	options.addOption("--Zinsitu",null,0,null,null,"Use (don't compute) insitu depth Z if supplied");
	options.addOption("--merge",null,0,null,null,
					  "Don't do any of the usual. Take the main argument as "+
					  "file1.txt,file2,txt...(split on commas). "+
					  "Just merge the simple ascii numbers in %nfile1,file2 etc "+
					  "into a big model.su file. One layer-property per file, in the order supplied. All files must "+
					  "have the same number count");
	options.addOption(null,"--chunk",1,"int","chunk-size","set mark header word to chunkify");

	options.parseOptions(argv,true);

	MM.filename = options.getMandatoryArg(1);

	//  stuff commandline options into fields
	if (options.isOptionSet("-b")) 
	  {
		MM.dumpDZlayer = ((Integer)options.getOptionValues("-b")[0]).intValue();
		if (MM.dumpDZlayer<1) MM.EH.die_err("Dump -b layer must be > 0"); 
	  }
	if (options.isOptionSet("--chunk")) 
	  {
		MM.chunk = options.getOptionInt("--chunk",0);
		if (MM.chunk<1) MM.EH.die_err("chunk must be > 0"); 
	  }
	if (options.isOptionSet("--TWT-top")) 
	  {
		MM.T_TOP_DEFAULT = ((Float)options.getOptionValues("--TWT-top")[0]).floatValue();
	  }
	if (options.isOptionSet("--min-Sgas")) 
	  {
		MM.SgasThresh = ((Float)options.getOptionValues("--min-Sgas")[0]).floatValue();
	  }
	if (options.isOptionSet("-o")) 
	  {
		MM.oFile = ((String)options.getOptionValues("-o")[0]);
	  }
	if (options.isOptionSet("--NP")) MM.noPad=true; 
	if (options.isOptionSet("--Zinsitu")) MM.useInsituDepth=true; 
	if (options.isOptionSet("--merge")) 
	  {
		MM.mergeFiles = MM.filename;
		MM.mergeMode=true;
	  }

	
	MM.run();
  }

  void run()
  {
	if (!mergeMode)
	  {
		runStandard();
	  }
	else
	  {
		runMerge();
	  }
  }

  void runStandard()
  {
	int i,j,k,Ncols,Nlayer_cols,Nprops,i_ecl=1,j_ecl=1,m,mD;
	int Necl_layers=0,Ndelivery_layers=0,Nlayers_above,Nlayers_below,Necl_props=0,Ndelivery_props=0;
	int m_phi=-1,m_pri=-1,m_NG=-1,m_Z=-1,m_T=-1,j_TWT=-1,j_X_grid=-1,j_X_map=-1,j_Y_grid=-1,j_Y_map=-1;
	int j_Ecl_i=-1,j_Ecl_j=-1,j_zbot=-1,j_t_base=-1;
	int m_S=-1,m_vpg=-1,m_rhog=-1,m_rocktype=-1,m_dfc=-1;
	int active_trace_cntr;
	double x=0,y=0,z,x_utm=0,y_utm=0,val,z_last=0,rho=1,dt,t_top,t=0,T_TOP=0,dz;
	double VP_APPROX, T_PAD;
	boolean write_Z=false, write_phi=false, write_pr=false, write_NG=false,write_Sg=false;
	boolean write_vpg=false, write_rhog=false, write_rock=false,write_dfc=false,write_t=false,have_t=false;
	ArrayList<String> outP=null;
	boolean parsedHeader=false, scalcoSet=false,getTfromDepthAndVelocities=false;
	String s,sc;
	String[] st;
	SuSegyFormat tr;
	SuWriter suwriter;
	FileOutputStream out_strm;
	String fileoutname=(oFile==null) ? "model.su" : oFile;
	long scalco=0;
	double scaleXY=1;

	if (noPad)
	  {
		Nlayers_above = Nlayers_below = 0;
	  }
	else
	  {
		Nlayers_above = Nlayers_below = 1;
	  }
	VP_APPROX=3000; // guess of velocities
	T_PAD=20.0;
	BufferedReader in=null;

	try
	  {
		if (filename.equals("-"))
		  {
			// use stdin
			in = new BufferedReader(new FileReader(FileDescriptor.in));
		  }
		else
		  {
			in = new BufferedReader(new FileReader(filename));
		  }
	  }
	catch (Exception e)
	  {
		EH.die_err("Can't find or use file "+filename);
	  }
		
	tr = new SuSegyFormat();
	suwriter = new SuWriter(); // Open stream.
	try 
	  {
		out_strm = (FileOutputStream) new FileOutputStream(fileoutname);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open "+fileoutname);
	  }

	k=1;
	Nprops=0;
	Ncols=0;
	Necl_props=0;
	mD=0;
	active_trace_cntr=0;
	outP=new ArrayList<String>();
	try
	  {
		while ((s = in.readLine()) != null)
		  {
			s = s.trim();
			if (s.startsWith("#")) continue;
			
			// remove {}
			sc = s.replaceAll("[{}]",""); 
			st = sc.split("\\s+");

			// parse header
			if (st.length==1 && !parsedHeader)
			  {
				sc = st[0];
				if (sc.startsWith("X") || sc.startsWith("Y") || sc.startsWith("Ecl") || 
					sc.startsWith("TWT") || sc.startsWith("zbot")|| sc.startsWith("t_base")) 
				  {
					// 1-Per-trace properties
					if (sc.startsWith("TWT")) j_TWT = Ncols;
					else if (sc.startsWith("X_grid")||sc.equals("X")) j_X_grid = Ncols;
					else if (sc.startsWith("Y_grid")||sc.equals("Y")) j_Y_grid = Ncols;
					else if (sc.startsWith("Ecl_i")) j_Ecl_i = Ncols;
					else if (sc.startsWith("Ecl_j")) j_Ecl_j = Ncols;
					else if (sc.startsWith("zbot")) j_zbot = Ncols;
					else if (sc.startsWith("t_base")) j_t_base = Ncols;
					else if (sc.startsWith("X_map")) j_X_map = Ncols;
					else if (sc.startsWith("Y_map")) j_Y_map = Ncols;
					mD=0;
				  }
				else 
				  {
					// 1-per gridblock properties
					if      (sc.startsWith("Z")) {m_Z = mD; Necl_props++; write_Z=true; outP.add("depth");}
					else if (sc.startsWith("t")) {m_T = mD; Necl_props++;have_t=true;} // we'll add t at end
					else if (sc.startsWith("phi")) {m_phi = mD; Necl_props++;write_phi=true; outP.add("phi");}
					else if (sc.startsWith("pr_i")) {m_pri = mD; Necl_props++;write_pr=true; outP.add("pressure");}
					else if (sc.startsWith("NG")) {m_NG = mD; Necl_props++;write_NG=true; outP.add("NG");}
					else if (sc.startsWith("S")) {m_S = mD; Necl_props++;write_Sg=true; outP.add("Sg");}
					else if (sc.startsWith("vpg")) {m_vpg = mD; Necl_props++;write_vpg=true; outP.add("vpg");}
					else if (sc.startsWith("rho")) {m_rhog = mD; Necl_props++;write_rhog=true; outP.add("rhog");}
					else if (sc.startsWith("rock")) {m_rocktype = mD; Necl_props++;write_rock=true; outP.add("rock_type");}	
					else if (sc.startsWith("depth_fc")) {m_dfc = mD; Necl_props++;write_dfc=true; outP.add("depth_fc");}	
					mD++;
				  }
				Nprops++;
				Ncols++;
			  }
			else if (!parsedHeader)
			  {
				Ncols = st.length;
				Nlayer_cols = Ncols;
				// Determine number of eclipse/model layers
				// probably lose 3 from i,j, zbot (last col)
				if (j_Ecl_i >= 0) {Nprops--;Nlayer_cols--;}
				if (j_Ecl_j >= 0) {Nprops--;Nlayer_cols--;}
				if (j_zbot >= 0) {Nprops--;Nlayer_cols--;}
				if (j_t_base >= 0) {Nprops--;Nlayer_cols--;}
				if (j_X_grid >= 0) {Nprops--; Nlayer_cols--;}
				if (j_Y_grid >= 0) {Nprops--; Nlayer_cols--;}
				if (j_X_map >= 0) {Nprops--; Nlayer_cols--;} 
				if (j_Y_map >= 0) {Nprops--; Nlayer_cols--;} 
				if (j_TWT >= 0) {Nprops--; Nlayer_cols--;}
				Necl_layers = Nlayer_cols/Nprops;
				// checks
				// here, parse header X Y ecl_i, ecl_j /tuples of z, phi, Sgas, vp, rho rock/ z_bottom
				// Necl_props = 6; // infer from header! z, phi, Sgas, vp, rho
				Ndelivery_props = Necl_props;   // z,phi,Sgas,vp_gas,rho_gas,rock..... 
				if (!have_t) Ndelivery_props++; //  then t{t_base}
				if (!have_t && write_Z) Ndelivery_props++;   // add thickness too

				//if (Nprops<6 || Nprops>7) EH.die_err("Nprops seems wrong. Expect "+Necl_props);

				Necl_layers = Ncols;    // Don't count X*,Y*, ecl_i,ecl_j, Zlast
				if (j_Ecl_i >= 0) Necl_layers--;
				if (j_Ecl_j >= 0) Necl_layers--;
				if (j_zbot >= 0) Necl_layers--;
				if (j_t_base >= 0) Necl_layers--;
				if (j_X_grid >= 0) Necl_layers--;
				if (j_X_map >= 0) Necl_layers--;
				if (j_Y_grid >= 0) Necl_layers--;
				if (j_Y_map >= 0) Necl_layers--;
				if (j_TWT >= 0) Necl_layers--;

				Necl_layers /= Necl_props;

				Ndelivery_layers = Necl_layers + Nlayers_below + Nlayers_above;

				EH.warn("Inferred %d eclipse properties, %d eclipse layers, => %d delivery layers%n",
						Necl_props,Necl_layers,Ndelivery_layers);
				if (!have_t) // not present on input
				  {
					if (write_Z) 
					  {
						// input model in depth. Get t from z, velocities
						getTfromDepthAndVelocities = true;
						write_t=true; // given depth, we'll attempt t :-)
						// last 2 props
						outP.add("thickness");	
						outP.add("t");	
					  }
					else EH.die_err("Can't form time samples: neither t not Z available!");
				  }
				else
				  {
					write_t=true; 
					outP.add("t");	
				  }
				tr.set("ns",(short)(Ndelivery_layers*Ndelivery_props + 1));
				tr.set("dt",(short)4000); // token
				if (j_TWT==-1)
				  {
					T_TOP=T_TOP_DEFAULT; // FIXME top time surface 
				  }
				parsedHeader = true;
			  }
			if (!parsedHeader) continue;

			boolean isNotActive=true;
			// read & fill values
			try
			  {
				i=0;
				if (j_X_grid  >=  0) {x = Double.parseDouble(st[j_X_grid]); i++;}
				if (j_Y_grid >=  0) {y = Double.parseDouble(st[j_Y_grid]); i++;}
				if (j_X_map >=  0) {x_utm = Double.parseDouble(st[j_X_map]); i++;}
				if (j_Y_map >=  0) {y_utm = Double.parseDouble(st[j_Y_map]); i++;}
				if (j_Ecl_i >= 0) {i_ecl = Integer.parseInt(st[j_Ecl_i]); i++;}
				if (j_Ecl_j >= 0) {j_ecl = Integer.parseInt(st[j_Ecl_j]); i++;}

				if (j_TWT==-1) // time not available: use default 
				  {
					// top depth
					t_top = T_TOP;
					t = t_top;
				  }
				else if (i==j_TWT)
				  {
					t_top = Double.parseDouble(st[j_TWT]); i++;
					t = t_top;
				  }
				
				for(j=0;i<Ncols;i++,j++)
				  {
					val = abs(Double.parseDouble(st[i]));
					// Delivery (zero-offset) layer no.
					k = j/Necl_props + Nlayers_above;
					// property indx
					m = j%Necl_props;
					// 0,1,2,3... = z, phi, Sgas, vp, rho
					if (m==m_phi) val = max(val,PHI_LWR_BND); // lower trun on phi FIXME
					else if (m==m_rhog) val*=0.001; // density to gm/m^3
					else if (m==m_S)
					  {
						if (val<SgasThresh) val=0.0;
					  }
					else if (m==m_T)
					  {
						// t available & input
						t = val;
						if (t!=0) 
						  {
							isNotActive = false;
						  }
					  }
					// delivery property index
					// z, phi, Sgas, vp, rho....
					mD = m*Ndelivery_layers + k;
					// special exception for t_base
					if (j_t_base>=0 && (i==(Ncols-1))) mD = Ndelivery_props*Ndelivery_layers;
					tr.set_data_value((float)val,mD);
					if (m==m_Z && !useInsituDepth) 
					  {
						// depth
						z = val;
						if (j>0)
						  {
							// not-the-top layer
							// add approximate two-way time (ms)
							dz = (z-z_last);
							dt = 2000*dz/VP_APPROX;
							t += dt;

							// write thickness of layer above
							mD = (Ndelivery_props-2)*Ndelivery_layers + (k-1);
							tr.set_data_value((float)dz,mD);
							if ((k-1)==dumpDZlayer) EH.printf("%d %d %g%n",i_ecl,j_ecl,dz);
						  }
						z_last = z;
						if (getTfromDepthAndVelocities)
						  {
							// write computed times into time block if needed
							mD = (Ndelivery_props-1)*Ndelivery_layers + k;
							tr.set_data_value((float)t,mD);
							if (j==0)
							  {
								// padding at start
								tr.set_data_value((float)(t-T_PAD),mD-1);
							  }
							if (i==Ncols-1) 
							  {
								// end time
								tr.set_data_value((float)(t+T_PAD),mD+1);
							  }
						  }
					  }
					//showTraceState(tr);					  
				  }

				if (getTfromDepthAndVelocities)
				  { 
					mD = (Ndelivery_props-1)*Ndelivery_layers + k;
					tr.set_data_value((float)t,mD);
				  }
				//showTraceState(tr);
	
				tr.set("ep",(int)i_ecl);
				tr.set("cdp",(int)j_ecl);
				if (!scalcoSet)
				  {
					double mx,pwr=0;
					mx = max(abs(x),abs(y));
					if (j_X_map>0)
					  {
						mx = max(max(abs(x_utm),abs(y_utm)),mx);
					  }
					if (mx!=0)
					  {
						scalco = 1+(int)floor(log10(mx/Integer.MAX_VALUE));
					  } 
					else scalco = 0; // guess
					scalco = min(scalco,4);
					scalco = max(scalco,-4);
					scaleXY = pow(10,scalco);
					scalcoSet = true;
				  }
				x/=scaleXY;
				y/=scaleXY;
				tr.set("scalco",(short)((scaleXY>1) ? scaleXY : -(1.0/scaleXY)));
				tr.set("gx",(int)round(x));
				tr.set("gy",(int)round(y));
				if (j_X_map>0)
				  {
					tr.set("sx",(int)round(x_utm/scaleXY));
					tr.set("sy",(int)round(y_utm/scaleXY));
 				  }
				if(chunk>0)
				  {
					tr.set("mark",(short)(active_trace_cntr/chunk));
				  }
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("Line containing "+s+" does not have numbers??");
			  }

			suwriter.fputtr(out_strm, tr);
			if (!isNotActive) 
			  {
				active_trace_cntr++;
			  }
		  }
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("Reading exception "+e);
	  }

	try
	  {
		suwriter.flush(out_strm);
		out_strm.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " + fileoutname);
	  }
	try
	  {
		in.close();
	  }
	catch (Exception e)
	  {
		EH.die_err("Couldn't close "+filename);
	  }
	 	
	EH.println("Warning: set scalco at "+scalco);
	EH.println("Number of active traces detected: "+active_trace_cntr);
	EH.println("Wrote "+fileoutname);
	EH.println("Written properties, in order:");
	EH.printf("%-12s %-15s %-11s%n","block_number","name","num_samples");

	j = outP.size();
	for(i=0;i<j;i++) 
	  {
		s=outP.get(i);
		EH.printf("%-12d %-15s %-11d%n",i+1,s,(i==(j-1)) ? (Ndelivery_layers+1) : Ndelivery_layers);
	  }
	EH.printf("Total samples = %d%n",j*Ndelivery_layers+1); 
  }

  void runMerge()
  {
	SuSegyFormat tr;
	SuWriter suwriter;
	FileOutputStream out_strm;
	String fileoutname=(oFile==null) ? "model.su" : oFile;
	BufferedReader in;
	String s=null;
	String[] st,sc;
	int i,j,N,nFiles;
	double x;
	String[] fArr = mergeFiles.split(",");
	nFiles = fArr.length;
	DoubleArrayList[] vals = new DoubleArrayList[nFiles];
	DoubleArrayList _val;
	for(i=0;i<nFiles;i++)
	  {
		try
		  {
			in = new BufferedReader(new FileReader(fArr[i]));
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't open "+fArr[i]);
		  }
		_val = vals[i] = new DoubleArrayList();
		try
		  {
			while ((s = in.readLine()) != null)
			  {
				s = s.trim();
				if (s.startsWith("#")) continue;
				// split on white space
				st = s.split("\\s+");
				for (j=0;j<st.length;j++)
				  {
					x = Double.parseDouble(st[j]);
					_val.add(x);
				  }
			  }
			in.close();
		  }
		catch (NumberFormatException e)
		  {
			EH.die_err("File "+fArr[i]+": Line containing "+s+" does not have numbers??");
		  }
		catch (IOException e) 
		  {
			throw new RuntimeException("couldn't close " + fArr[i]);
		  }
		if ((i>0) && (_val.size() != vals[0].size()))
		  {
			throw new RuntimeException("File " + fArr[i]+" has different total number-count to "+fArr[0]);
		  }
	  }
	tr = new SuSegyFormat();
	suwriter = new SuWriter();
	try 
	  {
		out_strm = (FileOutputStream) new FileOutputStream(fileoutname);
	  } 
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open "+fileoutname);
	  }
	N = vals[0].size();
	tr.set("ns",(short)nFiles);
	tr.set("dt",(short)4000); // token
	tr.set("cdp",0);
	for(j=0;j<N;j++)
	  {
		tr.set("ep",(int)(j+1));
		for(i=0;i<nFiles;i++)
		  {
			tr.set_data_value((float)(vals[i].get(j)), i);
		  }
		suwriter.fputtr(out_strm, tr);
	  }
	try
	  {
		suwriter.flush(out_strm);
		out_strm.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close " + fileoutname);
	  }
	EH.println("Wrote "+fileoutname);
  }


  // debugging method
  void showTraceState(SuSegyFormat tr)
  {
	int i,ns;
	ns = tr.get_short("ns");
	float values[] = tr.get_data_array();
	for(i=0;i<ns;i++)
	  {
		EH.printf("%-8g ",values[i]);
	  }
	EH.printf("%n");
  }
}