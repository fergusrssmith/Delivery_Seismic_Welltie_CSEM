package au.csiro.JamesGunning.seismicinverter;

import java.lang.*;
import java.text.*;
import java.io.*;
import java.util.*;

import au.csiro.JamesGunning.lib.*;
import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.util.*;


public final class ModelGluer
{
  Model m;
  ErrorHandler EH;
  int Nl;
  String filename1,filename2;
  String[] blocks = new String[2];

  ModelGluer()
  {
	EH = new ErrorHandler();
  }


  /** 
	  ModelGluer
	  Code to merge chunks of 2 delivery model.su type files
  */
  public static void main(String argv[])
  {
	ModelGluer MG = new ModelGluer();

	// Parse basic requirements and extra options
	OptionsParser options = new OptionsParser("ModelGluer",MG.EH);


   	options.addMandatoryArg(1,"N","Number of layers");
   	options.addMandatoryArg(2,"model1.su","Big endian delivery model file 1");
   	options.addMandatoryArg(3,"block-numbers1","blocks from model file 1 (e.g. 1,2)");
   	options.addMandatoryArg(4,"model2.su","Big endian delivery model file 2");
   	options.addMandatoryArg(5,"block-numbers2","blocks from model file 2 (e.g. 1)");



	options.parseOptions(argv,true);


	// stuff commandline options into Seismicinverter fields
	MG.Nl = Integer.parseInt(options.getMandatoryArg(1));
	MG.filename1 = options.getMandatoryArg(2);
	MG.blocks[0] = options.getMandatoryArg(3);
	MG.filename2 = options.getMandatoryArg(4);
	MG.blocks[1] = options.getMandatoryArg(5);
	
	MG.run();
  }

  void run()
  {
	int i,j,k,nTraces;
	int blocksize,Nprops1,Nprops2;
	short ns1,ns2;
	String s;
	String[] sa;
	SuSegyFormat tr1,tr2,trout;
	SuWriter suwriter;
	SuReader sureader;
	FileOutputStream out_strm;
	FileInputStream in_strm1,in_strm2;
	BufferedReader in;
	boolean gotTrace;

	if (Nl<1) EH.die_err("Number of layers must be >=1");
	
	tr1 = new SuSegyFormat();
	tr2 = new SuSegyFormat();
	trout = new SuSegyFormat();
	sureader = new SuReader();


	try
	  {
		in_strm1 = (FileInputStream) new FileInputStream(filename1);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open "+filename1);
	  }


	try
	  {
		in_strm2 = (FileInputStream) new FileInputStream(filename2);
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't open "+filename2);
	  }

	// layer numbers
	int[][] B = new int[2][];
	for(i=0;i<2;i++)
	  {
		sa = blocks[i].split(",");
		if (sa.length<1) EH.die_err("No blocks specified for model file 1");
	
		B[i] = new int[sa.length];
		for (j=0;j<sa.length;j++)
		  {
			try 
			  {
				B[i][j] = Integer.parseInt(sa[j]);
				if (B[i][j]<1) EH.die_err("Block numbers must be 1,2,3...");
			  }
			catch (NumberFormatException e)
			  {
				EH.die_err("block numbers must be numbers only");
			  }
		  }
	  }


	suwriter = new SuWriter(); // Open out stream.
	out_strm = (FileOutputStream) new FileOutputStream(FileDescriptor.out);

	int Nprops_out = B[0].length + B[1].length;
	short ns_out = (short)(Nl*Nprops_out+1);
	trout.set("ns",ns_out);
	trout.set("dt",(short)4000); // for surange???

	nTraces=0;
	for(;;)
	  {
		gotTrace = (sureader.fgettr(in_strm1,tr1)!=0);
		if (gotTrace)
		  {
			ns1 = tr1.get_short("ns");
			if (ns1<1) EH.die_err("ns <1 in "+filename1);
			if (sureader.fgettr(in_strm2,tr2)==0) EH.die_err("Coudn't get (matching) trace from "+filename2);
			else
			  {
				ns2 = tr2.get_short("ns");
				if (ns2<1) EH.die_err("ns < 1 in "+filename2);
				Nprops1 = (ns1-1)/Nl;
				Nprops2 = (ns2-1)/Nl;

				// inherit x,y s etc
				trout.copyHeader(tr1);
				trout.set("ns",ns_out);
				trout.set("dt",(short)4000); // for surange???

				// insert named blocks from tr1 successively into output trace
				float values[] = tr1.get_data_array();
				for (i=0,k=0;i<B[0].length;i++)
				  {
					blocksize = Nl;
					if (B[0][i] == Nprops1)
					  {
						// must be asking for a time from file 1. Not possble: will overwrite properties to come
						EH.die_err("Yuo're not allowed to  use last property (times) in "+filename1);
					  }
					int offset = Nl*(B[0][i]-1);
					for (j=0;j<blocksize;j++)
					  {
						trout.set_data_value(values[j+offset],k++);
					  }
				  }
				// then concatenate ones from tr2
				values = tr2.get_data_array();
				for (i=0;i<B[1].length;i++)
				  {
					blocksize = Nl;
					if (B[1][i] == Nprops2)
					  {
						// must be asking for a time from file 2. Be sure to add extra sample
						blocksize++;
					  }
					int offset = Nl*(B[1][i]-1);
					for (j=0;j<blocksize;j++)
					  {
						trout.set_data_value(values[j+offset],k++);
					  }
				  }
				// write out
				suwriter.fputtr(out_strm, trout);
			  }
			if (!gotTrace) EH.die_err("Coudn't get trace from "+filename2);
		  }
		else
		  {
			if (nTraces==0) EH.die_err("Coudn't get trace from "+filename1);
			else break;
		  }
		nTraces++;
	  }
	
	try
	  {
		in_strm1.close();
		in_strm2.close();
		suwriter.flush(out_strm);
		out_strm.close();
	  }
	catch (IOException e) 
	  {
		throw new RuntimeException("couldn't close streams");
	  }

  }
}