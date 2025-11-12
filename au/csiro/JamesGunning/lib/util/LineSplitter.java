package au.csiro.JamesGunning.lib.util;

/**
   Utility routine for parsing string sLine: 
   splits into 2 chunks separated by colons.
   code : value
   If code starts with % or !, sets bComment = true
   
*/

public class LineSplitter
{
  public boolean bComment=false;
  public String sCode=null;
  public String sValue=null;


  public void parse(String sLine)
  {
	bComment=false;
	if (sLine.length() == 0) 
	  {
		bComment=true;
		return;
	  }
	// split at colon
	String[] st = sLine.split(":");
	sCode = st[0].trim();
	if (sCode.startsWith("%") || sCode.startsWith("!")) 
	  {
		bComment=true;
		return;
	  }
	if (sCode.endsWith(":")) sCode = sCode.substring(0,sCode.length()-2) ;
	sCode = sCode.toLowerCase();	  
	if (st.length>1) 
	  {
		int i;
		sValue = st[1];
		if ((i=sValue.indexOf("!")) != -1) 
		  // test for & strip trailing comment, like "blah: foo !stuff"
		  {
			sValue = sValue.substring(0,i-1);
		  }
		sValue = sValue.trim();
	  }
  }
}
