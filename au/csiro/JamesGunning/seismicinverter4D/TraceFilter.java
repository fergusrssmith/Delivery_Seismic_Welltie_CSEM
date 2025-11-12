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



package au.csiro.JamesGunning.seismicinverter4D;

import java.util.*;
import java.util.regex.*;

import au.csiro.JamesGunning.sulib.*;
import au.csiro.JamesGunning.lib.ErrorHandler;


/**
   Basic class for filtering trace inputs by applying sequence of tests
   
   @author James.Gunning@csiro.au
   @version 1.0, 9/10/2001 
*/

public final class TraceFilter
{
  boolean applyTraceFilter = false;
  String  cmdLineTraceFilterString = null;
  ArrayList traceFilterList = null;
  ErrorHandler _EH;
  int verbosity;

  /**
	 Constructor: defaults to no trace filter
  */
  TraceFilter(ErrorHandler EH, int verbosity)
  {
	this._EH = EH;
	this.verbosity = verbosity;
  }


  /**
	 Sets up structure required to filter input traces by keyword if desired.
	 At this point, cmdLineTraceFilterString will have been set to a 
	 comma-separated list of logical requirements
	 for the trace header words, e.g. 'ep=2345,tracl<789,mark>34'. Also
	 applyTraceFilter will be set to true, indicating that trace filering is being applied.
	 Must have access to a working SuSegyFormat object tr
  */
  void setUpTraceInputFilter(SuSegyFormat tr)
  {
	if (tr==null)
	  {
		throw new RuntimeException("null SuSegyFormat object tr in TraceFilterSpecifications.setUpTraceInputFilter()");
	  }
	if (applyTraceFilter)
	  {
		String headerWord, operator, number;
		int comparison;
		int matches;
		// parse string into sequence of blocks separated by commas or colons. First, strip white space
		String cmd = cmdLineTraceFilterString.replaceAll(" ",""); 
		StringTokenizer args = new  StringTokenizer(cmd,",:");
		if (args!=null)
		  {
			traceFilterList = new ArrayList();
			// Use reluctant quantifiers X+? is 1 or more matches of X
			Pattern pattern = Pattern.compile("(\\w+?)(>|<|=|!=|<=|>=|%)(\\w+?)");  // "word operator word" , e.g. fldr=4
			while (args.hasMoreTokens())
			  {
				// break up logical expressions into bits with <,=,>,!=,>=,<= operators in between
				String token = args.nextToken();
				Matcher matcher = pattern.matcher(token);
				if (matcher.matches() && ((matches = matcher.groupCount())==3))
				  {
					// must match "headerword operator number" 
					headerWord =  matcher.group(1);
					operator =  matcher.group(2);
					number = matcher.group(3);
				  }
				else
				  {
					throw new RuntimeException("trace--filter expression "+token+" is not a sequence of header word comparisons, e.g. ep=4");
				  }

				try
				  {
					comparison = Integer.parseInt(number);
				  }
				catch (NumberFormatException e)
				  {
					throw new RuntimeException("trace--filter expression comparison is not a number: "+number+", in "+
											   cmdLineTraceFilterString);
				  }
				TraceFilterSpecifications tFS = new TraceFilterSpecifications(headerWord,operator,comparison);
				// verify keywords etc are sane
				if (!tFS.verifyKeyword(tr))
				  {
					throw new RuntimeException("trace--filter expression has illegal keyword: "+headerWord+", in "+
											   cmdLineTraceFilterString);
				  }
				if (!tFS.verifyOperator())
				  {
					throw new RuntimeException("trace--filter expression has illegal operator: "+operator+", in "+
											   cmdLineTraceFilterString);
				  }
				// things are ok
				traceFilterList.add(tFS);
			  }
		  }
	  }
  }


  /**
	 Method that returns a flag indicating whether to skip trace, able
	 to be applied once the array of tests in traceFilterList has been
	 set up correctly in the setUpTraceInputFilter() call (usually at
	 initialisation).  Requires a valid SuSegyFormat object to extract
	 header words from etc, for comparison to logical requirements in
	 the arrays of conditionals embedded in traceFilterList
  */
  boolean skipTraceIfRequired(SuSegyFormat tr)
  {
	if (tr==null)
	  {
		throw new RuntimeException("null SuSegyFormat object tr in TraceFilterSpecifications.skipTraceIfRequired()");
	  }

	boolean skip=false;
	if (!applyTraceFilter) skip=false;
	else
	  {
		// test words in trace etc
		// comparisons are made against seismic near_tr field
		// loop over tests. First failure signifies trace must be skipped
		for (int i=0;i<traceFilterList.size();i++)
		  {
			TraceFilterSpecifications tFS = (TraceFilterSpecifications)traceFilterList.get(i);
			// get the current value of the desired header word. A float's good enough
			double trace_value = tr.getAsFloat(tFS.headerWord);
			if (tFS.operator.equals(">"))
			  {
				if (trace_value <= tFS.comparison) 
				  {
					skip=true;
					break;
				  }
			  }
			else if (tFS.operator.equals("<"))
			  {
				if (trace_value >= tFS.comparison) 
				  {
					skip=true;
					break;
				  }
			  }
			else if (tFS.operator.equals("="))
			  {
				if (trace_value != tFS.comparison)
				  {
					skip=true;
					break;
				  }
			  }
			else if (tFS.operator.equals("!="))
			  {
				if (trace_value == tFS.comparison)
				  {
					skip=true;
					break;
				  }
			  }
			else if (tFS.operator.equals(">="))
			  {
				if (!(trace_value >= tFS.comparison))
				  {
					skip=true;
					break;
				  }
			  }
			else if (tFS.operator.equals("<="))
			  {
				if (!(trace_value <= tFS.comparison))
				  {
					skip=true;
					break;
				  }
			  }
			else if (tFS.operator.equals("%"))
			  {
				if (!(trace_value % tFS.comparison == 0))
				  {
					skip=true;
					break;
				  }
			  }
		  }
	  }
	if (verbosity>2 && skip) _EH.warn("TRACE skipped..");
	return skip;
  }


  /**
	 Small class to encapsulate & verify trace filter expressions
  */
  private static class TraceFilterSpecifications
  {
	String headerWord; // the SU header word used in the logical required e.g. fldr
	String operator;   // the operator, one of <,>,=,!=,%
	int comparison;    // the value to compare headerWord against, e.g. fldr>345
	
	public TraceFilterSpecifications(String headerWord,String operator,int comparison)
	  {
		this.headerWord = headerWord;
		this.operator = operator;
		this.comparison = comparison;
	  }

	// Check if header word seems valid in Su keywords list
	// must be passed any valid SuSegyFormat object with an intact headerWords HashMap
	boolean verifyKeyword(SuSegyFormat tr)
	{
	  if (tr==null)
		{
		  throw new RuntimeException("null SuSegyFormat object tr in TraceFilterSpecifications.verifyKeyword()");
		}
	  if (tr.headerMap.get(headerWord) == null) return false;
	  return true;
	}

	// operator must be one of <,>,=,!=,>=,<=
	boolean verifyOperator()
	{
	  boolean legal=false;
	  if (operator.equals(">")) legal=true;
	  else if (operator.equals("<")) legal=true;
	  else if (operator.equals("=")) legal=true;
	  else if (operator.equals("!=")) legal=true;
	  else if (operator.equals(">=")) legal=true;
	  else if (operator.equals("<=")) legal=true;
	  else if (operator.equals("%")) legal=true;
	  return legal;
	}	
  }
}
