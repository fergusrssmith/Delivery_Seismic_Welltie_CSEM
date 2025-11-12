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
package model;

/**
  * Simple debugging class whose only current functionality is to 
  * perform a println or print operation which has a cutoff priority so that 
  * the user can select the verbosity of debug output.
  */
public class DEditDebug {
    
    /* any DEditDebug.println call with a priority >= cutoff_priority will
     * be diplayed, unless cutoff_priority == -1 in which case none will */
    private static int cutoff_priority = -1;

    /**
      * Sets the cutoff priority for displaying debug messages.
      * @param priority the new cutoff priority
      */
    public static void setPriority(int priority) {
	cutoff_priority = priority;
    }
    
    /**
      * Returns the current cutoff priority.
      * @return the current cutoff priority
      */
    public static int getPriority() {
	return cutoff_priority;
    }

    /**
      * Replaces System.out.println() for debugging purposes so that 
      * output can be easily supressed.
      * @param caller the class instance in which this println call was made
      * @param priority the priority of the message to be printed
      * @param message the message to print
      */
    public static void println(Object caller, int priority, Object message) {
	if (cutoff_priority == -1 || cutoff_priority > priority) {
	    return;
	}
	System.out.println(caller.getClass().getName() + ": "
			   + message.toString());
    }

    /**
      * Replaces System.out.print() for debugging purposes so that 
      * output can be easily supressed.
      * @param caller the class instance in which this println call was made
      * @param priority the priority of the message to be printed
      * @param message the message to print
      */
    public static void print(Object caller, int priority, Object message) {
	if (cutoff_priority == -1 || cutoff_priority > priority) {
	    return;
	}
	System.out.print(caller.getClass().getName() + ": "
			   + message.toString());
    }
}
