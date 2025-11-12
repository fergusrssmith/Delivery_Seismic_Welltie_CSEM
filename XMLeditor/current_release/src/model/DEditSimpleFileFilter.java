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

import java.io.File;
import javax.swing.filechooser.FileFilter;

/**
  * This is a really simple file filter that accepts only one extension and
  * directories.  Add more extensions???
  */
public class DEditSimpleFileFilter extends FileFilter {
    
    private String extension; // extension to accept
    private String description; // description of file filter

    /**
      * Creates a simple file filter.
      * @param extension the extension to allow (verified by String.endsWith())
      * @param description the description of the FileFilter to be returned
      * by getDescription().
      */
    public DEditSimpleFileFilter(String extension, String description) {
	this.extension = extension;
	this.description = description;
    }

    /**
      * Checks to see if FileChooser should show this File to be shown
      * @param f the file to be checked
      * @return true if f is a directory or its absolute path .endsWith the 
      * extension string given to the constructor
      */
    public boolean accept(File f) {
	return f.isDirectory() || f.getAbsolutePath().endsWith(extension);
    }

    /**
      * Returns a string description of files allowed by FileFilter.
      * @return a string description of files allowed by FileFilter
      */
    public String getDescription() {
	return description;
    }
}
