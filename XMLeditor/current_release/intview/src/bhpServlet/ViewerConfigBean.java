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
package bhpServlet;

import java.util.*;
import java.io.*;
import javax.servlet.http.*;
import javax.servlet.ServletContext;

/**
 * Title:        BHP Viewer <br><br>
 * Description:  This class is used in jsp file for configuration. <br><br>
 *
 * Copyright:    Copyright (c) 2001 <br>
 * Company:      BHP INT <br>
 * @author Synthia Kong
 * @version 1.0
 */

public class ViewerConfigBean {
    private static final String DEFAULT_URI_HEADER = "http://";

    // need to find the webapp name dynamically
    private static final String READ_DIRECTORY = "/BHP/servlet/readDirectory";
    private static final String EXECUTE_SCRIPT = "/BHP/servlet/executeScript";
    private static final String SUMMARY_DATA   = "/BHP/servlet/summaryData";
    private static final String READ_DATA      = "/BHP/servlet/getData";
    private static final String SAVE_FILECFG   = "/BHP/servlet/saveConfig";
    private static final String READ_FILECFG   = "/BHP/servlet/openConfig";

    private ServletContext _servletContext;
    private HttpServletRequest _servletRequest;

    private HashMap _appletAttributeValues;
    private HashMap _serverContextValues;

    public ViewerConfigBean() {
        _serverContextValues = new HashMap();
        _appletAttributeValues = new HashMap();
    }

    /**
     * Usage: <jsp:param name = "HelpURL"           value = '<%= param_bean.getValue("HelpURL") %>'/>
     *  SPEC: This function is called in the JSP to pass the value of the attribute to
     *  the jsp:plugin statement. Applet values are values that are used as parameters
     *  when invoking the applet.
     */
    public String getAppletValue(String key) {
        if (!getAppletValueExistence(key)) {
            // Write error to screen
            System.out.println("ERROR: The applet attribute " + key + " Has no corresponding value. Check the configuration file." );
            return "";
        }
        return (String)_appletAttributeValues.get(key);
    }

    /**
     * Usage: getServerValue("HostName");
     *  SPEC: This function is called by the config bean sub-class.
     *  Server Values are ususally server conext parameters,
     *  or other server initialization parameters.
     */
    public String getServerValue(String key) {
        if (!getServerValueExistence(key)) {
            // Write error to screen
            System.out.println("ERROR: The Server parameter " + key + " Has no corresponding value." );
            return "";
        }
        return (String)_serverContextValues.get(key);
    }

    /**
     * Checks for the existance of the specified Applet key
     * @param key String identifier for the key
     */
    public boolean getAppletValueExistence(String key) {
        return _appletAttributeValues.containsKey(key);
    }

    /**
     * Checks for the existance of the specified Server key
     * @param key String identifier for the key
     */
    public boolean getServerValueExistence(String key) {
        return _serverContextValues.containsKey(key);
    }

    /**
     *  Sets the value of the specified Applet key. Usually called from viewerInitBean().
     *  @param key String identifier for the key
     *  @param value New value for the key
     */
    public void setAppletValue(String key, String value) {
        //System.out.println("ViewConfigBean.setAppletValue : " + key + "--" + value);
        _appletAttributeValues.put(key, value);
    }

    /**
     *  Sets the value of the specified Server key
     *  @param key String identifier for the key
     *  @param kalue New value for the key
     */
    public void setServerValue(String key, String value) {
        //System.out.println("ViewConfigBean.setServerValue : " + key + "--" + value);
        _serverContextValues.put(key, value);
    }

    /**
     *  This method loads the properties from a property file into the
     *  _appletAttributeValues hash table. This allows the Admin of the
     *  server to change the values passed to the applet without
     *  recompiling the code.
     *  This function must be called in the .jsp BEFORE veiwerBeanInit()
     *  @param filename name of .prop file.
     */
    public void setAppletPropertyFile(String filename) {
        System.out.println("ViewerConfigBean.setAppletPropertyFile : " + filename);
        Properties configFileValues = new Properties();
        try {
            configFileValues.load(new FileInputStream(filename));
            Enumeration filePropNames = configFileValues.propertyNames();
            while (filePropNames.hasMoreElements()) {
                String keyName = (String)filePropNames.nextElement();
                //if (_appletAttributeValues.containsKey(keyName)) {
                _appletAttributeValues.put(keyName, configFileValues.getProperty(keyName));
                //System.out.println("Attribue " + keyName + " has been read from the config file");
            }
        }
        catch (IOException iexcp) {
            System.out.println("NOTE: No config file found; using defaults." );
            File curDir = new File(".");
            System.out.println("Current directory : " + curDir.getAbsolutePath());
        }
        catch (SecurityException sexcp) {
            System.out.println("WARNING: Config file could not be opened for read" );
        }
        finally {
            //configFileValues.close()???
        }
    }

    /**
    * Save the servlet's context
    * @param svCtx the servlet context of the JSP
    */
    public void setServletContext(ServletContext svrCtx) {
         _servletContext = svrCtx;
    }

    /**
    * Save the servlet's http request, used to automagically snag the
    * page's host and port info
    * @param page_req the http request of the JSP
    */
    public void setRequest(HttpServletRequest pageReq) {
        System.out.println("ViewerConfigBean.setRequest");
         _servletRequest = pageReq;

         // !!! HACK ALERT !!!
         // do this first to force whatever magik needs to happen to
         // the the proper port number, else you will allways get 80 back
         // (at least on Tomcat 3.2.3 anyway)
         // !!! HACK ALERT !!!
         StringBuffer reqInfo = HttpUtils.getRequestURL(pageReq);

         // set the hostname and port as one string to HostName
         String hostname = _servletRequest.getServerName() + ":" +
                            _servletRequest.getServerPort();
         setServerValue("HostName", hostname);
    }

    /**
     * This method is called right before the applet is called.
     * This occurs after the server properties have been passed to the bean.
     * You should use this method to populate the applet parameters values
     * using the setAppletValue method.
     */
    public void viewerBeanInit() {
        String prefix = DEFAULT_URI_HEADER + getServerValue("HostName");
        setAppletValue("servletReadDirectoryURL", prefix + READ_DIRECTORY);
        setAppletValue("servletExecuteScriptURL", prefix + EXECUTE_SCRIPT);
        setAppletValue("servletSummaryDataURL", prefix + SUMMARY_DATA);
        setAppletValue("servletReadDataURL", prefix + READ_DATA);
        setAppletValue("servletSaveFilecfgURL", prefix + SAVE_FILECFG);
        setAppletValue("servletReadFilecfgURL", prefix + READ_FILECFG);
    }

}