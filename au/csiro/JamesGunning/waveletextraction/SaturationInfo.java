

package au.csiro.JamesGunning.waveletextraction;
/**
   Data class for storing saturation-related well info.
   
   @author James.Gunning@csiro.au
   @version 1.0, 25/11/2011
*/

import au.csiro.JamesGunning.lib.Numerical.PropertyMapper;

public class SaturationInfo
{
  // input values of "mean sat" after-production, and "before", in-log
  double Sa=.999; // Chevron = 0.15
  double Sb=0.001; // Chevron =.9
  // std devs after/before
  double sig_Sa=0.5; // Chevron = 0.5;
  //double sig_Sb=1;   // not used presently

  // working vars...
  PropertyMapper satPropMapper;
  // property exists in logit domain (-inf,inf). Mapped to within [0,1] by mapper
  Property satProp = new Property(Sb,false,false);  
}
