package au.csiro.JamesGunning.lib.util;

import static java.lang.Math.*;

import au.csiro.JamesGunning.lib.ErrorHandler;
import au.csiro.JamesGunning.lib.OptionsParser;
import static  au.csiro.JamesGunning.lib.util.MathUtil.sqr;

/**
   Basic class for Geophysics function, such as Zoeppritz
   reflectivities

@author james.gunning@dpr.csiro.au
@version 1.0, 26/7/2008 */

public class GeophysicsUtils
{
  /**
	 Zoeppritz Rpp OR Rps of welded interface between 2 isotropic layers, incident wave in
	 1 impinging on layer 2, from Aki & Richards. See diagram etc
	 Mavko et al, Rock Physics Handbook, section 3.4.
	 Rock properties are in
	 vp1,vp2,vs1,vs2,rho1,rho2, 
	 with obvious meanings.
	 Incident wave has angle sin^2(theta_1) = sinsq_th1.

	 If doRps==false, returns Rpp, else returns Rps

	 Returns -Double.MAX_VALUE for a critical reflection, else the documented Rpp value.

	 Visually checked with Crewes applet on 26/7/2008
	 http://www.crewes.org/ResearchLinks/ExplorerPrograms/ZE/ZEcrewes.html

	 Other possible sources are in the public domain CREWES matlab files, e.g.
	 zoeppritz.m
*/
  public static final double zoeppritz_Rpp_Rps(boolean doRps,double sinsq_th1,double vp1,double vp2,double vs1,double vs2,double rho1,double rho2)
  {
	double sin_th1;
	double cos_th1,cos_th2;
	double sinsq_phi1;
	double sin_phi1;
	double sinsq_phi2;
	double sin_phi2;
	double cos_phi1;
	double cos_phi2;

	double a,b,c,d,D,E,F,G,H,p,R;
	double tmp1,tmp2,tmp3;

	sin_th1 = sqrt(sinsq_th1);
	cos_th1 = sqrt(1.0 - sinsq_th1);
	p = sin_th1 / vp1;
	sin_phi1 = vs1 * p; // expect sin_phi1 = (vs1/vp1) * sin_th1 \in [-1,1] by safety of vs1/vp1
	sin_phi2 = vs2 * p; // Not guaranteed here, but violation practically very unlikely
	if (abs(sin_phi2) > 1.0) 
	  {
		return -Double.MAX_VALUE; // critical
	  }
	cos_th2 =  vp2 * p; // at this point, is sin_th2
	if (abs(cos_th2) > 1.0) 
	  {
		return -Double.MAX_VALUE;  // critical
	  }
	cos_th2 = sqrt(1.0 - cos_th2*cos_th2);


	E = sin_phi1*sin_phi1;
	F = sin_phi2*sin_phi2;
	cos_phi1 = sqrt(1.0 - E);
	cos_phi2 = sqrt(1.0 - F);

	b = rho2*(1.0 - 2.0*F);
	c = rho1*(1.0 - 2.0*E);
	a = b - c;
	b += 2.0*rho1*E;
	c += 2.0*rho2*F;
	d = 2.0*(rho2*vs2*vs2 - rho1*vs1*vs1);

	tmp1 = b * cos_th1/vp1;
	tmp2 = c * cos_th2/vp2;
	E = tmp1 + tmp2;    // b * cos_th1/vp1 + c * cos_th2/vp2;
	F = b * cos_phi1/vs1 + c * cos_phi2/vs2;
	tmp3 = d * cos_th1 * cos_phi2/(vp1*vs2);
	G = a - tmp3;       // a - d * cos_th1 * cos_phi2/(vp1*vs2);
	H = a - d * cos_th2 * cos_phi1/(vp2*vs1);

	H *= (p*p);      // H <- H*p*p in two equations below

	D = E*F + G*H;  // E*F + G*H*p*p;

	if (doRps) 
	  {
		// Rps formula, Mavko p62
		R = (-2.0*(cos_th1/vp1)*( a*b + tmp2*d*cos_phi2/vs2)*p*vp1)/(vs1*D);
	  }
	else
	  {
		// Rpp formula, Mavko p62
		R = ((tmp1 - tmp2) * F - (a + tmp3) * H)/D;
	  }
	return R;
  }

  /** 
	  Using supplied permeable rock properties
	  {vp_s,vs_s,phi_s}
	  and impermeable properties {vp_m,vs_m,rho_m}, with
	  net-to-gross NG, hydrocarbon props {S_h,vp_h,rho_h}
	  brine props {vp_b,rho_b}
	  grain properties {K_g,rho_g}
	  compute effective-medium triple {vp,vs,rho}=  vp_vs_rho[0,1,2]
	  Expects sane values. Returns false if physics goes feral.
	  Answer vector must be of size 3 or better
  */
  
  public static boolean GassmannBackus(double vp_s,double vs_s,double phi_s,
									   double vp_m,double vs_m,double rho_m,
									   double NG,
									   double S_h,double vp_h, double rho_h,
									   double vp_b,double rho_b,
									   double K_g,double rho_g,
									   double[] vp_vs_rho)
  {
	// Trivial pure shale case
	if (NG<=0)
	  {
		vp_vs_rho[0]=vp_m;
		vp_vs_rho[1]=vs_m;
		vp_vs_rho[2]=rho_m;
		return true;
	  }

	// (1) Fluid substitution

	double vp_eff,vs_eff,rho_eff,M_eff,mu_eff,mu_m,M_m;
	double rho_satb,M_satb,mu_satb,K_satb,K_satfl,rho_satfl,M_satfl,mu_satfl;
	double K_h,K_fl;

	// ref brine + matrix density		
	rho_satb = phi_s * rho_b + (1.0-phi_s) * rho_g;
	// moduli of ref-fluid saturated sand-member rock
	M_satb = rho_satb *  sqr(vp_s);
	mu_satb =  rho_satb *  sqr(vs_s);
	K_satb = M_satb - 1.3333333333333 * mu_satb;

	if (K_satb<0 && phi_s != 0.0) 
	  {			
		return false;
	  }

	mu_satfl = mu_satb;

	// brine bulk mod
	double K_b = rho_b * sqr(vp_b);
		
	// ref fluid bulk mod
	if (K_satb < K_g/(phi_s * (K_g/K_b - 1.0) +1.0) && phi_s!=0.0)
	  {
		System.err.println("K_satb must be at least= " + (K_g/(phi_s * (K_g/K_b - 1.0) +1.0)) + 
						   " under these conditions (N.B. can be large at small porosities)." +
						   " Rock speeds may be too low in trend curves");
		return false;
	  }
		
	// hydrocarbon bulk mod.		
	K_h = rho_h * sqr(vp_h);
		
	// New effective fluid bulk modulus & density...
	// A mixture of hydrocarbon and new brine
	K_fl = 1.0/(S_h/K_h + (1.0-S_h)/K_b);
			
	// do fluid substitution
	if (phi_s !=0)
	  {
		K_satfl = K_g/(1.0+1.0/((1.0/((K_g/K_fl)-1.0) - 1.0/((K_g/K_b)-1.0))/phi_s
								+ 1.0/((K_g/K_satb)-1.0)));
	  }
	else
	  {
		K_satfl = K_g;
	  }
	// new effective density
	rho_satfl = (1.0-phi_s) * rho_g + phi_s * (S_h*rho_h+(1.0-S_h)*rho_b); 
	// new effective p-wave modulus
	M_satfl = K_satfl + 1.3333333333333 * mu_satfl; 
	  

	// (2) Rock mixing
	
	// shale moduli
	mu_m = rho_m * sqr(vs_m);
	M_m = rho_m * sqr(vp_m);
	// NG weighted average for densities
	rho_eff = NG * rho_satfl + (1.0-NG)*rho_m;
	// And moduli are averaged using Backus
	if (NG<1)
	  {
		// Backus average of p-wave modulus.
		M_eff = 1.0/(NG/M_satfl + (1.0-NG)/M_m);
		// Backus average of shear modulus
		mu_eff = 1.0/(NG/mu_satfl + (1.0-NG)/mu_m);
	  }
	else 
	  {
		M_eff = M_satfl;
		mu_eff = mu_satfl;
	  }
	// and resulting equivalent velocities
	vp_eff = Math.sqrt(M_eff/rho_eff);
	vs_eff = Math.sqrt(mu_eff/rho_eff);

	// fill answer vector
	vp_vs_rho[0]=vp_eff;
	vp_vs_rho[1]=vs_eff;
	vp_vs_rho[2]=rho_eff;

  	return true;
  }


  /** 
	  Using supplied permeable rock properties
	  {vp_s,vs_s,phi_s}
	  hydrocarbon props {S_h,vp_h,rho_h}
	  brine props {vp_b,rho_b}
	  grain properties {K_g,rho_g}
	  compute effective-medium 4tuple {vp,vs,rho,E}=  vp_vs_rho_E[0,1,2,3]
      Using Gassmann and Wood's average for fluids
	  Expects sane values. Returns false if physics goes feral.
	  Answer vector must be of size 3 or better
  */
  
  public static boolean GassmannWood(double vp_s,double vs_s,double phi_s,
                                     double S_h,double vp_h, double rho_h,
                                     double vp_b,double rho_b,
                                     double K_g,double rho_g,
                                     double[] vp_vs_rho_E)
  {
	// (1) Fluid substitution

	double vp_eff,vs_eff,rho_eff,M_eff,mu_eff,mu_m,M_m;
	double rho_satb,M_satb,mu_satb,K_satb,K_satfl,rho_satfl,M_satfl,mu_satfl;
	double K_h,K_fl;

	// ref brine + matrix density		
	rho_satb = phi_s * rho_b + (1.0-phi_s) * rho_g;
	// moduli of ref-fluid saturated sand-member rock
	M_satb = rho_satb *  sqr(vp_s);
	mu_satb =  rho_satb *  sqr(vs_s);
	K_satb = M_satb - 1.3333333333333 * mu_satb;

	if (K_satb<0 && phi_s != 0.0) 
	  {			
		return false;
	  }

	mu_satfl = mu_satb;

	// brine bulk mod
	double K_b = rho_b * sqr(vp_b);
		
	// ref fluid bulk mod
	if (K_satb < K_g/(phi_s * (K_g/K_b - 1.0) +1.0) && phi_s!=0.0)
	  {
		System.err.println("K_satb must be at least= " + (K_g/(phi_s * (K_g/K_b - 1.0) +1.0)) + 
						   " under these conditions (N.B. can be large at small porosities)." +
						   " Rock speeds may be too low in trend curves");
		return false;
	  }
		
	// hydrocarbon bulk mod.		
	K_h = rho_h * sqr(vp_h);
		
	// New effective fluid bulk modulus & density...
	// A mixture of hydrocarbon and new brine
	K_fl = 1.0/(S_h/K_h + (1.0-S_h)/K_b);
			
	// do fluid substitution
	if (phi_s !=0)
	  {
		K_satfl = K_g/(1.0+1.0/((1.0/((K_g/K_fl)-1.0) - 1.0/((K_g/K_b)-1.0))/phi_s
								+ 1.0/((K_g/K_satb)-1.0)));
	  }
	else
	  {
		K_satfl = K_g;
	  }
	// new effective density
	rho_satfl = (1.0-phi_s) * rho_g + phi_s * (S_h*rho_h+(1.0-S_h)*rho_b); 
	// new effective p-wave modulus
	M_satfl = K_satfl + 1.3333333333333 * mu_satfl; 
	  
	// and resulting equivalent velocities
	vp_eff = Math.sqrt(M_satfl/rho_satfl);
	vs_eff = Math.sqrt(mu_satfl/rho_satfl);

	// fill answer vector
	vp_vs_rho_E[0]=vp_eff;
	vp_vs_rho_E[1]=vs_eff;
	vp_vs_rho_E[2]=rho_satfl;

	// Add also Young modulus: E_satfl Mavko table 2.1.1
	vp_vs_rho_E[3] = mu_satfl*(3.0*M_satfl - 4.0*mu_satfl)/(M_satfl - mu_satfl);

	return true;
  }



  /** 
	Sample driver for Zoeppritz calculation: can plot AVO reponses to stdout
  */
  public static void main(String argv[])
  {
	int i,testNo;
	double vp1=1,vp2=1,vs1=1,vs2=1,rho1=1,rho2=1;
	double rhoBar,vpBar ,vsBar ,delta_rho ,delta_vp ,delta_vs ,R_near , ratio;
	double th, th1,sinsq_th1,phiAv,tmp,p,thAv;
	double Rz,R;

    double vp_h,rho_h,S_h,NG,vp_b,rho_b,rho_g, K_g,phi_s,Z1,Z2;
    double vp_vs_rho[] = new double[4]; // May gave E in GassmannWood

    boolean showARPS=false;

	ErrorHandler EH = new ErrorHandler();
	OptionsParser options = new OptionsParser("GeophysicsUtils",EH);
	
	options.addSpielText("Computes reflection coefficients for chosen cases. "
                         +"-p: prints ascii stream of AVO response, PP,PS reflectivity, and Zoeppritz+Linearised-Zoeppritz");
    options.addOption("-c",null,1,"int","case","Do test case (case 1,2, or 3) showing Zoeppritz vs approximations");
	options.addOption("-A",null,0,null,null,"Show Aki-Richards --like Rps");
    options.addOption("-p",null,1,"String","vp1,vs1,rho1,vp2,vs2,rho2","use cmdline properties vp1,vs1,rho1,vp2,vs2,rho2 (m/s, kg/m3) for the two layers");
    options.addOption("-g",null,1,"String","vp1,vs1,rho1,vp2,vs2,rho2,vp_h,rho_h,S_h","use cmdline properties vp1,vs1,rho1,vp2,vs2,rho2,vp_h,rho_h,S_h (m/s, kg/m3), for fluid substitution in layer 2 below layer 1 (assumes generic quartz and brine properties). E.g. typical methane sub: -g 3000,1500,2300,3000,1500,2300,500,100,0.8");
    options.addOption("-f",null,1,"String","vp1,vs1,rho1,vp_h,rho_h","use cmdline properties vp1,vs1,rho1,vp_h,rho_h (m/s, kg/m3), for fluid substitution in layer 1, produce table of fluid response with S_h (assumes generic quartz and brine properties). E.g. typical methane sub: -f 3000,1500,2300,500,100");

	options.parseOptions(argv,true);

	testNo = 1;

    if (options.isOptionSet("-A")) showARPS=true;
	if (options.isOptionSet("-c")) 
	  {
		testNo = options.getOptionInt("-c",0);
	  }
	if (options.isOptionSet("-p")) 
	  {
		String s = options.getOptionString("-p",0);
        String[] a = s.split(",");
        vp1 = Double.parseDouble(a[0]);
        vs1 = Double.parseDouble(a[1]);
        rho1 = Double.parseDouble(a[2]);
        vp2 = Double.parseDouble(a[3]);
        vs2 = Double.parseDouble(a[4]);
        rho2 = Double.parseDouble(a[5]);
      }
    else if (options.isOptionSet("-g")) 
      {
        // Two layers with fluid substitution
        String s = options.getOptionString("-g",0);
        String[] a = s.split(",");
        // Examples
        // vp_h=350; rho_h=370; // co2 at ~15Mpa        
        // vp_h=1400; rho_h=12; //H2
        // vp_h=500; rho_h=103; //CH4

        vp1 = Double.parseDouble(a[0]);
        vs1 = Double.parseDouble(a[1]);
        rho1 = Double.parseDouble(a[2]);

        vp2 = Double.parseDouble(a[3]);
        vs2 = Double.parseDouble(a[4]);
        rho2 = Double.parseDouble(a[5]);

        vp_h = Double.parseDouble(a[6]);
        rho_h = Double.parseDouble(a[7]);
        S_h = Double.parseDouble(a[8]);

        NG=1;
        
        // brine
        vp_b = 1500;
        rho_b = 1000;
        // quartz grain
        rho_g = 2650;
        K_g = 38e9; //38GPa
        
        phi_s = (rho_g - rho2)/(rho_g-rho_b);
     
        Z1 = vp1*rho1;
        GassmannBackus(vp2,vs2,phi_s, vp1,vs1,rho1,NG,S_h/*S_h*/,vp_h, rho_h,vp_b,rho_b, K_g,rho_g, vp_vs_rho);
        Z2 = vp_vs_rho[0]*vp_vs_rho[2];
        System.out.printf("Layer2: phi=%f vp=%f rho=%f AI=%f %f\n",phi_s,vp_vs_rho[0],vp_vs_rho[2],Z2,K_g*1e-9);
        System.out.printf("Reff=%f\n",(Z2-Z1)/(Z1+Z2));

        System.exit(0);
      }
    
    else if (options.isOptionSet("-f")) 
      {
        // Sh dependence layer1 with fluid substitution
        String s = options.getOptionString("-f",0);
        String[] a = s.split(",");
        // Examples
        // vp_h=350; rho_h=370; // co2 at ~15Mpa        
        // vp_h=1400; rho_h=12; //H2
        // vp_h=500; rho_h=103; //CH4

        vp1 = Double.parseDouble(a[0]);
        vs1 = Double.parseDouble(a[1]);
        rho1 = Double.parseDouble(a[2]);
        vp_h = Double.parseDouble(a[3]);
        rho_h = Double.parseDouble(a[4]);
        NG=1;
        
        // brine
        vp_b = 1500;
        rho_b = 1000;
        // quartz grain
        rho_g = 2650;
        K_g = 38e9; //38GPa        
        phi_s = (rho_g - rho1)/(rho_g-rho_b);

        System.out.printf("#S_h  phi vp  rho  AI E\n");
        for(S_h=0;S_h<=1;S_h+=0.01)
          {
            GassmannWood(vp1,vs1,phi_s,S_h,vp_h, rho_h,vp_b,rho_b, K_g,rho_g, vp_vs_rho);
            Z1 = vp_vs_rho[0]*vp_vs_rho[2];
            System.out.printf("%f %f %f %f %f %f\n",S_h,phi_s,vp_vs_rho[0],vp_vs_rho[2],Z1,vp_vs_rho[3]);

          }
        System.exit(0);
      }
    else
      {
        if (testNo==1)
          {
            vp1=3300;
            vp2=4200;
            vs1=1700;
            vs2=2700;
            rho1=2.35;
            rho2=2.49;
          }
        else if (testNo==2)
          {
            vp1=3000;
            vp2=4000;
            vs1=1500;
            vs2=2000;
            rho1=2000;
            rho2=2200;	
          }
        else if (testNo==3)
          {
            vp1=2600;
            vp2=3000;
            vs1=1500;
            vs2=1700;
            rho1=2000;
            rho2=2200;	
          }
        else
          {
            EH.die_err("Illegal test case "+testNo);
          }
      }
	System.err.printf("Layer 1 props(vp vs rho): %10.5g %10.5g %10.5g%n",vp1,vs1,rho1);
	System.err.printf("Layer 2 props(vp vs rho): %10.5g %10.5g %10.5g%n",vp2,vs2,rho2);
	
	for(int pass=0;pass<2;pass++)
	  {
		// PP case, then PS
		boolean doPS = (pass==0);
		for (i=0;i<2;i++)
		  {
			// Zoeppritz, then approx...
			for (th=0; th<=65; th+=1.0)
			  {
				th1 = toRadians(th);  // = theta_1
				p = sin(th1)/vp1;
				sinsq_th1 = sqr(sin(th1));
				thAv = 0.5*(th1+asin(vp2*p));
				// (phi1+phi2)/2
				phiAv = (asin(vs1*p) + asin(vs2*p))/2.0;
				if (i==0)
				  {
					R = zoeppritz_Rpp_Rps(doPS,sinsq_th1,vp1,vp2,vs1,vs2,rho1,rho2);
				  }
				else
				  {
					rhoBar = 0.5*(rho2+rho1);
					vpBar = 0.5*(vp2+vp1);
					vsBar = 0.5*(vs2+vs1);
					delta_rho = rho2-rho1;
					delta_vp = vp2-vp1;
					delta_vs = vs2-vs1;
					double Dp = delta_vp/vpBar;
					double Ds = delta_vs/vsBar;
					double Dr = delta_rho/rhoBar;
					double lam = vsBar/vpBar;
					double ls=sqr(lam);
					if (!doPS)
					  {
						// PP case
						R_near = 0.5*(Dr + Dp);
						// Linearized Zoepritz equation
						ratio = Dp/2.0 - (2.0*ls*(Dr+2.0*Ds));
						R = R_near + th1*th1*ratio;
					  }
					else
					  {
						if (showARPS)
						  {
							double cc = lam*cos(thAv)*cos(phiAv);
							if (true)							  
							  {
								// everything in terms of averages..(CREWES)								
								tmp = sqr(sin(phiAv));
								R = -sin(thAv)/(2*cos(phiAv));
								R *= ((1 - 2*tmp + 2*cc)*Dr
									  -4*(tmp-cc)*Ds);
							  }
							else
							  {
								// Aki-Richards literally
								tmp = sqr(p*vsBar);
								R = -p*vpBar/(2*cos(phiAv));
								R *= ((1 - 2*tmp + 2*cc)*Dr
									  -4*(tmp-cc)*Ds);
							  }
						  }
						else
						  {
							// linear+cubic
						    tmp=Dr+2*Ds;
							R = th1*(-0.5*(Dr + 2.0*lam*tmp) 
									 + sqr(th1) * ((1.0 + lam*(8 + 9.0*lam))*Dr  +  lam*(16.0 + 24.0*lam)*Ds)/12.0);
							
						  }	
					  }				  
				  }
				System.out.printf("%8g %8g%n",th,(R>-10)?R:0);
			  }
			System.out.println("");
		  }
	  }
  }
}
