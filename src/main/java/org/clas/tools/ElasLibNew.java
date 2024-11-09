package org.clas.tools;

public class ElasLibNew {
	
	double gep;            // electric form factor of proton
	double gmp;            // magnetic form factor of proton
	double  pi = 3.1415926f;
	double rmp = 0.93827f;
	
	public boolean debug=false;
	public double z0sum=0, z1sum=0, z2sum=0;
	public SpenceFunction spencefunction = new SpenceFunction();
 
    public ElasLibNew() {
		
    }
        
    public class SpenceFunction {

        public double spence(double z) {
            if (z <= 0) {
                throw new IllegalArgumentException("Spence function is not defined for non-positive values.");
            }

            if (z == 1) {
                return 0;
            }

            double sum = 0;
            double term = z;
            int n = 1;

            while (Math.abs(term) > 1e-10) {
                sum += term / (n * n);
                term *= z / (n + 1);
                n++;
            }

            return sum;
        }

        public void main(String[] args) {
            double z = 0.5;
            double result = spence(z);
            System.out.println("Spence(" + z + ") = " + result);
        }
    }
    
    public double newspence(double z) {
    	return spencefunction.spence(z);
    }
	
    public double elas(double eb, double theta, int param) {

//  Returns the e-p differential cross section (ub/sr)
//  eb: beam energy (GeV)
//  theta: e- scattering angle (deg)
//  param: 1=dipole ff 2=Bosted ff 3=Brash ff 7=A1-Mainz

        double fsc,fsc2,theta2,s2,s22,c2,c22,s24,t2,t22,recoil,e_prime,q2,rmott,tau,eps,sig2;

    	fsc  = 1./137; fsc2=fsc*fsc; 

    	theta2 = Math.toRadians(theta)/2.;
    	s2 = Math.sin(theta2); s22=s2*s2; s24=s2*s2*s2*s2;
    	c2 = Math.cos(theta2); c22=c2*c2;
    	t2 = Math.tan(theta2); t22=t2*t2;
    
    	recoil = 1.+(2.*eb/rmp)*s22;
    	e_prime = eb/recoil;
    	q2 = 4.*eb*e_prime*s22;
    
        geteff(q2,param);

    	rmott = 389*fsc2*c22/(4.*eb*eb*s24);   
    	rmott = rmott/recoil;
    
    	tau = q2/(4.*rmp*rmp);
    	eps = 1+2*(1+tau)*t22;
    	eps = 1./eps;;
    
//    	sig1 = (gep*gep+tau*gmp*gmp)/(1.+tau)+2.*tau*gmp*gmp*t22;
    	sig2 = (tau*gmp*gmp+eps*gep*gep)/eps/(1+tau);
     
    	return rmott*sig2;
    
    }
    
    public double recoil(double eb, double theta) {   
    	double theta2 = Math.toRadians(theta)/2.;
    	double s2 = Math.sin(theta2);
    	double rec = 1.+(2.*eb/rmp)*s2*s2;    
    	return eb/rec;    
    }
    
    public double elasq(double eb, double q2, int param) {  
    	double cx = 1-q2/(2*eb*(eb-q2/2/rmp));
//    	double ep = eb-q2/2/rmp;
//      float jac = pi*(1+q2/(2*ep*rmp))/(ep*eb);
    	double epr = eb/(1+eb*(1-cx));
    	double jac = epr*epr/pi;
    	double th  = Math.toDegrees(Math.acos(cx));    
    	return elas(eb,th,param)*jac;
    }
    
    public double elasradq(double eb, double q2, double t, double wcut, int param) {
    	double cx = 1-q2/(2*eb*(eb-q2/2/rmp));
//    	double ep = eb-q2/2/rmp;
//      float jac = pi*(1+q2/(2*ep*0.938))/(ep*eb);
    	double epr=eb/(1+eb*(1-cx));
    	double jac=epr*epr/pi;
    	double th=Math.toDegrees(Math.acos(cx));    
    	return elasrad(eb,th,0,t,0.,0.,wcut,param,0)*jac;
    }

    
    public double asym(double e0, double q2, int param, int out) {
    	double[] aout = new double[3];
    	double nu 	= q2/2/rmp;
    	double e_pr	= e0-nu;
    	if (e_pr<=0.) return 0;
    	double tau	= q2/4/rmp/rmp;
    	double v_l	= 1/Math.pow(1+tau,2);
    	double snth2	= q2/4/e0/e_pr;
    	double csth2	= 1-snth2;
    	if (csth2<=0 || snth2<=0) return 0;
    	snth2	= Math.sqrt(snth2);
    	csth2	= Math.sqrt(csth2);
    	double tnth2	= snth2/csth2;
//    	double theta_e	= Math.toDegrees(Math.asin(snth2));
    	double x		= tnth2*tnth2;
    	double y		= 1./(1.+tau);
    	double v_t	= y/2 + x;
    	double v_tpr	=  tnth2*Math.sqrt(x+y);
    	double v_tlpr= -tnth2*y/Math.sqrt(2.);
    
    	geteff(q2,param);
    
    	double ge2	= gep*gep;
    	double gm2	= gmp*gmp;
    	double den	= v_l*(1.+tau)*ge2/gm2 + 2.*tau*v_t;
    	double qvec	= Math.sqrt((1.+tau)*q2);
    	double  snth_gam	= 2*snth2*csth2*e_pr/qvec;
    	double  csth_gam	= Math.sqrt(1. - snth_gam*snth_gam);
//    	double theta_gam	= Math.toDegrees(Math.asin(snth_gam));
    	aout[0]	=  2*tau*v_tpr/den;
    	aout[1]	= -2*Math.sqrt(2.*tau*(1.+tau))*gep/gmp*v_tlpr/den;
    	aout[2]	= aout[0]*csth_gam+aout[1]*snth_gam;    
    	return aout[out];
    }

/*    
    real function getff2(q2,opt,nf)
    integer nf
    real q2,opt
    
    call ff(q2,opt,gmp,gep)
    if (nf.eq.1) getff2=gep
    if (nf.eq.2) getff2=gmp 
    
    return
    end
*/    
    public void geteff(double q2, int opt) {
    	double corr1, corr2, q1 = Math.sqrt(q2);
    
    	switch (opt) {    
    	case 1: // dipole:
    		gep	= 1./Math.pow((1.+q2/0.71),2);
    		gmp	= 2.7928*gep;
    		break;
    	case 2:// Bosted parameterization of form factors: Phys. Rev. C 51, 409
    		corr1 = 1 + 0.62*q1 + 0.68*q2 + 2.8*q2*q1 + 0.83*q2*q2;
    		corr2 = 1 + 0.35*q1 + 2.44*q2 + 0.5*q2*q1 + 1.04*q2*q2 + 0.34*q2*q2*q1;
    		gep   = 1./corr1;
    		gmp   = 2.7928/corr2;
            break;
    	case 3: //Brash parameterization based on Hall A GeP/GmP measurement: hep-ex\0111038 PRC 65, 051001
    		corr1 = 1 + 0.116*q1 + 2.874*q2 + 0.241*q2*q1 + 1.006*q2*q2 + 0.345*q2*q2*q1;
    		corr2 = 1 - 0.13*(q2-0.04);
    		gep   = corr2/corr1;
    		gmp   = 2.7928/corr1;
            break;
    	case 4:     
    		corr2 = 1 + 0.35*q1 + 2.44*q2 + 0.5*q2*q1 + 1.04*q2*q2 + 0.34*q2*q2*q1;
    		gmp   = 2.7928/corr2;
    		gep   = 0.0;
    		break;
    	case 5: 
    		gep = 1.041/Math.pow(1+q2/0.765,2) - 0.041/Math.pow(1+q2/6.2,2);
    		gmp = 1.002/Math.pow(1+q2/0.749,2) - 0.002/Math.pow(1+q2/6.0,2);
    		gmp = 2.7928*gmp;
    		break;
    	case 6:
    		gep = 1.041/Math.pow(1+q2/0.765,2) - 0.041/Math.pow(1+q2/6.2,2);
    		gmp = 1.002/Math.pow(1+q2/0.749,2) - 0.002/Math.pow(1+q2/6.0,2);
    		double gepb = Math.exp(-0.5*Math.pow((q1-0.07)/0.27,2)) + Math.exp(-0.5*Math.pow((q1+0.07)/0.27,2));
    		gep = gep - 0.3*q2*gepb;    
    		double gmpb = Math.exp(-0.5*Math.pow((q1-0.35)/0.21,2)) + Math.exp(-0.5*Math.pow((q1+0.35)/0.21,2));
    		gmp = 2.7928*(gmp - 0.16*q2*gmpb);
    		break;
//    	case 7: //Direct extraction of Ge and Gm from A1-Mainz data (Bernauer, arXiv:1307.6227v1)          
//      gep = divdif(spige,spix,1000,q2,2)
//      gmp = 2.7928*divdif(spigm,spix,1000,q2,2)
    	}

    	return;
    }
    
    public double spence(double x) {
    	double spence;
    	if (Math.abs(x)<0.1) {
    		spence = x+x*x/4.;
    	} else if (x>0.99 && x<1.01) {
    		spence = pi*pi/6.;
    	} else if (x>-1.01 && x<-0.99) {
    		spence = -pi*pi/12.;
    	} else if (x>0) {
    		spence =  0.1025+sintp(x);
    	} else {
    		spence = -0.0975+sintn(x);
    	}    
    	return spence;
    }
    
    public double sintp(double x) {
    	double   arg,sum = 0.;
    	double xstep = (x-0.1)/100.;
    	double     y = 0.1-xstep/2.;
    	for (int i=1; i<101; i++) {
    		y = y+xstep;
    		arg = Math.abs(1.-y);
    		sum = sum - Math.log(arg)/y;
    	}
    	return sum*xstep;
    }
    
    public double sintn(double x) {
    	double   sum = 0.;
    	double    xa = Math.abs(x);
    	double ystep =(xa-0.1)/100.;
    	double     y = 0.1-ystep/2.;
    	for (int i=1; i<101; i++) {
    		y 	= y+ystep;
    		sum = sum - Math.log(1+y)/y;
    	}
    	return sum*ystep;
    }
    
    public double bfunc(double z) {
    	double xi1 = Math.log(1440.) - 2.*Math.log(z)/3.;
    	double xi2 = Math.log(183.)  -    Math.log(z)/3.;
    	double xi  = xi1/xi2;
    	return (4./3.)*(1.+(z+1.)/(z+xi)/xi2/9.);
    }
    
    public double elasrad(double es, double theta_d, int izn, double t1, double t2, double t3, double wcut, int param, int rc) {
    
//  Returns radiated cross section using equation II.6 of Mo and Tsai, Rev. Mod. Phys. 41, 205-235 (1969).

//  es: incident electron energy (GeV)
//  theta_d: scattered electron angle (deg)
//  izn: index for target nucleus (1,2,3=H,C,Ta)
//  t1: target thickness (radiation lengths)
//  t2: outgoing e- pathlength (r.l.)
//	t3: windows (has b factor included)
//  wcut: used to calculate maximum radiated photon energy. wcut should be >  W resolution (GeV).

//  This program calls function elas defined above.
    
    	double mp,wc,cst1,eel,epr,epcut,gamma4,beta4,e1,e3,e4,eta,qs,theta;
    	double delta,del_mo,delta_t1,delta_t2,delta_t3,delta_t;
    	double radcor,radcor1,sigunp;
    	double den,arg,arg11,arg15,arg19,arg23;
    	double[] deltac = new double[28];
    	double[] tmass = {1.007276470,12.00,180.94788};
    	int znuc, znuc2; int[] znucc = {1,6,73};
    	double me=0.000511, me2=me*me, alpha=7.299e-3, amu=0.931494061;
        int[] iz0 = {0,1,3,24}, iz1 = {2,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}, iz2 = {4,5,6,7,25,26,27};
        
    	theta = Math.toRadians(theta_d);

    	znuc  = znucc[izn]; znuc2 = znuc*znuc;
    	mp    = tmass[izn]*amu;
    	wc    = mp+wcut; // wc should be equal to W cut used in elastic analysis
    
    	cst1  = 1.-Math.cos(theta);
    	eel	  = es/(1.+es*cst1/mp);
    
    	epcut = (es+(mp*mp-wc*wc)/2/mp)*eel/es;
    	delta = wcut==0 ? 0.01*eel : eel-epcut; //maximum radiated energy by photon or dE/E=1% if wcut=0
    
    	epr	  = es+mp-eel;

    	gamma4 = epr/mp;
    	beta4  = Math.sqrt(1.-1/Math.pow(gamma4,2));
    	e1	   = es;
    	e3	   = eel;
    	e4	   = epr;
    	eta	   = es/eel;
    	qs	   = 2.*es*eel*cst1;
    
    	deltac[0] = 28./9.-13./6.*Math.log(qs/me2);
    	arg=2*Math.log(e1/delta)-3*Math.log(eta);
    	deltac[1] = (Math.log(qs/me2) - 1.)*arg; 
    	deltac[2] = 2*znuc*Math.log(eta)   *arg;
    	arg=(e3-e1)/e3;
    	deltac[3]=-spence(arg);
    	
    	deltac[4]=-znuc2 * Math.log(e4/mp);
    	deltac[5]= znuc2 * Math.log(mp/eta/delta) * (Math.log((1.+beta4)/(1.-beta4))/beta4-2.);
    	deltac[6]= znuc2 * 0.5/beta4 * (Math.log((1+beta4)/(1.-beta4)) * Math.log((e4+mp)/2/mp));
    	arg=Math.sqrt((e4-mp)*(1.+beta4)/(e4+mp)/(1.-beta4));
    	deltac[7]=-znuc2*spence(-arg)/beta4;
    	
    	arg=(mp-e3)/e1;
    	deltac[8]=znuc*spence(-arg);
    	den=2*e3*e4-mp*e1;
    	arg=  mp*(mp-e3)/den;
    	deltac[9]=-znuc*spence(arg);
    	arg=2*e3*(mp-e3)/den;
    	deltac[10]=znuc*spence(arg);
    	arg11=Math.abs(den/e1/(mp-2*e3));
    	deltac[11]=znuc*Math.log(arg11)*Math.log(mp/2/e3);
    	
    	arg=(e4-e3)/e3;
    	deltac[12]=-znuc*spence(arg);
    	den=2*e1*e4-mp*e3;
    	arg=(e4-e3)*mp/den;
    	deltac[13]=znuc*spence(arg);
    	arg=2*e1*(e4-e3)/den;
    	deltac[14]=-znuc*spence(arg);
    	arg15=Math.abs(den/e3/(mp-2*e1));
    	deltac[15]=-znuc*Math.log(arg15)*Math.log(mp/2/e1);
    	
    	arg=(mp-e1)/e1;
    	deltac[16]=-znuc*spence(-arg);
    	deltac[17]= znuc*spence( arg);
    	arg=2.*(mp-e1)/mp;
    	deltac[18]=-znuc*spence(arg);
    	arg19=Math.abs(mp/(2*e1-mp));
    	deltac[19]=-znuc*Math.log(arg19)*Math.log(mp/2/e1);
    	
    	arg=(mp-e3)/e3;
    	deltac[20]= znuc*spence(-arg);
    	deltac[21]=-znuc*spence( arg);
    	arg=2*(mp-e3)/mp;
    	deltac[22]=znuc*spence(arg);
    	arg23=Math.abs(mp/(2*e3-mp));
    	deltac[23]=znuc*Math.log(arg23)*Math.log(mp/2/e3);
    	
    	arg=(e1-e3)/e1;
    	deltac[24]=-spence(arg);
    	arg=(e4-mp)*(1-beta4)/(e4+mp)/(1+beta4);
    	arg=Math.sqrt(arg);
    	deltac[25]=znuc2*spence(arg)/beta4;
    	arg=(e4-mp)/(e4+mp);
    	arg=Math.sqrt(arg);
    	deltac[26]=-znuc2*spence( arg)/beta4;
    	deltac[27]= znuc2*spence(-arg)/beta4;

    	z0sum=0; z1sum=0; z2sum=0;
    	for (int i=0; i<iz0.length; i++) z0sum=z0sum+deltac[iz0[i]];
    	for (int i=0; i<iz1.length; i++) z1sum=z1sum+deltac[iz1[i]];
    	for (int i=0; i<iz2.length; i++) z2sum=z2sum+deltac[iz2[i]];
    	z0sum = -alpha*z0sum/pi; z1sum = -alpha*z1sum/pi; z2sum = -alpha*z2sum/pi;
   	
    	del_mo  = 0; 
    	for (int idel=0; idel<deltac.length; idel++) del_mo = del_mo + deltac[idel];
    	del_mo  = -alpha*del_mo/pi;

    	sigunp = izn==1 && rc==0 ? elas(es,theta_d,param) : 1.0;

    	delta_t1 = -0.5*bfunc(znuc)*t1*(Math.log(es/eta/eta/delta));
    	delta_t2 = -0.5*bfunc(znuc)*t2*(Math.log(eel/delta));
    	delta_t3 = -                t3*(Math.log(eel/delta));
    	delta_t  = delta_t1+delta_t2+delta_t3; // Straggling correction 
    
    	radcor1  = 1. + del_mo+delta_t;
    	radcor   = Math.exp(del_mo + delta_t);
    	
    	if(debug) {
        System.out.println(" ");
    	System.out.println(del_mo+" "+delta_t);
    	System.out.println(es+" "+eel+" "+del_mo+" "+delta_t1+" "+delta_t2+" "+delta_t3);
    	System.out.println(1/radcor+" "+1/radcor1);
    	}

    	return radcor*sigunp;
    }
    
    public void test1() {
    	
    	debug=true;
    	
    	double eb = 7.546; 
    	
    	for (int i=1; i<30; i++) {
    		double theta = (i-1)*1+6.5;
    		System.out.println(theta+" "+1e3*elas(eb, theta, 2)+" "+elasrad(eb,theta,1,0.0058,0.0058,0.0,0.0617,2,1));
    	}    	
    }
    
    public void test2() {
    	
    	debug=false;

        double[] eb = {17.314, 15.999, 14.649, 13.329, 11.999, 10.723, 6.032, 2.201, 2.206, 1.645};
    	double[] angle = {35.1, 19.7, 18.8, 17.6, 16.082, 14, 17.186, 38.601, 15.999, 12.0};
    	
    	for (int i=0; i<eb.length; i++) {
    		elasrad(eb[i],angle[i],0,0.0058,0.0058,0.0,0,2,1); //wcut=0.1117 for W=1.05 GeV
    		System.out.println("energy="+eb[i]+" "+angle[i]+" "+z0sum+" "+z1sum+" "+z2sum);
    	}
    }
    
    public void test3() {
    	double[] x = {0.1,0.2,0.3,0.8,1.0,2.0,3.0,10,20,30,40};
    	for (int i=0; i<x.length; i++) System.out.println(x[i]+" "+spence(x[i])+" "+spence(-x[i]));
    }
    
    public static void main(String[] args) {    	
    	ElasLibNew elib = new ElasLibNew();
    	elib.test2();    	
    }

}
