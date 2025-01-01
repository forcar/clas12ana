package org.clas.tools;

import java.text.DecimalFormat;
import java.util.Formatter;

public class MoTsai {
	
	double gep;             
	double gmp;            
	double  pi = 3.141592654;
	double rmp = 0.93827;
	double rme = 0.00051099895069;
	double fsc = 1./137.01;
	
	public double z0sum=0, z1sum=0, z2sum=0, del_mo, delta_t, xsect_raw, xsect_rad, radcor;
	public double q2, e_prime;
	public int param;
 
    public MoTsai() {
		
    }

    public double elas(double eb, double theta) {

        double fsc2,theta2,s2,s22,c2,c22,s24,t2,t22,recoil,rmott,tau,eps,sig;

    	fsc2 = fsc*fsc; 

    	theta2 = Math.toRadians(theta)/2.;
    	s2 = Math.sin(theta2); s22=s2*s2; s24=s2*s2*s2*s2;
    	c2 = Math.cos(theta2); c22=c2*c2;
    	t2 = Math.tan(theta2); t22=t2*t2;
    
    	recoil = 1.+(2.*eb/rmp)*s22;
    	e_prime = eb/recoil;
    	
    	q2 = 4.*eb*e_prime*s22;    
        getFF(q2,param);

    	rmott = 389*fsc2*c22/(4.*eb*eb*s24);   
    	rmott = rmott/recoil;
    
    	tau = q2/(4.*rmp*rmp);
    	eps = 1/(1+2*(1+tau)*t22);
    
    	sig = (tau*gmp*gmp+eps*gep*gep)/eps/(1+tau);
     
    	return rmott*sig;
    
    }    
    
    public void setFF(int param) {
    	this.param = param;
    }
    
    public double getBCC(double eb, double theta, double dtheta, int nsamples) { //bin centering correction
    	double xsect = elas(eb,theta);
    	double xsum = 0;
    	for (int i=0; i<nsamples; i++) xsum = xsum + elas(eb,theta + dtheta*(Math.random()-0.5));
        return nsamples*xsect/xsum;
    }
    
    public double elasq(double eb, double q2) {  
    	double cx = 1-q2/(2*eb*(eb-q2/2/rmp));
//    	double ep = eb-q2/2/rmp;
//      float jac = pi*(1+q2/(2*ep*rmp))/(ep*eb);
    	double epr = eb/(1+eb*(1-cx));
    	double jac = epr*epr/pi;
    	double th  = Math.toDegrees(Math.acos(cx));    
    	return elas(eb,th)*jac;
    }
    
    public double elasradq(double eb, double q2, double t, double wcut) {
    	double cx = 1-q2/(2*eb*(eb-q2/2/rmp));
//    	double ep = eb-q2/2/rmp;
//      float jac = pi*(1+q2/(2*ep*0.938))/(ep*eb);
    	double epr=eb/(1+eb*(1-cx));
    	double jac=epr*epr/pi;
    	double th=Math.toDegrees(Math.acos(cx));    
    	return radcor(eb,th,0,t,0.,0.,wcut,0)*jac;
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
    
    	getFF(q2,param);
    
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
    

    public void getFF(double q2, int opt) {
    	
    	// gep,gmn form factor parameterizations 1=dipole ff 2=Bosted ff 3=Brash ff 7=A1-Mainz
   	
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
    	case 3: //Brash parameterization of Hall A GeP/GmP measurement: hep-ex\0111038 PRC 65, 051001
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
//      	gep = divdif(spige,spix,1000,q2,2)
//      	gmp = 2.7928*divdif(spigm,spix,1000,q2,2)
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
    
    public double radcor(double es, double theta_d, int izn, double t1, double t2, double t3, double wcut, int rc) {
    
//  	Radiated elastic cross section from equation II.6 of Mo and Tsai, Rev. Mod. Phys. 41, 205-235 (1969).

//  	es: incident electron energy (GeV)
//  	theta_d: scattered electron angle (deg)
//  	izn: index for target nucleus (1,2,3=H,C,Ta)
//  	t1: target thickness (r.l.)
//  	t2: outgoing e- path length (r.l.)
//		t3: windows (b factor included)
//  	wcut: used to calculate maximum radiated photon energy. wcut should be >  W resolution (GeV).

//  	calls function elas defined above.
    
    	double mp,wc,cst1,eel,epr,epcut,gamma4,beta4,e1,e3,e4,eta,qs,theta;
    	double delta,delta_t1,delta_t2,delta_t3;
    	double radcor1,sigunp;
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
    	e_prime = eel;
    	eta	   = es/eel;
    	q2	   = 2.*es*eel*cst1;
    
    	deltac[0] = 28./9.-13./6.*Math.log(q2/me2);
    	arg=2*Math.log(e1/delta)-3*Math.log(eta);
    	deltac[1] = (Math.log(q2/me2) - 1.)*arg; 
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


    	delta_t1 = -0.5*bfunc(znuc)*t1*(Math.log(es/eta/eta/delta));
    	delta_t2 = -0.5*bfunc(znuc)*t2*(Math.log(eel/delta));
    	delta_t3 = -                t3*(Math.log(eel/delta));
    	delta_t  = delta_t1+delta_t2+delta_t3; // Straggling correction 
    
    	radcor1  = 1. + del_mo+delta_t;
    	radcor   = Math.exp(del_mo + delta_t);
    	
    	xsect_raw = elas(es,theta_d);    	
    	sigunp = izn==1 && rc==0 ? xsect_raw : 1.0;
    	xsect_rad = xsect_raw*radcor;

    	return radcor*sigunp;
    }

    public double radtail(double es_4, double ep_4, double thetae_4, double csthk_4, double phik_4, double egam_4, double d3sig_4) {
        	     
    // Calculates the exact one photon radiative cross section for elastic ep scattering. 
    // Value returned is integrand of equation B-4 in Mo and Tsai, Rev. Mod. Phys. 41, 205-235 (1969)

    // inputs: es = incident electron energy
    //     thetae = electron scattering angle
    //         ep = scattered electron energy
    //     thetak = photon angle with respect to the q vector

    // output: distributions of various quantities calculated on the assumption that both the proton 
    // and the scattered electron are detected, and that the interaction proceeded through a resonance.
    // These can be compared to distributions obtained from AO in which the interaction was actually 
    // through resonance production.

    // The calculations are made for a fixed incident electron energy and a series of scattered electron energies.

    // Units are GeV
    	
    	double es,ep,thetae,cthe,csths,snths,csthp,snthp,thetas,thetap;
    	double pvec, svec;
    	double csthk,snthk,egam;
        double qvec,q0,q2,sp,u2;
        double trm1,trm2,trm3,trm4,trm5,trm6,trm7,trm8,trm9,trm10;
        double trmf,trmg,Ge,Gm,tau,F0,G0,fac1,d3sig;
        double alpha, me,me2,mp,mp2;
        double w2,w;
        double q2e;
        double phik,sdotk,pdotk;

        es	= es_4;
        ep	= ep_4;
        
        cthe  = Math.cos(Math.toRadians(thetae_4));
        csthk = csthk_4;
        phik  = phik_4;

        me	= rme; me2=me*me;
        mp	= rmp; mp2=mp*mp;
        alpha	= 1/137.01;
        d3sig	= 0.;

        pvec = ep - 0.5*me2/ep;
        svec = es - 0.5*me2/es;
        q2e	 = -2*me2 + 2*es*ep - 2*svec*pvec*cthe;
        w2	 = mp2 + 2*mp*(es-ep) - q2e;
        w	 = Math.sqrt(w2);
        	      
        if (w<mp) return 0.0;

        qvec = Math.sqrt(svec*svec + pvec*pvec - 2*svec*pvec*cthe);
        q0	 = es-ep;

        sp	= es*ep - svec*pvec*cthe;
        u2	= Math.pow((es-ep+mp),2) - qvec*qvec;
        csths  = (svec*svec - pvec*svec*cthe)/svec/qvec;
        csthp  = (svec*pvec*cthe - pvec*pvec)/pvec/qvec;
        thetas = Math.acos(csths);
        thetap = Math.acos(csthp);
        snths  = Math.sin(thetas);
        snthp  = Math.sin(thetap);
        snthk  = Math.sqrt(1-csthk*csthk);

        egam   = 0.5*(u2-mp2)/(es - ep + mp - qvec*csthk);
        if (egam<0.)  egam=0.1E-6;
        sdotk  = es*egam-svec*egam*(csthk*csths+snthk*snths*Math.cos(phik));
        pdotk  = ep*egam-pvec*egam*(csthk*csthp+snthk*snthp*Math.cos(phik));

        q2 = 2*me2 - 2*ep*es + 2*pvec*svec*cthe - 2*egam*(es-ep) + 2*egam*qvec*csthk;

        trm1	= -Math.pow(me/pdotk,2)*(2*es*(ep+egam)+q2/2);
        trm2	= -Math.pow(me/sdotk,2)*(2*ep*(es-egam)+q2/2);
        trm3	= -2;
        trm4	= 2/sdotk/pdotk;
        trm4	= trm4*(me2*(sp-egam*egam)+sp*(2*es*ep-sp+egam*(es-ep)));
        trm5	= (1./pdotk)*(2.*(es*ep+es*egam+ep*ep)+q2/2-sp-me2);
        trm6	=-(1./sdotk)*(2.*(es*ep-ep*egam+es*es)+q2/2-sp-me2);
        
        trmf	= trm1+trm2+trm3+trm4+trm5+trm6;
        
        trm7	= me2*(2*me2+q2)*(1./pdotk/pdotk + 1./sdotk/sdotk);
        trm8	= 4.;
        trm9	= 4.*sp*(sp - 2.*me2)/pdotk/sdotk;
        trm10	= (2.*sp + 2*me2 - q2)*(1./pdotk - 1./sdotk);
        
        trmg	= trm7+trm8+trm9+trm10;
        
        Ge		= 1./Math.pow((1.-(q2/0.71)),2);
        Gm	 	= 2.793*Ge;
        tau	 	= -q2/4./mp2;
        F0		= 4.*(Ge*Ge + tau*Gm*Gm)/(1+tau);
        G0		= -q2*Gm*Gm;
        fac1	= Math.pow(alpha,3)/Math.pow(2.*pi,2)*(ep/es)*egam/mp/2./(q2*q2);
        fac1	= fac1/(q0 + mp - qvec*csthk);

        d3sig	= fac1*(mp2*F0*trmf + G0*trmg);
        d3sig_4 = d3sig*389.37966; // convert to ub
        egam_4  = egam;

        return d3sig_4;
   	
    }
    
    public void test1() {
    	
    	setFF(2);
    	
    	double eb = 7.546; 
    	
    	for (int i=1; i<30; i++) {
    		double theta = (i-1)*1+6.5;
    		System.out.println(theta+" "+1e3*elas(eb, theta)+" "+radcor(eb,theta,1,0.0058,0.0058,0.0,0.0617,2));
    	}    	
    }
    
    public void tsai_table1() {
    	
    	Formatter fmt = new Formatter();
    	DecimalFormat df = new DecimalFormat("#.####");
    	
    	setFF(1);

        double[] eb = {17.314, 15.999, 14.649, 13.329, 11.999, 10.723, 6.032, 2.201, 2.206};
    	double[] angle = {35.1, 19.7, 18.8, 17.6, 16.082, 14, 17.186, 38.601, 15.999};
    	
    	fmt.format("%10s %10s %8s %8s %10s %10s %10s\n","Ebeam","Angle","Eelec","-q2","Z0","Z1","Z2");
    	for (int i=0; i<eb.length; i++) {
    		radcor(eb[i],angle[i],0,0.0058,0.0058,0.0,0,1); //wcut=0.1117 for W=1.05 GeV, wcut=0 for Mo&Tsai Table 1   		
    		fmt.format("%10s %10s %8s %8s %10s %10s %10s\n",
    				   eb[i],angle[i],df.format(e_prime),df.format(q2),df.format(z0sum),df.format(z1sum),df.format(z2sum));
    	}
    	System.out.println(fmt);
    }
    
    public void rc(double ebeam, double thmin, double thmax, double bw) {
    	
    	Formatter fmt = new Formatter();
    	DecimalFormat df = new DecimalFormat("#.###");
    	DecimalFormat dfx = new DecimalFormat("####.###");
    	
    	String format = "%8s %5s %5s %5s %9s %9s %7s %7s %7s %7s\n";
    	
    	setFF(2);

    	fmt.format(format,"Ebeam","Angle","Eelec","-q2","xsraw(nb)","xsrad(nb)","rc_int","rc_ext","rc","bcc");

    	for (double theta=thmin/bw; theta<thmax/bw; theta++) {
    		double the=theta*bw;
    		radcor(ebeam,the,0,0.0058,0.0058,0.0,0,1);
    		fmt.format(format,
    		ebeam,df.format(the),df.format(e_prime),
    		            df.format(q2),
    		            dfx.format(xsect_raw*1e3),
    		            dfx.format(xsect_rad*1e3),
    		            df.format(del_mo),
    		            df.format(delta_t),
    		            df.format(radcor),
    		            df.format(getBCC(ebeam,the,bw,10000)));
    	}
    	System.out.println(fmt);
    }
    
    public void test3() {
    	double[] x = {0.1,0.2,0.3,0.8,1.0,2.0,3.0,10,20,30,40};
    	for (int i=0; i<x.length; i++) System.out.println(x[i]+" "+spence(x[i])+" "+spence(-x[i]));
    }
    
    public static void main(String[] args) {    	
    	MoTsai elib = new MoTsai();
//    	elib.tsai_table1(); 
    	elib.rc(7.546,5,40,2);
    }

}
