package org.clas.tools;

public class KinLib {
	
	public double    mp = 0.93828f;
	public double  m3He = 2.8094132f;
	public double melec = 0.000511f;
	public double  mpip = 0.13957f;
	public double  mpim = 0.13957f;
	public double  mpi0 = 0.13485f;
	public double    pi = 3.1415926f;
	public double   d2r = pi/180;
	
	double p0,e0,w0,p1,e1,p2,e2,th2,p2a,p2b,e2a,e2b,p3a,p3b,e3a,e3b,th3a,th3b,ph3a,ph3b,mma,mmb;
	double[] p, e, ct, st;
		
	public KinLib() {
		
	}
	
	public double ep_from_w(double es, double ethe, double wcut) {				
		if (wcut==0) wcut=mp;
		double cst1 = 1-Math.cos(Math.toRadians(ethe));
		double  eel = es/(1+es/mp*cst1);
		return (es+(mp*mp-wcut*wcut)/2/mp)*eel/es;		
	}
	
	public double[] getQW(double es, double epmin, double thmin) {
        double    nu = es - epmin;
        double q2min = 2*es*epmin*(1-Math.cos(thmin*d2r));
        double  wmax = -q2min + mp*mp + 2*mp*nu;
        double out[] = {q2min,Math.sqrt(wmax)};
        return out;        
	}
	
	public double p1vsth1(double x, double eb) {
		double rm[] = {melec,mp,melec,mp};
		twobod(eb,x,rm);
		return p1;
	}
	
	public double p2vsth2(double x, double eb) {
		double rm[] = {melec,mp,mp,melec};
		twobod(eb,x,rm);
		return p1;		
	}
	
	public double p2vsth1(double x, double eb) {
		double rm[] = {melec,mp,melec,mp};
		twobod(eb,x,rm);
		return p2;	
	}
	
	public double th1vsp1(double x, double eb) {
		return (Math.acos(1-mp*(eb-x)/eb/x)/d2r);
	}
	
	public double th2vsp2(double x, double eb) {
		double w0 = eb + mp;
		double e2 = Math.sqrt(x*x + mp*mp);
		double p1 = w0 - e2;
		double th1 = th1vsp1(p1,eb)*d2r;
		double p2 = p2vsp1(p1,eb);
		double c2 = Math.min(1.0, Math.max(-1.0, (eb-p1*Math.cos(th1))/p2));
		double s2 = Math.sqrt(1-c2*c2);
		double th2 = Math.atan(s2/c2)/d2r;
		if(th2<0) th2 = th2 + 180;
		return th2;
	}
	
	public double p2vsp1(double x, double eb) {
		double e2 = eb + mp - Math.sqrt(x*x+melec*melec);
		double p2 = Math.sqrt(e2*e2 - mp*mp);
		return p2;
	}
	
	public double p1vsp2(double x, double eb) {
		double e1 = eb + mp - Math.sqrt(x*x+mp*mp);
		double p1 = Math.sqrt(e1*e1 - melec*melec);
		return p1;
	}	
	
	public double th12vsth1(double x, double eb) {
		double rm[] = {melec,mp,melec,mp};
		twobod(eb,x,rm);
		return th2+x;
	}
	
	public double th2vsth1(double x, double eb) {
		double rm[] = {melec,mp,melec,mp};
		twobod(eb,x,rm);
		return th2;
	}
	
	public double th1vsth2(double x, double eb) {
		double rm[] = {melec,mp,mp,melec};
		twobod(eb,x,rm);
		return th2;
	}
	
	public double th2vsp1(double x, double eb) {
		p2 = p2vsp1(x,eb);
		return th2vsp2(p2,eb);
	}
		
	public double th1vsp2(double x, double eb) {
		p1 = p1vsp2(x,eb);
		return th1vsp1(p1,eb);
	}

	public void twobod(double p0, double th1, double[] rm) {
    
    	double rm0 = rm[0];
    	double rmn = rm[1];
    	double rm1 = rm[2];
    	double rm2 = rm[3];
    
        double theta1 = th1*d2r;
    
        double e0 = Math.sqrt(p0*p0 + rm0*rm0);
        double w0 = e0 + rmn;

        double ws = w0*w0 - p0*p0;
        double w12 = (rm1 - rm2)*(rm1 + rm2);

        double cs0 = Math.cos(theta1);

        double   d = w0*w0 - (p0*cs0)*(p0*cs0);
        double   a = (ws  + w12)*0.5;

        e1 = (a*w0 + p0*cs0*Math.sqrt(a*a - rm1*rm1*d))/d;
        p1 = Math.sqrt(e1*e1 - rm1*rm1);

        e2 = w0 - e1;
        p2 = Math.sqrt(e2*e2 - rm2*rm2);

        double c2 = Math.min(1.0,Math.max(-1.0,(p0 - p1*cs0)/p2));
        double s2 = Math.sqrt(1 - c2*c2);

        th2 = Math.atan(s2/c2)/d2r;
        
        if(th2<0) th2=th2 + 180;

        return;
    
	}
		
	public void test3b() {
		double rm[] = {mpip,m3He,mp,mp,mp};
		System.out.print("p1    p2a   p2b   p3a   p3b   th3a   th3b    ph3a   ph3b  mma   mmb\n");
		for (double p1=0; p1<0.88; p1+=0.01) {
			System.out.println(" ");
//			threebod(2,3,0.35+mpip,p1,40,108,0,180,rm);
			threebod(2,3,0.35+mpip,p1,60,82,1.01,179,rm);
			System.out.printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",p1,p2a,p2b,p3a,p3b,th3a,th3b,ph3a,ph3b,mma,mmb);
		}
	}
	
	public void threebod(int i, int j, double e0, double p1, double th1, double th2, double ph1, double ph2, double[] rm) {
		
		// 0 -> 1 + 2 + 3 (decay) or 0 + N -> 1 + 2 + 3 (collision)
		// P0 and E0 are mom. and energy for particle 0.
		// RM0 and RMN are the rest masses for 2-body initial state.
		// 3-body final state is assumed. Kinematics are determined if 5 quantities are specified.  
		// This method calculates P2 if P1 and cos12 are specified, (cos12 involves 4 quantities).  
		// **CAUTION** there are in general two roots for P2, however sometimes only 1 root is physical.  
		// This method calculates both roots. If root is not physical we set P2=E2=0
		
		double p2,p3,e2,e3,ct3,sp3;
		
        double rm0 = rm[0];
        double rmn = rm[1];
        double rm1 = rm[i];
        double rm2 = rm[j];
        double rm3 = rm[4];
        
        double[]  p2r = {0,0},  e2r = {0,0};
        double[]  p3r = {0,0},  e3r = {0,0};
        double[] ct3r = {0,0}, sp3r = {0,0};
        double[]  mmr = {0,0};

        int[] sign = {1,-1};

		double  et = e0 + rmn;
		double  p0 = Math.sqrt(e0*e0 - rm0*rm0);
		double  e1 = Math.sqrt(p1*p1 + rm1*rm1);
		double ct1 = Math.cos(Math.toRadians(th1));
		double ct2 = Math.cos(Math.toRadians(th2));
		double st1 = Math.sin(Math.toRadians(th1));
		double st2 = Math.sin(Math.toRadians(th2));
		double sp2 = Math.sin(Math.toRadians(ph2));
		double c12 = ct1*ct2 + st1*st2 * Math.cos(Math.toRadians(ph1)-Math.toRadians(ph2));

		double a = 0.5*(et*et - p0*p0 + rm1*rm1 + rm2*rm2 - rm3*rm3) - et*e1 + p0*p1*ct1;
		double b = e1 - et;
		double c = p0*ct1 - p1*c12;
		double bc = b*b - c*c;
		double d  = a*a - rm2*rm2*bc;

		if(d<-0.02) return;
		      
		d = Math.sqrt(Math.abs(d));     

		double ac = a*c;
		double bd = b*d;

		double check = 0.5*(et*et - p0*p0 + rm1*rm1 + rm2*rm2) - et*e1 + p0*p1*ct1;

		for (int k=0; k<2; k++) {

		      p2 = Math.abs((ac + sign[k]*bd)/bc);
		      e2 = Math.sqrt(p2*p2 + rm2*rm2);

		      double chk = check - et*e2 + p0*p2*ct2 + e1*e2 - p1*p2*c12;
		      double  ck = Math.sqrt(Math.abs(2*chk));

		      double em = Math.pow((et - e1 - e2),2);
		      double pm = p0*p0 + p1*p1 + p2*p2 - 2*p0*(p1*ct1+ p2*ct2) + 2*p1*p2*c12;
		      double mm = Math.abs(Math.sqrt(Math.abs(em-pm)) - rm3);
		      
		      e3 = et - e1 - e2;
		      p3 = Math.sqrt(e3*e3-rm3*rm3);
		      ct3 = (p0 - p1*ct1 - p2*ct2)/p3; 
		      System.out.println(e3+" "+rm3+" "+ct3);
		      sp3 = -p2*st2*sp2/p3/Math.sqrt(1-ct3*ct3);

		      p2r[k] = mm<=1 ?  p2: 0; 
		      e2r[k] = mm<=1 ?  e2: 0; 
		      p3r[k] = mm<=1 ?  p3: 0; 
		      e3r[k] = mm<=1 ?  e3: 0; 
		     ct3r[k] = mm<=1 ? ct3: 0; 
		     sp3r[k] = mm<=1 ? sp3: 0; 
		      mmr[k] = mm;

		}

		p2a = p2r[0]; p3a = p3r[0]; th3a = Math.toDegrees(Math.acos(ct3r[0])); 
		p2b = p2r[1]; p3b = p3r[1]; th3b = Math.toDegrees(Math.acos(ct3r[1]));
		e2a = e2r[0]; e3a = e3r[0]; ph3a = Math.toDegrees(Math.asin(sp3r[0]));
		e2b = e2r[1]; e3b = e3r[1]; ph3b = Math.toDegrees(Math.asin(sp3r[1]));
		mma = mmr[0]; mmb = mmr[1];
 		      
        return;
		
	}
	
    public static void main(String[] args) { 
    	
    	KinLib klib = new KinLib();
    	klib.test3b();

    }

}
