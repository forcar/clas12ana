package org.clas.tools;

public class KinLib {
	
	public double    mp = 0.93828f;
	public double melec = 0.000511f;
	public double mpip  = 0.13957f;
	public double mpim  = 0.13957f;
	public double mpi0  = 0.13485f;
	public double   pi  = 3.1415926f;
	public double d2r   = pi/180;
	
	double e0,w0,p1,e1,p2,e2,th2;
		
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
        double   a = (ws  + w12)*0.5f;

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

}
