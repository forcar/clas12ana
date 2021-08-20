package org.clas.tools;

public class KinLib {
	
	float    mp = 0.93828f;
	float melec = 0.000511f;
	float mpip  = 0.13957f;
	float mpim  = 0.13957f;
	float mpi0  = 0.13485f;
	float   pi  = 3.1415926f;
		
	public KinLib() {
		
	}
	
	public float ep_from_w(float es, float ethe, float wcut) {
		
		if (wcut==0) wcut=mp;
		float cst1 = (float) (1-Math.cos(Math.toRadians((float)ethe)));
		float  eel = es/(1+es/mp*cst1);
		return (es+(mp*mp-wcut*wcut)/2/mp)*eel/es;
		
	}

}
