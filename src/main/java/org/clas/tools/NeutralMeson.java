package org.clas.tools;

import java.util.ArrayList;
import java.util.List;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.geom.prim.Point3D;
import org.jlab.utils.groups.IndexedList;

public class NeutralMeson {

	Event    ev = null;
	
	public float mass;
	public float e;
	public float mom;
	public float the;
	public float phi;
	public float opa;
	public float X;
	public float ggp;
	
	public  LorentzVector VPI0;
	public  LorentzVector VG1;
	public  LorentzVector VG2;
	
	private boolean tag = false;
	private int mcsec = 2;
	private double thresh=0.2;
	private double mpi0 = 0.1349764;
		
	public List<Particle>  plist = new ArrayList<Particle>(); 
	public IndexedList<List<Particle>> ilist = new IndexedList<List<Particle>>(1);
	public List<Float> out = new ArrayList<Float>();	
	
	public NeutralMeson(boolean val) {
		this.tag = val;
	}
	
	public NeutralMeson(List<Particle> val) {
		this.plist = val;
	}
	
	public NeutralMeson(IndexedList<List<Particle>> val) {
		this.ilist = val;
	}
	
	public void addEvent(Event val) {
		ev = new Event();
		ev = val;
	}
	
	public void addPhoton(Particle p) {
		plist.add(p);
	}
	
	public void setThresh(double val) {
		thresh=val;
	}
	
	public void setPhotonSector(int val) {
		mcsec = val;
	}
	
	public Particle getPhoton(int n) {
		return plist.get(n);    		
	}
	
	public int getPINDEX(int n) {
		return (int) getPhoton(n).getProperty("pindex");
	}
	
	public int getPhotonSector(int n) {
		return (int) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("sector");
	}
	
	public int getPhotonLayer(int n) {
		return (int) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("layer");
	}
	
	public float getPhotonEnergy(int n) {
		return (float) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("energy");
	}
	
	public float getParticleEnergy(int n) {
		return (float) ev.part.get(getPINDEX(n)).p();
	}
	
	public float getBeta(int n) {
		return (float) getPhoton(n).getProperty("beta");    		
	}
	
	public float getEnergy(int n) {
		return (float) getPhoton(n).p();    		
	}  
	
	public Boolean filter(boolean val, int i1, int i2) {
		if (val) return  sameSector(i1,i2);
		return !sameSector(i1,i2);
	}
	
	public boolean sameSector(int i1, int i2) {
		return getPhotonSector(i1)==getPhotonSector(i2);		
	}
	
    public double Vangle(Vector3 v1, Vector3 v2){
        double res = 0;
        double l1 = v1.mag();
        double l2 = v2.mag();
        double prod = v1.dot(v2);
        if( l1 * l2 !=0 && Math.abs(prod)<l1*l2 )res = Math.toDegrees( Math.acos(prod/(l1*l2) ) );
        return res; 
    }
    
    public double getGGphi(int i1, int i2) {
        Vector3 n1 = new Vector3(), n2 = new Vector3();
        n1.copy(plist.get(i1).vector().vect()); n2.copy(plist.get(i2).vector().vect());
        n1.unit(); n2.unit();
        Point3D point1 = new Point3D(n1.x(),n1.y(),n1.z());
        Point3D point2 = new Point3D(n2.x(),n2.y(),n2.z());
    	point1.rotateZ(Math.toRadians(-60*(mcsec-1))); point2.rotateZ(Math.toRadians(-60*(mcsec-1)));
        point1.rotateY(Math.toRadians(-25));         point2.rotateY(Math.toRadians(-25));      
        Vector3 vv1 = new Vector3(point1.x(),point1.y(),point1.z()); 
        Vector3 vv2 = new Vector3(point2.x(),point2.y(),point2.z());
        Vector3 vv12 = vv1.cross(vv2);
        double ggp = Math.toDegrees(Math.atan2(vv12.y(),vv12.x()));
        if(ggp<0) ggp=-ggp;
        ggp=ggp-90;
        if(ggp<0) ggp=ggp+180;   
        return ggp;
    }
	
	public void getMeson(int i1, int i2) {
		Particle p1 = plist.get(i1);
		Particle p2 = plist.get(i2);
		VG1 = new LorentzVector(p1.px(),p1.py(),p1.pz(),p1.p());
		VG2 = new LorentzVector(p2.px(),p2.py(),p2.pz(),p2.p());		
		VPI0 = new LorentzVector(0,0,0,0);
		VPI0.add(VG1); 
		VPI0.add(VG2);
	}
	
	public boolean test(int i1, int i2) {
		if (ev==null) return true;
		return this.mass > 0.08 && filter(tag,i1,i2);
	}
	
	public boolean getPizeroKinematics(Particle p1, Particle p2) {
		
		if (p1.e()<thresh && p2.e()<thresh) return false; 
				
		double e1 = p1.e();
		double e2 = p2.e(); 
		
        Vector3 n1 = new Vector3(); Vector3 n2 = new Vector3();
        n1.copy(p1.vector().vect()); n2.copy( p2.vector().vect());
        n1.unit(); n2.unit();
       
        Particle g1 = new Particle(22,n1.x()*e1,n1.y()*e1,n1.z()*e1);                        
        Particle g2 = new Particle(22,n2.x()*e2,n2.y()*e2,n2.z()*e2);
       
        double cth1 = Math.cos(g1.theta());
        double cth2 = Math.cos(g2.theta());
        double  cth = g1.cosTheta(g2);         
        double    x = (e1-e2)/(e1+e2);
        double tpi2 = 2*mpi0*mpi0/(1-cth)/(1-X*X);
        double cpi0 = (e1*cth1+e2*cth2)/Math.sqrt(e1*e1+e2*e2+2*e1*e2*cth);
        
        g1.combine(g2, +1); double invm = Math.sqrt(g1.mass2());
       
        int n=0;
        out.add(n++,(float) Math.sqrt(tpi2));
        out.add(n++,(float) Math.toDegrees(Math.acos(cpi0)));
        out.add(n++,(float) ((invm-mpi0)/mpi0));
        out.add(n++,(float) Math.toDegrees(Math.acos(cth)));
        out.add(n++,(float) Math.abs(x));
        out.add(n++,(float) Math.sqrt(e1*e2));
        out.add(n++,(float) phi);
        out.add(n++,(float) (float)(Math.toDegrees(p1.theta())-Math.toDegrees(p2.theta())));
        out.add(n++,(float) (float)(Math.toDegrees(p1.phi())-Math.toDegrees(p2.phi())));
        
        return true;
	}
	
	public boolean getMesonKin(int i1, int i2) {
		this.mass = (float)VPI0.mass();
		this.e    = (float)VPI0.e();
		this.mom  = (float)VPI0.p();
		this.the  = (float)Math.toDegrees(VPI0.theta());
		this.phi  = (float)Math.toDegrees(VPI0.phi());
		this.opa  = (float)Vangle(VG1.vect(),VG2.vect());
		this.X    = (float)((VG1.e()-VG2.e())/(VG1.e()+VG2.e()));
		this.ggp  = (float) getGGphi(i1,i2);
		return test(i1,i2);			
	}
	
	public String toString(int n) {
		return getPINDEX(n)+" "+getPhotonSector(n)+" "+getPhotonLayer(n)+" "+getPhotonEnergy(n)+" "+getParticleEnergy(n);
	}
	    	    	
}

