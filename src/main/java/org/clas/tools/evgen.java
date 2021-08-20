package org.clas.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.ParticleGenerator;
import org.jlab.clas.physics.PhysicsEvent;

public class evgen {
	
	List<String> out = new ArrayList<String>();
	int   nc;
	
	GenEvent ge ;
	
	public evgen() {
		
	}
	
	public void genEvents(String str, int val) {
		if(str=="pim_pi0") ge = new gen_pim_pi0();
		if(str=="pim_gam") ge = new gen_pim_gam();
		if(str=="e_pim")   ge = new gen_e_pim();
		if(str=="e_mum")   ge = new gen_e_mum();
		if(str=="e_all")   ge = new gen_e_all();
		writeFile(str,val);
	}
		
	abstract class GenEvent {
		public abstract List<String> getev();		
	}
	
	public class gen_pim_g_g extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(3,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(getParticle(  22,4,4,15,15, 59, 61),2));
			out.add(getString(getParticle(  22,4,4,10,20, 50, 70),3));
			return out;			
		}		
	}
	
	public class gen_pim_gam extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			float vtx = (float) (-5.5+Math.random()*5);
			out.add(getString(getParticle(-211,2,2,18,18,240,260,vtx),1));
			out.add(getString(getParticle(  22,0.1,4,5.5,25.5,45,75,vtx),2));
			return out;			
		}		
	}	
	
	public class gen_pim_pi0 extends GenEvent {		
		public List<String> getev() {
			out.clear();
			List<Particle> plist = getDecay(111,4,8,14,16,-5,5);
			out.add(getHeader(3,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(plist.get(0),2));
			out.add(getString(plist.get(1),3));
			return out;			
		}		
	}
	
	public class gen_e_pim extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(  11,7,9,10,25,240,260),1));
			out.add(getString(getParticle(-211,2,5,15,30, 59, 61),2));
			return out;			
		}		
	}
	
	public class gen_e_pip extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(  11,7,9,10,25,240,260),1));
			out.add(getString(getParticle(+211,2,5,15,30, 59, 61),2));
			return out;			
		}		
	}
	
	public class gen_e_all extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(4,nc));
			out.add(getString(getParticle(   11,0.80,2.0, 5,25,250,270),1));
			out.add(getString(getParticle(+2212,0.55,1.4,28,30, 12, 52),2));
			out.add(getString(getParticle( +211,0.40,1.4,28,30, 72,112),3));
			out.add(getString(getParticle( -211,0.40,1.4,28,30,132,172),4));
			return out;			
		}		
	}
	
	public class gen_e_mum extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(  11,7,9,10,25,240,260),1));
			out.add(getString(getParticle(  13,2,5,15,30, 59, 61),2));
			return out;			
		}		
	}
	
	public String getHeader(int np, int n) {
		return np+" "+n+" 0. 0. 0. 0. 0. 0. 0. 0.";		
	}
	
	public Particle getParticle(double... pargs) {
		Particle p = new Particle();
		int   pid = (int) pargs[0];
		double pm = pargs[1]+Math.random()*(pargs[2]-pargs[1]);
	    double th = Math.toRadians((pargs[3]+Math.random()*(pargs[4]-pargs[3])));
	    double ph = Math.toRadians((pargs[5]+Math.random()*(pargs[6]-pargs[5])));
	    double px = pm*Math.sin(th)*Math.cos(ph);
	    double py = pm*Math.sin(th)*Math.sin(ph);
	    double pz = pm*Math.cos(th);		
	    p.setVector(pid,px,py,pz,0,0,pargs.length>7?pargs[7]:0);
	    return p;
	}
	
	public List<Particle> getDecay(double... pargs) {		
		List<Particle> plist = new ArrayList<Particle>();		
		TwoBodyDecay decay = new TwoBodyDecay();
    	ParticleGenerator pg = new ParticleGenerator((int)pargs[0]);    	
    	pg.setRange(pargs[1],pargs[2],pargs[3],pargs[4],pargs[5],pargs[6]); 
    	PhysicsEvent event = new PhysicsEvent(); event.clear();
    	event.addParticle(pg.getParticle());
    	decay.decayParticles(event);
    	plist.add(event.getGeneratedParticle(0)); 
    	plist.add(event.getGeneratedParticle(1));    	
		return plist;
	}
	
	public String getString(Particle p, int np) {
		return np+" "+p.charge()+" 1 "+p.pid()+" 0 0 "+String.format("%.4f",p.px())+" "
                                                      +String.format("%.4f",p.py())+" "
                                                      +String.format("%.4f",p.pz())+" "
                                                      +"1 "+"1"+" "
                                                      +String.format("%.4f",p.vx())+" "
                                                      +String.format("%.4f",p.vy())+" "
                                                      +String.format("%.4f",p.vz())+" "
                                                      +String.format("%.4f",-3.);			    
	}
	
	public void writeFile(String str, int val) {
		
		List<String> list = new ArrayList<String>();
		
		try { 
			File outputFile = new File("/Users/colesmith/CLAS12/sim/photon/"+str+".lund");
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("evgen writing to "+outputFile.getName());
			
			for (nc=1; nc<val+1; nc++) {				
				list=ge.getev();
				for (String line : list) {
					if (line!=null) {
//						System.out.println(line);
					    outputBw.write(line);
					    outputBw.newLine();
					}					
				}				
			}

			outputBw.close();
			outputFw.close();
		}
		catch(IOException ex) {
			System.out.println("Error writing file '" );                   
			ex.printStackTrace();
		}
		System.out.println("Finished");

	}
	
	public static void main(String[] args) {
		evgen gen = new evgen();
//		gen.genEvents("e_pim",10000);
//		gen.genEvents("e_mum", 10000);
//		gen.genEvents("e_all", 50000);
//		gen.genEvents("pim_g_g",1000);
//		gen.genEvents("pim_pi0",10000);
		gen.genEvents("pim_gam",50000);
	}

}
