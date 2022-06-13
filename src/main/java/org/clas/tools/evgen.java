package org.clas.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.ParticleGenerator;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;

public class evgen {
	
	GenEvent ge;
	List<String> out = new ArrayList<String>();
	int   nc, nev;
	boolean debug=false;
	
	public evgen(boolean val) {
		debug = val;
	}
	
	public void genEvents(String str, int nev) {
		this.nev = nev;
		if(str=="pim_g_g") ge = new gen_pim_g_g();
		if(str=="pim_g")   ge = new gen_pim_g();
		if(str=="pim_mu_mu") ge = new gen_pim_mu_mu();
		if(str=="pim_mup") ge = new gen_pim_mup();
		if(str=="pim_mum") ge = new gen_pim_mum();
		if(str=="pim_pi0") ge = new gen_pim_pi0();
		if(str=="pim_gam") ge = new gen_pim_gam(); 
		if(str=="pim_n")   ge = new gen_pim_n();
		if(str=="e_pim")   ge = new gen_e_pim();
		if(str=="e_mum")   ge = new gen_e_mum();
		if(str=="e_all")   ge = new gen_e_all();
		writeFile(ge,str,nev);
	}
			
	abstract class GenEvent {		
		public abstract List<String> getev();		
	}
	
	public class gen_pim_g_g extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(3,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(getParticle(  22,4,4,15,15, 60, 60),2));
			out.add(getString(getParticle(  22,4,4,10,20, 50, 70),3));
			return out;			
		}		
	}
	
	public class gen_pim_g extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(getParticle(  22,4,4,10,20, 50, 70),2));
			return out;			
		}		
	}
	
	public class gen_pim_mu_mu extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(3,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(getParticle(  13,4,4, 4,36, 45, 95),2));
			out.add(getString(getParticle(  13,4,4, 4,36, 45, 95),3));
			return out;			
		}		
	}
	
	public class gen_pim_mup extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(getParticle( -13,0.5,3.5, 5,25,  10, 70),2));
			return out;			
		}		
	}	
	
	public class gen_pim_mum extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(getParticle(  13,0.7,3.7, 4,40, 50, 110),2));
			return out;			
		}		
	}
	
	public class gen_pim_gam extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			float vtx = 0;
//			float vtx = (float) (-5.5+Math.random()*5);
			out.add(getString(getParticle(-211,2,2,18,18,240,260,vtx),1));
			out.add(getString(getParticle(  22, 0.05, 4, 5.5, 25.5, 45, 75,vtx),2));
			return out;			
		}		
	}	
	public class gen_pim_n extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
//			float vtx = (float) (-5.5+Math.random()*5);
			float vtx = 0;
			out.add(getString(getParticle(-211,2,2,18,18,240,260,vtx),1));
			out.add(getString(getParticle(2112,0.3,6,5.5,25.5,45,75,vtx),2));
			return out;			
		}		
	}
	
	public class gen_pim_pi0 extends GenEvent {	
		Decay decay = new Decay(111,1,8,24,26,55,65);
		public List<String> getev() {
			out.clear();
			List<Particle> plist = decay.doDecay();
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
	
	public class Decay {
		
		ParticleGenerator pg;
		TwoBodyDecay   decay = new TwoBodyDecay(); ;
		PhysicsEvent   event = new PhysicsEvent();
		List<Particle> plist = new ArrayList<Particle>();
		
    	H2F h;
		
		public Decay(double... pargs) { 
			initDecay(pargs);
			if (debug) h = new H2F("decay",50, -10, 10, 50,-10,10);			
		}
		
		public void initDecay(double... pargs) {			 
	    	pg    = new ParticleGenerator((int)pargs[0]);    	
	    	pg.setRange(pargs[1],pargs[2],pargs[3],pargs[4],pargs[5],pargs[6]); 
		}
		
		public List<Particle> doDecay() {	
			plist.clear(); event.clear();	    	
	    	event.addParticle(pg.getParticle());
	    	decay.decayParticles(event);
	    	plist.add(event.getGeneratedParticle(0)); 
	    	plist.add(event.getGeneratedParticle(1));  
	    	if(debug) store(plist); 	    	
			return plist;
		}
		
		public void store(List<Particle> list) {
    		Particle p1 = plist.get(0), p2 = plist.get(1);
        	float dthe = (float)(Math.toDegrees(p1.theta())-Math.toDegrees(p2.theta()));
        	float dphi = (float)(Math.toDegrees(p1.phi())  -Math.toDegrees(p2.phi()));
        	h.fill(dphi, dthe);	
        	if(nc==nev) display();
		}
		
		public void display() {
	        JFrame frame = new JFrame("Pizero");
	        frame.setSize(800,800);    	        
	        EmbeddedCanvas c = new EmbeddedCanvas(); 
	    	c.cd(0); c.draw(h);     
	        frame.add(c);
	        frame.setLocationRelativeTo(null);
	        frame.setVisible(true);				
		}
	}
	
	public void decaytest(int val) {
		nev = val;   	
		Decay decay = new Decay(111,1,8,24,26,55,65); 
    	for (nc=1; nc<val+1; nc++) decay.doDecay();		
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
	
	public void writeFile(GenEvent ge, String str, int val) {
		
		List<String> list = new ArrayList<String>();
		
		try { 
			File outputFile = new File("/Users/colesmith/CLAS12/sim/neutron/"+str+".lund");
			FileWriter     outputFw = new FileWriter(outputFile.getAbsoluteFile());
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
		evgen gen = new evgen(true);
//		gen.genEvents("e_pim",1000);
//		gen.genEvents("e_mum", 10000);
//		gen.genEvents("e_all", 50000);
//		gen.genEvents("pim_g",100000);
//		gen.genEvents("pim_pi0",10000);
//		gen.genEvents("pim_gam",50000);
		gen.genEvents("pim_n",80000);
//		gen.genEvents("pim_mu_mu", 100000);
//		gen.genEvents("pim_mup", 100000);
//		gen.genEvents("pim_mum", 100000);
//		gen.decaytest(10000);
	}

}
