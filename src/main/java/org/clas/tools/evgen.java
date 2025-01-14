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
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataSource;

public class evgen {
	
	GenEvent ge;
	String mcfile, inpath="/Users/colesmith/CLAS12/sim/",outpath=inpath+"gen/lund/";
	HipoDataSource reader = new HipoDataSource();
	List<String> out = new ArrayList<String>();
	int   nc, nev;
	boolean debug=false;
	
	public evgen(boolean val) {
		debug = val;
	}
	
	public void genEvents(int nev, String...str) {
		this.nev = nev; int n=0;
		if(str[n]=="pim_g_g") ge = new gen_pim_g_g();
		if(str[n]=="pim_g")   ge = new gen_pim_g();
		if(str[n]=="mum")     ge = new gen_mum();
		if(str[n]=="pim_mu_mu") ge = new gen_pim_mu_mu();
		if(str[n]=="pim_mup") ge = new gen_pim_mup();
		if(str[n]=="pim_mum") ge = new gen_pim_mum();
		if(str[n]=="pim_pi0") ge = new gen_pim_pi0();
		if(str[n]=="pim_gam") ge = new gen_pim_gam(); 
		if(str[n]=="pim_n")   ge = new gen_pim_n();
		if(str[n]=="e_pim")   ge = new gen_e_pim();
		if(str[n]=="e_mum")   ge = new gen_e_mum();
		if(str[n]=="e_g")     ge = new gen_e_g();
		if(str[n]=="e_all")   ge = new gen_e_all();
		if(str[n]=="file")   {ge = new gen_from_file(); init_gen_from_file(str); n=2;} 
		writeFile(ge,str[n],nev);
	}
			
	abstract class GenEvent {		
		public abstract List<String> getev();		
	}
	
	public void init_gen_from_file(String[] str) {
		reader.open(inpath+str[1]+"/"+str[2]+".hipo");
		reader.getNextEvent();
		reader.getNextEvent();
	}
	
	public class gen_from_file extends GenEvent {
		public List<String> getev() {
			out.clear();
			out.add(getHeader(1,nc));	
			out.add(getString(getParticle(reader.getNextEvent()),1));
			return out;
		}
	}
	
	public class gen_pim_g_g extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(3,nc));
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1));
			out.add(getString(getParticle(  22,4,4,20,20, 60, 60),2));
			out.add(getString(getParticle(  22,1,4, 8,32, 55, 70),3));
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
			out.add(getString(getParticle(-211,2,2,18,18,240,260),1)); //PDG +pi=211 -pi=-211
			out.add(getString(getParticle(  13,1,9, 4,36, 35, 85),2)); //PDG -mu=13  +mu=-13
			out.add(getString(getParticle( -13,1,9, 4,36, 95, 145),3));
			return out;			
		}		
	}
	
	public class gen_mum extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(1,nc));			
			out.add(getString(getParticle(13,2,2,11,35,30,90),1)); //PDG -mu=13  
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
	
	public class gen_e_g extends GenEvent {	//uses coatjava validation kinematics	
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc)); 
//			out.add(getString(getParticle(11,2.0,10,6,33,220,270,0,10),1));
//			out.add(getString(getParticle(22,0.5, 9,6,33, 30, 90),2));
			out.add(getString(getParticle(11,1.0,9.0,15,35,-10,10),1));
			out.add(getString(getParticle(22,0.5,4.5,20,35,110,130),2));
			return out;			
		}		
	}
	
	public class gen_e_all extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(5,nc));
			float vtx = (float) (-5.5+Math.random()*5);
			out.add(getString(getParticle(   11,1.00,9.0, 5,25,-10, 10,vtx),1));
			out.add(getString(getParticle(+2212,2.00,5.0,12,30, 50, 70,vtx),2));
			out.add(getString(getParticle( +211,1.00,4.0,12,30,110,130,vtx),3));
			out.add(getString(getParticle(   22,0.50,4.5,10,35,170,190,vtx),4));			
			out.add(getString(getParticle( 2112,0.30,6.0,10,35,230,250,vtx),5));
			return out;			
		}		
	}
	
	public class gen_e_mum extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(11,2,9,12,35,230,290),1));
			out.add(getString(getParticle(13,1,7,11,35, 40,100),2)); //PDG -mu=13  
			return out;			
		}		
	}
	
	public String getHeader(int np, int n) {
		return np+" "+n+" 0. 0. 0. 0. 0. 0. 0. 0.";		
	}
	
	public Particle getParticle(double... pargs) {
		Particle p = new Particle();
	    double rand = Math.random();
		double dphi = pargs.length==9 ? pargs[8]*(1-rand):0;
		pargs[5]=pargs[5]+dphi; pargs[6]=pargs[6]+dphi;
		int   pid = (int) pargs[0];
		double pm = pargs[1]+Math.random()*(pargs[2]-pargs[1]);
		double c2 = Math.cos(Math.toRadians(pargs[3]));
		double c1 = Math.cos(Math.toRadians(pargs[4]));
	    double cth = c1+Math.random()*(c2-c1);
	    double sth = Math.sqrt(1-cth*cth);
	    double ph = Math.toRadians((pargs[5]+Math.random()*(pargs[6]-pargs[5])));
	    double px = pm*sth*Math.cos(ph);
	    double py = pm*sth*Math.sin(ph);
	    double pz = pm*cth;		
	    p.setVector(pid,px,py,pz,0,0,pargs.length==8?pargs[7]:0);
	    return p;
	}
	
	public Particle getParticle(DataEvent ev) {
		Particle p = new Particle();
		if(ev.hasBank("MC::Particle")) {
			DataBank mc = ev.getBank("MC::Particle");
			int  pid = mc.getInt("pid",0);
			float px = mc.getFloat("px",0);
			float py = mc.getFloat("py",0);
			float pz = mc.getFloat("pz",0);
			float vx = mc.getFloat("vx",0);
			float vy = mc.getFloat("vy",0);
			float vz = mc.getFloat("vz",0);
			p.setVector(pid,px,py,pz,vx,vy,vz);
		}
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
			File outputFile = new File(outpath+str+".lund");
			FileWriter     outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("evgen writing to "+outputFile.getName());
			
			for (nc=1; nc<val+1; nc++) {				
				list=ge.getev();
				for (String line : list) {
					if (line!=null) {
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
//		gen.genEvents(1000000, "file","muon","muon_1GeV_376");
//		gen.genEvents(1000000,"mum");
//		gen.genEvents(100000,"e_pim");
//		gen.genEvents(50000, "e_mum");
		gen.genEvents(10000,  "e_all");
//		gen.genEvents(1000,  "pim_pi0");
//		gen.genEvents(100000,"pim_g_g");
//		gen.genEvents(50000, "pim_gam");
//		gen.genEvents(80000, "pim_n");
//		gen.genEvents(100000,"pim_mu_mu");
//		gen.genEvents(1000,"e_g");
//		gen.genEvents(100000,"pim_mup");
//		gen.genEvents(100000,"pim_mum");
//		gen.decaytest(10000);
	}

}
