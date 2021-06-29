package org.clas.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jlab.clas.physics.Particle;

public class evgen {
	
	Particle p = new Particle();
	
	List<String> out = new ArrayList<String>();
	int   nc;
	
	GenEvent ge ;
	
	public evgen() {
		
	}
	
	public void genEvents(String str, int val) {
		if(str=="pim_g_g") ge = new gen_pim_g_g();
		if(str=="e_pim")   ge = new gen_e_pim();
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
	
	public class gen_e_pim extends GenEvent {		
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nc));
			out.add(getString(getParticle(  11,7,9,10,25,240,260),1));
			out.add(getString(getParticle(-211,2.0,5,15,30, 59, 61),2));
			return out;			
		}		
	}
	
	public String getHeader(int np, int n) {
		return np+" "+n+" 0. 0. 0. 0. 0. 0. 0. 0.";		
	}
	
	public Particle getParticle(double... pargs) {
		int   pid = (int) pargs[0];
		double pm = pargs[1]+Math.random()*(pargs[2]-pargs[1]);
	    double th = Math.toRadians((pargs[3]+Math.random()*(pargs[4]-pargs[3])));
	    double ph = Math.toRadians((pargs[5]+Math.random()*(pargs[6]-pargs[5])));
	    double px = pm*Math.sin(th)*Math.cos(ph);
	    double py = pm*Math.sin(th)*Math.sin(ph);
	    double pz = pm*Math.cos(th);		
	    p.setVector(pid,px,py,pz,0,0,0);
	    return p;
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
			File outputFile = new File("/Users/colesmith/CLAS12/sim/2gamma/"+str+".lund");
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
		gen.genEvents("e_pim",10000);
//		gen.genEvents("pim_g_g",1000);
	}

}
