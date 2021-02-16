package org.clas.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class evgen {
	
	
	float[]  e = {2,4,4};
	float[]th1 = {18,15,0};
	float thmin=10;
	float thmax=20;
	float phmin=59;
	float phmax=61;
	float p=0f;
	float pmin=0.2f;
	float pmax=2.0f;
			
	
	public evgen() {
		
	}
	
	public void genevents() {
		writeFile();
	}
	
	public void writeFile() {
		
		String line = new String();
		
		try { 
			File outputFile = new File("/Users/colesmith/CLAS12/sim/2gamma/evgen.lund");
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("evgen writing to "+outputFile.getName());
			
			int nevents=100000;
			
			for (int n=1; n<nevents+1; n++) {
				for (int l=0; l<4; l++) {
					line =  getev(l,n);
					if (line!=null) {
//				      System.out.println(line);
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
	
	public String getev(int l, int n) {
	
		float theta,phi,px,py,pz;
		
		switch (l) {
		case 0:
			return "3 "+n+" 0. 0. 0. 0. 0. 0. 0. 0.";
		case 1:
			  theta = (float) Math.toRadians(th1[l-1]);
			    phi = (float) Math.toRadians((phmin+Math.random()*(phmax-phmin))+190);
			     px = (float) (e[l-1]*Math.sin(theta)*Math.cos(phi));
			     py = (float) (e[l-1]*Math.sin(theta)*Math.sin(phi));
			     pz = (float) (e[l-1]*Math.cos(theta));			    
			    return "1 -1. 1 -211 0 0 "+String.format("%.4f",px)+" "
			    		               +String.format("%.4f",py)+" "
			                           +String.format("%.4f",pz)+" "
			    		               +String.format("%.4f",e[l-1])+" "
					                   +String.format("%.4f",0.0)+" "
					                   +String.format("%.4f",0.0)+" "
					                   +String.format("%.4f",0.0)+" "
					                   +String.format("%.4f",-3.); 
			    
		case 2:
			  theta = (float) Math.toRadians(th1[l-1]);
			    phi = (float) Math.toRadians((phmin+Math.random()*(phmax-phmin)));
			     px = (float) (e[l-1]*Math.sin(theta)*Math.cos(phi));
			     py = (float) (e[l-1]*Math.sin(theta)*Math.sin(phi));
			     pz = (float) (e[l-1]*Math.cos(theta));			    
			    return "2 0. 1 22 0 0 "+String.format("%.4f",px)+" "
			    		               +String.format("%.4f",py)+" "
			                           +String.format("%.4f",pz)+" "
			    		               +String.format("%.4f",e[l-1])+" "
					                   +String.format("%.4f",0.0)+" "
					                   +String.format("%.4f",0.0)+" "
					                   +String.format("%.4f",0.0)+" "
					                   +String.format("%.4f",-3.); 
		case 3:
			      p = (float) (pmin+Math.random()*(pmax-pmin));
			  theta = (float) Math.toRadians((thmin+Math.random()*(thmax-thmin)));
			    phi = (float) Math.toRadians((phmin+Math.random()*(phmax-phmin)));
			     px = (float) (p*Math.sin(theta)*Math.cos(phi));
			     py = (float) (p*Math.sin(theta)*Math.sin(phi));
			     pz = (float) (p*Math.cos(theta));
			    return "3 0. 1 22 0 0 "+String.format("%.4f",px)+" "
  		                               +String.format("%.4f",py)+" "
                                       +String.format("%.4f",pz)+" "
		                               +String.format("%.4f",e[l-1])+" "
	                                   +String.format("%.4f",0.0)+" "
	                                   +String.format("%.4f",0.0)+" "
	                                   +String.format("%.4f",0.0)+" "
	                                   +String.format("%.4f",-3.);
		}
		return null;
	}
	
	public static void main(String[] args) {
		evgen gen = new evgen();
		gen.genevents();
	}

}
