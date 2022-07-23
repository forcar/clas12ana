package org.clas.analysis;

import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.Particle;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import java.util.ArrayList;
import java.util.HashMap;

public class ECmon extends DetectorMonitor {
	
    public ECmon(String name) {
        super(name);
        dgmActive=true; 
        this.setDetectorTabNames("MODE1",
        		                 "ADC",
        		                 "TDC",
        		                 "PEDS");
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();       
    }
    
    public void localinit() {
        System.out.println(getDetectorName()+".localinit()"); 
        tl.setFitData(Fits);
    } 
    
    public void localclear() {
    	System.out.println(getDetectorName()+".localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
    public void init(int run) {
	    setRunNumber(run);
	    runlist.add(run);     	
    }
    
    @Override  
    public void createHistos(int run) {
	    System.out.println(getDetectorName()+".createHistos("+run+")");
	    init(run);
	    createMODE1(0);
    }
    
    public void createMODE1(int st) {
    	String[] det = {"PCAL","ECIN","ECOU"}; 
    	String[] uvw = {"U","V","W"};
    	
    	for (int is=1; is<7; is++) {
    	for (int id=0; id<3; id++) {int np = id==0?68:36; dgm.add("MODE1",3,4,10*is+id,st,getRunNumber());
    	for (int il=0; il<3; il++) {String tit = "Sector "+is+" "+det[id]+" "+uvw[il]+" STRIPS";
    		dgm.makeH2("m1"+is+id+il,100, 0,250,100,0,300,-1,tit,"FADC","TDC");
    		dgm.makeH2("m2"+is+id+il,100, 0, 99,100,0,300,-1,"","FADC","TDC");
    		dgm.makeGraph("m3"+is+id+il,0,"","FADC","TDC EFFICIENCY");
    		dgm.makeH2("m4"+is+id+il,100, 0, 99,np,1,np+1,-1,"","FADC","TDC");
    	}
    	}
    	}
    	
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plot("MODE1");
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,10*getActiveSector()+getActiveLayer(),0,getRunNumber());
    }   
       

}
