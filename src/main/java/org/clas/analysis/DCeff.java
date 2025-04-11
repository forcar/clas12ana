package org.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.geom.prim.Point3D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

public class DCeff extends DetectorMonitor {
	
    Event ev = new Event();
	IndexedList<List<Particle>> ecpart = new IndexedList<List<Particle>>(1);
	List<Particle>[]  part = (ArrayList<Particle>[]) new ArrayList[4];

    public DCeff(String name) {
        super(name);
    	
        dgmActive=true; 
        this.setDetectorTabNames("RAW",
        		                 "EFF");

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
    	FitSummary.clear();
    	Fits.clear();
    	tl.fitData.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
    public void init(int run) {
	    setRunNumber(run);
	    runlist.add(run);     	
    }
    
    @Override  
    public void createHistos(int run) {
    	histosExist = true;
    	System.out.println(getDetectorName()+".createHistos("+run+")");
	    init(run);
	    createRAW(0);
    }
    
    public void createRAW(int st) {
    	
    	String[] lab = {"TRIGGERS","ELECTRONS","PROTONS"};
    	dgm.add("RAW",6,3,0,st,getRunNumber()); 
		for(int pid=0; pid<3; pid++){
			for(int is=1; is<7; is++) {
			dgm.makeH1("raw"+10*is+pid,20,5,25,50,lab[pid],"THETA(deg) ","COUNTS"); 			 
			}
    	}
    }
        
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);   	
    }
    
    public void plotSummary(int run) {
        setRunNumber(run);
        plot("RAW"); 
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,10*getActiveSector()+getActivePC(),0,getRunNumber());
    }
    
    public List<Particle> makePART(int pid, int stat, int chrg, float pmin) {
    	List<Particle> olist = new ArrayList<Particle>();        
        for (Particle p : ev.getPART(pmin,pid,stat)) olist.add(p);
        return olist;    	
    }
  
    @Override    
    public void processEvent(DataEvent event) {  
    	
    	ev.init(event);   	
    	ev.setHipoEvent(isHipo3Event);
    	ev.setEventNumber(getEventNumber());
    	ev.requireOneElectron(false);
    	ev.setElecTriggerSector(ev.getElecTriggerSector(shiftTrigBits(getRunNumber())));
	    
    	if(!ev.procEvent(event)) return;
    	
    	System.out.println("Trigger Sector = "+ev.trigger_sect);
    	
    	for (int i=0; i<4; i++) part[i].clear(); 
    	
     	part[0]  = makePART(0,   4,-1,0.1f); //
    	part[1]  = makePART(212, 2,-1,0.1f);
    	part[2]  = makePART(11,  2,-1,0.1f);
    	part[3]  = makePART(2212,2,+1,0.1f);
    	
    	for (int i=0; i<4; i++) {
    		if(part[i].size()>0) output(i);
    	}
    }
    
    public void output(int i) {
    	for (Particle p : part[i]) {
    		System.out.println("i:"+i+" "+p.pid()+" "+p.p()+" "+p.charge()+" "+p.getProperty("status"));
    	}
    	
    }

}
