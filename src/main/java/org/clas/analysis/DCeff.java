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

  
    @Override    
    public void processEvent(DataEvent event) {  
    	
    	ev.init(event);   	
    	ev.setHipoEvent(isHipo3Event);
    	ev.setEventNumber(getEventNumber());
    	ev.requireOneElectron(false);
    	ev.setElecTriggerSector(ev.getElecTriggerSector(shiftTrigBits(getRunNumber())));
	    
    	if(!ev.procEvent(event)) return;    	
    	
    	System.out.println("Trigger Sector = "+ev.trigger_sect+" Event = "+getEventNumber());
    	
//   	for (int i=0; i<4; i++) part[i].clear(); 
    	
     	part[0]  = ev.getPART(   0, 1.0f, 2, -1); 
    	part[1]  = ev.getPART( 212, 1.0f, 2, -1);
    	part[2]  = ev.getPART(  11, 1.0f, 2, -1);
    	part[3]  = ev.getPART(  22, 0.1f, 2,  0);
    	
    	for (int i=0; i<4; i++) {
    		if(part[i].size()>0) output(i);
    	}
    	
    	System.out.println(" ");
    }
    
    public void output(int i) {
    	for (Particle p : part[i]) {
    		List<Particle> e = ev.getECAL((int)p.getProperty("pindex"));
    		for (Particle ecal : e) {
    			System.out.println("i:"+i+" "+(int)ecal.getProperty("sector")+" "+(int)ecal.getProperty("layer")+" "+ecal.pid()+" "+(float)Math.toDegrees(ecal.theta())+" "+(float)ecal.p()+" "+ecal.charge());
    		}
    	}  	
    }

}
