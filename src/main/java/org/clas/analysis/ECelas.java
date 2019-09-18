package org.clas.analysis;

import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

//import org.jlab.jnp.hipo.io.*;
//import org.jlab.jnp.hipo.data.*;
//import org.jlab.jnp.reader.*;
//import org.jlab.jnp.physics.*;

import java.util.ArrayList;
import java.util.HashMap;


public class ECelas extends DetectorMonitor {
	
//	EventFilter filter = null;
	double beamEnergy = 2.22193;
	double mp = 0.93828;
	
    public ECelas(String name) {
        super(name);
        this.setDetectorTabNames("WvTh",
        		                 "W");
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
//        filter = new EventFilter("11");
    }
    
    public void localinit() {
        configEngine("elec");     	
    } 
    
    @Override  
    public void createHistos(int run) {
	     setRunNumber(run);
	     runlist.add(run);  
	     createSector2D(0);
	     createSector1D(1);
    }
    
    public void createSector1D(int k) {

	    int run = getRunNumber();
        H1F h;

        DataGroup dg = new DataGroup(3,2); 
    
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H1F("hi_w_"+tag,"hi_w_"+tag,100,0.8,1.2);
            h.setTitleX("Sector "+is+" W (GeV)");
            dg.addDataSet(h, is-1);          
        }
        
        this.getDataGroup().add(dg,0,0,k,run);    
    }
    
    public void createSector2D(int k) {

	    int run = getRunNumber();
        H2F h;

        DataGroup dg = new DataGroup(3,2); 
    
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_w_"+tag,"hi_w_"+tag,50,0.8,1.2,28,6.,21.);
            h.setTitleY("Theta (deg)");
            h.setTitleX("Sector "+is+" W (GeV)");
            dg.addDataSet(h, is-1);          
        }
        
        this.getDataGroup().add(dg,0,0,k,run);    
    }    
    
    @Override       
    public void plotHistos(int run) {
        setRunNumber(run);
        plotSummary(0);    	
        plotSummary(1);    	
    }

    public void processEvent(DataEvent event) { 
    	
    	if(!event.hasBank("REC::Calorimeter")) return;
        if(!event.hasBank("REC::Particle"))    return;
    	
    	DataBank bank = event.getBank("REC::Particle");
    	DataBank calb = event.getBank("REC::Calorimeter");
    	
    	HashMap<Integer,ArrayList<Integer>> part2calo = mapByIndex(calb);
    	
        for(int loop = 0; loop < bank.rows(); loop++){
        	if(bank.getInt("pid",loop)!=0) {
        		Particle p = new Particle(bank.getInt("pid", loop),
                                          bank.getFloat("px", loop),
                                          bank.getFloat("py", loop),
                                          bank.getFloat("pz", loop),
                                          bank.getFloat("vx", loop),
                                          bank.getFloat("vy", loop),
                                          bank.getFloat("vz", loop));
                    p.setProperty("beta", bank.getFloat("beta", loop));
                   
                 if (p.pid()==11&&part2calo.containsKey(loop)) {
                    int s = calb.getInt("sector", part2calo.get(loop).get(0));
               
                    double nu = beamEnergy - p.e();
                    double q2 = 2*beamEnergy*p.p()*(1-Math.cos(p.theta()));
                    double w2 = -q2 + mp*mp + 2*mp*nu;
                    if (w2>0) {
                    	double w = Math.sqrt(w2);
                    	((H1F) this.getDataGroup().getItem(0,0,1,getRunNumber()).getData(s-1).get(0)).fill(w);
                        ((H2F) this.getDataGroup().getItem(0,0,0,getRunNumber()).getData(s-1).get(0)).fill(w,p.theta()*180/Math.PI);
                    }
                 }
        	}
        } 
    }
    
    public void plotSummary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));
    }    
    
    
}
