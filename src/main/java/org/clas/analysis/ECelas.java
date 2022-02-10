package org.clas.analysis;

import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.jnp.hipo4.data.Bank;

//import org.jlab.jnp.hipo.io.*;
//import org.jlab.jnp.hipo.data.*;
//import org.jlab.jnp.reader.*;
//import org.jlab.jnp.physics.*;

import java.util.ArrayList;
import java.util.HashMap;


public class ECelas extends DetectorMonitor {
	
	boolean processWagon = true;
	boolean processEvent = true;
	double beamEnergy;
	double mp = 0.93828;
	
    public ECelas(String name) {
        super(name);
        dgmActive=true; 
        this.setDetectorTabNames("WAGON",
        		                 "EVENT");
        this.use123Buttons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();       
    }
    
    public void localinit() {
        System.out.println("ECelas.localinit()");           	
    } 
    
    public void localclear() {
    	System.out.println("ECelas.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    } 
    
    @Override  
    public void createHistos(int run) {
	    System.out.println("ECmc:createHistos("+run+")");
	    setRunNumber(run);
	    runlist.add(run);  
	    beamEnergy = getBeamEnergy(run);
    	histosExist = true;
    	createWAGON(0);
    	createEVENT(0);
    }
    
    public void createWAGON(int st) {
    	switch (st) {        
        case 0: 
        	dgm.add("WAGON",6,3,0,st,getRunNumber());         	
        	for(int is=1; is<7; is++) dgm.makeH1("WA0"+is,180,  0,360, -1,"SECTOR "+is,"#Delta#phi (DEG)","");
         	for(int is=1; is<7; is++) dgm.makeH1("WA1"+is,100, -1,  1, -1,"","#Delta EBEAM (GEV)","");
        	for(int is=1; is<7; is++) dgm.makeH1("WA2"+is,100,-20, 20, -1,"","#Delta#theta (DEG)","");
    	}
    }
    
    public void createEVENT(int st) {
    	switch (st) {        
        case 0: 
        	dgm.add("EVENT",6,2,0,st,getRunNumber());
        	for(int is=1; is<7; is++) {dgm.makeH2("EV1"+is,100,0.5,3.5,60,6,40,-1,"SECTOR "+is,"W (GEV)","#theta (DEG)");dgm.cc("EV"+is, false, true, 0, 0, 0, 0);}
        	for(int is=1; is<7; is++)  dgm.makeH1("EV2"+is,100,0.5,3.5,-1,"","W (GEV)","");	
        }
    }    
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plot("WAGON");
    	plot("EVENT");    	 	
    }
    
    public void processEvent(DataEvent event) {
    	
    	if(processWagon) processWag(event);
    	if(processEvent) processEvt(event);
    }
    
    public boolean processWag(DataEvent event) {
  
        DataBank RecPart = null, RecCal = null;
        
        if (event.hasBank("REC::Particle"))   RecPart = event.getBank("REC::Particle");
        if (event.hasBank("REC::Calorimeter")) RecCal = event.getBank("REC::Calorimeter");
        
    	HashMap<Integer,ArrayList<Integer>> part2calo = mapByIndex(RecCal);

        if (RecPart==null || RecPart.rows()==0) return false;

        ArrayList<Integer> eleCandi = new ArrayList<>();
        ArrayList<Integer> proCandi = new ArrayList<>();

        for (int ipart=0; ipart<RecPart.rows(); ipart++) {
            final int    pid = RecPart.getInt("pid",ipart);
            final int charge = RecPart.getInt("charge",ipart);
            final int status = RecPart.getInt("status",ipart);       

            final boolean isFD = (int)(Math.abs(status)/1000) == 2;

            if (isFD && charge < 0) eleCandi.add(ipart);
            if (pid==2212)          proCandi.add(ipart);
        }

        if (eleCandi.isEmpty() || proCandi.isEmpty()) return false;

        if (eleCandi.size()==1 && proCandi.size()==1) {
        	double epx  = RecPart.getFloat("px",eleCandi.get(0));
        	double epy  = RecPart.getFloat("py",eleCandi.get(0));
        	double epz  = RecPart.getFloat("pz",eleCandi.get(0));
        	double ep   = Math.sqrt(epx*epx + epy*epy + epz*epz);
        	double cthe = epz/ep, sthe=Math.sqrt(1-cthe*cthe);
        	Particle neg = new Particle(211,epx,epy,epx);

        	double ppx  = RecPart.getFloat("px",proCandi.get(0));
        	double ppy  = RecPart.getFloat("py",proCandi.get(0));
        	double ppz  = RecPart.getFloat("pz",proCandi.get(0));
        	double pp   = Math.sqrt(ppx*ppx + ppy*ppy + ppz*ppz);
        	double cthp = ppz/pp, sthp=Math.sqrt(1-cthp*cthp);
        	Particle pos = new Particle(211,ppx,ppy,ppz);

        	double ebeam  = (mp*(cthe + sthe*cthp/sthp)-mp)/(1-cthe); 
        	double delthe = Math.acos(cthe) - Math.atan(pp*sthp/(beamEnergy-pp*cthp));
        	double delphi = (neg.phi()-pos.phi())*180/3.14159;
        	if (delphi<0) delphi = 360+delphi;
        	
        	if (RecCal!=null) {
        		int s = part2calo.containsKey(eleCandi.get(0)) ? RecCal.getInt("sector", part2calo.get(eleCandi.get(0)).get(0)):0;           
        		if(s>0) {dgm.fill("WA0"+s, delphi);
        		         if(Math.abs(180-delphi)<5) dgm.fill("WA1"+s, beamEnergy-ebeam);
        		         dgm.fill("WA2"+s, delthe*180/3.14159);
        		}
        	}

        	return delthe > -0.026 ? true : false;
        } 

        return false;
  	
    }

    public void processEvt(DataEvent event) { 
    	
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
                    	dgm.fill("EV1"+s,w,p.theta()*180/Math.PI); dgm.fill("EV2"+s,w);
                    }
                 }
        	}
        } 
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    }   
        
}
