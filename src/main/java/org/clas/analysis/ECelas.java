package org.clas.analysis;

import org.clas.tools.KinLib;
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
	
	KinLib kin = new KinLib();
	float kinqw[] = null;
	boolean doWA = true, doEV = true;
	float beamEnergy;
	double mp = kin.mp;
	DataBank RecPart = null, RecCal = null;
	HashMap<Integer,ArrayList<Integer>> part2calo  = null;
	
    public ECelas(String name) {
        super(name);
        dgmActive=true; 
        this.setDetectorTabNames("WAGON",
        		                 "EVENT",
        		                 "KCOR",
        		                 "RCOR");
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
    
    public void init(int run) {
	    setRunNumber(run);
	    runlist.add(run);  
	    beamEnergy = getBeamEnergy(run);
	    kinqw = kin.getQW(beamEnergy, 0.2f*beamEnergy, 6.0f);   	
    }
    
    @Override  
    public void createHistos(int run) {
	    System.out.println("ECmc:createHistos("+run+")");
	    init(run);
    	createWAGON(0);
    	createEVENT(0);
    	createKCOR(0);
    	histosExist = true;
    }
    
    public void createWAGON(int st) {
    	switch (st) {        
        case 0: 
        	dgm.add("WAGON",6,3,0,st,getRunNumber());         	
        	for(int is=1; is<7; is++) dgm.makeH1("wa0"+is,180, 120,240,-1,"SECTOR "+is,"#Delta#phi (DEG)",  "",1,0,"11001100");
         	for(int is=1; is<7; is++) dgm.makeH1("wa1"+is,100, -1,  1, -1,"",          "#Delta EBEAM (GEV)","",1,0,"11001100");
        	for(int is=1; is<7; is++) dgm.makeH1("wa2"+is,100,-20, 20, -1,"",          "#Delta#theta (DEG)","",1,0,"11001100");
    	}
    }
    
    public void createEVENT(int st) {
    	switch (st) {        
        case 0: 
        	dgm.add("EVENT",6,3,0,st,getRunNumber());
        	for(int is=1; is<7; is++) {dgm.makeH2("ev0"+is, 100,0.7,kinqw[1],60,4,40,-1,"SECTOR "+is,"W (GEV)","#theta (DEG)"); dgm.cc("eV0"+is,false,true,0,0,0,0);}
        	for(int is=1; is<7; is++) {dgm.makeH2("ev1"+is, 100,0.7,kinqw[1],60,4,40,-1,"",          "W (GEV)","#theta (DEG)"); dgm.cc("ev1"+is,false,true,0,0,0,0);}
        	for(int is=1; is<7; is++) {dgm.makeH1("ev2a"+is,100,0.7,kinqw[1],-1,"","W (GEV)","",1,0,"1000000"); 
        	                           dgm.makeH1("ev2b"+is,100,0.7,kinqw[1],-2,"","W (GEV)","",1,1,"1000000");}	
        }
    }
    
    public void createKCOR(int st) {
    	switch (st) {
    	case 0:
    		dgm.add("KCOR",6,2,0,st,getRunNumber());
    		for(int is=1; is<7; is++) {dgm.makeH2("kc0"+is,100,40,80,50,-5,5, -1,"SECTOR "+is,"#theta p+ (DEG)","#Delta#phi (DEG)"); dgm.cc("kc0"+is, false,true,0,0,0,0);}
    		for(int is=1; is<7; is++) {dgm.makeH2("kc1"+is,100, 8,31,50,-5,5, -1,"SECTOR "+is,"#theta e- (DEG)","#Delta#phi (DEG)"); dgm.cc("kc1"+is, false,true,0,0,0,0);}
    	}
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plot("WAGON");
    	plot("EVENT");    
    	plot("KCOR");
    }
    
    @Override    
    public void processEvent(DataEvent event) {   	
    	if(!processFilter(event)) return;
    	if(doWA) processWA(event);
    	if(doEV) processEV(event);
    }
    
    public boolean processFilter(DataEvent event) {
        RecCal   = event.hasBank("REC::Calorimeter") ? event.getBank("REC::Calorimeter"):null;
        RecPart  = event.hasBank("REC::Particle")    ? event.getBank("REC::Particle"):null;       
    	boolean test1 = RecPart!=null && RecPart.rows()!=0;
    	boolean test2 = RecCal !=null && RecCal.rows()!=0;
    	if(test2) part2calo = mapByIndex(RecCal,"pindex"); 
    	return test1 && test2;
    }
     
    public boolean processWA(DataEvent event) {
         
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
            int    pid  = RecPart.getInt("pid", eleCandi.get(0));
        	double epx  = RecPart.getFloat("px",eleCandi.get(0));
        	double epy  = RecPart.getFloat("py",eleCandi.get(0));
        	double epz  = RecPart.getFloat("pz",eleCandi.get(0));
        	double ep   = Math.sqrt(epx*epx + epy*epy + epz*epz);
        	double cthe = epz/ep, sthe=Math.sqrt(1-cthe*cthe);
        	Particle neg = new Particle(11,epx,epy,epx);

        	double ppx  = RecPart.getFloat("px",proCandi.get(0));
        	double ppy  = RecPart.getFloat("py",proCandi.get(0));
        	double ppz  = RecPart.getFloat("pz",proCandi.get(0));
        	double pp   = Math.sqrt(ppx*ppx + ppy*ppy + ppz*ppz);
        	double cthp = ppz/pp, sthp=Math.sqrt(1-cthp*cthp);
        	Particle pos = new Particle(2212,ppx,ppy,ppz);

        	double ebeam  = (mp*(cthe + sthe*cthp/sthp)-mp)/(1-cthe); 
        	double delthe = Math.acos(cthe) - Math.atan(pp*sthp/(beamEnergy-pp*cthp));
        	double delphi = (neg.phi()-pos.phi())/kin.d2r;
        	if (delphi<0) delphi = 360+delphi;
        	
        	if (RecCal!=null) {
        		int s = part2calo.containsKey(eleCandi.get(0)) ? RecCal.getInt("sector", part2calo.get(eleCandi.get(0)).get(0)):0;           
        		if(s>0) {dgm.fill("wa0"+s, delphi);
        		         if(delphi>171 && delphi<194) {
        		        	 dgm.fill("wa1"+s, beamEnergy-ebeam);
        		         	 dgm.fill("wa2"+s, delthe/kin.d2r);	 
        		         	 dgm.fill("kc0"+s, Math.acos(cthp)/kin.d2r,delphi-180);
        		         	 dgm.fill("kc1"+s, Math.acos(cthe)/kin.d2r,delphi-180);
        		         }
                         double nu = beamEnergy - ep;
                         double q2 = 2*beamEnergy*ep*(1-cthe);
                         double w2 = -q2 + mp*mp + 2*mp*nu;
                         if (w2>0 && pid==11) {
                         	double w = Math.sqrt(w2);
                         	dgm.fill("ev1"+s,w,Math.acos(cthe)/kin.d2r); dgm.fill("ev2b"+s,w);
                         }
        		         
        		}
        	}

        	return delthe > -0.026 ? true : false;
        } 

        return false;
  	
    }

    public void processEV(DataEvent event) { 
    	    	
        for(int loop = 0; loop < RecPart.rows(); loop++){
        	if(RecPart.getInt("pid",loop)!=0) {
        		Particle p = new Particle(RecPart.getInt("pid", loop),
        		RecPart.getFloat("px", loop),
        		RecPart.getFloat("py", loop),
        		RecPart.getFloat("pz", loop),
        		RecPart.getFloat("vx", loop),
        		RecPart.getFloat("vy", loop),
        		RecPart.getFloat("vz", loop));
                p.setProperty("beta", RecPart.getFloat("beta", loop));
                   
                if (p.pid()==11 && part2calo.containsKey(loop)) {
                    double nu = beamEnergy - p.e();
                    double q2 = 2*beamEnergy*p.p()*(1-Math.cos(p.theta()));
                    double w2 = -q2 + mp*mp + 2*mp*nu;
                    if (w2>0) {
                    	int s = RecCal.getInt("sector", part2calo.get(loop).get(0));
                    	double w = Math.sqrt(w2);
                    	dgm.fill("ev0"+s,w,p.theta()*180/Math.PI); dgm.fill("ev2a"+s,w);
                    }
                    
                }
        	}
        } 
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    }   
        
}
