package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.tools.EBMC;
import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc1 extends DetectorMonitor {
	
	Event ev = new Event();
	EBMC  eb = new EBMC();
    List<DetectorParticle> np = new ArrayList<DetectorParticle>();
    HipoDataSync  writer = null;
	  
    List<Particle> phot = null;
    List<Particle> neut = null;
    
    public ECmc1(String name) {
        super(name);
        
        dgm.setDetectorTabNames("MAIN","STATUS","EFFICIENCY");

        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);

        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECmc2.localinit()");
        
        engine.init();
        engine.isMC = true;
        engine.setVariation("default");
        engine.setPCALTrackingPlane(9);
        engine.setCalRun(2);  
        
        eb.getCCDB(2);
        eb.setThresholds("Pizero",engine);
        eb.setGeom("2.5");
        eb.setGoodPhotons(12);
        eb.isMC = true;
        
        tl.setFitData(Fits);
        writer = new HipoDataSync();
        writer.open("/Users/colesmith/CLAS12ANA/ECmc1/photon_demo.hipo");
    }
    
    public void localclear() {
    	System.out.println("ECmc2.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	slider.setValue(0);
    }  
    
    @Override
    public void createHistos(int run) {  
	    System.out.println("ECmc2:createHistos("+run+")");
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;    	
    	createMAIN(0);
       	createEFFICIENCY(0);
       	createEFFICIENCY(1);
    }
        
    public void createMAIN(int st) {
    	       
    	switch (st) {        
        case 0:    
    		dgm.add("MAIN",3,7,0,st,getRunNumber());    
    		dgm.makeH1("h211", 5,0.5,5.5,-1,           "","PCAL Clusters",1,2,2);                              dgm.cc("h211", true,false,0,0,0,0);  
    		dgm.makeH1("h221", 5,0.5,5.5,-1,           "","ECIN Clusters",1,5,2);                              dgm.cc("h221", true,false,0,0,0,0);
    		dgm.makeH1("h231", 5,0.5,5.5,-1,           "","ECOU Clusters",1,4,2);                              dgm.cc("h231", true,false,0,0,0,0);
    		dgm.makeH1("h211g",5,0.5,5.5,-1,           "","PCAL Clusters",1,4,2);                              dgm.cc("h211g",true,false,0,0,0,0);  
    		dgm.makeH1("h211n",5,0.5,5.5,-2,           "","PCAL Clusters",1,43,2);
    		dgm.makeH1("h221g",5,0.5,5.5,-1,           "","ECIN Clusters",1,4,2);                              dgm.cc("h221g",true,false,0,0,0,0);  
    		dgm.makeH1("h221n",5,0.5,5.5,-2,           "","ECIN Clusters",1,43,2);
    		dgm.makeH1("h231g",5,0.5,5.5,-1,           "","ECOU Clusters",1,4,2);                              dgm.cc("h231g",true,false,0,0,0,0);  
    		dgm.makeH1("h231n",5,0.5,5.5,-2,           "","ECOU Clusters",1,43,2);
    		dgm.makeH2("st1",  5,0.5,5.5,3,0,3,-1,     "","PCAL Clusters","STATUS");                           dgm.cc("st1", false,true,0,0,1,100);
    		dgm.makeH2("st4",  5,0.5,5.5,3,0,3,-1,     "","ECIN Clusters","STATUS");                           dgm.cc("st4", false,true,0,0,1,100);
    		dgm.makeH2("st7",  5,0.5,5.5,3,0,3,-1,     "","ECOU Clusters","STATUS");                           dgm.cc("st7", false,true,0,0,1,100);                    
    		dgm.makeH2("set1", 50,0,1,4,1.5,5.5,-1,    "","PCAL Energy Fraction","PCAL Clusters");             dgm.cc("set1",false,true,0,0,1,100); 
    		dgm.makeH2("set4", 50,0,1,4,1.5,5.5,-1,    "","ECIN Energy Fraction","ECIN Clusters");             dgm.cc("set4",false,true,0,0,1,100); 
    		dgm.makeH2("set7", 50,0,1,4,1.5,5.5,-1,    "","ECOU Energy Fraction","ECOU Clusters");             dgm.cc("set7",false,true,0,0,1,100); 
    		dgm.makeH2("dee1", 50,0,3.8,20,-0.2,0.2,-1,"N#gamma=1","Photon Energy (GeV)","#DeltaE/E");         dgm.cc("dee1",false,true,0,0,1,100);
    		dgm.makeH2("dee2", 40,0,3.8,15,-0.2,0.2,-1,"N#gamma=2","Photon Energy (GeV)","#DeltaE/E");         dgm.cc("dee2",false,true,0,0,1,100); 
    		dgm.makeH2("dee3", 30,0,3.8,10,-0.2,0.2,-1,"N#gamma=3","Photon Energy (GeV)","#DeltaE/E");         dgm.cc("dee3",false,true,0,0,1,100); 
    		dgm.makeH1("hr21", 50,0,3.8,-1,            "","Photon Energy (GeV)","Avg. No. PCAL Clusters",1,2); dgm.cc("hr21",false,false,1,1.2f,0,0);
    		dgm.makeH1("hr22", 50,0,3.8,-1,            "","Photon Energy (GeV)","Avg. No. ECIN Clusters",1,5); dgm.cc("hr22",false,false,1,1.2f,0,0);
    		dgm.makeH1("hr23", 50,0,3.8,-1,            "","Photon Energy (GeV)","Avg. No. ECOU Clusters",1,4); dgm.cc("hr23",false,false,1,1.2f,0,0);            
    		dgm.makeH1("h212", 50,0.,1,-1,             "","PCAL Energy Fraction",1,2);
    		dgm.makeH1("h222", 50,0.,1,-1,             "","ECIN Energy Fraction",1,5);
    		dgm.makeH1("h232", 50,0.,1,-1,             "","ECOU Energy Fraction",1,4);
    	}

    }
    	
    public void createEFFICIENCY(int st) {
        
    	switch (st) {        
        case 0:                        
    		dgm.add("EFFICIENCY",2,2,0,st,getRunNumber());
        	dgm.makeH1("ef11",50,0,3.8,-1,"n>0 any layer",  "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef13",50,0,3.8,-2,"n>0 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef14",50,0,3.8,-2,"n>0 layer 1,4",  "Photon Energy (GeV)",1,5);
        	dgm.makeH1("ef15",50,0,3.8,-2,"n>0 layer 1,4,7","Photon Energy (GeV)",1,4);
        	dgm.makeH1("ef21",50,0,3.8,-1,"n=1 each layer", "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef22",50,0,3.8,-2,"n=1 any layer",  "Photon Energy (GeV)",1,1);
        	dgm.makeH1("ef23",50,0,3.8,-2,"n=1 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef24",50,0,3.8,-2,"n=1 layer 1,4",  "Photon Energy (GeV)",1,5);
        	dgm.makeH1("ef25",50,0,3.8,-2,"n=1 layer 1,4,7","Photon Energy (GeV)",1,4);
            break;
        case 1:            
    		dgm.add("EFFICIENCY",4,3,0,st,getRunNumber());
        	dgm.makeH1("h10", 50,0,3.8,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h11", 50,0,3.8,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h12", 50,0,3.8,-1,"PCAL","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h13", 50,0,3.8,-1,"PCAL","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h14", 50,0,3.8,-1,"ECIN","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h15", 50,0,3.8,-1,"ECOU","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h21", 50,0,3.8,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h22", 50,0,3.8,-1,"ECIN","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h23", 50,0,3.8,-1,"ECOU","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h210",50,0,3.8,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h220",50,0,3.8,-1,"ECIN","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h230",50,0,3.8,-1,"ECOU","Photon Energy (GeV)",1,5);
    	}

    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary();
    	plotAnalysis(run);
    }
    
    public void plotSummary() {
    	plotMCHistos();
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;
    	plotMCHistos();
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
    	geteff();
    	writer.close();
    	isAnalyzeDone = true;
    }
    
    public List<Particle> getPART(double thr, int pid) {   	
    	List<Particle> olist = new ArrayList<Particle>();    
    	for (Particle p : ev.getParticle(pid)) {
    		short status = (short) p.getProperty("status");
    		if(status>=2000 && status<3000 && p.p()>=thr) olist.add(p); 
    	}          	
       return olist;    	
    }
    
    @Override
    public void processEvent(DataEvent de) {
    	
		boolean goodev = eb.readMC(de) && eb.pmc.size()==1;
		
		int n=0,ng=0,npp=0,nphot=0,nneut=0;
		float pthresh=0.01f;	
		
    	List<Integer> s1 = new ArrayList<Integer>();
    	List<Integer> s4 = new ArrayList<Integer>();
    	List<Integer> s7 = new ArrayList<Integer>();	
    	
        if (goodev) {                   
            engine.processDataEvent(de); 
            
        	eb.readEC(de,"ECAL::clusters");
        	np.clear(); np = eb.getNeutralPart(); 
        	
        	eb.getRECBanks(de,eb.eb); writer.writeEvent(de);
       	
        	ev.init(de);        	
        	ev.setEventNumber(2);
        	ev.setStartTimeCut(-2000);
        	ev.requireOneElectron(false);
       	    ev.setElecTriggerSector(2);
       	    
       	    if(ev.procEvent(de)) {     		 
     		  phot = getPART(0.0, 22);      nphot = phot.size();
     		  neut = getPART(0.0,2112);     nneut = neut.size();
       	    }
       	    
        	float refP = (float) eb.pmc.get(0).p();        	
        	dgm.fill("h10",refP);
        	        	
        	int [][][] pid = new int[6][3][6]; pid = ev.getECALPID(2);

        	dgm.fill("h211g",pid[1][0][1]);		
        	dgm.fill("h211n",pid[0][0][1]);
        	dgm.fill("h221g",pid[1][1][1]);		
        	dgm.fill("h221n",pid[0][1][1]);		
        	dgm.fill("h231g",pid[1][2][1]);		   		
        	dgm.fill("h231n",pid[0][2][1]);		
    		
        	npp=0;        	        	
        	
        	if(np.size()>0) { // 1 or more neutral particles id=22,2112
        		 
        		double  sf = SamplingFractions.getMean(22, np.get(0), eb.ccdb);
        		double en1 = np.get(0).getEnergy(DetectorType.ECAL)/sf;             		
        		double dee = (refP-en1)/refP;
        		
        		dgm.fill("h11",refP); n++;
        		
    	    	int n1=0,n4=0,n7=0,sum=0;
    	    	s1.clear(); s4.clear(); s7.clear();   	    	
    	    	double e1=0, e4=0, e7=0, e1p=0, e4p=0, e7p=0;
    	    	
        	    for (DetectorParticle phot : np) {
        	    	double p1 = phot.getEnergy(DetectorType.ECAL);
        	    	if(p1>pthresh) {
        	    	npp++;
        	    	for (DetectorResponse dr : phot.getDetectorResponses()) {
        	    		int lay = dr.getDescriptor().getLayer();
        	    		double e = dr.getEnergy(); double t = dr.getTime(); 
        	    		int h1 = dr.getHitIndex(); int h2 = dr.getAssociation(); float pat = (float) dr.getPath();       	    		
        	    		if(lay==1) {n1++; e1p+=e; s1.add(dr.getStatus()); if(n1==1) {e1=e; dgm.fill("h210",refP);}}
        	    		if(lay==4) {n4++; e4p+=e; s4.add(dr.getStatus()); if(n4==1) {e4=e; dgm.fill("h220",refP);}}
        	    		if(lay==7) {n7++; e7p+=e; s7.add(dr.getStatus()); if(n7==1) {e7=e; dgm.fill("h230",refP);}}
        	    	}
        	    	}
        	    }
        	    
        		if(n1==1)       dgm.fill("dee1",refP,dee);
        		if(n1==1&&n4>1) dgm.fill("dee2",refP,dee);
        		if(n1==1&&n4>2) dgm.fill("dee3",refP,dee);   

        		dgm.fill("h211",n1); for (Integer i : s1) dgm.fill("st1",n1,i);
        		dgm.fill("h221",n4); for (Integer i : s4) dgm.fill("st4",n4,i);
        		dgm.fill("h231",n7); for (Integer i : s7) dgm.fill("st7",n7,i);
        	    
        		dgm.fill("h21",refP,n1);dgm.fill("h22",refP,n4); dgm.fill("h23",refP,n7);
        	    
                if(n1>1) dgm.fill("h212",e1/e1p); dgm.fill("set1",e1/e1p,n1);
                if(n4>1) dgm.fill("h222",e4/e4p); dgm.fill("set4",e4/e4p,n4);
                if(n7>1) dgm.fill("h232",e7/e7p); dgm.fill("set7",e7/e7p,n7);
        	}
        	if(npp==1) { // only 1 neutral particle id=22,2112            		
        		dgm.fill("h12",refP); ng++;
        	    for (DetectorParticle phot : np) {
        	    	double p1 = phot.getEnergy(DetectorType.ECAL);
        	    	int sum=0;
//           	    	if(p1>pthresh) {
        	    	for (DetectorResponse dr : phot.getDetectorResponses()) {
        	    		int lay = dr.getDescriptor().getLayer();
        	    		if(lay==1) sum+= 100;
        	    		if(lay==4) sum+=  40;
        	    		if(lay==7) sum+=   7;
        	    	}
        	    	if(sum>=100) dgm.fill("h13",refP);
        	    	if(sum>=140) dgm.fill("h14",refP);
        	    	if(sum==147) dgm.fill("h15",refP);
//        	    	}
        		}
        	}
        }  		
      
    }
   
    public void plotMCHistos() {      
        plot("MAIN");
        plot("EFFICIENCY");   
    }
    
    @Override
    public void plotEvent(DataEvent de) {
//        analyze();         
    }
    
    public void geteff() {
    	dgm.geteff("ef11","h11", "h10");
    	dgm.geteff("ef13","h210","h10");
    	dgm.geteff("ef14","h220","h10");
    	dgm.geteff("ef15","h230","h10");
    	dgm.geteff("ef21","h11", "h10");
    	dgm.geteff("ef22","h12", "h10");
    	dgm.geteff("ef23","h13", "h10");
    	dgm.geteff("ef24","h14", "h10");
    	dgm.geteff("ef25","h15", "h10");
    	dgm.geteff("hr21","h21", "h210");
    	dgm.geteff("hr22","h22", "h220");
    	dgm.geteff("hr23","h23", "h230");	
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,0, getRunNumber());
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
