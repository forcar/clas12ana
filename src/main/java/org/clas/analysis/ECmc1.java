package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.analysis.ECPart.SFFunction;
import org.clas.tools.DataGroupManager;
import org.clas.tools.EBMC;
import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc1 extends DetectorMonitor {
	
	Event ev = new Event();
	EBMC  eb = new EBMC();
    List<DetectorParticle> np = new ArrayList<DetectorParticle>();
    HipoDataSync  writer = null;
	DataGroup dg = null;
  
    List<Particle> phot = null;
    List<Particle> neut = null;
    
    public ECmc1(String name) {
        super(name);
        setDetectorTabNames("MAIN","STATUS","EFFICIENCY");

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
    	
    	String tab = "MAIN", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0;
        tag = st+"_"+k+"_"+run;    
        
    	switch (st) {        
        case 0:    
        	dg = new DataGroup(3,7);
            dg.addDataSet(makeH1("h211",5,0.5,5.5,"","PCAL Clusters",1,2,2),n++); cc("h211",true,false,0,0,0,0);  
            dg.addDataSet(makeH1("h221",5,0.5,5.5,"","ECIN Clusters",1,5,2),n++); cc("h221",true,false,0,0,0,0);
            dg.addDataSet(makeH1("h231",5,0.5,5.5,"","ECOU Clusters",1,4,2),n++); cc("h231",true,false,0,0,0,0);
            dg.addDataSet(makeH1("h211g",5,0.5,5.5,"","PCAL Clusters",1,4,2),n++); cc("h211g",true,false,0,0,0,0);  
            dg.addDataSet(makeH1("h211n",5,0.5,5.5,"","PCAL Clusters",1,43,2),n-1);
            dg.addDataSet(makeH1("h221g",5,0.5,5.5,"","ECIN Clusters",1,4,2),n++); cc("h221g",true,false,0,0,0,0);  
            dg.addDataSet(makeH1("h221n",5,0.5,5.5,"","ECIN Clusters",1,43,2),n-1);
            dg.addDataSet(makeH1("h231g",5,0.5,5.5,"","ECOU Clusters",1,4,2),n++); cc("h231g",true,false,0,0,0,0);  
            dg.addDataSet(makeH1("h231n",5,0.5,5.5,"","ECOU Clusters",1,43,2),n-1);
            dg.addDataSet(makeH2("st1",  5,0.5,5.5,3,0,3,"","PCAL Clusters","STATUS"),n++); cc("st1",false,true,0,0,1,100);
            dg.addDataSet(makeH2("st4",  5,0.5,5.5,3,0,3,"","ECIN Clusters","STATUS"),n++); cc("st4",false,true,0,0,1,100);
            dg.addDataSet(makeH2("st7",  5,0.5,5.5,3,0,3,"","ECOU Clusters","STATUS"),n++); cc("st7",false,true,0,0,1,100);                    
            dg.addDataSet(makeH2("set1", 50,0,1,4,1.5,5.5,"","PCAL Energy Fraction","PCAL Clusters"),n++);cc("set1",false,true,0,0,1,100); 
            dg.addDataSet(makeH2("set4", 50,0,1,4,1.5,5.5,"","ECIN Energy Fraction","ECIN Clusters"),n++);cc("set4",false,true,0,0,1,100); 
            dg.addDataSet(makeH2("set7", 50,0,1,4,1.5,5.5,"","ECOU Energy Fraction","ECOU Clusters"),n++);cc("set7",false,true,0,0,1,100); 
            dg.addDataSet(makeH2("dee1", 50,0,3.8,20,-0.2,0.2,"N#gamma=1","Photon Energy (GeV)","#DeltaE/E"),n++);cc("dee1",false,true,0,0,1,100);
            dg.addDataSet(makeH2("dee2", 40,0,3.8,15,-0.2,0.2,"N#gamma=2","Photon Energy (GeV)","#DeltaE/E"),n++);cc("dee4",false,true,0,0,1,100); 
            dg.addDataSet(makeH2("dee3", 30,0,3.8,10,-0.2,0.2,"N#gamma=3","Photon Energy (GeV)","#DeltaE/E"),n++);cc("dee7",false,true,0,0,1,100); 
            dg.addDataSet(makeH1("hr21",50,0,3.8,"","Photon Energy (GeV)","Avg. No. PCAL Clusters",1,2),n++);cc("hr21",false,false,1,1.2f,0,0);
            dg.addDataSet(makeH1("hr22",50,0,3.8,"","Photon Energy (GeV)","Avg. No. ECIN Clusters",1,5),n++);cc("hr22",false,false,1,1.2f,0,0);
            dg.addDataSet(makeH1("hr23",50,0,3.8,"","Photon Energy (GeV)","Avg. No. ECOU Clusters",1,4),n++);cc("hr23",false,false,1,1.2f,0,0);            
            dg.addDataSet(makeH1("h212",50,0.,1, "","PCAL Energy Fraction",1,2),n++);
            dg.addDataSet(makeH1("h222",50,0.,1, "","ECIN Energy Fraction",1,5),n++);
            dg.addDataSet(makeH1("h232",50,0.,1, "","ECOU Energy Fraction",1,4),n++);
    	}
    	
        getDataGroup().add(dg,0,st,k,run); 
    }
    public void createGENREC(int st) {
    	
    	String tab = "GENREC", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0;
        tag = st+"_"+k+"_"+run;    
        
    	switch (st) {        
        case 0: 
    	}
       	
    }
    	
    public void createEFFICIENCY(int st) {
    	
    	String tab = "EFFICIENCY", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0;
        tag = st+"_"+k+"_"+run;    
        
    	switch (st) {        
        case 0:                        
        	dg = new DataGroup(2,2);
            dg.addDataSet(makeH1("ef11",50,0,3.8,"n>0 any layer",  "Photon Energy (GeV)",1,3),n++);
            dg.addDataSet(makeH1("ef13",50,0,3.8,"n>0 layer 1",    "Photon Energy (GeV)",1,2),n-1);
            dg.addDataSet(makeH1("ef14",50,0,3.8,"n>0 layer 1,4",  "Photon Energy (GeV)",1,5),n-1);
            dg.addDataSet(makeH1("ef15",50,0,3.8,"n>0 layer 1,4,7","Photon Energy (GeV)",1,4),n-1);
            dg.addDataSet(makeH1("ef21",50,0,3.8,"n=1 each layer", "Photon Energy (GeV)",1,3),n++);
            dg.addDataSet(makeH1("ef22",50,0,3.8,"n=1 any layer",  "Photon Energy (GeV)",1,1),n-1);
            dg.addDataSet(makeH1("ef23",50,0,3.8,"n=1 layer 1",    "Photon Energy (GeV)",1,2),n-1);
            dg.addDataSet(makeH1("ef24",50,0,3.8,"n=1 layer 1,4",  "Photon Energy (GeV)",1,5),n-1);
            dg.addDataSet(makeH1("ef25",50,0,3.8,"n=1 layer 1,4,7","Photon Energy (GeV)",1,4),n-1);
            break;
        case 1:            
        	dg = new DataGroup(4,3);
            dg.addDataSet(makeH1("h10", 50,0,3.8,"PCAL","Photon Energy (GeV)",1,3),n++);
            dg.addDataSet(makeH1("h11", 50,0,3.8,"PCAL","Photon Energy (GeV)",1,3),n++);
            dg.addDataSet(makeH1("h12", 50,0,3.8,"PCAL","Photon Energy (GeV)",1,2),n++);
            dg.addDataSet(makeH1("h13", 50,0,3.8,"PCAL","Photon Energy (GeV)",1,5),n++);
            dg.addDataSet(makeH1("h14", 50,0,3.8,"ECIN","Photon Energy (GeV)",1,5),n++);
            dg.addDataSet(makeH1("h15", 50,0,3.8,"ECOU","Photon Energy (GeV)",1,5),n++);
            dg.addDataSet(makeH1("h21", 50,0,3.8,"PCAL","Photon Energy (GeV)",1,3),n++);
            dg.addDataSet(makeH1("h22", 50,0,3.8,"ECIN","Photon Energy (GeV)",1,2),n++);
            dg.addDataSet(makeH1("h23", 50,0,3.8,"ECOU","Photon Energy (GeV)",1,5),n++);
            dg.addDataSet(makeH1("h210",50,0,3.8,"PCAL","Photon Energy (GeV)",1,3),n++);
            dg.addDataSet(makeH1("h220",50,0,3.8,"ECIN","Photon Energy (GeV)",1,2),n++);
            dg.addDataSet(makeH1("h230",50,0,3.8,"ECOU","Photon Energy (GeV)",1,5),n++);

    	}
    	
        getDataGroup().add(dg,0,st,k,run); 
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
    	
		int run = getRunNumber();
		DataGroup  dg0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("MAIN"),run);
		DataGroup dg10 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg11 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		
		boolean goodev = eb.readMC(de) && eb.pmc.size()==1;
		
		int n=0,ng=0, npp=0;
		float pthresh=0.01f;	
		int nphot=0,nneut=0;
    	List<Integer> s1 = new ArrayList<Integer>();
    	List<Integer> s4 = new ArrayList<Integer>();
    	List<Integer> s7 = new ArrayList<Integer>();	
    	
        if (goodev) {                   
            engine.processDataEvent(de); 
            
        	eb.readEC(de,"ECAL::clusters");
        	np.clear(); np = eb.getNeutralPart(); 
        	eb.getRECBanks(de,eb.eb);
        	writer.writeEvent(de);
       	
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
        	dg11.getH1F("h10").fill(refP);
        	
        	int [][][] pid = new int[6][3][6]; pid = ev.getECALPID(2);

    		dg0.getH1F("h211g").fill(pid[1][0][1]);		
    		dg0.getH1F("h211n").fill(pid[0][0][1]);
    		dg0.getH1F("h221g").fill(pid[1][1][1]);		
    		dg0.getH1F("h221n").fill(pid[0][1][1]);		
    		dg0.getH1F("h231g").fill(pid[1][2][1]);		
    		dg0.getH1F("h231n").fill(pid[0][2][1]);		
    		
        	npp=0;
        	
        	
        	
        	if(np.size()>0) { // 1 or more neutral particles id=22,2112
        		 
        		double sf = SamplingFractions.getMean(22, np.get(0), eb.ccdb);
        		double en1= np.get(0).getEnergy(DetectorType.ECAL)/sf;             		
        		double dee=(refP-en1)/refP;
        		
        		dg11.getH1F("h11").fill(refP); n++;
        		
    	    	int n1=0,n4=0,n7=0,sum=0;
    	    	s1.clear(); s4.clear(); s7.clear();   	    	
    	    	double e1=0, e4=0, e7=0, e1p=0, e4p=0, e7p=0;
    	    	
        	    for (DetectorParticle phot : np) {
        	    	double p1 = phot.getEnergy(DetectorType.ECAL);
        	    	if(p1>pthresh) {
        	    	npp++;
//        	    	System.out.println(refP+" "+p1/SamplingFractions.getMean(22, phot, eb.ccdb));
        	    	for (DetectorResponse dr : phot.getDetectorResponses()) {
        	    		int lay = dr.getDescriptor().getLayer();
        	    		double e = dr.getEnergy(); double t = dr.getTime(); 
        	    		int h1 = dr.getHitIndex(); int h2 = dr.getAssociation(); float pat = (float) dr.getPath();       	    		
        	    		if(lay==1) {n1++; e1p+=e; s1.add(dr.getStatus()); if(n1==1) {e1=e; dg11.getH1F("h210").fill(refP);}}
        	    		if(lay==4) {n4++; e4p+=e; s4.add(dr.getStatus()); if(n4==1) {e4=e; dg11.getH1F("h220").fill(refP);}}
        	    		if(lay==7) {n7++; e7p+=e; s7.add(dr.getStatus()); if(n7==1) {e7=e; dg11.getH1F("h230").fill(refP);}}
//        	    		System.out.println(getEventNumber()+" "+n1+" "+n4+" "+n7+" "+lay+" "+e+" "+t+" "+pat+" "+h1+" "+h2+" "+stat);
        	    	}
        	    	}
        	    }
        	    
        		if(n1==1)       dg0.getH2F("dee1").fill(refP,dee);
        		if(n1==1&&n4>1) dg0.getH2F("dee2").fill(refP,dee);
        		if(n1==1&&n4>2) dg0.getH2F("dee3").fill(refP,dee);   

        	    dg0.getH1F("h211").fill(n1); for (Integer i : s1) dg0.getH2F("st1").fill(n1,i);
        	    dg0.getH1F("h221").fill(n4); for (Integer i : s4) dg0.getH2F("st4").fill(n4,i);
        	    dg0.getH1F("h231").fill(n7); for (Integer i : s7) dg0.getH2F("st7").fill(n7,i);
        	    
        	    dg11.getH1F("h21").fill(refP,n1);dg11.getH1F("h22").fill(refP,n4); dg11.getH1F("h23").fill(refP,n7);
        	    
                if(n1>1) dg0.getH1F("h212").fill(e1/e1p); dg0.getH2F("set1").fill(e1/e1p,n1);
                if(n4>1) dg0.getH1F("h222").fill(e4/e4p); dg0.getH2F("set4").fill(e4/e4p,n4);
                if(n7>1) dg0.getH1F("h232").fill(e7/e7p); dg0.getH2F("set7").fill(e7/e7p,n7);
        	}
        	if(npp==1) { // only 1 neutral particle id=22,2112            		
        		dg11.getH1F("h12").fill(refP); ng++;
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
        	    	if(sum>=100) dg11.getH1F("h13").fill(refP);
        	    	if(sum>=140) dg11.getH1F("h14").fill(refP);
        	    	if(sum==147) dg11.getH1F("h15").fill(refP);
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
		int run = getRunNumber();
		DataGroup  dg0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("MAIN"),run);
		DataGroup dg10 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg11 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		dg10.getH1F("ef11").add(H1F.divide(dg11.getH1F("h11"), dg11.getH1F("h10")));	
		dg10.getH1F("ef13").add(H1F.divide(dg11.getH1F("h210"),dg11.getH1F("h10")));	
		dg10.getH1F("ef14").add(H1F.divide(dg11.getH1F("h220"),dg11.getH1F("h10")));	
		dg10.getH1F("ef15").add(H1F.divide(dg11.getH1F("h230"),dg11.getH1F("h10")));	
		dg10.getH1F("ef21").add(H1F.divide(dg11.getH1F("h11"), dg11.getH1F("h10")));	
		dg10.getH1F("ef22").add(H1F.divide(dg11.getH1F("h12"), dg11.getH1F("h10")));	
		dg10.getH1F("ef23").add(H1F.divide(dg11.getH1F("h13"), dg11.getH1F("h10")));	
		dg10.getH1F("ef24").add(H1F.divide(dg11.getH1F("h14"), dg11.getH1F("h10")));	
		dg10.getH1F("ef25").add(H1F.divide(dg11.getH1F("h15"), dg11.getH1F("h10")));	
		dg0.getH1F("hr21").add(H1F.divide( dg11.getH1F("h21"), dg11.getH1F("h210")));	
		dg0.getH1F("hr22").add(H1F.divide( dg11.getH1F("h22"), dg11.getH1F("h220")));	
		dg0.getH1F("hr23").add(H1F.divide( dg11.getH1F("h23"), dg11.getH1F("h230")));	
    }
    
    public void plot(String tabname) {      	
    	int index = getDetectorTabNames().indexOf(tabname);
    	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));   	
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
