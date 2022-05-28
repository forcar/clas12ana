package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.service.ec.ECCluster;
import org.clas.service.ec.ECPeak;
import org.clas.service.ec.ECStrip;
import org.clas.tools.EBMCEngine;
import org.clas.tools.Event;
import org.clas.tools.ParallelSliceFitter;
import org.clas.tools.SFFunction;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc1 extends DetectorMonitor {
	
	Event            ev = new Event();
	EBMCEngine    ebmce = new EBMCEngine();
	ProcessECEngine pec = new ProcessECEngine();
	
	List<Float>           GEN = new ArrayList<Float>();
	List<Float>           REC = new ArrayList<Float>();   
    HipoDataSync       writer = null;		
	List<Particle>       phot = new ArrayList<Particle>();
	List<Particle>       neut = new ArrayList<Particle>();
    SFFit                  sf = null;
    
	String tit = null;
	double ethresh = 0.01;
	public float e1=0,e2=2.5f,p2=e2*4;
	static float refP, refTH;
	static int trSEC=5, trPID=-211, mcSEC=2, mcPID= 22;
   
    public ECmc1(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("CLUSTERS","RECGEN","EFFICIENCY","SF");

        this.use123Buttons(true);
        this.useZSliderPane(true);
        useECEnginePane(true);
        useEBCCDB(true);
        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECmc1.localinit()");          
        initEBMCE();
        tl.setFitData(Fits);
        sf = new SFFit(e1,e2);
    }
    
    public void initEBMCE() {        
        System.out.println("ECmc1.initEBMCE()");          
        ebmce.getCCDB(11);
        ebmce.setGeom("2.5");
        ebmce.setGoodPhotons(12);
        ebmce.setMCpid(mcPID);
        ebmce.setMCsec(mcSEC);
        ebmce.isMC = true;    	
    }
    
    public void localclear() {
    	System.out.println("ECmc1.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    } 
    
    public void openOutput(String file) {
    	System.out.println("ECmc1:openOutput("+file+")");
    	writer = new HipoDataSync();
        writer.open(file);    	
    } 
    
    @Override
    public void createHistos(int run) {  
	    System.out.println("ECmc1:createHistos("+run+")");
        if(dumpFiles) openOutput("/Users/colesmith/CLAS12ANA/ECmc1/photon_mc1.hipo");
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;    	
    	createCLUSTERS(0);
    	createRECGEN(0);
    	createRECGEN(1);
       	createSF(0);
       	createEFFICIENCY(0);
       	createEFFICIENCY(1);
    }
        
    public void createCLUSTERS(int st) {
    	       
    	switch (st) {        
        case 0:    
    		dgm.add("CLUSTERS",6,5,0,st,getRunNumber());    
    		dgm.makeH1("h211", 5,0.5,5.5,-1,"BLK: N#gamma > 0 RED: N#gamma(PC>0)" ,"PCAL Clusters",1,2,"100"); dgm.cc("h211", true,false,0,0,0,0);  
    		dgm.makeH1("h211t",5,0.5,5.5,-2,                                    "","PCAL Clusters",1,0,"100");                              
    		dgm.makeH1("h221", 5,0.5,5.5,-1,"BLK: N#gamma > 0 YEL: N#gamma(ECI>0)","ECIN Clusters",1,25,"100");dgm.cc("h221", true,false,0,0,0,0);
    		dgm.makeH1("h221t",5,0.5,5.5,-2,                                    "","ECIN Clusters",1,0,"100");                             
    		dgm.makeH1("h231", 5,0.5,5.5,-1,"BLK: N#gamma > 0 BLU: N#gamma(ECO>0)","ECOU Clusters",1,24,"100");dgm.cc("h231", true,false,0,0,0,0);
    		dgm.makeH1("h231t",5,0.5,5.5,-2,                                    "","ECOU Clusters",1,0,"100");                              

    		dgm.makeH2("st1",  5,0.5,5.5,7,0,7,-1,     "","PCAL Clusters","STATUS");                          dgm.cc("st1", false,true,0,0,0,0);
    		dgm.makeH2("st4",  5,0.5,5.5,7,0,7,-1,     "","ECIN Clusters","STATUS");                          dgm.cc("st4", false,true,0,0,0,0);
    		dgm.makeH2("st7",  5,0.5,5.5,7,0,7,-1,     "","ECOU Clusters","STATUS");                          dgm.cc("st7", false,true,0,0,0,0);                     

    		dgm.makeH2("set1", 50,0,1,4,1.5,5.5,-1,    "","PCAL Energy Deficit","PCAL Clusters");             dgm.cc("set1",false,true,0,0,0,0); 
    		dgm.makeH2("set4", 50,0,1,4,1.5,5.5,-1,    "","ECIN Energy Deficit","ECIN Clusters");             dgm.cc("set4",false,true,0,0,0,0); 
    		dgm.makeH2("set7", 50,0,1,4,1.5,5.5,-1,    "","ECOU Energy Deficit","ECOU Clusters");             dgm.cc("set7",false,true,0,0,0,0);

    		dgm.makeH2("sg1",  5,0.5,5.5,5,0.5,5.5,-1,     "","PCAL Clusters","N#gamma");                     dgm.cc("sg1", false,true,0,0,0,0);
    		dgm.makeH2("sg4",  5,0.5,5.5,5,0.5,5.5,-1,     "","ECIN Clusters","N#gamma");                     dgm.cc("sg4", false,true,0,0,0,0);
    		dgm.makeH2("sg7",  5,0.5,5.5,5,0.5,5.5,-1,     "","ECOU Clusters","N#gamma");                     dgm.cc("sg7", false,true,0,0,0,0);                     

    		dgm.makeH1("h212", 50,0.,1,-1,             "","PCAL Energy Deficit",1,2);
    		dgm.makeH1("h222", 50,0.,1,-1,             "","ECIN Energy Deficit",1,25);
    		dgm.makeH1("h232", 50,0.,1,-1,             "","ECOU Energy Deficit",1,24);

    		dgm.makeH1("hr21", 50,0,p2,-1,            "","Photon Energy (GeV)","Avg. No. PCAL Clusters",1,2);  dgm.cc("hr21",false,false,1,1.1f,0,0);
    		dgm.makeH1("hr22", 50,0,p2,-1,            "","Photon Energy (GeV)","Avg. No. ECIN Clusters",1,25); dgm.cc("hr22",false,false,1,1.1f,0,0);
    		dgm.makeH1("hr23", 50,0,p2,-1,            "","Photon Energy (GeV)","Avg. No. ECOU Clusters",1,24); dgm.cc("hr23",false,false,1,1.1f,0,0);            

            int nb=100; float x1=-180, x2=-20, y1=90, y2=250;
        	dgm.makeH2("cxy1", nb,x1,x2,nb,y1,y2, 0, "PCAL RAW","X ", "Y ");
        	dgm.makeH2("cxy4", nb,x1,x2,nb,y1,y2, 0, "ECIN RAW","X ", "Y ");
        	dgm.makeH2("cxy7", nb,x1,x2,nb,y1,y2, 0, "ECOU RAW","X ", "Y "); 
        	dgm.makeH1R("pef1",3,1, 4, -1,"","PCAL U,V,W","Peak/Cluster",1,2,"");
        	dgm.makeH1R("pef4",3,1, 4, -1,"","ECIN U,V,W","Peak/Cluster",1,25,"");
        	dgm.makeH1R("pef7",3,1, 4, -1,"","ECOU U,V,W","Peak/Cluster",1,24,"");
        	dgm.makeH1R("peft",9,1,10, -1,"","ECAL LAYER","Peak/Photon",1,0,"");
        	dgm.makeH2("peakE1", 9,1,10,50,0,0.6,  -1,"","ECAL LAYER","Peak/Cluster");
        	dgm.makeH2("peakE2",10,0,p2,270,1,10,  -1,"","Photon Energy (GeV)","Peak/Cluster");
        	dgm.makeH2("peakE3", 9,1,10,50,0.5,1.5,-1,"UV,UW,VW","ECAL LAYER","Peak/Peak");
        	dgm.makeH2("peakE4",10,0,p2,270,1,10,  -1,"UV,UW,VW","Photon Energy (GeV)","Peak/Peak");
    	}
    }
    
    public void createRECGEN(int st) {
    	
    	switch (st) {        
        case 0: 
    		dgm.add("RECGEN",5,4,0,st,getRunNumber());    
        	dgm.makeH2("dee0", 50,0,p2,20,-0.4,0.4,0,"N1>0",     "Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee0",false,true,0,0,0,0);
        	dgm.makeH2("dee1", 50,0,p2,20,-0.4,0.4,0,"N1=1",     "Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee1",false,true,0,0,0,0);
        	dgm.makeH2("dee2", 40,0,p2,15,-0.4,0.4,0,"N1=1 N4=2","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee2",false,true,0,0,0,0); 
        	dgm.makeH2("dee3", 30,0,p2,10,-0.4,0.4,0,"N1=1 N4=3","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee3",false,true,0,0,0,0);
        	dgm.makeH2("dee4", 30,0,p2,10,-0.4,0.4,0,"N1=1 N4=4","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee4",false,true,0,0,0,0);
        	dgm.makeH2("det0", 50,0,p2,20,-0.6,0.6,0,"N1>0",     "Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det0",false,true,0,0,0,0);
        	dgm.makeH2("det1", 50,0,p2,20,-0.6,0.6,0,"N1=1",     "Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det1",false,true,0,0,0,0);
        	dgm.makeH2("det2", 50,0,p2,20,-0.6,0.6,0,"N1=1 N4=2","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det2",false,true,0,0,0,0);
        	dgm.makeH2("det3", 50,0,p2,20,-0.6,0.6,0,"N1=1 N4=3","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det3",false,true,0,0,0,0);
        	dgm.makeH2("det4", 50,0,p2,20,-0.6,0.6,0,"N1=1 N4=4","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det4",false,true,0,0,0,0);
        	dgm.makeH2("dep0", 50,0,p2,20,-0.6,0.6,0,"N1>0",     "Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep0",false,true,0,0,0,0);
        	dgm.makeH2("dep1", 50,0,p2,20,-0.6,0.6,0,"N1=1",     "Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep1",false,true,0,0,0,0);
        	dgm.makeH2("dep2", 50,0,p2,20,-0.6,0.6,0,"N1=1 N4=2","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep2",false,true,0,0,0,0);
        	dgm.makeH2("dep3", 50,0,p2,20,-0.6,0.6,0,"N1=1 N4=3","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep3",false,true,0,0,0,0);
        	dgm.makeH2("dep4", 50,0,p2,20,-0.6,0.6,0,"N1=1 N4=4","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep4",false,true,0,0,0,0);
        	dgm.makeH2("deb0", 50,0,p2,20,-0.05,0.05,0,"N1>0",     "Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb0",false,true,0,0,0,0);
        	dgm.makeH2("deb1", 50,0,p2,20,-0.05,0.05,0,"N1=1",     "Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb1",false,true,0,0,0,0);
        	dgm.makeH2("deb2", 50,0,p2,20,-0.05,0.05,0,"N1=1 N4=2","Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb2",false,true,0,0,0,0);
        	dgm.makeH2("deb3", 50,0,p2,20,-0.05,0.05,0,"N1=1 N4=3","Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb3",false,true,0,0,0,0);
        	dgm.makeH2("deb4", 50,0,p2,20,-0.05,0.05,0,"N1=1 N4=4","Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb4",false,true,0,0,0,0);
        	break;
        case 1: 
    		dgm.add("RECGEN",5,4,0,st,getRunNumber());    
        	dgm.makeH2("dte0", 50,5,26,20,-0.4,0.4,0,"N1>0",     "Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte0",false,true,0,0,0,0);
        	dgm.makeH2("dte1", 50,5,26,20,-0.4,0.4,0,"N1=1",     "Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte1",false,true,0,0,0,0);
        	dgm.makeH2("dte2", 40,5,26,15,-0.4,0.4,0,"N1=1 N4=2","Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte2",false,true,0,0,0,0); 
        	dgm.makeH2("dte3", 30,5,26,10,-0.4,0.4,0,"N1=1 N4=3","Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte3",false,true,0,0,0,0);
        	dgm.makeH2("dte4", 30,5,26,10,-0.4,0.4,0,"N1=1 N4=4","Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte4",false,true,0,0,0,0);
        	dgm.makeH2("dtt0", 50,5,26,20,-0.6,0.6,0,"N1>0",     "Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt0",false,true,0,0,0,0);
        	dgm.makeH2("dtt1", 50,5,26,20,-0.6,0.6,0,"N1=1",     "Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt1",false,true,0,0,0,0);
        	dgm.makeH2("dtt2", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=2","Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt2",false,true,0,0,0,0);
        	dgm.makeH2("dtt3", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=3","Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt3",false,true,0,0,0,0);
        	dgm.makeH2("dtt4", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=4","Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt4",false,true,0,0,0,0);
           	dgm.makeH2("dtp0", 50,5,26,20,-0.6,0.6,0,"N1>0",     "Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp0",false,true,0,0,0,0);
           	dgm.makeH2("dtp1", 50,5,26,20,-0.6,0.6,0,"N1=1",     "Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp1",false,true,0,0,0,0);
        	dgm.makeH2("dtp2", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=2","Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp2",false,true,0,0,0,0);
        	dgm.makeH2("dtp3", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=3","Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp3",false,true,0,0,0,0);
        	dgm.makeH2("dtp4", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=4","Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp4",false,true,0,0,0,0);
        	dgm.makeH2("dtb0", 50,5,26,20,-0.05,0.05,0,"N1>0",     "Photon Theta (deg)","#Delta#beta");      dgm.cc("dtb0",false,true,0,0,0,0);
        	dgm.makeH2("dtb1", 50,5,26,20,-0.05,0.05,0,"N1=1",     "Photon Theta (deg)","#Delta#beta");      dgm.cc("dtb1",false,true,0,0,0,0);
        	dgm.makeH2("dtb2", 50,5,26,20,-0.05,0.05,0,"N1=1 N4=2","Photon Theta (deg)","#Delta#beta");      dgm.cc("dtb2",false,true,0,0,0,0);
        	dgm.makeH2("dtb3", 50,5,26,20,-0.05,0.05,0,"N1=1 N4=3","Photon Theta (deg)","#Delta#beta");      dgm.cc("dtb3",false,true,0,0,0,0);
        	dgm.makeH2("dtb4", 50,5,26,20,-0.05,0.05,0,"N1=1 N4=4","Photon Theta (deg)","#Delta#beta");      dgm.cc("dtb4",false,true,0,0,0,0);
    	}
    }
    	
    public void createEFFICIENCY(int st) {
        
    	switch (st) {        
        case 0:         	
    		dgm.add("EFFICIENCY",2,2,0,st,getRunNumber()); 
    		tit = "N#gamma>0, thresh = "+ethresh*1e3+" MeV R:PC>0 Y:ECI>0 B:ECO>0";
        	dgm.makeH1("ef11",50,0,p2,-1,tit,              "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef13",50,0,p2,-2,"n>0 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef14",50,0,p2,-2,"n>0 layer 1,4",  "Photon Energy (GeV)",1,25);
        	dgm.makeH1("ef15",50,0,p2,-2,"n>0 layer 1,4,7","Photon Energy (GeV)",1,24);
        	tit = "N#gamma=1, thresh = "+ethresh*1e3+" MeV  MeV R:PC=1 Y:ECI=1 B:ECO=1";
        	dgm.makeH1("ef21",50,0,p2,-1,tit,              "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef22",50,0,p2,-2,"n=1 any layer",  "Photon Energy (GeV)",1,1);
        	dgm.makeH1("ef23",50,0,p2,-2,"n=1 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef24",50,0,p2,-2,"n=1 layer 1,4",  "Photon Energy (GeV)",1,25);
        	dgm.makeH1("ef25",50,0,p2,-2,"n=1 layer 1,4,7","Photon Energy (GeV)",1,24);        	
        	break;
        case 1:            
    		dgm.add("EFFICIENCY",3,5,0,st,getRunNumber()); tit = "n>0 any layer, thresh = "+ethresh+" MeV";
        	dgm.makeH1("h10", 50,0,p2,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h11", 50,0,p2,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h12", 50,0,p2,-1,"PCAL","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h13", 50,0,p2,-1,"PCAL","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h14", 50,0,p2,-1,"ECIN","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h15", 50,0,p2,-1,"ECOU","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h21", 50,0,p2,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h22", 50,0,p2,-1,"ECIN","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h23", 50,0,p2,-1,"ECOU","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h210",50,0,p2,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h220",50,0,p2,-1,"ECIN","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h230",50,0,p2,-1,"ECOU","Photon Energy (GeV)",1,5);
    	}
    }
        
    public void createSF(int st) {
    	
    	switch (st) {
    	case 0: 
    		dgm.add("SF", 2, 2, 0, st, getRunNumber());
    		dgm.makeH2("sf1",80,e1,e2,50,0.15,0.35,-1,"N#gamma>0","Measured Photon Energy (GeV)","Sampling Fraction");  dgm.cc("sf1",false,true,0,0,0,0);
    		dgm.makeH2("sf2",80,e1,e2,50,0.15,0.35,-1,"N#gamma=1","Measured Photon Energy (GeV)","Sampling Fraction");  dgm.cc("sf2",false,true,0,0,0,0);
    		dgm.makeH2("sf3",80,5, 26,50,0.15,0.35,-1,"N#gamma>0","Photon Theta  (deg)","Sampling Fraction");           dgm.cc("sf3",false,true,0,0,0,0);
    		dgm.makeH2R("sf4",20,e1,e2,40,5.5,25.5,-1,"SF: N#gamma>0","Measured Photon Energy (GeV)","Photon Theta (deg)"); dgm.cc("sf4",false,false,0,0,0.18f,0.265f);		
    	}
    }
    
    @Override       
    public void plotHistos(int run) {    	
    	if(!histosExist) return;
    	plotSummary(run);
    	plotAnalysis(run);
    }
    
    public void plotSummary(int run) {
    	setRunNumber(run);
    	plotMCHistos();
    }

    public void plotMCHistos() {
        plot("CLUSTERS");
        plot("RECGEN");
        plot("EFFICIENCY");  
        plot("SF");
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;
     	sf.plot("SF");
    }
        
    @Override
    public void plotEvent(DataEvent de) {
    	System.out.println(getDetectorName()+".plotEvent() ");  
        analyze();         
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
    	geteff();
    	sf.fit("sf1","sf2");
    	if(dumpFiles) writer.close();
    	isAnalyzeDone = true;
    }
    
    public List<Float> getkin (List<Particle> list) {
    	List<Float> out = new ArrayList<Float>();
    	int n=0;
    	for (Particle p : list) {
		    out.add(n++,(float) p.e()); 
		    out.add(n++,(float) Math.toDegrees(p.theta())); 
		    out.add(n++,(float) Math.toDegrees(p.phi()));
		    if(p.hasProperty("beta")) out.add(n++,(float) p.getProperty("beta"));
    	}
		return out;
    }
    
    public class ProcessECEngine {
    	
    	int s,l;
    	double e,em,x,y;
    	ECCluster c = null; 
    	public int[] status = {-1,1};
    	
    	public ProcessECEngine() {
    		
    	}

        public void processEnergy(DetectorParticle dp) {

        	List<ECCluster> cl = eng.engine.getClusters();
    		em = (float) dp.getEnergy(DetectorType.ECAL);
        	for(DetectorResponse dr : dp.getDetectorResponses()) {
        		c = cl.get(dr.getHitIndex());
        		s = c.clusterPeaks.get(0).getDescriptor().getSector();
        		l = c.clusterPeaks.get(0).getDescriptor().getLayer();
        		e = c.getEnergy();
        		x = c.getHitPosition().x();
        		y = c.getHitPosition().y();
        		if(s==mcSEC && c.getStatus()>status[0]&&c.getStatus()<status[1]) {
        			if(l==1) {dgm.fill("cxy"+l,-x,y); for (int i=0;i<3;i++) processHist(i,l);}
        			if(l==4) {dgm.fill("cxy"+l,-x,y); for (int i=0;i<3;i++) processHist(i,l);}
        			if(l==7) {dgm.fill("cxy"+l,-x,y); for (int i=0;i<3;i++) processHist(i,l);}
        		}
        	}    	
        } 
        
        public void processHist(int i, int l) { 
        	double rat=0;
        	double ec = c.getEnergy(i);
        	double rat1 = ec/e, rat2=ec/em;
        	dgm.fill("pef"+l, i+1,rat1);
        	dgm.fill("peft",  i+l,rat2);
        	dgm.fill("peakE1",i+l,rat1);           	
        	if(rat1<0.6) dgm.fill("peakE2",refP,i+l+rat1/0.6);        	
        	if(i==0) {rat = c.getEnergy(0)/c.getEnergy(1); dgm.fill("peakE3",i+l,rat);}
        	if(i==0 && rat>0.5 && rat<1.5) dgm.fill("peakE4", refP, i+l+rat-0.5);
        	if(i==1) {rat = c.getEnergy(0)/c.getEnergy(2); dgm.fill("peakE3",i+l,rat);}
        	if(i==1 && rat>0.5 && rat<1.5) dgm.fill("peakE4", refP, i+l+rat-0.5);
        	if(i==2) {rat = c.getEnergy(1)/c.getEnergy(2); dgm.fill("peakE3",i+l,rat);}
        	if(i==2 && rat>0.5 && rat<1.5) dgm.fill("peakE4", refP, i+l+rat-0.5);
        	
        }
    		   
    }

    @Override
    public void processEvent(DataEvent de) {
    	
    	GEN.clear();  REC.clear();
    	
    	DetectorParticle p1 = new DetectorParticle();
    	DetectorParticle p2 = new DetectorParticle();  
    	
    	List<Integer> s1 = new ArrayList<Integer>();
    	List<Integer> s4 = new ArrayList<Integer>();
    	List<Integer> s7 = new ArrayList<Integer>();
    	
    	int n1=0, n4=0, n7=0, npart=0;
     	double e1=0, e4=0, e7=0, e1p=0, e4p=0, e7p=0;     
     	
		boolean goodev = ebmce.readMC(de) && ebmce.pmc.size()==1;
		
        if (goodev) { 
        	
        	GEN = getkin(ebmce.pmc); refP = (float) GEN.get(0); refTH = (float) GEN.get(1); 
        	
//        	processECEngine(de);        	
	    	if(dropBanks) dropBanks(de);  //drop ECAL banks and re-run ECEngine  
	    	
        	if(!ebmce.processDataEvent(de)) return;
        	
        	float stt = ebmce.starttime;  
            
            List<DetectorParticle> par = ebmce.eb.getEvent().getParticles(); //REC::Particle
        	List<DetectorResponse> cal = ebmce.eb.getEvent().getCalorimeterResponseList();  //REC::Calorimeter
        	
        	if(dumpFiles) {ebmce.getRECBanks(de,ebmce.eb); writer.writeEvent(de);}
        	
        	int trsec = -1; 
        	for (DetectorParticle dp : par) { //find trigger particle and sector
        		int pid = dp.getPid(); int sec = dp.getSector(DetectorType.ECAL);
        		if(trsec==-1 && sec==trSEC && pid==trPID) trsec=sec;
        	}
        	
        	if(trsec==-1) return; // reject events with missing trigger particle  
        	
        	if (dbgAnalyzer) {
        		System.out.println(" ");
        		System.out.println(getEventNumber()+" "+cal.size()+" "+par.size());
        	
        		for (DetectorResponse drr : cal) {
        			CalorimeterResponse dr = (CalorimeterResponse) drr;
        			System.out.println(dr.getAssociation()+" "
        		                      +dr.getDescriptor().getType()  +" "
        		                      +dr.getDescriptor().getSector()+" "
        				              +dr.getDescriptor().getLayer() +" "
        				              +dr.getStatus()                +" "
        				              +dr.getEnergy()                +" "
        				              +dr.getHitIndex()              +" "
        	                          +par.get(dr.getAssociation()).getPid());
        		}
        	
        		int nnn=0;
        		for (DetectorParticle dp : par) {
        			System.out.println("Particle "+nnn+"  "+dp.getSector(DetectorType.ECAL)+" "
        	                                               +dp.getEnergy(DetectorType.ECAL)+" "
        	                                               +dp.getPid()                    +" " 
        	                                               +dp.getBeta()                   +" "
        	                                               +ebmce.hasStartTime);nnn++;
        	        for (DetectorResponse dr : dp.getDetectorResponses()) {
        	        	System.out.println(dr.getAssociation()           +" "
        		                          +dr.getDescriptor().getType()  +" "
            		                      +dr.getDescriptor().getSector()+" "
              			                  +dr.getDescriptor().getLayer() +" "
              			                  +dr.getStatus()                +" "
              			                  +dr.getEnergy()                +" "
        	        					  +dr.getHitIndex());
        	        }
        		}	
        	}
        	       	
        	List<Particle> plist = new ArrayList<Particle>(); 
      	
        	int dpind = -1, cind=-1;
        	
        	for (DetectorParticle dp : par) { // make list of neutral Particle objects 
        		dpind++;
    		    if(dp.getSector(DetectorType.ECAL)==mcSEC && dp.getPid()==mcPID) { npart++; cind=dpind;
			    if(ebmce.hasStartTime && dp.getPid()==2112) {// / for photon MC you may want to recover these clusters
	 				double e = dp.getEnergy(DetectorType.ECAL)/ebmce.getSF(dp); 		
			    	Vector3D vec = new Vector3D() ; vec.copy(dp.getHit(DetectorType.ECAL).getPosition()); vec.unit(); 			    		
 			    	dp.vector().add(new Vector3(e*vec.x(),e*vec.y(),e*vec.z())); //track energy for neutrals in DetectorParticle
 			    	dp.setPid(mcPID);
 			    }
		    	//Particle momentum calculated in EBAnalyzer.assignNeutralMomenta 
			    //p.vector().setMag(p.getEnergy(DetectorType.ECAL)/SF 
			    Particle p = dp.getPhysicsParticle(mcPID); p.setProperty("beta",dp.getBeta()); plist.add(p);
		        }
 			}
        	
        	if(npart==1) pec.processEnergy(par.get(cind));
//        	if(npart==2) {pec.status[0]=0; pec.status[1]=4; pec.processEnergy(par.get(cind));}
        	
        	int ipp=0, npp=0;
        	
        	for (DetectorParticle dp : par) { 
        		if(dp.getSector(DetectorType.ECAL)==mcSEC && dp.getPid()==mcPID) { 
        			float em = (float) (dp.getEnergy(DetectorType.ECAL)); 
        			float  e = (float) (em/ebmce.getSF(dp)); 				  
        			if(e>ethresh) {
        				dgm.fill("sf1",em,em/refP); dgm.fill("sf3",refTH,em/refP);
        				for (DetectorResponse dr : cal){ CalorimeterResponse r = (CalorimeterResponse) dr;          					
        					if (r.getAssociation()==ipp) { 
        						int lay = r.getDescriptor().getLayer(); 	
        						if(lay==1) {n1++; e1p+=e; s1.add(r.getStatus()); if(n1==1) {e1=e; dgm.fill("h210",refP);}}
        						if(lay==4) {n4++; e4p+=e; s4.add(r.getStatus()); if(n4==1) {e4=e; dgm.fill("h220",refP);}}
        						if(lay==7) {n7++; e7p+=e; s7.add(r.getStatus()); if(n7==1) {e7=e; dgm.fill("h230",refP);}}
        					}
        				}
        				npp++;
        			}		        		
        		}
        		ipp++;
 			}  
       	    
        	dgm.fill("h10",refP);	
        	        	
        	if(npart>0) { // 1 or more neutral particles 
         		
           		REC = getkin(plist);
        		double  delE1 = (REC.get(0) - GEN.get(0))/refP;
        		double delTH1 = (REC.get(1) - GEN.get(1));
        		double delPH1 = (REC.get(2) - GEN.get(2));
        		double delBET = (REC.get(3) - 1);
        		
        		dgm.fill("h11",refP);       	   
        	    
        		//RECGEN.0
        		if(n1>0)        {dgm.fill("dee0",refP,delE1);dgm.fill("det0",refP,delTH1); dgm.fill("dep0",refP,delPH1);dgm.fill("deb0",refP,delBET);}
        		if(n1==1)       {dgm.fill("dee1",refP,delE1);dgm.fill("det1",refP,delTH1); dgm.fill("dep1",refP,delPH1);dgm.fill("deb1",refP,delBET);}
        		if(n1==1&&n4>1) {dgm.fill("dee2",refP,delE1);dgm.fill("det2",refP,delTH1); dgm.fill("dep2",refP,delPH1);dgm.fill("deb2",refP,delBET);}
        		if(n1==1&&n4>2) {dgm.fill("dee3",refP,delE1);dgm.fill("det3",refP,delTH1); dgm.fill("dep3",refP,delPH1);dgm.fill("deb3",refP,delBET);}   
        		if(n1==1&&n4>3) {dgm.fill("dee4",refP,delE1);dgm.fill("det4",refP,delTH1); dgm.fill("dep4",refP,delPH1);dgm.fill("deb4",refP,delBET);} 
        		
        		//RECGEN.1
        		if(refP>0) {
            	if(n1>0)        {dgm.fill("dte0",refTH,delE1);dgm.fill("dtt0",refTH,delTH1); dgm.fill("dtp0",refTH,delPH1);dgm.fill("dtb0",refTH,delBET);}
            	if(n1==1)       {dgm.fill("dte1",refTH,delE1);dgm.fill("dtt1",refTH,delTH1); dgm.fill("dtp1",refTH,delPH1);dgm.fill("dtb1",refTH,delBET);}
        		if(n1==1&&n4>1) {dgm.fill("dte2",refTH,delE1);dgm.fill("dtt2",refTH,delTH1); dgm.fill("dtp2",refTH,delPH1);dgm.fill("dtb2",refTH,delBET);}
        		if(n1==1&&n4>2) {dgm.fill("dte3",refTH,delE1);dgm.fill("dtt3",refTH,delTH1); dgm.fill("dtp3",refTH,delPH1);dgm.fill("dtb3",refTH,delBET);}   
        		if(n1==1&&n4>3) {dgm.fill("dte4",refTH,delE1);dgm.fill("dtt4",refTH,delTH1); dgm.fill("dtp4",refTH,delPH1);dgm.fill("dtb4",refTH,delBET);}  
        		}
        		//CLUSTERS.0.1
        		dgm.fill("h211",n1); if(n1>0) dgm.fill("h211t", npart); 
        		dgm.fill("h221",n4); if(n4>0) dgm.fill("h221t", npart);  
        		dgm.fill("h231",n7); if(n7>0) dgm.fill("h231t", npart);
        		//CLUSTERS.0.2
        		for (Integer i : s1) dgm.fill("st1",n1,i); 
       		    for (Integer i : s4) dgm.fill("st4",n4,i); 
        		for (Integer i : s7) dgm.fill("st7",n7,i); 
        		//CLUSTERS.0.4
        		dgm.fill("sg1",n1,npart); 
       		    dgm.fill("sg4",n4,npart); 
        		dgm.fill("sg7",n7,npart); 
        		//CLUSTERS.0.6
        		dgm.fill("h21",refP,n1); 
        		dgm.fill("h22",refP,n4); 
        		dgm.fill("h23",refP,n7);        		
        		//CLUSTERS.0.4
        		dgm.fill("set1",e1/e1p,n1);
        		dgm.fill("set4",e4/e4p,n4); 
        		dgm.fill("set7",e7/e7p,n7);
        		//CLUSTERS.0.5
        		if(n1>1) dgm.fill("h212",e1/e1p); 
        		if(n4>1) dgm.fill("h222",e4/e4p); 
        		if(n7>1) dgm.fill("h232",e7/e7p);
        	}
        	
        	if(npp==1) { // only 1 neutral particle above threshold           		
    	    	dgm.fill("h12",refP); 
    	    	for (DetectorParticle dp : par) { 
    	    		float em = (float) dp.getEnergy(DetectorType.ECAL);
	        		double e = em/ebmce.getSF(dp); 
        			int sum=0;
        			if(dp.getSector(DetectorType.ECAL)==mcSEC && e>ethresh && dp.getPid()==mcPID) {         				
	        			for (DetectorResponse dr : dp.getDetectorResponses()) {
	        				if(dr.getDescriptor().getType()==DetectorType.ECAL) {	        					
		        				int lay = dr.getDescriptor().getLayer();
	        					if(lay==1) sum+= 100;
	        					if(lay==4) sum+=  40;
	        					if(lay==7) sum+=   7;
	        				}
	        			}	        			
	        			if(sum>=100){dgm.fill("h13",refP); dgm.fill("sf2", em, em/refP); dgm.fill("sf4", em, refTH, em/refP);}
	        			if(sum>=140) dgm.fill("h14",refP);
	        			if(sum==147) dgm.fill("h15",refP);	        					        			
	        		}
        	    }
        	}
        } 
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
    	dgm.geteff("pef1");
    	dgm.geteff("pef4");
    	dgm.geteff("pef7");	
    	dgm.geteff("peft");
    	dgm.geteff("sf4");
    }
    
    public class SFFit {
    	
        float dmin=0.045f, dmax=0.1f, min=0, max=0;

        public SFFit(float e1, float e2) {
        	min = e1+dmin;
        	max = e2-dmax;
        }
    	
        public void fit(String...  h2) {
          	
           	int n=0; 
                   
        	for (String h : h2) {    
        		ParallelSliceFitter fitter = new ParallelSliceFitter(dgm.getH2F(h));
        		fitter.setBackgroundOrder(0); fitter.setMin(0.16); fitter.setMax(0.30); fitter.fitSlicesX(); 
        		GraphErrors MeanGraph = fitter.getMeanSlices();      
        		FitSummary.add(MeanGraph,n, 0, 7, getRunNumber());  
        		tl.fitData.add(fitEngine(MeanGraph,15,min,max,min,max),n++,0,5,getRunNumber());
        	}    	    
        }
        
        public void plot(String tab) {
        	
            int run = getRunNumber();
            
       	    int index = getDetectorTabNames().indexOf(tab);    	
            EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));

            SFFunction sf = new SFFunction("gsf",22,0,ebccdb,min,max); 
            sf.setLineStyle(2); sf.setLineWidth(3) ; sf.setLineColor(1); sf.setOptStat("");

       		for (int n=0; n<4; n++) {
                tl.fitData.getItem(n,0,5,run).graph.getFunction().setLineColor(1); 
                tl.fitData.getItem(n,0,5,run).graph.getFunction().setLineWidth(2);
                tl.fitData.getItem(n,0,5,run).graph.getFunction().setOptStat("1110");
                tl.fitData.getItem(n,0,5,run).graph.getFunction().setRange(min,max);
       			if (FitSummary.hasItem(n,0,7,run)) {       				
       				GraphPlot((GraphErrors)FitSummary.getItem(n,0,7,run),c,n,e1,e2,0.15f,0.35f,1,4,1,""," E/P","same"); 
       				c.draw(sf,"same");
       			} 
       		}
         	
        } 
    	
    }     
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(),getRunNumber()); //here and only here is run number used
    }  
        
    @Override
    public void timerUpdate() {  

    }


}
