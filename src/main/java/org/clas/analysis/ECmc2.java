package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.analysis.ECPart.SFFunction;
import org.clas.tools.DataGroupManager;
import org.clas.tools.EBMCEngine;
import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.CalorimeterResponse;
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
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;
import org.jlab.clas.physics.Vector3;

public class ECmc2 extends DetectorMonitor {
	
	Event          ev = new Event();
	EBMCEngine  ebmce = new EBMCEngine();
	
	List<DetectorParticle> np = new ArrayList<DetectorParticle>();    
	List<Float>           GEN = new ArrayList<Float>();
	List<Float>           REC = new ArrayList<Float>();   
	HipoDataSync       writer = null;		
	List<Particle>       phot = new ArrayList<Particle>();
	List<Particle>       neut = new ArrayList<Particle>();
	
	String tit = null;
	double ethresh = 0.3;
	
    public ECmc2(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("PHOTONS","CLUSTERS","GENREC","EFFICIENCY");

        this.use123Buttons(true);
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
        engine.setCalRun(10);   
        engine.debug = false;
//        engine.setLogParam(0);
        
        ebmce.getCCDB(10);
//        eb.setThresholds("Pizero",engine);
        ebmce.setThresholds("Pizero",engine);
        ebmce.setGeom("2.5");
        ebmce.setGoodPhotons(12);
        ebmce.setMCpid(22);
        ebmce.isMC = true;
        
        tl.setFitData(Fits);
        writer = new HipoDataSync();
        writer.open("/Users/colesmith/CLAS12ANA/ECmc2/photon_mc2.hipo");
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
    	createPHOTONS(0);
    	createCLUSTERS(0);
    	createCLUSTERS(1);
    	createGENREC(0);
    	createGENREC(1);
    	createGENREC(2);
       	createEFFICIENCY(0);
       	createEFFICIENCY(1);
    }
    
    public void createPHOTONS(int st) {
    	
    	switch (st) {
    	case 0:
    		dgm.add("PHOTONS", 4,2,0,st,getRunNumber());
    		dgm.makeH1("p00",100,0.8,1.3,-1,             "N#gamma>1","#gamma1 #beta","",1,1);
    		dgm.makeH1("p01",100,0.8,1.3,-1,             "N#gamma>1","#gamma2 #beta","",1,1);
    		dgm.makeH1("p02",100,0.8,1.3,-1,             "N#gamma=2","#gamma1 #beta","",1,1);
    		dgm.makeH1("p03",100,0.8,1.3,-1,             "N#gamma=2","#gamma2 #beta","",1,1);
    		dgm.makeH2("p04",100,0.8,1.3,16,-0.5,15.5,-1,"N#gamma>1","#gamma1 #beta","Photons");
    		dgm.makeH2("p05",100,0.8,1.3,16,-0.5,15.5,-1,"N#gamma>1","#gamma2 #beta","Photons");
    		dgm.makeH2("p06",100,0.8,1.3,50,0., 5.0,-1,  "N#gamma>1","#gamma1 #beta","Energy (GeV)");
    		dgm.makeH2("p07",100,0.8,1.3,50,0., 5.0,-1,  "N#gamma>1","#gamma2 #beta","Energy (GeV)");
    		break;
    	}
    }
    
    public void createCLUSTERS(int st) {
    	      
    	switch (st) {        
        case 0:  
        	dgm.add("CLUSTERS", 5,3,0,st,getRunNumber());
            dgm.makeH2("c00",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>1.95","PC Clusters" ,"Photons"); 
            dgm.makeH2("c01",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>1.95","ECi Clusters","Photons"); 
            dgm.makeH2("c02",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>1.95","ECo Clusters","Photons"); 
            dgm.makeH2("c03",16,-0.5,15.5, 80,0.4,1.6,1,"",          "Photons",     "REC (E1+E2)/2E");
            dgm.makeH2("c04",16,-0.5,15.5, 80,0.4,1.6,1,"OPA>1.95",  "Photons",     "REC (E1+E2)/2E");
            dgm.makeH2("c20",50,0,5.2,     16,-0.5,15.5,-1,"",         "GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c21",50,0,5.2,     16,-0.5,15.5,-1,"pc=2 ec=0","GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c22",50,0,5.2,     16,-0.5,15.5,-1,"pc=2 ec=1","GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c23",50,0,5.2,     16,-0.5,15.5,-1,"pc=2 ec>1","GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c24",50,0,5.2,     16,-0.5,15.5,-1,"pc=1 ec>1","GEN #theta12 (deg))","Photons");  
            break;
        case 1:
        	dgm.add("CLUSTERS", 5,3,0,st,getRunNumber());
            dgm.makeH2("c05",16,-0.5,15.5, 16,-0.5,15.5,-1,"",       "PC Clusters" ,"Photons"); 
            dgm.makeH2("c06",16,-0.5,15.5, 16,-0.5,15.5,-1,"",       "ECi Clusters","Photons"); 
            dgm.makeH2("c07",16,-0.5,15.5, 16,-0.5,15.5,-1,"",       "ECo Clusters","Photons");             
    	}       
    }
    
    public void createGENREC(int st) {
	  
    	switch (st) {        
    	case 0:
    		dgm.add("GENREC",6,7,0,st,getRunNumber());
                                                 tit = "pc>1 ec=0";
            dgm.makeH2("d00",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1");
            dgm.makeH2("d10",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2");
            dgm.makeH2("d02",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1"); 
            dgm.makeH2("d12",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2"); 
            dgm.makeH2("d03",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "E1/E #gamma1");  
            dgm.makeH2("d13",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "E2/E #gamma2");  
                                                 tit = "pc>1 ec=1";
            dgm.makeH2("d20",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1"); 
            dgm.makeH2("d30",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2");  
            dgm.makeH2("d22",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1");   
            dgm.makeH2("d32",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2");  
            dgm.makeH2("d23",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "E1/E #gamma1");   
            dgm.makeH2("d33",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "E2/E #gamma2");
                                                 tit = "pc>1 ec>1";
            dgm.makeH2("d40",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1"); 
            dgm.makeH2("d50",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2"); 
            dgm.makeH2("d42",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1"); 
            dgm.makeH2("d52",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2");  
            dgm.makeH2("d43",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "E1/E #gamma1");  
            dgm.makeH2("d53",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "E2/E #gamma2");              
                                                 tit = "pc=2 ec=1"; 
            dgm.makeH2("d60",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d70",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d62",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d72",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d63",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "(E1+E2)/2E");  
            dgm.makeH2("d73",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "SQRT(E1*E2)/E");  
                                                 tit = "pc>1 ec>1";
            dgm.makeH2("d80",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d90",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d82",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d92",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d83",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "(E1+E2)/2E");  
            dgm.makeH2("d93",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg))",     "SQRT(E1*E2)/E");  
                                                 tit = "pc=2 ec=2";
            dgm.makeH2("d100",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d110",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d102",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d112",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d103",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg))",     "(E1+E2)/2E");  
            dgm.makeH2("d113",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg))",     "SQRT(E1*E2)/E");  
            									 tit = "pc=2 ec=2";
            dgm.makeH2("d120",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1"); 
            dgm.makeH2("d130",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2"); 
            dgm.makeH2("d122",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1"); 
            dgm.makeH2("d132",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2");  
            dgm.makeH2("d123",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg))",     "E1/E #gamma1");  
            dgm.makeH2("d133",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg))",     "E2/E #gamma2");              

            break;
    	case 1:
    		dgm.add("GENREC",6,4,0,st,getRunNumber()); int y1=-2, y2=2;
                                                   tit = "pc>1 ec=0";
            dgm.makeH2("dt00",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1");
            dgm.makeH2("dt10",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2");
            dgm.makeH2("dt02",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1"); 
            dgm.makeH2("dt12",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2"); 
            dgm.makeH2("dt03",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma1");  
            dgm.makeH2("dt13",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma2");  
                                                   tit = "pc>1 ec=1";
            dgm.makeH2("dt20",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1"); 
            dgm.makeH2("dt30",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2");  
            dgm.makeH2("dt22",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1");   
            dgm.makeH2("dt32",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2");  
            dgm.makeH2("dt23",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma1");   
            dgm.makeH2("dt33",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma2");
                                                   tit = "pc>1 ec>1";
            dgm.makeH2("dt40",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1"); 
            dgm.makeH2("dt50",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2"); 
            dgm.makeH2("dt42",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1"); 
            dgm.makeH2("dt52",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2");  
            dgm.makeH2("dt43",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma1");  
            dgm.makeH2("dt53",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma2"); 
                                                    tit = "pc>1 ec>0";            
            dgm.makeH2("dt60",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta12 (deg)"); 
            dgm.makeH2("dt62",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta12 (deg)");    
            dgm.makeH2("dt63",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta12 (deg)");  
                                                   tit = "pc=1 ec>1"; 
            dgm.makeH2("dt70",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta12 (deg)");   
            dgm.makeH2("dt72",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta12 (deg)");  
            dgm.makeH2("dt73",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta12 (deg)");  
            break;
    	case 2:
    		dgm.add("GENREC",6,4,0,st,getRunNumber()); y1=-3; y2=3;
                                                   tit = "pc>1 ec=0";
            dgm.makeH2("df00",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1");
            dgm.makeH2("df10",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2");
            dgm.makeH2("df02",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1"); 
            dgm.makeH2("df12",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2"); 
            dgm.makeH2("df03",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma1");  
            dgm.makeH2("df13",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma2");  
                                                   tit = "pc>1 ec=1";
            dgm.makeH2("df20",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1"); 
            dgm.makeH2("df30",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2");  
            dgm.makeH2("df22",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1");   
            dgm.makeH2("df32",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2");  
            dgm.makeH2("df23",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma1");   
            dgm.makeH2("df33",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma2");
                                                   tit = "pc>1 ec>1";
            dgm.makeH2("df40",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1"); 
            dgm.makeH2("df50",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2"); 
            dgm.makeH2("df42",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1"); 
            dgm.makeH2("df52",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2");  
            dgm.makeH2("df43",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma1");  
            dgm.makeH2("df53",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma2");    		            
    	}        
   	   
    }   
    
    public void createEFFICIENCY(int st) {
  
    	switch (st) {        
        case 0:          
        	dgm.add("EFFICIENCY",6,5,0,st,getRunNumber());
            dgm.makeH1("eff1",50,0,5.2,-1,"n#gamma>1",          "Opening Angle (deg)","Efficiency",1,3);
            dgm.makeH1("eff1a",50,0,5.2,-2,"pc>0 ec>1",         "Opening Angle (deg)","Efficiency",1,4);
            dgm.makeH1("eff3",50,0,5.2,-2,"pc>1 ec=1 pc>1 ec>1","Opening Angle (deg)","Efficiency",1,5);
            dgm.makeH1("eff4",50,0,5.2,-2,"pc>1 ec=2",          "Opening Angle (deg)","Efficiency",1,2);
            dgm.makeH1("eff5",50,0,5.2,-2,"pc=2 ec=2",          "Opening Angle (deg)","Efficiency",1,1);
            dgm.makeH1("eff1a"); dgm.makeH1("eff3"); 
            dgm.makeH1("eff2",50,0,5.2,-2,"pc>1 ec=1",          "Opening Angle (deg)","Efficiency",1,0);                      
            dgm.makeH1("eff4");  dgm.makeH1("eff5"); dgm.makeH1("eff2");
            dgm.makeH2("ef10",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","PC Clusters" ,"Photons"); dgm.cc("ef10",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH2("ef11",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","ECi Clusters","Photons"); dgm.cc("ef11",false,false, 0,0,0.8f,1.2f);  
            dgm.makeH2("ef12",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","ECo Clusters","Photons"); dgm.cc("ef12",false,false, 0,0,0.8f,1.2f);              
            dgm.makeH2("ef13",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","PC Clusters" ,"Photons");  dgm.cc("ef13",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH2("ef14",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","ECi Clusters","Photons");  dgm.cc("ef14",false,false, 0,0,0.8f,1.2f);  
            dgm.makeH2("ef15",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","ECo Clusters","Photons");  dgm.cc("ef15",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH1("eff111",50,0,5.2,-1,"pc>1 ec>0",        "Opening Angle (deg)","Efficiency",1,2);
            dgm.makeH1("ef16",180,0,180,-1,"","#gamma #gamma #phi (deg)","");             
        	dgm.makeH1("ef18",160,-80,80,-1,"","dU","",1,4);        	       	
           break;
        case 1:
        	dgm.add("EFFICIENCY",4,4,0,st,getRunNumber());          
            dgm.makeH1("h10",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
            dgm.makeH1("h1a",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
        	dgm.makeH1("h11",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,5);
        	dgm.makeH1("h111",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,5);
        	dgm.makeH1("h12",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,2);
        	dgm.makeH1("h13",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,1);
        	dgm.makeH1("h14",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
        	dgm.makeH1("h15",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);        	
        	dgm.makeH1("h16",180,0,180,-1,"","#gamma #gamma #phi (deg)","",1,4);        	
        	dgm.makeH1("h17",180,0,180,-1,"","#gamma #gamma #phi (deg)","",1,4);        	
            dgm.makeH2("e10",16,-0.5,15.5, 16,-0.5,15.5,-1,"","PC Clusters" ,"Photons");  
            dgm.makeH2("e11",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECi Clusters","Photons");  
            dgm.makeH2("e12",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECo Clusters","Photons");             
            dgm.makeH2("e13",16,-0.5,15.5, 16,-0.5,15.5,-1,"","PC Clusters" ,"Photons");  
            dgm.makeH2("e14",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECi Clusters","Photons");  
            dgm.makeH2("e15",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECo Clusters","Photons");             
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
    
    public List<Float> getkin(List<Particle> list) {
    	List<Float> out = new ArrayList<Float>();
    	int n=0;
    	for (Particle p : list) {
		    out.add(n++,(float) p.e());
		    out.add(n++,(float) Math.toDegrees(p.theta())); 
		    out.add(n++,(float) Math.toDegrees(p.phi()));
    	}
		return out;
    }
    
    public double getSF(DetectorParticle dp) {
    	return SamplingFractions.getMean(22, dp, ebmce.ccdb);
    }
    
    @Override
    public void processEvent(DataEvent de) {
		
    	GEN.clear();  REC.clear();
    	
    	DetectorParticle p1 = new DetectorParticle();
    	DetectorParticle p2 = new DetectorParticle();
		
		int npart=0;
		double e1,e2;
		
		boolean debug = false, correct=true;
						
		boolean goodev = ebmce.readMC(de) && ebmce.pmc.size()==2;
				
        if (goodev) { 
        	
        	GEN = getkin(ebmce.pmc);   
        	
	    	double  opa = ebmce.pmv.get(0).theta(ebmce.pmv.get(1));
	    	Vector3 ggc = ebmce.pmv.get(0).cross(ebmce.pmv.get(1));
	    	double ggp = Math.toDegrees(Math.atan2(ggc.y(),ggc.x()));
	    	if(ggp<0) ggp=-ggp;
	    	ggp=ggp-90;
	    	if(ggp<0) ggp=ggp+180;	  
	    	
	    	if(!ebmce.hasStartTime) engine.processDataEvent(de);  
            
        	if (!ebmce.processDataEvent(de)) return;
        	
        	np.clear(); np = ebmce.eb.getEvent().getParticles();  
        	
        	List<DetectorResponse> cal = ebmce.eb.getEvent().getCalorimeterResponseList(); 
/*        	
        	System.out.println(" ");        	
            for(int iresp = 0; iresp < cal.size(); iresp++){
               CalorimeterResponse r = (CalorimeterResponse)cal.get(iresp);
                for(int iass = 0; iass < r.getNAssociations(); iass++) {
                    System.out.println(iresp+" "+r.getHitIndex()+" "+r.getAssociation(iass)+" "+
                    +r.getDescriptor().getSector()+" "+r.getDescriptor().getLayer());
                }            
            }
*/        	
        	if(dumpFiles) {ebmce.getRECBanks(de,ebmce.eb); writer.writeEvent(de);}
       
        	List<Particle> plist = new ArrayList<Particle>(); 
        	
//        	System.out.println("PART ");
        	
        	int trsec = -1; int trpid = -211;
        	for (DetectorParticle dp: np) {
        		DetectorResponse dr = dp.getHit(DetectorType.ECAL);
//        		System.out.println(dp.getPid()+" "+dp.getSector(DetectorType.ECAL));
        		int pid = dp.getPid(); int sec = dp.getSector(DetectorType.ECAL);
        		if(trsec==-1 && pid==trpid) trsec=sec;
        	}
        	
        	for (DetectorParticle dp : np) { // make list of neutral Particle objects 
		        if(dp.getSector(DetectorType.ECAL)!=trsec && dp.getPid()==22) { npart++;
			    if(!ebmce.hasStartTime && dp.getPid()==2112) {// this repairs zero momentum neutrons from non-PCAL seeded neutrals
	 				double e = dp.getEnergy(DetectorType.ECAL)/getSF(dp); 		
			    	Vector3D vec = dp.getHit(DetectorType.ECAL).getPosition(); vec.unit(); 			    		
 			    	dp.vector().add(new Vector3(e*vec.x(),e*vec.y(),e*vec.z())); //track energy for neutrals in DetectorParticle
 			    	dp.setPid(22);
 			    }
			    plist.add(dp.getPhysicsParticle(22)); //only this gives you SF corrected Particle energy from DetectorParticle
		        }
 			}
        	
       	    dgm.fill("c20",opa,npart); 
       	    dgm.fill("h10",opa);
       	    dgm.fill("h16",ggp);
       	    
     		if (npart>=2) {
     			
 			    	double dist=0, du=0;
 					int[]     npc = new int[50];       int[] neci = new int[50];       int[] neco = new int[50]; 
 					double[]  epc = new double[50]; double[]  eci = new double[50]; double[]  eco = new double[50]; 
 					 
					Vector3D[] r1 = new Vector3D[50]; Vector3[] c1 = new Vector3[50];
					Vector3D[] r4 = new Vector3D[50]; Vector3[] c4 = new Vector3[50]; 
					Vector3D[] r7 = new Vector3D[50]; Vector3[] c7 = new Vector3[50]; 
               	
 			    	double nopa = Math.toDegrees(Math.acos(plist.get(0).cosTheta(plist.get(1))));
                    double dopa = opa-nopa;  //GEN-REC               
 			    	
 			    	int npp = 0, ipp = 0;
 			    	for (DetectorParticle dp : np) {
 			    		if(dp.getSector(DetectorType.ECAL)==2 && dp.getPid()==22) {
 			    		if(npp==0) p1 = dp; //Photon 1
 			    		if(npp==1) p2 = dp; //Photon 2
 			            for(int iresp = 0; iresp < cal.size(); iresp++){
 			               CalorimeterResponse dr = (CalorimeterResponse)cal.get(iresp); 			              
 			    		   int lay = dr.getDescriptor().getLayer(); int index = dr.getAssociation(0); 			    		  
 			    		   if (index==ipp && dr.getDescriptor().getType()==DetectorType.ECAL) {
// 			            for (DetectorResponse dr : dp.getDetectorResponses()) { 			            	
// 			            	int lay = dr.getDescriptor().getType()==DetectorType.ECAL ? dr.getDescriptor().getLayer():0; 			           
 			            	if(lay==1) { npc[0]++ ; npc[npp+1]++ ; epc[npp+1]=dr.getEnergy() ; r1[npp+1]=dr.getPosition(); c1[npp+1]=dr.getCoordUVW();}    					
 			            	if(lay==4) {neci[0]++ ;neci[npp+1]++ ; eci[npp+1]=dr.getEnergy() ; r4[npp+1]=dr.getPosition(); c4[npp+1]=dr.getCoordUVW();}    					
 			            	if(lay==7) {neco[0]++ ;neco[npp+1]++ ; eco[npp+1]=dr.getEnergy() ; r7[npp+1]=dr.getPosition(); c7[npp+1]=dr.getCoordUVW();}
 			                }
 			            } 			    		   
 			            npp++;
 			    		}
 			    		ipp++;
			    	}
 			    	
 			    	e1=p1.getEnergy(DetectorType.ECAL)/getSF(p1); 
 			    	e2=p2.getEnergy(DetectorType.ECAL)/getSF(p2); 
 			    	
 			    	if (opa>1.95) {
 			    	float b1 = (float) p1.getBeta(); float b2 = (float) p2.getBeta();
 			    	if(npart>1)  {dgm.fill("p00",b1);       dgm.fill("p01",b2); 
 			    	              dgm.fill("p04",b1,npart); dgm.fill("p05",b2,npart);
 			    	              dgm.fill("p06",b1,e1);    dgm.fill("p07",b2,e2);}
 			    	if(npart==2) {dgm.fill("p02",b1);dgm.fill("p03",b2);}
 			    	}
		    	
 			    	if(npc[0]>=2)                            {r1[2].sub(r1[1]); dist=r1[2].mag(); c1[2].sub(c1[1]); du=c1[2].x();}
 			    	if(npc[1]==1 && npc[2]==0 && neci[2]>=1) {r4[2].sub(r1[1]); dist=r4[2].mag(); c4[2].sub(c1[1]); du=c4[2].x();}
 			    	if(npc[2]==1 && npc[1]==0 && neci[1]>=1) {r4[1].sub(r1[2]); dist=r4[1].mag(); c4[1].sub(c1[2]); du=c4[1].x();}
 			    	 			    	
 			    	if (correct && p1.countResponses(DetectorType.ECAL,1)==1 && p2.countResponses(DetectorType.ECAL,1)==1) {
 			    		//Merged cluster associated with photon 2 	
 	 			    	if (p1.countResponses(DetectorType.ECAL,4)==0 && p2.countResponses(DetectorType.ECAL,4)==1) {
 	 			    		double rat12 = 0.5*epc[1]/epc[2]; double rat21 = 0.5*epc[2]/epc[1]; rat12=0.5;rat21=0.5;
 	 			    		plist.get(0).setP((epc[1]+eci[2]*rat21+eco[1])/getSF(p2));
 	 			    		plist.get(1).setP((epc[2]+eci[2]*rat12+eco[2])/getSF(p2)); 	 			    		
 	 			    	}
 	 			        //Merged cluster associated with photon 1 	
 	 	 			    if (p1.countResponses(DetectorType.ECAL,4)==1 && p2.countResponses(DetectorType.ECAL,4)==0) {
	 			    		double rat12 = 0.5*epc[1]/epc[2]; double rat21 = 0.5*epc[2]/epc[1]; rat12=0.5; rat21=0.5;
 	 	 			    	plist.get(0).setP((epc[1]+eci[1]*rat21+eco[1])/getSF(p1));
 	 	 			    	plist.get(1).setP((epc[2]+eci[1]*rat12+eco[2])/getSF(p1));	 	 			    	
 	 	 			    }	 	 			    
 			    	} 	
 			    	
 	 			    REC=getkin(plist); 
 	 			    
     				dgm.fill("h11",opa);
     	       	    dgm.fill("h17",ggp);
     	       	    dgm.fill("ef18",du);
     			 
     	       	    //Truth matching based on fixed angle of gamma 1
      				Boolean swap = Math.abs(GEN.get(1)- REC.get(1))<0.17 ? false:true;

     				double  delE1 = REC.get(swap?3:0)/GEN.get(0);
     				double delTH1 = REC.get(swap?4:1)-GEN.get(1);
     				double delPH1 = REC.get(swap?5:2)-GEN.get(2);
     				double  delE2 = REC.get(swap?0:3)/GEN.get(3);
     				double delTH2 = REC.get(swap?1:4)-GEN.get(4);
     				double delPH2 = REC.get(swap?2:5)-GEN.get(5);
     				
 	 			    double delesum =          (REC.get(0)+REC.get(3))/(GEN.get(0)+GEN.get(3));    	
 	 			    double delepro = Math.sqrt(REC.get(0)*REC.get(3))/Math.sqrt(GEN.get(0)*GEN.get(3)); 	 			    
 	 			    
 	 			   //Debug 
     				if(debug && neci[0]==1 && opa>1.8 && opa<2.5) {
     				System.out.println(" ");
     				System.out.println(getEventNumber());
     				
     				double sf = getSF(p1);
     				
			        System.out.println(p1.getEnergy(DetectorType.ECAL)+" "+epc[1]+" "+eci[1]+" "+eco[1]);
			        System.out.println(p2.getEnergy(DetectorType.ECAL)+" "+epc[2]+" "+eci[2]+" "+eco[2]);
     				
     				System.out.println(p1.getEnergy(DetectorType.ECAL)+" "+(epc[1]+eci[1]/2+eco[1])/sf);
     				System.out.println(p2.getEnergy(DetectorType.ECAL)+" "+(epc[2]+eci[1]/2+eco[2])/sf);
			    	System.out.println(npart+" "+npc[0]+" "+neci[0]+" "+neco[0]);
 	 			    System.out.println("Photon 1: "+p1.getMass()+" "+e1+" "+npc[1]+" "+neci[1]); 
 	 			    System.out.println("Photon 2: "+p2.getMass()+" "+e2+" "+npc[2]+" "+neci[2]); 
 	 			    
     				System.out.println("GEN,REC EN1,EN2 "+GEN.get(0)+" "+ REC.get(0)+" "+GEN.get(3)+" "+ REC.get(3));
     				System.out.println("GEN,REC TH1,TH2 "+GEN.get(1)+" "+ REC.get(1)+" "+GEN.get(4)+" "+ REC.get(4));
     				System.out.println("GEN,REC PH1,PH2 "+GEN.get(2)+" "+ REC.get(2)+" "+GEN.get(5)+" "+ REC.get(5));
 			    	System.out.println("GEN,REC opa "+opa+" "+nopa + " "+dist);
 		
    				System.out.println("Phot 1 "+swap+" "+delE1+" "+delTH1+" "+delPH1);
     				System.out.println("Phot 2 "+swap+" "+delE2+" "+delTH2+" "+delPH2);   
     				}
     				
     				if (opa>1.95) {
     				dgm.fill("c00", npc[0],npart);        dgm.fill("c01", neci[0],npart);        dgm.fill("c02", neco[0],npart);
     				dgm.fill("e10", npc[0],npart,delesum);dgm.fill("e11", neci[0],npart,delesum);dgm.fill("e12", neco[0],npart,delesum);
     				}
     				if (opa>2.40) {
     				dgm.fill("c05", npc[0],npart);        dgm.fill("c06", neci[0],npart);        dgm.fill("c07", neco[0],npart);
     				dgm.fill("e13", npc[0],npart,delesum);dgm.fill("e14", neci[0],npart,delesum);dgm.fill("e15", neco[0],npart,delesum);
     				}
     				
     			    if (npc[0]>=2 && neci[0]==0) {
     			    	dgm.fill("c21",opa,npart);
     			    	dgm.fill("d00",GEN.get(4),delE1);
     			    	dgm.fill("d02",dist,      delE1);        			
     			    	dgm.fill("d03",opa,       delE1);        			
     			    	dgm.fill("d10",GEN.get(4),delE2);
     			    	dgm.fill("d12",dist,      delE2);        			
     			    	dgm.fill("d13",opa,       delE2); 
     			    	
     			    	dgm.fill("dt00",GEN.get(4),delTH1);
     			    	dgm.fill("dt02",dist,      delTH1);        			
     			    	dgm.fill("dt03",opa,       delTH1);        			
     			    	dgm.fill("dt10",GEN.get(4),delTH2);
     			    	dgm.fill("dt12",dist,      delTH2);        			
     			    	dgm.fill("dt13",opa,       delTH2); 
     			    	
     			    	dgm.fill("df00",GEN.get(4),delPH1);
     			    	dgm.fill("df02",dist,      delPH1);        			
     			    	dgm.fill("df03",opa,       delPH1);        			
     			    	dgm.fill("df10",GEN.get(4),delPH2);
     			    	dgm.fill("df12",dist,      delPH2);        			
     			    	dgm.fill("df13",opa,       delPH2); 
     			    }
     				
     				if (npc[0]>=2 && neci[0]==1) {
     			    	dgm.fill("c22",opa,npart);
     					dgm.fill("d20",GEN.get(4),delE1);
     					dgm.fill("d22",dist,      delE1);        			
     					dgm.fill("d23",opa,       delE1);        			
     					dgm.fill("d30",GEN.get(4),delE2);
     					dgm.fill("d32",dist,      delE2);        			     				
     					dgm.fill("d33",opa,       delE2);  
     					
     			    	dgm.fill("dt20",GEN.get(4),delTH1);
     			    	dgm.fill("dt22",dist,      delTH1);        			
     			    	dgm.fill("dt23",opa,       delTH1);        			
     			    	dgm.fill("dt30",GEN.get(4),delTH2);
     			    	dgm.fill("dt32",dist,      delTH2);        			
     			    	dgm.fill("dt33",opa,       delTH2); 
     			    	
     			    	dgm.fill("df20",GEN.get(4),delPH1);
     			    	dgm.fill("df22",dist,      delPH1);        			
     			    	dgm.fill("df23",opa,       delPH1);        			
     			    	dgm.fill("df30",GEN.get(4),delPH2);
     			    	dgm.fill("df32",dist,      delPH2);        			
     			    	dgm.fill("df33",opa,       delPH2); 
     			    	
     					dgm.fill("d60",GEN.get(4),delesum);
     					dgm.fill("d62",dist,      delesum);        			
     					dgm.fill("d63",opa,       delesum);        			
     					dgm.fill("d70",GEN.get(4),delepro);
     					dgm.fill("d72",dist,      delepro);        			     				
     					dgm.fill("d73",opa,       delepro);   
     					dgm.fill("h12",opa);
     			    }
     				
     				              dgm.fill("c03",npart,delesum);
     				if (opa>1.95) dgm.fill("c04",npart,delesum);
     				
     				if (npc[0]>=2 && neci[0]>=2) {
     			    	dgm.fill("c23",opa,npart);
     					dgm.fill("d40",GEN.get(4),delE1);
     					dgm.fill("d42",dist,      delE1);        			
     					dgm.fill("d43",opa,       delE1);        			
     					dgm.fill("d50",GEN.get(4),delE2);
     					dgm.fill("d52",dist,      delE2);        			     				
     					dgm.fill("d53",opa,       delE2);
     					
     					dgm.fill("d80",GEN.get(4),delesum);
     					dgm.fill("d82",dist,      delesum);        			
     					dgm.fill("d83",opa,       delesum);        			
     					dgm.fill("d90",GEN.get(4),delepro);
     					dgm.fill("d92",dist,      delepro);        			     				
     					dgm.fill("d93",opa,       delepro);   
     					
     			    	dgm.fill("dt40",GEN.get(4),delTH1);
     			    	dgm.fill("dt42",dist,      delTH1);        			
     			    	dgm.fill("dt43",opa,       delTH1);        			
     			    	dgm.fill("dt50",GEN.get(4),delTH2);
     			    	dgm.fill("dt52",dist,      delTH2);        			
     			    	dgm.fill("dt53",opa,       delTH2); 
     			    	
     			    	dgm.fill("df40",GEN.get(4),delPH1);
     			    	dgm.fill("df42",dist,      delPH1);        			
     			    	dgm.fill("df43",opa,       delPH1);        			
     			    	dgm.fill("df50",GEN.get(4),delPH2);
     			    	dgm.fill("df52",dist,      delPH2);        			
     			    	dgm.fill("df53",opa,       delPH2); 
     					dgm.fill("h13", opa);
     				}
     				
     				if (npc[0]>=2 && neci[0]>0) {
     			    	dgm.fill("dt60",GEN.get(4),dopa);
     					dgm.fill("dt62",dist,      dopa);        			
     					dgm.fill("dt63",opa,       dopa); 
     					dgm.fill("h111",opa);
     				}
     				
     				if (npc[0]>0  && neci[0]>=2) dgm.fill("h1a", opa);
     				if (npc[0]>=2 && neci[0]==2) dgm.fill("h14", opa);
     				if (npc[0]==2 && neci[0]==2) {
     					dgm.fill("h15", opa);
     					dgm.fill("d100",GEN.get(4),delesum);
     					dgm.fill("d102",dist,      delesum);        			
     					dgm.fill("d103",opa,       delesum);        			
     					dgm.fill("d110",GEN.get(4),delepro);
     					dgm.fill("d112",dist,      delepro);        			     				
     					dgm.fill("d113",opa,       delepro);   
     					dgm.fill("d120",GEN.get(4),delE1);
     					dgm.fill("d122",dist,      delE1);        			
     					dgm.fill("d123",opa,       delE1);        			
     					dgm.fill("d130",GEN.get(4),delE2);
     					dgm.fill("d132",dist,      delE2);        			     				
     					dgm.fill("d133",opa,       delE2);     					
     				}
   				
     				if (npc[0]==1 && neci[0]>=2) {
     			    	dgm.fill("c24",opa,npart);
     					dgm.fill("dt70",GEN.get(4),dopa);
     					dgm.fill("dt72",dist,      dopa);        			
     					dgm.fill("dt73",opa,       dopa);  
     				}
     			}    			     			
     		     		
       	    }   
    }
  
   
    public void plotMCHistos() {  
    	plot("PHOTONS");
        plot("CLUSTERS");
        plot("GENREC");
        plot("EFFICIENCY");   
    }
    
    @Override
    public void plotEvent(DataEvent de) {
       // analyze();         
    }
    
    public void geteff() {
		dgm.geteff("eff1",  "h11", "h10");
		dgm.geteff("eff1a", "h1a", "h10");
		dgm.geteff("eff2",  "h12", "h10");
		dgm.geteff("eff3",  "h13", "h10");
		dgm.geteff("eff4",  "h14", "h10");
		dgm.geteff("eff5",  "h15", "h10");
		dgm.geteff("eff111","h111","h10");
		dgm.geteff("ef10",  "e10", "c00");
		dgm.geteff("ef11",  "e11", "c01");
		dgm.geteff("ef12",  "e12", "c02");
		dgm.geteff("ef13",  "e13", "c05");
		dgm.geteff("ef14",  "e14", "c06");
		dgm.geteff("ef15",  "e15", "c07");
		dgm.geteff("ef16",  "h17", "h16");
		
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    } 
  
    @Override
    public void timerUpdate() {  

    }


}
