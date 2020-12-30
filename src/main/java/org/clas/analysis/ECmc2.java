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
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc2 extends DetectorMonitor {
	
	Event ev = new Event();
	EBMC  eb = new EBMC();
    List<DetectorParticle> np = new ArrayList<DetectorParticle>();    
	List<Float> GEN =  new ArrayList<Float>();
	List<Float> REC1 = new ArrayList<Float>();   
    HipoDataSync  writer = null;		
	List<Particle> phot = new ArrayList<Particle>();
	List<Particle> neut = new ArrayList<Particle>();
	String tit = null;
	double ethresh = 0.3;
	
    public ECmc2(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("CLUSTERS","GENREC","EFFICIENCY");

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
        
        eb.getCCDB(10);
        eb.setThresholds("Pizero",engine);
        eb.setGeom("2.5");
        eb.setGoodPhotons(12);
        eb.isMC = true;
        
        tl.setFitData(Fits);
        writer = new HipoDataSync();
        writer.open("/Users/colesmith/CLAS12ANA/ECmc2/photon_dist.hipo");
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
    	createCLUSTERS(0);
    	createCLUSTERS(1);
    	createGENREC(0);
    	createGENREC(1);
    	createGENREC(2);
       	createEFFICIENCY(0);
       	createEFFICIENCY(1);
    }
    
    public void createCLUSTERS(int st) {
    	      
    	switch (st) {        
        case 0:  
        	dgm.add("CLUSTERS", 5,3,0,st,getRunNumber());
            dgm.makeH2("c00",16,-0.5,15.5,16,-0.5,15.5,-1,"opa>1.95","PC Clusters" ,"Photons"); 
            dgm.makeH2("c01",16,-0.5,15.5,16,-0.5,15.5,-1,"opa>1.95","ECi Clusters","Photons"); 
            dgm.makeH2("c02",16,-0.5,15.5,16,-0.5,15.5,-1,"opa>1.95","ECo Clusters","Photons"); 
            dgm.makeH2("c03",16,-0.5,15.5,80,0.4,1.6,1,"",          "Photons",     "REC (E1+E2)/2E");
            dgm.makeH2("c04",16,-0.5,15.5,80,0.4,1.6,1,"OPA>1.95",  "Photons",     "REC (E1+E2)/2E");
            dgm.makeH2("c20",50,0,5.2, 16,-0.5,15.5,-1,"",         "GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c21",50,0,5.2, 16,-0.5,15.5,-1,"pc=2 ec=0","GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c22",50,0,5.2, 16,-0.5,15.5,-1,"pc=2 ec=1","GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c23",50,0,5.2, 16,-0.5,15.5,-1,"pc=2 ec>1","GEN #theta12 (deg))","Photons"); 
            dgm.makeH2("c24",50,0,5.2, 16,-0.5,15.5,-1,"pc=1 ec>1","GEN #theta12 (deg))","Photons");  
            break;
        case 1:
        	dgm.add("CLUSTERS", 5,3,0,st,getRunNumber());
            dgm.makeH2("c05",16,-0.5,15.5,16,-0.5,15.5,-1,"",       "PC Clusters" ,"Photons"); 
            dgm.makeH2("c06",16,-0.5,15.5,16,-0.5,15.5,-1,"",       "ECi Clusters","Photons"); 
            dgm.makeH2("c07",16,-0.5,15.5,16,-0.5,15.5,-1,"",       "ECo Clusters","Photons");             
    	}       
    }
    
    public void createGENREC(int st) {
	  
    	switch (st) {        
    	case 0:
    		dgm.add("GENREC",6,6,0,st,getRunNumber());
        	                                      tit = "pc>1 ec=0";
            dgm.makeH2("d00",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma1");
            dgm.makeH2("d10",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma2");
            dgm.makeH2("d02",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma1"); 
            dgm.makeH2("d12",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma2"); 
            dgm.makeH2("d03",50,0,5.2,  50,-1,1,0,tit,"GEN #theta12 (deg))",     "#DeltaE/E #gamma1");  
            dgm.makeH2("d13",50,0,5.2,  50,-1,1,0,tit,"GEN #theta12 (deg))",     "#DeltaE/E #gamma2");  
                                                  tit = "pc>1 ec=1";
            dgm.makeH2("d20",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma1"); 
            dgm.makeH2("d30",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma2");  
            dgm.makeH2("d22",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma1");   
            dgm.makeH2("d32",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma2");  
            dgm.makeH2("d23",50,0,5.2,  50,-1,1,0,tit,"GEN #theta12 (deg))",     "#DeltaE/E #gamma1");   
            dgm.makeH2("d33",50,0,5.2,  50,-1,1,0,tit,"GEN #theta12 (deg))",     "#DeltaE/E #gamma2");
                                                  tit = "pc>1 ec>1";
            dgm.makeH2("d40",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma1"); 
            dgm.makeH2("d50",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma2"); 
            dgm.makeH2("d42",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma1"); 
            dgm.makeH2("d52",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma2");  
            dgm.makeH2("d43",50,0,5.2,  50,-1,1,0,tit,"GEN #theta12 (deg))",     "#DeltaE/E #gamma1");  
            dgm.makeH2("d53",50,0,5.2,  50,-1,1,0,tit,"GEN #theta12 (deg))",     "#DeltaE/E #gamma2");              
                                                  tit = "pc=2 ec=1"; 
            dgm.makeH2("d60",50,10,20,  50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d70",50,10,20,  50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d62",50,0,70,   50, 0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d72",50,0,70,   50, 0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d63",50,0,5.2,  50, 0,2,1,tit,"GEN #theta12 (deg))",     "(E1+E2)/2E");  
            dgm.makeH2("d73",50,0,5.2,  50, 0,2,1,tit,"GEN #theta12 (deg))",     "SQRT(E1*E2)/E");  
                                                  tit = "pc>1 ec>1";
            dgm.makeH2("d80",50,10,20,  50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d90",50,10,20,  50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d82",50,0,70,   50, 0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d92",50,0,70,   50, 0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d83",50,0,5.2,  50, 0,2,1,tit,"GEN #theta12 (deg))",     "(E1+E2)/2E");  
            dgm.makeH2("d93",50,0,5.2,  50, 0,2,1,tit,"GEN #theta12 (deg))",     "SQRT(E1*E2)/E");  
                                                  tit = "pc=2 ec=2";
            dgm.makeH2("d100",50,10,20, 50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d110",50,10,20, 50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d102",50,0,70,  50, 0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d112",50,0,70,  50, 0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d103",50,0,5.2, 50, 0,2,1,tit,"GEN #theta12 (deg))",     "(E1+E2)/2E");  
            dgm.makeH2("d113",50,0,5.2, 50, 0,2,1,tit,"GEN #theta12 (deg))",     "SQRT(E1*E2)/E");  

            break;
    	case 1:
    		dgm.add("GENREC",6,4,0,st,getRunNumber()); int y1=-2, y2=2;
                                                   tit = "pc>1 ec=0";
            dgm.makeH2("dt00",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1");
            dgm.makeH2("dt10",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2");
            dgm.makeH2("dt02",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1"); 
            dgm.makeH2("dt12",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2"); 
            dgm.makeH2("dt03",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma1");  
            dgm.makeH2("dt13",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma2");  
                                                   tit = "pc>1 ec=1";
            dgm.makeH2("dt20",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1"); 
            dgm.makeH2("dt30",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2");  
            dgm.makeH2("dt22",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1");   
            dgm.makeH2("dt32",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2");  
            dgm.makeH2("dt23",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma1");   
            dgm.makeH2("dt33",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma2");
                                                   tit = "pc>1 ec>1";
            dgm.makeH2("dt40",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1"); 
            dgm.makeH2("dt50",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2"); 
            dgm.makeH2("dt42",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1"); 
            dgm.makeH2("dt52",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2");  
            dgm.makeH2("dt43",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma1");  
            dgm.makeH2("dt53",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta #gamma2");    		                                                  
            dgm.makeH2("dt60",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta12 (deg)"); 
            dgm.makeH2("dt62",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta12 (deg)");    
            dgm.makeH2("dt63",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta12 (deg)");  
                                                   tit = "pc=1 ec>1"; 
            dgm.makeH2("dt70",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta12 (deg)");   
            dgm.makeH2("dt72",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta12 (deg)");  
            dgm.makeH2("dt73",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#theta12 (deg)");  
            break;
    	case 2:
    		dgm.add("GENREC",6,4,0,st,getRunNumber()); y1=-3; y2=3;
                                                   tit = "pc>1 ec=0";
            dgm.makeH2("df00",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1");
            dgm.makeH2("df10",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2");
            dgm.makeH2("df02",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1"); 
            dgm.makeH2("df12",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2"); 
            dgm.makeH2("df03",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma1");  
            dgm.makeH2("df13",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma2");  
                                                   tit = "pc>1 ec=1";
            dgm.makeH2("df20",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1"); 
            dgm.makeH2("df30",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2");  
            dgm.makeH2("df22",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1");   
            dgm.makeH2("df32",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2");  
            dgm.makeH2("df23",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma1");   
            dgm.makeH2("df33",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma2");
                                                   tit = "pc>1 ec>1";
            dgm.makeH2("df40",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1"); 
            dgm.makeH2("df50",50,10,20,  50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2"); 
            dgm.makeH2("df42",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1"); 
            dgm.makeH2("df52",50,0,70,   50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2");  
            dgm.makeH2("df43",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma1");  
            dgm.makeH2("df53",50,0,5.2,  50,y1,y2,0,tit,"GEN #theta12 (deg))",     "#Delta#phi #gamma2");    		            
    	}        
   	   
    }   
    
    public void createEFFICIENCY(int st) {
  
    	switch (st) {        
        case 0:          
        	dgm.add("EFFICIENCY",6,5,0,st,getRunNumber());
            dgm.makeH1("eff1",50,0,5.2,-1,"n#gamma>1","Opening Angle (deg)","Efficiency",1,3);
            dgm.makeH1("eff1a",50,0,5.2,-2,"pc>0 ec>1","Opening Angle (deg)","Efficiency",1,4);
            dgm.makeH1("eff3",50,0,5.2,-2,"pc>1 ec>1","Opening Angle (deg)","Efficiency",1,5);
            dgm.makeH1("eff4",50,0,5.2,-2,"pc>1 ec=2","Opening Angle (deg)","Efficiency",1,2);
            dgm.makeH1("eff5",50,0,5.2,-2,"pc=2 ec=2","Opening Angle (deg)","Efficiency",1,1);
            dgm.makeH1("eff1a"); dgm.makeH1("eff3"); 
            dgm.makeH1("eff2",50,0,5.2,-2,"pc>1 ec=1","Opening Angle (deg)","Efficiency",1,0);                      
            dgm.makeH1("eff4");  dgm.makeH1("eff5"); dgm.makeH1("eff2");
            dgm.makeH2("ef10",16,-0.5,15.5,16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","PC Clusters" ,"Photons"); dgm.cc("ef10",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH2("ef11",16,-0.5,15.5,16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","ECi Clusters","Photons"); dgm.cc("ef11",false,false, 0,0,0.8f,1.2f);  
            dgm.makeH2("ef12",16,-0.5,15.5,16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","ECo Clusters","Photons"); dgm.cc("ef12",false,false, 0,0,0.8f,1.2f);              
            dgm.makeH2("ef13",16,-0.5,15.5,16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","PC Clusters" ,"Photons");  dgm.cc("ef13",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH2("ef14",16,-0.5,15.5,16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","ECi Clusters","Photons");  dgm.cc("ef14",false,false, 0,0,0.8f,1.2f);  
            dgm.makeH2("ef15",16,-0.5,15.5,16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","ECo Clusters","Photons");  dgm.cc("ef15",false,false, 0,0,0.8f,1.2f);              
           break;
        case 1:
        	dgm.add("EFFICIENCY",4,3,0,st,getRunNumber());          
            dgm.makeH1("h10",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
            dgm.makeH1("h1a",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
        	dgm.makeH1("h11",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,5);
        	dgm.makeH1("h12",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,2);
        	dgm.makeH1("h13",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,1);
        	dgm.makeH1("h14",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
        	dgm.makeH1("h15",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);        	
            dgm.makeH2("e10",16,-0.5,15.5,16,-0.5,15.5,-1,"","PC Clusters" ,"Photons");  
            dgm.makeH2("e11",16,-0.5,15.5,16,-0.5,15.5,-1,"","ECi Clusters","Photons");  
            dgm.makeH2("e12",16,-0.5,15.5,16,-0.5,15.5,-1,"","ECo Clusters","Photons");             
            dgm.makeH2("e13",16,-0.5,15.5,16,-0.5,15.5,-1,"","PC Clusters" ,"Photons");  
            dgm.makeH2("e14",16,-0.5,15.5,16,-0.5,15.5,-1,"","ECi Clusters","Photons");  
            dgm.makeH2("e15",16,-0.5,15.5,16,-0.5,15.5,-1,"","ECo Clusters","Photons");             
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
    
    public List<Float> getkin (List<Particle> list) {
    	List<Float> out = new ArrayList<Float>();
    	int n=0;
    	for (Particle p : list) {
		    out.add(n++,(float) p.e()); 
		    out.add(n++,(float) Math.toDegrees(p.theta())); 
		    out.add(n++,(float) Math.toDegrees(p.phi()));
    	}
		return out;
    }
    
    @Override
    public void processEvent(DataEvent de) {
		
    	DetectorParticle p1 = new DetectorParticle();
    	DetectorParticle p2 = new DetectorParticle();
		GEN.clear(); REC1.clear();
		
		int indx=-1; float pthresh=1f;
		Boolean debug = false;
		
		int n,np2,npp,sec,nphot=0, nneut=0, npart=0, etot=8;
		double e1,e2,oparec=-10,the2,epc1,epc2,eeci1, eeci2,eeco1,eeco2;
		Vector3D r11=null,r12=null, r41=null, r42=null;
				
		boolean goodev = eb.readMC(de) && eb.pmc.size()==2;
				
        if (goodev) { 
        	
        	GEN = getkin(eb.pmc);   
        	
	    	double ct1 = Math.cos(Math.toRadians(GEN.get(1))); double st1=Math.sqrt(1-ct1*ct1);
	    	double ct2 = Math.cos(Math.toRadians(GEN.get(4))); double st2=Math.sqrt(1-ct2*ct2);
	    	double opa = Math.toDegrees(Math.acos(ct1*ct2+st1*st2*Math.cos(Math.toRadians(GEN.get(2)-GEN.get(5))))); 
	    	
            engine.processDataEvent(de);  
            
        	eb.readEC(de,"ECAL::clusters");
        	np.clear(); np = eb.getNeutralPart(); npart = np.size(); 
        	eb.getRECBanks(de,eb.eb); writer.writeEvent(de);
       
        	List<Particle> plist = new ArrayList<Particle>(); 
        	
        	for (DetectorParticle dp : np) { 				
 				double e = dp.getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, dp, eb.ccdb); 				
			    if(dp.getPid()==2112) {// this repairs zero momentum neutrons from non-PCAL seeded neutrals
 			    	Vector3D vec = dp.getHit(DetectorType.ECAL).getPosition(); vec.unit(); 			    		
 			    	dp.vector().add(new Vector3(e*vec.x(),e*vec.y(),e*vec.z())); 
 			    	dp.setPid(22);
 			    }
			    plist.add(dp.getPhysicsParticle(22));			    
 			}
        	
       	    dgm.fill("c20",opa,npart); 
       	    dgm.fill("h10",opa);
     		
     		if (npart>=2) {
     			
 			    	double dist=0;
 					int[]     npc = new int[50];       int[] neci = new int[50];       int[] neco = new int[50]; 
 					double[]  epc = new double[50]; double[]  eci = new double[50]; double[]  eco = new double[50]; 
 					Vector3D[] r1 = new Vector3D[50]; 
 					Vector3D[] r4 = new Vector3D[50];
					Vector3D[] r7 = new Vector3D[50];
               	
 			    	e1=np.get(0).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p1, eb.ccdb);
 			    	e2=np.get(1).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p2, eb.ccdb);
 			    	
 			    	double nopa = Math.toDegrees(Math.acos(plist.get(0).cosTheta(plist.get(1))));
                    double dopa = opa-nopa;                  
 			    	
 			    	npp = 0;
 			    	for (DetectorParticle dp : np) {
 			    		if(npp==0) p1 = dp; //Photon 1
 			    		if(npp==1) p2 = dp; //Photon 2			    		
 			            for (DetectorResponse dr : dp.getDetectorResponses()) {
 			            	int lay = dr.getDescriptor().getLayer();
 			            	if(lay==1) { npc[0]++ ; npc[npp+1]++ ; epc[npp+1]=dr.getEnergy() ; r1[npp+1]=dr.getPosition();}    					
 			            	if(lay==4) {neci[0]++ ;neci[npp+1]++ ; eci[npp+1]=dr.getEnergy() ; r4[npp+1]=dr.getPosition();}    					
 			            	if(lay==7) {neco[0]++ ;neco[npp+1]++ ; eco[npp+1]=dr.getEnergy() ; r7[npp+1]=dr.getPosition();}
 			            }
 			            npp++;
			    	}
		    	
 			    	if(npc[0]>=2)                            {r1[2].sub(r1[1]); dist=r1[2].mag();}
 			    	if(npc[1]==1 && npc[2]==0 && neci[2]>=1) {r4[2].sub(r1[1]); dist=r4[2].mag();}
 			    	if(npc[2]==1 && npc[1]==0 && neci[1]>=1) {r4[1].sub(r1[2]); dist=r4[1].mag();}
 			    	
 	 			    REC1=getkin(plist); 
     				dgm.fill("h11",opa);
     			
     				double dth11 = GEN.get(1)-REC1.get(1);  
     				double dth14 = Math.abs(GEN.get(1)-REC1.get(4));
     				Boolean swap = Math.abs(dth11)<0.17 ? false:true;
     				
//     				System.out.println(dth11+" "+dth14+" "+swap);
     				
     				double  delE1 = swap ? GEN.get(0)-REC1.get(3):GEN.get(0)-REC1.get(0);
     				double delTH1 = swap ? GEN.get(1)-REC1.get(4):GEN.get(1)-REC1.get(1);
     				double delPH1 = swap ? GEN.get(2)-REC1.get(5):GEN.get(2)-REC1.get(2);
     				double  delE2 = swap ? GEN.get(3)-REC1.get(0):GEN.get(3)-REC1.get(3);
     				double delTH2 = swap ? GEN.get(4)-REC1.get(1):GEN.get(4)-REC1.get(4);
     				double delPH2 = swap ? GEN.get(5)-REC1.get(2):GEN.get(5)-REC1.get(5);
 	 			    double delesum =      0.5*(REC1.get(0)+REC1.get(3))/GEN.get(0);    	
 	 			    double delepro = Math.sqrt(REC1.get(0)*REC1.get(3))/GEN.get(0);
 	 			    
     				if(debug && dopa>0.5) {
     				System.out.println(" ");
     				System.out.println(getEventNumber());
     				
			    	System.out.println(npart+" "+nphot+" "+npc[0]+" "+neci[0]);
 	 			    System.out.println("Photon 1: "+p1.getMass()+" "+e1+" "+npc[1]+" "+neci[1]); 
 	 			    System.out.println("Photon 2: "+p2.getMass()+" "+e2+" "+npc[2]+" "+neci[2]); 
 	 			    
     				System.out.println("GEN,REC EN1,EN2 "+GEN.get(0)+" "+REC1.get(0)+" "+GEN.get(3)+" "+REC1.get(3));
     				System.out.println("GEN,REC TH1,TH2 "+GEN.get(1)+" "+REC1.get(1)+" "+GEN.get(4)+" "+REC1.get(4));
     				System.out.println("GEN,REC PH1,PH2 "+GEN.get(2)+" "+REC1.get(2)+" "+GEN.get(5)+" "+REC1.get(5));
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
     			    	dgm.fill("d00",GEN.get(4),delE1/GEN.get(0));
     			    	dgm.fill("d02",dist,      delE1/GEN.get(0));        			
     			    	dgm.fill("d03",opa,       delE1/GEN.get(0));        			
     			    	dgm.fill("d10",GEN.get(4),delE2/GEN.get(3));
     			    	dgm.fill("d12",dist,      delE2/GEN.get(3));        			
     			    	dgm.fill("d13",opa,       delE2/GEN.get(3)); 
     			    	
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
     					dgm.fill("d20",GEN.get(4),delE1/GEN.get(0));
     					dgm.fill("d22",dist,      delE1/GEN.get(0));        			
     					dgm.fill("d23",opa,       delE1/GEN.get(0));        			
     					dgm.fill("d30",GEN.get(4),delE2/GEN.get(3));
     					dgm.fill("d32",dist,      delE2/GEN.get(3));        			     				
     					dgm.fill("d33",opa,       delE2/GEN.get(3));  
     					
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
     					dgm.fill("h12", opa);
     			    }
     				
     				              dgm.fill("c03",npart,delesum);
     				if (opa>1.95) dgm.fill("c04",npart,delesum);
     				
     				if (npc[0]>=2 && neci[0]>=2) {
     			    	dgm.fill("c23",opa,npart);
     					dgm.fill("d40",GEN.get(4),delE1/GEN.get(0));
     					dgm.fill("d42",dist,      delE1/GEN.get(0));        			
     					dgm.fill("d43",opa,       delE1/GEN.get(0));        			
     					dgm.fill("d50",GEN.get(4),delE2/GEN.get(3));
     					dgm.fill("d52",dist,      delE2/GEN.get(3));        			     				
     					dgm.fill("d53",opa,       delE2/GEN.get(3));
     					
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

     			    	dgm.fill("dt60",GEN.get(4),dopa);
     					dgm.fill("dt62",dist,      dopa);        			
     					dgm.fill("dt63",opa,       dopa);        			     					
     					dgm.fill("h13", opa);
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
 
/*    			        	
     			npp=0;indx=-1;
        	
     			dg11.getH1F("h11").fill(GEN.get(6));
     			for (DetectorParticle phot : np) {
     				indx++;
     				double ep = phot.getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, phot, eb.ccdb);
     				if(ep>pthresh) {
//        				System.out.println(npp+" "+indx+" "+ep);
     					npp++;           				
     					for (DetectorResponse dr : phot.getDetectorResponses()) {
        					
     					}
     				}
     			}

     			dg0.getH2F("hnpp").fill(GEN.get(6),npp);
        	
     			if (npp>1)  dg11.getH1F("h12").fill(GEN.get(6));
     			if (npp==2) dg11.getH1F("h13").fill(GEN.get(6));
     			if (npp>=2) dg11.getH1F("h14").fill(GEN.get(6));
     			if (npart==2) {
     				double e1=0,e2=0,oparec=0,the2=0,dist=0;
            	
     				DetectorParticle p1 = np.get(0);  //Photon 1
     				DetectorParticle p2 = np.get(1);  //Photon 2  
               
     				Vector3 n1 = p1.vector(); n1.unit();
     				Vector3 n2 = p2.vector(); n2.unit();
                
     				e1=np.get(0).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p1, eb.ccdb);
     				e2=np.get(1).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p2, eb.ccdb);
                
     				List<Particle> list = new ArrayList<Particle>();
     				list.add(new Particle(22,n1.x()*e1,n1.y()*e1,n1.z()*e1)); 
     				list.add(new Particle(22,n2.x()*e2,n2.y()*e2,n2.z()*e2)); 
                             	
 //    				System.out.println(GEN.get(0)+" "+REC1.get(0)+" "+GEN.get(3)+" "+REC1.get(3));
 //    				System.out.println(GEN.get(1)+" "+REC1.get(1)+" "+GEN.get(4)+" "+REC1.get(4));
 //    				System.out.println(GEN.get(2)+" "+REC1.get(2)+" "+GEN.get(5)+" "+REC1.get(5));
                                
     				Vector3D r1=null,r2=null;
     				npc=0; nec=0;
     				double epc1=0,epc2=0,eeci1=0,eeci2=0;
     				for (DetectorResponse dr : p1.getDetectorResponses()) {
     					int lay = dr.getDescriptor().getLayer();
     					if(lay==1) {r1=dr.getPosition();npc++;epc1=dr.getEnergy();}    					
     					if(lay==4) {nec++;eeci1=dr.getEnergy();}    					
     				}
     				for (DetectorResponse dr : p2.getDetectorResponses()) {
     					int lay = dr.getDescriptor().getLayer();
     					if(lay==1) {r2=dr.getPosition();npc++;epc2=dr.getEnergy();}    					
     					if(lay==4) {nec++;eeci2=dr.getEnergy();}    					
     				}
     				dg0.getH2F("h2f").fill(GEN.get(6), (e1+e2)/etot-1);
     				if(npc==2) {
     					r2.sub(r1); dist=r2.mag();
     					dg0.getH2F("h2a").fill(dist, e2>0?e1/e2:1000);
     					dg0.getH2F("h2b").fill(dist, (e1+e2)/etot-1);
                		dg0.getH2F("h2e").fill(dist,GEN.get(6)-REC1.get(6));
                		if(nec==2)              dg0.getH2F("h2c" ).fill(epc1/epc2,eeci1/eeci2);
                		if(nec==2 && e1/e2>1.3) dg0.getH2F("h2cc").fill(epc1/epc2,eeci1/eeci2);
                	
                		dg0.getH2F("h2g").fill(npp,(e1+e2)/etot-1);
            			dgr.getH2F("d00").fill(GEN.get(3),1-delE1/GEN.get(0));
            			dgr.getH2F("d01").fill(GEN.get(2),1-delE1/GEN.get(0));
            			dgr.getH2F("d02").fill(dist,      1-delE1/GEN.get(0));
            			dgr.getH2F("d03").fill(dopa,      1-delE1/GEN.get(0));
            			dgr.getH2F("d10").fill(GEN.get(3),1-delE2/GEN.get(0));
            			dgr.getH2F("d11").fill(GEN.get(2),1-delE2/GEN.get(0));
            			dgr.getH2F("d12").fill(dist,      1-delE2/GEN.get(0));
            			dgr.getH2F("d13").fill(dopa,      1-delE2/GEN.get(0));
            			
     				}
     				
     			}
     			*/
        
    
   
    public void plotMCHistos() {      
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
		dgm.geteff("eff1a", "h1a", "h11");
		dgm.geteff("eff2",  "h12", "h11");
		dgm.geteff("eff3",  "h13", "h11");
		dgm.geteff("eff4",  "h14", "h11");
		dgm.geteff("eff5",  "h15", "h11");
		dgm.geteff("ef10",  "e10", "c00");
		dgm.geteff("ef11",  "e11", "c01");
		dgm.geteff("ef12",  "e12", "c02");
		dgm.geteff("ef13",  "e13", "c05");
		dgm.geteff("ef14",  "e14", "c06");
		dgm.geteff("ef15",  "e15", "c07");
//		dg0.getH1F("h1gx").add(dg0.getH2F("h2g").projectionY());
//		dg0.getH1F("h1gy").add(dg0.getH2F("h2g").projectionX());		
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    } 
  
    @Override
    public void timerUpdate() {  

    }


}
