package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.analysis.ECsf.SFFunction;
import org.clas.tools.EBMCEngine;
import org.clas.tools.Event;
import org.clas.tools.ParallelSliceFitter;
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
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc1 extends DetectorMonitor {
	
	Event          ev = new Event();
	EBMCEngine  ebmce = new EBMCEngine();
	
	List<Float>           GEN = new ArrayList<Float>();
	List<Float>           REC = new ArrayList<Float>();   
    HipoDataSync       writer = null;		
	List<Particle>       phot = new ArrayList<Particle>();
	List<Particle>       neut = new ArrayList<Particle>();
       
	String tit = null;
	double ethresh = 0.01;
    
    public ECmc1(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("CLUSTERS","RECGEN","EFFICIENCY","SF");

        this.use123Buttons(true);
        this.useZSliderPane(true);

        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECmc1.localinit()");
        
        engine.init();
        engine.isMC = true;
        engine.setVariation("default");
        engine.setPCALTrackingPlane(9);
        engine.setCalRun(11);  
        
        ebmce.getCCDB(11);
        ebmce.setGeom("2.5");
        ebmce.setGoodPhotons(12);
        ebmce.setMCpid(22);
        ebmce.setMCsec(2);
        ebmce.isMC = true;
        
        tl.setFitData(Fits);
        writer = new HipoDataSync();
        writer.open("/Users/colesmith/CLAS12ANA/ECmc1/photon_mc1.hipo");
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
    
    @Override
    public void createHistos(int run) {  
	    System.out.println("ECmc1:createHistos("+run+")");
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
    		dgm.add("CLUSTERS",3,5,0,st,getRunNumber());    
    		dgm.makeH1("h211", 5,0.5,5.5,-1,           "","PCAL Clusters",1,2,2);                              dgm.cc("h211", true,false,0,0,0,0);  
    		dgm.makeH1("h211t",5,0.5,5.5,-2,           "","PCAL Clusters",1,0,2);                             
    		dgm.makeH1("h221", 5,0.5,5.5,-1,           "","ECIN Clusters",1,5,2);                              dgm.cc("h221", true,false,0,0,0,0);
    		dgm.makeH1("h221t",5,0.5,5.5,-2,           "","ECIN Clusters",1,0,2);                              
    		dgm.makeH1("h231", 5,0.5,5.5,-1,           "","ECOU Clusters",1,4,2);                              dgm.cc("h231", true,false,0,0,0,0);
    		dgm.makeH1("h231t",5,0.5,5.5,-2,           "","ECOU Clusters",1,0,2);                             
    		dgm.makeH2("st1",  5,0.5,5.5,4,0,4,-1,     "","PCAL Clusters","STATUS");                           dgm.cc("st1", false,true,0,0,0,0);
    		dgm.makeH2("st4",  5,0.5,5.5,4,0,4,-1,     "","ECIN Clusters","STATUS");                           dgm.cc("st4", false,true,0,0,0,0);
    		dgm.makeH2("st7",  5,0.5,5.5,4,0,4,-1,     "","ECOU Clusters","STATUS");                           dgm.cc("st7", false,true,0,0,0,0);                     
    		dgm.makeH1("hr21", 50,0,10,-1,            "","Photon Energy (GeV)","Avg. No. PCAL Clusters",1,2);  dgm.cc("hr21",false,false,1,1.4f,0,0);
    		dgm.makeH1("hr22", 50,0,10,-1,            "","Photon Energy (GeV)","Avg. No. ECIN Clusters",1,5);  dgm.cc("hr22",false,false,1,1.4f,0,0);
    		dgm.makeH1("hr23", 50,0,10,-1,            "","Photon Energy (GeV)","Avg. No. ECOU Clusters",1,4);  dgm.cc("hr23",false,false,1,1.4f,0,0);            
    		dgm.makeH2("set1", 50,0,1,4,1.5,5.5,-1,    "","PCAL Energy Fraction","PCAL Clusters");             dgm.cc("set1",false,true,0,0,0,0); 
    		dgm.makeH2("set4", 50,0,1,4,1.5,5.5,-1,    "","ECIN Energy Fraction","ECIN Clusters");             dgm.cc("set4",false,true,0,0,0,0); 
    		dgm.makeH2("set7", 50,0,1,4,1.5,5.5,-1,    "","ECOU Energy Fraction","ECOU Clusters");             dgm.cc("set7",false,true,0,0,0,0);
    		dgm.makeH1("h212", 50,0.,1,-1,             "","PCAL Energy Fraction",1,2);
    		dgm.makeH1("h222", 50,0.,1,-1,             "","ECIN Energy Fraction",1,5);
    		dgm.makeH1("h232", 50,0.,1,-1,             "","ECOU Energy Fraction",1,4);
    	}

    }
    
    public void createRECGEN(int st) {
    	
    	switch (st) {        
        case 0: 
    		dgm.add("RECGEN",4,4,0,st,getRunNumber());    
        	dgm.makeH2("dee1", 50,0,10,20,-0.4,0.4,0,"N1=1",     "Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee1",false,true,0,0,0,0);
        	dgm.makeH2("dee2", 40,0,10,15,-0.4,0.4,0,"N1=1 N4=2","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee2",false,true,0,0,0,0); 
        	dgm.makeH2("dee3", 30,0,10,10,-0.4,0.4,0,"N1=1 N4=3","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee3",false,true,0,0,0,0);
        	dgm.makeH2("dee4", 30,0,10,10,-0.4,0.4,0,"N1=1 N4=4","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee4",false,true,0,0,0,0);
        	dgm.makeH2("det1", 50,0,10,20,-0.6,0.6,0,"N1=1",     "Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det1",false,true,0,0,0,0);
        	dgm.makeH2("det2", 50,0,10,20,-0.6,0.6,0,"N1=1 N4=2","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det2",false,true,0,0,0,0);
        	dgm.makeH2("det3", 50,0,10,20,-0.6,0.6,0,"N1=1 N4=3","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det3",false,true,0,0,0,0);
        	dgm.makeH2("det4", 50,0,10,20,-0.6,0.6,0,"N1=1 N4=4","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det4",false,true,0,0,0,0);
        	dgm.makeH2("dep1", 50,0,10,20,-0.6,0.6,0,"N1=1",     "Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep1",false,true,0,0,0,0);
        	dgm.makeH2("dep2", 50,0,10,20,-0.6,0.6,0,"N1=1 N4=2","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep2",false,true,0,0,0,0);
        	dgm.makeH2("dep3", 50,0,10,20,-0.6,0.6,0,"N1=1 N4=3","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep3",false,true,0,0,0,0);
        	dgm.makeH2("dep4", 50,0,10,20,-0.6,0.6,0,"N1=1 N4=4","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep4",false,true,0,0,0,0);
        	dgm.makeH2("deb1", 50,0,10,20,-0.05,0.05,0,"N1=1",     "Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb1",false,true,0,0,0,0);
        	dgm.makeH2("deb2", 50,0,10,20,-0.05,0.05,0,"N1=1 N4=2","Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb2",false,true,0,0,0,0);
        	dgm.makeH2("deb3", 50,0,10,20,-0.05,0.05,0,"N1=1 N4=3","Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb3",false,true,0,0,0,0);
        	dgm.makeH2("deb4", 50,0,10,20,-0.05,0.05,0,"N1=1 N4=4","Photon Energy (GeV)","#Delta#beta");      dgm.cc("deb4",false,true,0,0,0,0);
        	break;
        case 1: 
    		dgm.add("RECGEN",4,4,0,st,getRunNumber());    
        	dgm.makeH2("dte1", 50,5,26,20,-0.4,0.4,0,"N1=1",     "Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte1",false,true,0,0,0,0);
        	dgm.makeH2("dte2", 40,5,26,15,-0.4,0.4,0,"N1=1 N4=2","Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte2",false,true,0,0,0,0); 
        	dgm.makeH2("dte3", 30,5,26,10,-0.4,0.4,0,"N1=1 N4=3","Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte3",false,true,0,0,0,0);
        	dgm.makeH2("dte4", 30,5,26,10,-0.4,0.4,0,"N1=1 N4=4","Photon Theta (deg)","#DeltaE/E");          dgm.cc("dte4",false,true,0,0,0,0);
        	dgm.makeH2("dtt1", 50,5,26,20,-0.6,0.6,0,"N1=1",     "Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt1",false,true,0,0,0,0);
        	dgm.makeH2("dtt2", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=2","Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt2",false,true,0,0,0,0);
        	dgm.makeH2("dtt3", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=3","Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt3",false,true,0,0,0,0);
        	dgm.makeH2("dtt4", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=4","Photon Theta (deg)","#Delta#theta (deg)"); dgm.cc("dtt4",false,true,0,0,0,0);
        	dgm.makeH2("dtp1", 50,5,26,20,-0.6,0.6,0,"N1=1",     "Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp1",false,true,0,0,0,0);
        	dgm.makeH2("dtp2", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=2","Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp2",false,true,0,0,0,0);
        	dgm.makeH2("dtp3", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=3","Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp3",false,true,0,0,0,0);
        	dgm.makeH2("dtp4", 50,5,26,20,-0.6,0.6,0,"N1=1 N4=4","Photon Theta (deg)","#Delta#phi (deg)");   dgm.cc("dtp4",false,true,0,0,0,0);
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
    		tit = "N#gamma>0, thresh = "+ethresh*1e3+" MeV";
        	dgm.makeH1("ef11",50,0,10,-1,tit,              "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef13",50,0,10,-2,"n>0 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef14",50,0,10,-2,"n>0 layer 1,4",  "Photon Energy (GeV)",1,5);
        	dgm.makeH1("ef15",50,0,10,-2,"n>0 layer 1,4,7","Photon Energy (GeV)",1,4);
        	tit = "N#gamma=1, thresh = "+ethresh*1e3+" MeV";
        	dgm.makeH1("ef21",50,0,10,-1,tit,              "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef22",50,0,10,-2,"n=1 any layer",  "Photon Energy (GeV)",1,1);
        	dgm.makeH1("ef23",50,0,10,-2,"n=1 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef24",50,0,10,-2,"n=1 layer 1,4",  "Photon Energy (GeV)",1,5);
        	dgm.makeH1("ef25",50,0,10,-2,"n=1 layer 1,4,7","Photon Energy (GeV)",1,4);
            break;
        case 1:            
    		dgm.add("EFFICIENCY",3,4,0,st,getRunNumber()); tit = "n>0 any layer, thresh = "+ethresh+" MeV";
        	dgm.makeH1("h10", 50,0,10,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h11", 50,0,10,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h12", 50,0,10,-1,"PCAL","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h13", 50,0,10,-1,"PCAL","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h14", 50,0,10,-1,"ECIN","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h15", 50,0,10,-1,"ECOU","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h21", 50,0,10,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h22", 50,0,10,-1,"ECIN","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h23", 50,0,10,-1,"ECOU","Photon Energy (GeV)",1,5);
        	dgm.makeH1("h210",50,0,10,-1,"PCAL","Photon Energy (GeV)",1,3);
        	dgm.makeH1("h220",50,0,10,-1,"ECIN","Photon Energy (GeV)",1,2);
        	dgm.makeH1("h230",50,0,10,-1,"ECOU","Photon Energy (GeV)",1,5);
    	}

    }
        
    public void createSF(int st) {
    	
    	switch (st) {
    	case 0: 
    		dgm.add("SF", 2, 2, 0, st, getRunNumber());
    		dgm.makeH2("sf1",80,0,2.8,50,0.15,0.35,-1,"N#gamma>0","Reconstructed Photon Energy (GeV)","Sampling Fraction");           dgm.cc("sf1",false,true,0,0,0,0);
    		dgm.makeH2("sf2",80,0,2.8,50,0.15,0.35,-1,"N#gamma>0 #theta<10","Reconstructed Photon Energy (GeV)","Sampling Fraction"); dgm.cc("sf2",false,true,0,0,0,0);
    		dgm.makeH2("sf3",80,0,2.8,50,0.15,0.35,-1,"N#gamma=1","Reconstructed Photon Energy (GeV)","Sampling Fraction");           dgm.cc("sf3",false,true,0,0,0,0);
    		dgm.makeH2("sf4",80,0,2.8,50,0.15,0.35,-1,"N#gamma=1 #theta<10","Reconstructed Photon Energy (GeV)","Sampling Fraction"); dgm.cc("sf4",false,true,0,0,0,0);  		
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
    	plotSFfit("SF","sf1","sf2","sf3","sf4");
    }
    
    public void plotMCHistos() {      
        plot("CLUSTERS");
        plot("RECGEN");
        plot("EFFICIENCY");  
        plot("SF");
    }
    
    @Override
    public void plotEvent(DataEvent de) {
        analyze();         
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
    	geteff();
    	getSFfit("sf1","sf2","sf3","sf4");
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
		    if(p.hasProperty("beta")) out.add(n++,(float) p.getProperty("beta"));
    	}
		return out;
    }
    
    @Override
    public void processEvent(DataEvent de) {
    	
    	GEN.clear();  REC.clear();
    	
    	DetectorParticle p1 = new DetectorParticle();
    	DetectorParticle p2 = new DetectorParticle();  
    	
    	List<Integer> s1 = new ArrayList<Integer>();
    	List<Integer> s4 = new ArrayList<Integer>();
    	List<Integer> s7 = new ArrayList<Integer>();
    	
    	int n1=0, n4=0, n7=0, npp=0, npart=0;
     	double e1=0, e4=0, e7=0, e1p=0, e4p=0, e7p=0;     
     	
		boolean goodev = ebmce.readMC(de) && ebmce.pmc.size()==1;
		
        if (goodev) { 
        	
        	GEN = getkin(ebmce.pmc); float refP = (float) GEN.get(0); float refTH = (float) GEN.get(1);        	    	
 	
//	    	if(!ebmce.hasStartTime) dropBanks(de);  //drop ECAL banks and re-run ECEngine            
        	if(!ebmce.processDataEvent(de)) return;
        	
        	float stt = ebmce.starttime;  
        	
            List<DetectorParticle> par = ebmce.eb.getEvent().getParticles(); 
        	List<DetectorResponse> cal = ebmce.eb.getEvent().getCalorimeterResponseList(); 
        	
        	if(dumpFiles) {ebmce.getRECBanks(de,ebmce.eb); writer.writeEvent(de);}
        	
        	int trsec = -1; int trpid = -211;
        	for (DetectorParticle dp: par) {
        		int pid = dp.getPid(); int sec = dp.getSector(DetectorType.ECAL);
        		if(trsec==-1 && sec>0 && pid==trpid) trsec=sec;
        	}
        	
        	if(trsec==-1) return;   
        	
        	if (dbgAnalyzer) {
        	System.out.println(" ");
        	System.out.println(getEventNumber()+" "+cal.size()+" "+par.size());
        	
        	for (DetectorResponse drr : cal) {
        		CalorimeterResponse dr = (CalorimeterResponse) drr;
        		System.out.println(dr.getAssociation()+" "+dr.getDescriptor().getType()+" "+dr.getDescriptor().getSector()+" "+dr.getDescriptor().getLayer()+" "
        	                      +par.get(dr.getAssociation()).getPid());
        	}
        	
        	int nnn=0;
        	for (DetectorParticle dp : par) {
        		System.out.println("Particle "+nnn+"  "+dp.getSector(DetectorType.ECAL)+" "+dp.getEnergy(DetectorType.ECAL));nnn++;
        		for (DetectorResponse dr : dp.getDetectorResponses()) {
              	  System.out.println(dr.getAssociation()+" "+dr.getDescriptor().getType()+" "+dr.getDescriptor().getLayer());       			
        		}
        	}
        	}
        	       	
        	List<Particle> plist = new ArrayList<Particle>(); 
      	
        	for (DetectorParticle dp : par) { // make list of neutral Particle objects 
//        		System.out.println(trsec+" "+dp.getSector(DetectorType.ECAL)+" "+ebmce.hasStartTime+" "+dp.getPid()+" "+dp.getBeta()+" "+dp.getEnergy(DetectorType.ECAL));
		        if(dp.getSector(DetectorType.ECAL)!=trsec && dp.getPid()==22) { npart++;
			    if(!ebmce.hasStartTime && dp.getPid()==2112) {// this repairs zero momentum neutrons from non-PCAL seeded neutrals
	 				double e = dp.getEnergy(DetectorType.ECAL)/ebmce.getSF(dp); 		
			    	Vector3D vec = new Vector3D() ; vec.copy(dp.getHit(DetectorType.ECAL).getPosition()); vec.unit(); 			    		
 			    	dp.vector().add(new Vector3(e*vec.x(),e*vec.y(),e*vec.z())); //track energy for neutrals in DetectorParticle
 			    	dp.setPid(22);
 			    }
		    	//SF corrected Particle energy from DetectorParticle
			    Particle p = dp.getPhysicsParticle(22); p.setProperty("beta",dp.getBeta()); plist.add(p);
		        }
 			}
        	
        	for (DetectorParticle dp : par) { // make list of neutral Particle objects 
        		float em = (float) (dp.getEnergy(DetectorType.ECAL));
        		float  e = (float) (em/ebmce.getSF(dp));
        		if(dp.getSector(DetectorType.ECAL)==2 && e>ethresh && dp.getPid()==22) {
        			npp++; dgm.fill("sf1", em, em/refP); if(refTH<8) {dgm.fill("sf2", em, em/refP);
        		}
        			for (DetectorResponse dr : dp.getDetectorResponses()) {
        				if(dr.getDescriptor().getType()==DetectorType.ECAL && dr.getDescriptor().getSector()==2) {		        		
        					if(e>ethresh) {
        						for (DetectorResponse drr : cal){
        							CalorimeterResponse dcr = (CalorimeterResponse) drr; 			              		        			
        							boolean go = dcr.getAssociation()==dr.getAssociation() && 
        										 dcr.getDescriptor().getLayer()==dr.getDescriptor().getLayer() &&
        										 dcr.getDescriptor().getSector()==2;
        							int lay = go ? dr.getDescriptor().getLayer():0; 	
        							if(lay==1) {n1++; e1p+=e; s1.add(dr.getStatus()); if(n1==1) {e1=e; dgm.fill("h210",refP);}}
        							if(lay==4) {n4++; e4p+=e; s4.add(dr.getStatus()); if(n4==1) {e4=e; dgm.fill("h220",refP);}}
        							if(lay==7) {n7++; e7p+=e; s7.add(dr.getStatus()); if(n7==1) {e7=e; dgm.fill("h230",refP);}}
        						}
        					}		        		
        				}
        			}
        		}
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
        		if(n1==1)       {dgm.fill("dee1",refP,delE1);dgm.fill("det1",refP,delTH1); dgm.fill("dep1",refP,delPH1);dgm.fill("deb1",refP,delBET);}
        		if(n1==1&&n4>1) {dgm.fill("dee2",refP,delE1);dgm.fill("det2",refP,delTH1); dgm.fill("dep2",refP,delPH1);dgm.fill("deb2",refP,delBET);}
        		if(n1==1&&n4>2) {dgm.fill("dee3",refP,delE1);dgm.fill("det3",refP,delTH1); dgm.fill("dep3",refP,delPH1);dgm.fill("deb3",refP,delBET);}   
        		if(n1==1&&n4>3) {dgm.fill("dee4",refP,delE1);dgm.fill("det4",refP,delTH1); dgm.fill("dep4",refP,delPH1);dgm.fill("deb4",refP,delBET);} 
        		
        		//RECGEN.1
        		if(refP>5) {
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
        		//CLUSTERS.0.3
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
        			if(dp.getSector(DetectorType.ECAL)==2 && e>ethresh && dp.getPid()==22) {         				
	        			for (DetectorResponse dr : dp.getDetectorResponses()) {
	        				int lay = dr.getDescriptor().getLayer();
	        				if(dr.getDescriptor().getType()==DetectorType.ECAL) {	        					
	        					if(lay==1) sum+= 100;
	        					if(lay==4) sum+=  40;
	        					if(lay==7) sum+=   7;
	        				}
	        			}	        			
	        			if(sum>=100) {dgm.fill("h13",refP); dgm.fill("sf3", em, em/refP); if(refTH<8) dgm.fill("sf4", em, em/refP);}
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
    }
    
    public void getSFfit(String...  h2) {
    	
    	int run = getRunNumber(), n=0;
    	ParallelSliceFitter fitter;
    	
        float min=0.01f, max=2.4f;
        
    	for (String h : h2) {
    		fitter = new ParallelSliceFitter(dgm.getH2F(h));
    		fitter.setBackgroundOrder(0); fitter.setMin(0.16); fitter.setMax(0.30); fitter.fitSlicesX(); 
    		GraphErrors MeanGraph = fitter.getMeanSlices();      
    		FitSummary.add(MeanGraph,n, 0, 7, run);      
    		tl.fitData.add(fitEngine(MeanGraph,15,min,max,min,max),n,0,5,run);
    		n++;
    	}
	    
    }
    
    public void plotSFfit(String tab, String... h2) {
    	
        int run = getRunNumber(), n=0;
        
   	    int index = getDetectorTabNames().indexOf(tab);    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); c.divide(2, 2);

        SFFunction sf = new SFFunction("gsf",22,0,ebccdb,0.01,2.5);  sf.setLineWidth(2) ; sf.setLineColor(2); sf.setOptStat("");
    	
   		for (String h : h2) {
        if (FitSummary.hasItem(n,0,7,run)) {
        	c.cd(n);c.getPad().getAxisZ().setLog(true); c.draw(dgm.getH2F(h));
            GraphPlot((GraphErrors)FitSummary.getItem(n,0,7,run),c,n,0.0f,2.8f,0.15f,0.35f,4,6,1,""," E/P","same"); c.draw(sf,"same");
            tl.fitData.getItem(n,0,5,run).graph.getFunction().setLineColor(1); 
            tl.fitData.getItem(n,0,5,run).graph.getFunction().setLineWidth(2);
            tl.fitData.getItem(n,0,5,run).graph.getFunction().setOptStat("1110");
            tl.fitData.getItem(n,0,5,run).graph.getFunction().setRange(0.01,2.8);
//            c.draw(tl.fitData.getItem(n,0,5,run).graph.getFunction(),"same");
            n++;
        } 
   		}
     	
    }      
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(),getRunNumber());
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
