package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.tools.EBMC;
import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc1 extends DetectorMonitor {
	
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
    
    public ECmc1(String name) {
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
        System.out.println("ECmc1.localinit()");
        
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
    	System.out.println("ECmc1.localclear()");
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
	    System.out.println("ECmc1:createHistos("+run+")");
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;    	
    	createCLUSTERS(0);
    	createGENREC(0);
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
    		dgm.makeH2("st1",  5,0.5,5.5,4,0,4,-1,     "","PCAL Clusters","STATUS");                           dgm.cc("st1", false,true,0,0,1,2000);
    		dgm.makeH2("st4",  5,0.5,5.5,4,0,4,-1,     "","ECIN Clusters","STATUS");                           dgm.cc("st4", false,true,0,0,1,2000);
    		dgm.makeH2("st7",  5,0.5,5.5,4,0,4,-1,     "","ECOU Clusters","STATUS");                           dgm.cc("st7", false,true,0,0,1,2000);                     
    		dgm.makeH1("hr21", 50,0,3.8,-1,            "","Photon Energy (GeV)","Avg. No. PCAL Clusters",1,2); dgm.cc("hr21",false,false,1,1.2f,0,0);
    		dgm.makeH1("hr22", 50,0,3.8,-1,            "","Photon Energy (GeV)","Avg. No. ECIN Clusters",1,5); dgm.cc("hr22",false,false,1,1.2f,0,0);
    		dgm.makeH1("hr23", 50,0,3.8,-1,            "","Photon Energy (GeV)","Avg. No. ECOU Clusters",1,4); dgm.cc("hr23",false,false,1,1.2f,0,0);            
    		dgm.makeH2("set1", 50,0,1,4,1.5,5.5,-1,    "","PCAL Energy Fraction","PCAL Clusters");             dgm.cc("set1",false,true,0,0,1,100); 
    		dgm.makeH2("set4", 50,0,1,4,1.5,5.5,-1,    "","ECIN Energy Fraction","ECIN Clusters");             dgm.cc("set4",false,true,0,0,1,100); 
    		dgm.makeH2("set7", 50,0,1,4,1.5,5.5,-1,    "","ECOU Energy Fraction","ECOU Clusters");             dgm.cc("set7",false,true,0,0,1,100);
    		dgm.makeH1("h212", 50,0.,1,-1,             "","PCAL Energy Fraction",1,2);
    		dgm.makeH1("h222", 50,0.,1,-1,             "","ECIN Energy Fraction",1,5);
    		dgm.makeH1("h232", 50,0.,1,-1,             "","ECOU Energy Fraction",1,4);
    	}

    }
    
    public void createGENREC(int st) {
    	
    	switch (st) {        
        case 0: 
    		dgm.add("GENREC",4,3,0,st,getRunNumber());    
        	dgm.makeH2("dee1", 50,0,3.8,20,-0.4,0.4,0,"N#gamma=1","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee1",false,true,0,0,1,100);
        	dgm.makeH2("dee2", 40,0,3.8,15,-0.4,0.4,0,"N#gamma=2","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee2",false,true,0,0,1,100); 
        	dgm.makeH2("dee3", 30,0,3.8,10,-0.4,0.4,0,"N#gamma=3","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee3",false,true,0,0,1,100);
        	dgm.makeH2("dee4", 30,0,3.8,10,-0.4,0.4,0,"N#gamma=4","Photon Energy (GeV)","#DeltaE/E");          dgm.cc("dee4",false,true,0,0,1,100);
        	dgm.makeH2("det1", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=1","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det1",false,true,0,0,1,100);
        	dgm.makeH2("det2", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=2","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det2",false,true,0,0,1,100);
        	dgm.makeH2("det3", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=3","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det3",false,true,0,0,1,100);
        	dgm.makeH2("det4", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=4","Photon Energy (GeV)","#Delta#theta (deg)"); dgm.cc("det4",false,true,0,0,1,100);
        	dgm.makeH2("dep1", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=1","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep1",false,true,0,0,1,100);
        	dgm.makeH2("dep2", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=2","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep2",false,true,0,0,1,100);
        	dgm.makeH2("dep3", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=3","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep3",false,true,0,0,1,100);
        	dgm.makeH2("dep4", 50,0,3.8,20,-0.6,0.6,0,"N#gamma=4","Photon Energy (GeV)","#Delta#phi (deg)");   dgm.cc("dep4",false,true,0,0,1,100);
    	}
    }
    	
    public void createEFFICIENCY(int st) {
        
    	switch (st) {        
        case 0:                        
    		dgm.add("EFFICIENCY",2,2,0,st,getRunNumber()); 
    		tit = "n>0 any layer, thresh = "+ethresh*1e3+" MeV";
        	dgm.makeH1("ef11",50,0,3.8,-1,tit,              "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef13",50,0,3.8,-2,"n>0 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef14",50,0,3.8,-2,"n>0 layer 1,4",  "Photon Energy (GeV)",1,5);
        	dgm.makeH1("ef15",50,0,3.8,-2,"n>0 layer 1,4,7","Photon Energy (GeV)",1,4);
        	tit = "n=1 each layer, thresh = "+ethresh*1e3+" MeV";
        	dgm.makeH1("ef21",50,0,3.8,-1,tit,              "Photon Energy (GeV)",1,3);
        	dgm.makeH1("ef22",50,0,3.8,-2,"n=1 any layer",  "Photon Energy (GeV)",1,1);
        	dgm.makeH1("ef23",50,0,3.8,-2,"n=1 layer 1",    "Photon Energy (GeV)",1,2);
        	dgm.makeH1("ef24",50,0,3.8,-2,"n=1 layer 1,4",  "Photon Energy (GeV)",1,5);
        	dgm.makeH1("ef25",50,0,3.8,-2,"n=1 layer 1,4,7","Photon Energy (GeV)",1,4);
            break;
        case 1:            
    		dgm.add("EFFICIENCY",4,3,0,st,getRunNumber()); tit = "n>0 any layer, thresh = "+ethresh+" MeV";
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
    	
		GEN.clear(); REC1.clear(); 
		
    	List<Integer> s1 = new ArrayList<Integer>();
    	List<Integer> s4 = new ArrayList<Integer>();
    	List<Integer> s7 = new ArrayList<Integer>();
    	
    	int n1=0, n4=0, n7=0, npp=0, npart=0;
     	double e1=0, e4=0, e7=0, e1p=0, e4p=0, e7p=0;
		
		boolean debug = false;
		if(debug) {System.out.println(" ") ; System.out.println(getEventNumber());}
    	
		boolean goodev = eb.readMC(de) && eb.pmc.size()==1;
			    	
        if (goodev) {    
        	
        	GEN = getkin(eb.pmc);   
        	float refP = (float) GEN.get(0);        	

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
   	    	    if(e>ethresh) {
   	    	    	npp++;
   	    	    	for (DetectorResponse dr : dp.getDetectorResponses()) {
   	    	    		int lay = dr.getDescriptor().getLayer();    		
   	    	    		if(lay==1) {n1++; e1p+=e; s1.add(dr.getStatus()); if(n1==1) {e1=e; dgm.fill("h210",refP);}}
   	    	    		if(lay==4) {n4++; e4p+=e; s4.add(dr.getStatus()); if(n4==1) {e4=e; dgm.fill("h220",refP);}}
   	    	    		if(lay==7) {n7++; e7p+=e; s7.add(dr.getStatus()); if(n7==1) {e7=e; dgm.fill("h230",refP);}}
   	    	    	}
    	    	}			    
 			}  
       	    
        	dgm.fill("h10",refP);	      	        	
        	
        	if(npart>0) { // 1 or more neutral particles 
         		
           		REC1 = getkin(plist);
        		double  delE1 = (GEN.get(0)-REC1.get(0))/refP;
        		double delTH1 = (GEN.get(1)-REC1.get(1));
        		double delPH1 = (GEN.get(2)-REC1.get(2));
        		
        		dgm.fill("h11",refP);       	   
        	    
        		if(n1==1)       {dgm.fill("dee1",refP,delE1);dgm.fill("det1",refP,delTH1); dgm.fill("dep1",refP,delPH1);}
        		if(n1==1&&n4>1) {dgm.fill("dee2",refP,delE1);dgm.fill("det2",refP,delTH1); dgm.fill("dep2",refP,delPH1);}
        		if(n1==1&&n4>2) {dgm.fill("dee3",refP,delE1);dgm.fill("det3",refP,delTH1); dgm.fill("dep3",refP,delPH1);}   
        		if(n1==1&&n4>3) {dgm.fill("dee4",refP,delE1);dgm.fill("det4",refP,delTH1); dgm.fill("dep4",refP,delPH1);}   
        		dgm.fill("h211",n1); dgm.fill("h21",refP,n1); for (Integer i : s1) dgm.fill("st1",n1,i);
        		dgm.fill("h221",n4); dgm.fill("h22",refP,n4); for (Integer i : s4) dgm.fill("st4",n4,i);
        		dgm.fill("h231",n7); dgm.fill("h23",refP,n7); for (Integer i : s7) dgm.fill("st7",n7,i);
        		if(n1>0) dgm.fill("h211t", npart); 
        		if(n4>0) dgm.fill("h221t", npart);  
        		if(n7>0) dgm.fill("h231t", npart);         		        	    
                if(n1>1) dgm.fill("h212",e1/e1p); dgm.fill("set1",e1/e1p,n1);
                if(n4>1) dgm.fill("h222",e4/e4p); dgm.fill("set4",e4/e4p,n4);
                if(n7>1) dgm.fill("h232",e7/e7p); dgm.fill("set7",e7/e7p,n7);
        	}
        	
        	if(npp==1) { // only 1 neutral particle above threshold           		
        		dgm.fill("h12",refP); 
        	    for (DetectorParticle dp : np) {
        	    	double e = dp.getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, dp, eb.ccdb);
        	    	int sum=0;
           	    	if(e>ethresh) {
           	    		for (DetectorResponse dr : dp.getDetectorResponses()) {
           	    			int lay = dr.getDescriptor().getLayer();
           	    			if(lay==1) sum+= 100;
           	    			if(lay==4) sum+=  40;
           	    			if(lay==7) sum+=   7;
           	    		}
           	    		if(sum>=100) dgm.fill("h13",refP);
           	    		if(sum>=140) dgm.fill("h14",refP);
           	    		if(sum==147) dgm.fill("h15",refP);
        	    	}
        		}
        	}
        }  		
      
    }
   
    public void plotMCHistos() {      
        plot("CLUSTERS");
        plot("GENREC");
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
    	dgm.drawGroup(tabname,0,getActive123(),getRunNumber());
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
