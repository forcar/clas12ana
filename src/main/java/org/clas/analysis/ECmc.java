package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.clas.analysis.ECPart.SFFunction;
import org.clas.tools.EBMCEngine;
import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECmc extends DetectorMonitor {
	
	Event          ev = new Event();
    IndexGenerator ig = new IndexGenerator();
	   
    float dp1=0.1f, dp2=0.1f; 
    
    Boolean isGEMC = true, isGENREC = false;
		
    public ECmc(String name) {
        super(name);
        
        dgmActive=true; 
        
        if(isGENREC) setDetectorTabNames("GENREC","KINEMATICS","EFFICIENCY");
        if(isGEMC)   setDetectorTabNames("GEMC");
        
        use123Buttons(true);
        useSliderPane(true);
        useECEnginePane(true);

        init();
        localinit("rga_fall2018");
        localclear();
    }
    
    @Override
    public void localinit(String variation) {
        System.out.println("ECmc.localinit("+variation+")");
        eng.engine.init();
        eng.engine.setIsMC(true);
        eng.engine.setGeomVariation(variation);
        eng.engine.setVariation("default");  
        eng.setEngineConfig("phot");
      
        tl.setFitData(Fits);
    }
    
    @Override
    public void localclear() {
    	System.out.println("ECmc.localclear()");
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
    	histosExist = true; 
    	if(isGEMC)    createGEMC(0);
    	if(isGENREC) {createGENREC(0); createKINEMATICS(0); for (int i=0; i<3; i++) createEFFICIENCY(i);}
    }
    
    public void createGEMC(int st) {
    	
        String[] partname = {"ELECTRON", "PROTON", "PION+", "PHOTON", "NEUTRON"};
    	switch (st) {        
        case 0: 
        	dgm.add("GEMC",5,3,0,st,getRunNumber());
        	for (int i=1; i<6; i++) {
        		dgm.makeH1("g0"+i,50,0,9,-1,partname[i-1],"MC Momentum (GeV)",1,4,"1000000");
        		dgm.makeH1("g1"+i,50,0,9,-2,partname[i-1]+" PCAL","MC Momentum (GeV)",1,1,"1000000");
        		dgm.makeH1("g2"+i,50,0,9,-2,partname[i-1]+" ECIN","MC Momentum (GeV)",1,3,"1000000");
        		dgm.makeH1("g3"+i,50,0,9,-2,partname[i-1]+" ECOU","MC Momentum (GeV)",1,5,"1000000");
//        		dgm.makeGE("e1"+i,-1,"","EFF",)
        	}
    	}
    }
    
    public void createGENREC(int st) {

    	switch (st) {        
        case 0: 
        	dgm.add("GENREC",4,5,0,st,getRunNumber());
        	dgm.makeH1("h00",50,-0.1,0.1,-1,"E#gamma<10 MeV","#DeltaPx GEN-REC (GeV)",1,4,"1000100"); 
        	dgm.makeH1("h01",50,-0.1,0.1,-1,"E#gamma<10 MeV","#DeltaPy GEN-REC (GeV)",1,4,"1000100"); 
        	dgm.makeH1("h02",50,-0.1,0.1,-1,"E#gamma<10 MeV","#DeltaPz GEN-REC (GeV)",1,4,"1000100"); 
        	dgm.makeH1("h03",50,-0.1,0.1,-1,"E#gamma<10 MeV","#DeltaP  GEN-REC (GeV)",1,4,"1000100");
        	
        	dgm.makeH2("h10",50,0,5.0, 50,-0.1,1.0,0,"E#gamma>10 MeV","ECAL #gamma E (GeV)","#Deltapx (GeV)");
        	dgm.makeH2("h11",50,0,5.0, 50,-0.1,0.1,0,"E#gamma>10 MeV","ECAL #gamma E (GeV)","#Deltapy (GeV)");        	
        	dgm.makeH2("h12",50,0,5.0, 50,-0.1,2.2,0,"E#gamma>10 MeV","ECAL #gamma E (GeV)","#Deltapz (GeV)");
        	dgm.makeH2("h13",50,0,5.0, 50,-0.1,2.2,0,"E#gamma>10 MeV","ECAL #gamma E (GeV)","#Deltap  (GeV)");
        
        	dgm.makeH2("hv0",50,0,4.0, 50,-3,3,    0,"E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVx (cm)");
        	dgm.makeH2("hv1",50,0,4.0, 50,-3,3,    0,"E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVy (cm)");
        	dgm.makeH2("hv2",50,0,4.0, 50,-3,3,    0,"E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVz (cm)");
        	dgm.makeH2("h23",50,0,  9,100, 0.0,5.0,0,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","ECAL #gamma E (GeV)");   

        	
        	dgm.makeH2("hv3",50,-6, 0, 60,-20,20,  0,"#Deltap>"+dp1,"Vz GEN (cm)","z0 (cm)");        	
        	dgm.makeH2("hv4",50, 5,35, 60,-20,20,  0,"#Deltap>"+dp1,"True Electron Theta (deg)","z0 (cm)");
        	dgm.makeH2("hv5",50, 0, 9, 60,-20,20,  0,"#Deltap>"+dp1,"True Electron Momentum (GeV)","z0 (cm)");
        	
        	dgm.makeH2("h20",50,0,100,100,-0.2,5.0,0,"",            "Distance e-#gamma (cm)","#Deltap (GeV)");   
        	dgm.makeH2("h21",50,0,100,100,-0.2,5.0,0,"#Deltap>"+dp1,"Distance e-#gamma (cm)","#Deltap (GeV)");         	
        	dgm.makeH2("h22",50,0,100,100,-0.2,5.0,0,"#Deltap<"+dp2,"Distance e-#gamma (cm)","#Deltap (GeV)");  
        	
    	}        
    }
    
    public void createKINEMATICS(int st) {
     
    	switch (st) {        
        case 0: 
        	dgm.add("KINEMATICS",4,2,0,st,getRunNumber());      
        	String lab[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};
        	dgm.makeH2("h30",50,0,9, 70,-2,12,0,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","PCAL "+lab[0]);   
        	dgm.makeH2("h31",50,0,9, 70,-5,30,0,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","PCAL "+lab[1]);   
        	dgm.makeH2("h32",50,0,9, 70,-2,12,0,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","ECIN "+lab[0]);   
        	dgm.makeH2("h33",50,0,9, 70,-5,30,0,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","ECIN "+lab[1]);
        	
        	dgm.makeH2("h40",50,0,9, 70,-2,12,0,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","PCAL "+lab[0]);   
        	dgm.makeH2("h41",50,0,9, 70,-5,30,0,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","PCAL "+lab[1]);   
        	dgm.makeH2("h42",50,0,9, 70,-2,12,0,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","ECIN "+lab[0]);   
        	dgm.makeH2("h43",50,0,9, 70,-5,30,0,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","ECIN "+lab[1]);
    	}      
    }
    
    public void createEFFICIENCY(int st) {
    	
    	switch (st) {        
        case 0: 
        dgm.add("EFFICIENCY",4,3,0,st,getRunNumber());               
        dgm.makeH2("h50",50,0.0,10., 60,0.0,0.31,-1,"","True Electron Momentum (GeV)","Sampling Fraction E/P");   
        dgm.makeH2("h51",50,0.0,10., 60,0.0,0.31,-1,"","Track Electron Momentum (GeV)","Sampling Fraction E/P");   
        dgm.makeH2("h52",50,0.0,10., 60,0.0,0.31,-1,"","Track Electron Momentum (GeV)","Sampling Fraction E/P");   
        dgm.makeH2("h53",50,0.0,2.5, 60,0.0,0.31,-1,"","ECAL Electron Energy (GeV)","Sampling Fraction E/P");
//        dgm.makeGE("g53",-2,"","","",1); SFFunction sf = new SFFunction("esf",-11,ebmce.ccdb,0.1,2.5); dgm.addDataSet(sf,-2); cannot save to histos
        dgm.makeH1("h60a",50,0,10,-1,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Energy (GeV)","Efficiency #theta>15",1,3);
        dgm.makeH1("h60b",50,0,10,-2,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Energy (GeV)","Efficiency #theta>15",1,5);
        dgm.makeH1("h60c",50,0,10,-2,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Energy (GeV)","Efficiency #theta>15",1,2);
        dgm.makeH1("h61a",50,5,35,-1,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Theta (deg)","Efficiency",1,3);
        dgm.makeH1("h61b",50,5,35,-2,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Theta (deg)","Efficiency",1,5);
        dgm.makeH1("h61c",50,5,35,-2,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Theta (deg)","Efficiency",1,2);
        dgm.makeH1("h62a",50,0,10,-1,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,3);
        dgm.makeH1("h62b",50,0,10,-2,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,5);
        dgm.makeH1("h62c",50,0,10,-2,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,2);
        dgm.makeH1("h62d",50,0,10,-2,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,1);
        dgm.makeH1("h63a",50,5,35,-1,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Theta (deg)","Efficiency",1,3);
        dgm.makeH1("h63b",50,5,35,-2,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Theta (deg)","Efficiency",1,5);
        dgm.makeH1("h63c",50,5,35,-2,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Theta (deg)","Efficiency",1,2);
        dgm.makeH2("h70",30,0,10, 30,5,30,-1,"e-,#gamma,#pi-","True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("h71",30,0,10, 30,5,30,-1,"#chiPID<3.5",   "True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("h72",30,0,10, 30,5,30,-1,"#DeltaP<"+dp2,  "True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH1("h73a",10,0.5,10.5,-1,"","Multiplicity",1, 1,"1000000");
        dgm.makeH1("h73b",10,0.5,10.5,-2,"","Multiplicity",1,34,"1000000");
        dgm.makeH1("h73c",10,0.5,10.5,-2,"","Multiplicity",1, 2,"1000000");
        break;
        case 1:
        dgm.add("EFFICIENCY",4,3,0,st,getRunNumber());               
        dgm.makeH2("eff1", 30,0,10, 30,5,30,-1,"","True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("eff1a",30,0,10, 30,5,30,-1,"","True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("eff1b",30,0,10, 30,5,30,-1,"","True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("eff2a",30,0,10, 30,5,30,-1,"","True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("eff2b",30,0,10, 30,5,30,-1,"","True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("eff3a",30,0,10, 30,5,30,-1,"","True Electron Momentum (GeV)","True Electron Theta (deg)");
        dgm.makeH2("eff3b",30,0,10, 30,5,30,-1,"","True Electron Momentum (GeV)","True Electron Theta (deg)");
        break;
        case 2:
        dgm.add("EFFICIENCY",6,6,0,st,getRunNumber());               
        for (int i=0;i<8;i++) dgm.makeH1("effmom"+i,50,0,10,1,"","MC Electron P (GeV)",4,4);
        for (int i=0;i<8;i++) dgm.makeH1("effthe"+i,50,5,35,1,"","MC Electron #theta (deg)",4,4);
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
  
    public void plotMCHistos() {
    	if(isGEMC)   plot("GEMC");
        if(isGENREC) plot("GENREC"); plot("KINEMATICS"); plot("EFFICIENCY");   
    }
      
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;
    	plotMCHistos();
    }

    @Override
    public void plotEvent(DataEvent de) {
        analyze();         
    }
 
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
//    	fith52(); 
    	if(isGEMC)   geteffGEMC();
    	if(isGENREC) geteffGENREC();
    	isAnalyzeDone = true;
    }
   
    @Override
    public void processEvent(DataEvent event) {     	
    	if(processEV(event)) {
    	    if(isGEMC)   processGEMC(event);  
    	    if(isGENREC) processGENREC(event); 
    	}
    }

    public boolean processEV(DataEvent event) {
    	ev.init(event);        	
    	ev.setEventNumber(getEventNumber());
    	ev.requireOneElectron(false);
   	    ev.setElecTriggerSector(trigFD);
   	    ev.setMCpid(11);
  	    
        return ev.procEvent(event);
    }
 
    public void processGEMC(DataEvent event) {
    	
    	int[] pid = {11,2212,211,22,2112};
        String[] partname = {"ELECTRON", "PROTON", "PION+", "PHOTON", "NEUTRON"};

    	for (int i=0; i<pid.length; i++) {int is = i+1, ip=pid[i];
    		Particle mc = ev.pmc.get(i);
    		dgm.fill("g0"+is, mc.p());
    		for (Map.Entry<Long,List<Particle>>  entry :  ev.filterECALClusters(ev.getECALClusters(ev.getPART(ip,0.0)),2,ip).getMap().entrySet()) {
    			int isp = ig.getIndex(entry.getKey(), 0);
    			 for (Particle rp : entry.getValue()) {	
                     int il = (int) rp.getProperty("layer");                     
                     if(isp==is && il==1) dgm.fill("g1"+is, mc.p());
                     if(isp==is && il==4) dgm.fill("g2"+is, mc.p());
                     if(isp==is && il==7) dgm.fill("g3"+is, mc.p());
    			 }
    		}
    	}   	
    }

    public void processGENREC(DataEvent event) {
 
   	    float refE=0,refP=0,refTH=0,refPH=0;
   	    int nelec=0, nphot=0, npim=0;
        float esum=0,dp=0,chi2pid=1000;
          
        List<Particle> elec = null;
        List<Particle> phot = null;
        List<Particle>  pim = null;
        
    	IndexedList<List<Particle>> ecphot = new IndexedList<List<Particle>>(1);

        ecphot.clear();
        	            
        refE  = (float) ev.pmc.get(0).e(); 
        refP  = (float) ev.pmc.get(0).p(); 
        refTH = (float) Math.toDegrees(ev.pmc.get(0).theta());
        refPH = (float) Math.toDegrees(ev.pmc.get(0).phi());
            
        elec = ev.getPART(11,0,0);  nelec = elec.size(); dgm.fill("h73a",nelec);     
        phot = ev.getPART(22,0.0);  nphot = phot.size(); dgm.fill("h73b",nphot);    
        pim  = ev.getPART(212,0.0); npim  = pim.size();  dgm.fill("h73c",npim);   		
        int npart = nelec+nphot+npim;        		        		
		      
        ecphot = ev.filterECALClusters(ev.getECALClusters(phot),2,22); //used in fillECelec for e-g residuals
	
        dgm.fill("eff1",refE, refTH);
        dgm.fill("effthe0",refTH);  

        if(           refTH>15)  dgm.fill("effmom0",refE);
        if(npart>0 && refTH>15) {dgm.fill("effmom4",refE);dgm.fill("effthe4",refTH); dgm.fill("eff1a",refE,refTH);}
        if(npim==0 && refTH>15) {dgm.fill("effmom5",refE);dgm.fill("effthe5",refTH);}
   		
        if(nelec==1) {
        	
            for (Map.Entry<Long,List<Particle>>  entry : ecphot.getMap().entrySet()) {
                if(ig.getIndex(entry.getKey(), 0)==trigFD) for (Particle rp : entry.getValue()) esum+=rp.e(); 
            }

            chi2pid = (float) Math.abs(elec.get(0).getProperty("chi2pid"));
            dp = (float)(ev.pmc.get(0).p()-elec.get(0).p());
                
            if(npim==0 && chi2pid<3.5)  {
                if(refTH>15)                     dgm.fill("effmom6",refE); 
                if(refTH>15 && Math.abs(dp)<dp2) dgm.fill("effmom7",refE); 
                if(Math.abs(dp)<dp2)       dgm.fill("eff3a",refE,refTH);
                dgm.fill("effthe6",refTH); dgm.fill("eff2a",refE,refTH); 	                    
            }
        		
            if(esum<=0.0) {
                dgm.fill("h00",ev.pmc.get(0).px()-elec.get(0).px());
                dgm.fill("h01",ev.pmc.get(0).py()-elec.get(0).py());
                dgm.fill("h02",ev.pmc.get(0).pz()-elec.get(0).pz());
                dgm.fill("h03",dp);
            }
    			
            if(esum>0.01) {
                dgm.fill("h10",esum,ev.pmc.get(0).px()-elec.get(0).px());
                dgm.fill("h11",esum,ev.pmc.get(0).py()-elec.get(0).py());
                dgm.fill("h12",esum,ev.pmc.get(0).pz()-elec.get(0).pz());
                dgm.fill("h13",esum,dp);
                	
                dgm.fill("hv0",esum,ev.pmc.get(0).vx()-elec.get(0).vx());
                dgm.fill("hv1",esum,ev.pmc.get(0).vy()-elec.get(0).vy());
                dgm.fill("hv2",esum,ev.pmc.get(0).vz()-elec.get(0).vz());
            }
                              		
            float e_mom = (float) elec.get(0).p();
            float e_the = (float) Math.toDegrees(elec.get(0).theta());
            float e_phi = (float) Math.toDegrees(elec.get(0).phi());
                
            double    pt = Math.sqrt(e_mom*e_mom-elec.get(0).pz()*elec.get(0).pz());
            double ptref = Math.sqrt(refP*refP-ev.pmc.get(0).pz()*ev.pmc.get(0).pz());
                
            List<Particle> elecECAL = ev.getECAL((int)elec.get(0).getProperty("pindex"));
       	    float ex=ev.getVar(elecECAL, "x", 1), ey=ev.getVar(elecECAL, "y", 1), ez=ev.getVar(elecECAL, "z", 1);
       			
            for (Map.Entry<Long,List<Particle>>  entry : ecphot.getMap().entrySet()) {
                int is = ig.getIndex(entry.getKey(), 0);
                for (Particle rp : entry.getValue()) {				
                    int il = (int) rp.getProperty("layer");
                    if(rp.getProperty("sector")!=trigFD) continue; //photon must be in e- sector and PCAL
                    float ec_the = (float) Math.toDegrees(rp.theta());
                    float ec_phi = (float) Math.toDegrees(rp.phi()); 
                    float dthe = e_the-ec_the; 
                    float dphi = e_phi-ec_phi; 
                    float dphi_ref = refPH-ec_phi;
           			  
                    if(il==1) {  //use only PCAL and ECIN entries                        
//                      double z0  =    pt*0.667*Math.sin(dphi*3.14159/180.)*100; 
                        double z0  =    ptref*0.667*Math.sin(dphi_ref*3.14159/180.)*100; 
                      
                        if(Math.abs(dp)>dp1) {dgm.fill("hv3",ev.pmc.get(0).vz(),z0);dgm.fill("hv4",refTH,z0);dgm.fill("hv5",refP,z0);}	
              	      
           			    float gx=(float)rp.getProperty("x");float gy=(float)rp.getProperty("y");float gz=(float)rp.getProperty("z");
           			    float dist = (float) Math.sqrt((ex-gx)*(ex-gx)+(ey-gy)*(ey-gy)+(ez-gz)*(ez-gz));
           			  
                        dgm.fill("h20",dist, dp);
                        
           			    if(Math.abs(dp)>dp1) {dgm.fill("h30",e_mom,dthe); dgm.fill("h31",e_mom,dphi); dgm.fill("h21",dist,esum);}
           			    if(Math.abs(dp)<dp2) {dgm.fill("h40",e_mom,dthe); dgm.fill("h41",e_mom,dphi); dgm.fill("h22",dist,esum);dgm.fill("h23",e_mom,esum);}
                    }
                        
                    if(il==4 && Math.abs(dp)>dp1) {dgm.fill("h32",e_mom,dthe); dgm.fill("h33",e_mom,dphi);}
                    if(il==4 && Math.abs(dp)<dp2) {dgm.fill("h42",e_mom,dthe); dgm.fill("h43",e_mom,dphi);}
                }
            }
            
            float en=0; for (Particle p : elecECAL) en += (float) p.getProperty("energy");   
            
                            dgm.fill("h51",e_mom, en*1e-3/e_mom);
            if(chi2pid<3.5) dgm.fill("h52",e_mom, en*1e-3/e_mom);
    	}
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
    
    public void processEBMCE(DataEvent de) { // still debugging this
    	
    	EBMCEngine  ebmce = new EBMCEngine();
    	List<Float>   GEN = new ArrayList<Float>();
    	List<Float>   REC = new ArrayList<Float>(); 
       
        ebmce.getCCDB(10);  ebmce.setGeom("2.5"); ebmce.isMC = true;
	
    	GEN.clear();  REC.clear();
    	
		boolean goodev = ebmce.readMC(de) && ebmce.pmc.size()==1;   
		
	    if (goodev) { 
		
        	GEN = getkin(ebmce.pmc); float refP = (float) GEN.get(0); float refTH = (float) GEN.get(1); float refE = refP;       	    	

	    	dropBanks(de);
	    	
        	if(!ebmce.processDataEvent(de)) return;  
        	
        	ebmce.processDataEvent(de);  
        	ebmce.getUnmatchedResponses();
        
        	List<DetectorResponse> rPC = ebmce.getPCResponses(trigFD);
        
        	double ecalE = ebmce.getEcalEnergy(trigFD);

//      Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.04 && energy>0.12;
//      Boolean trig2 = good_pcal && !good_ecal && part.epc>0.04;
//      Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.06 && energy>0.15;
//      Boolean trig2 = good_pcal && !good_ecal && part.epc>0.15;
        
//        	System.out.println("energy,1,2,3 = "+ecalE+" "+ebmce.epc+" "+ebmce.eec1+" "+ebmce.eec2);
        
        	Boolean good_pcal = ebmce.epc>0.00;                 //VTP reported cluster
        	Boolean good_ecal = (ebmce.eec1+ebmce.eec2)>0.001 ; //VTP reported cluster
        	Boolean trig1 = good_pcal &&  good_ecal && ebmce.epc>0.04 && ecalE>0.12;
        	Boolean trig2 = good_pcal && !good_ecal && ebmce.epc>0.12;
        		
//        trig1=true; trig2=true;
        
        	if(rPC.size()>0)                   {if(refTH>15)dgm.fill("effmom1",refE); dgm.fill("effthe1",refTH);}
        	if(rPC.size()==1 || rPC.size()==2) {if(refTH>15)dgm.fill("effmom2",refE); dgm.fill("effthe2",refTH);}
        	if(rPC.size()==1)                  {if(refTH>15)dgm.fill("effmom3",refE); dgm.fill("effthe3",refTH);
        	dgm.fill("h50",refP, ecalE/refP); 
        	dgm.fill("h53",ecalE,ecalE/refP);}  
        	
	    }
              
    }

    public void fith52() {
    	
		ParallelSliceFitter fitter = new ParallelSliceFitter(dgm.getH2F("h53"));
		fitter.setMin(0.2);
        fitter.fitSlicesX();
        
        GraphErrors sigGraph = fitter.getSigmaSlices();   
        sigGraph.setTitleX("ECAL Electron Energy (GeV)");  sigGraph.setTitleY("SIGMA E/P"); 	
        sigGraph.setMarkerSize(4); sigGraph.setMarkerStyle(1);
        
        dgm.getGE("g53").copy(fitter.getMeanSlices());  
        dgm.getGE("g53").setTitleX("ECAL Electron Energy (GeV)");  dgm.getGE("g53").setTitleY("Sampling Fraction E/P"); 	
        dgm.getGE("g53").setMarkerSize(4);  dgm.getGE("g53").setMarkerStyle(1); 
          
        int npts =  dgm.getGE("g53").getDataSize(0)  ;
        double[] xm  = new double[npts];
        double[] ym  = new double[npts];
        double[] yme = new double[npts];
        double[] xs  = new double[npts];
        double[] ys  = new double[npts];
        double[] yse = new double[npts];
        
        for (int i=0; i<npts; i++) {
      	   xm[i] = 1/Math.sqrt(dgm.getGE("g53").getDataX(i));
      	   ym[i] = dgm.getGE("g53").getDataY(i);
      	  yme[i] = dgm.getGE("g53").getDataEY(i);     
        } 
        
        GraphErrors resGraph = new GraphErrors();
        resGraph.setTitleX("1/E (GeV^-^1)");  resGraph.setTitleY("(#sigma / E)^2"); 
        resGraph.setMarkerSize(4); resGraph.setMarkerStyle(1);
        
        int n=0;
        for (int i=0; i<sigGraph.getDataSize(0); i++) {
            double y = sigGraph.getDataY(i); double ye = sigGraph.getDataEY(i);
            if(ym[i]>0&&y>0) {
                xs[n] = 1/sigGraph.getDataX(i)/ym[i];  //sig(E)/E vs True Energy
        	    ys[n] = y*y/ym[i]/ym[i]; //sigma(E)/E = sigma(E/P)*(P/E)
        	   yse[n] = 2*ys[n]*Math.sqrt(Math.pow(ye/y,2)+Math.pow(yme[i]/ym[i],2));
//        	   System.out.println(xs[n]+" "+ys[n]+" "+yse[n]);
                xs[n] = 1/Math.sqrt(sigGraph.getDataX(i)/ym[i]);  //sig(E)/E vs True Energy
       	        ys[n] = y*y/ym[i]/ym[i]; //sigma(E)/E = sigma(E/P)*(P/E)
       	       yse[n] = 2*ys[n]*Math.sqrt(Math.pow(ye/y,2)+Math.pow(yme[i]/ym[i],2));
//        	   System.out.println(xs[n]+" "+ys[n]+" "+yse[n]);
        	   resGraph.addPoint(xs[n], ys[n], 0., yse[n]);
        	}
        } 
    }
    
    public void geteffGENREC() {
		dgm.geteff("h60a","effmom1","effmom0");
		dgm.geteff("h60b","effmom2","effmom0");
		dgm.geteff("h60c","effmom3","effmom0");
		dgm.geteff("h62a","effmom4","effmom0");
		dgm.geteff("h62b","effmom5","effmom0");
		dgm.geteff("h62c","effmom6","effmom0");
		dgm.geteff("h62d","effmom7","effmom0");
		dgm.geteff("h61a","effthe1","effthe0");
		dgm.geteff("h61b","effthe2","effthe0");
		dgm.geteff("h61c","effthe3","effthe0");
		dgm.geteff("h63a","effthe4","effthe0");
		dgm.geteff("h63b","effthe5","effthe0");
		dgm.geteff("h63c","effthe6","effthe0");
		dgm.geteff("h70","eff1a","eff1");
		dgm.geteff("h71","eff2a","eff1");
		dgm.geteff("h72","eff3a","eff1");    
    }
    
    public void geteffGEMC() {
    	
    }    
    
    public GraphErrors getEff(int index, int i2, int icol, int ind) {
		DataGroup dg = this.getDataGroup().getItem(0,1,index,getRunNumber());    	
    	GraphErrors geff = new GraphErrors(); geff.setMarkerColor(icol); geff.setLineColor(icol);
    	geff.getAttributes().setTitleX("p_m_m (GeV)");
    	geff.getAttributes().setTitleY("Efficiency");
    	GraphErrors g1 = ((H1F)dg.getData(ind).get(0)).getGraph();
    	GraphErrors g2 = ((H1F)dg.getData(ind).get(i2)).getGraph();
    	for (int ix=0; ix<g1.getDataSize(0); ix++) {
    		float x1 = (float) g1.getDataX(ix); float y1 = (float) g1.getDataY(ix);
    		if(y1>0) {
        		float   y2 = (float)g2.getDataY(ix);
    			float  eff = y2/y1;
//    			float effe = (float) Math.sqrt((y1-1)*eff*(1-eff))/(y1-1); 
    			float effe = (float) Math.sqrt(eff*(1-eff)/(y1-1)); 
    			geff.addPoint(x1, eff, 0f, effe);
    		}
    	}
    	return geff;
    } 
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
