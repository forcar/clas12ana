package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

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

public class ECmc extends DetectorMonitor {
	
	Event       ev = new Event();
	EBMCEngine  eb = new EBMCEngine();
	
	DataGroup dg = null;
	   
    float dp1=0.1f, dp2=0.1f; int mcsec=1;
    
    List<Particle> elec = null;
    List<Particle> phot = null;
    List<Particle>  pim = null;
		
    public ECmc(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("GENREC","KINEMATICS","EFFICIENCY");
        
        this.use123Buttons(true);
        this.useSliderPane(true);

        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECmc.localinit()");
        
        engine.init();
        engine.isMC = true;
        engine.setVariation("default");        
        engine.setCalRun(10);  
        
        eb.getCCDB(10);
        eb.setGeom("2.5");
        eb.isMC = true;
        
        tl.setFitData(Fits);
    }
    
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
    	createGENREC(0);
    	createKINEMATICS(0);
       	createEFFICIENCY(0);
       	createEFFICIENCY(1);
       	createEFFICIENCY(2);
    }
    
    public void createGENREC(int st) {

    	switch (st) {        
        case 0: 
        	dgm.add("GENREC",4,5,0,st,getRunNumber());
        	dgm.makeH1("h00",50,-0.1,0.1,-1,"","#Deltapx GEN-REC (GeV)",1,4); 
        	dgm.makeH1("h01",50,-0.1,0.1,-1,"","#Deltapy GEN-REC (GeV)",1,4); 
        	dgm.makeH1("h02",50,-0.1,0.1,-1,"","#Deltapz GEN-REC (GeV)",1,4); 
        	dgm.makeH1("h03",50,-0.1,0.1,-1,"","#Deltap  GEN-REC (GeV)",1,4);
        
        	dgm.makeH2("hv0",50,0,4.0, 50,-5,5,    0,"ECAL E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVx  GEN-REC (cm)");
        	dgm.makeH2("hv1",50,0,4.0, 50,-5,5,    0,"ECAL E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVy  GEN-REC (cm)");
        	dgm.makeH2("hv2",50,0,4.0, 50,-5,5,    0,"ECAL E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVz  GEN-REC (cm)");
        	dgm.makeH2("hv3",50,-6,0,  60,-10,30,  0,"#Deltap>"+dp1,"Vz GEN (cm)","z0 (cm)");
        	dgm.makeH2("hv4",50,5,35,  60,-10,30,  0,"#Deltap>"+dp1,"True Electron Theta (deg)","z0 (cm)");
        	dgm.makeH2("hv5",50,0, 9,  60,-10,30,  0,"#Deltap>"+dp1,"True Electron Momentum (GeV)","z0 (cm)");
        	dgm.makeH2("h10",50,0,4.0, 50,-0.1,1.0,0,"","ECAL #gamma E (GeV)","#Deltapx (GeV)");
        	dgm.makeH2("h11",50,0,1.0, 50,-0.1,0.1,0,"","ECAL #gamma E (GeV)","#Deltapy (GeV)");
        	dgm.makeH2("h12",50,0,1.0, 50,-0.1,1.0,0,"","ECAL #gamma E (GeV)","#Deltapz (GeV)");
        	dgm.makeH2("h13",50,0,5.0, 50,-0.1,5.0,0,"","ECAL #gamma E (GeV)","#Deltap  (GeV)");
        	dgm.makeH2("h20",50,0,100,100,-0.1,5.0,0,"","Distance e-#gamma (cm)","#Deltap  GEN-REC (GeV)");   
        	dgm.makeH2("h21",50,0,100,100, 0.0,5.0,0,"#Deltap>"+dp1,"Distance e-#gamma (cm)","#Deltap  GEN-REC (GeV)");   
        	dgm.makeH2("h22",50,0,100,100, 0.0,5.0,0,"#Deltap<"+dp2,"Distance e-#gamma (cm)","#Deltap  GEN-REC (GeV)");   
        	dgm.makeH2("h23",50,0,  9,100, 0.0,5.0,0,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","ECAL #gamma E (GeV)");   
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
        dgm.makeGraph("g53",-2,1); SFFunction sf = new SFFunction("esf",-11,eb.ccdb,0.1,2.5); dgm.addDataSet(sf,-2);
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
        dgm.makeH1("h73a",10,0.5,10.5,-1,"","Multiplicity",1,1,2);
        dgm.makeH1("h73b",10,0.5,10.5,-2,"","Multiplicity",1,34,2);
        dgm.makeH1("h73c",10,0.5,10.5,-2,"","Multiplicity",1,2,2);
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
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;
    	plotMCHistos();
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
    	fith52(); geteff();
    	isAnalyzeDone = true;
    }
    
    @Override
    public void processEvent(DataEvent event) {
    	
    	ev.init(event);        	
    	ev.setEventNumber(10);
    	ev.requireOneElectron(false);
   	    ev.setElecTriggerSector(mcsec);
   	    ev.setMCpid(11);
   	    
   	    float refE=0,refP=0,refTH=0,refPH=0;
   	    		
        if(ev.procEvent(event)) {
        	            
            refE  = (float) ev.pmc.get(0).e(); 
            refP  = (float) ev.pmc.get(0).p(); 
            refTH = (float) Math.toDegrees(ev.pmc.get(0).theta());
            refPH = (float) Math.toDegrees(ev.pmc.get(0).phi());
            
            dgm.fill("eff1",refE, refTH);
            if(refTH>15)dgm.fill("effmom0",refE);
            dgm.fill("effthe0",refTH);  
            
    		elec = ev.getPART(0.0, 11);      int nelec = elec.size(); dgm.fill("h73a",nelec);
    		phot = ev.getPART(0.0, 22);      int nphot = phot.size(); dgm.fill("h73b",nphot);
    		pim  = ev.getPART(0.0, 212);     int  npim = pim.size();  dgm.fill("h73c",npim);
    		
    		int npart = nelec+nphot+npim;        		        		
    		float esum=0,dp=0,chi2pid=1000;
    		
    		if(npart>0) {if(refTH>15)dgm.fill("effmom4",refE);dgm.fill("effthe4",refTH); dgm.fill("eff1a",refE, refTH);}
    		if(npim==0) {if(refTH>15)dgm.fill("effmom5",refE);dgm.fill("effthe5",refTH);}
   		
    		for (Particle p : phot) {esum+=p.e();}
    		        		
    		if(nelec==1) {
        		chi2pid = (float) Math.abs(elec.get(0).getProperty("chi2pid"));
        		dp = (float)(ev.pmc.get(0).p()-elec.get(0).p());
        		
        		if(npim==0&&chi2pid<3.5)  {
        			if(refTH>15)                     dgm.fill("effmom6",refE); 
        			if(refTH>15 && Math.abs(dp)<dp2) dgm.fill("effmom7",refE); 
        			if(Math.abs(dp)<dp2)  dgm.fill("eff3a",refE,refTH);
        			dgm.fill("effthe6",refTH); dgm.fill("eff2a",refE,refTH); 	                    
        		}
        		
        		if(esum<0.01) {
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
       			float ex=ev.getVar(elecECAL, "x", 1);float ey=ev.getVar(elecECAL, "y", 1);float ez=ev.getVar(elecECAL, "z", 1);
       			
       			for (Particle rp : phot) { //loop over photons from REC::Particle
         			List<Particle> photECAL = ev.getECAL((int)rp.getProperty("pindex"));
           			for (Particle rc : photECAL) { //loop over photon entries in REC::Calorimeter
           			  int il = (int)rc.getProperty("layer");
         		      float ec_the = (float) Math.toDegrees(rc.theta());
           		      float ec_phi = (float) Math.toDegrees(rc.phi()); 
           			  float dthe = e_the-ec_the; 
           			  float dphi = e_phi-ec_phi; 
           			  float dphi_ref = refPH-ec_phi;
           			  
                      if(il==1) {  //use only PCAL and ECIN entries                        
//                    double z0     =    pt*0.667*Math.sin(dphi*3.14159/180.)*100; 
                      double z0  =    ptref*0.667*Math.sin(dphi_ref*3.14159/180.)*100; 
               	      if(Math.abs(dp)>dp1) {dgm.fill("hv3",ev.pmc.get(0).vz(),z0);dgm.fill("hv4",refTH,z0);dgm.fill("hv5",refP,z0);}	
           			  float gx=(float)rc.getProperty("x");float gy=(float)rc.getProperty("y");float gz=(float)rc.getProperty("z");
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
//    		System.out.println(" "+nevent+" "+nev+" "+nelec+" "+nphot+" "+npim);
        }
        
        dropBanks(event);
        engine.processDataEvent(event);  
        eb.processDataEvent(event);  
        eb.getUnmatchedResponses();
        List<DetectorResponse> rPC = eb.getPCResponses(mcsec);
        
        double ecalE = eb.getEcalEnergy(mcsec);

//      Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.04 && energy>0.12;
//      Boolean trig2 = good_pcal && !good_ecal && part.epc>0.04;
//      Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.06 && energy>0.15;
//      Boolean trig2 = good_pcal && !good_ecal && part.epc>0.15;
        
//        System.out.println("energy,1,2,3 = "+energy+" "+part.epc+" "+part.eec1+" "+part.eec2);
        
        Boolean good_pcal = eb.epc>0.00;              //VTP reported cluster
        Boolean good_ecal = (eb.eec1+eb.eec2)>0.001 ; //VTP reported cluster
        Boolean trig1 = good_pcal &&  good_ecal && eb.epc>0.04 && ecalE>0.12;
        Boolean trig2 = good_pcal && !good_ecal && eb.epc>0.12;
        		
//        trig1=true; trig2=true;
        
        if(rPC.size()>0)                   {if(refTH>15)dgm.fill("effmom1",refE); dgm.fill("effthe1",refTH);}
        if(rPC.size()==1 || rPC.size()==2) {if(refTH>15)dgm.fill("effmom2",refE); dgm.fill("effthe2",refTH);}
        if(rPC.size()==1)                  {if(refTH>15)dgm.fill("effmom3",refE); dgm.fill("effthe3",refTH);
        dgm.fill("h50",refP, ecalE/refP); 
        dgm.fill("h53",ecalE,ecalE/refP);}                       
              
    }
    
    public void plotMCHistos() {      
        plot("GENREC");
        plot("KINEMATICS");
        plot("EFFICIENCY");   
    }
    
    @Override
    public void plotEvent(DataEvent de) {
//        analyze();         
    }
    
    public void fith52() {
    	
		ParallelSliceFitter fitter = new ParallelSliceFitter(dgm.getH2F("h53"));
		fitter.setMin(0.2);
        fitter.fitSlicesX();
        
        GraphErrors sigGraph = fitter.getSigmaSlices();   
        sigGraph.setTitleX("ECAL Electron Energy (GeV)");  sigGraph.setTitleY("SIGMA E/P"); 	
        sigGraph.setMarkerSize(4); sigGraph.setMarkerStyle(1);
        
        dgm.getGraph("g53").copy(fitter.getMeanSlices());  
        dgm.getGraph("g53").setTitleX("ECAL Electron Energy (GeV)");  dgm.getGraph("g53").setTitleY("Sampling Fraction E/P"); 	
        dgm.getGraph("g53").setMarkerSize(4);  dgm.getGraph("g53").setMarkerStyle(1); 
          
        int npts =  dgm.getGraph("g53").getDataSize(0)  ;
        double[] xm  = new double[npts];
        double[] ym  = new double[npts];
        double[] yme = new double[npts];
        double[] xs  = new double[npts];
        double[] ys  = new double[npts];
        double[] yse = new double[npts];
        
        for (int i=0; i<npts; i++) {
      	   xm[i] = 1/Math.sqrt(dgm.getGraph("g53").getDataX(i));
      	   ym[i] = dgm.getGraph("g53").getDataY(i);
      	  yme[i] = dgm.getGraph("g53").getDataEY(i);     
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
    
    public void geteff() {
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
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
