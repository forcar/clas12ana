package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.analysis.ECPart.SFFunction;
import org.clas.tools.EBMC;
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
	
	Event ev = new Event();
	EBMC  eb = new EBMC();
	DataGroup dg = null;
	   
    float dp1=0.1f, dp2=0.1f; int mcsec=1;
    
    List<Particle> elec = null;
    List<Particle> phot = null;
    List<Particle>  pim = null;
		
    public ECmc(String name) {
        super(name);
        setDetectorTabNames("GENREC",
        		            "KINEMATICS",
        		            "EFFICIENCY");
        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
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
        eb.setThresholds("Pizero",engine);
        eb.setGeom("2.5");
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
    	slider.setValue(0);
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
        
    	String tab = "GENREC", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0;
        tag = st+"_"+k+"_"+run;
 
    	switch (st) {        
        case 0: 
        dg = new DataGroup(4,5);
        dg.addDataSet(makeH1("h00",50,-0.1,0.1,"","#Deltapx GEN-REC (GeV)",1,4),n++); 
        dg.addDataSet(makeH1("h01",50,-0.1,0.1,"","#Deltapy GEN-REC (GeV)",1,4),n++); 
        dg.addDataSet(makeH1("h02",50,-0.1,0.1,"","#Deltapz GEN-REC (GeV)",1,4),n++); 
        dg.addDataSet(makeH1("h03",50,-0.1,0.1,"","#Deltap  GEN-REC (GeV)",1,4),n++);
        dg.addDataSet(makeH2("hv0",50,0,4.0, 50,-5,5,"ECAL E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVx  GEN-REC (cm)"),n++);
        dg.addDataSet(makeH2("hv1",50,0,4.0, 50,-5,5,"ECAL E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVy  GEN-REC (cm)"),n++);
        dg.addDataSet(makeH2("hv2",50,0,4.0, 50,-5,5,"ECAL E#gamma>10 MeV","ECAL #gamma E (GeV)","#DeltaVz  GEN-REC (cm)"),n++);
        dg.addDataSet(makeH2("hv3",50,-6,0,60,-10,30,"#Deltap>"+dp1,"Vz GEN (cm)","z0 (cm)"),n++);
        dg.addDataSet(makeH2("hv4",50,5,35,60,-10,30,"#Deltap>"+dp1,"True Electron Theta (deg)","z0 (cm)"),n++);
        dg.addDataSet(makeH2("hv5",50,0, 9,60,-10,30,"#Deltap>"+dp1,"True Electron Momentum (GeV)","z0 (cm)"),n++);
        dg.addDataSet(makeH2("h10",50,0,4.0, 50,-0.1,1.0,"","ECAL #gamma E (GeV)","#Deltapx (GeV)"),n++);
        dg.addDataSet(makeH2("h11",50,0,1.0, 50,-0.1,0.1,"","ECAL #gamma E (GeV)","#Deltapy (GeV)"),n++);
        dg.addDataSet(makeH2("h12",50,0,1.0, 50,-0.1,1.0,"","ECAL #gamma E (GeV)","#Deltapz (GeV)"),n++);
        dg.addDataSet(makeH2("h13",50,0,5.0, 50,-0.1,5.0,"","ECAL #gamma E (GeV)","#Deltap  (GeV)"),n++);
        dg.addDataSet(makeH2("h20",50,0,100,100,-0.1,5.0,"","Distance e-#gamma (cm)","#Deltap  GEN-REC (GeV)"),n++);   
        dg.addDataSet(makeH2("h21",50,0,100,100, 0.0,5.0,"#Deltap>"+dp1,"Distance e-#gamma (cm)","#Deltap  GEN-REC (GeV)"),n++);   
        dg.addDataSet(makeH2("h22",50,0,100,100, 0.0,5.0,"#Deltap<"+dp2,"Distance e-#gamma (cm)","#Deltap  GEN-REC (GeV)"),n++);   
        dg.addDataSet(makeH2("h23",50,0,  9,100, 0.0,5.0,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","ECAL #gamma E (GeV)"),n++);   
    	}
        this.getDataGroup().add(dg,0,st,k,run);        
    }
    
    public void createKINEMATICS(int st) {
    	
    	String tab = "KINEMATICS", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0;
        tag = st+"_"+k+"_"+run;
      
    	switch (st) {        
        case 0: 
        dg = new DataGroup(4,2);      
        String lab[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};
        dg.addDataSet(makeH2("h30",50,0,9,70,-2,12,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","PCAL "+lab[0]),n++);   
        dg.addDataSet(makeH2("h31",50,0,9,70,-5,30,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","PCAL "+lab[1]),n++);   
        dg.addDataSet(makeH2("h32",50,0,9,70,-2,12,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","ECIN "+lab[0]),n++);   
        dg.addDataSet(makeH2("h33",50,0,9,70,-5,30,"#Deltap>"+dp1,"Track Electron Momentum (GeV)","ECIN "+lab[1]),n++);
        dg.addDataSet(makeH2("h40",50,0,9,70,-2,12,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","PCAL "+lab[0]),n++);   
        dg.addDataSet(makeH2("h41",50,0,9,70,-5,30,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","PCAL "+lab[1]),n++);   
        dg.addDataSet(makeH2("h42",50,0,9,70,-2,12,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","ECIN "+lab[0]),n++);   
        dg.addDataSet(makeH2("h43",50,0,9,70,-5,30,"#Deltap<"+dp2,"Track Electron Momentum (GeV)","ECIN "+lab[1]),n++);
    	}
        this.getDataGroup().add(dg,0,st,k,run);        
    }
    
    public void createEFFICIENCY(int st) {
    	
    	String tab = "EFFICIENCY", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0, nh=8;;
    	
    	switch (st) {        
        case 0: 
        dg = new DataGroup(4,3);         
        dg.addDataSet(makeH2("h50",50,0.0,10.,60,0.0,0.31,"","True Electron Momentum (GeV)","Sampling Fraction E/P"),n++);   
        dg.addDataSet(makeH2("h51",50,0.0,10.,60,0.0,0.31,"","Track Electron Momentum (GeV)","Sampling Fraction E/P"),n++);   
        dg.addDataSet(makeH2("h52",50,0.0,10.,60,0.0,0.31,"","Track Electron Momentum (GeV)","Sampling Fraction E/P"),n++);   
        dg.addDataSet(makeH2("h53",50,0.0,2.5,60,0.0,0.31,"","ECAL Electron Energy (GeV)","Sampling Fraction E/P"),n++);
        dg.addDataSet(makeGraph("g53",1),n-1); SFFunction sf = new SFFunction("esf",-11,eb.ccdb,0.1,2.5); dg.addDataSet(sf,n-1);
        dg.addDataSet(makeH1("h60a",50,0,10,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Energy (GeV)","Efficiency #theta>15",1,3),n++);
        dg.addDataSet(makeH1("h60b",50,0,10,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Energy (GeV)","Efficiency #theta>15",1,5),n-1);
        dg.addDataSet(makeH1("h60c",50,0,10,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Energy (GeV)","Efficiency #theta>15",1,2),n-1);
        dg.addDataSet(makeH1("h61a",50,5,35,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Theta (deg)","Efficiency",1,3),n++);
        dg.addDataSet(makeH1("h61b",50,5,35,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Theta (deg)","Efficiency",1,5),n-1);
        dg.addDataSet(makeH1("h61c",50,5,35,"G: PC > 0  Y: PC = 1 or 2  R: PC = 1","True Electron Theta (deg)","Efficiency",1,2),n-1);
        dg.addDataSet(makeH1("h62a",50,0,10,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,3),n++);
        dg.addDataSet(makeH1("h62b",50,0,10,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,5),n-1);
        dg.addDataSet(makeH1("h62c",50,0,10,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,2),n-1);
        dg.addDataSet(makeH1("h62d",50,0,10,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Energy (GeV)","Efficiency #theta>15",1,1),n-1);
        dg.addDataSet(makeH1("h63a",50,5,35,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Theta (deg)","Efficiency",1,3),n++);
        dg.addDataSet(makeH1("h63b",50,5,35,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Theta (deg)","Efficiency",1,5),n-1);
        dg.addDataSet(makeH1("h63c",50,5,35,"G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2,"True Electron Theta (deg)","Efficiency",1,2),n-1);
    	dg.addDataSet(makeH2("h70",30,0,10,30,5,30,"e-,#gamma,#pi-","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
    	dg.addDataSet(makeH2("h71",30,0,10,30,5,30,"#chiPID<3.5","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
    	dg.addDataSet(makeH2("h72",30,0,10,30,5,30,"#DeltaP<"+dp2,"True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
    	dg.addDataSet(makeH1("h73a",10,1,11,"","Multiplicity",1,1), n++);
    	dg.addDataSet(makeH1("h73b",10,1,11,"","Multiplicity",1,34), n-1);
    	dg.addDataSet(makeH1("h73c",10,1,11,"","Multiplicity",1,2), n-1);
        break;
        case 1:
        dg = new DataGroup(4,3);
        dg.addDataSet(makeH2("eff1", 30,0,10,30,5,30,"","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
        dg.addDataSet(makeH2("eff1a",30,0,10,30,5,30,"","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
        dg.addDataSet(makeH2("eff1b",30,0,10,30,5,30,"","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
        dg.addDataSet(makeH2("eff2a",30,0,10,30,5,30,"","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
        dg.addDataSet(makeH2("eff2b",30,0,10,30,5,30,"","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
        dg.addDataSet(makeH2("eff3a",30,0,10,30,5,30,"","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
        dg.addDataSet(makeH2("eff3b",30,0,10,30,5,30,"","True Electron Momentum (GeV)","True Electron Theta (deg)"),n++);
        break;
        case 2:
        dg = new DataGroup(6,6);
        for (int i=0;i<nh;i++) dg.addDataSet(makeH1("effmom"+i,50,0.,10.,"","MC Electron P (GeV)",4,4),n++);
        for (int i=0;i<nh;i++) dg.addDataSet(makeH1("effthe"+i,50,5.,35.,"","MC Electron #theta (deg)",4,4),n++);
    	}
    	
        this.getDataGroup().add(dg,0,st,k,run);
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
    
    
    public List<Particle> getPART(double thr, int pid) {   	
    	List<Particle> olist = new ArrayList<Particle>();    
    	for (Particle p : ev.getParticle(pid)) {
    		short status = (short) p.getProperty("status");
    		if(status>=2000 && status<3000 && p.p()>thr) olist.add(p); 
    	}          	
       return olist;    	
    }
    
    @Override
    public void processEvent(DataEvent event) {
    	
		int run = getRunNumber();
		DataGroup  dg0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("GENREC"),run);
		DataGroup  dg1 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("KINEMATICS"),run);
		DataGroup dg20 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg21 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg22 = this.getDataGroup().getItem(0,2,getDetectorTabNames().indexOf("EFFICIENCY"),run);

    	if(eb.readMC(event)) { 
            
        float refE  = (float) eb.pmc.get(0).e(); 
        float refP  = (float) eb.pmc.get(0).p(); 
        float refTH = (float) Math.toDegrees(eb.pmc.get(0).theta());
        float refPH = (float) Math.toDegrees(eb.pmc.get(0).phi());
        
        dg21.getH2F("eff1").fill(refE, refTH);
        if(refTH>15)dg22.getH1F("effmom0").fill(refE);
        dg22.getH1F("effthe0").fill(refTH);    
        
    	ev.init(event);        	
    	ev.setEventNumber(10);
    	ev.requireOneElectron(false);
   	    ev.setElecTriggerSector(mcsec);
   	    		
        if(ev.procEvent(event)) {
    		elec = getPART(0.0, 11);      int nelec = elec.size(); dg20.getH1F("h73a").fill(nelec);
    		phot = getPART(0.0, 22);      int nphot = phot.size(); dg20.getH1F("h73b").fill(nphot);
    		pim  = getPART(0.0, 212);     int  npim = pim.size();  dg20.getH1F("h73c").fill(npim);
    		
    		int npart = nelec+nphot+npim;        		        		
    		float esum=0,dp=0,chi2pid=1000;
    		
    		if(npart>0) {if(refTH>15)dg22.getH1F("effmom4").fill(refE); dg22.getH1F("effthe4").fill(refTH); dg21.getH2F("eff1a").fill(refE, refTH);}
    		if(npim==0) {if(refTH>15)dg22.getH1F("effmom5").fill(refE); dg22.getH1F("effthe5").fill(refTH);}
   		
    		for (Particle p : phot) {esum+=p.e();}
    		        		
    		if(nelec==1) {
        		chi2pid = (float) Math.abs(elec.get(0).getProperty("chi2pid"));
        		dp = (float)(eb.pmc.get(0).p()-elec.get(0).p());
        		
        		if(npim==0&&chi2pid<3.5)  {
        			if(refTH>15)                      dg22.getH1F("effmom6").fill(refE); 
        			if(refTH>15 && Math.abs(dp)<dp2)  dg22.getH1F("effmom7").fill(refE); 
        			if(Math.abs(dp)<dp2)   dg21.getH2F("eff3a").fill(refE,refTH);
        			dg22.getH1F("effthe6").fill(refTH); dg21.getH2F("eff2a").fill(refE,refTH); 	                    
        		}
        		
        		if(esum<0.01) {
        			dg0.getH1F("h00").fill(eb.pmc.get(0).px()-elec.get(0).px());
        			dg0.getH1F("h01").fill(eb.pmc.get(0).py()-elec.get(0).py());
        			dg0.getH1F("h02").fill(eb.pmc.get(0).pz()-elec.get(0).pz());
        			dg0.getH1F("h03").fill(dp);
        		}
    			
                if(esum>0.01) {
        			dg0.getH2F("h10").fill(esum,eb.pmc.get(0).px()-elec.get(0).px());
        			dg0.getH2F("h11").fill(esum,eb.pmc.get(0).py()-elec.get(0).py());
        			dg0.getH2F("h12").fill(esum,eb.pmc.get(0).pz()-elec.get(0).pz());
        			dg0.getH2F("h13").fill(esum,dp);
        			dg0.getH2F("hv0").fill(esum,eb.pmc.get(0).vx()-elec.get(0).vx());
        			dg0.getH2F("hv1").fill(esum,eb.pmc.get(0).vy()-elec.get(0).vy());
        			dg0.getH2F("hv2").fill(esum,eb.pmc.get(0).vz()-elec.get(0).vz());
                }
                              		
          		float e_mom = (float) elec.get(0).p();
                float e_the = (float) Math.toDegrees(elec.get(0).theta());
                float e_phi = (float) Math.toDegrees(elec.get(0).phi());
                
                double    pt = Math.sqrt(e_mom*e_mom-elec.get(0).pz()*elec.get(0).pz());
                double ptref = Math.sqrt(refP*refP-eb.pmc.get(0).pz()*eb.pmc.get(0).pz());
                
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
               	      if(Math.abs(dp)>dp1) {dg0.getH2F("hv3").fill(eb.pmc.get(0).vz(),z0);dg0.getH2F("hv4").fill(refTH,z0);dg0.getH2F("hv5").fill(refP,z0);}	
           			  float gx=(float)rc.getProperty("x");float gy=(float)rc.getProperty("y");float gz=(float)rc.getProperty("z");
           			  float dist = (float) Math.sqrt((ex-gx)*(ex-gx)+(ey-gy)*(ey-gy)+(ez-gz)*(ez-gz));
           			  dg0.getH2F("h20").fill(dist, dp); 
           			  if(Math.abs(dp)>dp1) {dg1.getH2F("h30").fill(e_mom,dthe); dg1.getH2F("h31").fill(e_mom,dphi); dg0.getH2F("h21").fill(dist,esum);}
           			  if(Math.abs(dp)<dp2) {dg1.getH2F("h40").fill(e_mom,dthe); dg1.getH2F("h41").fill(e_mom,dphi); dg0.getH2F("h22").fill(dist,esum);dg0.getH2F("h23").fill(e_mom,esum);}
                      }
                      if(il==4 && Math.abs(dp)>dp1) {dg1.getH2F("h32").fill(e_mom,dthe); dg1.getH2F("h33").fill(e_mom,dphi);}
                      if(il==4 && Math.abs(dp)<dp2) {dg1.getH2F("h42").fill(e_mom,dthe); dg1.getH2F("h43").fill(e_mom,dphi);}
           			}
           		}
           		float en=0; for (Particle p : elecECAL) en += (float) p.getProperty("energy");    
           		dg20.getH2F("h51").fill(e_mom, en*1e-3/e_mom);
           		if(chi2pid<3.5) dg20.getH2F("h52").fill(e_mom, en*1e-3/e_mom);
    		}
//    		System.out.println(" "+nevent+" "+nev+" "+nelec+" "+nphot+" "+npim);
        }
        
        dropBanks(event);
        engine.processDataEvent(event);  
        eb.readEC(event,"ECAL::clusters");  
        eb.getUnmatchedResponses();
        List<DetectorResponse> rPC = eb.getPCResponses(mcsec);
        
        double ecalE = eb.getEcalEnergy(mcsec);

//      Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.04 && energy>0.12;
//      Boolean trig2 = good_pcal && !good_ecal && part.epc>0.04;
//      Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.06 && energy>0.15;
//      Boolean trig2 = good_pcal && !good_ecal && part.epc>0.15;
        
//        System.out.println("energy,1,2,3 = "+energy+" "+part.epc+" "+part.eec1+" "+part.eec2);
        
        Boolean good_pcal = eb.epc>0.00;               //VTP reported cluster
        Boolean good_ecal = (eb.eec1+eb.eec2)>0.001 ; //VTP reported cluster
        Boolean trig1 = good_pcal &&  good_ecal && eb.epc>0.04 && ecalE>0.12;
        Boolean trig2 = good_pcal && !good_ecal && eb.epc>0.12;
        		
//        trig1=true; trig2=true;
        
        if(rPC.size()>0)                   {if(refTH>15)dg22.getH1F("effmom1").fill(refE); dg22.getH1F("effthe1").fill(refTH);}
        if(rPC.size()==1 || rPC.size()==2) {if(refTH>15)dg22.getH1F("effmom2").fill(refE); dg22.getH1F("effthe2").fill(refTH);}
        if(rPC.size()==1)                  {if(refTH>15)dg22.getH1F("effmom3").fill(refE); dg22.getH1F("effthe3").fill(refTH);
        dg20.getH2F("h50").fill(refP, ecalE/refP); 
        dg20.getH2F("h53").fill(ecalE,ecalE/refP);}                       
    }          
    }
    
    public void plotMCHistos() {      
        plot("GENREC");
        plot("KINEMATICS");
        plot("EFFICIENCY");   
    }
    
    public void plotEFFICIENCY() {
    	
        EmbeddedCanvas c = null; int ncd;

        int index=getDetectorTabNames().indexOf("EFFICIENCY");
        c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); c.divide(4,3); ncd=3;           
		c.setAxisTitleSize(18);
		c.setAxisLabelSize(18);
		
		int run = getRunNumber();
		DataGroup dg2 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
          
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.05);  c.getPad(ncd).setTitleFontSize(18);
        c.draw(dg2.getH2F("h60a"));  c.draw(dg2.getH2F("h60b"),"same"); c.draw(dg2.getH2F("h60c"),"same");
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.05);  c.getPad(ncd).setTitleFontSize(18);
        c.draw(dg2.getH2F("h61a"));  c.draw(dg2.getH2F("h61b"),"same"); c.draw(dg2.getH2F("h61c"),"same");
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.015); c.getPad(ncd).setTitleFontSize(18);
                          c.getPad(ncd).getAxisX().setGrid(true);       c.getPad(ncd).getAxisY().setGrid(true);                               
        c.draw(dg2.getH2F("h62a"));  c.draw(dg2.getH2F("h62b"),"same"); c.draw(dg2.getH2F("h62c"),"same"); c.draw(dg2.getH2F("h62d"),"same"); 
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.015); c.getPad(ncd).setTitleFontSize(18);
        c.draw(dg2.getH2F("h63a"));  c.draw(dg2.getH2F("h63b"),"same"); c.draw(dg2.getH2F("h63c"),"same"); 
        
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setRange(0.9,1.01); c.getPad(ncd).setTitleFontSize(18); c.draw(dg2.getH2F("h70"));
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setRange(0.9,1.01); c.getPad(ncd).setTitleFontSize(18); c.draw(dg2.getH2F("h71"));
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setRange(0.5,1.01); c.getPad(ncd).setTitleFontSize(18); c.draw(dg2.getH2F("h72"));
    }
    
    @Override
    public void plotEvent(DataEvent de) {
//        analyze();         
    }
    
    public void fith52() {
    	
		int run = getRunNumber();
		DataGroup dg2 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
       
		ParallelSliceFitter fitter = new ParallelSliceFitter(dg2.getH2F("h53"));
		fitter.setMin(0.2);
        fitter.fitSlicesX();
        
        GraphErrors sigGraph = fitter.getSigmaSlices();   
        sigGraph.setTitleX("ECAL Electron Energy (GeV)");  sigGraph.setTitleY("SIGMA E/P"); 	
        sigGraph.setMarkerSize(4); sigGraph.setMarkerStyle(1);
        
        dg2.getGraph("g53").copy(fitter.getMeanSlices());  
        dg2.getGraph("g53").setTitleX("ECAL Electron Energy (GeV)");  dg2.getGraph("g53").setTitleY("Sampling Fraction E/P"); 	
        dg2.getGraph("g53").setMarkerSize(4);  dg2.getGraph("g53").setMarkerStyle(1); 
          
        int npts =  dg2.getGraph("g53").getDataSize(0)  ;
        double[] xm  = new double[npts];
        double[] ym  = new double[npts];
        double[] yme = new double[npts];
        double[] xs  = new double[npts];
        double[] ys  = new double[npts];
        double[] yse = new double[npts];
        
        for (int i=0; i< dg2.getGraph("g53").getDataSize(0); i++) {
      	   xm[i] = 1/Math.sqrt(dg2.getGraph("g53").getDataX(i));
      	   ym[i] = dg2.getGraph("g53").getDataY(i);
      	  yme[i] = dg2.getGraph("g53").getDataEY(i);     
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
		int run = getRunNumber();
		DataGroup dg20 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg21 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg22 = this.getDataGroup().getItem(0,2,getDetectorTabNames().indexOf("EFFICIENCY"),run);
        dg20.getH1F("h60a").add(H1F.divide(dg22.getH1F("effmom1"),dg22.getH1F("effmom0")));
        dg20.getH1F("h60b").add(H1F.divide(dg22.getH1F("effmom2"),dg22.getH1F("effmom0"))); 
        dg20.getH1F("h60c").add(H1F.divide(dg22.getH1F("effmom3"),dg22.getH1F("effmom0"))); 
        dg20.getH1F("h62a").add(H1F.divide(dg22.getH1F("effmom4"),dg22.getH1F("effmom0"))); 
        dg20.getH1F("h62b").add(H1F.divide(dg22.getH1F("effmom5"),dg22.getH1F("effmom0"))); 
        dg20.getH1F("h62c").add(H1F.divide(dg22.getH1F("effmom6"),dg22.getH1F("effmom0")));
        dg20.getH1F("h62d").add(H1F.divide(dg22.getH1F("effmom7"),dg22.getH1F("effmom0")));        
        dg20.getH1F("h61a").add(H1F.divide(dg22.getH1F("effthe1"),dg22.getH1F("effthe0"))); 
        dg20.getH1F("h61b").add(H1F.divide(dg22.getH1F("effthe2"),dg22.getH1F("effthe0"))); 
        dg20.getH1F("h61c").add(H1F.divide(dg22.getH1F("effthe3"),dg22.getH1F("effthe0"))); 
        dg20.getH1F("h63a").add(H1F.divide(dg22.getH1F("effthe4"),dg22.getH1F("effthe0"))); 
        dg20.getH1F("h63b").add(H1F.divide(dg22.getH1F("effthe5"),dg22.getH1F("effthe0"))); 
        dg20.getH1F("h63c").add(H1F.divide(dg22.getH1F("effthe6"),dg22.getH1F("effthe0")));                
        dg20.getH2F("h70").add(H2F.divide(dg21.getH2F("eff1a"),dg21.getH2F("eff1")));
        dg20.getH2F("h71").add(H2F.divide(dg21.getH2F("eff2a"),dg21.getH2F("eff1")));
        dg20.getH2F("h72").add(H2F.divide(dg21.getH2F("eff3a"),dg21.getH2F("eff1")));     
    }
    
    public void plot(String tabname) {      	
    	int index = getDetectorTabNames().indexOf(tabname);
    	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));   	
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
