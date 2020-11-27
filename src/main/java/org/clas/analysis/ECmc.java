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
	
    H1F h00,h01,h02,h03;
    H2F h10,h11,h12,h13,h20,h21,h22,h23,h30,h31,h32,h33,h40,h41,h42,h43,h50,h51,h52,eff1,eff1a,eff1b,eff2a,eff2b,eff3a,eff3b;
	H1F effmom1,effmom2,effmom3,effmom4,effmom5,effmom6,effmom7,effthe1,effthe2,effthe3,effthe4,effthe5,effthe6;
	GraphErrors meanGraph = null;
	
    int nh=10;
    H1F effmom[] = new H1F[nh]; H1F effthe[] = new H1F[nh];
    float dp1=0.1f, dp2=0.1f; int mcsec=1;
    
    List<Particle> elec = null;
    List<Particle> phot = null;
    List<Particle>  pim = null;
		
    public ECmc(String name) {
        super(name);
        setDetectorTabNames("GEN-REC",
        		            "DTHE,DPHI",
        		            "EFF");

        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECmc.localinit()");
        configEngine("muon");  
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
    	histosExist = true;    	
    	DataGroup dg = null;
        int n;
        
        dg = new DataGroup(4,3); n=0;
        h00 = new H1F("h00",50,-0.1,0.1);           h00.setTitleX("#Deltapx GEN-REC (GeV)"); dg.addDataSet(h00, n); n++;
        h01 = new H1F("h01",50,-0.1,0.1);           h01.setTitleX("#Deltapy GEN-REC (GeV)"); dg.addDataSet(h01, n); n++;
        h02 = new H1F("h02",50,-0.1,0.1);           h02.setTitleX("#Deltapz GEN-REC (GeV)"); dg.addDataSet(h02, n); n++;
        h03 = new H1F("h03",50,-0.1,0.1);           h03.setTitleX("#Deltap  GEN-REC (GeV)"); dg.addDataSet(h03, n); n++;       
        h10 = new H2F("h10",50,0,4.0,50,-0.1,1.0);  h10.setTitleX("ECAL #gamma E (GeV)"); h10.setTitleY("#Deltapx  (GeV)"); dg.addDataSet(h10, n); n++;
        h11 = new H2F("h11",50,0,1.0,50,-0.1,0.1);  h11.setTitleX("ECAL #gamma E (GeV)"); h11.setTitleY("#Deltapy  (GeV)"); dg.addDataSet(h11, n); n++;
        h12 = new H2F("h12",50,0,1.0,50,-0.1,1.0);  h12.setTitleX("ECAL #gamma E (GeV)"); h12.setTitleY("#Deltapz  (GeV)"); dg.addDataSet(h12, n); n++;
        h13 = new H2F("h13",50,0,5.0,100,-0.1,5.0); h13.setTitleX("ECAL #gamma E (GeV)"); h13.setTitleY("#Deltap   (GeV)"); dg.addDataSet(h13, n); n++;        
        h20 = new H2F("h20",50,0,100.0,100,-1,5);   h20.setTitleX("Distance e-#gamma (cm)");        h20.setTitleY("#Deltap  GEN-REC (GeV)"); dg.addDataSet(h20, n); n++;
        h21 = new H2F("h21",50,0,100.0,100,0,5);    h21.setTitleX("Distance e-#gamma (cm)");        h21.setTitleY("ECAL #gamma E (GeV)");    dg.addDataSet(h21, n); n++;
        h22 = new H2F("h23",50,0,100.0,100,0,5);    h22.setTitleX("Distance e-#gamma (cm)");        h22.setTitleY("ECAL #gamma E (GeV)");    dg.addDataSet(h22, n); n++; 
        h23 = new H2F("h23",50,0,9,50,0,5);         h23.setTitleX("Track Electron Momentum (GeV)"); h23.setTitleY("ECAL #gamma E (GeV)");    dg.addDataSet(h23, n); n++;
        h21.setTitle("#Deltap>"+dp1); h22.setTitle("#Deltap<"+dp2);  h23.setTitle("#Deltap<"+dp2); 
        this.getDataGroup().add(dg,0,0,0,run);
        
        dg = new DataGroup(4,2); n=0;       
        String lab[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};
        h30 = new H2F("h30",50,0,9,70,-2,12); h30.setTitleX("Track Electron Momentum (GeV)");h30.setTitleY("PCAL "+lab[0]); dg.addDataSet(h30, n); n++;
        h31 = new H2F("h31",50,0,9,70,-5,30); h31.setTitleX("Track Electron Momentum (GeV)");h31.setTitleY("PCAL "+lab[1]); dg.addDataSet(h31, n); n++;
        h32 = new H2F("h32",50,0,9,70,-2,12); h32.setTitleX("Track Electron Momentum (GeV)");h32.setTitleY("ECIN "+lab[0]); dg.addDataSet(h32, n); n++;
        h33 = new H2F("h33",50,0,9,70,-5,30); h33.setTitleX("Track Electron Momentum (GeV)");h33.setTitleY("ECIN "+lab[1]); dg.addDataSet(h33, n); n++;
        h30.setTitle("#Deltap>"+dp1);h31.setTitle("#Deltap>"+dp1);h32.setTitle("#Deltap>"+dp1);h33.setTitle("#Deltap>"+dp1);
        h40 = new H2F("h40",50,0,9,70,-2,12); h40.setTitleX("Track Electron Momentum (GeV)");h40.setTitleY("PCAL "+lab[0]); dg.addDataSet(h40, n); n++;
        h41 = new H2F("h41",50,0,9,70,-5,30); h41.setTitleX("Track Electron Momentum (GeV)");h41.setTitleY("PCAL "+lab[1]); dg.addDataSet(h41, n); n++;
        h42 = new H2F("h42",50,0,9,70,-2,12); h42.setTitleX("Track Electron Momentum (GeV)");h42.setTitleY("ECIN "+lab[0]); dg.addDataSet(h42, n); n++;
        h43 = new H2F("h43",50,0,9,70,-5,30); h43.setTitleX("Track Electron Momentum (GeV)");h43.setTitleY("ECIN "+lab[1]); dg.addDataSet(h43, n); n++;
        h40.setTitle("#Deltap<"+dp2);h41.setTitle("#Deltap<"+dp2);h42.setTitle("#Deltap<"+dp2);h43.setTitle("#Deltap<"+dp2);
        this.getDataGroup().add(dg,0,0,1,run);
        
        
        dg = new DataGroup(4,3); n=0;       
        h50 = new H2F("h40",50,0.0,10.,50,0.15,0.31); dg.addDataSet(h50, n); n++;      
        h50.setTitleX("True Electron Momentum (GeV)");
        h50.setTitleY("Sampling Fraction E/P");
        
        h51 = new H2F("h41",50,0.0,10.,50,0.15,0.31); dg.addDataSet(h51, n); n++;      
        h51.setTitleX("Track Electron Momentum (GeV)");
        h51.setTitleY("Sampling Fraction E/P");
        
        h52 = new H2F("h42",50,0.0,2.5,50,0.15,0.31); dg.addDataSet(h52, n); n++;      
        h52.setTitleX("ECAL Electron Energy (GeV)");
        h52.setTitleY("Sampling Fraction E/P");
      
        eff1 = new H2F("eff1",30,0,10,30,5,30); eff1.setTitleX("True Electron Momentum (GeV)"); eff1.setTitleY("True Electron Theta (deg)");
        
        eff1a = new H2F("eff1a",30,0,10,30,5,30);
        eff1b = new H2F("eff1b",30,0,10,30,5,30); 
        
        eff2a = new H2F("eff2a",30,0,10,30,5,30);
        eff2b = new H2F("eff2b",30,0,10,30,5,30); 
        
        eff3a = new H2F("eff3a",30,0,10,30,5,30);
        eff3b = new H2F("eff3b",30,0,10,30,5,30);
            
        for (int i=0;i<nh;i++) {effmom[i] = new H1F("effmom",50,0.,10.); effmom[i].setTitleX("MC Electron P (GeV)");}
        for (int i=0;i<nh;i++) {effthe[i] = new H1F("effthe",50,5.,35.); effthe[i].setTitleX("MC Electron #theta (deg)");}
        
        eff1b.setTitle("e-,#gamma,#pi-"); eff2b.setTitle("#chiPID<3.5"); eff3b.setTitle("#DeltaP<"+dp2);       
        eff1b.setTitleX("True Electron Energy (GeV)" );eff1b.setTitleY("True Electron Theta (deg)");
        eff2b.setTitleX("True Electron Energy (GeV)" );eff2b.setTitleY("True Electron Theta (deg)");               
        eff3b.setTitleX("True Electron Energy (GeV)" );eff3b.setTitleY("True Electron Theta (deg)");               
        effmom1.setTitle("G: PC > 0  Y: PC = 1 or 2  R: PC = 1");        
        effmom1.setTitleX("True Electron Energy (GeV)"); effmom1.setTitleY("Efficiency #theta>15");
        effmom2.setTitleX("True Electron Energy (GeV)"); effmom2.setTitleY("Efficiency #theta>15");
        effmom3.setTitleX("True Electron Energy (GeV)"); effmom3.setTitleY("Efficiency #theta>15");                
        effthe1.setTitleX("True Electron Theta (deg)");  effthe1.setTitleY("Efficiency");
        effthe2.setTitleX("True Electron Theta (deg)");  effthe2.setTitleY("Efficiency");
        effthe3.setTitleX("True Electron Theta (deg)");  effthe3.setTitleY("Efficiency");        
        effmom4.setTitle("G: e-,#gamma,#pi-  Y: NO #pi-  R: #chiPID<3.5 B: #DeltaP<"+dp2);        
        effmom4.setTitleX("True Electron Energy (GeV)"); effmom4.setTitleY("Efficiency #theta>15");        
        effthe4.setTitleX("True Electron Theta (deg)");  effthe4.setTitleY("Efficiency");    	
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
    	
        if(eb.readMC(event)) { 
            
        float refE = (float) eb.pmc.get(0).e(); float refP = (float) eb.pmc.get(0).p(); float refTH = (float) (eb.pmc.get(0).theta()*180/Math.PI);
        
        eff1.fill(refE, refTH);
        if(refTH>15)effmom[0].fill(refE);
        effthe[0].fill(refTH);    
        
    	ev.init(event);        	
    	ev.setEventNumber(10);
    	ev.requireOneElectron(false);
   	    ev.setElecTriggerSector(mcsec);
   	    		
        if(ev.procEvent(event)) {
    		elec = getPART(0.0, 11);      int nelec = elec.size();
    		phot = getPART(0.0, 22);      int nphot = phot.size();
    		pim  = getPART(0.0, 212);     int  npim = pim.size();
    		
    		int npart = nelec+nphot+npim;        		        		
    		float esum=0,dp=0,chi2pid=1000;
    		
    		if(npart>0) {if(refTH>15)effmom[4].fill(refE); effthe[4].fill(refTH); eff1a.fill(refE, refTH);}
    		if(npim==0) {if(refTH>15)effmom[5].fill(refE); effthe[5].fill(refTH);}
   		
    		for (Particle p : phot) {esum+=p.e();}
    		        		
    		if(nelec==1) {
        		chi2pid = (float) Math.abs(elec.get(0).getProperty("chi2pid"));
        		dp = (float)(eb.pmc.get(0).p()-elec.get(0).p());
        		
        		if(npim==0&&chi2pid<3.5)  {
        			if(refTH>15)                      effmom[6].fill(refE); 
        			if(refTH>15 && Math.abs(dp)<dp2)  effmom[7].fill(refE); 
        			if(Math.abs(dp)<dp2)   eff3a.fill(refE,refTH);
        			effthe[6].fill(refTH); eff2a.fill(refE,refTH); 	                    
        		}
        		
    			if(esum<0.01) {
    			h00.fill(eb.pmc.get(0).px()-elec.get(0).px());
           		h01.fill(eb.pmc.get(0).py()-elec.get(0).py());
           		h02.fill(eb.pmc.get(0).pz()-elec.get(0).pz());
           		h03.fill(dp);
    			}
    			
                if(esum>0.01) {
    			h10.fill(esum,eb.pmc.get(0).px()-elec.get(0).px());
           		h11.fill(esum,eb.pmc.get(0).py()-elec.get(0).py());
           		h12.fill(esum,eb.pmc.get(0).pz()-elec.get(0).pz());
           		h13.fill(esum,dp);
                }
                              		
          		float e_mom = (float) elec.get(0).p();
                float e_the = (float) Math.toDegrees(elec.get(0).theta());
                float e_phi = (float) Math.toDegrees(elec.get(0).phi());
                
          		List<Particle> elecECAL = ev.getECAL((int)elec.get(0).getProperty("pindex"));
       			float ex=ev.getVar(elecECAL, "x", 1);float ey=ev.getVar(elecECAL, "y", 1);float ez=ev.getVar(elecECAL, "z", 1);
       			
       			for (Particle rp : phot) { //loop over photons from REC::Particle
         			List<Particle> photECAL = ev.getECAL((int)rp.getProperty("pindex"));
           			for (Particle rc : photECAL) { //loop over photon entries in REC::Calorimeter
           			  int il = (int)rc.getProperty("layer");
         		      float ec_the = (float) Math.toDegrees(rc.theta());
           		      float ec_phi = (float) Math.toDegrees(rc.phi());               			
           			  float dthe = e_the-ec_the; float dphi = e_phi-ec_phi;
                      if(il==1) {  //use only PCAL and ECIN entries                        
           			  float gx=(float)rc.getProperty("x");float gy=(float)rc.getProperty("y");float gz=(float)rc.getProperty("z");
           			  float dist = (float) Math.sqrt((ex-gx)*(ex-gx)+(ey-gy)*(ey-gy)+(ez-gz)*(ez-gz));
           			  h20.fill(dist, dp); 
           			  if(Math.abs(dp)>dp1) {h30.fill(e_mom,dthe); h31.fill(e_mom,dphi); h21.fill(dist,esum);}
           			  if(Math.abs(dp)<dp2) {h40.fill(e_mom,dthe); h41.fill(e_mom,dphi); h22.fill(dist,esum); h23.fill(e_mom,esum);}
                      }
                      if(il==4 && Math.abs(dp)>dp1) {h32.fill(e_mom,dthe); h33.fill(e_mom,dphi);}
                      if(il==4 && Math.abs(dp)<dp2) {h42.fill(e_mom,dthe); h43.fill(e_mom,dphi);}
           			}
           		}
           		float en=0; for (Particle p : elecECAL) en += (float) p.getProperty("energy");         		  
           		if(chi2pid<3.5) h51.fill(e_mom, en*1e-3/e_mom);
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
        
        if(rPC.size()>0)                   {if(refTH>15)effmom[1].fill(refE);effthe[1].fill(refTH);}
        if(rPC.size()==1 || rPC.size()==2) {if(refTH>15)effmom[2].fill(refE);effthe[2].fill(refTH);}
        if(rPC.size()==1)                  {if(refTH>15)effmom[3].fill(refE);effthe[3].fill(refTH);
                                           h50.fill(refP, ecalE/refP); 
                                           h52.fill(ecalE,ecalE/refP);}                       
    }  

        
    }
    
    public void plotMCHistos() {
    	int index = 0;
      
        EmbeddedCanvas c = null; int ncd;
       
        index=getDetectorTabNames().indexOf("GEN-REC");
        c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(0)); c.divide(4,3); ncd=-1;           
		c.setAxisTitleSize(18);
		c.setAxisLabelSize(18);
		
		ncd++; c.cd(ncd) ; c.draw(h00);  
		ncd++; c.cd(ncd) ; c.draw(h01);  
		ncd++; c.cd(ncd) ; c.draw(h02);  
		ncd++; c.cd(ncd) ; c.draw(h03);  
		
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.draw(h10);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.draw(h11);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.draw(h12);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.draw(h13);  

        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setLog(true); c.draw(h20);
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setLog(true); c.draw(h21);
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setLog(true); c.draw(h22);
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setLog(true); c.draw(h23);

        index=getDetectorTabNames().indexOf("DTHE,DPHI");
        c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); c.divide(4,2); ncd=-1;           
		c.setAxisTitleSize(18);
		c.setAxisLabelSize(18);

		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h30);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h31);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h32);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h33);  
		
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h40);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h41);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h42);  
		ncd++; c.cd(ncd) ; c.getPad(ncd).getAxisZ().setLog(true); c.getPad(ncd).setTitleFontSize(18); c.draw(h43);  

        index=getDetectorTabNames().indexOf("EFF");
        c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); c.divide(4,3); ncd=-1;           
		c.setAxisTitleSize(18);
		c.setAxisLabelSize(18);

        ncd++; c.cd(ncd); c.draw(h50);          
        ncd++; c.cd(ncd); c.draw(h51);  
        ncd++; c.cd(ncd); c.draw(h52);  
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.15,0.31); c.draw(meanGraph);          
        SFFunction sf = new SFFunction("esf",-11,eb.ccdb,0.1,2.5);  c.draw(sf,"same");        
    
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.05);  c.getPad(ncd).setTitleFontSize(18);
        c.draw(effmom1);  c.draw(effmom2,"same"); c.draw(effmom3,"same");
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.05);  c.getPad(ncd).setTitleFontSize(18);
        c.draw(effthe1);  c.draw(effthe2,"same"); c.draw(effthe3,"same");
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.015); c.getPad(ncd).setTitleFontSize(18);
                               c.getPad(ncd).getAxisX().setGrid(true);       c.getPad(ncd).getAxisY().setGrid(true);                               
        c.draw(effmom4);  c.draw(effmom5,"same"); c.draw(effmom6,"same"); c.draw(effmom7,"same"); 
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisY().setRange(0.5,1.015); c.getPad(ncd).setTitleFontSize(18);
        c.draw(effthe4);  c.draw(effthe5,"same"); c.draw(effthe6,"same"); 

        
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setRange(0.9,1.01); c.getPad(ncd).setTitleFontSize(18); c.draw(eff1b);
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setRange(0.9,1.01); c.getPad(ncd).setTitleFontSize(18); c.draw(eff2b);
        ncd++; c.cd(ncd); c.getPad(ncd).getAxisZ().setRange(0.5,1.01); c.getPad(ncd).setTitleFontSize(18); c.draw(eff3b);
    }
    
    @Override
    public void plotEvent(DataEvent de) {
        analyze();         
    }
    
    public GraphErrors fith52() {
	    ParallelSliceFitter fitter = new ParallelSliceFitter(h52);
        fitter.fitSlicesX();
        
        GraphErrors sigGraph = fitter.getSigmaSlices();   
        sigGraph.setTitleX("ECAL Electron Energy (GeV)");  sigGraph.setTitleY("SIGMA E/P"); 	
        sigGraph.setMarkerSize(4); sigGraph.setMarkerStyle(1);
        
        GraphErrors meanGraph = fitter.getMeanSlices();  
        meanGraph.setTitleX("ECAL Electron Energy (GeV)"); meanGraph.setTitleY("Sampling Fraction E/P"); 	
        meanGraph.setMarkerSize(4); meanGraph.setMarkerStyle(1); 
          
        int npts = meanGraph.getDataSize(0)  ;
        double[] xm  = new double[npts];
        double[] ym  = new double[npts];
        double[] yme = new double[npts];
        double[] xs  = new double[npts];
        double[] ys  = new double[npts];
        double[] yse = new double[npts];
        
        for (int i=0; i<meanGraph.getDataSize(0); i++) {
      	   xm[i] = 1/Math.sqrt(meanGraph.getDataX(i));
      	   ym[i] = meanGraph.getDataY(i);
      	  yme[i] = meanGraph.getDataEY(i);      			  
        } 
        
        GraphErrors resGraph = new GraphErrors();
        resGraph.setTitleX("1/E (GeV^-^1)");  resGraph.setTitleY("(#sigma / E)^2"); 
        resGraph.setMarkerSize(4); meanGraph.setMarkerStyle(1);
        
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
        return meanGraph;
    }
    
    public void geteff() {
        H1F effmom1 = H1F.divide(effmom[1],effmom[0]); effmom1.setFillColor(3);
        H1F effmom2 = H1F.divide(effmom[2],effmom[0]); effmom2.setFillColor(5);
        H1F effmom3 = H1F.divide(effmom[3],effmom[0]); effmom3.setFillColor(2);
        H1F effmom4 = H1F.divide(effmom[4],effmom[0]); effmom4.setFillColor(3);
        H1F effmom5 = H1F.divide(effmom[5],effmom[0]); effmom5.setFillColor(5);
        H1F effmom6 = H1F.divide(effmom[6],effmom[0]); effmom6.setFillColor(2);
        H1F effmom7 = H1F.divide(effmom[7],effmom[0]); effmom7.setFillColor(1);        
        H1F effthe1 = H1F.divide(effthe[1],effthe[0]); effthe1.setFillColor(3);
        H1F effthe2 = H1F.divide(effthe[2],effthe[0]); effthe2.setFillColor(5);
        H1F effthe3 = H1F.divide(effthe[3],effthe[0]); effthe3.setFillColor(2);
        H1F effthe4 = H1F.divide(effthe[4],effthe[0]); effthe4.setFillColor(3);
        H1F effthe5 = H1F.divide(effthe[5],effthe[0]); effthe5.setFillColor(5);
        H1F effthe6 = H1F.divide(effthe[6],effthe[0]); effthe6.setFillColor(2);        
        eff1b=eff1b.divide(eff1a, eff1); eff2b=eff2b.divide(eff2a, eff1); eff3b=eff3b.divide(eff3a, eff1);    	
    }
        
    @Override
    public void timerUpdate() {  
    	meanGraph = fith52();
    	geteff();
    }


}
