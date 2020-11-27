package org.clas.tools;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorData;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.math.Func1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.eb.EBCCDBConstants;
import org.jlab.rec.eb.EBCCDBEnum;
import org.jlab.rec.eb.EBRadioFrequency;
import org.jlab.rec.eb.SamplingFractions;
import org.jlab.service.eb.EBAnalyzer;
import org.jlab.service.eb.EBEngine;
import org.jlab.service.eb.EBMatching;
import org.jlab.service.eb.EBTBEngine;
import org.jlab.service.eb.EventBuilder;
import org.jlab.service.ec.ECEngine;
import org.jlab.utils.groups.IndexedList;

public class EBMC extends EBTBEngine {
	
	public EventBuilder eb = null;
	EBEngine         ebe = new EBEngine("EBMC");
	public EBCCDBConstants  ccdb;
	EBMatching       ebm;
	EBRadioFrequency rf;
	
    public List<List<DetectorResponse>>     unmatchedResponses = new ArrayList<>(); 
    
    IndexedList<List<DetectorResponse>>  singleNeutrals = new IndexedList<>(1);
    IndexedList<List<DetectorResponse>>      singleMIPs = new IndexedList<>(1);
    
    DetectorParticle p1 = new DetectorParticle();
    DetectorParticle p2 = new DetectorParticle();
    
    public List<Particle> pmc = new ArrayList<>();;
    
    public double distance11,distance12,distance21,distance22;
    public double e1,e2,e1c,e2c,cth,cth1,cth2;
    public double X,tpi2,cpi0,ppx1,ppy1,ppz1,refE,refP,refTH,refPH;
    public double x1,y1,x2,y2;
    public double epc,eec1,eec2;
    public LorentzVector VG1,VG2;
//    public static int[] ip1,ip2,is1,is2;
    public int[]    iis = new int[2];
    public int[][]  iip = new int[2][3];
    public double[][] x = new double[2][3];
    public double[][] y = new double[2][3];
    public double[][] t = new double[2][3];
    public double[] distance1 = new double[2];
    public double[] distance2 = new double[2];
    public double mpi0 = 0.1349764;
    public double melec = 0.000511;
    public double mneut = 0.93957;
    public float c = 29.98f; 
    public String geom = "2.4";
    public String config = null;
    public double SF1 = 0.27, SF1db=0.27;
    public double SF2 = 0.27, SF2db=0.27;
    public boolean isMC = true;
    public H1F h5,h6,h7,h8;
  
    public int n2mc=0;

    public int[] mip = {0,0,0,0,0,0};
    public int runNumber=11;
    
    int photonMult = 12;
    
    Vector3 vtx = new Vector3(0,0,0);
    double opa = 0;
    
    String eventBank        = "REC::Event";    
    String particleBank     = "REC::Particle";
    String calorimeterBank  = "REC::Calorimeter";  
    
    private EmbeddedCanvasTabbed      detectorCanvas = null;  
    private ArrayList<String>       detectorTabNames = new ArrayList();
    
    public EBMC() {  
    	initBankNames();
    	detectorCanvas = new EmbeddedCanvasTabbed();
    }
    
    public void setDetectorCanvas(EmbeddedCanvasTabbed canvas) {
        detectorCanvas = canvas;
    }
    
    public void setDetectorTabNames(String... names) {
        for(String name : names) {
            detectorTabNames.add(name);
        }
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed(names);
        setDetectorCanvas(canvas);
    }
    
    public EmbeddedCanvasTabbed getDetectorCanvas() {
        return detectorCanvas;
    }
    
    public ArrayList<String> getDetectorTabNames() {
        return detectorTabNames;
    }
    
    public void getCCDB(int runno) {
    	System.out.println("ECpart.setCCDB("+runno+")");
    	ebe.init();
    	this.ccdb = new EBCCDBConstants(runno,ebe.getConstantsManager());    	
    }
    
    public boolean readMC(DataEvent event) {
    	
        pmc.clear();
        if(event.hasBank("MC::Particle")) {
            DataBank bank = event.getBank("MC::Particle");
            for (int i=0; i<bank.rows(); i++) {
            	double px = bank.getFloat("px",i);
            	double py = bank.getFloat("py",i);
            	double pz = bank.getFloat("pz",i);
            	double vx = bank.getFloat("vx",i);
            	double vy = bank.getFloat("vy",i);
            	double vz = bank.getFloat("vz",i);
            	int   pid = bank.getInt("pid",i);  
                pmc.add(new Particle(pid, px, py, pz, vx, vy, vz));
                if(i==0) vtx.setXYZ(vx, vy, vz);
            }
            n2mc++;
            return true;
        }                  
        return false;
    }
    
	// readEC: Copies relevant parts of EBEngine.processDataEvent    
    public void  readEC(DataEvent de, String bank){    	    	
        rf = new EBRadioFrequency(ccdb);    	
        eb = new EventBuilder(ccdb);    	   	
        eb.initEvent(); //don't bother with event header    	
        eb.getEvent().getEventHeader().setRfTime(rf.getTime(de)+ccdb.getDouble(EBCCDBEnum.RF_OFFSET));
        eb.addDetectorResponses(CalorimeterResponse.readHipoEvent(de, "ECAL::clusters", DetectorType.ECAL,null));  	
        eb.getPindexMap().put(0, 0); 
        eb.getPindexMap().put(1, 0);         
    } 
    
	// getNeutralPart: Copies relevant parts of EBEngine.processDataEvent 
    // Note for MC w/o charged trigger particles EBAnalyzer.foundTriggerTime will be false and SF not applied to PID=22
    public List<DetectorParticle> getNeutralPart() {
    	eb.processNeutralTracks();    	
    	EBAnalyzer analyzer = new EBAnalyzer(ccdb, rf);
        analyzer.processEvent(eb.getEvent());
        if(eb.getEvent().getParticles().size()>0) {
            Collections.sort(eb.getEvent().getParticles());
            eb.setParticleStatuses();  
            return eb.getEvent().getParticles();
        } 
     	return new ArrayList<DetectorParticle>() ;     	
    }
    
    public void getRECBanks(DataEvent de, EventBuilder eb) {
        DataBank bankP = DetectorData.getDetectorParticleBank(eb.getEvent().getParticles(), de, particleBank);
        de.appendBanks(bankP);
      
        DataBank bankEve = DetectorData.getEventBank(eb.getEvent(), de, eventBank);
        de.appendBanks(bankEve);

        List<DetectorResponse> calorimeters = eb.getEvent().getCalorimeterResponseList();
        if(calorimeterBank!=null && calorimeters.size()>0) {
            DataBank bankCal = DetectorData.getCalorimeterResponseBank(calorimeters, de, calorimeterBank);
            de.appendBanks(bankCal);
        }
    	
    }
    
    public void getUnmatchedResponses() {        
        unmatchedResponses.clear();
        unmatchedResponses.add(eb.getUnmatchedResponses(null, DetectorType.ECAL,1)); //PCAL
        unmatchedResponses.add(eb.getUnmatchedResponses(null, DetectorType.ECAL,4)); //ECIN
        unmatchedResponses.add(eb.getUnmatchedResponses(null, DetectorType.ECAL,7)); //ECOU
    }    
    
    public void getMIPResponses() {        
        getUnmatchedResponses();
     	getSingleMIPResponses(); 
    }
    
    public void getSingleMIPResponses() {
        List<DetectorResponse> rPC = new ArrayList<>();
        singleMIPs.clear();
        for (int is=1; is<7; is++) {
            rPC = DetectorResponse.getListBySector(unmatchedResponses.get(0),  DetectorType.ECAL, is);  //look in PCAL only          
            if(rPC.size()==1&&mip[is-1]==1) singleMIPs.add(rPC,is);
        }     	
    }   
    
    public double getEcalEnergy(int sector){
        iis[0]=0;epc=0;eec1=0;eec2=0;
        return processSingleMIP(doHitMatching(getMIParticles(sector)));
    }   
    
    public double processSingleMIP(List<DetectorParticle> particles) {
        if (particles.size()==0) return 0.0;
      	return particles.get(0).getEnergy(DetectorType.ECAL);
    }  
    
    public  List<DetectorResponse> getPCResponses(int sector) {
    	return DetectorResponse.getListBySector(unmatchedResponses.get(0), DetectorType.ECAL, sector);   	
    }
    
    public List<DetectorParticle> getMIParticles(int sector) {
    	
        List<DetectorParticle> particles = new ArrayList<>();          
        List<DetectorResponse>      rPC  = new ArrayList<>();        
        rPC = getPCResponses(sector);
        if (rPC.size()==1) particles.add(DetectorParticle.createNeutral(rPC.get(0),vtx)); // Zero torus/solenoid fields
        return particles;
    }
    
    public void getNeutralResponses() {        
        getUnmatchedResponses();
       	getSingleNeutralResponses(); // For two-photon decays in different sectors     	
    }
          
    public void getSingleNeutralResponses() {
        List<DetectorResponse> rPC = new ArrayList<>();
        singleNeutrals.clear();
        for (int is=1; is<7; is++) {
            rPC = DetectorResponse.getListBySector(unmatchedResponses.get(0),  DetectorType.ECAL, is); //look in PCAL only
            if(rPC.size()==1&&mip[is-1]!=1) singleNeutrals.add(rPC,is);
        } 
    }
    
    public List<DetectorResponse> findSecondPhoton(int sector) {
        int neut=0, isave=0;
        for (int is=sector+1; is<7; is++) {
            if(singleNeutrals.hasItem(is)) {neut++; isave=is;}
        }
        return (neut==1) ? singleNeutrals.getItem(isave):new ArrayList<DetectorResponse>();
    }
     
    public double getTwoPhotonInvMass(int sector){
        iis[0]=0;iis[1]=-1;
        return processTwoPhotons(doHitMatching(getNeutralParticles(sector)));
//        return processTwoPhotons(getNeutralPart());
    }   
    
    // getNeutralParticles: Similar to EBMatching.findNeutrals(1)
    // Note here I am requiring 2 PCAL responses, NOT 2 complete ECAL particles
    public List<DetectorParticle> getNeutralParticles(int sector) {
              
        List<DetectorResponse>      rEC  = new ArrayList<>();        
        List<DetectorParticle> particles = new ArrayList<>();          

        rEC = DetectorResponse.getListBySector(unmatchedResponses.get(0), DetectorType.ECAL, sector); //get PCAL responses

        switch (rEC.size()) {
        case 1:  List<DetectorResponse> rEC2 = findSecondPhoton(sector);
                if (rEC2.size()>0) {
                   particles.add(DetectorParticle.createNeutral(rEC.get(0),vtx));                    // make neutral particle 1 from PCAL sector                   
                   particles.add(DetectorParticle.createNeutral(rEC2.get(0),vtx)); return particles; // make neutral particle 2 from other PCAL sector
                }
                break;
        case 2: particles.add(DetectorParticle.createNeutral(rEC.get(0),vtx));                   // make neutral particle 1 from PCAL sector
                particles.add(DetectorParticle.createNeutral(rEC.get(1),vtx)); return particles; // make neutral particle 2 from PCAL sector
//        case 3: particles.add(DetectorParticle.createNeutral(rEC.get(0)));                   // make neutral particle 1 from PCAL sector
//        		particles.add(DetectorParticle.createNeutral(rEC.get(1))); return particles; // make neutral particle 2 from PCAL sector
        }
        
       return particles;
    }
    
    public List<DetectorParticle> doHitMatching(List<DetectorParticle> particles) {
    	
        int ii=0;
        for (DetectorParticle p: particles) {
            DetectorResponse rPC = p.getDetectorResponses().get(0); //get 1st PCAL responses
                   epc = rPC.getEnergy();
            iip[ii][0] = rPC.getHitIndex();  
            iis[ii]    = rPC.getDescriptor().getSector();
              x[ii][0] = rPC.getPosition().x();
              y[ii][0] = rPC.getPosition().y();
              t[ii][0] = rPC.getTime()-rPC.getPath()/c;;
             
            distance1[ii] = doPCECMatch(p,ii,"Inner");
            distance2[ii] = doPCECMatch(p,ii,"Outer");
            ii++;            
        }
   
        return particles;
    }  
    
    // doPCECMatch: Similar to EBMatching.addResponsesECAL
    public double doPCECMatch(DetectorParticle p, int ii, String io) {
        
        int index = 0;
        double distance = -10;
        List<DetectorResponse> rEC = new ArrayList<>();        
        
        int is = p.getDetectorResponses().get(0).getDescriptor().getSector();
       
        switch (io) {   
        
        case "Inner": rEC = DetectorResponse.getListBySector(unmatchedResponses.get(1), DetectorType.ECAL, is);
                      index  = p.getDetectorHit(rEC,DetectorType.ECAL,4,eb.ccdb.getDouble(EBCCDBEnum.ECIN_MATCHING));
                      if(index>=0){p.addResponse(rEC.get(index),true); rEC.get(index).setAssociation(0);
                      iip[ii][1] = rEC.get(index).getHitIndex(); 
                      x[ii][1]   = rEC.get(index).getPosition().x(); y[ii][1] = rEC.get(index).getPosition().y();                      
                      t[ii][1]   = rEC.get(index).getTime()-rEC.get(index).getPath()/c;
                      distance = p.getDistance(rEC.get(index)).length(); eec1=rEC.get(index).getEnergy();} break;
                      
        case "Outer": rEC = DetectorResponse.getListBySector(unmatchedResponses.get(2), DetectorType.ECAL, is); 
                      index  = p.getDetectorHit(rEC,DetectorType.ECAL,7,eb.ccdb.getDouble(EBCCDBEnum.ECOUT_MATCHING));
                      if(index>=0){p.addResponse(rEC.get(index),true); rEC.get(index).setAssociation(0);
                      iip[ii][2] = rEC.get(index).getHitIndex();
                      x[ii][2]   = rEC.get(index).getPosition().x(); y[ii][2] = rEC.get(index).getPosition().y(); 
                      t[ii][2]   = rEC.get(index).getTime()-rEC.get(index).getPath()/c;                      
                      distance = p.getDistance(rEC.get(index)).length(); eec2=rEC.get(index).getEnergy();}
                      
        }
        
        return distance;        
    }
        
    public double processTwoPhotons(List<DetectorParticle> particles) {
        
        if (particles.size()<2) return 0.0;
        
        p1 = particles.get(0);  //Photon 1
        p2 = particles.get(1);  //Photon 2
       
        Vector3 n1 = p1.vector(); n1.unit();
        Vector3 n2 = p2.vector(); n2.unit();
                
        e1 = p1.getEnergy(DetectorType.ECAL);
        e2 = p2.getEnergy(DetectorType.ECAL);
        
        SF1db = SamplingFractions.getMean(22, p1, this.ccdb);
        e1c = e1/SF1db;
        Particle g1 = new Particle(22,n1.x()*e1c,n1.y()*e1c,n1.z()*e1c);        
        VG1 = new LorentzVector(n1.x()*e1c,n1.y()*e1c,n1.z()*e1c,e1c);
        
        SF2db = SamplingFractions.getMean(22, p2, this.ccdb);
        e2c = e2/SF2db;
        Particle g2 = new Particle(22,n2.x()*e2c,n2.y()*e2c,n2.z()*e2c);
        VG2 = new LorentzVector(n2.x()*e2c,n2.y()*e2c,n2.z()*e2c,e2c);
       
        cth1 = Math.cos(g1.theta());
        cth2 = Math.cos(g2.theta());
         cth = g1.cosTheta(g2);         
           X = (e1c-e2c)/(e1c+e2c);
        tpi2 = 2*mpi0*mpi0/(1-cth)/(1-X*X);
        cpi0 = (e1c*cth1+e2c*cth2)/Math.sqrt(e1c*e1c+e2c*e2c+2*e1c*e2c*cth);
        
        g1.combine(g2, +1);
                               
        return goodPhotons(photonMult,p1,p2) ? g1.mass2():0.0;
    }
    
    public boolean goodPhotons(int test, DetectorParticle pp1, DetectorParticle pp2) {
    	
       // Require two photons in PCAL  
       boolean pc12 = DetectorResponse.getListByLayer(pp1.getDetectorResponses(),DetectorType.ECAL, 1).size()!=0  &&
                      DetectorResponse.getListByLayer(pp2.getDetectorResponses(),DetectorType.ECAL, 1).size()!=0;
        
       // Require photon 1 cluster in PCAL and ECin
       boolean pcec1 = DetectorResponse.getListByLayer(pp1.getDetectorResponses(),DetectorType.ECAL, 1).size()!=0 &&           
                       DetectorResponse.getListByLayer(pp1.getDetectorResponses(),DetectorType.ECAL, 4).size()!=0; 

       // Require photon 2 cluster in PCAL and ECin
       boolean pcec2 = DetectorResponse.getListByLayer(pp2.getDetectorResponses(),DetectorType.ECAL, 1).size()!=0 &&
                       DetectorResponse.getListByLayer(pp2.getDetectorResponses(),DetectorType.ECAL, 4).size()!=0;  
       
       switch (test) {
           case 121:  return pc12 && pcec1;          
           case 122:  return pc12 && pcec2;        
           case 120:  return pc12 && (pcec1 || pcec2); 
           case 1212: return pc12 && (pcec1 && pcec2);
       }
       
       return pc12;
       
    }
    
    public void setGoodPhotons(int num) {
    	this.photonMult = num;
    }
    
    public void setGeom(String geom) {
        this.geom = geom;
    }
    
    public void setConfig(String config) {
        this.config = config;
    }
    
    private int getRunNumber(DataEvent event) {
        DataBank bank = event.getBank("RUN::config"); 
        return (bank!=null) ? bank.getInt("run",0):this.runNumber;
    }
    
    public static double getSF(String geom, double e) {
        switch (geom) {
        case "2.4": return 0.268*(1.0510 - 0.0104/e - 0.00008/e/e); 
        case "2.5": return 0.250*(1.0286 - 0.0150/e + 0.00012/e/e);
        }
        return Double.parseDouble(geom);
    } 
    
    public static class SFFunction extends Func1D{
    	
        EBCCDBConstants ccdb = new EBCCDBConstants();
   	    DetectorParticle p = new DetectorParticle(); 
        int pid;
  	    
        public SFFunction(String name, int pid, EBCCDBConstants ccdb, double min, double max) {
            super(name, min, max);
            this.ccdb = ccdb;
            this.pid  = pid;
            
            p.addResponse(new CalorimeterResponse(1,1,0));
            p.getDetectorResponses().get(0).getDescriptor().setType(DetectorType.ECAL);
        }
        @Override
        public double evaluate(double x){        	 
        	 p.getDetectorResponses().get(0).setEnergy(x);
       	     return  SamplingFractions.getMean(pid, p, ccdb);
        }
    }
    
    public void setThresholds(String part, ECEngine engine) {
    	switch (part) {
    	case "Electron_lo":engine.setStripThresholds(10,10,10);
                           engine.setPeakThresholds(30,30,30);
                           engine.setClusterCuts(7,15,20); break;    		
    	case "Electron_hi":engine.setStripThresholds(20,50,50);
                           engine.setPeakThresholds(40,80,80);
                           engine.setClusterCuts(7,15,20); break;    		
    	case      "Pizero":engine.setStripThresholds(10,9,8);  
                           engine.setPeakThresholds(18,20,15); 
                           engine.setClusterCuts(7,15,20); break;
    	case        "Test":engine.setStripThresholds(10,200,8);  
                           engine.setPeakThresholds(18,500,15); 
                           engine.setClusterCuts(7,15,20); 
    	}
    }
    
    public void dumpGraph(String filename, GraphErrors graph) {
    	PrintWriter writer = null;
		try 
		{
			writer = new PrintWriter(filename);
		} 
		catch (FileNotFoundException e) 
		{			// TODO Auto-generated catch block
			e.printStackTrace();
		}				
    	for (int i=0; i<graph.getDataSize(0); i++) writer.println(String.format("%1$.3f %2$.3f %3$.3f",graph.getDataX(i),graph.getDataY(i),graph.getDataEY(i))); 
    	writer.close();    	
    }
    
    public void dropBanks(DataEvent event) {    	
        if(event.hasBank("ECAL::clusters")) event.removeBanks("ECAL::hits","ECAL::peaks","ECAL::clusters","ECAL::calib","ECAL::moments");
        if(event.hasBank("ECAL::clusters")) event.removeBank("ECAL::clusters");
        if(event.hasBank("ECAL::hits"))     event.removeBank("ECAL::hits");
        if(event.hasBank("ECAL::peaks"))    event.removeBank("ECAL::peaks");
        if(event.hasBank("ECAL::calib"))    event.removeBank("ECAL::calib");
        if(event.hasBank("ECAL::moments"))  event.removeBank("ECAL::moments");  	
    }   
    
    public void initGraphics() {
        GStyle.getAxisAttributesX().setTitleFontSize(14);
        GStyle.getAxisAttributesX().setLabelFontSize(14);
        GStyle.getAxisAttributesY().setTitleFontSize(14);
        GStyle.getAxisAttributesY().setLabelFontSize(14);
        GStyle.getAxisAttributesZ().setLabelFontSize(14); 
        GStyle.getAxisAttributesX().setAxisGrid(false);
        GStyle.getAxisAttributesY().setAxisGrid(false);
        GStyle.getAxisAttributesX().setLabelFontName("Avenir");
        GStyle.getAxisAttributesY().setLabelFontName("Avenir");
        GStyle.getAxisAttributesZ().setLabelFontName("Avenir");
        GStyle.getAxisAttributesX().setTitleFontName("Avenir");
        GStyle.getAxisAttributesY().setTitleFontName("Avenir");
        GStyle.getAxisAttributesZ().setTitleFontName("Avenir");
    }

}
    
    

