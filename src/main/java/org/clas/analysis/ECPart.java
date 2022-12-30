package org.clas.analysis;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.JFrame;

import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorData;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.EventFilter;
import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.evio.EvioDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;


//import org.jlab.io.hipo3.Hipo3DataSource;
import org.jlab.service.eb.EBAnalyzer;
//import org.jlab.service.eb.EBConstants;
import org.jlab.service.eb.EBEngine;
import org.jlab.service.eb.EBMatching;
import org.jlab.service.eb.EBTBEngine;
import org.clas.service.ec.ECEngine;
import org.jlab.service.eb.EventBuilder;
import org.jlab.utils.groups.IndexedList;
import org.jlab.rec.eb.EBCCDBConstants;
import org.jlab.rec.eb.EBCCDBEnum;
import org.jlab.rec.eb.EBRadioFrequency;
import org.jlab.rec.eb.EBUtil;
import org.jlab.rec.eb.SamplingFractions;
import org.clas.tools.EmbeddedCanvasTabbed;
import org.clas.tools.Event;

public class ECPart extends EBEngine {
	
	public EventBuilder eb = null;
	EBEngine         ebe = new EBEngine("ECPART");
	
	Event             ev = new Event();
	
	EBCCDBConstants  ccdb;
	EBMatching       ebm;
	EBRadioFrequency rf;
	
    public List<List<DetectorResponse>>     unmatchedResponses = new ArrayList<>();     
    IndexedList<List<DetectorResponse>>         singleNeutrals = new IndexedList<>(1);
    IndexedList<List<DetectorResponse>>             singleMIPs = new IndexedList<>(1);
    
    DetectorParticle p1 = new DetectorParticle();
    DetectorParticle p2 = new DetectorParticle();
    
    public List<Particle> pmc = new ArrayList<>();;
    public List<Vector3>  pmv = new ArrayList<>();
    
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
    public double[][] z = new double[2][3];
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
    public boolean isMC = true, hasStartTime = false;;
    public H1F h5,h6,h7,h8;
  
    public int n2mc=0, MCpid=11, MCsec=2;
    
    public float mip[][] = {{0,0,0,0,0,0},{0,0,0,0,0,0}};
    public int runNumber=11;
    public float starttime=0;
    
    public boolean FTOFveto;
    
    int photonMult = 12;
    
    Vector3 vtx = new Vector3(0,0,0);
    double opa = 0;
    
    String eventBank        = "REC::Event";    
    String particleBank     = "REC::Particle";
    String calorimeterBank  = "REC::Calorimeter";  
    
    private EmbeddedCanvasTabbed      detectorCanvas = new EmbeddedCanvasTabbed();  
    private ArrayList<String>       detectorTabNames = new ArrayList();
    
    public ECPart() {  
    	super("EBMC");
    	initBankNames();
    }
    
	public void setMCpid(int val) {
		MCpid = val;
	}
	
	public void setMCsec(int val) {
		MCsec = val;
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
    
    public void setCCDB(int runno) {
    	System.out.println("ECpart.setCCDB("+runno+")");
    	ebe.init();
    	this.ccdb = new EBCCDBConstants(runno,ebe.getConstantsManager());    	
    }   
    
    public boolean hasStartTime(int pid) {
    	return (pid==11 || pid==-11 || pid==-211 || pid==211);
    }
    
    public float getStartTime(DataEvent event) { 
    	float time=0;
    	if(event.hasBank("REC::Event")) {
    		DataBank bank = event.getBank("REC::Event");
    		time = bank.getFloat("startTime", 0);	
    	}
    	return time;
    }
    
    public boolean readMC(DataEvent event) {    	
        pmc.clear(); pmv.clear(); hasStartTime=false;
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
                if(pid==MCpid) {
                	pmc.add(new Particle(pid, px, py, pz, vx, vy, vz)); pmv.add(new Vector3(px,py,pz));
                	refE = pmc.get(0).e(); 
                	refP = pmc.get(0).p();
                	refTH = Math.toDegrees(pmc.get(0).theta());
                	refPH = Math.toDegrees(pmc.get(0).phi());	
                }
                if(i==0) {vtx.setXYZ(vx, vy, vz); hasStartTime = hasStartTime(pid); starttime = getStartTime(event);}
            }
            n2mc++;
            return true;
        }                  
        return false;
    }
    
	// Copies relevant parts of EBEngine.processDataEvent    
    public boolean  processDataEvent(DataEvent de){ 
    	
        eb = new EventBuilder(ccdb);    	   	
        eb.initEvent(); //don't bother with event header  
       
        rf = new EBRadioFrequency(ccdb);    	
        eb.getEvent().getEventHeader().setRfTime(rf.getTime(de)+ccdb.getDouble(EBCCDBEnum.RF_OFFSET));
       
        eb.addDetectorResponses(CalorimeterResponse.readHipoEvent(de, "ECAL::clusters", DetectorType.ECAL,null));
       
        eb.getPindexMap().put(0, 0); 
        eb.getPindexMap().put(1, 0); 
        
//        eb.processHitMatching();
//        eb.assignTrigger(); 
        
        return true;
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
            if(rPC.size()==1&&(FTOFveto?mip[0][is-1]>0&&mip[1][is-1]>0:true)) singleMIPs.add(rPC,is);
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
            if(rPC.size()==1&&(FTOFveto?mip[0][is-1]==0&&mip[1][is-1]==0:true)) singleNeutrals.add(rPC,is);
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
    	iis = new int[2]; for (int i=0; i<2; i++) { for (int j=0; j<3; j++) {iip[i][j]=-1;}};
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
              z[ii][0] = rPC.getPosition().z();
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
                      x[ii][1]   = rEC.get(index).getPosition().x(); y[ii][1] = rEC.get(index).getPosition().y(); z[ii][1] = rEC.get(index).getPosition().z();                      
                      t[ii][1]   = rEC.get(index).getTime()-rEC.get(index).getPath()/c;
                      distance = p.getDistance(rEC.get(index)).length(); eec1=rEC.get(index).getEnergy();} break;
                      
        case "Outer": rEC = DetectorResponse.getListBySector(unmatchedResponses.get(2), DetectorType.ECAL, is); 
                      index  = p.getDetectorHit(rEC,DetectorType.ECAL,7,eb.ccdb.getDouble(EBCCDBEnum.ECOUT_MATCHING));
                      if(index>=0){p.addResponse(rEC.get(index),true); rEC.get(index).setAssociation(0);
                      iip[ii][2] = rEC.get(index).getHitIndex();
                      x[ii][2]   = rEC.get(index).getPosition().x(); y[ii][2] = rEC.get(index).getPosition().y(); z[ii][2] = rEC.get(index).getPosition().z();  
                      t[ii][2]   = rEC.get(index).getTime()-rEC.get(index).getPath()/c;                      
                      distance = p.getDistance(rEC.get(index)).length(); eec2=rEC.get(index).getEnergy();}
                      
        }
        
        return distance;        
    }
        
    public double processTwoPhotons(List<DetectorParticle> particles) {
        
        if (particles.size()<2) return 0.0;
        
        p1 = particles.get(0);  //Photon 1
        p2 = particles.get(1);  //Photon 2
        
        Vector3 n1 = new Vector3(); Vector3 n2 = new Vector3();
        n1.copy(p1.vector()); n2.copy( p2.vector());
        n1.unit(); n2.unit();        
                
        e1 = p1.getEnergy(DetectorType.ECAL);
        e2 = p2.getEnergy(DetectorType.ECAL);
        
        SF1db = SamplingFractions.getMean(22, p1, this.ccdb);
        e1c = e1/SF1db;
        Particle g1 = new Particle(22,n1.x()*e1c,n1.y()*e1c,n1.z()*e1c);        
//        VG1 = new LorentzVector(n1.x()*e1c,n1.y()*e1c,n1.z()*e1c,e1c);
        
        SF2db = SamplingFractions.getMean(22, p2, this.ccdb);
        e2c = e2/SF2db;
        Particle g2 = new Particle(22,n2.x()*e2c,n2.y()*e2c,n2.z()*e2c);
//        VG2 = new LorentzVector(n2.x()*e2c,n2.y()*e2c,n2.z()*e2c,e2c);
       
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
    	case        "Test":engine.setStripThresholds(20,20,20);  
                           engine.setPeakThresholds(18,20,15); 
                           engine.setClusterCuts(3,7,7); 
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
    
    public List<Particle> getPART(double thr, int pid) {   	
    	List<Particle> olist = new ArrayList<Particle>();    
    	for (Particle p : ev.getParticle(pid)) {
    		short status = (short) p.getProperty("status");
    		if(status>=2000 && status<3000 && p.p()>thr) olist.add(p); 
    	}          	
       return olist;    	
    }
    
    
    public void electronDemo(String[] args) {
    	
        HipoDataSource reader = new HipoDataSource();
        ECEngine       engine = new ECEngine();
        List<DetectorResponse> rPC = null;
        
        setDetectorTabNames("GEN-REC","DTHE,DPHI","EFF");
        int mcsec=2;
     
//        String hipoPath = "/Users/colesmith/imac/clas12/sim/radstudy/";
        String hipoPath = "/Users/colesmith/clas12/sim/";
        String hipo3File0 = "fc-elec-40k-s5-r2.hipo";
        String hipo3File3 = "fc-elec-40k-s5-r2.hipo";
        String hipo3File1 = "fc-ecpcsc-elec-s5-20k.hipo";
        String hipo4File1 = "fc-ecpcsc-elec-s5-20k.hipo";
        String hipo3File2 = "clasdispr-large.hipo";
        String hipo4File4 = "clas12-radstudy-5gev-22deg-pm8.hipo";
        
        String hipoFile = "out-rga-fall2018-torsol-100k.hipo"; mcsec=1;
//        String hipoFile = "out_rga_fall2018_s6_9m.hipo"; mcsec=6;
///        String hipoFile = "fc-ecpcsc-elec-s1-eff-lh2.hipo"; mcsec=1;
        
        reader.open(hipoPath+hipoFile);
        
        engine.init();
        engine.setIsMC(true);
        engine.setVariation("default");
                
        setCCDB(10);
        setThresholds("Pizero",engine);
        setGeom("2.5");
        
        float dp1=0.1f, dp2=0.1f;
        
        H1F h00 = new H1F("h00",50,-0.1,0.1); h00.setTitleX("#Deltapx GEN-REC (GeV)");
        H1F h01 = new H1F("h01",50,-0.1,0.1); h01.setTitleX("#Deltapy GEN-REC (GeV)");
        H1F h02 = new H1F("h02",50,-0.1,0.1); h02.setTitleX("#Deltapz GEN-REC (GeV)");
        H1F h03 = new H1F("h03",50,-0.1,0.1); h03.setTitleX("#Deltap  GEN-REC (GeV)");
        
        H2F h10 = new H2F("h10",50,0,4.0,50,-0.1,1.0);  h10.setTitleX("ECAL #gamma E (GeV)"); h10.setTitleY("#Deltapx  (GeV)");
        H2F h11 = new H2F("h11",50,0,1.0,50,-0.1,0.1);  h11.setTitleX("ECAL #gamma E (GeV)"); h11.setTitleY("#Deltapy  (GeV)");
        H2F h12 = new H2F("h12",50,0,1.0,50,-0.1,1.0);  h12.setTitleX("ECAL #gamma E (GeV)"); h12.setTitleY("#Deltapz  (GeV)");
        H2F h13 = new H2F("h13",50,0,5.0,100,-0.1,5.0); h13.setTitleX("ECAL #gamma E (GeV)"); h13.setTitleY("#Deltap   (GeV)");
        
        H2F h20 = new H2F("h20",50,0,100.0,100,-1,5); h20.setTitleX("Distance e-#gamma (cm)");        h20.setTitleY("#Deltap  GEN-REC (GeV)");
        H2F h21 = new H2F("h21",50,0,100.0,100,0,5);  h21.setTitleX("Distance e-#gamma (cm)");        h21.setTitleY("ECAL #gamma E (GeV)");
        H2F h22 = new H2F("h23",50,0,100.0,100,0,5);  h22.setTitleX("Distance e-#gamma (cm)");        h22.setTitleY("ECAL #gamma E (GeV)"); 
        H2F h23 = new H2F("h23",50,0,9,50,0,5);       h23.setTitleX("Track Electron Momentum (GeV)"); h23.setTitleY("ECAL #gamma E (GeV)");
        h21.setTitle("#Deltap>"+dp1); h22.setTitle("#Deltap<"+dp2);  h23.setTitle("#Deltap<"+dp2); 
        
        String lab[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};
        H2F h30 = new H2F("h20",50,0,9,70,-2,12); h30.setTitleX("Track Electron Momentum (GeV)");h30.setTitleY("PCAL "+lab[0]);
        H2F h31 = new H2F("h21",50,0,9,70,-5,30); h31.setTitleX("Track Electron Momentum (GeV)");h31.setTitleY("PCAL "+lab[1]);
        H2F h32 = new H2F("h22",50,0,9,70,-2,12); h32.setTitleX("Track Electron Momentum (GeV)");h32.setTitleY("ECIN "+lab[0]);
        H2F h33 = new H2F("h23",50,0,9,70,-5,30); h33.setTitleX("Track Electron Momentum (GeV)");h33.setTitleY("ECIN "+lab[1]);
        h30.setTitle("#Deltap>"+dp1);h31.setTitle("#Deltap>"+dp1);h32.setTitle("#Deltap>"+dp1);h33.setTitle("#Deltap>"+dp1);
        H2F h40 = new H2F("h30",50,0,9,70,-2,12); h40.setTitleX("Track Electron Momentum (GeV)");h40.setTitleY("PCAL "+lab[0]);
        H2F h41 = new H2F("h31",50,0,9,70,-5,30); h41.setTitleX("Track Electron Momentum (GeV)");h41.setTitleY("PCAL "+lab[1]);
        H2F h42 = new H2F("h32",50,0,9,70,-2,12); h42.setTitleX("Track Electron Momentum (GeV)");h42.setTitleY("ECIN "+lab[0]);
        H2F h43 = new H2F("h33",50,0,9,70,-5,30); h43.setTitleX("Track Electron Momentum (GeV)");h43.setTitleY("ECIN "+lab[1]);
        h40.setTitle("#Deltap<"+dp2);h41.setTitle("#Deltap<"+dp2);h42.setTitle("#Deltap<"+dp2);h43.setTitle("#Deltap<"+dp2);
               
        H2F h50 = new H2F("h40",50,0.0,10.,50,0.15,0.31);      
        h50.setTitleX("True Electron Momentum (GeV)");
        h50.setTitleY("Sampling Fraction E/P");
        
        H2F h51 = new H2F("h41",50,0.0,10.,50,0.15,0.31);      
        h51.setTitleX("Track Electron Momentum (GeV)");
        h51.setTitleY("Sampling Fraction E/P");
        
        H2F h52 = new H2F("h42",50,0.0,2.5,50,0.15,0.31);      
        h52.setTitleX("ECAL Electron Energy (GeV)");
        h52.setTitleY("Sampling Fraction E/P");
      
        H2F eff1 = new H2F("eff1",30,0,10,30,5,30); eff1.setTitleX("True Electron Momentum (GeV)"); eff1.setTitleY("True Electron Theta (deg)");
        
        H2F eff1a = new H2F("eff1a",30,0,10,30,5,30);
        H2F eff1b = new H2F("eff1b",30,0,10,30,5,30); 
        
        H2F eff2a = new H2F("eff2a",30,0,10,30,5,30);
        H2F eff2b = new H2F("eff2b",30,0,10,30,5,30); 
        
        H2F eff3a = new H2F("eff3a",30,0,10,30,5,30);
        H2F eff3b = new H2F("eff3b",30,0,10,30,5,30);
            
        int nh=10;
        H1F effmom[] = new H1F[nh]; H1F effthe[] = new H1F[nh]; 
        for (int i=0;i<nh;i++) {effmom[i] = new H1F("effmom",50,0.,10.); effmom[i].setTitleX("MC Electron P (GeV)");}
        for (int i=0;i<nh;i++) {effthe[i] = new H1F("effthe",50,5.,35.); effthe[i].setTitleX("MC Electron #theta (deg)");}
        
        int nevent = 0, nev=0, nelec,nphot,npim;
        
        List<Particle> elec = null;
        List<Particle> phot = null;
        List<Particle>  pim = null;
        
        while(reader.hasEvent()){
            
        	nevent++;
            DataEvent event = reader.getNextEvent();
            
            if(readMC(event)) { 
           
            refE = pmc.get(0).e(); refP = pmc.get(0).p(); refTH = pmc.get(0).theta()*180/Math.PI;
            
            eff1.fill(refE, refTH);
            if(refTH>15)effmom[0].fill(refE);
            effthe[0].fill(refTH);    
            
        	ev.init(event);        	
        	ev.setEventNumber(10);
        	ev.requireOneElectron(false);
       	    ev.setElecTriggerSector(mcsec);
       	    		
            if(ev.procEvent(event)) {nev++;
        		elec = getPART(0.0, 11);      nelec = elec.size();
        		phot = getPART(0.0, 22);      nphot = phot.size();
        		pim  = getPART(0.0, 212);      npim = pim.size();
        		
        		int npart = nelec+nphot+npim;        		        		
        		float esum=0,dp=0,chi2pid=1000;
        		
        		if(npart>0) {if(refTH>15)effmom[4].fill(refE); effthe[4].fill(refTH); eff1a.fill(refE, refTH);}
        		if(npim==0) {if(refTH>15)effmom[5].fill(refE); effthe[5].fill(refTH);}
       		
        		for (Particle p : phot) {esum+=p.e();}
        		        		
        		if(nelec==1) {
            		chi2pid = (float) Math.abs(elec.get(0).getProperty("chi2pid"));
            		dp = (float)(pmc.get(0).p()-elec.get(0).p());
            		
            		if(npim==0&&chi2pid<3.5)  {
            			if(refTH>15)                      effmom[6].fill(refE); 
            			if(refTH>15 && Math.abs(dp)<dp2)  effmom[7].fill(refE); 
            			if(Math.abs(dp)<dp2)   eff3a.fill(refE,refTH);
            			effthe[6].fill(refTH); eff2a.fill(refE,refTH); 	                    
            		}
            		
        			if(esum<0.01) {
        			h00.fill(pmc.get(0).px()-elec.get(0).px());
               		h01.fill(pmc.get(0).py()-elec.get(0).py());
               		h02.fill(pmc.get(0).pz()-elec.get(0).pz());
               		h03.fill(dp);
        			}
        			
                    if(esum>0.01) {
        			h10.fill(esum,pmc.get(0).px()-elec.get(0).px());
               		h11.fill(esum,pmc.get(0).py()-elec.get(0).py());
               		h12.fill(esum,pmc.get(0).pz()-elec.get(0).pz());
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
//        		System.out.println(" "+nevent+" "+nev+" "+nelec+" "+nphot+" "+npim);
            }
            
            dropBanks(event);
            engine.processDataEvent(event);  
            processDataEvent(event);  
            getUnmatchedResponses();
            rPC = getPCResponses(mcsec);
            
            double ecalE = getEcalEnergy(mcsec);

//          Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.04 && energy>0.12;
//          Boolean trig2 = good_pcal && !good_ecal && part.epc>0.04;
//          Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.06 && energy>0.15;
//          Boolean trig2 = good_pcal && !good_ecal && part.epc>0.15;
            
//            System.out.println("energy,1,2,3 = "+energy+" "+part.epc+" "+part.eec1+" "+part.eec2);
            
            Boolean good_pcal = epc>0.00;               //VTP reported cluster
            Boolean good_ecal = (eec1+eec2)>0.001 ; //VTP reported cluster
            Boolean trig1 = good_pcal &&  good_ecal && epc>0.04 && ecalE>0.12;
            Boolean trig2 = good_pcal && !good_ecal && epc>0.12;
            		
//            trig1=true; trig2=true;
            
            if(rPC.size()>0)                   {if(refTH>15)effmom[1].fill(refE);effthe[1].fill(refTH);}
            if(rPC.size()==1 || rPC.size()==2) {if(refTH>15)effmom[2].fill(refE);effthe[2].fill(refTH);}
            if(rPC.size()==1)                  {if(refTH>15)effmom[3].fill(refE);effthe[3].fill(refTH);
                                               h50.fill(refP, ecalE/refP); 
                                               h52.fill(ecalE,ecalE/refP);}                       
            }
        }
        
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
        	   System.out.println(xs[n]+" "+ys[n]+" "+yse[n]);
                xs[n] = 1/Math.sqrt(sigGraph.getDataX(i)/ym[i]);  //sig(E)/E vs True Energy
       	        ys[n] = y*y/ym[i]/ym[i]; //sigma(E)/E = sigma(E/P)*(P/E)
       	       yse[n] = 2*ys[n]*Math.sqrt(Math.pow(ye/y,2)+Math.pow(yme[i]/ym[i],2));
        	   System.out.println(xs[n]+" "+ys[n]+" "+yse[n]);
        	   resGraph.addPoint(xs[n], ys[n], 0., yse[n]);
        	}
        }
        
        
        JFrame frame = new JFrame("Electron Reconstruction");
        frame.setSize(1700,1900);
        
        GStyle.getH1FAttributes().setMarkerSize(4);
        GStyle.getH1FAttributes().setMarkerStyle(1);
        GStyle.getH1FAttributes().setLineWidth(1);
        

        /*
        canvas.cd(3); canvas.getPad(3).getAxisY().setRange(0.0,0.028) ; canvas.draw(sigGraph);
        canvas.cd(4); canvas.getPad(4).getAxisX().setRange(0.1,0.51) ; 
                      canvas.getPad(4).getAxisY().setRange(0.0030,0.0075) ; canvas.draw(resGraph);
        canvas.cd(4); canvas.getPad(4).getAxisX().setRange(0.00,1.5) ; 
                      canvas.getPad(4).getAxisY().setRange(0.00,0.18) ; canvas.draw(resGraph);
      */
    
//        H2F ratio = H2F.divide(eff2, eff1);
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
        
        EmbeddedCanvas c = null; int ncd;
               
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

        c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(1)); c.divide(4,2); ncd=-1;           
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

        c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(2)); c.divide(4,3); ncd=-1;           
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
        
        
        /*
        canvas.cd(3); canvas.getPad(3).getAxisY().setRange(0.5,1.015);
        canvas.draw(hrat1); 
        canvas.draw(hrat1); canvas.draw(hrat2,"same"); canvas.draw(hrat3,"same");
        canvas.cd(4); canvas.getPad(4).getAxisY().setRange(0.5,1.015);
        canvas.draw(hrat2,"same"); 
        canvas.cd(5); canvas.getPad(5).getAxisY().setRange(0.5,1.015);
        canvas.draw(hrat3,"same");
        */
        frame.add(c);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
        /*
        H2_a_Hist.add(5, 0, 0, h1);
        H2_a_Hist.add(5, 0, 1, h2);
        H2_a_Hist.add(5, 0, 2, h3);
        String hipoFileName = "/Users/colesmith/test.hipo";
        HipoFile histofile = new HipoFile(hipoFileName);
        histofile.addToMap("H2_a_Hist", H2_a_Hist); 
        histofile.writeHipoFile(hipoFileName);     
        */   
    }
      
    public void neutronDemo(String[] args) {
    	
        HipoDataSource reader = new HipoDataSource();
        ECEngine       engine = new ECEngine();
        
        JFrame frame = new JFrame("Neutron Reconstruction");
        frame.setSize(800,800);
        
        EmbeddedCanvas c = new EmbeddedCanvas(); 
        
    	int trSEC=5, trPID=-211, mcSEC=2, mcPID= 2112;
    	
        setMCpid(mcPID);
        setMCsec(mcSEC);
        
        String evioPath = "/Users/colesmith/clas12/sim/neutron/";
        String evioFile = "fc-neut-80k-s2-r5424.hipo";
//        String evioFile = "out_pim_n.hipo";
        
        if (args.length == 0) { 
            reader.open(evioPath+evioFile);
        } else {
            String inputFile = args[0];
            reader.open(inputFile);
        } 
        
        H1F  h1 = new H1F("h1",50,0.,3);   
        H1F  h2 = new H1F("h2",50,0.,3);  
        H1F  h3 = new H1F("h3",50,0.,3);   
        H1F  h4 = new H1F("h4",50,0.,3);   
        H1F  h5 = new H1F("h5",10,0,10);          
        H2F  h6 = new H2F("h6",50,0,3,100,0,1000); 
        H2F  h7 = new H2F("h7",50,0,3,100,0,1000); 
        H2F  h8 = new H2F("h8",50,0,3,10,1,11);
        H2F  h9 = new H2F("h9",50,0,3,10,1,11);
        H1F h10 = new H1F("h10",50,0,2200);
        H1F h11 = new H1F("h11",50,0,2200);
            
        engine.init();
//        engine.setOccMax(100000);
        engine.setIsMC(true);
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 

        setCCDB(2);
        
        setThresholds("Test",engine);
        setGeom("2.5");
        
        int run = 0;
        
        while(reader.hasEvent()) {
            DataEvent event = reader.getNextEvent();
            run = getRunNumber(event); 
            
            if(!readMC(event)) continue; 
            
            engine.processDataEvent(event); processDataEvent(event);
            
       		h1.fill(refP);
       		
   			boolean n1=true, n2=true; float en=0; int np=0;
       		boolean goodpc=false, goodeci=false, goodeco=false;
       		
       		for (DetectorParticle neut : getNeutralPart()) {         		
       		if(neut.getSector(DetectorType.ECAL)==mcSEC) {np++;
       			for (DetectorResponse ds : neut.getDetectorResponses()) { 
      				if (ds.getDescriptor().getType()==DetectorType.ECAL) {      					
       					if(n1 && ds.getDescriptor().getLayer()>0) {h2.fill(refP); n1=false;} //PCAL+ECAL
       					if(n2 && ds.getDescriptor().getLayer()>1) {h3.fill(refP); n2=false;} //ECAL ONLY
       					en+=ds.getEnergy()*1e3;
       					if(ds.getDescriptor().getLayer()==1) goodpc=true;
       					if(ds.getDescriptor().getLayer()==4) goodeci=true;
       					if(ds.getDescriptor().getLayer()==7) goodeco=true;       					
        			}        	       		
      			}
       		}
       		}       			
            if(en>0)                     {h6.fill(refP,en);h8.fill(refP, np);h10.fill(1e3*(refE-0.93957));}
            if(goodpc&&goodeci&&goodeco) {h7.fill(refP,en);h9.fill(refP, np);h11.fill(1e3*(refE-0.93957));}
        }
        
        H1F h21 = H1F.divide(h2,h1); H1F h31 = H1F.divide(h3,h1);  

        h1.setTitle("GRN:GEN    BLK:REC    RED:ECAL ONLY"); h1.setTitleX("GEN P (GeV)"); h1.setTitleY("COUNTS");
        h1.setFillColor(3);h1.setLineColor(3); h2.setFillColor(1); h3.setFillColor(2); h3.setLineColor(2);        
        
        h21.setTitle("BLK:REC    RED:ECAL ONLY");h21.setTitleX("GEN P (GeV)"); h21.setTitleY("EFFICIENCY"); 
        h6.setTitle("ANY LAYER");h6.setTitleX("GEN P (GeV)");h6.setTitleY("REC E (MeV)");
        h7.setTitle("PCAL+ECIN+ECOU");h7.setTitleX("GEN P (GeV)");h7.setTitleY("REC E (MeV)");
        h8.setTitleX("GEN P (GeV)"); h8.setTitleY("REC NEUTRONS");
        h9.setTitleX("GEN P (GeV)"); h9.setTitleY("REC NEUTRONS");
        h10.setTitleX("GEN K.E. (MeV)"); h11.setTitleX("GEN K.E. (MeV)");
        GraphErrors g21 = h21.getGraph(); GraphErrors g31 = h31.getGraph();
        g21.setTitle(h21.getTitle());g21.setTitleX(h21.getTitleX());g21.setTitleY(h21.getTitleY());
        g21.setMarkerColor(1); g31.setMarkerColor(2); g31.setLineColor(2);

        c.divide(2,4);
        c.cd(0); c.draw(h1); c.draw(h2,"same"); c.draw(h3,"same"); 
        c.cd(1); c.getPad().getAxisY().setRange(0, 0.85); c.draw(g21); c.draw(g31,"same");
        c.cd(2); c.getPad().getAxisZ().setLog(true); c.draw(h6);
        c.cd(3); c.getPad().getAxisZ().setLog(true); c.draw(h7);
        c.cd(4); c.getPad().getAxisZ().setRange(1,1000); c.getPad().getAxisZ().setLog(true); c.draw(h8);
        c.cd(5); c.getPad().getAxisZ().setRange(1,1000); c.getPad().getAxisZ().setLog(true); c.draw(h9);
        c.cd(6); c.draw(h10);
        c.cd(7); c.draw(h11);
              
//        dumpGraph("/Users/colesmith/CLAS12ANA/ECpart/files/neutronDemo_r"+run+"_h21_nopc_noftof.vec",h21.getGraph());
//        dumpGraph("/Users/colesmith/CLAS12ANA/ECpart/files/neutronDemo_r"+run+"_h31_nopc_noftof.vec",h31.getGraph());

        frame.add(c);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);  
        System.out.println("FINISHED");
    }
    
    public void photonDist(String[] args) {
    	
        HipoDataSource reader = new HipoDataSource();
        ECEngine       engine = new ECEngine();
        List<DetectorParticle> np = new ArrayList<DetectorParticle>();
        int run = this.runNumber;    
        
        H2F h2a = new H2F("Dist Error 1",50,0,70,50,0.95,4);         
        h2a.setTitleX("Distance (cm)"); h2a.setTitleY("E1 / E2");
        
        H2F h2b = new H2F("Dist Error 2",50,0,70,50,-0.4,0.4);         
        h2b.setTitleX("Distance (cm)"); h2b.setTitleY("#Delta E / E");
        
        H2F h2c = new H2F("Dist Error 3",50,0,4,50,0,4);         
        h2c.setTitleX("E #gamma 1 / E #gamma 2 (PCAL)"); h2c.setTitleY("E #gamma 1 / E #gamma 2 (ECIN)");
        
        H2F h2cc = new H2F("Dist Error 3c",50,0,4,50,0,4);         
        h2cc.setTitleX("E #gamma 1 / E #gamma 2 (PCAL)"); h2cc.setTitleY("E #gamma 1 / E #gamma 2 (ECIN)");
        
        H2F h2d = new H2F("Dist Error 4",50,0,4,50,0,4);         
        h2d.setTitleX("E #gamma 1 / E #gamma 2 (PCAL)"); h2d.setTitleY("E #gamma 2 / E #gamma 1 (ECIN)");

        H2F h2e = new H2F("Dist Error 5",80,0,70,80,-1.5,1.5);         
        h2e.setTitleX("Distance (cm)"); h2e.setTitleY("#Delta#theta 12");
        
        H2F h2f = new H2F("Dist Error 6",50,0,5.2,50,-0.4,0.4);         
        h2f.setTitleX("#theta 12"); h2f.setTitleY("#Delta E / E");
        
        H2F h2g = new H2F("Dist Error 7",10,1,11,50,-1,1);
        h2g.setTitleX("Number of photons"); h2g.setTitleY("#Delta E / E");
        
        H1F h11 = new H1F("Eff1",50,0,5.2);
        H1F h12 = new H1F("Eff2",50,0,5.2);
        H1F h13 = new H1F("Eff3",50,0,5.2);
        H1F h14 = new H1F("Eff4",50,0,5.2);
        H1F h15 = new H1F("Eff5",50,0,5.2);
        H1F h16 = new H1F("Eff6",50,0,5.2);
        
        H2F hnpp = new H2F("npp",50,0,5.2,10,1,11);
                
        String evioPath = "/Users/colesmith/clas12/sim/2gamma/";
        String evioFile = "fc4gev-100k.hipo"; int sec=2; int etot=8;
        
        if (args.length == 0) { 
            reader.open(evioPath+evioFile);
        } else {
            String inputFile = args[0];
            reader.open(inputFile);
        } 
        
        HipoDataSync  writer = new HipoDataSync();
        writer.open("/Users/colesmith/photon_dist.hipo");
        
        engine.init();
        engine.setIsMC(true);
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        
        setCCDB(2);
        setThresholds("Test",engine);
        setGeom("2.5");
        setGoodPhotons(12);
        
        int n=0,np2=0,npc=0,nec=0;
        List<DetectorResponse> rPC  = new ArrayList<>();        
        List<DetectorResponse> rECi  = new ArrayList<>();        
        List<DetectorResponse> rECo  = new ArrayList<>();        
        List<DetectorResponse> rrPC  = new ArrayList<>();        
        List<DetectorResponse> rrECi1  = new ArrayList<>();        
        List<DetectorResponse> rrECi2  = new ArrayList<>();        
        List<DetectorResponse> rrECo  = new ArrayList<>();     
        
        int npp=0,indx=-1;
        float pthresh=1f;
       
        while(reader.hasEvent()){
            DataEvent de = reader.getNextEvent();
            engine.processDataEvent(de);
            run = getRunNumber(de);            
            if (readMC(de)) {
            	processDataEvent(de);
            	np.clear(); np = getNeutralPart(); 
            	getRECBanks(de,eb);
            	writer.writeEvent(de);
        		float      opa =(float) Math.toDegrees(Math.acos(pmc.get(0).cosTheta(pmc.get(1))));
            	h11.fill(opa);
            	npp=0;indx=-1;
            	if (np.size()>0) {
            		h12.fill(opa);
            		for (DetectorParticle phot : np) {
            			indx++;
                		double ep = phot.getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, phot, eb.ccdb);
            			if(ep>pthresh) {
//            				System.out.println(npp+" "+indx+" "+ep);
            				npp++;           				
            				for (DetectorResponse dr : phot.getDetectorResponses()) {
            					
            				}
            			}
            		}
            	}
            	hnpp.fill(opa,npp);
            	if (npp>1)  h13.fill(opa);
            	if (npp==2) h14.fill(opa);
            	if (npp>=2) h15.fill(opa);
            	if (npp>=2) {
                	double e1=0,e2=0,the1=0,the2=0,dist=0;
                    p1 = np.get(0);  //Photon 1
                    p2 = np.get(1);  //Photon 2   
                    Vector3 n1 = p1.vector(); n1.unit();
                    Vector3 n2 = p2.vector(); n2.unit();
                    e1=np.get(0).getEnergy(DetectorType.ECAL);
                    e2=np.get(1).getEnergy(DetectorType.ECAL);
                    SF1db = SamplingFractions.getMean(22, np.get(0), eb.ccdb);
                    e1 = e1/SF1db;
                    SF1db = SamplingFractions.getMean(22, np.get(1), eb.ccdb);
                    e2 = e2/SF1db;
                    Particle g1 = new Particle(22,n1.x()*e1,n1.y()*e1,n1.z()*e1);        
                    Particle g2 = new Particle(22,n2.x()*e2,n2.y()*e2,n2.z()*e2); 
                    the1 = Math.acos(g1.cosTheta(g2))*180/3.14159; 
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
                    h2f.fill(opa, (e1+e2)/etot-1);
                    if(npc==2) {
                    	r2.sub(r1); dist=r2.mag();
                    	h2a.fill(dist, e2>0?e1/e2:1000);
                    	h2b.fill(dist, (e1+e2)/etot-1);
                    	h2e.fill(dist,opa-the1);
                    	if(nec==2) h2c.fill(epc1/epc2,eeci1/eeci2);
                    	if(nec==2 && e1/e2>1.3) h2cc.fill(epc1/epc2,eeci1/eeci2);
                    	h2g.fill(npp,(e1+e2)/etot-1);
            	    }
            	}           		
            }
        }
        
        writer.close();
            	
            	
/*
                n++; 
            	if(np.size()==2)   np2++; 
            	double e1=0,e2=0,the1=0,the2=0,dist=0;
//            	System.out.println(n+" "+np.size()+" "+refE+" "+" "+rPC.size()+" "+rECi.size()+" "+rECo.size());
            
            	
            	double esum=0; int npp=0;
            	for (DetectorParticle p : np) {
            		double ep = p.getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p, eb.ccdb);
            		if(ep>1.4) {esum+= ep; npp++;}
            	}
            	h2g.fill(npp,esum/etot-1);
            	
            	if (npp>0)  h12.fill(opa);
            	if (npp>1)  h13.fill(opa);
            	if (npp>=2) h14.fill(opa);
            	if (npp==2) h16.fill(opa);
            	if (false) {
            		h15.fill(opa);
                	p1 = np.get(0);  //Photon 1
                	p2 = np.get(1);  //Photon 2   
                	Vector3 n1 = p1.vector(); n1.unit();
                	Vector3 n2 = p2.vector(); n2.unit();
                	e1=np.get(0).getEnergy(DetectorType.ECAL);
                	e2=np.get(1).getEnergy(DetectorType.ECAL);
                	SF1db = SamplingFractions.getMean(22, np.get(0), eb.ccdb);
                	e1 = e1/SF1db;
                	SF1db = SamplingFractions.getMean(22, np.get(1), eb.ccdb);
                	e2 = e2/SF1db;
                    Particle g1 = new Particle(22,n1.x()*e1,n1.y()*e1,n1.z()*e1);        
                    Particle g2 = new Particle(22,n2.x()*e2,n2.y()*e2,n2.z()*e2);   
                    the1 = Math.acos(g1.cosTheta(g2))*180/3.14159;      
                    
                	Vector3D r1 = rPC.get(0).getPosition();
                	Vector3D r2 = rPC.get(1).getPosition();
                	r2.sub(r1); dist=r2.mag();  
//                	float epc1  = (float) rPC.get(0).getEnergy();
//                	float epc2  = (float) rPC.get(1).getEnergy();
//                	float eeci1 = (float) rECi.get(0).getEnergy();
//                	float eeci2 = (float) rECi.get(1).getEnergy();
//                	float eeco1 = (float) rPC.get(0).getEnergy();
//                	float eeco2 = (float) rPC.get(1).getEnergy();
                	
                	double epc1 = DetectorResponse.getListByLayer(np.get(0).getDetectorResponses(),DetectorType.ECAL, 1).get(0).getEnergy();
                	double epc2 = DetectorResponse.getListByLayer(np.get(1).getDetectorResponses(),DetectorType.ECAL, 1).get(0).getEnergy();
//                	double eeci1 = DetectorResponse.getListByLayer(np.get(0).getDetectorResponses(),DetectorType.ECAL, 4).get(0).getEnergy();
//                	double eeci2 = DetectorResponse.getListByLayer(np.get(1).getDetectorResponses(),DetectorType.ECAL, 4).get(0).getEnergy();
                	rrECi1 = DetectorResponse.getListByLayer(np.get(0).getDetectorResponses(),DetectorType.ECAL, 4);
                	rrECi2 = DetectorResponse.getListByLayer(np.get(1).getDetectorResponses(),DetectorType.ECAL, 4);
                	double eeci1 = rrECi1.size()>0 ? rrECi1.get(0).getEnergy():0;
                	double eeci2 = rrECi2.size()>0 ? rrECi2.get(0).getEnergy():0;
                	h2a.fill(dist, e2>0?e1/e2:1000);
                	h2b.fill(dist, (e1+e2)/etot-1); h2f.fill(opa, (e1+e2)/etot-1);
                	if(eeci2>0 && epc2>0) h2c.fill(epc1/epc2,eeci1/eeci2);
                	if(eeci1>0 && epc2>0) h2d.fill(epc1/epc2,eeci2/eeci1);
                	h2e.fill(dist,opa-the1);
                	
            	}
            	*/
//            	System.out.println(dist+" "+the1+" "+opa+" "+(e1+e2)+" "+(e2>0?e1/e2:0));


        System.out.println("Total events= "+n+" np2= "+np2+" npc= "+npc);
        
        JFrame frame = new JFrame("Photon Efficieny");
        frame.setSize(800,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas(); canvas.divide(4,4);
        
        canvas.cd(0); canvas.draw(h2a);
        canvas.cd(1); canvas.draw(h2b);
        canvas.cd(2); canvas.draw(h2f);
        canvas.cd(3); canvas.draw(h2e);
        
        canvas.cd(4); canvas.draw(h2c);
        canvas.cd(5); canvas.draw(h2cc);
         
        canvas.cd(6); canvas.draw(hnpp);
        canvas.cd(7); canvas.draw(h2g);
        canvas.cd(8); canvas.draw(h2g.projectionY());
        canvas.cd(9); canvas.draw(h2g.projectionX());
        
        H1F hr12 = H1F.divide(h12,  h11); hr12.setFillColor(4); hr12.setTitleY("Eff n>0");   hr12.setTitleX("Opening Angle (deg)");
        H1F hr13 = H1F.divide(h13,  h11); hr13.setFillColor(3); hr13.setTitleY("Eff n>1");   hr13.setTitleX("Opening Angle (deg)");
        H1F hr14 = H1F.divide(h14,  h11); hr14.setFillColor(2); hr14.setTitleY("Eff n==2");  hr14.setTitleX("Opening Angle (deg)");
        H1F hr15 = H1F.divide(h15,  h11); hr15.setFillColor(5); hr15.setTitleY("Eff n>=2");  hr15.setTitleX("Opening Angle (deg)");
//        H1F hr16 = H1F.divide(h16,  h11); hr16.setFillColor(6); hr16.setTitleY("EFF n==2");  hr16.setTitleX("Opening Angle (deg)");
        canvas.cd(10);  canvas.draw(h11); canvas.draw(h12,"same"); 
        canvas.cd(11);  canvas.draw(hr13); 
        canvas.cd(12); canvas.draw(hr14); 
        canvas.cd(13); canvas.draw(hr15); 
//        canvas.cd(12); canvas.draw(hr16); 
        
        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);          
    }
    
    public void photonDemo(String[] args) {
    	
        HipoDataSource reader = new HipoDataSource();
        ECEngine       engine = new ECEngine();
        List<DetectorParticle> np = new ArrayList<DetectorParticle>();
        int run = this.runNumber;
        
        String evioPath = "/Users/colesmith/clas12/sim/photon/";
        String evioFile = "fc-phot-20k-25deg-r2-0.2-3.8-newgeom-gemc2.6.hipo"; int sec=2;
        
        if (args.length == 0) { 
            reader.open(evioPath+evioFile);
        } else {
            String inputFile = args[0];
            reader.open(inputFile);
        } 
        
        h5 = new H1F("Thrown",50,0.,3.5); h5.setTitleX("MC Photon E (MeV)");
        
        H1F  h11 = new H1F("n>0 any layer",  50,0.,3.5); h11.setTitleX("MC Photon E (MeV)"); 
        H1F  h12 = new H1F("n=1 any layer",  50,0.,3.5); h12.setTitleX("MC Photon E (MeV)");  
        H1F  h13 = new H1F("n=1 layer 1",    50,0.,3.5); h13.setTitleX("MC Photon E (MeV)");
        H1F  h14 = new H1F("n=1 layer 1,4",  50,0.,3.5); h14.setTitleX("MC Photon E (MeV)");
        H1F  h15 = new H1F("n=1 layer 1,4,7",50,0.,3.5); h15.setTitleX("MC Photon E (MeV)");
        
        H1F  h21  = new H1F("ndist21",50,0.,3.5);  h21.setTitleX("MC Photon E (MeV)");  
        H1F  h22  = new H1F("ndist22",50,0.,3.5);  h22.setTitleX("MC Photon E (MeV)");  
        H1F  h23  = new H1F("ndist23",50,0.,3.5);  h23.setTitleX("MC Photon E (MeV)");  
        H1F  h210 = new H1F("ndist210",50,0.,3.5); h210.setTitleX("MC Photon E (MeV)");  
        H1F  h220 = new H1F("ndist220",50,0.,3.5); h220.setTitleX("MC Photon E (MeV)");  
        H1F  h230 = new H1F("ndist230",50,0.,3.5); h230.setTitleX("MC Photon E (MeV)"); 
        H1F  h211 = new H1F("nmult211",5,1,6);     h211.setTitleX("PCAL Clusters");  
        H1F  h221 = new H1F("nmult221",5,1,6);     h221.setTitleX("ECIN Clusters");  
        H1F  h231 = new H1F("nmult231",5,1,6);     h231.setTitleX("ECOU Clusters"); 
        H1F  h212 = new H1F("nmult212",50,0.,1);   h212.setTitleX("PCAL Energy Fraction");  
        H1F  h222 = new H1F("nmult222",50,0.,1);   h222.setTitleX("ECIN Energy Fraction");  
        H1F  h232 = new H1F("nmult232",50,0.,1);   h232.setTitleX("ECOU Energy Fraction"); 
      
        engine.init();
        engine.setIsMC(true);
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        
        setCCDB(2);
        setThresholds("Photon",engine);
        setGeom("2.5");
        setGoodPhotons(12);
      
        int n=0, ng=0, npp=0;
        float pthresh=0.01f;
        
        while(reader.hasEvent()){
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event);
            if (readMC(event)) {
            	processDataEvent(event);
            	float refP = (float) pmc.get(0).p();
            	np.clear(); np = getNeutralPart();
            	h5.fill(refP);
            	npp=0;
            	if(np.size()>0) { // 1 or more neutral particles id=22,2112
            		h11.fill(refP); n++;
        	    	int n1=0,n4=0,n7=0,sum=0;
        	    	double e1=0, e4=0, e7=0, e1p=0, e4p=0, e7p=0;
            	    for (DetectorParticle phot : np) {
            	    	double p1 = phot.getEnergy(DetectorType.ECAL);
            	    	if(p1>pthresh) {
            	    	npp++;
//            	    	System.out.println(refP+" "+p1/SamplingFractions.getMean(22, phot, eb.ccdb));
            	    	for (DetectorResponse dr : phot.getDetectorResponses()) {
            	    		int lay = dr.getDescriptor().getLayer();
            	    		double e = dr.getEnergy();
            	    		if(lay==1) {n1++; e1p+=e; h211.fill(n1); if(n1==1) {e1=e; h210.fill(refP);}}
            	    		if(lay==4) {n4++; e4p+=e; h221.fill(n4); if(n4==1) {e4=e; h220.fill(refP);}}
            	    		if(lay==7) {n7++; e7p+=e; h231.fill(n7); if(n7==1) {e7=e; h230.fill(refP);}}
            	    	}
            	    	}
            	    }
            	    h21.fill(refP,n1); h22.fill(refP,n4); h23.fill(refP,n7);
                    if(n1>1) h212.fill(e1/e1p);
                    if(n4>1) h222.fill(e4/e4p);
                    if(n7>1) h232.fill(e7/e7p);
            	}
            	if(npp==1) { // only 1 neutral particle id=22,2112            		
            		h12.fill(refP); ng++;
            	    for (DetectorParticle phot : np) {
            	    	double p1 = phot.getEnergy(DetectorType.ECAL);
            	    	int sum=0;
 //           	    	if(p1>pthresh) {
            	    	for (DetectorResponse dr : phot.getDetectorResponses()) {
            	    		int lay = dr.getDescriptor().getLayer();
            	    		if(lay==1) sum+= 100;
            	    		if(lay==4) sum+=  40;
            	    		if(lay==7) sum+=   7;
            	    	}
            	    	if(sum>=100) h13.fill(refP);
            	    	if(sum>=140) h14.fill(refP);
            	    	if(sum==147) h15.fill(refP);
//            	    	}
            		}
            	}
            }
        }
        System.out.println("Total events= "+n+" ng= "+ng);

        JFrame frame = new JFrame("Photon Reconstruction");
        frame.setSize(800,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas(); canvas.divide(4,3);
        
        
        
        H1F hr11 = H1F.divide(h11,  h5); hr11.setFillColor(3); hr11.setTitleY("Photon Efficiency n>0"); hr11.setTitleX("Photon Momentum (GeV)");
        H1F hr12 = H1F.divide(h12,  h5); hr12.setFillColor(1); hr12.setTitleY("Photon Efficiency n=1"); hr12.setTitleX("Photon Momentum (GeV)");
        H1F hr13 = H1F.divide(h13,  h5); hr13.setFillColor(2); hr13.setTitleY("Photon Eff 1");          hr13.setTitleX("Photon Momentum (GeV)");
        H1F hr14 = H1F.divide(h14,  h5); hr14.setFillColor(5); hr14.setTitleY("Photon Eff 14");         hr14.setTitleX("Photon Momentum (GeV)");
        H1F hr15 = H1F.divide(h15,  h5); hr15.setFillColor(4); hr15.setTitleY("Photon Eff 147");        hr15.setTitleX("Photon Momentum (GeV)");
        canvas.cd(4);  
        canvas.draw(hr11); canvas.draw(hr12,"same"); canvas.draw(hr13,"same"); canvas.draw(hr14,"same"); canvas.draw(hr15,"same");
        
        H1F hr21 = H1F.divide(h21,  h210); hr21.setFillColor(2); hr21.setTitleY("Avg. No. PCAL Clusters"); hr21.setTitleX("Photon Momentum (GeV)");
        H1F hr22 = H1F.divide(h22,  h220); hr22.setFillColor(5); hr22.setTitleY("Avg. No. ECIN Clusters"); hr22.setTitleX("Photon Momentum (GeV)");
        H1F hr23 = H1F.divide(h23,  h230); hr23.setFillColor(4); hr23.setTitleY("Avg. No. ECOU Clusters"); hr23.setTitleX("Photon Momentum (GeV)");
        canvas.cd(5); canvas.getPad().getAxisY().setRange(1, 1.25); canvas.draw(hr21); 
        canvas.cd(6); canvas.getPad().getAxisY().setRange(1, 1.25); canvas.draw(hr22); 
        canvas.cd(7); canvas.getPad().getAxisY().setRange(1, 1.25); canvas.draw(hr23);
        
        H1F hr210 = H1F.divide(h210,  h5); hr210.setFillColor(2); hr210.setTitleY("Photon Eff 1");   hr210.setTitleX("Photon Momentum (GeV)");
        H1F hr220 = H1F.divide(h220,  h5); hr220.setFillColor(5); hr220.setTitleY("Photon Eff 14");  hr220.setTitleX("Photon Momentum (GeV)");
        H1F hr230 = H1F.divide(h230,  h5); hr230.setFillColor(4); hr230.setTitleY("Photon Eff 147"); hr230.setTitleX("Photon Momentum (GeV)");
        canvas.cd(0); canvas.draw(hr11); canvas.draw(hr210,"same"); canvas.draw(hr220,"same"); canvas.draw(hr230,"same");
        canvas.cd(1); canvas.getPad().getAxisY().setLog(true); h211.setFillColor(2); canvas.draw(h211);
        canvas.cd(2); canvas.getPad().getAxisY().setLog(true); h221.setFillColor(5); canvas.draw(h221);
        canvas.cd(3); canvas.getPad().getAxisY().setLog(true); h231.setFillColor(4); canvas.draw(h231);

        canvas.cd(9);  h212.setFillColor(2); canvas.draw(h212);
        canvas.cd(10); h222.setFillColor(5); canvas.draw(h222);
        canvas.cd(11); h232.setFillColor(4); canvas.draw(h232);      
  
        System.out.println("Done");
        
//        dumpGraph("/Users/colesmith/CLAS12ANA/ECpart/files/photeff_r"+run+".vec",hrat1.getGraph());
//        dumpGraph("/Users/colesmith/CLAS12ANA/ECpart/files/photeff_r"+run+".vec",hrat2.getGraph());
//        dumpGraph("/Users/colesmith/CLAS12ANA/ECpart/files/photeff_r"+run+".vec",hrat3.getGraph());
//        dumpGraph("/Users/colesmith/CLAS12ANA/ECpart/files/photeff_r"+run+".vec",hrat4.getGraph());

        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);  
    }

    public void pizeroDemo(String[] args) {
    	
        HipoDataSource reader = new HipoDataSource();
//     	HipoReader reader = new HipoReader();
        ECEngine       engine = new ECEngine();
        List<DetectorParticle> np = new ArrayList<DetectorParticle>();
        
        H2F h2a,h2b,h2c,h2d,h2e,h2f;
        H1F h5a,h5b,h5c,h6,h7a,h7b,h8;
        H2F hmp;
        
        int n2hit=0,n2rec1=0,n2rec2=0,nimcut=0;        
        double emax = 12.5;
        
        h2a = new H2F("Invariant Mass",50,0.,emax,50,100.,200);         
        h2a.setTitleX("Pizero Energy (GeV)"); h2a.setTitleY("Two-Photon Invariant Mass (MeV)");
        
        h2b = new H2F("Energy Asymmetry",50,0.,emax,50,-1.0,1.0);      
        h2b.setTitleX("Pizero Energy (GeV)"); h2b.setTitleY("X:(E1-E2)/(E1+E2)");
       
        h2c = new H2F("Pizero Energy Error",50,0.,emax,50,-1,1); 
        h2c.setTitleX("Pizero Energy (GeV)"); h2c.setTitleY("Pizero Energy Error (GeV)");
       
        h2d = new H2F("Pizero Theta Error",50,0.,emax,50,-1.,1.);      
        h2d.setTitleX("Pizero Energy (GeV)") ; h2d.setTitleY("Pizero Theta Error (deg)");
        
        h2e = new H2F("Pizero Energy Error",50,0.,10,50,-1,1.);      
        h2e.setTitleX("Opening Angle (deg)") ; h2e.setTitleY("Pizero Energy Error (GeV)");
        
        h2f = new H2F("Pizero Energy Error",50,0.,10,50,-15.,15.);      
        h2f.setTitleX("Opening Angle (deg)") ; h2f.setTitleY("Invariant Mass Error (MeV)");
        
        h5 = new H1F("Thrown",50,0.,emax); h5.setTitleX("MC Pizero E (MeV)");
        
        hmp = new H2F("Mult",50,0,12,10,1,11);
        
//        h5a = new H1F("2Gamma",50,0.,20.);  h5a.setTitleX("Opening Angle (deg)");
//        h5b = new H1F("2Gamma",50,0.,20.);  h5b.setTitleX("Opening Angle (deg)");
//        h5c = new H1F("2Gamma",50,0.,20.);  h5c.setTitleX("Opening Angle (deg)");
        
        h6 = new H1F("2Gamma",50,0.,emax);  h6.setTitleX("MC Pizero E (MeV)"); 
        h7a = new H1F("PC/EC", 50,0.,emax); h7a.setTitleX("MC Pizero E (MeV)");
        h7b = new H1F("PC/EC", 50,0.,emax); h7b.setTitleX("MC Pizero E (MeV)");
        h8 = new H1F("Mcut",  50,0.,emax);  h8.setTitleY("ECIN 2 Photon Eff"); h8.setTitle("80<InvMass<200");
        H2F h9 = new H2F("Beta1", 50,0.75, 1.1,50,650.,750.); H2F h10 = new H2F("Beta2",50,0.75,1.1,50,650.,750.);
        
        H1F[] hview = new H1F[4];
        hview[0] = new H1F("Shared View U",50,100,200);         
        hview[0].setTitleX("Invariant Mass (MeV)");
        hview[1] = new H1F("Shared View V",50,100,200);         
        hview[1].setTitleX("Invariant Mass (MeV)");
        hview[2] = new H1F("Shared View W",50,100,200);         
        hview[2].setTitleX("Invariant Mass (MeV)");
        hview[3] = new H1F("Shared View UVW",50,100,200);         
        hview[3].setTitleX("Invariant Mass (MeV)");
                
        String evioPath = "/Users/colesmith/clas12/sim/2gamma/aug.1.2021/";
//        String evioPath = "/Users/colesmith/clas12/sim/";
        
//        String evioFile = "fc-pizero-50k-s2-newgeom-0.35-8.35.hipo"; int sec=2;
//      String evioFile = "fc-pizero-50k-s2-newgeom-15-0.35-8.35.hipo4"; int sec=2;
//    String evioFile = "fc-pizero-50k-s2-newgeom-15-0.35-8.35-r5716.hipo"; int sec=2;
//        String evioFile = "fc-pizero-50k-s2-newgeom-0.35-8.35-r5716.hipo"; int sec=2;
        
//        String evioFile = "fc-pizero-100k-s2-newgeom-15-1.0-12.0.hipo"; int sec=2;
        String evioFile = "out-pim-pi0.hipo"; int sec=2;
//        String evioFile = "rga_fall2018_s2-pizero.hipo"; int sec=2;
//        String evioFile = "fc-pizero-s2-new-4.4.0.hipo"; int sec=2;
//        evioFile = "pi0_hi.hipo";
        
//         evioPath = "/Users/colesmith/ECMON/HIPO4/fbossu/3deg/";
//         evioFile = "out_6b.4.0_photons.hipo";sec=1;
         isMC = true;
        
        // GEMC file: 10k 2.0 GeV pizeros thrown at 25 deg into Sector 2 using GEMC 2.4 geometry
        // JLAB: evioPath = "/lustre/expphy/work/hallb/clas12/lcsmith/clas12/forcar/gemc/evio/";
         
        
        if (args.length == 0) { 
            reader.open(evioPath+evioFile);
        } else {
            String inputFile = args[0];
            reader.open(inputFile);
        }
     
        HipoDataSync  writer = new HipoDataSync();
        writer.open("/Users/colesmith/pizero_demo.hipo");
             
//        engine.setGeomVariation("rga_fall2018"); //this must be set before engine.init()
        engine.init();
        engine.setIsMC(true);
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        engine.setVeff(16); //GEMC default
        engine.setPCTrackingPlane(9);
        
        setCCDB(10);
        setThresholds("Pizero",engine);
        setGeom("2.5");
        setGoodPhotons(12);
                
        while(reader.hasEvent()){
            DataEvent de = reader.getNextEvent();
            dropBanks(de);
            if (engine.processDataEvent(de)) {   
//            int iview = engine.getClusters().get(0).getStatus();
            
            if(readMC(de)) {  
            	processDataEvent(de);
            	getNeutralResponses();
            	getRECBanks(de,eb);
            	writer.writeEvent(de);
 //           	np.clear(); np = getNeutralPart(); hmp.fill(refE,np.size());            	
            	double   invmass = 1e3*Math.sqrt(getTwoPhotonInvMass(sec));            
            	boolean goodmass = invmass>0 && invmass<200;            
            	boolean    pcec1 = goodPhotons(12,p1,p2);
            	boolean    pcec2 = goodPhotons(12,p1,p2);
           
//            h9.fill(part.p1.getBeta(DetectorType.ECAL,1,0.),part.p1.getHit(DetectorType.ECAL).getPosition().z());  
//            h10.fill(part.p2.getBeta(DetectorType.ECAL,1,0.),part.p2.getHit(DetectorType.ECAL).getPosition().z());

            	h5.fill(refE);
            	if (goodmass) {
                                      n2hit++;   h6.fill(refE);  	            
                  if(pcec1 || pcec2) {n2rec1++; h7a.fill(refE);}
                  if(pcec1 && pcec2) {n2rec2++; h7b.fill(refE);}
          
                  if (pcec1&&pcec2) {
//                	  if(engine.hasSharedView()) {hview[iview].fill(invmass);hview[3].fill(invmass);}
                	  if(true) {
                	  h2a.fill(refE, invmass);                                        //Two-photon invariant mass                
                	  h2b.fill(refE, X);                                     	      //Pizero energy asymmetry
                	  h2c.fill(refE,(Math.sqrt(tpi2)-refE));            			  //Pizero total energy error
                	  h2e.fill(Math.acos(cth)*180/Math.PI,(Math.sqrt(tpi2)-refE));    //Pizero total energy error
                	  h2f.fill(Math.acos(cth)*180/Math.PI,invmass-mpi0*1e3);          //Pizero total energy error
                	  h2d.fill(refE,Math.acos(cpi0)*180/Math.PI-refTH); 			  //Pizero theta angle error
                	  nimcut++; h8.fill(refE);
                	  }
                  }
            	}
            }
            }
        }
        
        writer.close();
        
        H1F hrat1 = H1F.divide(h6,  h5); hrat1.setFillColor(2); hrat1.setTitleY("PC 2 Photon Eff");    hrat1.setTitleX("Pizero Energy (GeV)");
        H1F hrat2 = H1F.divide(h7a, h5); hrat2.setFillColor(2); hrat2.setTitleY("PC*EC 1 Photon Eff"); hrat2.setTitleX("Pizero Energy (GeV)");
        H1F hrat3 = H1F.divide(h7b, h5); hrat3.setFillColor(2); hrat3.setTitleY("PC*EC 2 Photon Eff"); hrat3.setTitleX("Pizero Energy (GeV)");
        H1F hrat4 = H1F.divide(h8, h6);  hrat4.setFillColor(2); hrat4.setTitleY("PC*EC NIMCUT Eff");   hrat4.setTitleX("Pizero Energy (GeV)");
        
        H1F h1 = h2a.projectionY();  h1.setOptStat("1100") ; h1.setFillColor(4); h1.setTitleX("Two-Photon Invariant Mass (MeV)");
        H1F h2 = h2b.projectionY();  h2.setOptStat("1100") ; h2.setFillColor(4); h2.setTitleX("X:(E1-E2)/(E1+E2)");
        H1F h3 = h2c.projectionY();  h3.setOptStat("1100") ; h3.setFillColor(4); h3.setTitleX("Pizero Energy Error");
        H1F h4 = h2d.projectionY();  h4.setOptStat("1100") ; h4.setFillColor(4); h4.setTitleX("Pizero Theta Error (deg)");
       
        System.out.println("THROWN TWOPHOTONS PCEC MATCH1 PCEC MATCH2 INVMASS CUT");
        System.out.println(n2mc+"     "+n2hit+"         "+n2rec1+"     "+n2rec2+"   "+nimcut);
        System.out.println("Eff1 = "+(float)n2hit/(float)n2mc+
                          " Eff2 = "+(float)n2rec1/(float)n2mc+
                          " Eff3 = "+(float)n2rec2/(float)n2mc+
                          " Eff4 = "+(float)nimcut/(float)n2hit);
        
        JFrame frame = new JFrame("Pizero Reconstruction");
        frame.setSize(800,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        canvas.divide(4,6);
		canvas.setAxisTitleSize(18);
		canvas.setAxisLabelSize(18);
		
        canvas.cd(0); canvas.draw(h2a);
        canvas.cd(1); canvas.draw(h2b);
        canvas.cd(2); canvas.draw(h2c);
        canvas.cd(3); canvas.draw(h2d);
        
        canvas.cd(4); canvas.draw(h1);
        canvas.cd(5); canvas.draw(h2);
        canvas.cd(6); canvas.draw(h3);
        canvas.cd(7); canvas.draw(h4);
        
        canvas.cd(8);  canvas.draw(hrat1); 
        canvas.cd(9);  canvas.draw(hrat2);
        canvas.cd(10); canvas.draw(hrat3);
        canvas.cd(11); canvas.draw(hrat4);
        
        canvas.cd(12); canvas.draw(h9);
        canvas.cd(13); canvas.draw(h10);
        canvas.cd(14); canvas.draw(h2e);
        canvas.cd(15); canvas.draw(h2f);
        
        canvas.cd(16); hview[0].setOptStat("1100"); canvas.draw(hview[0]);
        canvas.cd(17); hview[1].setOptStat("1100"); canvas.draw(hview[1]);
        canvas.cd(18); hview[2].setOptStat("1100"); canvas.draw(hview[2]);
        canvas.cd(19); hview[3].setOptStat("1100"); canvas.draw(hview[3]);
        
        canvas.cd(20); canvas.draw(hmp);
        
	    ParallelSliceFitter fitter = new ParallelSliceFitter(h2a);
        fitter.fitSlicesX();  
        dumpGraph("/Users/colesmith/pi0fit.vec",fitter.getSigmaSlices());
        
        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);     
    }
    
    public void scalerdemo(String[] args) {
        HipoDataSource reader = new HipoDataSource();
        ECEngine       engine = new ECEngine();    	
        String evioPath = "/Users/colesmith/clas12/sim/neutron/";
        String evioFile = "junk.hipo"; 
     
        if (args.length == 0) { 
        	reader.open(evioPath+evioFile);
        } else {
        	String inputFile = args[0];
        	reader.open(inputFile);
        } 
        H2F[] ha = new H2F[3];
        H2F[] ht = new H2F[3];
        
        ha[0] = new H2F("pcsca", 9,1,10,1152,1,1153);
        ha[1] = new H2F("ecisca",9,1,10,648,1,649);
        ha[2] = new H2F("ecosca",9,1,10,648,1,649);
        ht[0] = new H2F("pcsct", 9,1,10,1152,1,1153);
        ht[1] = new H2F("ecisct",9,1,10,648,1,649);
        ht[2] = new H2F("ecosct",9,1,10,648,1,649);
        
        engine.init();
//        engine.setOccMax(10000);
        engine.setIsMC(true);
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 

        setCCDB(2);
        
        int nev=1;
        
        while(reader.hasEvent()) {
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event); 
            
            if(event.hasBank("ECAL::scaler")) {
            	DataBank  bank = event.getBank("ECAL::scaler");
            	System.out.println(bank.rows());
            	for(int i=0; i<bank.rows(); i++) {
            		int is = bank.getByte("sector", i);
            		int il = bank.getByte("layer", i);
            		int ip = bank.getShort("component", i);
            		int ie = bank.getInt("event", i);
            		int ic = bank.getInt("acount", i);
            		int it = bank.getInt("tcount", i);
            		if(il>0) {
            			ha[getLayer(il)].fill(nev, getIndex(is,il,ip),ic);
            			ht[getLayer(il)].fill(nev, getIndex(is,il,ip),it);
            		}
            	}
                nev++;
           }
            
        }
        
        reader.close();
        
        JFrame frame = new JFrame("Scaler Bank Demo");
        frame.setSize(800,800);
        
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        
        canvas.divide(3,2);
        for (int i=0; i<3; i++) {canvas.cd(i);   canvas.draw(ha[i]);}
        for (int i=0; i<3; i++) {canvas.cd(i+3); canvas.draw(ht[i]);}
        
        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);          
    }
    
    public int getLayer(int il) {
    	int[] lay = {0,0,0,1,1,1,2,2,2};
    	return lay[il-1];
    }
    
    public int getIndex(int is, int il, int ip) {  
        int off[] = {0,68,130,0,36,72,0,36,72};
        int sca[] = {192,192,192,216,216,216,216,216,216};
        return is>0 ? (is-1)*sca[il-1]+off[il-1]+ip:0;
    }  
    
    public static void main(String[] args){
        ECPart part = new ECPart();  
        part.initGraphics();
        String env = System.getenv("CLAS12DIR");
 //    	part.pizeroDemo(args);
 //    	part.photonDist(args);
//     	part.photonDemo(args);
    	part.neutronDemo(args);
//        part.electronDemo(args);
//        part.scalerdemo(args);
    }
    
}
