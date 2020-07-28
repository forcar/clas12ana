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
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.Func1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.evio.EvioDataEvent;
import org.jlab.io.hipo.HipoDataSource;
//import org.jlab.io.hipo3.Hipo3DataSource;
import org.jlab.service.eb.EBAnalyzer;
//import org.jlab.service.eb.EBConstants;
import org.jlab.service.eb.EBEngine;
import org.jlab.service.eb.EBMatching;
import org.jlab.service.ec.ECEngine;
import org.jlab.service.eb.EventBuilder;
import org.jlab.utils.groups.IndexedList;
import org.jlab.rec.eb.EBCCDBConstants;
import org.jlab.rec.eb.EBCCDBEnum;
import org.jlab.rec.eb.EBRadioFrequency;
import org.jlab.rec.eb.EBUtil;
import org.jlab.rec.eb.SamplingFractions;

public class ECPart  {
	
	public EventBuilder eb = null;
	EBEngine         ebe = new EBEngine("ECPART");
	EBCCDBConstants  ccdb;
	EBMatching       ebm;
	EBRadioFrequency rf;
	
    public List<List<DetectorResponse>>     unmatchedResponses = new ArrayList<>(); 
    
    IndexedList<List<DetectorResponse>>  singleNeutrals = new IndexedList<>(1);
    IndexedList<List<DetectorResponse>>      singleMIPs = new IndexedList<>(1);
    DetectorParticle p1 = new DetectorParticle();
    DetectorParticle p2 = new DetectorParticle();
    
    public double distance11,distance12,distance21,distance22;
    public double e1,e2,e1c,e2c,cth,cth1,cth2;
    public double X,tpi2,cpi0,refE,refP,refTH;
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
    
    public ECPart() {  
    	
    }
    
    public void getCCDB(int runno) {
    	System.out.println("ECpart.setCCDB()");
    	ebe.init();
    	this.ccdb = new EBCCDBConstants(runno,ebe.getConstantsManager());    	
    }
    
    public boolean readMC(DataEvent event) {
        int pid1=0, pid2=0, off=0;
        double ppx1=0,ppy1=0,ppz1=0;
        double ppx2=0,ppy2=0,ppz2=0;
        double rm = 0.;
        
        Boolean isEvio = event instanceof EvioDataEvent;      
        
        if (!isEvio && !event.hasBank("MC::Particle")) return false;
        
        if(isEvio&&event.hasBank("GenPart::true")) {
            EvioDataBank bank = (EvioDataBank) event.getBank("GenPart::true");
            ppx1 = bank.getDouble("px",0);
            ppy1 = bank.getDouble("py",0);
            ppz1 = bank.getDouble("pz",0);
            pid1 = bank.getInt("pid",0);
        }
        
        if(!isEvio&&event.hasBank("MC::Particle")) {
            DataBank bank = event.getBank("MC::Particle");
            
           if(bank.rows()==3) off=1;

            ppx1 = bank.getFloat("px",0+off);
            ppy1 = bank.getFloat("py",0+off);
            ppz1 = bank.getFloat("pz",0+off);
            pid1 = bank.getInt("pid",0+off);  
            n2mc++;
            
            if (pid1==22&&bank.rows()>1) {  // Two-photon runs
                ppx2 = bank.getFloat("px",1+off);
                ppy2 = bank.getFloat("py",1+off);
                ppz2 = bank.getFloat("pz",1+off);
                pid2 = bank.getInt("pid",1+off);                   
            }
        }
        
        if (pid1==11)           rm=melec;
        if (pid1==111||off==1)  rm=mpi0;
        if (pid1==2112)         rm=mneut;
        
        double ppx=ppx1+ppx2; double ppy=ppy1+ppy2; double ppz=ppz1+ppz2;
        
        refP  = Math.sqrt(ppx*ppx+ppy*ppy+ppz*ppz);  
        refE  = Math.sqrt(refP*refP+rm*rm);            
        refTH = Math.acos(ppz/refP)*180/Math.PI; 
        
        return true;
    }
    
    public static List<DetectorResponse>  getResponses(DataEvent event, String bankName, DetectorType type){        
            List<DetectorResponse> responseList = new ArrayList<>();
            if(event.hasBank(bankName)==true){
                DataBank bank = event.getBank(bankName);
                for(int row = 0; row < bank.rows(); row++){
                    int sector = bank.getByte("sector", row);
                    int  layer = bank.getByte("layer",  row);
                    DetectorResponse  response = new DetectorResponse(sector,layer,0);
                    response.getDescriptor().setType(type);
                    float x = bank.getFloat("x", row);
                    float y = bank.getFloat("y", row);
                    float z = bank.getFloat("z", row);
                    response.setHitIndex(bankName.equals("ECAL::clusters")?row:bank.getInt("index",row));
                    response.setPosition(x, y, z);
                    response.setEnergy(bank.getFloat("energy", row));
                    response.setTime(bank.getFloat("time", row));
                    responseList.add(response);
                }
            }
            return responseList;
    } 
    
	// readEC: Copies relevant parts of EBEngine.processDataEvent    
    public void  readEC(DataEvent event, String bank){    	
        rf = new EBRadioFrequency(ccdb);    	
        eb = new EventBuilder(ccdb);    	
        eb.initEvent(); //don't bother with event header
        eb.getEvent().getEventHeader().setRfTime(rf.getTime(event)+ccdb.getDouble(EBCCDBEnum.RF_OFFSET));        
        eb.addDetectorResponses(getResponses(event, bank, DetectorType.ECAL)); 
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
    
    public void getNeutralResponses() {        
        getUnmatchedResponses();
       	getSingleNeutralResponses(); // For two-photon decays in different sectors     	
    }
    
    public void getSingleMIPResponses() {
        List<DetectorResponse> rEC = new ArrayList<>();
        singleMIPs.clear();
        for (int is=1; is<7; is++) {
            rEC = DetectorResponse.getListBySector(unmatchedResponses.get(0),  DetectorType.ECAL, is);  //look in PCAL only          
            if(rEC.size()==1&&mip[is-1]==1) singleMIPs.add(rEC,is);
        }     	
    }
        
    public void getSingleNeutralResponses() {
        List<DetectorResponse> rEC = new ArrayList<>();
        singleNeutrals.clear();
        for (int is=1; is<7; is++) {
            rEC = DetectorResponse.getListBySector(unmatchedResponses.get(0),  DetectorType.ECAL, is); //look in PCAL only
            if(rEC.size()==1&&mip[is-1]!=1) singleNeutrals.add(rEC,is);
        } 
    }
    
    public List<DetectorResponse> findSecondPhoton(int sector) {
        int neut=0, isave=0;
        for (int is=sector+1; is<7; is++) {
            if(singleNeutrals.hasItem(is)) {neut++; isave=is;}
        }
        return (neut==1) ? singleNeutrals.getItem(isave):new ArrayList<DetectorResponse>();
    }
    
    public double getEcalEnergy(int sector){
        iis[0]=0;epc=0;eec1=0;eec2=0;
        return processSingleMIP(doHitMatching(getMIParticles(sector)));
    }   
    
    public double processSingleMIP(List<DetectorParticle> particles) {
        if (particles.size()==0) return 0.0;
      	return particles.get(0).getEnergy(DetectorType.ECAL);
    }  
    
    public List<DetectorParticle> getMIParticles(int sector) {
    	
        List<DetectorParticle> particles = new ArrayList<>();          
        List<DetectorResponse>      rPC  = new ArrayList<>();        
        rPC = DetectorResponse.getListBySector(unmatchedResponses.get(0), DetectorType.ECAL, sector);
        if (rPC.size()==1) particles.add(DetectorParticle.createNeutral(rPC.get(0))); // Zero torus/solenoid fields
        return particles;
    }
     
    public double getTwoPhotonInvMass(int sector){
        iis[0]=0;iis[1]=-1;
        return processTwoPhotons(doHitMatching(getNeutralParticles(sector)));
//        return processTwoPhotons(getNeutralPart());
    }   
    
    // getNeutralParticles: Similar to EBMatching.findNeutrals(1)
    public List<DetectorParticle> getNeutralParticles(int sector) {
              
        List<DetectorResponse>      rEC  = new ArrayList<>();        
        List<DetectorParticle> particles = new ArrayList<>();          
        
        rEC = DetectorResponse.getListBySector(unmatchedResponses.get(0), DetectorType.ECAL, sector); //get PCAL responses

        switch (rEC.size()) {
        case 1:  List<DetectorResponse> rEC2 = findSecondPhoton(sector);
                if (rEC2.size()>0) {
                   particles.add(DetectorParticle.createNeutral(rEC.get(0)));                    // make neutral particle 1 from PCAL sector                   
                   particles.add(DetectorParticle.createNeutral(rEC2.get(0))); return particles; // make neutral particle 2 from other PCAL sector
                }
                break;
        case 2: particles.add(DetectorParticle.createNeutral(rEC.get(0)));                   // make neutral particle 1 from PCAL sector
                particles.add(DetectorParticle.createNeutral(rEC.get(1))); return particles; // make neutral particle 2 from PCAL sector
//        case 3: particles.add(DetectorParticle.createNeutral(rEC.get(0)));                   // make neutral particle 1 from PCAL sector
//        		particles.add(DetectorParticle.createNeutral(rEC.get(1))); return particles; // make neutral particle 2 from PCAL sector
        }
        
       return particles;
    }
    
    public List<DetectorParticle> doHitMatching(List<DetectorParticle> particles) {
    	
        int ii=0;
        for (DetectorParticle p: particles) {
            DetectorResponse rPC = p.getDetectorResponses().get(0); //get 1st PCAL responses
//                   epc = rPC.getEnergy();
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
        
        SF1db = SamplingFractions.getMean(22, p1, eb.ccdb);
        e1c = e1/SF1db;
        Particle g1 = new Particle(22,n1.x()*e1c,n1.y()*e1c,n1.z()*e1c);        
        VG1 = new LorentzVector(n1.x()*e1c,n1.y()*e1c,n1.z()*e1c,e1c);
        
        SF2db = SamplingFractions.getMean(22, p2, eb.ccdb);
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
       
       boolean goodPhotons = pc12;
       
       switch (test) {
           case 121:  goodPhotons = pc12 && pcec1;          break;
           case 122:  goodPhotons = pc12 && pcec2;          break;
           case 120:  goodPhotons = pc12 && (pcec1||pcec2); break;
           case 1212: goodPhotons = pc12 && (pcec1&&pcec2);
       }
       
       return goodPhotons;
       
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
//    					   engine.setStripThresholds(20,19,18);
                           engine.setPeakThresholds(18,20,15); 
//                           engine.setPeakThresholds(28,30,25);
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
    
    public void electronDemo(String[] args) {
    	
        HipoDataSource reader = new HipoDataSource();
        ECEngine       engine = new ECEngine();
        
        String id ;
     
        String hipoPath = "/Users/colesmith/clas12/sim/";
//        String hipoPath = "/Users/colesmith/clas12/sim/radstudy/";
        String hipo3File0 = "fc-elec-40k-s5-r2.hipo";
        String hipo3File3 = "fc-elec-40k-s5-r2.hipo";
        String hipo3File1 = "fc-ecpcsc-elec-s5-20k.hipo";
        String hipo4File1 = "fc-ecpcsc-elec-s5-20k.hipo";
        String hipo3File2 = "clasdispr-large.hipo";
        String hipo4File4 = "clas12-radstudy-5gev-22deg-pm8.hipo";
        String hipo4File3 = "rga-fall2018-1.hipo";
        String  hipoFile  = "radstudy/"+hipo4File4;
        		
        reader.open(hipoPath+hipoFile);
        
        engine.init();
        engine.isMC = true;
        engine.setVariation("default");
        engine.setCalRun(10);                
//.setThresholds("Electron_lo",engine);
        getCCDB(10);
        setThresholds("Pizero",engine);
        setGeom("2.5");
        
        id="_s"+Integer.toString(5)+"_l"+Integer.toString(0)+"_c";
        H2F h1 = new H2F("E over P"+id+0,50,0.0,2.5,50,0.15,0.31);      
        h1.setTitleX("Measured Electron Energy (GeV))");
        h1.setTitleY("Sampling Fraction (E/P)");
        
        id="_s"+Integer.toString(5)+"_l"+Integer.toString(0)+"_c";
        H2F h2 = new H2F("E over P"+id+1,50,0.0,10.,50,0.15,0.31);      
        h2.setTitleX("True Electron Energy (GeV))");
        h2.setTitleY("Sampling Fraction (E/P)");
        
        id="_s"+Integer.toString(5)+"_l"+Integer.toString(0)+"_c";
        H2F h3 = new H2F("E vs. P"+id+2,100,0.0,10.,100,0.,2.7);      
        h3.setTitleX("True Electron Energy (GeV))");
        h3.setTitleY("Measured Electron Energy (GeV)");
        
        h5 = new H1F("Throw",50,0.,10.); h5.setTitleX("MC Electron E (MeV)");
        h6 = new H1F("Recon",50,0.,10.); h6.setTitleX("Efficiency");
        h6.setTitleX("Measured Electron Energy (GeV))");
        
        int nevent = 0;
        
        while(reader.hasEvent()&&nevent<40000){
            nevent++;
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event);   
            readMC(event);
            readEC(event,"ECAL::clusters") ;           
            getMIPResponses();
            double energy = getEcalEnergy(5);
//          Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.04 && energy>0.12;
//          Boolean trig2 = good_pcal && !good_ecal && part.epc>0.04;
//          Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.06 && energy>0.15;
//          Boolean trig2 = good_pcal && !good_ecal && part.epc>0.15;
            
//            System.out.println("energy,1,2,3 = "+energy+" "+part.epc+" "+part.eec1+" "+part.eec2);
            
            Boolean good_pcal = epc>0.001;               //VTP reported cluster
            Boolean good_ecal = (eec1+eec2)>0.001 ; //VTP reported cluster
            Boolean trig1 = good_pcal &&  good_ecal && epc>0.04 && energy>0.12;
            Boolean trig2 = good_pcal && !good_ecal && epc>0.12;
            		
            trig1=true; trig2=true;
            if (trig1||trig2) {
               h6.fill(refE);
         	   h1.fill(energy,energy/refP);
               h2.fill(refE,energy/refP);
               h3.fill(refP,energy);
            }
        }
        
        JFrame frame = new JFrame("Electron Reconstruction");
        frame.setSize(800,800);
        GStyle.getH1FAttributes().setMarkerSize(4);
        GStyle.getH1FAttributes().setMarkerStyle(1);
        GStyle.getH1FAttributes().setLineWidth(1);
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        
	    ParallelSliceFitter fitter = new ParallelSliceFitter(h1);
        fitter.fitSlicesX();
        
        GraphErrors sigGraph = fitter.getSigmaSlices();   
        sigGraph.setTitleX("Measured Electron Energy (GeV)");  sigGraph.setTitleY("SIGMA(E/P)"); 	
        sigGraph.setMarkerSize(4); sigGraph.setMarkerStyle(1);
        
        GraphErrors meanGraph = fitter.getMeanSlices();  
        meanGraph.setTitleX("Measured Electron Energy (GeV)"); meanGraph.setTitleY("Sampling Fraction (E/P)"); 	
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
               
        canvas.divide(3,3);
		canvas.setAxisTitleSize(18);
		canvas.setAxisLabelSize(18);
        canvas.cd(0); canvas.draw(h1);        
        canvas.cd(1); canvas.draw(h2);
        canvas.cd(2); canvas.draw(h3);
        canvas.cd(3); canvas.getPad(3).getAxisY().setRange(0.0,0.028) ; canvas.draw(sigGraph);
        canvas.cd(4); canvas.getPad(4).getAxisX().setRange(0.1,0.51) ; 
                      canvas.getPad(4).getAxisY().setRange(0.0030,0.0075) ; canvas.draw(resGraph);
//        canvas.cd(4); canvas.getPad(4).getAxisX().setRange(0.00,1.5) ; 
//                      canvas.getPad(4).getAxisY().setRange(0.00,0.18) ; canvas.draw(resGraph);
        canvas.cd(5); canvas.getPad(5).getAxisY().setRange(0.18,0.26) ; canvas.draw(meanGraph);
        SFFunction sf = new SFFunction("esf",-11,eb.ccdb,0.1,2.5); 
        canvas.draw(sf,"same");
        H1F hrat1 = H1F.divide(h6,h5); hrat1.setFillColor(2);
        hrat1.setTitleX("True Electron Energy (GeV))"); hrat1.setTitleY("Efficiency");
        canvas.cd(6); canvas.getPad(6).getAxisY().setRange(0.85,1.02);canvas.draw(hrat1);
        
        frame.add(canvas);
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
        List<DetectorParticle> np = new ArrayList<DetectorParticle>();
        int run = this.runNumber;
        
        String evioPath = "/Users/colesmith/clas12/sim/neutron/hipo/";
        String evioFile = "fc-neut-80k-s2-r5424.hipo"; int sec=2;
        
        if (args.length == 0) { 
            reader.open(evioPath+evioFile);
        } else {
            String inputFile = args[0];
            reader.open(inputFile);
        } 
        
        h5 = new H1F("Thrown",50,0.,3); h5.setTitleX("MC Neutron E (MeV)");
        H1F  h1 = new H1F("Thrown",50,0.,3);      h1.setTitleX("MC Neutron E (MeV)");
       
        engine.init();
        engine.isMC = true;
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        engine.setCalRun(2);
        getCCDB(2);
        
        setThresholds("Pizero",engine);
        setGeom("2.5");
        
        while(reader.hasEvent()){
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event);
            run = getRunNumber(event);
            readMC(event); readEC(event,"ECAL::clusters");
            np.clear();
            np= getNeutralPart();
       		int n=0;
//       		System.out.println(" ");
       		h5.fill(refP);
       		for (DetectorParticle neut : np) {
       			if(n==0) h1.fill(refP);
       			
/*       		    System.out.println(n+" "+neut.getSector(DetectorType.ECAL,1)+" "
       		                            +neut.getTime(DetectorType.ECAL)+" "
       		    		                +neut.getEnergy(DetectorType.ECAL)+" "
       		    		                +neut.getBeta(DetectorType.ECAL)+" "
       		                            +refP);
*/
       		    n++;
       		}	
        }
        
        JFrame frame = new JFrame("Neutron Reconstruction");
        frame.setSize(800,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        canvas.divide(2,2);
        H1F hrat1 = H1F.divide(h1,  h5); hrat1.setFillColor(2); hrat1.setTitleY("Neutron Eff");    hrat1.setTitleX("Neutron Momentum (GeV)");
        canvas.cd(0);  canvas.draw(h5);  canvas.draw(h1,"same");       
        canvas.cd(1);  canvas.draw(hrat1); 
        dumpGraph("/Users/colesmith/CLAS12ANA/ECpart/files/neuteff_r"+run+".vec",hrat1.getGraph());

        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);  
    }
    
    public void photonDemo(String[] args) {
    	
        HipoDataSource reader = new HipoDataSource();
        ECEngine       engine = new ECEngine();
        List<DetectorParticle> np = new ArrayList<DetectorParticle>();
        int run = this.runNumber;
        
        String evioPath = "/Users/colesmith/clas12/sim/photon/hipo/";
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
        engine.isMC = true;
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        engine.setPCALTrackingPlane(9);
        engine.setCalRun(2);
        
        getCCDB(2);
        setThresholds("Pizero",engine);
        setGeom("2.5");
        setGoodPhotons(12);
      
        while(reader.hasEvent()){
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event);
            run = getRunNumber(event);
            if (readMC(event)) {
            	readEC(event,"ECAL::clusters");
            	np.clear(); np = getNeutralPart();
            	h5.fill(refP);
            	
            	if(np.size()>0) {
            		h11.fill(refP);
        	    	int n1=0,n4=0,n7=0,sum=0;
        	    	double e1=0, e4=0, e7=0, e1p=0, e4p=0, e7p=0;
            	    for (DetectorParticle phot : np) {
            	    	for (DetectorResponse dr : phot.getDetectorResponses()) {
            	    		int lay = dr.getDescriptor().getLayer();
            	    		double e = dr.getEnergy();
            	    		if(lay==1) {n1++; e1p+=e; h211.fill(n1); if(n1==1) {e1=e; h210.fill(refP);}}
            	    		if(lay==4) {n4++; e4p+=e; h221.fill(n4); if(n4==1) {e4=e; h220.fill(refP);}}
            	    		if(lay==7) {n7++; e7p+=e; h231.fill(n7); if(n7==1) {e7=e; h230.fill(refP);}}
            	    	}
            	    }
            	    h21.fill(refP,n1); h22.fill(refP,n4); h23.fill(refP,n7);
                    if(n1>1) h212.fill(e1/e1p);
                    if(n4>1) h222.fill(e4/e4p);
                    if(n7>1) h232.fill(e7/e7p);
            	}
            	if(np.size()==1) {
            		h12.fill(refP);
            	    for (DetectorParticle phot : np) {
            	    	int sum=0;
            	    	for (DetectorResponse dr : phot.getDetectorResponses()) {
            	    		int lay = dr.getDescriptor().getLayer();
            	    		if(lay==1) sum+= 100;
            	    		if(lay==4) sum+=  40;
            	    		if(lay==7) sum+=   7;
            	    	}
            	    	if(sum>=100) h13.fill(refP);
            	    	if(sum>=140) h14.fill(refP);
            	    	if(sum==147) h15.fill(refP);
            		}
            	}
            }
        }
        
        JFrame frame = new JFrame("Photon Reconstruction");
        frame.setSize(800,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas(); canvas.divide(4,3);
        
        H1F hr11 = H1F.divide(h11,  h5); hr11.setFillColor(3); hr11.setTitleY("Photon Efficiency n>0"); hr11.setTitleX("Photon Momentum (GeV)");
        H1F hr12 = H1F.divide(h12,  h5); hr12.setFillColor(1); hr12.setTitleY("Photon Efficiency n=1"); hr12.setTitleX("Photon Momentum (GeV)");
        H1F hr13 = H1F.divide(h13,  h5); hr13.setFillColor(2); hr13.setTitleY("Photon Eff 1");   hr13.setTitleX("Photon Momentum (GeV)");
        H1F hr14 = H1F.divide(h14,  h5); hr14.setFillColor(5); hr14.setTitleY("Photon Eff 14");  hr14.setTitleX("Photon Momentum (GeV)");
        H1F hr15 = H1F.divide(h15,  h5); hr15.setFillColor(4); hr15.setTitleY("Photon Eff 147"); hr15.setTitleX("Photon Momentum (GeV)");
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
        ECEngine       engine = new ECEngine();
        
        H2F h2a,h2b,h2c,h2d,h2e,h2f;
        H1F h5a,h5b,h5c,h6,h7a,h7b,h8;
        
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
                
        String evioPath = "/Users/colesmith/clas12/sim/pizero/hipo/";
//        String evioFile = "fc-pizero-50k-s2-newgeom-0.35-8.35.hipo"; int sec=2;
//      String evioFile = "fc-pizero-50k-s2-newgeom-15-0.35-8.35.hipo4"; int sec=2;
//    String evioFile = "fc-pizero-50k-s2-newgeom-15-0.35-8.35-r5716.hipo"; int sec=2;
//        String evioFile = "fc-pizero-50k-s2-newgeom-0.35-8.35-r5716.hipo"; int sec=2;
        String evioFile = "fc-pizero-100k-s2-newgeom-15-1.0-12.0.hipo"; int sec=2;
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
              
//        engine.setGeomVariation("rga_fall2018");
        engine.init();
        engine.isMC = isMC;
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        engine.setVeff(16); //GEMC default
        engine.setPCALTrackingPlane(9);
        
        getCCDB(10);
        setThresholds("Pizero",engine);
        setGeom("2.5");
        setGoodPhotons(12);
                
        while(reader.hasEvent()){
            DataEvent event = reader.getNextEvent();
            dropBanks(event);
            if (engine.processDataEvent(event)) {   
//            int iview = engine.getSharedView();
            
            if(readMC(event)) {  
            	readEC(event,"ECAL::clusters");
            	getNeutralResponses();
            	double   invmass = 1e3*Math.sqrt(getTwoPhotonInvMass(sec));            
            	boolean goodmass = invmass>0 && invmass<200;            
            	boolean    pcec1 = goodPhotons(121,p1,p2);
            	boolean    pcec2 = goodPhotons(122,p1,p2);
           
//            h9.fill(part.p1.getBeta(DetectorType.ECAL,1,0.),part.p1.getHit(DetectorType.ECAL).getPosition().z());  
//            h10.fill(part.p2.getBeta(DetectorType.ECAL,1,0.),part.p2.getHit(DetectorType.ECAL).getPosition().z());

            	h5.fill(refE);
            	if (goodmass) {
                                  	n2hit++;   h6.fill(refE);  	            
                  if(pcec1||pcec2) {n2rec1++; h7a.fill(refE);}
                  if(pcec1&&pcec2) {n2rec2++; h7b.fill(refE);}
          
                  if (pcec1&&pcec2) {
//                	  if( engine.hasSharedView()) {hview[iview].fill(invmass);hview[3].fill(invmass);}
                	  if(true) {
                	  h2a.fill(refE, invmass);                                    			  //Two-photon invariant mass                
                	  h2b.fill(refE, X);                                     			  //Pizero energy asymmetry
                	  h2c.fill(refE,(Math.sqrt(tpi2)-refE));            			  //Pizero total energy error
                	  h2e.fill(Math.acos(cth)*180/Math.PI,(Math.sqrt(tpi2)-refE)); //Pizero total energy error
                	  h2f.fill(Math.acos(cth)*180/Math.PI,invmass-mpi0*1e3);            	  //Pizero total energy error
                	  h2d.fill(refE,Math.acos(cpi0)*180/Math.PI-refTH); 			  //Pizero theta angle error
                	  nimcut++; h8.fill(refE);
                	  }
                  }
            	}
            }
            }
        }
        
        H1F hrat1 = H1F.divide(h6,  h5); hrat1.setFillColor(2); hrat1.setTitleY("PC 2 Photon Eff");    hrat1.setTitleX("Pizero Energy (GeV)");
        H1F hrat2 = H1F.divide(h7a, h5); hrat2.setFillColor(2); hrat2.setTitleY("PC*EC 1 Photon Eff"); hrat2.setTitleX("Pizero Energy (GeV)");
        H1F hrat3 = H1F.divide(h7b, h5); hrat3.setFillColor(2); hrat3.setTitleY("PC*EC 2 Photon Eff"); hrat3.setTitleX("Pizero Energy (GeV)");
        H1F hrat4 = H1F.divide(h8, h6);  hrat4.setFillColor(2); hrat4.setTitleY("PC*EC NIMCUT Eff");        hrat4.setTitleX("Pizero Energy (GeV)");
        
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
        canvas.divide(4,5);
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
        
	    ParallelSliceFitter fitter = new ParallelSliceFitter(h2a);
        fitter.fitSlicesX();  
        dumpGraph("/Users/colesmith/pi0fit.vec",fitter.getSigmaSlices());
        
        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);     
    }
  	 
    public static void main(String[] args){
        ECPart part = new ECPart();  
        part.initGraphics();
//     	part.pizeroDemo(args);
     	part.photonDemo(args);
//     	part.neutronDemo(args);
//        part.electronDemo(args);
    }
    
}
