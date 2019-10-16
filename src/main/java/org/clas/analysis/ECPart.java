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
import org.jlab.io.hipo3.Hipo3DataSource;
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
	
    List<List<DetectorResponse>>     unmatchedResponses = new ArrayList<List<DetectorResponse>>(); 
    IndexedList<List<DetectorResponse>>  singleNeutrals = new IndexedList<List<DetectorResponse>>(1);
    IndexedList<List<DetectorResponse>>      singleMIPs = new IndexedList<List<DetectorResponse>>(1);
    DetectorParticle p1 = new DetectorParticle();
    DetectorParticle p2 = new DetectorParticle();
    
    public double distance11,distance12,distance21,distance22;
    public double e1,e2,e1c,e2c,cth,cth1,cth2;
    public double X,tpi2,cpi0,refE,refP,refTH;
    public double x1,y1,x2,y2;
    public double epc,eec1,eec2;
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
    public double SF1 = 0.27;
    public double SF2 = 0.27;
    
    public H1F h5,h6,h7,h8;
  
    public int n2mc=0;

    public int[] mip = {0,0,0,0,0,0};
    
    int photonMult = 12;
    
    public ECPart() {  
        initccdb();
    	ccdb = new EBCCDBConstants(10,ebe.getConstantsManager());
    }
    
    public void initccdb() {	   
	    ebe.init();
    }
    
    public void readMC(DataEvent event) {
        int pid1=0; int pid2=0;
        double ppx1=0,ppy1=0,ppz1=0;
        double ppx2=0,ppy2=0,ppz2=0;
        double rm = 0.;
        
        Boolean isEvio = event instanceof EvioDataEvent;      
        
        if(isEvio&&event.hasBank("GenPart::true")) {
            EvioDataBank bank = (EvioDataBank) event.getBank("GenPart::true");
            ppx1 = bank.getDouble("px",0);
            ppy1 = bank.getDouble("py",0);
            ppz1 = bank.getDouble("pz",0);
            pid1 = bank.getInt("pid",0);
        }
        
        if(!isEvio&&event.hasBank("MC::Particle")) {
        	   
            DataBank bank = event.getBank("MC::Particle");
            ppx1 = bank.getFloat("px",0);
            ppy1 = bank.getFloat("py",0);
            ppz1 = bank.getFloat("pz",0);
            pid1 = bank.getInt("pid",0);  
            n2mc++;
            
            if (pid1==22&&bank.rows()==2) {  // FX two-photon runs
                ppx2 = bank.getFloat("px",1);
                ppy2 = bank.getFloat("py",1);
                ppz2 = bank.getFloat("pz",1);
                pid2 = bank.getInt("pid",1);                   
            }
        }
        
        if (pid1==11)   rm=melec;
        if (pid1==111)  rm=mpi0;
        if (pid1==2112) rm=mneut;
        
        double ppx=ppx1+ppx2; double ppy=ppy1+ppy2; double ppz=ppz1+ppz2;
        
        refP  = Math.sqrt(ppx*ppx+ppy*ppy+ppz*ppz);  
        refE  = Math.sqrt(refP*refP+rm*rm);            
        refTH = Math.acos(ppz/refP)*180/Math.PI; 
        
        h5.fill(refP);
    }
    
    public static List<DetectorResponse>  readEvent(DataEvent event, String bankName, DetectorType type){        
            List<DetectorResponse> responseList = new ArrayList<DetectorResponse>();
            if(event.hasBank(bankName)==true){
                DataBank bank = event.getBank(bankName);
                int nrows = bank.rows();
                for(int row = 0; row < nrows; row++){
                    int sector = bank.getByte("sector", row);
                    int  layer = bank.getByte("layer",  row);
                    DetectorResponse  response = new DetectorResponse(sector,layer,0);
                    response.getDescriptor().setType(type);
                    float x = bank.getFloat("x", row);
                    float y = bank.getFloat("y", row);
                    float z = bank.getFloat("z", row);
                    response.setHitIndex(row);
                    response.setPosition(x, y, z);
                    response.setEnergy(bank.getFloat("energy", row));
                    response.setTime(bank.getFloat("time", row));
                    responseList.add(response);
                }
            }
            return responseList;
    }   
    
    public List<DetectorResponse>  readEC(DataEvent event){
        eb   = new EventBuilder(ccdb);    	
        rf   = new EBRadioFrequency(ccdb);
        eb.initEvent();
        eb.getEvent().getEventHeader().setRfTime(rf.getTime(event)+ccdb.getDouble(EBCCDBEnum.RF_OFFSET));
        List<DetectorResponse> rEC = new ArrayList<DetectorResponse>();
        rEC = readEvent(event, "ECAL::clusters", DetectorType.ECAL);
        eb.addDetectorResponses(rEC); 
        eb.getPindexMap().put(0, 0); eb.getPindexMap().put(1, 0);
        return rEC;
    } 
    
    public List<DetectorParticle> getNeutralPart() {
        ebm  = new EBMatching(eb);
    	eb.processNeutralTracks();    	
    	EBAnalyzer analyzer = new EBAnalyzer(ccdb,rf);
        analyzer.processEvent(eb.getEvent());
        if(eb.getEvent().getParticles().size()>0) {
            Collections.sort(eb.getEvent().getParticles());
            eb.setParticleStatuses();  
            return eb.getEvent().getParticles();
        } 
     	return new ArrayList<DetectorParticle>() ;     	
    }
    
    public void getUnmatchedResponses(List<DetectorResponse> response) {        
        unmatchedResponses.clear();
        unmatchedResponses.add(eb.getUnmatchedResponses(response, DetectorType.ECAL,1));
        unmatchedResponses.add(eb.getUnmatchedResponses(response, DetectorType.ECAL,4));
        unmatchedResponses.add(eb.getUnmatchedResponses(response, DetectorType.ECAL,7));
    }    
        
    public void getNeutralResponses(List<DetectorResponse> response) {        
        getUnmatchedResponses(response);
       	getSingleNeutralResponses();
       	System.out.println(" ");
       	for(List<DetectorResponse> lresp : unmatchedResponses) {
       		int n=0;
       		for (DetectorResponse resp : lresp ) {
       		    System.out.println(n+" "+resp.getSector()+" "+resp.getTime()+" "+resp.getEnergy()+" "+refP);
       		    n++;
       		}
       	}
    }
    
    public void getMIPResponses(List<DetectorResponse> response) {        
        getUnmatchedResponses(response);
     	getSingleMIPResponses();
    }
    
    public void getSingleMIPResponses() {
        List<DetectorResponse> rEC = new ArrayList<DetectorResponse>();
        singleMIPs.clear();
        for (int is=1; is<7; is++) {
            rEC = DetectorResponse.getListBySector(unmatchedResponses.get(0),  DetectorType.ECAL, is);            
            if(rEC.size()==1&&mip[is-1]==1) singleMIPs.add(rEC,is);
        }     	
    }
        
    public void getSingleNeutralResponses() {
        List<DetectorResponse> rEC = new ArrayList<DetectorResponse>();
        singleNeutrals.clear();
        for (int is=1; is<7; is++) {
            rEC = DetectorResponse.getListBySector(unmatchedResponses.get(0),  DetectorType.ECAL, is);
            if(rEC.size()==1&&mip[is-1]!=1) singleNeutrals.add(rEC,is);
        } 
    }
    
    public double getTwoPhotonInvMass(int sector){
        iis[0]=0;iis[1]=-1;
        return processTwoPhotons(doHitMatching(getNeutralParticles(sector)));
    }   
    
    public double getEcalEnergy(int sector){
        iis[0]=0;epc=0;eec1=0;eec2=0;
        return processSingleMIP(doHitMatching(getMIParticles(sector)));
    }   
    
    public double processSingleMIP(List<DetectorParticle> particles) {
        if (particles.size()==0) return 0.0;
      	return particles.get(0).getEnergy(DetectorType.ECAL);
    }  
    
    public List<DetectorParticle> doHitMatching(List<DetectorParticle> particles) {
        
        for (int ii=0; ii<particles.size(); ii++) {
          	DetectorParticle      p = particles.get(ii);          
            DetectorResponse rPC = p.getDetectorResponses().get(0) ;
            epc = rPC.getEnergy();
//            System.out.println(Math.sqrt(rPC.getPosition().x()*rPC.getPosition().x()+
//                                        rPC.getPosition().y()*rPC.getPosition().y()+
//                                         rPC.getPosition().z()*rPC.getPosition().z()));
            iip[ii][0] = rPC.getHitIndex();  
            iis[ii]    = rPC.getDescriptor().getSector();
              x[ii][0] = rPC.getPosition().x();
              y[ii][0] = rPC.getPosition().y();
              t[ii][0] = rPC.getTime()-rPC.getPath()/c;;
             
            distance1[ii] = doPCECMatch(p,ii,"Inner");
            distance2[ii] = doPCECMatch(p,ii,"Outer");
            
        }
                
        return particles;
    }  
    
    public List<DetectorParticle> getMIParticles(int sector) {
    	
        List<DetectorParticle> particles = new ArrayList<DetectorParticle>();          
        List<DetectorResponse>   rPC  = new ArrayList<DetectorResponse>();        
        rPC = DetectorResponse.getListBySector(unmatchedResponses.get(0), DetectorType.ECAL, sector);
        if (rPC.size()==1) particles.add(DetectorParticle.createNeutral(rPC.get(0))); // Zero torus/solenoid fields
        return particles;
    }
     
    public List<DetectorParticle> getNeutralParticles(int sector) {
              
        List<DetectorParticle> particles = new ArrayList<DetectorParticle>();          
        List<DetectorResponse>      rEC  = new ArrayList<DetectorResponse>();        
        
        rEC = DetectorResponse.getListBySector(unmatchedResponses.get(0), DetectorType.ECAL, sector); //get PCAL responses
        
        switch (rEC.size()) {
        case 1:  List<DetectorResponse> rEC2 = findSecondPhoton(sector);
                if (rEC2.size()>0) {
                   particles.add(DetectorParticle.createNeutral(rEC.get(0))); // make neutral particle 1 from PCAL sector                   
                   particles.add(DetectorParticle.createNeutral(rEC2.get(0))); return particles; // make neutral particle 2 from other PCAL sector
                }
                break;
        case 2: particles.add(DetectorParticle.createNeutral(rEC.get(0)));                   // make neutral particle 1 from PCAL sector
                particles.add(DetectorParticle.createNeutral(rEC.get(1))); return particles; // make neutral particle 2 from PCAL sector
        }
       return particles;
    }
    
    public List<DetectorResponse> findSecondPhoton(int sector) {
        int neut=0, isave=0;
        List<DetectorResponse> rEC = new ArrayList<DetectorResponse>();        
        for (int is=sector+1; is<7; is++) {
            if(singleNeutrals.hasItem(is)) {neut++; isave=is;}
        }
        return (neut==1) ? singleNeutrals.getItem(isave):rEC;
    }
    
    public double doPCECMatch(DetectorParticle p, int ii, String io) {
        
        int index=0;
        double distance = -10;
        List<DetectorResponse> rEC = new ArrayList<DetectorResponse>();        
        
        int is = p.getDetectorResponses().get(0).getDescriptor().getSector();
       
        switch (io) {       
        case "Inner": rEC = DetectorResponse.getListBySector(unmatchedResponses.get(1), DetectorType.ECAL, is);
                      index  = p.getDetectorHit(rEC,DetectorType.ECAL,4,eb.ccdb.getDouble(EBCCDBEnum.ECIN_MATCHING));
                      if(index>=0){p.addResponse(rEC.get(index),true); rEC.get(index).setAssociation(0); 
                      iip[ii][1] = rEC.get(index).getHitIndex(); 
                      x[ii][1] = rEC.get(index).getPosition().x(); y[ii][1] = rEC.get(index).getPosition().y(); t[ii][1] = rEC.get(index).getTime()-rEC.get(index).getPath()/c;
                      distance = p.getDistance(rEC.get(index)).length();eec1=rEC.get(index).getEnergy();} break;
        case "Outer": rEC = DetectorResponse.getListBySector(unmatchedResponses.get(2), DetectorType.ECAL, is); 
                      index  = p.getDetectorHit(rEC,DetectorType.ECAL,7,eb.ccdb.getDouble(EBCCDBEnum.ECOUT_MATCHING));
                      if(index>=0){p.addResponse(rEC.get(index),true); rEC.get(index).setAssociation(0);
                      iip[ii][2] = rEC.get(index).getHitIndex();
                      x[ii][2] = rEC.get(index).getPosition().x(); y[ii][2] = rEC.get(index).getPosition().y(); t[ii][2] = rEC.get(index).getTime()-rEC.get(index).getPath()/c;
                      distance = p.getDistance(rEC.get(index)).length();eec2=rEC.get(index).getEnergy();}
        }
        
        return distance;        
    }
/*        
    public List<DetectorParticle> doHitMatching(List<DetectorParticle> particles) {
                       
        if (particles.size()==0) return particles;
        
        DetectorParticle p1 = particles.get(0);  //PCAL Photon 1
        DetectorParticle p2 = particles.get(1);  //PCAL Photon 2       
        
        CalorimeterResponse rPC1 = p1.getCalorimeterResponse().get(0) ;
        CalorimeterResponse rPC2 = p2.getCalorimeterResponse().get(0) ;
        
        ip1 = rPC1.getHitIndex();
        ip2 = rPC2.getHitIndex();
        
        is1 = rPC1.getDescriptor().getSector();
        is2 = rPC2.getDescriptor().getSector();

        x1 = rPC1.getPosition().x();
        y1 = rPC1.getPosition().y();
        x2 = rPC2.getPosition().x();
        y2 = rPC2.getPosition().y();
        
        distance11 = doHitMatch(p1,"Inner");
        distance12 = doHitMatch(p1,"Outer");
        distance21 = doHitMatch(p2,"Inner");
        distance22 = doHitMatch(p2,"Outer");
                
        return particles;
    }
    */ 

    
    public double processTwoPhotons(List<DetectorParticle> particles) {
        
        if (particles.size()==0) return 0.0;
        
        p1 = particles.get(0);  //Photon 1
        p2 = particles.get(1);  //Photon 2
       
        Vector3 n1 = p1.vector(); n1.unit();
        Vector3 n2 = p2.vector(); n2.unit();
                
        e1 = p1.getEnergy(DetectorType.ECAL);
        e2 = p2.getEnergy(DetectorType.ECAL);
        
        double SF1db = SamplingFractions.getMean(22, p1, eb.ccdb);
        e1c = e1/SF1db;
        Particle g1 = new Particle(22,n1.x()*e1c,n1.y()*e1c,n1.z()*e1c);
        
        double SF2db = SamplingFractions.getMean(22, p2, eb.ccdb);
        e2c = e2/SF2db;
        Particle g2 = new Particle(22,n2.x()*e2c,n2.y()*e2c,n2.z()*e2c);
        
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
    	
        //Require two photons in PCAL  
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
    	
        HipoDataSource  reader = new HipoDataSource();
//        Hipo3DataSource reader = new Hipo3DataSource();
        ECEngine       engine = new ECEngine();
        ECPart           part = new ECPart();  
        
        String id ;
     
        String hipoPath = "/Users/colesmith/clas12/sim/";
//        String hipoPath = "/Users/colesmith/clas12/sim/radstudy/";
        String hipo3File0 = "fc-elec-40k-s5-r2.hipo";
        String hipo3File3 = "fc-elec-40k-s5-r2.hipo";
        String hipo3File1 = "fc-ecpcsc-elec-s5-20k.hipo";
        String hipo4File1 = "fc-ecpcsc-elec-s5-20k.hipo";
        String hipo3File2 = "clasdispr-large.hipo";
        String hipo4File4 = "clas12-radstudy-5gev-22deg-pm8.hipo";
        String  hipoFile  = hipo4File1;
        		
        reader.open(hipoPath+hipoFile);
        
        engine.init();
        engine.isMC = true;
        engine.setVariation("default");
        engine.setCalRun(10);                
//.setThresholds("Electron_lo",engine);
        part.setThresholds("Pizero",engine);
        part.setGeom("2.5");
        
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
        
        part.h5 = new H1F("Throw",50,0.,10.); part.h5.setTitleX("MC Electron E (MeV)");
        part.h6 = new H1F("Recon",50,0.,10.); part.h6.setTitleX("Efficiency");
        part.h6.setTitleX("Measured Electron Energy (GeV))");
        
        int nevent = 0;
        
        while(reader.hasEvent()&&nevent<40000){
            nevent++;
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event);   
            part.readMC(event);
            part.getMIPResponses(part.readEC(event));
            double energy = part.getEcalEnergy(5);
//          Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.04 && energy>0.12;
//          Boolean trig2 = good_pcal && !good_ecal && part.epc>0.04;
//          Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.06 && energy>0.15;
//          Boolean trig2 = good_pcal && !good_ecal && part.epc>0.15;
            
//            System.out.println("energy,1,2,3 = "+energy+" "+part.epc+" "+part.eec1+" "+part.eec2);
            
            Boolean good_pcal = part.epc>0.001;               //VTP reported cluster
            Boolean good_ecal = (part.eec1+part.eec2)>0.001 ; //VTP reported cluster
            Boolean trig1 = good_pcal &&  good_ecal && part.epc>0.04 && energy>0.12;
            Boolean trig2 = good_pcal && !good_ecal && part.epc>0.12;
            		
            trig1=true; trig2=true;
            if (trig1||trig2) {
               part.h6.fill(part.refE);
         	   h1.fill(energy,energy/part.refP);
               h2.fill(part.refE,energy/part.refP);
               h3.fill(part.refP,energy);
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
        SFFunction sf = new SFFunction("esf",-11,part.eb.ccdb,0.1,2.5); 
        canvas.draw(sf,"same");
        H1F hrat1 = H1F.divide(part.h6, part.h5); hrat1.setFillColor(2);
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
        ECPart           part = new ECPart();      	
        List<DetectorParticle> np = new ArrayList<DetectorParticle>();
        
        String evioPath = "/Users/colesmith/clas12/sim/neutron/hipo/";
        String evioFile = "fc-neut-80k-s2.hipo"; int sec=2;
        
        if (args.length == 0) { 
            reader.open(evioPath+evioFile);
        } else {
            String inputFile = args[0];
            reader.open(inputFile);
        } 
        
        part.h5 = new H1F("Thrown",50,0.,3); part.h5.setTitleX("MC Neutron E (MeV)");
        H1F  h1 = new H1F("Thrown",50,0.,3);      h1.setTitleX("MC Neutron E (MeV)");
       
        engine.init();
        engine.isMC = true;
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        engine.setCalRun(2);
       
        part.setThresholds("Pizero",engine);
        part.setGeom("2.5");
        
        while(reader.hasEvent()){
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event);   
            part.readMC(event); part.readEC(event);
            np.clear();
            np= part.getNeutralPart();
       		int n=0;
//       		System.out.println(" ");
       		for (DetectorParticle neut : np) {
       			if(n==0) h1.fill(part.refP);
/*       			
       		    System.out.println(n+" "+neut.getSector(DetectorType.ECAL,1)+" "
       		                            +neut.getTime(DetectorType.ECAL)+" "
       		    		                +neut.getEnergy(DetectorType.ECAL)+" "
       		    		                +neut.getBeta(DetectorType.ECAL)+" "
       		                            +part.refP);
*/
       		    n++;
       		}	
        }
        
        JFrame frame = new JFrame("Pizero Reconstruction");
        frame.setSize(800,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        canvas.divide(2,2);
        H1F hrat1 = H1F.divide(h1,  part.h5); hrat1.setFillColor(2); hrat1.setTitleY("Neutron Eff");    hrat1.setTitleX("Neutron Momentum (GeV)");
        canvas.cd(0);  canvas.draw(part.h5);  canvas.draw(h1,"same");       
        canvas.cd(1);  canvas.draw(hrat1); 
        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);  
    }
    
    public void pizeroDemo(String[] args) {
    	
        Hipo3DataSource reader = new Hipo3DataSource();
        ECEngine        engine = new ECEngine();
        ECPart            part = new ECPart();  
        
        H2F h2a,h2b,h2c,h2d,h2e,h2f;
        H1F h5a,h5b,h5c,h6,h7a,h7b,h8;
        int n2hit=0;
        int n2rec1=0;
        int n2rec2=0;
        
        String evioPath = "/Users/colesmith/clas12/sim/pizero/hipo/";
//        String evioFile = "fc-pizero-50k-s2-newgeom-0.35-8.35.hipo"; int sec=2;
        String evioFile = "fc-pizero-50k-s2-newgeom-15-0.35-8.35.hipo"; int sec=2;
//        evioFile = "pi0_hi.hipo";
        
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
        engine.isMC = true;
        engine.setVariation("default"); // Use clas6 variation for legacy simulation 10k-s2-newgeom 
        engine.setCalRun(2);
        engine.setVeff(16); //GEMC default
        engine.setPCALTrackingPlane(9);
        
        part.setThresholds("Pizero",engine);
        part.setGeom("2.5");
        part.setGoodPhotons(12);
        
        double emax = 8.5;
        
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
        
        part.h5 = new H1F("Thrown",50,0.,emax); part.h5.setTitleX("MC Pizero E (MeV)");
        
//        h5a = new H1F("2Gamma",50,0.,20.);  h5a.setTitleX("Opening Angle (deg)");
//        h5b = new H1F("2Gamma",50,0.,20.);  h5b.setTitleX("Opening Angle (deg)");
//        h5c = new H1F("2Gamma",50,0.,20.);  h5c.setTitleX("Opening Angle (deg)");
        h6 = new H1F("2Gamma",50,0.,emax);  h6.setTitleX("MC Pizero E (MeV)"); 
        h7a = new H1F("PC/EC", 50,0.,emax); h7a.setTitleX("MC Pizero E (MeV)");
        h7b = new H1F("PC/EC", 50,0.,emax); h7b.setTitleX("MC Pizero E (MeV)");
        h8 = new H1F("Mcut",  50,0.,emax);  h8.setTitleY("ECIN 2 Photon Eff"); h8.setTitle("80<InvMass<200");
        H2F h9 = new H2F("Beta1", 50,0.75, 1.1,50,650.,750.); H2F h10 = new H2F("Beta2",50,0.75,1.1,50,650.,750.);
        
        int nimcut = 0;
        
        while(reader.hasEvent()){
            DataEvent event = reader.getNextEvent();
            engine.processDataEvent(event);   
            part.readMC(event);
            part.getNeutralResponses(part.readEC(event));
            double invmass = 1e3*Math.sqrt(part.getTwoPhotonInvMass(sec));
            
            boolean goodmass = invmass>0 && invmass<200;
            
            boolean pcec1 = goodPhotons(121,part.p1,part.p2);
            boolean pcec2 = goodPhotons(122,part.p1,part.p2);
           
            
            h9.fill(part.p1.getBeta(DetectorType.ECAL,1,0.),part.p1.getHit(DetectorType.ECAL).getPosition().z());  
            h10.fill(part.p2.getBeta(DetectorType.ECAL,1,0.),part.p2.getHit(DetectorType.ECAL).getPosition().z());

            if (goodmass) {
                                  n2hit++;   h6.fill(part.refE);  	            
                if(pcec1||pcec2) {n2rec1++; h7a.fill(part.refE);}
                if(pcec1&&pcec2) {n2rec2++; h7b.fill(part.refE);}
          
                if (pcec1&&pcec2) {
                    h2a.fill(part.refE, invmass);                                    //Two-photon invariant mass                
                    h2b.fill(part.refE, part.X);                                     //Pizero energy asymmetry
                    h2c.fill(part.refE,(Math.sqrt(part.tpi2)-part.refE));            //Pizero total energy error
                    h2e.fill(Math.acos(part.cth)*180/Math.PI,(Math.sqrt(part.tpi2)-part.refE));            //Pizero total energy error
                    h2f.fill(Math.acos(part.cth)*180/Math.PI,invmass-mpi0*1e3);            //Pizero total energy error
                    h2d.fill(part.refE,Math.acos(part.cpi0)*180/Math.PI-part.refTH); //Pizero theta angle error
                    nimcut++; h8.fill(part.refE);
                }
            }
        }
        
        H1F hrat1 = H1F.divide(h6,  part.h5); hrat1.setFillColor(2); hrat1.setTitleY("PC 2 Photon Eff");    hrat1.setTitleX("Pizero Energy (GeV)");
        H1F hrat2 = H1F.divide(h7a, part.h5); hrat2.setFillColor(2); hrat2.setTitleY("PC*EC 1 Photon Eff"); hrat2.setTitleX("Pizero Energy (GeV)");
        H1F hrat3 = H1F.divide(h7b, part.h5); hrat3.setFillColor(2); hrat3.setTitleY("PC*EC 2 Photon Eff"); hrat3.setTitleX("Pizero Energy (GeV)");
        H1F hrat4 = H1F.divide(h8, h6);  hrat4.setFillColor(2); hrat4.setTitleY("PC*EC NIMCUT Eff");        hrat4.setTitleX("Pizero Energy (GeV)");
        
        H1F h1 = h2a.projectionY();  h1.setOptStat("1100") ; h1.setFillColor(4); h1.setTitleX("Two-Photon Invariant Mass (MeV)");
        H1F h2 = h2b.projectionY();  h2.setOptStat("1100") ; h2.setFillColor(4); h2.setTitleX("X:(E1-E2)/(E1+E2)");
        H1F h3 = h2c.projectionY();  h3.setOptStat("1100") ; h3.setFillColor(4); h3.setTitleX("Pizero Energy Error");
        H1F h4 = h2d.projectionY();  h4.setOptStat("1100") ; h4.setFillColor(4); h4.setTitleX("Pizero Theta Error (deg)");
       
        System.out.println("THROWN TWOPHOTONS PCEC MATCH1 PCEC MATCH2 INVMASS CUT");
        System.out.println(part.n2mc+"     "+n2hit+"         "+n2rec1+"     "+n2rec2+"   "+nimcut);
        System.out.println("Eff1 = "+(float)n2hit/(float)part.n2mc+
                          " Eff2 = "+(float)n2rec1/(float)part.n2mc+
                          " Eff3 = "+(float)n2rec2/(float)part.n2mc+
                          " Eff4 = "+(float)nimcut/(float)n2hit);
        
        JFrame frame = new JFrame("Pizero Reconstruction");
        frame.setSize(800,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        canvas.divide(4,4);
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
     	part.neutronDemo(args);
//        part.electronDemo(args);
    }
    
}
