package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.service.ec.ECCommon;
import org.clas.tools.EBMCEngine;
import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.H1F;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;

public class ECmcn extends DetectorMonitor {
	
	Event          ev = new Event();
	EBMCEngine  ebmce = new EBMCEngine();	
	
	List<Float>           GEN = new ArrayList<Float>();
	List<Float>           REC = new ArrayList<Float>(); 
	
	List<Particle>       phot = new ArrayList<Particle>();
	List<Particle>       neut = new ArrayList<Particle>();
	
	HipoDataSync       writer = new HipoDataSync();	
	
	String tit = null;
	double ethresh = 0.3;	
	
    public ECmcn(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("DEMO","NEUTRONS","CLUSTERS","RECGEN","EFFICIENCY");

        this.use123Buttons(true);
        this.useZSliderPane(true);

        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println(getDetectorName()+".localinit()");
        
        engine.setIsMC(true);
        engine.setGeomVariation("rga_fall2018"); 
        engine.setDebug(false);
        setEngineConfig("phot");
        
        ebmce.getCCDB(10);
        ebmce.setGeom("2.5");
        ebmce.setMCpid(2112);
        ebmce.setMCsec(2);
        ebmce.isMC = true;
        
        tl.setFitData(Fits);
    }
    
    public void localclear() {
    	System.out.println(getDetectorName()+".localclear()");
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
	    System.out.println(getDetectorName()+".createHistos("+run+")");
        if(dumpFiles) {
        	System.out.println(getDetectorName()+" Writing REC banks");
        	writer.open("/Users/colesmith/CLAS12ANA/ECmcn/neutron_mc.hipo");
        }
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;    	
    	createDEMO(0);
    	createNEUTRONS(0);
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);
    	plotAnalysis(run);
    }
    
    public void plotSummary(int run) {
    	setRunNumber(run);
    	plotMCHistos();
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
   	    if(!isAnalyzeDone) return;
        plot("DEMO");  
    }
       
    public void plotMCHistos() {  
        plot("NEUTRONS");
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
    	geteff();
    	if(dumpFiles) writer.close();
    	isAnalyzeDone = true;
    }
    
	public void plot(String tabname) {    
    	drawGroup(tabname,0,getActive123(), getRunNumber());
    }    
    
    public void createDEMO(int st) {
    	switch (st) {
    	case 0:   
    		dgm.add("DEMO", 2,4,0,st,getRunNumber());
    		dgm.makeH1("h1",50,0,3,-1,"GRN:GEN    BLK:REC    RED:ECAL ONLY","GEN P (GeV)","COUNTS",3,3);
    		dgm.makeH1("h2",50,0,3,-2,"","","",1,1);
    		dgm.makeH1("h3",50,0,3,-2,"","","",2,2);
    		dgm.makeGraph("g21",-1,"BLK:REC  RED:ECAL ONLY","GEN P (GeV)","EFFICIENCY",1); dgm.cc("g21",false,false,0,0.85f,0,0); 
    		dgm.makeGraph("g31",-2," "," "," ",2);                                         dgm.cc("g31",false,false,0,0.85f,0,0); 
    		dgm.makeH2("h6", 50, 0, 3, 100, 0, 1000, -1, "ANY LAYER", "GEN P (GeV)", "REC E (MeV)");
    		dgm.cc("h6", false, true, 0, 0, 0, 0);
    		dgm.makeH2("h7", 50, 0, 3, 100, 0, 1000, -1, "PCAL+ECIN+ECOU", "GEN P (GeV)", "REC E (MeV)");
    		dgm.cc("h7", false, true, 0, 0, 0, 0);
    		dgm.makeH2("h8", 50, 0, 3, 10, 1, 11, -1, "", "GEN P (GeV)", "REC NEUTRONS");
    		dgm.cc("h8", false, true, 0, 0, 1, 1000);
    		dgm.makeH2("h9", 50, 0, 3, 10, 1, 11, -1, "", "GEN P (GeV)", "REC NEUTRONS");
    		dgm.cc("h9", false, true, 0, 0, 1, 1000);
    		dgm.makeH1("h10",50,0,2200,-1,"","GEN K.E. (MeV)","",1);
    		dgm.makeH1("h11",50,0,2200,-1,"","GEN K.E. (MeV)","",1);
    		break;  		
    	}
    }
    
    public void createNEUTRONS(int st) {
    	switch (st) {
    	case 0:
    		dgm.add("NEUTRONS", 3,2,0,st,getRunNumber());
//    		dgm.makeH2("n00",100,0.2,1.3,9,0.5,9.5,-1,"","Neutron #beta","Neutrons");
//    		dgm.makeH2("n01",100,0,2,    9,0.5,9.5,-1,"","Neutron REC Energy (GeV)","Neutrons");
    		dgm.makeH2("n100", 50,0,6,   50,0.2,1.3,-1,"","Neutron GEN Momentum (GeV)","PCAL #beta");
    		dgm.makeH2("n010", 50,0,6,   50,0.2,1.3,-1,"","Neutron GEN Momentum (GeV)","ECIN #beta");
    		dgm.makeH2("n001", 50,0,6,   50,0.2,1.3,-1,"","Neutron GEN Momentum (GeV)","ECOU #beta");
    	}
    }
    
    @Override
    public void processEvent(DataEvent de) {
		
    	GEN.clear();  
    	
		boolean goodev = ebmce.readMC(de) && ebmce.pmc.size()==1;
       
        if (goodev) {         	
        	GEN   = ebmce.getkin(ebmce.pmc); dgm.fill("h1",GEN.get(0));        	
	    	if (dropBanks) dropBanks(de);  //drop ECAL banks and re-run ECEngine            
        	if(!ebmce.processDataEvent(de)) return;        	
        	processDEMO(); 
        	processNEUTRONS();
            if(dumpFiles) {ebmce.getRECBanks(de,ebmce.eb); writer.writeEvent(de);}     
        }
    }
    
    public void processDEMO() {
    	
    	float refP = GEN.get(0), refE = (float)ebmce.pmc.get(0).e();
    	
	    boolean n1=true, n2=true; float en=0; int np=0;
	    boolean goodpc=false, goodeci=false, goodeco=false;
		
	    for (DetectorParticle neut : ebmce.eb.getEvent().getParticles()) { np++;          			
			    for (DetectorResponse ds : neut.getDetectorResponses()) { 
				    if (ds.getDescriptor().getType()==DetectorType.ECAL) {      					
					    if(n1 && ds.getDescriptor().getLayer()>0) {dgm.fill("h2",refP); n1=false;} //PCAL+ECAL
					    if(n2 && ds.getDescriptor().getLayer()>1) {dgm.fill("h3",refP); n2=false;} //ECAL ONLY
					    en+=ds.getEnergy()*1e3;
					    if(ds.getDescriptor().getLayer()==1) goodpc=true;
					    if(ds.getDescriptor().getLayer()==4) goodeci=true;
					    if(ds.getDescriptor().getLayer()==7) goodeco=true;       					
			    }          		
			    }
		    }       			
        if(en>0)                     {dgm.fill("h6", refP,en); dgm.fill("h8",refP, np); dgm.fill("h10",1e3*(refE-0.93957));}
        if(goodpc&&goodeci&&goodeco) {dgm.fill("h7", refP,en); dgm.fill("h9",refP, np); dgm.fill("h11",1e3*(refE-0.93957));}    	
    }
    
    public void processNEUTRONS() {
    	
    	int npart=0;	
        float stt = ebmce.starttime;
        
    	DetectorParticle p1 = new DetectorParticle();
    	DetectorParticle p2 = new DetectorParticle();
        
        List<DetectorParticle> par = ebmce.eb.getEvent().getParticles();  
        List<DetectorResponse> cal = ebmce.eb.getEvent().getCalorimeterResponseList();         
        
        int trsec = -1; int trpid = -211;
        for (DetectorParticle dp: par) { //find sector of trpid
        	int pid = dp.getPid(), sec = dp.getSector(DetectorType.ECAL);
        	if(trsec==-1 && sec>0 && pid==trpid) trsec=sec;
        }
        	
        if(trsec==-1) return;
        
        for (DetectorParticle dp : par) if(dp.getSector(DetectorType.ECAL)==2 && dp.getPid()==2112) npart++;

 		double dist=0, du=0;
 		
		int[]     npc = new int[50];       int[] neci = new int[50];       int[] neco = new int[50]; 
		int[]     spc = new int[50];       int[]  sci = new int[50];       int[]  sco = new int[50]; 
		double[]  epc = new double[50]; double[]  eci = new double[50]; double[]  eco = new double[50];  					
		double[]  bpc = new double[50]; double[]  bci = new double[50]; double[]  bco = new double[50];
					
		Vector3D[] r1 = new Vector3D[50]; Vector3[] c1 = new Vector3[50];
		Vector3D[] r4 = new Vector3D[50]; Vector3[] c4 = new Vector3[50]; 
		Vector3D[] r7 = new Vector3D[50]; Vector3[] c7 = new Vector3[50];  
		
		System.out.println(trsec+" "+npart+" "+GEN.get(0));
		
		if (npart>=1) {
			int npp = 0, ipp = 0;
			for (DetectorParticle dp : par) {
			    if(dp.getSector(DetectorType.ECAL)==2 && dp.getPid()==2112) {
			    	if(npp==0) p1 = dp; //Neutron 1
			    	if(npp==1) p2 = dp; //Neutron 2
			        for(int iresp = 0; iresp < cal.size(); iresp++){
			        	CalorimeterResponse dr = (CalorimeterResponse)cal.get(iresp); 			              
			    		int lay = dr.getDescriptor().getLayer(); 			    		  
			    		if (dr.getAssociation(0)==ipp && dr.getDescriptor().getType()==DetectorType.ECAL) {
//			            for (DetectorResponse dr : dp.getDetectorResponses()) { 			            	
//			            	int lay = dr.getDescriptor().getType()==DetectorType.ECAL ? dr.getDescriptor().getLayer():0; 	
			    			double dre=dr.getEnergy(),drt = dr.getTime()-stt, drb=dr.getPath()/drt/29.97; int drs=dr.getStatus(); 
			    			System.out.println(lay+" "+dr.getPath()+" "+dr.getTime()+" "+stt);
			    			Vector3D drp=dr.getPosition(); Vector3 drc=dr.getCoordUVW(); 
			            	if(lay==1) { npc[0]++ ; npc[npp+1]++;epc[npp+1]=dre;r1[npp+1]=drp;c1[npp+1]=drc;spc[npp+1]=drs;bpc[npp+1]=drb; }    					
			            	if(lay==4) {neci[0]++ ;neci[npp+1]++;eci[npp+1]=dre;r4[npp+1]=drp;c4[npp+1]=drc;sci[npp+1]=drs;bci[npp+1]=drb; }    					
			            	if(lay==7) {neco[0]++ ;neco[npp+1]++;eco[npp+1]=dre;r7[npp+1]=drp;c7[npp+1]=drc;sco[npp+1]=drs;bco[npp+1]=drb; }
			            }
			         } 			    		   
			         npp++;
			    }
			    ipp++;
			}
			
			System.out.println(npc[0]+" "+neci[0]+" "+neco[0]);
			
			boolean n100 = npc[0]==1 && neci[0]==0 && neco[0]==0; 
			boolean n010 = npc[0]==0 && neci[0]==1 && neco[0]==0;
			boolean n001 = npc[0]==0 && neci[0]==0 && neco[1]==1;
			
			if(n100) dgm.fill("n100",GEN.get(0), bpc[1]);
			if(n010) dgm.fill("n010",GEN.get(0), bci[1]);
			if(n001) dgm.fill("n001",GEN.get(0), bco[1]);
			
		}   	       
    }
    
    public void geteff() {
    	dgm.geteff("g21", "h2", "h1");
    	dgm.geteff("g31", "h3", "h1");
    }
   
}
