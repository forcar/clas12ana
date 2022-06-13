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
	
	List<DetectorParticle> par;
	List<DetectorResponse> cal;	
	int trsec = -1;
	String tit = null;
	double ethresh = 0.3;	
	
	static double refE, refP, refTH, refPH;
	static int trSEC=5, trPID=-211, mcSEC=2, mcPID= 2112;
	
    public ECmcn(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("DEMO","NEUTRONS","CLUSTERS","RECGEN","EFFICIENCY");

        this.use123Buttons(true);
        this.useZSliderPane(true);
        useECEnginePane(true);

        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println(getDetectorName()+".localinit()");
        initEBMCE();
        tl.setFitData(Fits);          
    }
        
   public void initEBMCE() {        
       System.out.println(getDetectorName()+".initEBMCE()");          
       ebmce.getCCDB(10);
       ebmce.setGeom("2.5");
       ebmce.setMCpid(mcPID);
       ebmce.setMCsec(mcSEC);
       ebmce.isMC = true;
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
    
    
    public void openOutput(String file) {
    	System.out.println(getDetectorName()+".openOutput("+file+")");
    	writer = new HipoDataSync();
        writer.open(file);    	
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
    		dgm.add("NEUTRONS", 3,4,0,st,getRunNumber());
//    		dgm.makeH2("n00",100,0.2,1.3,9,0.5,9.5,-1,"","Neutron #beta","Neutrons");
//    		dgm.makeH2("n01",100,0,2,    9,0.5,9.5,-1,"","Neutron REC Energy (GeV)","Neutrons");
    		dgm.makeH2("n100", 50,0,6,   50,-0.15,0.1,0,"pc=1 eci=0 eco=0","Neutron GEN Momentum (GeV)","PCAL #Delta#beta");
    		dgm.makeH2("n010", 50,0,6,   50,-0.15,0.1,0,"pc=0 eci=1 eco=0","Neutron GEN Momentum (GeV)","ECIN #Delta#beta");
    		dgm.makeH2("n001", 50,0,6,   50,-0.15,0.1,0,"pc=0 eci=0 eco=1","Neutron GEN Momentum (GeV)","ECOU #Delta#beta");
    		dgm.makeH2("n+00", 50,0,6,   50,-0.15,0.1,0,"pc>1 eci=0 eco=0","Neutron GEN Momentum (GeV)","PCAL #Delta#beta");
    		dgm.makeH2("n0+0", 50,0,6,   50,-0.15,0.1,0,"pc=0 eci>1 eco=0","Neutron GEN Momentum (GeV)","ECIN #Delta#beta");
    		dgm.makeH2("n00+", 50,0,6,   50,-0.15,0.1,0,"pc=0 eci=0 eco>1","Neutron GEN Momentum (GeV)","ECOU #Delta#beta");
    		dgm.makeH2("n1110",50,0,6,   50,-0.15,0.1,0,"pc=1 eci=1 eco=1","Neutron GEN Momentum (GeV)","PCAL #Delta#beta");
    		dgm.makeH2("n1111",50,0,6,   50,-0.15,0.1,0,"pc=1 eci=1 eco=1","Neutron GEN Momentum (GeV)","ECIN #Delta#beta");
    		dgm.makeH2("n1112",50,0,6,   50,-0.15,0.1,0,"pc=1 eci=1 eco=1","Neutron GEN Momentum (GeV)","ECOU #Delta#beta");
    		dgm.makeH2("pgr1", 50,0,6,   50,-0.6,0.6, 0,"pc=1 eci=1 eco=1","Neutron GEN Momentum (GeV)","PCAL #DeltaP/P" );
    		dgm.makeH2("pgr2", 50,0,6,   50,-0.6,0.6, 0,"pc=1 eci=1 eco=1","Neutron GEN Momentum (GeV)","ECIN #DeltaP/P" );
    		dgm.makeH2("pgr3", 50,0,6,   50,-0.6,0.6, 0,"pc=1 eci=1 eco=1","Neutron GEN Momentum (GeV)","ECOU #DeltaP/P" );
    	}
    }
    
    @Override
    public void processEvent(DataEvent de) {
		
		boolean goodev = ebmce.readMC(de) && ebmce.pmc.size()==1;
		
        if (goodev) { 
        	
        	refE  = ebmce.pmc.get(0).e(); 
        	refP  = ebmce.pmc.get(0).p();
        	refTH = Math.toDegrees(ebmce.pmc.get(0).theta());
        	refPH = Math.toDegrees(ebmce.pmc.get(0).phi());
        	
        	dgm.fill("h1",refP);  
        	
	    	if (dropBanks) dropBanks(de);  //drop ECAL banks and re-run ECEngine            
        	if(!ebmce.processDataEvent(de)) return;  
        	
            par = ebmce.eb.getEvent().getParticles(); //REC::Particle
        	cal = ebmce.eb.getEvent().getCalorimeterResponseList();  //REC::Calorimeter
       	
        	trsec = -1; 
        	for (DetectorParticle dp : par) { //find trigger particle and sector
        		int pid = dp.getPid(); int sec = dp.getSector(DetectorType.ECAL);
        		if(trsec==-1 && sec==trSEC && pid==trPID) trsec=sec;
        	}
        	
        	if(trsec==-1) return; // reject events with missing trigger particle  
        	
        	processDEMO(); 
        	processNEUTRONS();
        	
            if(dumpFiles) {ebmce.getRECBanks(de,ebmce.eb); writer.writeEvent(de);}     
        }
    }
    
    public void processDEMO() {
    	
	    boolean n1=true, n2=true; float en=0; int np=0;
	    boolean goodpc=false, goodeci=false, goodeco=false;
		
	    for (DetectorParticle dp : par) {         			
	    	if(dp.getSector(DetectorType.ECAL)==mcSEC) { np++; 
	    		for (DetectorResponse ds : dp.getDetectorResponses()) { 
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
		}
	    
        if(en>0)                     {dgm.fill("h6", refP,en); dgm.fill("h8",refP, np); dgm.fill("h10",1e3*(refE-0.93957));}
        if(goodpc&&goodeci&&goodeco) {dgm.fill("h7", refP,en); dgm.fill("h9",refP, np); dgm.fill("h11",1e3*(refE-0.93957));}    	
    }
    
    public boolean isNeutral(DetectorParticle dp) {
    	return dp.getPid()==2112 || dp.getPid()==22;
    }
    
    public float bp(double rm, double b) {
    	return (float) (rm*b/Math.sqrt(1-b*b));
    }
    
    public void processNEUTRONS() {
    	
    	int npart=0;	
        float stt = ebmce.starttime;
        
    	DetectorParticle p1 = new DetectorParticle();
        
        for (DetectorParticle dp : par) if(dp.getSector(DetectorType.ECAL)==mcSEC && isNeutral(dp)) npart++;

 		double dist=0, du=0;
 		
		int[]     npc = new int[50];       int[] neci = new int[50];       int[] neco = new int[50]; 
		int[]     spc = new int[50];       int[]  sci = new int[50];       int[]  sco = new int[50]; 
		double[]  epc = new double[50]; double[]  eci = new double[50]; double[]  eco = new double[50];  					
		double[]  bpc = new double[50]; double[]  bci = new double[50]; double[]  bco = new double[50];
					
		Vector3D[] r1 = new Vector3D[50]; Vector3[] c1 = new Vector3[50];
		Vector3D[] r4 = new Vector3D[50]; Vector3[] c4 = new Vector3[50]; 
		Vector3D[] r7 = new Vector3D[50]; Vector3[] c7 = new Vector3[50];  
		
		if (npart>=1) { // number of EB ID=22,2112 > 0
			
			int npp = 0, ipp = 0;
			
			for (DetectorParticle dp : par) {
			    if(dp.getSector(DetectorType.ECAL)==mcSEC && isNeutral(dp)) {
			    	if(npp==0) p1 = dp; //Neutron 1
 			        for(DetectorResponse dr : cal){ CalorimeterResponse r = (CalorimeterResponse) dr; 
			        	if (r.getAssociation(0)==ipp) { //assume r.getNAssociations=1			        				              
			        		int lay = r.getDescriptor().getLayer(); 			    		  
			    			double dre=dr.getEnergy(),drt = dr.getTime()-stt, drb=dr.getPath()/drt/29.97; 
			    			int drs=dr.getStatus(); 			    			
			    			Vector3D drp=dr.getPosition(); Vector3 drc=r.getCoordUVW(); 
			            	if(lay==1) { npc[0]++ ; npc[npp+1]++;epc[npp+1]=dre;r1[npp+1]=drp;c1[npp+1]=drc;spc[npp+1]=drs;bpc[npp+1]=drb; }    					
			            	if(lay==4) {neci[0]++ ;neci[npp+1]++;eci[npp+1]=dre;r4[npp+1]=drp;c4[npp+1]=drc;sci[npp+1]=drs;bci[npp+1]=drb; }    					
			            	if(lay==7) {neco[0]++ ;neco[npp+1]++;eco[npp+1]=dre;r7[npp+1]=drp;c7[npp+1]=drc;sco[npp+1]=drs;bco[npp+1]=drb; }
			            }
			         } 			    		   
			         npp++;
			    }
			    ipp++;
			}
			
			float bmc = (float) Math.sqrt(1/(1+0.93957*0.93957/refP/refP));
			
			boolean n100 = npc[0]==1 && neci[0]==0 && neco[0]==0; 
			boolean n010 = npc[0]==0 && neci[0]==1 && neco[0]==0;
			boolean n001 = npc[0]==0 && neci[0]==0 && neco[0]==1;
			boolean np00 = npc[0]>1  && neci[0]==0 && neco[0]==0; 
			boolean n0p0 = npc[0]==0 && neci[0]>1  && neco[0]==0;
			boolean n00p = npc[0]==0 && neci[0]==0 && neco[0]>1;
			boolean n111 = npc[0]==1 && neci[0]==1 && neco[0]==1;
			
			float ppc = bp(0.93957,bpc[1]);
			float pci = bp(0.93957,bci[1]);
			float pco = bp(0.93957,bco[1]);
			
			if(n100) dgm.fill("n100",refP, bpc[1]-bmc);
			if(n010) dgm.fill("n010",refP, bci[1]-bmc);
			if(n001) dgm.fill("n001",refP, bco[1]-bmc);
			if(np00) dgm.fill("n+00",refP, bpc[1]-bmc);
			if(n0p0) dgm.fill("n0+0",refP, bci[1]-bmc);
			if(n00p) dgm.fill("n00+",refP, bco[1]-bmc);
			if(n111) dgm.fill("n1110",refP,bpc[1]-bmc);
			if(n111) dgm.fill("n1111",refP,bci[1]-bmc);
			if(n111) dgm.fill("n1112",refP,bco[1]-bmc);
			if(n111) dgm.fill("pgr1", refP,(ppc-refP)/refP);
			if(n111) dgm.fill("pgr2", refP,(pci-refP)/refP);
			if(n111) dgm.fill("pgr3", refP,(pco-refP)/refP);
			
		}   	       
    }
    
    public void geteff() {
    	dgm.geteff("g21", "h2", "h1");
    	dgm.geteff("g31", "h3", "h1");
    }
   
}
