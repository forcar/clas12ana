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
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;

public class ECmc2 extends DetectorMonitor {
	
	Event          ev = new Event();
	EBMCEngine  ebmce = new EBMCEngine();
	
	List<Float>           GEN   = new ArrayList<Float>();
	List<Float>           REC   = new ArrayList<Float>();   
	List<Float>           GENPZ = new ArrayList<Float>();
	List<Float>           RECPZ = new ArrayList<Float>(); 
	
	HipoDataSync       writer = null;		
	List<Particle>       phot = new ArrayList<Particle>();
	List<Particle>       neut = new ArrayList<Particle>();
	
	String tit = null;
	double ethresh = 0.3;
	int MPpid = 22;
	
    public ECmc2(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("PHOTONS","CLUSTERS","GENREC","EFFICIENCY","PIZERO");

        this.use123Buttons(true);
        this.useZSliderPane(true);

        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECmc2.localinit()");
                
        engine.setIsMC(true);       
        engine.setGeomVariation("rga_fall2018");  
        engine.setDebug(false);
        setEngineConfig("test");
        
        ebmce.getCCDB(10); 
        ebmce.setGeom("2.5");
        ebmce.setGoodPhotons(12);
        ebmce.setMCpid(MPpid);
        ebmce.isMC = true;
        
        tl.setFitData(Fits);
    }
    
    public void localclear() {
    	System.out.println("ECmc2.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(1);
    }  
    
    @Override
    public void createHistos(int run) {
    	System.out.println("ECmc2:createHistos("+run+")");
        writer = new HipoDataSync();
        engine.setUnsharedEnergyCut(6);
        if(dumpFiles) writer.open("/Users/colesmith/CLAS12ANA/ECmc2/photon_mc2.hipo");
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;  
    	createPHOTONS(0);
    	createPHOTONS(1);
    	createPHOTONS(2);
    	createPHOTONS(3);
    	createCLUSTERS(0);
    	createCLUSTERS(1);
    	createCLUSTERS(2);
    	createGENREC(0);
    	createGENREC(1);
    	createGENREC(2);
       	createEFFICIENCY(0);
       	createEFFICIENCY(1);
       	createPIZERO(0);
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
        plot("EFFICIENCY");  
    }
       
    public void plotMCHistos() {  
    	plot("PHOTONS");
        plot("CLUSTERS");
        plot("GENREC");       
        plot("PIZERO");
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
	
    @Override
    public void plotEvent(DataEvent de) {
       // analyze();         
    }
    
    public void createPHOTONS(int st) {
    	
    	switch (st) {
    	case 0:
    		dgm.add("PHOTONS", 4,2,0,st,getRunNumber());
    		dgm.makeH2("p00",100,0.9,1.3,9,0.5,9.5,-1,"","#gamma1 #beta","Photons");
    		dgm.makeH2("p01",100,0.9,1.3,9,0.5,9.5,-1,"","#gamma2 #beta","Photons");
    		dgm.makeH2("p02",100,0,5,9,0.5,9.5,-1,    "","#gamma1 Energy (GeV)","Photons");
    		dgm.makeH2("p03",100,0,5,9,0.5,9.5,-1,    "","#gamma2 Energy (GeV)","Photons");
            dgm.makeH2("p04",50,0,4,50,0,4,    -1,    "","E1 / E2 (PCAL)","E1 / E2 (ECIN)");
            dgm.makeH2("p05",50,0,4,50,0,4,    -1,    "","E2 / E1 (PCAL)","E1 / E2 (ECIN)");
            dgm.makeH2("p06",80,0,180,50,0,2,   1,  "pc=2 ec=2","#phi_#gamma_#gamma (deg)","E2/E #gamma2");            
            dgm.makeH2("p07",80,0,180,50,0,2,   1,  "pc>1 ec>1","#phi_#gamma_#gamma (deg)","E2/E #gamma2");            
    		break;
    	case 1:
    		dgm.add("PHOTONS",2,2,0,st,getRunNumber());
    		dgm.makeH2("p10", 50, 2, 9, 50,  0, 1, -1, "GEN", "Pizero Energy (GeV)", "X");
    		dgm.makeH2("p11", 50, 0, 1, 50,  0, 10,-1, "GEN", "X", "OPA (DEG)");    		
    		dgm.makeH2("p12", 50, 2, 9, 50,  0, 1, -1, "REC", "Pizero Energy (GeV)", "X");
    		dgm.makeH2("p13", 50, 0, 1, 50,  0, 10,-1, "REC", "X", "OPA (DEG)");    
    		break;
    	case 2:
    		dgm.add("PHOTONS",7,3,0,st,getRunNumber());
    		for (int il=0; il<3; il++) {
    			for (int ic=0; ic<7; ic++) {
    				String txt = layer[il]+" cluster status = "+ic;
    	            dgm.makeH2("p2"+il+ic,80,0,180,50,0,2,   1,  "pc=2 ec=2 "+txt,"#phi_#gamma_#gamma (deg)","E2/E #gamma2");               				
    			}
    		} 
    		break;
    	case 3:
    		dgm.add("PHOTONS",7,3,0,st,getRunNumber());
    		for (int il=0; il<3; il++) {
    			for (int ic=0; ic<7; ic++) {
    				String txt = layer[il]+" cluster status = "+ic;
    	            dgm.makeH2("p3"+il+ic,80,0,180,80,-0.2,0.2,   0,  "pc=2 ec=2 "+txt,"#phi_#gamma_#gamma (deg)","#Delta #beta");               				
    			}
    		}    		
    	}
    }
    
    public void createCLUSTERS(int st) {
    	String test;      
    	switch (st) {        
        case 0:  
        	dgm.add("CLUSTERS", 5,3,0,st,getRunNumber()); 
            dgm.makeH2("c003",16,-0.5,15.5, 80,0.4,1.6,1,"",          "Photons",     "REC (E1+E2)/2E");
            dgm.makeH2("c004",16,-0.5,15.5, 80,0.4,1.6,1,"OPA>1.95",  "Photons",     "REC (E1+E2)/2E");
            dgm.makeH2("c020",50,0,5.2,     16,-0.5,15.5,-1,"",         "GEN #theta12 (deg)","Photons"); 
            dgm.makeH2("c021",50,0,5.2,     16,-0.5,15.5,-1,"pc=2 ec=0","GEN #theta12 (deg)","Photons"); 
            dgm.makeH2("c022",50,0,5.2,     16,-0.5,15.5,-1,"pc=2 ec=1","GEN #theta12 (deg)","Photons"); 
            dgm.makeH2("c023",50,0,5.2,     16,-0.5,15.5,-1,"pc=2 ec>1","GEN #theta12 (deg)","Photons"); 
            dgm.makeH2("c024",50,0,5.2,     16,-0.5,15.5,-1,"pc=1 ec>1","GEN #theta12 (deg)","Photons");  
            break;
        case 1:
        	dgm.add("CLUSTERS", 7,3,0,st,getRunNumber()); test = "pc=2 ec=2";
            dgm.makeH2("c100",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (PCAL)","#beta #gamma1");
            dgm.makeH2("c101",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (PCAL)","#beta #gamma2");
            dgm.makeH2("c102",80,0,4,  80,-0.2,0.2, 0,test,"E1 / E2 (PCAL)",        "PCAL #Delta #beta");
            dgm.makeH2("c103",80,0,5.4,80,-0.2,0.2, 0,test,"Opening Angle (deg)",   "PCAL #Delta #beta");
            dgm.makeH2("c104",80,0,180,80,-0.2,0.2, 0,test,"#phi_#gamma_#gamma (deg)","PCAL #Delta #beta");
            dgm.makeH2("c104a",7,0,7,  80,-0.2,0.2, 0,test,"Cluster1 Status","PCAL #Delta #beta");
            dgm.makeH2("c104b",80,0,180,7,0,7,      1,test,"#phi_#gamma_#gamma (deg)","Cluster1 Status");
            
            dgm.makeH2("c105",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECIN)","#beta #gamma1");
            dgm.makeH2("c106",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECIN)","#beta #gamma2");
            dgm.makeH2("c107",80,0,4,  80,-0.2,0.2, 0,test,"E1 / E2 (ECIN)",        "ECIN #Delta #beta");
            dgm.makeH2("c108",80,0,5.4,80,-0.2,0.2, 0,test,"Opening Angle (deg)",   "ECIN #Delta #beta");
            dgm.makeH2("c109",80,0,180,80,-0.2,0.2, 0,test,"#phi_#gamma_#gamma (deg)","ECIN #Delta #beta");
            dgm.makeH2("c109a",7,0,7,  80,-0.2,0.2, 0,test,"Cluster1 Status","ECIN #Delta #beta");
            dgm.makeH2("c109b",80,0,180,7,0,7,      1,test,"#phi_#gamma_#gamma (deg)","Cluster1 Status");

            dgm.makeH2("c110",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECOU)","#beta #gamma1 (ns)");
            dgm.makeH2("c111",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECOU)","#beta #gamma2 (ns)");
            dgm.makeH2("c112",80,0,4,  80,-0.2,0.2, 0,test,"E1 / E2 (ECOU)",        "ECOU #Delta #beta");
            dgm.makeH2("c113",80,0,5.4,80,-0.2,0.2, 0,test,"Opening Angle (deg)",   "ECOU #Delta #beta");
            dgm.makeH2("c114",80,0,180,80,-0.2,0.2, 0,test,"#phi_#gamma_#gamma (deg)","ECOU #Delta #beta");
            dgm.makeH2("c114a",7,0,7,  80,-0.2,0.2, 0,test,"Cluster1 Status","ECOU #Delta #beta");
            dgm.makeH2("c114b",80,0,180,7,0,7,      1,test,"#phi_#gamma_#gamma (deg)","Cluster1 Status");
            break;
        case 2:
        	dgm.add("CLUSTERS", 7,3,0,st,getRunNumber()); test = "pc>1 ec>1";
            dgm.makeH2("c200",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (PCAL)","#beta #gamma1");
            dgm.makeH2("c201",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (PCAL)","#beta #gamma2");
            dgm.makeH2("c202",80,0,4,  80,-0.2,0.2, 0,test,"E1 / E2 (PCAL)",        "PCAL #Delta #beta");
            dgm.makeH2("c203",80,0,5.4,80,-0.2,0.2, 0,test,"Opening Angle (deg)",   "PCAL #Delta #beta");
            dgm.makeH2("c204",80,0,180,80,-0.2,0.2, 0,test,"#phi_#gamma_#gamma (deg)","PCAL #Delta #beta");
            dgm.makeH2("c204a",7,0,7,  80,-0.2,0.2, 0,test,"Cluster1 Status","PCAL #Delta #beta");          
            dgm.makeH2("c204b",80,0,180,7,0,7,      1,test,"#phi_#gamma_#gamma (deg)","Cluster1 Status");
            
            dgm.makeH2("c205",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECIN)","#beta #gamma1");
            dgm.makeH2("c206",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECIN)","#beta #gamma2");
            dgm.makeH2("c207",80,0,4,  80,-0.2,0.2, 0,test,"E1 / E2 (ECIN)",        "ECIN #Delta #beta");
            dgm.makeH2("c208",80,0,5.4,80,-0.2,0.2, 0,test,"Opening Angle (deg)",   "ECIN #Delta #beta");
            dgm.makeH2("c209",80,0,180,80,-0.2,0.2, 0,test,"#phi_#gamma_#gamma (deg)","ECIN #Delta #beta");
            dgm.makeH2("c209a",7,0,7,  80,-0.2,0.2, 0,test,"Cluster1 Status","ECIN #Delta #beta");          
            dgm.makeH2("c209b",80,0,180,7,0,7,      1,test,"#phi_#gamma_#gamma (deg)","Cluster1 Status");

            dgm.makeH2("c210",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECOU)","#beta #gamma1 (ns)");
            dgm.makeH2("c211",80,0,1,  80, 0.8,1.2, 1,test,"E #gamma1 (ECOU)","#beta #gamma2 (ns)");
            dgm.makeH2("c212",80,0,4,  80,-0.2,0.2, 0,test,"E1 / E2 (ECOU)",        "ECOU #Delta #beta");
            dgm.makeH2("c213",80,0,5.4,80,-0.2,0.2, 0,test,"Opening Angle (deg)",   "ECOU #Delta #beta");
            dgm.makeH2("c214",80,0,180,80,-0.2,0.2, 0,test,"#phi_#gamma_#gamma (deg)","ECOU #Delta #beta");
            dgm.makeH2("c214a",7,0,7,  80,-0.2,0.2, 0,test,"Cluster1 Status","ECOU #Delta #beta");                   
            dgm.makeH2("c214b",80,0,180,7,0,7,      1,test,"#phi_#gamma_#gamma (deg)","Cluster1 Status");
    	}       
    }
    
    public void createGENREC(int st) {
	  
    	switch (st) {        
    	case 1:
    		dgm.add("GENREC",7,7,0,st,getRunNumber());
                                                 tit = "pc>1 ec=0";
            dgm.makeH2("d00",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1");
            dgm.makeH2("d10",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2");
            dgm.makeH2("d02",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1"); 
            dgm.makeH2("d12",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2"); 
            dgm.makeH2("d03",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "E1/E #gamma1");  
            dgm.makeH2("d13",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "E2/E #gamma2");  
            dgm.makeH2("d14",7,0,7,     50,0,2,1,tit+" opa>2.5","Cluster1 Status",        "E2/E #gamma2");  
                                                 tit = "pc>1 ec=1";
            dgm.makeH2("d20",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1"); 
            dgm.makeH2("d30",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2");  
            dgm.makeH2("d22",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1");   
            dgm.makeH2("d32",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2");  
            dgm.makeH2("d23",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "E1/E #gamma1");   
            dgm.makeH2("d33",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "E2/E #gamma2");
            dgm.makeH2("d34",7,0,7,     50,0,2,1,tit+" opa>2.5","Cluster1 Status",        "E2/E #gamma2");  
                                                 tit = "pc>1 ec>1";
            dgm.makeH2("d40",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1"); 
            dgm.makeH2("d50",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2"); 
            dgm.makeH2("d42",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1"); 
            dgm.makeH2("d52",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2");  
            dgm.makeH2("d43",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "E1/E #gamma1");  
            dgm.makeH2("d53",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "E2/E #gamma2");              
            dgm.makeH2("d54",7,0,7,     50,0,2,1,tit+" opa>2.5","Cluster1 Status",        "E2/E #gamma2");  
                                                 tit = "pc=2 ec=1"; 
            dgm.makeH2("d60",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d70",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d62",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d72",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d63",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "(E1+E2)/2E");  
            dgm.makeH2("d73",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "SQRT(E1*E2)/E");  
            dgm.makeH2("d74",7,0,7,     50,0,2,1,tit+" opa>2.5","Cluster1 Status",        "SQRT(E1*E2)/E");  
                                                 tit = "pc>1 ec>1";
            dgm.makeH2("d80",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d90",50,10,20,  50,0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d82",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d92",50,0,70,   50,0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d83",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "(E1+E2)/2E");  
            dgm.makeH2("d93",50,0,5.2,  50,0,2,1,tit,"GEN #theta12 (deg)",     "SQRT(E1*E2)/E");  
            dgm.makeH2("d94",7,0,7,     50,0,2,1,tit+" opa>2.5","Cluster1 Status",        "SQRT(E1*E2)/E");  
                                                 tit = "pc=2 ec=2";
            dgm.makeH2("d100",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d110",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d102",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d112",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d103",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg)",     "(E1+E2)/2E");  
            dgm.makeH2("d113",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg)",     "SQRT(E1*E2)/E");  
            dgm.makeH2("d114",7,0,7,    50,0,2,1,tit+" opa>2.5","Cluster1 Status",        "SQRT(E1*E2)/E");  
            									 tit = "pc=2 ec=2";
            dgm.makeH2("d120",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E1/E #gamma1"); 
            dgm.makeH2("d130",50,10,20, 50,0,2,1,tit,"GEN #theta #gamma2 (deg)","E2/E #gamma2"); 
            dgm.makeH2("d122",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "E1/E #gamma1"); 
            dgm.makeH2("d132",50,0,70,  50,0,2,1,tit,"Distance (cm)",           "E2/E #gamma2");  
            dgm.makeH2("d123",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg)",     "E1/E #gamma1");  
            dgm.makeH2("d133",50,0,5.2, 50,0,2,1,tit,"GEN #theta12 (deg)",     "E2/E #gamma2");              
            dgm.makeH2("d134",7,0,7,    50,0,2,1,tit+" opa>2.5","Cluster1 Status",        "E2/E #gamma2");  

            break;
    	case 0:
    		dgm.add("GENREC",6,4,0,st,getRunNumber()); int y1=-2, y2=2;
                                                   tit = "pc>1 ec=0";
            dgm.makeH2("dt00",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1");
            dgm.makeH2("dt10",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2");
            dgm.makeH2("dt02",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1"); 
            dgm.makeH2("dt12",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2"); 
            dgm.makeH2("dt03",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta #gamma1");  
            dgm.makeH2("dt13",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta #gamma2");  
                                                   tit = "pc>1 ec=1";
            dgm.makeH2("dt20",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1"); 
            dgm.makeH2("dt30",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2");  
            dgm.makeH2("dt22",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1");   
            dgm.makeH2("dt32",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2");  
            dgm.makeH2("dt23",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta #gamma1");   
            dgm.makeH2("dt33",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta #gamma2");
                                                   tit = "pc>1 ec>1";
            dgm.makeH2("dt40",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma1"); 
            dgm.makeH2("dt50",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta #gamma2"); 
            dgm.makeH2("dt42",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma1"); 
            dgm.makeH2("dt52",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta #gamma2");  
            dgm.makeH2("dt43",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta #gamma1");  
            dgm.makeH2("dt53",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta #gamma2");                                                    tit = "pc=1 ec>1"; 
                                                   tit = "pc>1 ec>0";            
            dgm.makeH2("dt60",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta12 (deg)"); 
            dgm.makeH2("dt62",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta12 (deg)");    
            dgm.makeH2("dt63",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta12 (deg)");   
                                                   tit = "pc=1 ec>1";             
            dgm.makeH2("dt70",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#theta12 (deg)");   
            dgm.makeH2("dt72",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#theta12 (deg)");  
            dgm.makeH2("dt73",50,0,5.5, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#theta12 (deg)");  
            break;
    	case 2:
    		dgm.add("GENREC",6,4,0,st,getRunNumber()); y1=-3; y2=3;
                                                   tit = "pc>1 ec=0";
            dgm.makeH2("df00",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1");
            dgm.makeH2("df10",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2");
            dgm.makeH2("df02",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1"); 
            dgm.makeH2("df12",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2"); 
            dgm.makeH2("df03",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#phi #gamma1");  
            dgm.makeH2("df13",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#phi #gamma2");  
                                                   tit = "pc>1 ec=1";
            dgm.makeH2("df20",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1"); 
            dgm.makeH2("df30",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2");  
            dgm.makeH2("df22",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1");   
            dgm.makeH2("df32",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2");  
            dgm.makeH2("df23",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#phi #gamma1");   
            dgm.makeH2("df33",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#phi #gamma2");
                                                   tit = "pc>1 ec>1";
            dgm.makeH2("df40",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma1"); 
            dgm.makeH2("df50",50,10,20, 50,y1,y2,0,tit,"GEN #theta #gamma2 (deg)","#Delta#phi #gamma2"); 
            dgm.makeH2("df42",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma1"); 
            dgm.makeH2("df52",50,0,70,  50,y1,y2,0,tit,"Distance (cm)",           "#Delta#phi #gamma2");  
            dgm.makeH2("df43",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#phi #gamma1");  
            dgm.makeH2("df53",50,0,5.2, 50,y1,y2,0,tit,"GEN #theta12 (deg)",     "#Delta#phi #gamma2");    		            
    	}        
   	   
    }   
    
    public void createEFFICIENCY(int st) {
  
    	switch (st) {        
        case 0:          
        	dgm.add("EFFICIENCY",6,5,0,st,getRunNumber());
            dgm.makeH1("eff1",50,0,5.2,-1,"n#gamma>1",          "Opening Angle (deg)","Efficiency",1,3);
            dgm.makeH1("eff1a",50,0,5.2,-2,"pc>0 ec>1",         "Opening Angle (deg)","Efficiency",1,4);
            dgm.makeH1("eff3",50,0,5.2,-2,"pc>1 ec=1 pc>1 ec>1","Opening Angle (deg)","Efficiency",1,5);
            dgm.makeH1("eff4",50,0,5.2,-2,"pc>1 ec=2",          "Opening Angle (deg)","Efficiency",1,2);
            dgm.makeH1("eff5",50,0,5.2,-2,"pc=2 ec=2",          "Opening Angle (deg)","Efficiency",1,1);
            dgm.makeH1("eff1a"); dgm.makeH1("eff3"); 
            dgm.makeH1("eff2",50,0,5.2,-2,"pc>1 ec=1",          "Opening Angle (deg)","Efficiency",1,0);                      
            dgm.makeH1("eff4");  dgm.makeH1("eff5"); dgm.makeH1("eff2");
            dgm.makeH1("eff3dpcb",50,0,5.2,-1,"pc>1 ec>1 #Delta#beta<0.1","Opening Angle (deg)","Efficiency",1,5);
            
            dgm.makeH2("ef10",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","PC Clusters" ,"Photons"); dgm.cc("ef10",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH2("ef11",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","ECi Clusters","Photons"); dgm.cc("ef11",false,false, 0,0,0.8f,1.2f);  
            dgm.makeH2("ef12",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>1.95","ECo Clusters","Photons"); dgm.cc("ef12",false,false, 0,0,0.8f,1.2f);              
            dgm.makeH2("ef13",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","PC Clusters" ,"Photons");  dgm.cc("ef13",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH2("ef14",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","ECi Clusters","Photons");  dgm.cc("ef14",false,false, 0,0,0.8f,1.2f);  
            dgm.makeH2("ef15",16,-0.5,15.5, 16,-0.5,15.5,-1,"(E1+E2)/2E opa>2.4","ECo Clusters","Photons");  dgm.cc("ef15",false,false, 0,0,0.8f,1.2f); 
            dgm.makeH1("eff111",50,0,5.2,-1,"pc>1 ec>0",        "Opening Angle (deg)","Efficiency",1,2);
            dgm.makeH1("ef17",180,0,180,-1,"pc>1 ec>1 #Delta#beta<0.1","#phi_#gamma_#gamma (deg)","");             
            dgm.makeH1("ef18",180,0,180,-1,"pc>1 ec>1","#phi_#gamma_#gamma (deg)","");             
            dgm.makeH1("ef19",180,0,180,-1,"pc=2 ec=2","#phi_#gamma_#gamma (deg)","");             
        	dgm.makeH1("ef118",160,-80,80,-1,"","dU","",1,4);        	       	
        	dgm.makeH1("ef119",160,-80,80,-1,"","dV","",1,4);        	       	
        	dgm.makeH1("ef120",160,-80,80,-1,"","dW","",1,4);        	       	
           break;
        case 1:
        	dgm.add("EFFICIENCY",6,6,0,st,getRunNumber());          
            dgm.makeH1("h10",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
            dgm.makeH1("h1a",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
        	dgm.makeH1("h11",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,5);
        	dgm.makeH1("h111",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,5);
        	dgm.makeH1("h12",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,2);
        	dgm.makeH1("h13",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,1);
        	dgm.makeH1("h13dpcb",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,1);
        	dgm.makeH1("h14",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);
        	dgm.makeH1("h15",50,0,5.2,-1,"","Opening Angle (deg)","Efficiency",1,4);        	
        	dgm.makeH1("h16",180,0,180,-1,"","#phi_#gamma_#gamma (deg)","",1,4);        	
        	dgm.makeH1("h17",180,0,180,-1,"","#phi_#gamma_#gamma (deg)","",1,4);   
        	dgm.makeH1("h18",180,0,180,-1,"","#phi_#gamma_#gamma (deg)","",1,4);   
        	dgm.makeH1("h19",180,0,180,-1,"","#phi_#gamma_#gamma (deg)","",1,4);   
            dgm.makeH2("e00",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>1.95", "PC Clusters","Photons"); 
            dgm.makeH2("e01",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>1.95","ECi Clusters","Photons"); 
            dgm.makeH2("e02",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>1.95","ECo Clusters","Photons");
            dgm.makeH2("e03",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>2.4",  "PC Clusters","Photons"); 
            dgm.makeH2("e04",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>2.4", "ECi Clusters","Photons"); 
            dgm.makeH2("e05",16,-0.5,15.5, 16,-0.5,15.5,-1,"opa>2.3", "ECo Clusters","Photons");         	
            dgm.makeH2("e10",16,-0.5,15.5, 16,-0.5,15.5,-1,"", "PC Clusters","Photons");  
            dgm.makeH2("e11",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECi Clusters","Photons");  
            dgm.makeH2("e12",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECo Clusters","Photons");             
            dgm.makeH2("e13",16,-0.5,15.5, 16,-0.5,15.5,-1,"", "PC Clusters","Photons");  
            dgm.makeH2("e14",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECi Clusters","Photons");  
            dgm.makeH2("e15",16,-0.5,15.5, 16,-0.5,15.5,-1,"","ECo Clusters","Photons");             
    	}
    	
     }
    
    public void createPIZERO(int st) {
    	
    	switch (st) {
    	case 0:
    		dgm.add("PIZERO",5,3,0,st,getRunNumber());
    		dgm.makeH2("pz0", 50,0,10,50,0,5, -1,   "","Opening Angle (deg)", "SQRT(E1E2)");
    		dgm.makeH2("pz1", 50,4, 8,50,100,200,-1,"","Pizero Energy (GeV)","#gamma#gamma Invariant Mass (MeV)");   
    		dgm.makeH2("pz2", 50,4, 8,50,-1,1,-1,   "","Pizero Energy (GeV)","Pizero Energy Error (GeV)"); 
    		dgm.makeH2("pz3", 50,4, 8,50,-1,1,-1,   "","Pizero Energy (GeV)","Pizero Theta Error (deg)"); 
    		dgm.makeH2("pz4", 50,0,10,50,-1,1,-1,   "","Opening Angle (deg)","Pizero Energy Error (GeV)"); 
    		dgm.makeH2("pz5", 50,0,10,50,-12,12,-1,   "","Opening Angle (deg)","Invariant Mass Error (MeV)"); 
    		break;
    	case 1:
    	}
    }
    
    public List<Float> getkin(List<Particle> list) {
    	List<Float> out = new ArrayList<Float>();
    	int n=0;
    	for (Particle p : list) {
    		out.add(n++,(float) p.e());
		    out.add(n++,(float) Math.toDegrees(p.theta())); 
		    out.add(n++,(float) Math.toDegrees(p.phi()));
    	}
		return out;
    }
    
    public int getStat(int pc, int ei, int eo) {
    	if(pc==0&&ei==0) return 0;
    	if(pc==1&&ei==1) return 1;
    	if(pc==2&&ei==2) return 2;
    	if(pc==3&&ei==3) return 3;
    	if(pc==4&&ei==4) return 4;
    	if(pc==5&&ei==5) return 5;
    	if(pc==6&&ei==6) return 6;
    	return -1;
    }
    
    @Override
    public void processEvent(DataEvent de) {
		
    	GEN.clear();  REC.clear(); GENPZ.clear(); RECPZ.clear();
   	
    	DetectorParticle p1 = new DetectorParticle();
    	DetectorParticle p2 = new DetectorParticle();
		
		int npart=0;		
		boolean correct=false;
		
		boolean goodev = ebmce.readMC(de) && ebmce.pmc.size()==2;
		
		if (!goodev) return;
        	
        GEN   = getkin(ebmce.pmc);	
        	        	
        GENPZ = ebmce.getPizeroKinematics(ebmce.pmc); float opa = GENPZ.get(3); float x = GENPZ.get(4);	    	
	    dgm.fill("p10",GENPZ.get(0),GENPZ.get(4)); dgm.fill("p11",GENPZ.get(4),GENPZ.get(3));	      
	    	
	    if (dropBanks) dropBanks(de);  //drop ECAL banks and re-run ECEngine   
        if(!ebmce.processDataEvent(de)) return;
        	
        float stt = ebmce.starttime;
        	
        List<DetectorParticle> par = ebmce.eb.getEvent().getParticles();  
        List<DetectorResponse> cal = ebmce.eb.getEvent().getCalorimeterResponseList(); 

        if(dumpFiles) {ebmce.getRECBanks(de,ebmce.eb); writer.writeEvent(de);}
       
        int trsec = -1; int trpid = -211;
        for (DetectorParticle dp: par) { //find sector of trpid
        	int pid = dp.getPid(), sec = dp.getSector(DetectorType.ECAL);
        	if(trsec==-1 && sec>0 && pid==trpid) trsec=sec;
        }
        
        if(trsec==-1) return;
        	
        if (dbgAnalyzer) {
        	System.out.println(" ");
        	System.out.println(getEventNumber()+" "+cal.size()+" "+par.size());
        	
        	for (DetectorResponse drr : cal) {
        		CalorimeterResponse dr = (CalorimeterResponse) drr;
        		System.out.println("Response "+dr.getAssociation()+" "+dr.getDescriptor().getType()+" "+dr.getDescriptor().getSector()+" "+dr.getDescriptor().getLayer()+" "
        	                      +par.get(dr.getAssociation()).getPid());
        	}
        	
        	int nnn=0;
        	for (DetectorParticle dp : par) {
        		System.out.println("Particle "+nnn+"  "+dp.getSector(DetectorType.ECAL)+" "+dp.getEnergy(DetectorType.ECAL));nnn++;
        		for (DetectorResponse dr : dp.getDetectorResponses()) {
              	  System.out.println(dr.getAssociation()+" "+dr.getDescriptor().getType()+" "+dr.getDescriptor().getLayer());       			
        		}
        	}
        } //dbgAnalyzer
        	
        List<Particle> plist = new ArrayList<Particle>(); 
        	
        for (DetectorParticle dp : par) { // make list of neutral Particle objects 
        	if(dbgAnalyzer) {
        		System.out.println("Plist "+trsec+" "+dp.getSector(DetectorType.ECAL)+" "+ebmce.hasStartTime+" "+dp.getPid()+" "+dp.getBeta()+" "+dp.getEnergy(DetectorType.ECAL));
        	}
		    if(dp.getSector(DetectorType.ECAL)!=trsec && dp.getPid()==MPpid) { npart++;
		    	if(!ebmce.hasStartTime && dp.getPid()==2112) {// this repairs zero momentum neutrons from non-PCAL seeded neutrals
	 				double e = dp.getEnergy(DetectorType.ECAL)/ebmce.getSF(dp); 		
			    	Vector3D vec = new Vector3D() ; vec.copy(dp.getHit(DetectorType.ECAL).getPosition()); vec.unit(); 			    		
 			    	dp.vector().add(new Vector3(e*vec.x(),e*vec.y(),e*vec.z())); //track energy for neutrals in DetectorParticle
 			    	dp.setPid(MPpid);
 			    }
		    	//SF corrected Particle energy from DetectorParticle
			    Particle p = dp.getPhysicsParticle(MPpid); p.setProperty("beta",dp.getBeta()); plist.add(p);
		    }
 		}
	    	
        //gamma-gamma phi angle
        int mcsec = 2; // sector of MC generated neutron
    	Point3D point1 = new Point3D((float)ebmce.pmv.get(0).x(),(float)ebmce.pmv.get(0).y(),(float)ebmce.pmv.get(0).z());
    	Point3D point2 = new Point3D((float)ebmce.pmv.get(1).x(),(float)ebmce.pmv.get(1).y(),(float)ebmce.pmv.get(1).z());
    	point1.rotateZ(Math.toRadians(-60*(mcsec-1))); point2.rotateZ(Math.toRadians(-60*(mcsec-1)));
        point1.rotateY(Math.toRadians(-25));           point2.rotateY(Math.toRadians(-25));      
        Vector3 vv1 = new Vector3(point1.x(),point1.y(),point1.z()); Vector3 vv2 = new Vector3(point2.x(),point2.y(),point2.z());
        Vector3 vv12 = vv1.cross(vv2);
        double ggp = Math.toDegrees(Math.atan2(vv12.y(),vv12.x()));
        if(ggp<0) ggp=-ggp;
        ggp=ggp-90;
        if(ggp<0) ggp=ggp+180;
/*        
        Vector3 ggc = ebmce.pmv.get(0).cross(ebmce.pmv.get(1));
        double ggp = Math.toDegrees(Math.atan2(ggc.y(),ggc.x()));
        if(ggp<0) ggp=-ggp;
        ggp=ggp-90;
        if(ggp<0) ggp=ggp+180;
*/        
       	dgm.fill("c020",opa,npart); 
        dgm.fill("h10",opa);
       	    
       	if(opa>=2) dgm.fill("h16",ggp);
   
     	if (npart>=2) {
     		double dist=0, du=0;
 			int[]     npc = new int[50];       int[] neci = new int[50];       int[] neco = new int[50]; 
 			int[]     spc = new int[50];       int[]  sci = new int[50];       int[]  sco = new int[50]; 
 			double[]  epc = new double[50]; double[]  eci = new double[50]; double[]  eco = new double[50];  					
 			double[]  bpc = new double[50]; double[]  bci = new double[50]; double[]  bco = new double[50];
 					
			Vector3D[] r1 = new Vector3D[50]; Vector3[] c1 = new Vector3[50];
			Vector3D[] r4 = new Vector3D[50]; Vector3[] c4 = new Vector3[50]; 
			Vector3D[] r7 = new Vector3D[50]; Vector3[] c7 = new Vector3[50]; 
               	
 			double nopa = Math.toDegrees(Math.acos(plist.get(0).cosTheta(plist.get(1))));
            double dopa = opa-nopa;  //GEN-REC               
 			    	
 			int npp = 0, ipp = 0;
 			for (DetectorParticle dp : par) {
 			    if(dp.getSector(DetectorType.ECAL)!=trsec && dp.getPid()==MPpid) {
 			    	if(npp==0) p1 = dp; //Photon 1
 			    	if(npp==1) p2 = dp; //Photon 2
 			        for(int iresp = 0; iresp < cal.size(); iresp++){
 			        	CalorimeterResponse dr = (CalorimeterResponse)cal.get(iresp); 			              
 			    		int lay = dr.getDescriptor().getLayer(); 			    		  
 			    		if (dr.getAssociation(0)==ipp && dr.getDescriptor().getType()==DetectorType.ECAL) {
// 			            for (DetectorResponse dr : dp.getDetectorResponses()) { 			            	
// 			            	int lay = dr.getDescriptor().getType()==DetectorType.ECAL ? dr.getDescriptor().getLayer():0; 	
			    			
 			    			double dre=dr.getEnergy(),drt = dr.getTime()-stt, drb=dr.getPath()/drt/29.97; int drs=dr.getStatus(); 
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
		    	
 			if(npc[0]>=2)                            {r1[2].sub(r1[1]); dist=r1[2].mag();} 
 			if(npc[1]==1 && npc[2]==0 && neci[2]>=1) {r4[2].sub(r1[1]); dist=r4[2].mag();} 
 			if(npc[2]==1 && npc[1]==0 && neci[1]>=1) {r4[1].sub(r1[2]); dist=r4[1].mag();} 
 			    	
 			float dpcb = (float) (bpc[1]-bpc[2]); float dcib = (float) (bci[1]-bci[2]); float dcob = (float) (bco[1]-bco[2]); 
 			    	 			    	
 			if (correct && p1.countResponses(DetectorType.ECAL,1)==1 && p2.countResponses(DetectorType.ECAL,1)==1) {
 			    //Merged cluster associated with photon 2 	
 	 			if (p1.countResponses(DetectorType.ECAL,4)==0 && p2.countResponses(DetectorType.ECAL,4)==1) {
 	 			    double rat12 = 0.5*epc[1]/epc[2]; double rat21 = 0.5*epc[2]/epc[1]; rat12=0.5;rat21=0.5;
 	 			    plist.get(0).setP((epc[1]+eci[2]*rat21+eco[1])/ebmce.getSF(p2));
 	 			    plist.get(1).setP((epc[2]+eci[2]*rat12+eco[2])/ebmce.getSF(p2)); 	 			    		
 	 			}
 	 			//Merged cluster associated with photon 1 	
 	 	        if (p1.countResponses(DetectorType.ECAL,4)==1 && p2.countResponses(DetectorType.ECAL,4)==0) {
	 			    double rat12 = 0.5*epc[1]/epc[2]; double rat21 = 0.5*epc[2]/epc[1]; rat12=0.5; rat21=0.5;
 	 	 			plist.get(0).setP((epc[1]+eci[1]*rat21+eco[1])/ebmce.getSF(p1));
 	 	 			plist.get(1).setP((epc[2]+eci[1]*rat12+eco[2])/ebmce.getSF(p1));	 	 			    	
 	 	 	    }	 	 			    
 			} 	
 			    	
 	 	    REC   = getkin(plist);
 	 	    RECPZ = ebmce.getPizeroKinematics(plist);
 	 			    
 	 		dgm.fill("p12",GENPZ.get(0),RECPZ.get(4)); dgm.fill("p13",RECPZ.get(4),RECPZ.get(3));
 	 		    	
 	 		dgm.fill("pz0",RECPZ.get(3),RECPZ.get(5)); 
 	 		dgm.fill("pz1",GENPZ.get(0),RECPZ.get(2));
 	 		dgm.fill("pz2",GENPZ.get(0),GENPZ.get(0)-RECPZ.get(0));
 	 		dgm.fill("pz3",GENPZ.get(0),GENPZ.get(1)-RECPZ.get(1));
 	 		dgm.fill("pz4",RECPZ.get(3),GENPZ.get(0)-RECPZ.get(0));
 	 		dgm.fill("pz5",RECPZ.get(3),GENPZ.get(2)-RECPZ.get(2));
 	 			    
     	    dgm.fill("h11",opa);
     	       	    
     	    if (npc[0]>=2 && opa>2.0) {
     			Vector3 v1 = new Vector3(c1[1]); Vector3 v2 = new Vector3(c1[2]); v1.sub(v2);
     			dgm.fill("ef118",Math.abs(v1.x())>0?v1.x():-10+Math.random()*20);  //dU
     			dgm.fill("ef119",Math.abs(v1.y())>0?v1.y():-10+Math.random()*20);  //dV
     			dgm.fill("ef120",Math.abs(v1.z())>0?v1.z():-10+Math.random()*20);  //dW  					     					
     		}
     			 
     	    //Truth matching based on fixed angle of gamma 1
      		Boolean swap = Math.abs(GEN.get(1)- REC.get(1))<0.17 ? false:true;

     		double  delE1 = REC.get(swap?3:0)/GEN.get(0);
     		double delTH1 = REC.get(swap?4:1)-GEN.get(1);
     		double delPH1 = REC.get(swap?5:2)-GEN.get(2);
     		double  delE2 = REC.get(swap?0:3)/GEN.get(3);
     		double delTH2 = REC.get(swap?1:4)-GEN.get(4);
     		double delPH2 = REC.get(swap?2:5)-GEN.get(5);
     				
 	 		double delesum =          (REC.get(0)+REC.get(3))/(GEN.get(0)+GEN.get(3));    	
 	 	    double delepro = Math.sqrt(REC.get(0)*REC.get(3))/Math.sqrt(GEN.get(0)*GEN.get(3)); 	 			    
 	 			    
 			float e1 = (float) (p1.getEnergy(DetectorType.ECAL)/ebmce.getSF(p1)), b1 = (float) p1.getBeta();
 			float e2 = (float) (p2.getEnergy(DetectorType.ECAL)/ebmce.getSF(p2)), b2 = (float) p2.getBeta();			    	 

// 			    	if(opa>2.5) {
// 			    	System.out.println(REC.get(swap?3:0)+" "+e1); 
// 			    	System.out.println(REC.get(swap?0:3)+" "+e2);
// 			    	}
 			    	
 	 			    //Debug 
//     				if(debug && neci[0]==1 && opa>1.8 && opa<2.5) {
     		if(dbgAnalyzer && npc[0]>1 && neci[0]>1 && opa>2.5) {
     				System.out.println(" "); int scaf = 2;
     				System.out.println(getEventNumber());
     				     				
			        System.out.println(p1.getEnergy(DetectorType.ECAL)+" "+epc[1]+" "+eci[1]+" "+eco[1]);
			        System.out.println(p2.getEnergy(DetectorType.ECAL)+" "+epc[2]+" "+eci[2]+" "+eco[2]);
     				
     				System.out.println(p1.getEnergy(DetectorType.ECAL)+" "+(epc[1]+eci[1]/scaf+eco[1])/ebmce.getSF(p1));
     				System.out.println(p2.getEnergy(DetectorType.ECAL)+" "+(epc[2]+eci[1]/scaf+eco[2])/ebmce.getSF(p2));
     				
			    	System.out.println(npart+" "+npc[0]+" "+neci[0]+" "+neco[0]);
 	 			    System.out.println("Photon 1: "+p1.getMass()+" "+e1+" "+npc[1]+" "+neci[1]+" "+b1+" "+bpc[1]); 
 	 			    System.out.println("Photon 2: "+p2.getMass()+" "+e2+" "+npc[2]+" "+neci[2]+" "+b2+" "+bpc[2]); 
 	 			    
     				System.out.println("GEN,REC EN1,EN2 "+GEN.get(0)+" "+ REC.get(swap?3:0)+" "+GEN.get(3)+" "+ REC.get(swap?0:3));
     				System.out.println("GEN,REC TH1,TH2 "+GEN.get(1)+" "+ REC.get(swap?4:1)+" "+GEN.get(4)+" "+ REC.get(swap?1:4));
     				System.out.println("GEN,REC PH1,PH2 "+GEN.get(2)+" "+ REC.get(swap?5:2)+" "+GEN.get(5)+" "+ REC.get(swap?2:5));
 			    	System.out.println("GEN,REC opa "+opa+" "+nopa + " "+dist);
 		
    				System.out.println("Phot 1 "+swap+" "+delE1+" "+delTH1+" "+delPH1);
     				System.out.println("Phot 2 "+swap+" "+delE2+" "+delTH2+" "+delPH2);   
     		}
 			    	 			    	
 			if (opa>1.95) {
     			dgm.fill("e00", npc[0],npart);        dgm.fill("e01", neci[0],npart);        dgm.fill("e02", neco[0],npart);
     			dgm.fill("e10", npc[0],npart,delesum);dgm.fill("e11", neci[0],npart,delesum);dgm.fill("e12", neco[0],npart,delesum);
 			}     				

     		if (opa>2.40) {
     			dgm.fill("e03", npc[0],npart);        dgm.fill("e04", neci[0],npart);        dgm.fill("e05", neco[0],npart);
     			dgm.fill("e13", npc[0],npart,delesum);dgm.fill("e14", neci[0],npart,delesum);dgm.fill("e15", neco[0],npart,delesum);
     		}
     				
     		int stat = getStat(spc[1],sci[1],sco[1]);
     				
     		boolean opat = opa>2.5;
     				
     	    if (npc[0]>=2 && neci[0]==0) {
     			    	dgm.fill("c021",opa,npart);
     			    	dgm.fill("d00",GEN.get(4),delE1);
     			    	dgm.fill("d02",dist,      delE1);        			
     			    	dgm.fill("d03",opa,       delE1);        			
     			    	dgm.fill("d10",GEN.get(4),delE2);
     			    	dgm.fill("d12",dist,      delE2);        			
     			    	dgm.fill("d13",opa,       delE2); 
     			    	if(opat) dgm.fill("d14",stat,    delE2); 
     			    	
     			    	dgm.fill("dt00",GEN.get(4),delTH1);
     			    	dgm.fill("dt02",dist,      delTH1);        			
     			    	dgm.fill("dt03",opa,       delTH1);        			
     			    	dgm.fill("dt10",GEN.get(4),delTH2);
     			    	dgm.fill("dt12",dist,      delTH2);        			
     			    	dgm.fill("dt13",opa,       delTH2); 
     			    	
     			    	dgm.fill("df00",GEN.get(4),delPH1);
     			    	dgm.fill("df02",dist,      delPH1);        			
     			    	dgm.fill("df03",opa,       delPH1);        			
     			    	dgm.fill("df10",GEN.get(4),delPH2);
     			    	dgm.fill("df12",dist,      delPH2);        			
     			    	dgm.fill("df13",opa,       delPH2); 
     	    }
     				
     		if (npc[0]>=2 && neci[0]==1) {
     			    	dgm.fill("c022",opa,npart);
     					dgm.fill("d20",GEN.get(4),delE1);
     					dgm.fill("d22",dist,      delE1);        			
     					dgm.fill("d23",opa,       delE1);        			
     					dgm.fill("d30",GEN.get(4),delE2);
     					dgm.fill("d32",dist,      delE2);        			     				
     					dgm.fill("d33",opa,       delE2);  
     					if(opat) dgm.fill("d34",stat,    delE2);  
     					
     			    	dgm.fill("dt20",GEN.get(4),delTH1);
     			    	dgm.fill("dt22",dist,      delTH1);        			
     			    	dgm.fill("dt23",opa,       delTH1);        			
     			    	dgm.fill("dt30",GEN.get(4),delTH2);
     			    	dgm.fill("dt32",dist,      delTH2);        			
     			    	dgm.fill("dt33",opa,       delTH2); 
     			    	
     			    	dgm.fill("df20",GEN.get(4),delPH1);
     			    	dgm.fill("df22",dist,      delPH1);        			
     			    	dgm.fill("df23",opa,       delPH1);        			
     			    	dgm.fill("df30",GEN.get(4),delPH2);
     			    	dgm.fill("df32",dist,      delPH2);        			
     			    	dgm.fill("df33",opa,       delPH2); 
     			    	
     					dgm.fill("d60",GEN.get(4),delesum);
     					dgm.fill("d62",dist,      delesum);        			
     					dgm.fill("d63",opa,       delesum);        			
     					dgm.fill("d70",GEN.get(4),delepro);
     					dgm.fill("d72",dist,      delepro);        			     				
     					dgm.fill("d73",opa,       delepro);   
     					if(opat) dgm.fill("d74",stat,    delepro);   
     					dgm.fill("h12",opa);
     	   }
     				
     	   dgm.fill("c003",npart,delesum);
     	   if (opa>1.95) dgm.fill("c004",npart,delesum);
     				
     	   if (npc[0]>=2 && neci[0]>=2) {
 			    	    if(opat) {dgm.fill("p00",b1,npart);dgm.fill("p01",b2,npart);dgm.fill("p02",e1,npart);dgm.fill("p03",e2,npart);}
     			    	dgm.fill("c023",opa,npart);
     					dgm.fill("d40",GEN.get(4),delE1);
     					dgm.fill("d42",dist,      delE1);        			
     					dgm.fill("d43",opa,       delE1);        			
     					dgm.fill("d50",GEN.get(4),delE2);
     					dgm.fill("d52",dist,      delE2);        			     				
     					dgm.fill("d53",opa,       delE2);
     					dgm.fill("p07",ggp,       delE2);
     					if(opat) dgm.fill("d54",stat,    delE2);
     					
     					dgm.fill("d80",GEN.get(4),delesum);
     					dgm.fill("d82",dist,      delesum);        			
     					dgm.fill("d83",opa,       delesum);        			
     					dgm.fill("d90",GEN.get(4),delepro);
     					dgm.fill("d92",dist,      delepro);        			     				
     					dgm.fill("d93",opa,       delepro);   
     					if(opat)dgm.fill("d94",stat,    delepro);   
     					
     			    	dgm.fill("dt40",GEN.get(4),delTH1);
     			    	dgm.fill("dt42",dist,      delTH1);        			
     			    	dgm.fill("dt43",opa,       delTH1);        			
     			    	dgm.fill("dt50",GEN.get(4),delTH2);
     			    	dgm.fill("dt52",dist,      delTH2);        			
     			    	dgm.fill("dt53",opa,       delTH2); 
     			    	
     			    	dgm.fill("df40",GEN.get(4),delPH1);
     			    	dgm.fill("df42",dist,      delPH1);        			
     			    	dgm.fill("df43",opa,       delPH1);        			
     			    	dgm.fill("df50",GEN.get(4),delPH2);
     			    	dgm.fill("df52",dist,      delPH2);        			
     			    	dgm.fill("df53",opa,       delPH2); 
     			    	dgm.fill("h13", opa); 
     			    	
     			    	if(Math.abs(dpcb)<0.1) dgm.fill("h13dpcb", opa);
     					if(opa>=2) dgm.fill("h18", ggp);
     					if(Math.abs(dpcb)<0.1 && opa>=2) dgm.fill("h17",ggp);     					
     					
     					dgm.fill("c200", epc[1],bpc[1]);
     					dgm.fill("c201", epc[2],bpc[2]);
     					dgm.fill("c202", epc[1]/epc[2],dpcb);
     					dgm.fill("c203", opa, dpcb);
     					dgm.fill("c204", ggp, dpcb);
     					dgm.fill("c204a",spc[1], dpcb);
     					dgm.fill("c204b",ggp,spc[1]);    					
     					
     					dgm.fill("c205", eci[1],bci[1]);
     					dgm.fill("c206", eci[2],bci[2]);
     					dgm.fill("c207", eci[1]/eci[2],dcib);
     					dgm.fill("c208", opa, dcib);
     					dgm.fill("c209", ggp, dcib);
     					dgm.fill("c209a",sci[1], dcib);
     					dgm.fill("c209b",ggp,sci[1]);
    					
     					if(neco[0]>=2) {
     					dgm.fill("c210", eco[1],bco[1]);
     					dgm.fill("c211", eco[2],bco[2]);
     					dgm.fill("c212", eco[1]/eco[2],dcob);
     					dgm.fill("c213", opa, dcob);
     					dgm.fill("c214", ggp, dcob);
     					dgm.fill("c214a",sco[1], dcob);
     					dgm.fill("c214b",ggp,sco[1]);
     					}     				
     					
     		}
     				
     		if (npc[0]>=2 && neci[0]>0) {
     			    	dgm.fill("dt60",GEN.get(4),dopa);
     					dgm.fill("dt62",dist,      dopa);        			
     					dgm.fill("dt63",opa,       dopa); 
     					dgm.fill("h111",opa);
     		}
     				
     		if (npc[0]>0  && neci[0]>=2) dgm.fill("h1a", opa);
     		if (npc[0]>=2 && neci[0]==2) dgm.fill("h14", opa);
     		
     		if (npc[0]==2 && neci[0]==2) {
     					dgm.fill("h15", opa);
     					dgm.fill("d100",GEN.get(4),delesum);
     					dgm.fill("d102",dist,      delesum);        			
     					dgm.fill("d103",opa,       delesum);        			
     					dgm.fill("d110",GEN.get(4),delepro);
     					dgm.fill("d112",dist,      delepro);        			     				
     					dgm.fill("d113",opa,       delepro);   
     					dgm.fill("d114",spc[1],    delepro);   
     					dgm.fill("d120",GEN.get(4),delE1);
     					dgm.fill("d122",dist,      delE1);        			
     					dgm.fill("d123",opa,       delE1);        			
     					dgm.fill("d130",GEN.get(4),delE2);
     					dgm.fill("d132",dist,      delE2);        			     				
    					dgm.fill("d133",opa,       delE2);  
    					dgm.fill("d134",spc[1],    delE2);  
     					dgm.fill("p06",ggp,        delE2);    					
     					
     					for (int ic=0; ic<7; ic++) {
     						if(spc[1]==ic) dgm.fill("p20"+ic, ggp, delE2);
     						if(sci[1]==ic) dgm.fill("p21"+ic, ggp, delE2);
     						if(sco[1]==ic) dgm.fill("p22"+ic, ggp, delE2);
     						if(spc[1]==ic) dgm.fill("p30"+ic, ggp, dpcb);
     						if(sci[1]==ic) dgm.fill("p31"+ic, ggp, dcib);
     						if(sco[1]==ic) dgm.fill("p32"+ic, ggp, dcob);
    					}     				
    					     					
     					dgm.fill("p04", epc[1]/epc[2],eci[1]/eci[2]);
     					dgm.fill("p05", epc[2]/epc[1],eci[1]/eci[2]);
     					
     					if(opa>=2) dgm.fill("h19", ggp);
     					
    					dgm.fill("c100", epc[1],bpc[1]);
     					dgm.fill("c101", epc[2],bpc[2]);
     					dgm.fill("c102", epc[1]/epc[2],dpcb);
     					dgm.fill("c103", opa, dpcb);
     					dgm.fill("c104", ggp, dpcb);
     					dgm.fill("c104a",spc[1], dpcb);
     					dgm.fill("c104b",ggp,spc[1]);
     					
     					dgm.fill("c105", eci[1],bci[1]);
     					dgm.fill("c106", eci[2],bci[2]);
     					dgm.fill("c107", eci[1]/eci[2],dcib);
     					dgm.fill("c108", opa, dcib);
     					dgm.fill("c109", ggp, dcib);
     					dgm.fill("c109a",sci[1], dcib);
     					dgm.fill("c109b",ggp,sci[1]);
    					
     					if(neco[0]==2) {
     					dgm.fill("c110", eco[1],bco[1]);
     					dgm.fill("c111", eco[2],bco[2]);
     					dgm.fill("c112", eco[1]/eco[2],dcob);
     					dgm.fill("c113", opa, dcob);
     					dgm.fill("c114", ggp, dcob);
     					dgm.fill("c114a",sco[1], dcob);
     					dgm.fill("c114b",ggp,sco[1]);
     					}

     		}
   				
     		if (npc[0]==1 && neci[0]>1) {
     			    	dgm.fill("c024",opa,npart);
     					dgm.fill("dt70",GEN.get(4),dopa);
     					dgm.fill("dt72",dist,      dopa);        			
     					dgm.fill("dt73",opa,       dopa);  
     		}
     	}    			     			
   
    } 
   
    public void geteff() {
		dgm.geteff("eff1",  "h11", "h10");
		dgm.geteff("eff1a", "h1a", "h10");
		dgm.geteff("eff2",  "h12", "h10");
		dgm.geteff("eff3",  "h13", "h10"); dgm.geteff("eff3dpcb",  "h13dpcb", "h10");
		dgm.geteff("eff4",  "h14", "h10");
		dgm.geteff("eff5",  "h15", "h10");
		dgm.geteff("eff111","h111","h10");
		dgm.geteff("ef10",  "e10", "e00");
		dgm.geteff("ef11",  "e11", "e01");
		dgm.geteff("ef12",  "e12", "e02");
		dgm.geteff("ef13",  "e13", "e03");
		dgm.geteff("ef14",  "e14", "e04");
		dgm.geteff("ef15",  "e15", "e05");
		dgm.geteff("ef17",  "h17", "h16");
		dgm.geteff("ef18",  "h18", "h16");
		dgm.geteff("ef19",  "h19", "h16");		
    }
  
    @Override
    public void timerUpdate() {  

    }


}
