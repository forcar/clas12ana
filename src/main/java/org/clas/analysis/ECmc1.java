package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.analysis.ECPart.SFFunction;
import org.clas.tools.EBMC;
import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc1 extends DetectorMonitor {
	
	Event ev = new Event();
	EBMC  eb = new EBMC();
    List<DetectorParticle> np = new ArrayList<DetectorParticle>();
    HipoDataSync  writer = null;
	DataGroup dg = null;

    public ECmc1(String name) {
        super(name);
        setDetectorTabNames("MAIN","EFFICIENCY");

        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);

        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECmc2.localinit()");
        engine.init();
        engine.isMC = true;
        engine.setVariation("default");
        engine.setPCALTrackingPlane(9);
        engine.setCalRun(10);                
        eb.getCCDB(10);
        eb.setThresholds("Test",engine);
        eb.setGeom("2.5");
        eb.setGoodPhotons(12);
        tl.setFitData(Fits);
        writer = new HipoDataSync();
        writer.open("/Users/colesmith/CLAS12ANA/ECmc2/photon_dist.hipo");
    }
    
    public void localclear() {
    	System.out.println("ECmc2.localclear()");
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
	    System.out.println("ECmc2:createHistos("+run+")");
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;    	
    	createMAIN(0);
       	createEFFICIENCY(0);
       	createEFFICIENCY(1);
    }
    
    public void createMAIN(int st) {
    	
    	String tab = "MAIN", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0;
        tag = st+"_"+k+"_"+run;    
        
    	switch (st) {        
        case 0:  
            dg = new DataGroup(4,4);
            dg.addDataSet(makeH2("h2a",50,0,70,50,0.95,4,"","Distance (cm)","E1 / E2"),n++);
            dg.addDataSet(makeH2("h2b",50,0,70,50,-0.4,0.4,"","Distance (cm)","#Delta E / E"),n++);
            dg.addDataSet(makeH2("h2f",50,0,5.2,50,-0.4,0.4,"","#theta 12","#Delta E / E"),n++);
            dg.addDataSet(makeH2("h2e",80,0,70,80,-1.5,1.5,"","Distance (cm)","#Delta#theta 12"),n++);
            dg.addDataSet(makeH2("h2c",50,0,4,50,0,4,"nec=2","E #gamma 1 / E #gamma 2 (PCAL)","E #gamma 1 / E #gamma 2 (ECIN)"),n++);
            dg.addDataSet(makeH2("h2cc",50,0,4,50,0,4,"nec=2 E#gamma1/E#gamma2>1.3","E #gamma 1 / E #gamma 2 (PCAL)","E #gamma 1 / E #gamma 2 (ECIN)"),n++);
            dg.addDataSet(makeH2("hnpp",50,0,5.2,10,1,11,"","Opening Angle (deg)","Number of photons"),n++);
            dg.addDataSet(makeH2("h2g",10,1,11,50,-1,1,"","Number of photons","#Delta E / E"),n++);
            dg.addDataSet(makeH1("h1gx",50,-1,1,"","#Delta E / E",1,4),n++);
            dg.addDataSet(makeH1("h1gy",10,1,11,"","Number of photons",1,4),n++);
            dg.addDataSet(makeH2("h2d",50,0,4,50,0,4,"","E #gamma 1 / E #gamma 2 (PCAL)","E #gamma 1 / E #gamma 2 (ECIN)"),n++);
    	}
    	
        this.getDataGroup().add(dg,0,st,k,run);         
    }
    
    public void createEFFICIENCY(int st) {
    	
    	String tab = "EFFICIENCY", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab), n=0;
        tag = st+"_"+k+"_"+run;    
        
    	switch (st) {        
        case 0:          
            dg = new DataGroup(5,2);
            dg.addDataSet(makeH1("h10", 50,0,5.2,"GEN","Opening Angle (deg)","Efficiency",1,4),n++);
            dg.addDataSet(makeH1("eff1",50,0,5.2,"n>0","Opening Angle (deg)","Efficiency",1,4),n++);
            dg.addDataSet(makeH1("eff2",50,0,5.2,"n>1","Opening Angle (deg)","Efficiency",1,3),n++);
            dg.addDataSet(makeH1("eff3",50,0,5.2,"n==2","Opening Angle (deg)","Efficiency",1,2),n++);
            dg.addDataSet(makeH1("eff4",50,0,5.2,"n>=2","Opening Angle (deg)","Efficiency",1,5),n++);
            break;
        case 1:
            dg = new DataGroup(4,2);
            dg.addDataSet(makeH1("h11",50,0,5.2,"n>0","Opening Angle (deg)","Efficiency",1,4),n++);
            dg.addDataSet(makeH1("h12",50,0,5.2,"n>1","Opening Angle (deg)","Efficiency",1,3),n++);
            dg.addDataSet(makeH1("h13",50,0,5.2,"n==2","Opening Angle (deg)","Efficiency",1,2),n++);
            dg.addDataSet(makeH1("h14",50,0,5.2,"n>=2","Opening Angle (deg)","Efficiency",1,5),n++);
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
    	geteff();
    	isAnalyzeDone = true;
    }
    
    @Override
    public void processEvent(DataEvent de) {
    	
		int run = getRunNumber();
		DataGroup  dg0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("MAIN"),run);
		DataGroup dg10 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg11 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		
		boolean goodev = eb.readMC(de) && eb.pmc.size()==2;
		float      opa =(float) Math.toDegrees(Math.acos(eb.pmc.get(0).cosTheta(eb.pmc.get(1))));
		
		int npp=0,indx=-1; float pthresh=1f;
		int n=0,np2=0,npc=0,nec=0, sec=2, etot=8;
		
        if (goodev) {                   
            engine.processDataEvent(de);                      
        	eb.readEC(de,"ECAL::clusters");
        	np.clear(); np = eb.getNeutralPart(); 
        	eb.getRECBanks(de,eb.eb);
        	writer.writeEvent(de);
        	dg10.getH1F("h10").fill(opa);
        	npp=0;indx=-1;
        	if (np.size()>0) {
        		dg11.getH1F("h11").fill(opa);
        		for (DetectorParticle phot : np) {
        			indx++;
            		double ep = phot.getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, phot, eb.ccdb);
        			if(ep>pthresh) {
//        				System.out.println(npp+" "+indx+" "+ep);
        				npp++;           				
        				for (DetectorResponse dr : phot.getDetectorResponses()) {
        					
        				}
        			}
        		}
        	}
        	dg0.getH2F("hnpp").fill(opa,npp);
        	if (npp>1)  dg11.getH1F("h12").fill(opa);
        	if (npp==2) dg11.getH1F("h13").fill(opa);
        	if (npp>=2) dg11.getH1F("h14").fill(opa);
        	if (npp>=2) {
            	double e1=0,e2=0,the1=0,the2=0,dist=0;
                DetectorParticle p1 = np.get(0);  //Photon 1
                DetectorParticle p2 = np.get(1);  //Photon 2   
                Vector3 n1 = p1.vector(); n1.unit();
                Vector3 n2 = p2.vector(); n2.unit();
                e1=np.get(0).getEnergy(DetectorType.ECAL);
                e2=np.get(1).getEnergy(DetectorType.ECAL);
                eb.SF1db = SamplingFractions.getMean(22, np.get(0), eb.ccdb);
                e1 = e1/eb.SF1db;
                eb.SF1db = SamplingFractions.getMean(22, np.get(1), eb.ccdb);
                e2 = e2/eb.SF1db;
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
				dg0.getH2F("h2f").fill(opa, (e1+e2)/etot-1);
                if(npc==2) {
                r2.sub(r1); dist=r2.mag();
                dg0.getH2F("h2a").fill(dist, e2>0?e1/e2:1000);
                dg0.getH2F("h2b").fill(dist, (e1+e2)/etot-1);
                dg0.getH2F("h2e").fill(dist,opa-the1);
            	if(nec==2)              dg0.getH2F("h2c" ).fill(epc1/epc2,eeci1/eeci2);
            	if(nec==2 && e1/e2>1.3) dg0.getH2F("h2cc").fill(epc1/epc2,eeci1/eeci2);
            	dg0.getH2F("h2g").fill(npp,(e1+e2)/etot-1);
        	    }
        	}           		
        }
        
    }
   
    public void plotMCHistos() {      
        plot("MAIN");
        plot("EFFICIENCY");   
    }
    
    @Override
    public void plotEvent(DataEvent de) {
//        analyze();         
    }
    
    public void geteff() {
		int run = getRunNumber();
		DataGroup  dg0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("MAIN"),run);
		DataGroup dg10 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg11 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		dg10.getH1F("eff1").add(H1F.divide(dg11.getH1F("h11"),dg10.getH1F("h10")));
		dg10.getH1F("eff2").add(H1F.divide(dg11.getH1F("h12"),dg10.getH1F("h10")));
		dg10.getH1F("eff3").add(H1F.divide(dg11.getH1F("h13"),dg10.getH1F("h10")));
		dg10.getH1F("eff4").add(H1F.divide(dg11.getH1F("h14"),dg10.getH1F("h10")));
		dg0.getH1F("h1gx").add(dg0.getH2F("h2g").projectionY());
		dg0.getH1F("h1gy").add(dg0.getH2F("h2g").projectionX());		
    }
    
    public void plot(String tabname) {      	
    	int index = getDetectorTabNames().indexOf(tabname);
    	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));   	
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
