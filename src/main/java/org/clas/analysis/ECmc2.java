package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.analysis.ECPart.SFFunction;
import org.clas.tools.DataGroupManager;
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
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.eb.SamplingFractions;

public class ECmc2 extends DetectorMonitor {
	
	Event ev = new Event();
	EBMC  eb = new EBMC();
    List<DetectorParticle> np = new ArrayList<DetectorParticle>();
   
    HipoDataSync  writer = null;
	DataGroup dg = null;
	String tit = null;
	List<Particle> phot = new ArrayList<Particle>();
	List<Particle> neut = new ArrayList<Particle>();
//	List<Float> GEN =  new ArrayList<Float>();
//	List<Float> REC1 = new ArrayList<Float>();
//	List<Float> REC2 = new ArrayList<Float>();

    public ECmc2(String name) {
        super(name);
        
        dgm.setDetectorTabNames("MAIN","GENREC","EFFICIENCY");

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
        eb.setThresholds("Pizero",engine);
        eb.setGeom("2.5");
        eb.setGoodPhotons(12);
        eb.isMC = true;
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
    	createGENREC(0);
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
            dg.addDataSet(makeH2("h2a", 50,0,70, 50,0.95,4,"","Distance (cm)","E1 / E2"),n++); 
            dg.addDataSet(makeH2("h2b", 50,0,70, 50,-0.4,0.4,"","Distance (cm)","#Delta E / E"),n++); 
            dg.addDataSet(makeH2("h2f", 50,0,5.2,50,-0.4,0.4,"","#theta 12","#Delta E / E"),n++); 
            dg.addDataSet(makeH2("h2e", 80,0,70, 80,-1.5,1.5,"","Distance (cm)","#Delta#theta 12"),n++); 
            dg.addDataSet(makeH2("h2c", 50,0,4,  50,0,4,"nec=2","E #gamma 1 / E #gamma 2 (PCAL)","E #gamma 1 / E #gamma 2 (ECIN)"),n++);
            dg.addDataSet(makeH2("h2cc",50,0,4,  50,0,4,"nec=2 E#gamma1/E#gamma2>1.3","E #gamma 1 / E #gamma 2 (PCAL)","E #gamma 1 / E #gamma 2 (ECIN)"),n++);
            dg.addDataSet(makeH2("hnpp",50,0,5.2,10,1,11,"","Opening Angle (deg)","Number of photons"),n++);
            dg.addDataSet(makeH2("h2g", 10,1,11, 50,-1,1,"","Number of photons","#Delta E / E"),n++);
            dg.addDataSet(makeH1("h1gx",50,-1,1,"","#Delta E / E",1,4),n++);
            dg.addDataSet(makeH1("h1gy",10,1,11,"","Number of photons",1,4),n++);
            dg.addDataSet(makeH2("h2d",50,0,4,50,0,4,"","E #gamma 1 / E #gamma 2 (PCAL)","E #gamma 1 / E #gamma 2 (ECIN)"),n++);
    	}
    	
        this.getDataGroup().add(dg,0,st,k,run);         
    }
    

    public void createGENREC(int st) {
	  
    	switch (st) {        
    	case 0:
    		dgm.add("GENREC",6,4,0,st,getRunNumber());
        	                                      tit = "pc=2 ec=0";
            dgm.makeH2("d00",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma1");
            dgm.makeH2("d10",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma2");
            dgm.makeH2("d02",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma1"); 
            dgm.makeH2("d12",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma2"); 
            dgm.makeH2("d03",50,0,5.2,  50,-1,1,0,tit,"#theta12 (deg))",         "#DeltaE/E #gamma1");  
            dgm.makeH2("d13",50,0,5.2,  50,-1,1,0,tit,"#theta12 (deg))",         "#DeltaE/E #gamma2");  
                                                  tit = "pc=2 ec=1";
            dgm.makeH2("d20",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma1"); 
            dgm.makeH2("d30",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma2");  
            dgm.makeH2("d22",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma1");   
            dgm.makeH2("d32",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma2");  
            dgm.makeH2("d23",50,0,5.2,  50,-1,1,0,tit,"#theta12 (deg))",         "#DeltaE/E #gamma1");   
            dgm.makeH2("d33",50,0,5.2,  50,-1,1,0,tit,"#theta12 (deg))",         "#DeltaE/E #gamma2");
                                                  tit = "pc=2 ec=2";
            dgm.makeH2("d40",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma1"); 
            dgm.makeH2("d50",50,10,20,  50,-1,1,0,tit,"GEN #theta #gamma2 (deg)","#DeltaE/E #gamma2"); 
            dgm.makeH2("d42",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma1"); 
            dgm.makeH2("d52",50,0,70,   50,-1,1,0,tit,"Distance (cm)",           "#DeltaE/E #gamma2");  
            dgm.makeH2("d43",50,0,5.2,  50,-1,1,0,tit,"#theta12 (deg))",         "#DeltaE/E #gamma1");  
            dgm.makeH2("d53",50,0,5.2,  50,-1,1,0,tit,"#theta12 (deg))",         "#DeltaE/E #gamma2");              
                                                  tit = "pc=2 ec=1"; 
            dgm.makeH2("d60",50,10,20,  50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","(E1+E2)/2E"); 
            dgm.makeH2("d70",50,10,20,  50, 0,2,1,tit,"GEN #theta #gamma2 (deg)","SQRT(E1*E2)/E");   
            dgm.makeH2("d62",50,0,70,   50, 0,2,1,tit,"Distance (cm)",           "(E1+E2)/2E");    
            dgm.makeH2("d72",50,0,70,   50, 0,2,1,tit,"Distance (cm)",           "SQRT(E1*E2)/E");  
            dgm.makeH2("d63",50,0,5.2,  50, 0,2,1,tit,"#theta12 (deg))",         "(E1+E2)/2E");  
            dgm.makeH2("d73",50,0,5.2,  50, 0,2,1,tit,"#theta12 (deg))",         "SQRT(E1*E2)/E");             
    	}        
   	   
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
//    	geteff();
    	writer.close();
    	isAnalyzeDone = true;
    }
    
    public List<Float> get2 (List<Particle> list) {
    	List<Float> out = new ArrayList<Float>();
    	int n=0;
    	for (Particle p : list) {
		    out.add(n++,(float) p.e()); 
		    out.add(n++,(float) Math.toDegrees(p.theta()));  
		    out.add(n++,(float) Math.toDegrees(p.phi()));
    	}		 
		out.add(n++,(float) Math.toDegrees(Math.acos(list.get(0).cosTheta(list.get(1)))));   
		return out;
    }
    
    @Override
    public void processEvent(DataEvent de) {
		
		List<Float> GEN =  new ArrayList<Float>();
		List<Float> REC1 = new ArrayList<Float>();	
		
		int npp=0,indx=-1; float pthresh=1f;
		int n=0,np2=0,npc=0,neci=0,neco=0,sec=2,nphot=0, nneut=0, npart=0, etot=8;
		double e1=0,e2=0,oparec=0,the2=0,dist=0,epc1=0,epc2=0,eeci1=0,eeci2=0,eeco1=0,eeco2=0;
		Vector3D r1=null,r2=null;
				        
//		System.out.println(getEventNumber());
		
		boolean goodev = eb.readMC(de) && eb.pmc.size()==2;
				
        if (goodev) {    
        	
            engine.processDataEvent(de);  
            
        	eb.readEC(de,"ECAL::clusters");
        	np.clear(); np = eb.getNeutralPart(); npart = np.size(); 
        	
        	eb.getRECBanks(de,eb.eb); writer.writeEvent(de);
        	
    		GEN = get2(eb.pmc); 
        	
        	ev.init(de);        	
        	ev.setEventNumber(10);
        	ev.setStartTimeCut(-2000);
        	ev.requireOneElectron(false);
       	    ev.setElecTriggerSector(10);
       	    
 //    		phot.clear(); REC1.clear();
     		
//        	dg10.getH1F("h10").fill(GEN.get(6)); 
//       	    dgm.fill("h10", GEN.get(6));

       	    if(ev.procEvent(de)) {
       	    	
       	    phot = ev.getPART(0.0, 22);  neut = ev.getPART(0.0,2112); 
       	    nphot = phot.size(); nneut = neut.size();
      	    
//       	    System.out.println(npart+" "+nphot+" "+nneut);
       	    
 			for (DetectorParticle dp : np) { 
 				double ep = dp.getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, dp, eb.ccdb);
// 				System.out.println(ep+" "+dp.countResponses(DetectorType.ECAL,1)+" "+
// 						                  dp.countResponses(DetectorType.ECAL,4)+" "+
// 						                  dp.countResponses(DetectorType.ECAL,7));
// 				if(ep>pthresh) {       				
// 					for (DetectorResponse dr : dp.getDetectorResponses()) {
// 			    		int lay = dr.getDescriptor().getLayer();
// 			    		if(lay==1) {r1=dr.getPosition();npc++;epc1=dr.getEnergy();}    					
// 			    		if(lay==4) {nec++;eeci1=dr.getEnergy();}     					
// 					}
// 				}
 			}      
     		
     		if (npart>0) {
     			
 			    if (npart>=2) {
                	
 			    	DetectorParticle p1 = np.get(0);  //Photon 1
 			    	DetectorParticle p2 = np.get(1);  //Photon 2  
           
 			    	Vector3 n1 = p1.vector(); n1.unit();
 			    	Vector3 n2 = p2.vector(); n2.unit();
            
 			    	e1=np.get(0).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p1, eb.ccdb);
 			    	e2=np.get(1).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p2, eb.ccdb);
            
 			    	List<Particle> list = new ArrayList<Particle>();
 			    	list.add(new Particle(22,n1.x()*e1,n1.y()*e1,n1.z()*e1)); 
 			    	list.add(new Particle(22,n2.x()*e2,n2.y()*e2,n2.z()*e2));
                            
 			    	for (DetectorResponse dr : p1.getDetectorResponses()) {
 			    		int lay = dr.getDescriptor().getLayer();
 			    		if(lay==1) {npc++  ;  epc1=dr.getEnergy() ; r1=dr.getPosition();}    					
 			    		if(lay==4) {neci++ ; eeci1=dr.getEnergy();}
 			    		if(lay==7) {neco++ ; eeco1=dr.getEnergy();}
 			    	}
 				
 			    	for (DetectorResponse dr : p2.getDetectorResponses()) {
 			    		int lay = dr.getDescriptor().getLayer();
 			    		if(lay==1) {npc++  ;  epc2=dr.getEnergy() ; r2=dr.getPosition();}    					
 			    		if(lay==4) {neci++ ; eeci2=dr.getEnergy();}    					
 			    		if(lay==7) {neco++ ; eeco2=dr.getEnergy();}
 			    	}
 			    	
 			    	if(npc==2) {r2.sub(r1); dist=r2.mag();}
// 	 			    System.out.println(npc+" "+nec); 
 			    	     				
     				REC1=get2(phot);
     				
//     				System.out.println(GEN.get(0)+" "+REC1.get(0)+" "+GEN.get(3)+" "+REC1.get(3));
//     				System.out.println(GEN.get(1)+" "+REC1.get(1)+" "+GEN.get(4)+" "+REC1.get(4));
//     				System.out.println(GEN.get(2)+" "+REC1.get(2)+" "+GEN.get(5)+" "+REC1.get(5));
       	    
//     				dg10.getH1F("h10").fill(GEN.get(6));
//     				dgm.fill("h10",GEN.get(6));
     			
     				double dth11 = GEN.get(1)-REC1.get(1);  
     				double dth14 = Math.abs(GEN.get(1)-REC1.get(4));
     				Boolean swap = Math.abs(dth11)<0.17 ? false:true;
     				
//     				System.out.println(dth11+" "+dth14+" "+swap);
     				
     				double  delE1 = swap ? GEN.get(0)-REC1.get(3):GEN.get(0)-REC1.get(0);
     				double delTH1 = swap ? GEN.get(1)-REC1.get(4):GEN.get(1)-REC1.get(1);
     				double delPH1 = swap ? GEN.get(2)-REC1.get(5):GEN.get(2)-REC1.get(2);
     				double  delE2 = swap ? GEN.get(3)-REC1.get(0):GEN.get(3)-REC1.get(3);
     				double delTH2 = swap ? GEN.get(4)-REC1.get(1):GEN.get(4)-REC1.get(4);
     				double delPH2 = swap ? GEN.get(5)-REC1.get(2):GEN.get(5)-REC1.get(5);
     				double dopa = GEN.get(6)-REC1.get(6);
     				
 //    				System.out.println(delE1+" "+delTH1+" "+delPH1);
 //    				System.out.println(delE2+" "+delTH2+" "+delPH2);   
     				
     			    if (npc==2 && neci==0) {
     			    	dgm.fill("d00",GEN.get(4),delE1/GEN.get(0));
     			    	dgm.fill("d02",dist,      delE1/GEN.get(0));        			
     			    	dgm.fill("d03",GEN.get(6),delE1/GEN.get(0));        			
     			    	dgm.fill("d10",GEN.get(4),delE2/GEN.get(3));
     			    	dgm.fill("d12",dist,      delE2/GEN.get(3));        			
     			    	dgm.fill("d13",GEN.get(6),delE2/GEN.get(3)); 
     			    }
     				
     				if (npc==2 && neci==1) {
     					dgm.fill("d20",GEN.get(4),delE1/GEN.get(0));
     					dgm.fill("d22",dist,      delE1/GEN.get(0));        			
     					dgm.fill("d23",GEN.get(6),delE1/GEN.get(0));        			
     					dgm.fill("d30",GEN.get(4),delE2/GEN.get(3));
     					dgm.fill("d32",dist,      delE2/GEN.get(3));        			     				
     					dgm.fill("d33",GEN.get(6),delE2/GEN.get(3));      
               			
     					dgm.fill("d60",GEN.get(4),(REC1.get(0)+REC1.get(3))/2/GEN.get(0));
     					dgm.fill("d62",dist,      (REC1.get(0)+REC1.get(3))/2/GEN.get(0));        			
     					dgm.fill("d63",GEN.get(6),(REC1.get(0)+REC1.get(3))/2/GEN.get(0));        			
     					dgm.fill("d70",GEN.get(4),Math.sqrt(REC1.get(0)*REC1.get(3))/GEN.get(3));
     					dgm.fill("d72",dist,      Math.sqrt(REC1.get(0)*REC1.get(3))/GEN.get(3));        			     				
     					dgm.fill("d73",GEN.get(6),Math.sqrt(REC1.get(0)*REC1.get(3))/GEN.get(3));               			
     			    }
     				
     				if (npc==2 && neci==2) {
     					dgm.fill("d40",GEN.get(4),delE1/GEN.get(0));
     					dgm.fill("d42",dist,      delE1/GEN.get(0));        			
     					dgm.fill("d43",GEN.get(6),delE1/GEN.get(0));        			
     					dgm.fill("d50",GEN.get(4),delE2/GEN.get(3));
     					dgm.fill("d52",dist,      delE2/GEN.get(3));        			     				
     					dgm.fill("d53",GEN.get(6),delE2/GEN.get(3));
     				}   			
     			}    			     			
     		}     		
       	    }   
        }
 
/*    			        	
     			npp=0;indx=-1;
        	
     			dg11.getH1F("h11").fill(GEN.get(6));
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

     			dg0.getH2F("hnpp").fill(GEN.get(6),npp);
        	
     			if (npp>1)  dg11.getH1F("h12").fill(GEN.get(6));
     			if (npp==2) dg11.getH1F("h13").fill(GEN.get(6));
     			if (npp>=2) dg11.getH1F("h14").fill(GEN.get(6));
     			if (npart==2) {
     				double e1=0,e2=0,oparec=0,the2=0,dist=0;
            	
     				DetectorParticle p1 = np.get(0);  //Photon 1
     				DetectorParticle p2 = np.get(1);  //Photon 2  
               
     				Vector3 n1 = p1.vector(); n1.unit();
     				Vector3 n2 = p2.vector(); n2.unit();
                
     				e1=np.get(0).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p1, eb.ccdb);
     				e2=np.get(1).getEnergy(DetectorType.ECAL)/SamplingFractions.getMean(22, p2, eb.ccdb);
                
     				List<Particle> list = new ArrayList<Particle>();
     				list.add(new Particle(22,n1.x()*e1,n1.y()*e1,n1.z()*e1)); 
     				list.add(new Particle(22,n2.x()*e2,n2.y()*e2,n2.z()*e2)); 
                             	
 //    				System.out.println(GEN.get(0)+" "+REC1.get(0)+" "+GEN.get(3)+" "+REC1.get(3));
 //    				System.out.println(GEN.get(1)+" "+REC1.get(1)+" "+GEN.get(4)+" "+REC1.get(4));
 //    				System.out.println(GEN.get(2)+" "+REC1.get(2)+" "+GEN.get(5)+" "+REC1.get(5));
                                
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
     				dg0.getH2F("h2f").fill(GEN.get(6), (e1+e2)/etot-1);
     				if(npc==2) {
     					r2.sub(r1); dist=r2.mag();
     					dg0.getH2F("h2a").fill(dist, e2>0?e1/e2:1000);
     					dg0.getH2F("h2b").fill(dist, (e1+e2)/etot-1);
                		dg0.getH2F("h2e").fill(dist,GEN.get(6)-REC1.get(6));
                		if(nec==2)              dg0.getH2F("h2c" ).fill(epc1/epc2,eeci1/eeci2);
                		if(nec==2 && e1/e2>1.3) dg0.getH2F("h2cc").fill(epc1/epc2,eeci1/eeci2);
                	
                		dg0.getH2F("h2g").fill(npp,(e1+e2)/etot-1);
            			dgr.getH2F("d00").fill(GEN.get(3),1-delE1/GEN.get(0));
            			dgr.getH2F("d01").fill(GEN.get(2),1-delE1/GEN.get(0));
            			dgr.getH2F("d02").fill(dist,      1-delE1/GEN.get(0));
            			dgr.getH2F("d03").fill(dopa,      1-delE1/GEN.get(0));
            			dgr.getH2F("d10").fill(GEN.get(3),1-delE2/GEN.get(0));
            			dgr.getH2F("d11").fill(GEN.get(2),1-delE2/GEN.get(0));
            			dgr.getH2F("d12").fill(dist,      1-delE2/GEN.get(0));
            			dgr.getH2F("d13").fill(dopa,      1-delE2/GEN.get(0));
            			
     				}
     				
     			}
     			*/
        
    }
   
    public void plotMCHistos() {      
//        plot("MAIN");
        plot("GENREC");
//        plot("EFFICIENCY");   
    }
    
    @Override
    public void plotEvent(DataEvent de) {
//        analyze();         
    }
    
    public void geteff() {
		int run = getRunNumber();
		DataGroup  dg0 = dgm.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("MAIN"),run);
		DataGroup dg10 = dgm.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		DataGroup dg11 = dgm.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("EFFICIENCY"),run);
		dg10.getH1F("eff1").add(H1F.divide(dg11.getH1F("h11"),dg10.getH1F("h10")));
		dg10.getH1F("eff2").add(H1F.divide(dg11.getH1F("h12"),dg10.getH1F("h10")));
		dg10.getH1F("eff3").add(H1F.divide(dg11.getH1F("h13"),dg10.getH1F("h10")));
		dg10.getH1F("eff4").add(H1F.divide(dg11.getH1F("h14"),dg10.getH1F("h10")));
		dg0.getH1F("h1gx").add(dg0.getH2F("h2g").projectionY());
		dg0.getH1F("h1gy").add(dg0.getH2F("h2g").projectionX());		
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,0, getRunNumber());
    } 
        
    @Override
    public void timerUpdate() {  

    }


}
