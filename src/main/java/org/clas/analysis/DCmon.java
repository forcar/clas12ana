package org.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;

import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.geom.prim.Point3D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

public class DCmon extends DetectorMonitor {
	
	boolean doDC = true;
	DataBank RecPart = null, RecTrk=null;
	HashMap<Integer,ArrayList<Integer>> part2trk  = null;
	
    public DCmon(String name) {
        super(name);
    	
        dgmActive=true; 
        this.setDetectorTabNames("DOCA",
        		                 "DTIM",
        		                 "WIRE");

        this.usePCCheckBox(true);
        this.useSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();       
    }
    
    public void localinit() {
        System.out.println(getDetectorName()+".localinit()"); 
        tl.setFitData(Fits);
    } 
    
    public void localclear() {
    	System.out.println(getDetectorName()+".localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	FitSummary.clear();
    	Fits.clear();
    	tl.fitData.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
    public void init(int run) {
	    setRunNumber(run);
	    runlist.add(run);     	
    }
    
    @Override  
    public void createHistos(int run) {
    	histosExist = true;
    	System.out.println(getDetectorName()+".createHistos("+run+")");
	    init(run);
	    createDOCA(0);
	    createDTIM(0);
	    createWIRE(0);

    }
    
    public void createDOCA(int st) {
    	
    	float max[] = {0.85f,0.85f,1.3f,1.3f,2f,2f};
    	for(int is=1; is<7; is++) {
    		for(int lr=1; lr<3; lr++){
            	dgm.add("DOCA",6,6,10*is+lr,st,getRunNumber()); 
   			 	for(int sl=1; sl<7; sl++)  {
   			 		for(int il=1; il<7; il++) {
   			 			String tit = "S"+is+" SL"+sl+" L"+il+" LR"+(lr==1?"+":"-");    	    	   	    	    	   	
   			 			dgm.makeH2("doca"+is+sl+il+lr,50,0,max[sl-1],50,-0.2,0.2,0,tit+" "," "," "); 
   			 		}
        		}
    		}
    	}
    }
    
    public void createDTIM(int st) {
    	
    	float max[] = {170f,170f,350f,350f,700f,700f};
    	for(int is=1; is<7; is++) {
    		for(int lr=1; lr<3; lr++){
            	dgm.add("DTIM",6,6,10*is+lr,st,getRunNumber()); 
   			 	for(int sl=1; sl<7; sl++)  {
   			 		for(int il=1; il<7; il++) {
   			 			String tit = "S"+is+" SL"+sl+" L"+il+" LR"+(lr==1?"+":"-");    	    	   	    	    	   	
   			 			dgm.makeH2("dtim"+is+sl+il+lr,50,-30,max[sl-1],50,-0.2,0.2,0,tit+" "," "," "); 
   			 		}
        		}
    		}
    	}
    }
    
    public void createWIRE(int st) {
    	
    	for(int is=1; is<7; is++) {
    		for(int lr=1; lr<3; lr++){
            	dgm.add("WIRE",6,6,10*is+lr,st,getRunNumber()); 
   			 	for(int sl=1; sl<7; sl++)  {
   			 		for(int il=1; il<7; il++) {
   			 			String tit = "S"+is+" SL"+sl+" L"+il+" LR"+(lr==1?"+":"-");    	    	   	    	    	   	
   			 			dgm.makeH2("wire"+is+sl+il+lr,112,0.5,112.5,50,-0.2,0.2,0,tit+" "," "," ");  
   			 		}
        		}
    		}
    	}
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);   	
    }
    
    public void plotSummary(int run) {
        setRunNumber(run);
        plot("DOCA"); plot("DTIM"); plot("WIRE");
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,10*getActiveSector()+getActivePC(),0,getRunNumber());
    }
    
    public boolean processFilter(DataEvent event) {
    	int sec = getElecTriggerSector(true);       
        RecPart  = event.hasBank("REC::Particle") ? event.getBank("REC::Particle"):null;       
        RecTrk   = event.hasBank("REC::Track")    ? event.getBank("REC::Track"):null;       
    	boolean test1 = RecPart!=null && RecPart.rows()!=0;
    	boolean test2 = RecTrk!=null  && RecTrk.rows()!=0;
    	boolean test3 = sec>0 && sec<7;
    	if(test2) part2trk = mapByIndex(RecTrk,"pindex");
    	return test1 && test2 && test3;
    }
    
    @Override    
    public void processEvent(DataEvent event) {  
    	if(!processFilter(event)) return;
    	if(doDC) processDC(event);
    }
    
    public Particle getPart(int ipart) {
    	int    pid  = RecPart.getInt("pid",ipart);
    	double epx  = RecPart.getFloat("px",ipart);
    	double epy  = RecPart.getFloat("py",ipart);
    	double epz  = RecPart.getFloat("pz",ipart);
    	return new Particle(pid>0?pid:-11,epx,epy,epz);    	    	
    }
 
    public void processDC(DataEvent event) {
    	
    	Particle neg=null,eneg=null;
    	
        ArrayList<Integer> eleCandi = new ArrayList<>();
        ArrayList<Integer>    eleEB = new ArrayList<>();
        ArrayList<Integer>    proEB = new ArrayList<>();

        for (int ipart=0; ipart<RecPart.rows(); ipart++) {        	
            int pid    = RecPart.getInt("pid",ipart);
            int charge = RecPart.getInt("charge",ipart);
            int status = RecPart.getInt("status",ipart);       
            final boolean isFD = (int)(Math.abs(status)/1000) == 2;

            if (isFD && charge < 0 && Math.abs(pid)>0) eleCandi.add(ipart);
            if (isFD && pid==11)                          eleEB.add(ipart);
            if (isFD && pid==2212)                        proEB.add(ipart);
        }
        
        if (!(proEB.size()==1)) return;

        int itr = -1, s=0;
        for (Integer ipart : proEB) {
        	neg = getPart(ipart);
    		s   = part2trk.containsKey(ipart) ? RecTrk.getInt("sector", part2trk.get(ipart).get(0)):0; 
    		itr = part2trk.containsKey(ipart) ? RecTrk.getInt("index",  part2trk.get(ipart).get(0)):0; 
        }
        

	   if (event.hasBank("TimeBasedTrkg::TBHits")==true){
            DataBank  bank = event.getBank("TimeBasedTrkg::TBHits");
            for(int i = 0; i < bank.rows(); i++){
                int   is = bank.getByte("sector",i);
                int   il = bank.getByte("layer",i);
                int   sl = bank.getByte("superlayer",i);
                int wire = bank.getInt("wire",i);               
                int  tdc = bank.getInt("TDC",i);
                int   lr = bank.getByte("LR",i);
                int trid = bank.getByte("trkID",i);
                
                float  time = bank.getFloat("time",i);               
                float  doca = bank.getFloat("doca",i); 
                float trdoc = bank.getFloat("trkDoca",i); 
                float  fitr = bank.getFloat("fitResidual",i);
                float  timr = bank.getFloat("timeResidual",i);
               
                if(itr==trid && tdc>0)  {
                	dgm.fill("doca"+is+sl+il+(lr==1?1:2),doca,timr);
                	dgm.fill("dtim"+is+sl+il+(lr==1?1:2),time,timr);
                	dgm.fill("wire"+is+sl+il+(lr==1?1:2),wire,timr);
                }
            }
       }
    }

}
