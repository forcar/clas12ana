package org.clas.analysis;

import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ECmon extends DetectorMonitor {
	
	IndexedList<List<Integer>> tdcs = new IndexedList<List<Integer>>(3);
    public static double[] SCALE  = {10,10,10,10,10,10,10,10,10}; // Fitter.ADC/SCALE is plotted and fitted in ECMon
    public static double[] SCALE5 = {10,10,10,5,5,5,5,5,5};       // Sector 5 ECAL uses EMI PMTs near max voltage
    IndexedTable time=null, ftime=null, gtw=null, fo=null, fgo=null, tmf=null, tmfcut=null, tgo=null, gain=null, veff=null, rfT=null;
	float FTOFFSET, TOFFSET, TMFCUT;
    static float tps = (float) 0.02345;
    
	
    public ECmon(String name) {
        super(name);
    	
        dgmActive=true; 
        this.setDetectorTabNames("AT",
        		                 "EFF");
        this.useCALUVWSECButtons(true);
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
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
    public void init(int run) {
	    setRunNumber(run);
	    runlist.add(run);     	
    }
    
    @Override  
    public void createHistos(int run) {
	    System.out.println(getDetectorName()+".createHistos("+run+")");
	    init(run);
	    createAT(0);
	    createEFF(0);
    	histosExist = true;
    }
    
    public void createAT(int st) {
    	
        int il;
    	String[] det = {"PCAL","ECIN","ECOU"}; 
    	String[] uvw = {"U","V","W"};
    	
    	for(int is=1; is<7; is++) {
    	for(int id=0; id<3; id++) {
    	int np = id==0?68:36; dgm.add("AT",3,4,10*is+id,st,getRunNumber());
    	String tit = "Sector "+is+" "+det[id]+" ";
        for(il=0; il<3; il++) dgm.makeH2("at1"+is+id+il,100, 0,150,100,0,240,-1,tit+uvw[il]+" STRIPS","FADC","TDC");
        for(il=0; il<3; il++) dgm.makeH2("at2"+is+id+il,100, 0,150,100,0,240,-1," ","FADC","TDC");
        for(il=0; il<3; il++) dgm.makeH2("at3"+is+id+il,100, 0,150,100,-11,11,-1," ","FADC","TDC-FADC");
        for(il=0; il<3; il++) dgm.makeH2("at4"+is+id+il,100, 0,150,np,1,np+1,-1," ","FADC","PMT");
    	}
    	}   	
    }
    
    public void createEFF(int st) {
        int il;
        String tit2 = "ABS(FADC-TDC)<10  55<TDC<120";
    	String[] det = {"PCAL","ECIN","ECOU"}; 
    	String[] uvw = {"U","V","W"};
    	for(int is=1; is<7; is++) {
    	for(int id=0; id<3; id++) {    		
        	int np = id==0?68:36; dgm.add("EFF",3,2,10*is+id,st,getRunNumber());
        	String tit1 = "Sector "+is+" "+det[id]+" ";
            for(il=0; il<3; il++) {dgm.makeGE("ef1"+is+id+il,0,tit1+uvw[il],"FADC","TDC EFFICIENCY");
                                   dgm.cc("ef1"+is+id+il,false,false,0,1.05f,0,0);}
            for(il=0; il<3; il++) {dgm.makeGE("ef2"+is+id+il,0,tit2,"FADC","TDC EFFICIENCY");
                                   dgm.cc("ef2"+is+id+il,false,false,0,1.05f,0,0);}           
    	}
    	}   	
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);
       	plotAnalysis(run);   	
    }
    
    public void plotSummary(int run) {
        setRunNumber(run);
        plot("AT");
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;  
    	plot("EFF");
    }
    
    @Override         
    public void plotEvent(DataEvent event) {
    	analyze();
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");   
    	getEff();
    	isAnalyzeDone = true;
    }
    
    @Override   
    public void initCCDB(int runno) {
    	System.out.println(getDetectorName()+".initCCDB("+runno+") "); 
    	fo      = cm.getConstants(runno, "/calibration/ec/fadc_offset");        //Crate/fiber FADC offsets 
        fgo     = cm.getConstants(runno, "/calibration/ec/fadc_global_offset"); //FADC capture window 
        gtw     = cm.getConstants(runno, "/calibration/ec/global_time_walk");   //Global time walk correction using raw ADC
        tmfcut  = cm.getConstants(runno, "/calibration/ec/tmf_window");
        tmf     = cm.getConstants(runno, "/calibration/ec/tmf_offset");         //TDC-FADC offsets
        time    = cm.getConstants(runno, "/calibration/ec/timing");
        ftime   = cm.getConstants(runno, "/calibration/ec/ftiming"); 
        
        FTOFFSET = (float) fgo.getDoubleValue("global_offset",0,0,0);
        TMFCUT   = (float) tmfcut.getDoubleValue("window", 0,0,0); //acceptance window for TDC-FADC cut
    }
    
    public void getEff() {
       
       for (int is=1; is<7; is++) {
       for (int id=0; id<3; id++) {
       for (int il=0; il<3; il++) {
    	   H2F h2at1 = dgm.getH2F("at1"+is+id+il);
    	   H2F h2at2 = dgm.getH2F("at2"+is+id+il);
    	   H2F h2at4 = dgm.getH2F("at4"+is+id+il);
    	   H1F h1at1 = projectionX(h2at1,1,240);
    	   H1F h1at2 = projectionX(h2at2,1,240);
    	   H1F h1at4 = h2at4.projectionX();
    	   dgm.getGE("ef1"+is+id+il).copy(H1F.divide(h1at1, h1at4).getGraph());
    	   dgm.getGE("ef1"+is+id+il).setMarkerColor(1); 
    	   dgm.getGE("ef1"+is+id+il).setLineColor(2);
    	   dgm.getGE("ef2"+is+id+il).copy(H1F.divide(h1at2, h1at4).getGraph());
    	   dgm.getGE("ef2"+is+id+il).setMarkerColor(1); 
    	   dgm.getGE("ef2"+is+id+il).setLineColor(2);
       }
       }
       }
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,10*getActiveSector()+getActiveLayer(),0,getRunNumber());
    }
    
    public int getDet(int layer) {
        int[] il = {0,0,0,1,1,1,2,2,2}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
     }
    
    public int getLay(int layer) {
        int[] il = {1,2,3,1,2,3,1,2,3}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
    }   
    
    @Override     
    public void processEvent(DataEvent event) {

        tdcs.clear();
        
        if(event.hasBank("ECAL::tdc")==true){
            DataBank  bank = event.getBank("ECAL::tdc");
            for(int i = 0; i < bank.rows(); i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  ip = bank.getShort("component",i);               
                int tdc = bank.getInt("TDC",i);
                if(tdc>0) {
                    if(!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Integer>(),is,il,ip);
                        tdcs.getItem(is,il,ip).add(tdc);       
                }
            }
        }
        
        if(event.hasBank("ECAL::adc")==true){
            DataBank  bank = event.getBank("ECAL::adc");
            for(int i = 0; i < bank.rows(); i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  ip = bank.getShort("component",i);
                int adc = Math.abs(bank.getInt("ADC",i));
                float t = bank.getFloat("time",i) + (float) tmf.getDoubleValue("offset",is,il,ip) 
                                                  + (float)  fo.getDoubleValue("offset",is,il,0);            
                
                String tag = ""+is+getDet(il)+(getLay(il)-1);
                
                float  sca = (float) ((is==5)?SCALE5[il-1]:SCALE[il-1]);
                float sadc = (float)(adc/sca); 
                
                float  tmax = 1000; float tdc = 0;
                
                double a0 = time.getDoubleValue("a0", is, il, ip); 
                double a2 = time.getDoubleValue("a2", is, il, ip);
                double a3 = time.getDoubleValue("a3", is, il, ip);
                double a4 = time.getDoubleValue("a4", is, il, ip);
                
                float ftdc_corr = t+FTOFFSET;
                
                if (tdcs.hasItem(is,il,ip)) {
                    float radc = (float)Math.sqrt(adc);
                    for (int tdcc : tdcs.getItem(is,il,ip)) {
                    	float tcor = tps*tdcc -  getTriggerPhase() - (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - FTOFFSET;
//                        float tdcmc = tdcm - a0 -  (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - a2/radc - a3 - a4/Math.sqrt(radc);
                        float tdif = tcor - t; 
                        if (Math.abs(tdif)<TMFCUT && tdif<tmax) {tmax = tdif; tdc = tcor;} 
                        dgm.fill("at1"+tag,sadc,tcor); 
                    }
                } 
                
                if(sadc>0) {
                	dgm.fill("at2"+tag,sadc,tdc);
                	dgm.fill("at3"+tag,sadc,tmax);
                	dgm.fill("at4"+tag,sadc,ip);                	
                }              
            }
        }
    }

}
