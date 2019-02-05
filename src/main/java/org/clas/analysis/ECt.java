package org.clas.analysis;

import java.awt.Point;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.swing.SwingUtilities;

import org.clas.tools.FitData;
import org.clas.tools.ParallelSliceFitter;
import org.clas.viewer.DetectorMonitor;
import org.dom4j.CDATA;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.data.DataVector;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
//import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.service.ec.ECCluster;
import org.jlab.service.ec.ECEngine;
import org.jlab.service.ec.ECPeak;
import org.jlab.service.ec.ECStrip;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;

public class ECt extends DetectorMonitor {

    IndexedList<List<Float>> tdcs = new IndexedList<List<Float>>(3);
    IndexedTable time=null, offset=null, goffset=null, gain=null, veff=null;

    int[]     npmt = {68,62,62,36,36,36,36,36,36};    
    int[]    npmts = new int[]{68,36,36};
    String[]   det = new String[]{"pcal","ecin","ecou"};
    String[]     v = new String[]{"u","v","w"};  
    
    IndexedList<GraphErrors>  TDCSummary = new IndexedList<GraphErrors>(4);
    IndexedList<FitData>         TDCFits = new IndexedList<FitData>(4);
    Boolean                isAnalyzeDone = false;
    
    int trigger_sect = 0;
    int        phase = 0;
    Boolean     isMC = false;
    
//  static float TOFFSET = 436; 
    static float   FTOFFSET = 0;
    static float    TOFFSET = 0;
    static float  TOFFSETMC = 180;
    static float          c = 29.98f;
    float               tps =  (float) 0.02345;
    float[] shiftTV = {0,40,0,0,0,0}; //Run 3050 t0 calibration
    float[] tw = {300,300,300,100,100,100,100,100,100};
    
    public ECt(String name) {
        super(name);
        this.setDetectorTabNames("Raw TDC",
                                 "PhaseCorr TDC",
                                 "Triggered TDC",
                                 "Matched TDC",
                                 "Calib TDC",
                                 "Hit Time",
                                 "Peak Time",
                                 "Cluster Time",
                                 "TIME-FADC",
                                 "TIME-TVERT",
                                 "RESID v STRIP",
                                 "RESID v PATH",
                                 "RESID v ENERGY",
                                 "RESID v LEFF",
                                 "RESID v TIME ",
                                 "RESID v ADC ",
                                 "ADC v TIME",
                                 "ADC v TTW",
                                 "ADCHI v TTW",
                                 "LEFF v TIME",
                                 "LEFF v TTW",
                                 "LEFF v TVERT",
                                 "TDIF",
                                 "MISC",
                                 "LTFITS",
                                 "TWFITS");

        
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.variation = "default";
        configEngine("muon");
        engine.setVeff(18.1f);
        engine.setNewTimeCal(true);
        this.init();
        this.localinit();
    }
    
    public void localinit() {
    	System.out.println("ECt.localinit()");
    	tl.setFitData(Fits); 
    }  
    
    public void localclear() {
    	System.out.println("ECt:localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	FitSummary.clear();
    	Fits.clear();
    	tl.Timeline.clear();    	   	
    	slider.setValue(0);
    } 
    
    @Override
    public void createHistos(int run) {  
        setRunNumber(run);
        runlist.add(run);        
        createTDCHistos(0,150,350,"TIME (ns)");
        createTDCHistos(1,150,350,"TIME (ns)");
        createTDCHistos(2,150,350,"TIME (ns)");
        createTDCHistos(3,150,350,"TIME (ns)");
        createTDCHistos(4,150,350,"TIME (ns)");    
        createTDCHistos(5,150,350,"TIME (ns)");    
        createTDCHistos(6,150,350,"TIME (ns)");    
        createTDCHistos(7,150,350,"TIME (ns)");    
        createTDCHistos(8,-50.,50.,"TIME-FADC (ns)");    
        createTDCHistos(9,  0.,50.,"T-TVERT (ns)");    
        createTDCHistos(10,-10.,10.,"T-TVERT-PATH/c (ns)"); 
        createUVWHistos(11,50,50,700,800,-5,5,"PATH","RESID ");
        createUVWHistos(12,50,50,0,1600,-5,5,"ENERGY","RESID ");
      
        createUVWHistos(13,50,50,0,430,-5,5,"LEFF","RESID ");
        createUVWHistos(14,50,50,190,220,-5,5,"T","RESID ");
        createUVWHistos(15,50,50,0,6000,-5,5,"ADC","RESID ");
//      createUVWHistos(16,0,50,0,6000,"T","ADC ");
//      createUVWHistos(17,0,50,0,6000,"TTW","ADC ");
        createUVWHistos(16,60,50,-10,50,0,200,"T","ADC ");
        createUVWHistos(17,60,50,-10,50,0,200,"TTW","ADC ");
        createUVWHistos(18,60,50,-10,50,0,12000,"TTW","ADCHI ");
        createUVWHistos(19,50,50,0,50,0,430,"T","LEFF ");
        createUVWHistos(20,50,50,0,50,0,430,"TTW","LEFF ");
        createUVWHistos(21,50,50,160,190,0,430,"TVERT","LEFF ");        
        createMISCHistos(23,0,300,"Start Time","Counts");
    }

    @Override        
    public void plotHistos(int run) {  
    	    setRunNumber(run);
    	    plotTDCHistos(0);
    	    plotTDCHistos(1);    	    	    
    	    plotTDCHistos(2);    	    	    
    	    plotTDCHistos(3);    	    	    
    	    plotTDCHistos(4);    	    	    
    	    plotTDCHistos(5);    	    	    
    	    plotTDCHistos(6);    	    	    
    	    plotTDCHistos(7);    	    	    
    	    plotTDCHistos(8);    	    	    
    	    plotTDCHistos(9);    	    	    
    	    plotTDCHistos(10); 
    	    plotUVWHistos(11);
    	    plotUVWHistos(12);
    	    plotUVWHistos(13);
    	    plotUVWHistos(14);
    	    plotUVWHistos(15);
    	    plotUVWHistos(16);
    	    plotUVWHistos(17);
    	    plotUVWHistos(18);
    	    plotUVWHistos(19);
    	    plotUVWHistos(20);
    	    plotUVWHistos(21);
    	    if(isAnalyzeDone) {/*updateUVW(22)*/; updateFITS(24);updateFITS(25);}
    	    plotMISCHistos(23);
    }
    
    public void createMISCHistos(int k, int xmin, int xmax, String xtxt, String ytxt) {
        H1F h;
        DataGroup dg1 = new DataGroup(1,1);
        int xbins = 100;
        int run = getRunNumber();
        h = new H1F("misc_"+k+"_"+run,"misc_"+k+"_"+run,xbins,xmin,xmax);
        h.setTitleX(xtxt); h.setTitleY(ytxt);       
        dg1.addDataSet(h,0);
        this.getDataGroup().add(dg1,0,0,k,run);        
    }
    
    public void createUVWHistos(int k, int xbins, int ybins, double xmin, double xmax, double ymin, double ymax, String xtxt, String ytxt) {
    	
        H2F h;  F1D f1; 
        
        double sca1=1.0,sca2=1.0; double xoff=0; boolean scaly=false; boolean scaly2=false;
        int run = getRunNumber();
        if (k==12) {sca1=0.7; sca2=0.6;}
        if (k==19||k==20) {scaly=true;}
        if (k==16||k==17) {scaly2=true;}
        
        for (int is=1; is<7; is++) {      
        	
            DataGroup dg1 = new DataGroup(9,8); DataGroup dg2 = new DataGroup(8,8); DataGroup dg3 = new DataGroup(8,8);        	    
            f1 = new F1D("p0"+is+1+k,"[a]",xmin,xmax); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);          
            for (int ip=1; ip<npmts[0]+1; ip++) {int uvw = (ip>52)?(52+(ip-52)*2):ip; double ymx=(scaly)?2*uvw*4.5*0.51:ymax;int ybns=(scaly)?((ip>6)?ybins*ip/npmts[0]:5):ybins;
                if(scaly2) {ymx=ymax*(1-0.4*ip/npmt[0]);}
                h = new H2F("uvw_pcal_u"+ip+"_s"+is+"_"+k+"_"+run,"uvw_pcal_u"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin,xmax,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"U"+ip);       
                dg1.addDataSet(h,ip-1); dg1.addDataSet(f1,ip-1);
                uvw = (ip>15)?(30+(ip-15)):2*ip;
                ymx=(scaly)?uvw*4.5*1.23:ymax;
                if(scaly2) {ymx=0.82*ymax*(0.75+0.25*ip/npmt[1]) ;}
                h = new H2F("uvw_pcal_v"+ip+"_s"+is+"_"+k+"_"+run,"uvw_pcal_v"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin,xmax,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"V"+ip);
                dg2.addDataSet(h,ip-1); dg2.addDataSet(f1,ip-1);
                h = new H2F("uvw_pcal_w"+ip+"_s"+is+"_"+k+"_"+run,"uvw_pcal_w"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin,xmax,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"W"+ip); 
                dg3.addDataSet(h,ip-1); dg3.addDataSet(f1,ip-1);
     	    }
            this.getDataGroup().add(dg1,is,1,k,run); this.getDataGroup().add(dg2,is,2,k,run); this.getDataGroup().add(dg3,is,3,k,run);
            
            DataGroup dg4 = new DataGroup(6,6); DataGroup dg5 = new DataGroup(6,6); DataGroup dg6 = new DataGroup(6,6);        	         	   
            f1 = new F1D("p0"+is+2+k,"[a]",xmin*sca1-xoff,xmax*sca1-xoff); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
     	    for (int ip=1; ip<npmts[1]+1; ip++) {double ymx=(scaly)?ymax*ip/npmts[1]:ymax;int ybns=(scaly)?((ip>4)?ybins*ip/npmts[1]:5):ybins;
                if(scaly2) {ymx=0.75*ymax*(1-0.47*ip/npmts[1]);}
                h = new H2F("uvw_ecin_u"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecin_u"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca1-xoff,xmax*sca1-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECIN "+xtxt);  h.setTitleY(ytxt+"U"+ip); 
                dg4.addDataSet(h,ip-1); dg4.addDataSet(f1,ip-1);
                if(scaly2) {ymx=0.55*ymax*(0.63+0.37*ip/npmts[1]) ;}
                h = new H2F("uvw_ecin_v"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecin_v"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca1-xoff,xmax*sca1-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"V"+ip); 
                dg5.addDataSet(h,ip-1); dg5.addDataSet(f1,ip-1);
                h = new H2F("uvw_ecin_w"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecin_w"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca1-xoff,xmax*sca1-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"W"+ip);
                dg6.addDataSet(h,ip-1); dg6.addDataSet(f1,ip-1);
                
     	    }
            this.getDataGroup().add(dg4,is,4,k,run); this.getDataGroup().add(dg5,is,5,k,run); this.getDataGroup().add(dg6,is,6,k,run);
     	   
            DataGroup dg7 = new DataGroup(6,6); DataGroup dg8 = new DataGroup(6,6); DataGroup dg9 = new DataGroup(6,6);        	         	   
            f1 = new F1D("p0"+is+3+k,"[a]",xmin*sca2-xoff,xmax*sca2-xoff); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
     	    for (int ip=1; ip<npmts[2]+1; ip++) {double ymx=(scaly)?ymax*ip/npmts[2]:ymax;int ybns=(scaly)?((ip>4)?ybins*ip/npmts[2]:5):ybins;
                if(scaly2) {ymx=0.75*ymax*(1-0.47*ip/npmts[2]);}
                h = new H2F("uvw_ecou_u"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecou_u"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca2-xoff,xmax*sca2-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECOU "+xtxt); h.setTitleY(ytxt+"U"+ip);
                dg7.addDataSet(h,ip-1); dg7.addDataSet(f1,ip-1);
                if(scaly2) {ymx=0.5*ymax*(0.63+0.37*ip/npmts[2]) ;}
                h = new H2F("uvw_ecou_v"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecou_v"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca2-xoff,xmax*sca2-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECOU "+xtxt);  h.setTitleY(ytxt+"V"+ip);
                dg8.addDataSet(h,ip-1); dg8.addDataSet(f1,ip-1);
                h = new H2F("uvw_ecou_w"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecou_w"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca2-xoff,xmax*sca2-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECOU "+xtxt); h.setTitleY(ytxt+"W"+ip);
                dg9.addDataSet(h,ip-1); dg9.addDataSet(f1,ip-1);    
     	    }
            this.getDataGroup().add(dg7,is,7,k,run); this.getDataGroup().add(dg8,is,8,k,run); this.getDataGroup().add(dg9,is,9,k,run);     	    
        }        
    }
  
    public void createTDCHistos(int k, double tmin, double tmax, String txt) {
    	
        int run = getRunNumber();
        H2F h;  
        
        for (int is=1; is<7; is++) {
            DataGroup dg = new DataGroup(3,3);
            h = new H2F("tdc_pcal_u_"+is+"_"+k+"_"+run,"tdc_pcal_u_"+is+"_"+k+"_"+run,100, tmin, tmax, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("tdc_pcal_v_"+is+"_"+k+"_"+run,"tdc_pcal_v_"+is+"_"+k+"_"+run,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("tdc_pcal_w_"+is+"_"+k+"_"+run,"tdc_pcal_w_"+is+"_"+k+"_"+run,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
            
            h = new H2F("tdc_ecin_u_"+is+"_"+k+"_"+run,"tdc_ecin_u_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("tdc_ecin_v_"+is+"_"+k+"_"+run,"tdc_ecin_v_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("tdc_ecin_w_"+is+"_"+k+"_"+run,"tdc_ecin_w_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
            
            h = new H2F("tdc_ecou_u_"+is+"_"+k+"_"+run,"tdc_ecou_u_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("tdc_ecou_v_"+is+"_"+k+"_"+run,"tdc_ecou_v_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("tdc_ecou_w_"+is+"_"+k+"_"+run,"tdc_ecou_w_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,0,k,run);
        }            

    } 
    
    public void initCCDB(int runno) {
        gain    = engine.getConstantsManager().getConstants(runno, "/calibration/ec/gain");
        time    = engine.getConstantsManager().getConstants(runno, "/calibration/ec/timing");
        veff    = engine.getConstantsManager().getConstants(runno, "/calibration/ec/effective_velocity");
        offset  = engine.getConstantsManager().getConstants(runno, "/calibration/ec/fadc_offset");
        goffset = engine.getConstantsManager().getConstants(runno, "/calibration/ec/fadc_global_offset");    	
    }
    
    @Override
    public void processEvent(DataEvent event) {   
       isMC = (getRunNumber()<100) ? true:false;
       trigger_sect = getElecTriggerSector(); 
       phase        = getTriggerPhase();
       processRaw(event);
       processRec(event);             
    }
    
    public void processRaw(DataEvent event) { //To cross-check ECengine for consistency
    	
 	   int run = getRunNumber();
       FTOFFSET = (float) goffset.getDoubleValue("global_offset",0,0,0);
        
       tdcs.clear();
       
       if(event.hasBank("ECAL::tdc")==true){
           DataBank  bank = event.getBank("ECAL::tdc");
           for(int i = 0; i < bank.rows(); i++){
               int  is = bank.getByte("sector",i);
               int  il = bank.getByte("layer",i);
               int  ip = bank.getShort("component",i);   
               float tdcd  = bank.getInt("TDC",i)*tps;
               float tdcdc = tdcd-phase;
               if(is>0&&is<7&&tdcd>0) {
                   if(!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Float>(),is,il,ip);    
              	       ((H2F) this.getDataGroup().getItem(is,0,0,run).getData(il-1).get(0)).fill(tdcd, ip); // raw time
              	       ((H2F) this.getDataGroup().getItem(is,0,1,run).getData(il-1).get(0)).fill(tdcdc,ip); // phase corrected time
                       if (true||is==trigger_sect||isMC) {
                    	       tdcs.getItem(is,il,ip).add((float)tdcdc);
                  	       ((H2F) this.getDataGroup().getItem(is,0,2,run).getData(il-1).get(0)).fill(tdcdc,ip); // triggered time
                       }
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
               float t = bank.getFloat("time",i) + (float) offset.getDoubleValue("offset",is,il,0);  // FADC-TDC offset (sector, UVW layer)              
               
               float tmax = 1000; float tdcm = 1000;
               
               if (tdcs.hasItem(is,il,ip)) { // sector,layer,component FADC/TDC match
                   for (float tdc : tdcs.getItem(is,il,ip)) {
                	    float tdif = tdc-FTOFFSET-t;
             	      ((H2F) this.getDataGroup().getItem(is,0,8,run).getData(il-1).get(0)).fill(tdif,ip); // FADC t - TDC time
            	          if (Math.abs(tdif)<20 && tdif<tmax) {tmax = tdif; tdcm = tdc;}                	    
                   }
                   float radc = (float)Math.sqrt(adc);
        	       double a0 = time.getDoubleValue("a0", is, il, ip); 
        	       double a1 = time.getDoubleValue("a1", is, il, ip);
        	       double a2 = time.getDoubleValue("a2", is, il, ip);
        	       double a3 = time.getDoubleValue("a3", is, il, ip);
        	       double a4 = time.getDoubleValue("a4", is, il, ip);
        	       double tdcmc = tdcm - a0 - tw[il-1]/radc - a2 - a3/radc - a4/Math.sqrt(radc);
//        	       System.out.println(a0+" "+a1+" "+a2+" "+a3+" "+a4);
//        	       System.out.println(tdcm+" "+adc+" "+tw[il-1]/radc+" "+tdcmc);
//        	       System.out.println(" ");
//                   double tdcmc = tdcm - a0 - a2 / Math.sqrt(adc);
          	       ((H2F) this.getDataGroup().getItem(is,0,3,run).getData(il-1).get(0)).fill(tdcm,ip);  // matched FADC/TDC
          	       ((H2F) this.getDataGroup().getItem(is,0,4,run).getData(il-1).get(0)).fill(tdcmc,ip); // calibrated time
               }
           }
       }    	
    }
    
    public void processRec(DataEvent event) {
  	   
       float Tvertex = event.hasBank("REC::Event") ? event.getBank("REC::Event").getFloat("STTime", 0):0;
        
       if(!(Tvertex>0)) return;
       
       int run = getRunNumber();
        
       if(dropBanks) dropBanks(event);
       
       List<ECStrip>     strips = engine.getStrips();
       List<ECPeak>       peaks = engine.getPeaks();
       List<ECCluster> clusters = engine.getClusters();
       
       if (run>=3818) {phase=0; shiftTV[1]=0;} // Corrections needed until runs<4013 are recooked
       
       if(event.hasBank("ECAL::hits")){
          	DataBank  bank = event.getBank("ECAL::hits");
            for(int loop = 0; loop < bank.rows(); loop++){
               int   is = bank.getByte("sector", loop);
               int   il = bank.getByte("layer", loop); 
               int   ip = bank.getByte("strip", loop);
               float  t = bank.getFloat("time", loop);
               if (true||is==trigger_sect||isMC) {
            	      ((H2F) this.getDataGroup().getItem(is,0,5,run).getData(il-1).get(0)).fill(t, ip); //calibrated triggered matched hits
               }
            }
       }
       
       boolean goodEvent = event.hasBank("REC::Particle")&&event.hasBank("REC::Calorimeter");
       
       if(!goodEvent) return;
       
       ((H1F) this.getDataGroup().getItem(0,0,23,run).getData(0).get(0)).fill(Tvertex);  
       
       IndexedList<Integer> pathlist = new IndexedList<Integer>(3);    
       
       DataBank bankc = event.getBank("REC::Calorimeter");
       DataBank bankp = event.getBank("REC::Particle");
      
       Map<Integer,List<Integer>> caloMap = loadMapByIndex(bankc,"pindex");
       Map<Integer,List<Integer>> partMap = loadMapByIndex(bankp,"pid");    
       
       trigger_sect = getElecTriggerSector(); 
//       if (!(trigger_sect>0)) return;
       
//       if(!partMap.containsKey(11)) return;
      
       for(int loop = 0; loop < bankc.rows(); loop++){
          int   is = bankc.getByte("sector", loop);
          int   il = bankc.getByte("layer", loop);
          int   in = bankc.getShort("index", loop);
          int  det = bankc.getByte("detector", loop);
          if (det==7 && !pathlist.hasItem(is,il,in)) pathlist.add(loop,is,il,in);                 
       }
       
       int[] iip = new int[3];int[] iid = new int[3];
       float[] tid = new float[3]; float[] lef = new float[3]; float[] add = new float[3];
       float x=0,y=0,z=0;
       
       if(event.hasBank("ECAL::clusters")){       
       DataBank  bank1 = event.getBank("ECAL::clusters");
       DataBank  bank2 = event.getBank("ECAL::calib");
       DataBank  bank3 = event.getBank("ECAL::peaks");
       for(int loop = 0; loop < bank1.rows(); loop++){
           int is = bank1.getByte("sector", loop);
           if (true||is==trigger_sect||isMC){
               int     il = bank1.getByte("layer", loop);
               float ener = bank1.getFloat("energy",loop)*1000;
               float    t = bank1.getFloat("time",loop);
               iip[0] = (bank1.getInt("coordU", loop)-4)/8+1;
               iip[1] = (bank1.getInt("coordV", loop)-4)/8+1;
               iip[2] = (bank1.getInt("coordW", loop)-4)/8+1;
               iid[0] = bank1.getInt("idU",loop)-1;
               iid[1] = bank1.getInt("idV",loop)-1;
               iid[2] = bank1.getInt("idW",loop)-1;	
               
               Point3D  pc = new Point3D(bank1.getFloat("x",loop),
            		                     bank1.getFloat("y",loop),
            		                     bank1.getFloat("z",loop));
//               System.out.println("Cluster "+loop+" "+pc.toString()+" "+clusters.get(loop).getHitPosition().toString());
               lef[0] = getLeff(pc,getPeakline(iid[0],pc,bank3));
               lef[1] = getLeff(pc,getPeakline(iid[1],pc,bank3));
               lef[2] = getLeff(pc,getPeakline(iid[2],pc,bank3));
               tid[0] = bank3.getFloat("time",iid[0])-lef[0]/(float)veff.getDoubleValue("veff", is,il+0,iip[0]);
               tid[1] = bank3.getFloat("time",iid[1])-lef[1]/(float)veff.getDoubleValue("veff", is,il+1,iip[1]);
               tid[2] = bank3.getFloat("time",iid[2])-lef[2]/(float)veff.getDoubleValue("veff", is,il+2,iip[2]);
               add[0] = bank3.getFloat("energy",iid[0]);
               add[1] = bank3.getFloat("energy",iid[1]);
               add[2] = bank3.getFloat("energy",iid[2]);
               if (pathlist.hasItem(is,il,loop)) {
                   int    pin = bankc.getShort("pindex", pathlist.getItem(is,il,loop));
                   float path = bankc.getFloat("path",   pathlist.getItem(is,il,loop));
                   int    pid = bankp.getInt("pid",pin);
                   float beta = bankp.getFloat("beta",pin);  
                   for (int i=0; i<3; i++) {  
                	   float tu=0,tdc=0,tdcc=0,tdccc=0,leff=0,adc=0; int ip=0;
                	   if (clusters.size()>0) {
                         tu    = (float) clusters.get(loop).getTime(i); 
                         ip    =         clusters.get(loop).getPeak(i).getMaxStrip();
                         adc   =         clusters.get(loop).getPeak(i).getMaxECStrip().getADC();
                         tdc   = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getRawTime(true);
                         tdcc  = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTWCTime();
//                         tdcc = tdc-tw[il+i-1]/(float)Math.sqrt(adc);
                         tdccc = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTime(); 
                         leff  = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTdist();
//                         System.out.println("C Layer "+il+" "+clusters.get(loop).getPeak(i).getMaxECStrip().getLine().toString());
//                         System.out.println("C Layer "+il+" "+clusters.get(loop).getPeak(i).getLine().toString());
//                         System.out.println("P Layer "+il+" "+peaks.get(iid[i]).getLine().toString());
//                         System.out.println("P Layer "+il+" "+getPeakline(iid[i],pc,bank3).toString());
//                         System.out.println("Layer "+il+" "+tu+" "+tid[i]+" "+leff+ " "+lef[i]);
                	   } else {
                		 tu    = tid[i];
                		 ip    = iip[i];
                		 leff  = lef[i];
                         adc   = 10000*add[i]/(float)gain.getDoubleValue("gain", is, il+i, ip);
                	   }
                	   
                	   float radc = (float) Math.sqrt(adc);
                	   
                       float vel=c; if(Math.abs(pid)==211) vel=Math.abs(beta*c);
               	   
                	   float vcorr = Tvertex - phase;
                	   float pcorr = path/vel;
                	   float lcorr = leff/(float)veff.getDoubleValue("veff", is, il+i, ip);
                       float tdif  = tu  - vcorr;
                       float resid = tdif - pcorr;
                       
                       ((H2F) this.getDataGroup().getItem(is,0,6,run).getData(il+i-1).get(0)).fill(tu, ip); //peak times
                       ((H2F) this.getDataGroup().getItem(is,0,7,run).getData(il+i-1).get(0)).fill(t,  ip); //cluster times
                       ((H2F) this.getDataGroup().getItem(is,0,9,run).getData(il+i-1).get(0)).fill(tdif, ip);
//                       ((H2F) this.getDataGroup().getItem(is,   0,10,run).getData(il+i-1).get(0)).fill(tdifp, ip);
//                       if (pid==22) {
                           if (true||pid==11||Math.abs(pid)==211) {
                           ((H2F) this.getDataGroup().getItem(is,   0,10,run).getData(il+i-1).get(0)).fill(resid, ip);
                           ((H2F) this.getDataGroup().getItem(is,il+i,11,run).getData(ip-1).get(0)).fill(path, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,12,run).getData(ip-1).get(0)).fill(ener, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,13,run).getData(ip-1).get(0)).fill(leff, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,14,run).getData(ip-1).get(0)).fill(tu,   resid);  
                           ((H2F) this.getDataGroup().getItem(is,il+i,15,run).getData(ip-1).get(0)).fill(adc,  resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,16,run).getData(ip-1).get(0)).fill(tdc- vcorr-pcorr-lcorr, radc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,17,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr-lcorr, radc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,18,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr-lcorr, adc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,19,run).getData(ip-1).get(0)).fill(tdc- vcorr-pcorr, leff);
                           ((H2F) this.getDataGroup().getItem(is,il+i,20,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr, leff);
                           ((H2F) this.getDataGroup().getItem(is,il+i,21,run).getData(ip-1).get(0)).fill(vcorr,leff);  	
                       } 	      
                   } 
               }
           }
       }
       }
   	
    }
    
    public Line3D getPeakline(int iid, Point3D point, DataBank bank) {    	
        Point3D  po = new Point3D(bank.getFloat("xo",iid),
                                  bank.getFloat("yo",iid),
                                  bank.getFloat("zo",iid));
        Point3D  pe = new Point3D(bank.getFloat("xe",iid),
                                  bank.getFloat("ye",iid),
                                  bank.getFloat("ze",iid));
        return  new Line3D(po,pe);
    }
    
    public float getLeff(Point3D point, Line3D peakline) {
    	return (float) point.distance(peakline.end());
    }

    private void updateUVW(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.divide(3, 3);
        
        int    is = getActiveSector(); 
               
        for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            c.cd(3*i+j); c.getPad(3*i+j).getAxisY().setLog(false); 
            c.draw(TDCFits.getItem(is,i,j,0).getGraph());
        }
        }
        
    }  
    
    public void updateFITS(int index) {
        
    	GraphErrors g = null;
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
       
        int    is = getActiveSector(); 
        int    il = getActiveLayer();
        int    iv = getActiveView();
        
        int    np = npmt[3*il+iv];
        int    off = (index-24)*100;
        
        c.clear();
        if (np==36) c.divide(6,6); if (np>36) c.divide(9,8);
        
        for (int ipp=0; ipp<np ; ipp++) {
        	int ip=ipp+off;
        	if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
        		g=tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).getGraph();        	
        		if(g.getDataSize(0)>0) {
                   c.cd(ipp); c.getPad(ipp); c.getPad(ipp).getAxisX().setRange(-1,g.getDataX(g.getDataSize(0)-1)*1.05);
                                  if(off!=100) {double min = tl.fitData.getItem(is,3*il+iv+1,ipp,getRunNumber()).p0;
                                                double max = g.getDataY(g.getDataSize(0)-1); double mn=Math.min(min, max); double mx=Math.max(min, max);                                                
                                                c.getPad(ipp).getAxisY().setRange(mn*0.95,mx*1.05);}
                                  if(off==100)  c.getPad(ipp).getAxisY().setRange(-5,5);
                   c.draw(g);
                   F1D f1 = new F1D("p0","[a]",0.,g.getDataX(g.getDataSize(0)-1)*1.05); f1.setParameter(0,0);
                   f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        		}
        	}
        }
    }
    
	public void writeFile(String table, int is1, int is2, int il1, int il2, int iv1, int iv2) {
		
		String path = "/Users/colesmith/CLAS12ANA/";
		String line = new String();
		
		try { 
			File outputFile = new File(path+table+"_"+getRunNumber());
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("ECt.writefile("+table+")");

			for (int is=is1; is<is2; is++) {
				for (int il=il1; il<il2; il++ ) {
					for (int iv=iv1; iv<iv2; iv++) {
						for (int ip=0; ip<npmt[3*il+iv]; ip++) {
							switch (table) {
							case "timing":             line =  getTW(is,il,iv,ip); break;
							case "effective_velocity": line =  getEV(is,il,iv,ip);  
							}
						    System.out.println(line);
						    outputBw.write(line);
						    outputBw.newLine();
						}
					}
				}
			}

			outputBw.close();
			outputFw.close();
		}
		catch(IOException ex) {
			System.out.println("Error writing file '" );                   
			ex.printStackTrace();
		}

	}
	
	public String getTW(int is, int il, int iv, int ip) {
		if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p0
				+" 0.02345 "
				+(time.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+100,getRunNumber()).p0)+" "
				+(time.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+100,getRunNumber()).p1)+" "
				+(time.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+100,getRunNumber()).p2)+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0";
		}
		
	}
	
	public String getEV(int is, int il, int iv, int ip) {
		if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
			double ev = tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p1;
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "+Math.max(10,Math.min(22,ev!=0?Math.abs(1/ev):17.0));
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 17.0";
		}		
	}
    
    @Override
    public void plotEvent(DataEvent de) {
       analyze();
    }

    public void analyze() {    
       System.out.println("I am in analyze()");
       if(!fitEnable) return;
       analyzeGraphs(1,7,0,3,0,3);
       System.out.println("Finished");
       writeFile("timing",1,7,0,3,0,3);
       writeFile("effective_velocity",1,7,0,3,0,3);
       isAnalyzeDone = true;
    }
    
    public GraphErrors graphShift(GraphErrors g, double shift) {
    	GraphErrors gnew = new GraphErrors();      	
    	for (int i=0; i<g.getDataSize(0); i++) {
           gnew.addPoint(g.getDataX(i),g.getDataY(i)+shift,g.getDataEX(i),g.getDataEY(i));
    	}
    	return gnew;
    }
    
    public void analyzeGraphs(int is1, int is2, int il1, int il2, int iv1, int iv2) {
        
       ParallelSliceFitter fitter1, fitter2;
       GraphErrors g;
       
       int run=getRunNumber();
       double min=0,max=420;
       for (int is=is1; is<is2; is++) {            
          for (int il=il1; il<il2; il++) {
             for (int iv=iv1; iv<iv2; iv++) {
                for (int ip=0; ip<npmt[3*il+iv]; ip++) {
                   System.out.println("Fitting Sector "+is+" Layer "+il+" View "+iv+" PMT "+ip+" "+this.getDataGroup().hasItem(is,3*il+iv+1,20,run)); 
                   
                   fitter1 = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(is,3*il+iv+1,20,run).getData(ip).get(0));
                   fitter1.setRange(0,50); fitter1.fitSlicesY();
                   g = fitter1.getMeanSlices(); 
                   g.getAttributes().setTitleX("Sector "+is+" "+det[il]+" "+v[iv]+(ip+1)); g.getAttributes().setTitleY("");
               	   tl.fitData.add(fitEngine(g,6,0),is,3*il+iv+1,ip,run); //LEFF fits
              	   
                   fitter2 = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(is,3*il+iv+1,17,run).getData(ip).get(0));
                   fitter2.setRange(-10,50); fitter2.fitSlicesY();
                   g = graphShift(fitter2.getMeanSlices(),-tl.fitData.getItem(is,3*il+iv+1,ip,run).p0);
                   g.getAttributes().setTitleX("Sector "+is+" "+det[il]+" "+v[iv]+(ip+1)); g.getAttributes().setTitleY("");  
            	   tl.fitData.add(fitEngine(g,13,20),is,3*il+iv+1,ip+100,run); //TW fits
            	   
                }
             }
          }
       } 
    }  
    
    public void plotUVWHistos(int index) {
       int run = getRunNumber();
       drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index,run));	    
    }

    
    public void plotTDCHistos(int index) {
       int run = getRunNumber();
       drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),0,index,run));	    
    }
    
    public void plotMISCHistos(int index) {
       int run = getRunNumber();
       drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,run));	    
    }  
    
    @Override
    public void timerUpdate() {
    	
    }
     
}
