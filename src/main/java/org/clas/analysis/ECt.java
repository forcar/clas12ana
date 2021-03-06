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
import org.jlab.groot.data.DataLine;
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
    IndexedTable time=null, gtw=null, fo=null, fgo=null, tmf=null, tgo=null, gain=null, veff=null, rfT=null;

    int[]     npmt = {68,62,62,36,36,36,36,36,36};    
    int[]    npmts = new int[]{68,36,36};
    String[]   det = new String[]{"pcal","ecin","ecou"};
    String[]   pid = new String[]{"elec","pion","neut"};
    String[]     v = new String[]{"u","v","w"};  
    
    IndexedList<GraphErrors>  TDCSummary = new IndexedList<GraphErrors>(4);
    IndexedList<FitData>         TDCFits = new IndexedList<FitData>(4);
    Boolean                isAnalyzeDone = false;
    Boolean               isResidualDone = false;
    Boolean                    isTMFDone = false;
    Boolean                   isGTMFDone = false;
    Boolean           isTimeLineFitsDone = false;
    
    int trigger_sect = 0;
    int        phase = 0;
    float        STT = 0;
    
    static int tlnum;
    
    Boolean       isMC = false;    
    Boolean   isGoodTL = false;
    
//    static float TOFFSET = 600; 
   
    static float   FTOFFSET = 0;
    static float    TOFFSET = 0;
    static float  TOFFSETMC = 180;
    static float  TOFFSETER = 305;
    static float   RFPERIOD = 4.008f;    
    static float          c = 29.98f;
    static float   A0offset = 0;
    static float   A0sector = 0;
    
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
                                 "BETA",
                                 "RESID",
                                 "TMF",
                                 "GTMF",
                                 "LTFITS",
                                 "TWFITS",
                                 "TL",
                                 "Timeline",
                                 "EV",
                                 "T0");
        

        
        this.useCALUVWSECButtons(true);
        this.usePCCheckBox(true);
        this.useSliderPane(true);
        this.variation = "default";
        engine.setVeff(18.1f);
        engine.setNewTimeCal(true);
        engine.setPCALTrackingPlane(0);
        tlnum = getDetectorTabNames().indexOf("TL");
        this.init();
        this.localinit();
        this.localclear();
    }
    
    public ECt(String name, int runno) {
    	super(name);
    	initCCDB(runno);
    	
    }
    
    public void localinit() {
    	System.out.println("ECt.localinit()");
    }  
    
    public void localclear() {
    	System.out.println("ECt:localclear()");
    	isAnalyzeDone = false;
        isTimeLineFitsDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();    	   	
    	slider.setValue(0);
    } 
    
    @Override
    public void createHistos(int run) {  
	    histosExist = true;
	    System.out.println("ECt.createHistos("+run+")");
        setRunNumber(run);
        runlist.add(run);        
        this.setNumberOfEvents(0);  
        
        int t1=-70, t2=150, t3=50, t4=250;
        createTDCHistos(10,-10.,10.,"T-TVERT-PATH/c (ns)"); 
        createBETAHistos(22);
        createTLHistos(tlnum,t1,t2,t3,t4);
        
        if(dropSummary) return;
        createTDCHistos(0,t3,t4,"TIME (ns)");
        createTDCHistos(1,t3,t4,"TIME (ns)");
        createTDCHistos(2,t3,t4,"TIME (ns)");
        createTDCHistos(3,t3,t4,"TIME (ns)");
        createTDCHistos(4,t3,t4,"TIME (ns)");    
        createTDCHistos(5,t3,t4,"TIME (ns)");    
        createTDCHistos(6,t3,t4,"TIME (ns)");    
        createTDCHistos(7,t3,t4,"TIME (ns)");    
        createTDCHistos(8,-30.,25.,"TIME-FADC (ns)");    
        createTDCHistos(9,  0.,50.,"T-TVERT (ns)");    
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
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);
    	plotAnalysis(run);
    }
      
    public void plotSummary(int run) {  
    	    setRunNumber(run);
    	    plotTLHistos(tlnum);  	
    	    plotTLHistos(22);
    	    plotTDCHistos(10); 
            if(dropSummary) return;
    	    plotTDCHistos(0);
    	    plotTDCHistos(1);    	    	    
    	    plotTDCHistos(2);    	    	    
    	    plotTDCHistos(3);    	    	    
    	    plotTDCHistos(4);    	    	    
    	    plotTDCHistos(8);    	    	    
    	    plotTDCHistos(5);    	    	    
    	    plotTDCHistos(6);    	    	    
    	    plotTDCHistos(7);    	    	    
    	    plotTDCHistos(8);    	    	    
    	    plotTDCHistos(9);    	    	    
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
    }
    
    public void plotAnalysis(int run) {
    	    setRunNumber(run);
    	    if(!isAnalyzeDone) return;
	    	if(isResidualDone) {plotResidualSummary(23);}
    	    if(!dropSummary) {
    	    	if(isAnalyzeDone) {/*updateUVW(22)*/; updateFITS(26);updateFITS(27);plotEVSummary(30);plotT0Summary(31);}
//    	    	if(isResidualDone) plotResidualSummary(23);
    	    	if(isTMFDone)      plotTMFSummary(24);
    	    	if(isGTMFDone)     plotGTMFSummary(25);
    	    }
    	    if(!isTimeLineFitsDone) return;
    	    if(!isMC) plotTimeLines(29);
    }
    
    public void createBETAHistos(int k) {
        H1F h;
        int run = getRunNumber();
   
        DataGroup dg = new DataGroup(3,3);
        for (int id=1; id<4; id++) {
    	for (int il=1; il<4; il++) {
        for (int is=1; is<7 ; is++) {
           h = new H1F("beta_"+k+"_"+is+"_"+il+"_"+id+"_"+run,"beta_"+k+"_"+is+"_"+il+"_"+id+"_"+run,100,0.4,1.5);
           h.setTitleX(pid[id-1]+" "+det[il-1]+" beta "); h.setTitleY("Counts"); h.setLineColor(is==6?9:is);
           dg.addDataSet(h,id*3+il-4);
        }
    	}
        }
        this.getDataGroup().add(dg,0,0,k,run);  
                
    }
    
    public void createTLHistos(int k, int t1, int t2, int t3, int t4) {
        H1F h;
        int run = getRunNumber(); 
        
        DataGroup dg = new DataGroup(6,2);
        
        for (int is=1; is<7; is++) {
        	h = new H1F("TL1_"+k+"_"+is+"_"+run,"TL1_"+k+"_"+is+"_"+run,200,t1,t2);
        	h.setTitleX("Sector "+is+" Start Time-fgo"); h.setTitleY("Counts"); h.setOptStat("1100");
            dg.addDataSet(h,is-1);    	
        	h = new H1F("TL2_"+k+"_"+is+"_"+run,"TL2_"+k+"_"+is+"_"+run,100,t3,t4);
        	h.setTitleX("Sector "+is+" PCAL U3 Cluster Time-fgo"); h.setTitleY("Counts"); h.setOptStat("1100");     
            dg.addDataSet(h,is+5);
        }
        this.getDataGroup().add(dg,0,0,k,run);  
        createTDC1DHistos(k,-10.,10.,"T-TVERT-PATH/c (ns)"); 
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

    public void createTDC1DHistos(int k, double tmin, double tmax, String txt) {
    	
        int run = getRunNumber();
        H1F h;  
        
        for (int is=1; is<7; is++) {
            DataGroup dg = new DataGroup(3,3);
            h = new H1F("tdc1_pcal_u_"+is+"_"+k+"_"+run,"tdc1_pcal_u_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" PCAL U "+txt);
            dg.addDataSet(h,0);  
            h = new H1F("tdc1_pcal_v_"+is+"_"+k+"_"+run,"tdc1_pcal_v_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" PCAL V "+txt);       
            dg.addDataSet(h,1);            
            h = new H1F("tdc1_pcal_w_"+is+"_"+k+"_"+run,"tdc1_pcal_w_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" PCAL W "+txt);  
            dg.addDataSet(h,2); 
            
            h = new H1F("tdc1_ecin_u_"+is+"_"+k+"_"+run,"tdc1_ecin_u_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECIN U "+txt);  
            dg.addDataSet(h,3);  
            h = new H1F("tdc1_ecin_v_"+is+"_"+k+"_"+run,"tdc1_ecin_v_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECIN V "+txt);        
            dg.addDataSet(h,4);            
            h = new H1F("tdc1_ecin_w_"+is+"_"+k+"_"+run,"tdc1_ecin_w_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECIN W "+txt);   
            dg.addDataSet(h,5); 
            
            h = new H1F("tdc1_ecou_u_"+is+"_"+k+"_"+run,"tdc_ecou_u_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECOU U "+txt);   
            dg.addDataSet(h,6);  
            h = new H1F("tdc1_ecou_v_"+is+"_"+k+"_"+run,"tdc1_ecou_v_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECOU V "+txt);      
            dg.addDataSet(h,7);            
            h = new H1F("tdc1_ecou_w_"+is+"_"+k+"_"+run,"tdc1_ecou_w_"+is+"_"+k+"_"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECOU W "+txt); 
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,0,k,run);
        }            

    } 
    public void initCCDB(int runno) {
        gain    = cm.getConstants(runno, "/calibration/ec/gain");
        time    = cm.getConstants(runno, "/calibration/ec/timing");
        veff    = cm.getConstants(runno, "/calibration/ec/effective_velocity");
        fo      = cm.getConstants(runno, "/calibration/ec/fadc_offset");        //Crate/fiber FADC offsets 
        fgo     = cm.getConstants(runno, "/calibration/ec/fadc_global_offset"); //FADC capture window 
        tgo     = cm.getConstants(runno, "/calibration/ec/tdc_global_offset");  //TDC capture window
        gtw     = cm.getConstants(runno, "/calibration/ec/global_time_walk");   //Global time walk correction using raw ADC
        tmf     = cm.getConstants(runno, "/calibration/ec/tmf_offset");         //TDC-FADC offsets
        rfT     = ebcm.getConstants(runno, "/calibration/eb/rf/config");  
        
        FTOFFSET = (float) fgo.getDoubleValue("global_offset",0,0,0);
        TOFFSET  = (float) tgo.getDoubleValue("offset", 0,0,0);
        RFPERIOD = (float) rfT.getDoubleValue("clock",1,1,1);     
    }
    
    @Override
    public void processEvent(DataEvent event) {   
    	
       isMC = (getRunNumber()<100) ? true:false;
       
       trigger_sect = isMC ? (event.hasBank("ECAL::adc") ? event.getBank("ECAL::adc").getByte("sector",0):5) : getElecTriggerSector(); 
       boolean goodSector = trigger_sect>0 && trigger_sect<7; 
       
       if(!goodSector) return;
       
       phase  = getTriggerPhase(); 
       STT = event.hasBank("REC::Event") ? 
    		 (isHipo3Event ? 
    		 event.getBank("REC::Event").getFloat("STTime", 0):
             event.getBank("REC::Event").getFloat("startTime", 0)):
             0; 
       isGoodTL = event.hasBank("REC::Event")&&
    		      event.hasBank("REC::Particle")&&
    		      event.hasBank("REC::Calorimeter")&&
    		      STT>0;
    		     	   
       if(isGoodTL)     processTL(event); 
       if(!dropSummary) processRaw(event);
       
       processRec(event);       
       
    }
    
    public void processTL(DataEvent event) {
    	
  	   int  run = getRunNumber();
  	   ((H1F) this.getDataGroup().getItem(0,0,tlnum,run).getData(trigger_sect-1).get(0)).fill(STT-FTOFFSET);
  	   

  	   if(event.hasBank("ECAL::clusters")){       
       DataBank  bank1 = event.getBank("ECAL::clusters");
       for(int loop = 0; loop < bank1.rows(); loop++){
           int is = bank1.getByte("sector", loop);
           int il = bank1.getByte("layer", loop);
           if (is==trigger_sect&&il==1){
               float    t = bank1.getFloat("time",loop);
               int iU = (bank1.getInt("coordU", loop)-4)/8+1;
               if(iU==3) ((H1F) this.getDataGroup().getItem(0,0,tlnum,run).getData(is+5).get(0)).fill(t+TOFFSET-FTOFFSET);
            }
       }
       }
      
    }
    
    public void processRaw(DataEvent event) { //To cross-check ECengine for consistency
    	
 	   int  run = getRunNumber();
 	  
       tdcs.clear();
       
       if(event.hasBank("ECAL::tdc")==true){
           DataBank  bank = event.getBank("ECAL::tdc");
           for(int i = 0; i < bank.rows(); i++){
               int  is = bank.getByte("sector",i);
               int  il = bank.getByte("layer",i);
               int  ip = bank.getShort("component",i);   
               float tdcd  = bank.getInt("TDC",i)*tps; // raw time
               float tdcdc = tdcd-phase; // phase corrected time
               if(is>0&&is<7&&tdcd>0) {
                   if(!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Float>(),is,il,ip);    
              	       ((H2F) this.getDataGroup().getItem(is,0,0,run).getData(il-1).get(0)).fill(tdcd-FTOFFSET, ip); 
              	       ((H2F) this.getDataGroup().getItem(is,0,1,run).getData(il-1).get(0)).fill(tdcdc-FTOFFSET,ip); 
                       if (is==trigger_sect||isMC) {
                    	       tdcs.getItem(is,il,ip).add((float)tdcdc);
                  	       ((H2F) this.getDataGroup().getItem(is,0,2,run).getData(il-1).get(0)).fill(tdcdc-FTOFFSET,ip); // triggered time
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
               float t = bank.getFloat("time",i) + (float) tmf.getDoubleValue("offset",is,il,ip)  // FADC-TDC offset (sector, layer, PMT) 
                                                 + (float)  fo.getDoubleValue("offset",is,il,0)   // FADC-TDC offset (sector, layer) 
                                                 + FTOFFSET;
               
               float tmax = 1000; float tdcm = 1000;
              
               if (tdcs.hasItem(is,il,ip)) { // sector,layer,component FADC/TDC match
                 float radc = (float)Math.sqrt(adc);
                 for (float tdc : tdcs.getItem(is,il,ip)) {
                	    float tdif = tdc - (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - t;
             	      ((H2F) this.getDataGroup().getItem(is,0,8,run).getData(il-1).get(0)).fill(tdif,ip); // FADC t - TDC time
            	          if (Math.abs(tdif)<10 && tdif<tmax) {tmax = tdif; tdcm = tdc;}                	    
                   }
        	       double a0 = time.getDoubleValue("a0", is, il, ip); 
        	       double a1 = time.getDoubleValue("a1", is, il, ip);
        	       double a2 = time.getDoubleValue("a2", is, il, ip);
        	       double a3 = time.getDoubleValue("a3", is, il, ip);
        	       double a4 = time.getDoubleValue("a4", is, il, ip);
        	       double tdcmc = tdcm - a0 -  (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - a2/radc - a3 - a4/Math.sqrt(radc);
          	       ((H2F) this.getDataGroup().getItem(is,0,3,run).getData(il-1).get(0)).fill(tdcm-FTOFFSET,ip);  // matched FADC/TDC
          	       ((H2F) this.getDataGroup().getItem(is,0,4,run).getData(il-1).get(0)).fill(tdcmc-FTOFFSET,ip); // calibrated time
          	       if(isGoodTL&&il==1&&ip==3) ((H1F) this.getDataGroup().getItem(0,0,tlnum,run).getData(is+5).get(0)).fill(tdcmc); 
               }
           }
       }    	
    }
    
    public void processRec(DataEvent event) {
  	   
       float rftime = event.hasBank("REC::Event") ? event.getBank("REC::Event").getFloat("RFTime",0):0;

       DataBank recRunRF  = null; float trf = 0;

       if(event.hasBank("RUN::rf"))  {
    	   recRunRF = event.getBank("RUN::rf");             
           for(int k = 0; k < recRunRF.rows(); k++){
        	   if(recRunRF.getInt("id", k)==1) trf = recRunRF.getFloat("time",k);
           }    
       }
       
       if(isMC && STT<0) STT=124.25f;
       
       if(!(STT>-1)) return;
       
       int run = getRunNumber(); 
        
       if(dropBanks) dropBanks(event);
       
       List<ECStrip>     strips = engine.getStrips();
       List<ECPeak>       peaks = engine.getPeaks();
       List<ECCluster> clusters = engine.getClusters();
       
       if (run>=2000) {phase=0; shiftTV[1]=0;} // Corrections needed until runs<4013 are recooked
       
       if(!dropSummary) {
       if(event.hasBank("ECAL::hits")){
          	DataBank  bank = event.getBank("ECAL::hits");
            for(int loop = 0; loop < bank.rows(); loop++){
               int   is = bank.getByte("sector", loop);
               int   il = bank.getByte("layer", loop); 
               int   ip = bank.getByte("strip", loop);
               float  t = bank.getFloat("time", loop);
               ((H2F) this.getDataGroup().getItem(is,0,5,run).getData(il-1).get(0)).fill(t+TOFFSET-FTOFFSET, ip); //calibrated triggered matched hits
            }
       }
       }
       
       if(!(event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter"))) return;
       
       IndexedList<Integer> pathlist = new IndexedList<Integer>(3);    
       
       DataBank bankc = event.getBank("REC::Calorimeter");
       DataBank bankp = event.getBank("REC::Particle");
      
       Map<Integer,List<Integer>> caloMap = loadMapByIndex(bankc,"pindex");
       Map<Integer,List<Integer>> partMap = loadMapByIndex(bankp,"pid");    
       
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
             if (true) {
               int     il =  bank1.getByte("layer", loop);
               float ener =  bank1.getFloat("energy",loop)*1000;
               float    t =  bank1.getFloat("time",loop);
               iip[0]     = (bank1.getInt("coordU", loop)-4)/8+1;
               iip[1]     = (bank1.getInt("coordV", loop)-4)/8+1;
               iip[2]     = (bank1.getInt("coordW", loop)-4)/8+1;
               iid[0]     =  bank1.getInt("idU",loop)-1;
               iid[1]     =  bank1.getInt("idV",loop)-1;
               iid[2]     =  bank1.getInt("idW",loop)-1;	
               
               Point3D  pc = new Point3D(bank1.getFloat("x",loop),
            		                     bank1.getFloat("y",loop),
            		                     bank1.getFloat("z",loop));
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
                   int   stat = Math.abs(bankp.getInt("status",pin));
                   for (int i=0; i<3; i++) {  
                	   float tu=0,tdc=0,tdcc=0,tdccc=0,leff=0,adc=0; int ip=0;
                	   if (clusters.size()>0) {
                         tu    = (float) clusters.get(loop).getTime(i); 
                         ip    =         clusters.get(loop).getPeak(i).getMaxStrip();
                         adc   =         clusters.get(loop).getPeak(i).getMaxECStrip().getADC();
                         tdc   = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getRawTime(true)-TOFFSET;
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
                	   
                       float vel=c; 
                       if(Math.abs(pid)==211) vel=Math.abs(beta*c); //use EB beta for pion calibration residuals
                       
                	   float vcorr = STT - phase + TVOffset;  // phase=0 unless STT is not phase corrected (early engineering runs)
                	   float pcorr = path/vel;
                       float tvcor = tu  - vcorr;
                       float resid = tvcor - pcorr; 
                       
                       float mybet = path/tvcor/c; // use ECAL beta for beta distribution plots
                       
                	   float  radc = (float) Math.sqrt(adc);                	   
                       float lcorr = leff/(float)veff.getDoubleValue("veff", is, il+i, ip); 
                                              
//                       System.out.println((tu-pcorr)+" "+STT+" "+phase+" "+TVOffset+" "+rftime);
//                       double dt = (time - path/(beta*29.97) - trf + 120.5*this.rfPeriod)%this.rfPeriod-this.rfPeriod/2;
                       
                       float dt = 0;
                       if(recRunRF!=null) dt = (tu-pcorr-trf+120.5f*RFPERIOD)%RFPERIOD-RFPERIOD/2;   
                       
                       boolean goodSector = dropEsect?is!=trigger_sect:is==trigger_sect;
                       boolean goodStatus = stat>=2000 && stat<3000; 
                       boolean goodHisto  = Math.abs(pid)==TRpid && goodSector && goodStatus;
                       
                       if (goodHisto) {
                    	   ((H1F) this.getDataGroup().getItem(is,0,tlnum,run).getData(il+i-1).get(0)).fill(resid); //timelines
                           ((H2F) this.getDataGroup().getItem(is,0,10,run).getData(il+i-1).get(0)).fill(resid, ip);
                       }
                       
//                       ((H2F) this.getDataGroup().getItem(is,   0,10,run).getData(il+i-1).get(0)).fill(tdifp, ip);
                       if (!dropSummary && goodHisto) {  //Menu selection
                           ((H2F) this.getDataGroup().getItem(is,0,6,run).getData(il+i-1).get(0)).fill(tu+TOFFSET-FTOFFSET, ip); //peak times
                           ((H2F) this.getDataGroup().getItem(is,0,7,run).getData(il+i-1).get(0)).fill(t+TOFFSET-FTOFFSET,  ip); //cluster times
                           ((H2F) this.getDataGroup().getItem(is,0,9,run).getData(il+i-1).get(0)).fill(tvcor, ip);      // TVertex corrected time
                           ((H2F) this.getDataGroup().getItem(is,il+i,11,run).getData(ip-1).get(0)).fill(path, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,12,run).getData(ip-1).get(0)).fill(ener, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,13,run).getData(ip-1).get(0)).fill(leff, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,14,run).getData(ip-1).get(0)).fill(tu+TOFFSET,   resid);  
                           ((H2F) this.getDataGroup().getItem(is,il+i,15,run).getData(ip-1).get(0)).fill(adc,  resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,16,run).getData(ip-1).get(0)).fill(tdc- vcorr-pcorr-lcorr, radc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,17,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr-lcorr, radc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,18,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr-lcorr, adc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,19,run).getData(ip-1).get(0)).fill(tdc- vcorr-pcorr, leff);
                           ((H2F) this.getDataGroup().getItem(is,il+i,20,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr, leff);
                           ((H2F) this.getDataGroup().getItem(is,il+i,21,run).getData(ip-1).get(0)).fill(vcorr+TOFFSET+(isMC?50:0),leff);  	
                       } 
                       if(pid==11)                     ((H1F) this.getDataGroup().getItem(0,0,22,run).getData(getDet(il)).get(is-1)).fill(mybet);  
                       if(Math.abs(pid)==211)          ((H1F) this.getDataGroup().getItem(0,0,22,run).getData(getDet(il)+3).get(is-1)).fill(mybet);  
                       if(pid==22||pid==2112)          ((H1F) this.getDataGroup().getItem(0,0,22,run).getData(getDet(il)+6).get(is-1)).fill(mybet);                         
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
        int    off = (index-26)*200;
        
        c.clear();
        if (np==36) c.divide(6,6); if (np>36) c.divide(9,8);
        
        for (int ipp=0; ipp<np ; ipp++) {
        	int ip=ipp+off;
        	if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
        		g=tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).getGraph(); 
        		if(g.getDataSize(0)>0) {
                   c.cd(ipp); c.getPad(ipp); c.getPad(ipp).getAxisX().setRange(-1,g.getDataX(g.getDataSize(0)-1)*1.05);
                                  if(off!=200) {double min = tl.fitData.getItem(is,3*il+iv+1,ipp,getRunNumber()).p0;
                                                double max = g.getDataY(g.getDataSize(0)-1); double mn=Math.min(min, max); double mx=Math.max(min, max);                                                
                                                c.getPad(ipp).getAxisY().setRange(mn*0.95,mx*1.05);}
                                  if(off==200)  c.getPad(ipp).getAxisY().setRange(-5,5);
                   c.draw(g);
                   F1D f1 = new F1D("p0","[a]",0.,g.getDataX(g.getDataSize(0)-1)*1.05); f1.setParameter(0,0);
                   f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        		}
        	}
        }
    }
    
    @Override
    public void plotEvent(DataEvent de) {
       analyze();
    }

    public void analyze() {    
       System.out.println(getDetectorName()+".Analyze() ");
	   if(sfitEnable)  analyzeResiduals();
    
       if(!dropSummary) {
    	   if(cfitEnable)  analyzeCalibration();
    	   if(dfitEnable)  analyzeTMF();
    	   if(gdfitEnable) analyzeGTMF();
       }
       
       analyzeTimeLineFits();
       isAnalyzeDone = true;
       System.out.println(getDetectorName()+".Analyze() Finished");
    }
    
    public void analyzeTimeLineFits() {
    	cfitEnable = true;
    	fitTLGraphs(1,7,0,0,0,0); 
        if(!isTimeLineFitsDone) createTimeLineHistos();
    	fillTimeLineHisto();   	
        System.out.println("analyzeTimeLineFits Finished");
        isTimeLineFitsDone = true;      
    }
    
    public void fitTLGraphs(int is1, int is2, int id1, int id2, int il1, int il2) {    	
        int run = getRunNumber();	 
        for (int is=is1; is<is2; is++) {
//        	tl.fitData.add(fitEngine(((H1F)this.getDataGroup().getItem(0,0,tlnum,run).getData(is-1).get(0)),0,70,105,75,101,1.5,2.5),is,0,100,run); 
        	tl.fitData.add(fitEngine(((H1F)this.getDataGroup().getItem(0,0,tlnum,run).getData(is-1).get(0)),0,70,200,75,200,1.5,2.5),is,0,100,run); 
        	tl.fitData.add(fitEngine(((H1F)this.getDataGroup().getItem(0,0,tlnum,run).getData(is+5).get(0)),0,190,230,190,130,2.0,2.0),is,1,100,run); 
        	for (int il=0; il<9; il++) {
        		double f1 = il<3 ? 2.2:2 ; double f2 = il<3 ? 1.4:1.1;
            	tl.fitData.add(fitEngine(((H1F)this.getDataGroup().getItem(is,0,tlnum,run).getData(il).get(0)),0,-2,4,-10,10,f1,f2),is,10+il,100,run);    
        	}
        }
    }
    
    public void analyzeCalibration() {
        analyzeGraphs(1,7,0,3,0,3);
        System.out.println("analyzeCalibration Finished");
        writeFile("timing",1,7,0,3,0,3);
        writeFile("effective_velocity",1,7,0,3,0,3);    
        isAnalyzeDone = true;
    }
    
    public void analyzeResiduals() {
        getResidualSummary(1,7,0,3,0,3);
        System.out.println("analyzeResiduals Finished");
        writeFile("timing_update",1,7,0,3,0,3);
        isResidualDone = true;
    }
    
    public void analyzeTMF() {
        getTMFSummary(1,7,0,3,0,3);
        System.out.println("analyzeTMF Finished");
        writeFile("tmf_offset",1,7,0,3,0,3);
        isTMFDone = true;
    }
    
    public void analyzeGTMF() {
        getGTMFSummary(1,7,0,3,0,3);
        System.out.println("analyzeGTMF Finished");
        writeFile("fadc_offset",1,7,0,3,0,3);
        isGTMFDone = true;
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
            	   tl.fitData.add(fitEngine(g,13,20),is,3*il+iv+1,ip+200,run); //TW fits            	   
                }
             }
          }
       } 
    }
    
    public void getResidualSummary(int is1, int is2, int il1, int il2, int iv1, int iv2) {
    	
        ParallelSliceFitter fitter;
        GraphErrors g1,g2;
        
        int run=getRunNumber();
        
        for (int is=is1; is<is2; is++) {            
           for (int il=il1; il<il2; il++) {
              for (int iv=iv1; iv<iv2; iv++) {
                  fitter = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(is,0,10,run).getData(3*il+iv).get(0));
                  fitter.setRange(-10,10); fitter.fitSlicesY(); 
                  g1 = sliceToGraph(fitter.getMeanSlices(),il,iv); 
                  g2 = sliceToGraph(fitter.getSigmaSlices(),il,iv);
                  GraphErrors mean = new GraphErrors("RESIDUAL_"+is+"_"+il+" "+iv,g1.getVectorX().getArray(),
                		                                                          g1.getVectorY().getArray(),
                		                                                          new double[g1.getDataSize(0)],
                		                                                          g2.getVectorY().getArray()); 
                  mean.getAttributes().setTitleX("Sector "+is+" "+det[il]+" "+v[iv]); mean.getAttributes().setTitleY("");  
             	  FitSummary.add(mean, is,il,iv,run);
              }
           }
        } 
    }
    
    public void plotResidualSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));     
        
        int n = 0;        
        double ymin=-3f, ymax=3f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
         
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int il=0; il<3; il++) {
        for (int iv=0; iv<3; iv++) {           	
        	if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
            F1D f1 = new F1D("p0","[a]",0.,npmt[3*il+iv]); f1.setParameter(0,0);
            GraphErrors plot1 = FitSummary.getItem(is,il,iv,getRunNumber());
            plot1.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) plot1.getAttributes().setTitleY("RESIDUAL");
            plot1.getAttributes().setTitleX("SECTOR "+is+" "+det[il]+" "+v[iv].toUpperCase()+" PMT");
            n++; c.draw(plot1); 
            f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        	}
        }
        }
        }        
    }
    
    public void plotEVSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));     
        
        int n = 0;        
        double ymin=10f, ymax=30f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
         
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int il=0; il<3; il++) {
        for (int iv=0; iv<3; iv++) {
        	GraphErrors plot = new GraphErrors(); GraphErrors plotdb = new GraphErrors();
        	for (int ip=0; ip<npmt[3*il+iv]; ip++) {
    			float  y = (float) (tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())?tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p1:18);
    			float ye = (float) (tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())?tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p1e:0);
    			float ev = (float) Math.max(10,Math.min(22,y!=0?Math.abs(1/y):17f));
    			float eve = y!=0 && y>=10 && y<=22 ? ye/y/y:0;
        		plot.addPoint(ip+1, ev, 0, eve); 
        		float ydb = (float) veff.getDoubleValue("veff", is, 3*il+iv+1, ip+1);
        		plotdb.addPoint(ip+1, ydb, 0, 0);
        	}            
            plot.setMarkerColor(1); plot.setLineColor(1); plotdb.setMarkerColor(4); plotdb.setLineColor(4);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) plot.getAttributes().setTitleY("EV (cm/ns)");
            plot.getAttributes().setTitleX("SECTOR "+is+" "+det[il]+" "+v[iv].toUpperCase()+" PMT");
            n++; c.draw(plot); c.draw(plotdb,"same"); 
        }
        }
        }
             
    } 
    
    public void plotT0Summary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));     
        
        int n = 0;        
        double ymin=0f, ymax=40f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
         
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int il=0; il<3; il++) {
        for (int iv=0; iv<3; iv++) {
        	GraphErrors plot1 = new GraphErrors(); GraphErrors plot2 = new GraphErrors(); GraphErrors plot3 = new GraphErrors();
        	for (int ip=0; ip<npmt[3*il+iv]; ip++) {
        		float ydb = (float) time.getDoubleValue("a0", is, 3*il+iv+1, ip+1);        		
        		float y = 0, ye = 0, y1 = 0, y1e = 0;
        		if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
        			y1 =  (float) tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p0;
        			y1e = (float) tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p0e;  			
        		}
        		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
        			 y = (float) (FitSummary.getItem(is,il,iv,getRunNumber()).getDataY(ip)+ydb);
        			ye = (float)  FitSummary.getItem(is,il,iv,getRunNumber()).getDataEY(ip);
        		}
        		plot1.addPoint(ip+1,   y, 0, ye); 
        		plot2.addPoint(ip+1, ydb, 0, 0);
        		plot3.addPoint(ip+1,  y1, 0, y1e);
        	}            
            plot1.setMarkerColor(1); plot1.setLineColor(1); 
            plot2.setMarkerColor(4); plot2.setLineColor(4);
            plot3.setMarkerColor(2); plot3.setLineColor(2);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) plot1.getAttributes().setTitleY("T0 (ns)");
            plot1.getAttributes().setTitleX("SECTOR "+is+" "+det[il]+" "+v[iv].toUpperCase()+" PMT");
            n++; c.draw(plot1); c.draw(plot2,"same"); c.draw(plot3,"same"); 
        }
        }
        }
             
    } 
    
    public void getTMFSummary(int is1, int is2, int il1, int il2, int iv1, int iv2) {
        ParallelSliceFitter fitter;
        GraphErrors g1;
        
        int run=getRunNumber();
        
        for (int is=is1; is<is2; is++) {            
           for (int il=il1; il<il2; il++) {
              for (int iv=iv1; iv<iv2; iv++) {
                  fitter = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(is,0,8,run).getData(3*il+iv).get(0));
                  fitter.setRange(-20,15); fitter.fitSlicesY();
                  g1 = sliceToGraph(fitter.getMeanSlices(),il,iv); 
                  GraphErrors mean = new GraphErrors("TMF_"+is+"_"+il+" "+iv,g1.getVectorX().getArray(),
                		                                                     g1.getVectorY().getArray(),
                		                                                     new double[g1.getDataSize(0)],
                		                                                     new double[g1.getDataSize(0)]); 
                  mean.getAttributes().setTitleX("Sector "+is+" "+det[il]+" "+v[iv]); mean.getAttributes().setTitleY("");  
             	  FitSummary.add(mean, is+10,il,iv,run);
              }
           }
        } 
    }
    
    public void plotTMFSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));      
        int            n = 0;
        
        double ymin=-30f, ymax=30f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
         
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int il=0; il<3; il++) {
        for (int iv=0; iv<3; iv++) {            
            F1D f1 = new F1D("p0","[a]",0.,npmt[3*il+iv]); f1.setParameter(0,0);
            GraphErrors plot1 = FitSummary.getItem(is+10,il,iv,getRunNumber());
            plot1.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) plot1.getAttributes().setTitleY("TDC-FADC (ns)");
            plot1.getAttributes().setTitleX("SECTOR "+is+" "+det[il]+" "+v[iv].toUpperCase()+" PMT");
            n++; c.draw(plot1); 
            f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        }
        }
        }        
    }
    
    public void getGTMFSummary(int is1, int is2, int il1, int il2, int iv1, int iv2) {
    	
        int run=getRunNumber();
        
        for (int is=is1; is<is2; is++) {            
           GraphErrors g1 = new GraphErrors();
           for (int il=il1; il<il2; il++) {
              for (int iv=iv1; iv<iv2; iv++) {            	  
                  tl.fitData.add(fitEngine(((H2F)this.getDataGroup().getItem(is,0,8,run).getData(3*il+iv).get(0)).projectionX(),0,-10,10,-10,10),is+20,3*il+iv+1,0,run);
                  g1.addPoint(3*il+iv+1, tl.fitData.getItem(is+20,3*il+iv+1,0,run).mean,0,tl.fitData.getItem(is+20,3*il+iv+1,0,run).meane);
              }
           } 
           FitSummary.add(g1,is+20,0,0,run);
        } 
      
    }  
    
    public void plotGTMFSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));      
        int            n = 0;
        
        double ymin=-10f, ymax=10f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
         
        c.clear(); c.divide(3, 2);
        
        for (int is=1; is<7; is++) {          	
            F1D f1 = new F1D("p0","[a]",0.,10.); f1.setParameter(0,0);
            GraphErrors plot1 = FitSummary.getItem(is+20,0,0,getRunNumber());
            plot1.setMarkerColor(1); plot1.setMarkerSize(4);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            plot1.getAttributes().setTitleY("Global TMF Offset (ns)");
            plot1.getAttributes().setTitleX("SECTOR "+is+" "+" LAYERS");
            n++; c.draw(plot1); 
            f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        }        
    }
    
	public void writeFile(String table, int is1, int is2, int il1, int il2, int iv1, int iv2) {
		
		String line = new String();
		
		try { 
			File outputFile = new File(outPath+table+"_"+getRunNumber());
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("ECt.writefile("+table+")");

			for (int is=is1; is<is2; is++) {
				for (int il=il1; il<il2; il++ ) {
					for (int iv=iv1; iv<iv2; iv++) {
						for (int ip=0; ip<npmt[3*il+iv]; ip++) {
							switch (table) {
							case "A0offset":           line =  getNewA0(is,il,iv,ip); break;
							case "timing_update":      line =  getA0(is,il,iv,ip);  break;
							case "timing":             line =  getTW(is,il,iv,ip);  break;
							case "effective_velocity": line =  getEV(is,il,iv,ip);  break;
							case "tmf_offset":         line =  getTMF(is,il,iv,ip); break;  
							case "fadc_offset":        line =  getGTMF(is,il,iv,ip);  
							}
							if (line!=null) {
						      System.out.println(line);
						      outputBw.write(line);
						      outputBw.newLine();
							}
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
				+(time.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p1)+" "
				+(time.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p0)+" "
				+(time.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p2)+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0";
		}
		
	}
	
	public String getA0(int is, int il, int iv, int ip) {
		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(time.getDoubleValue("a0", is, 3*il+iv+1, ip+1)
				+FitSummary.getItem(is,il,iv,getRunNumber()).getDataY(ip))+" "  //residuals
				+" 0.02345 "
				+time.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+" "
				+time.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+" "
				+time.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0";
		}
		
	}
	
	public void setA0offset(int sector, float offset) {
		A0sector = sector;
		A0offset = offset;
	}
	
	public String getNewA0(int is, int il, int iv, int ip) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(time.getDoubleValue("a0", is, 3*il+iv+1, ip+1)
				+(is==A0sector ? A0offset:0))+" "
				+" 0.02345 "
				+time.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+" "
				+time.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+" "
				+time.getDoubleValue("a4", is, 3*il+iv+1, ip+1);		
	}
	
	public String getTMF(int is, int il, int iv, int ip) {
		float off = 0;
		if(FitSummary.hasItem(is+10,il,iv,getRunNumber())) {
			double ev = FitSummary.getItem(is+10,il,iv,getRunNumber()).getDataY(ip);
			off = (float) tmf.getDoubleValue("offset",is,3*il+iv+1,ip+1);
//			float off = (float)  fo.getDoubleValue("offset",is,3*il+iv+1,0);
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "+ (ev + off);
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 0.0";
		}		
	}
	
	public String getGTMF(int is, int il, int iv, int ip) {		
		if (ip>0) return null;
		float off = 0;
		if(tl.fitData.hasItem(is+20,3*il+iv+1,0,getRunNumber())) {
			off = (float) fo.getDoubleValue("offset",is,3*il+iv+1,0);
			double ev = tl.fitData.getItem(is+20,3*il+iv+1,0,getRunNumber()).mean;
		    return is+" "+(3*il+iv+1)+" "+ip+" "+ (ev + off);
		} else {
			return is+" "+(3*il+iv+1)+" "+ip+" 0.0";
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
	
    public void plotUVWHistos(int index) {
       int run = getRunNumber();
       drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index,run));	    
    }

    
    public void plotTDCHistos(int index) {
       int run = getRunNumber();
       drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),0,index,run));	    
    }
    
    public void plotTLHistos(int index) {
       int run = getRunNumber();
       drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActivePC()==1?getActiveSector():0,0,index,run));	    
    }  
    
/*   TIMELINES */
    
    public void createTimeLineHistos() {   
    	System.out.println("Initializing "+TLname+" timeline"); 
    	runIndex = 0;
    	tl.createTimeLineHisto(10,"StartTime","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto(20,"PCALU3","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 1,"PCAL U Peak Resid Mean","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 2,"PCAL V Peak Resid Mean","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 3,"PCAL W Peak Resid Mean","Sector",TLmax,6,1,7);    	
    	tl.createTimeLineHisto( 4,"ECIN U Peak Resid Mean","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 5,"ECIN V Peak Resid Mean","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 6,"ECIN W Peak Resid Mean","Sector",TLmax,6,1,7);    	
    	tl.createTimeLineHisto( 7,"ECOU U Peak Resid Mean","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 8,"ECOU V Peak Resid Mean","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 9,"ECOU W Peak Resid Mean","Sector",TLmax,6,1,7);    
    	tl.createTimeLineHisto( 11,"PCAL U Peak Resid Sigma","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 12,"PCAL V Peak Resid Sigma","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 13,"PCAL W Peak Resid Sigma","Sector",TLmax,6,1,7);    	
    	tl.createTimeLineHisto( 14,"ECIN U Peak Resid Sigma","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 15,"ECIN V Peak Resid Sigma","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 16,"ECIN W Peak Resid Sigma","Sector",TLmax,6,1,7);    	
    	tl.createTimeLineHisto( 17,"ECOU U Peak Resid Sigma","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 18,"ECOU V Peak Resid Sigma","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 19,"ECOU W Peak Resid Sigma","Sector",TLmax,6,1,7);    
    	
    }   
    
    public void fillTimeLineHisto() {
    	System.out.println("Filling "+TLname+" timeline"); 
        for (int is=1; is<7; is++) {
            float   y = (float) tl.fitData.getItem(is,0,100,getRunNumber()).mean; 
            float  ye = (float) tl.fitData.getItem(is,0,100,getRunNumber()).meane;			 
            ((H2F)tl.Timeline.getItem(10,0)).fill(runIndex,is,y);	
            ((H2F)tl.Timeline.getItem(10,1)).fill(runIndex,is,Math.abs(ye));  
            
            float  y1 = (float) tl.fitData.getItem(is,1,100,getRunNumber()).mean; 
            float y1e = (float) tl.fitData.getItem(is,1,100,getRunNumber()).meane;			 
            ((H2F)tl.Timeline.getItem(20,0)).fill(runIndex,is,y1);	
            ((H2F)tl.Timeline.getItem(20,1)).fill(runIndex,is,Math.abs(y1e)); 
        }
            
        for (int is=1; is<7; is++) {
  		  for (int il=0; il<3; il++) {	
  	        for (int iv=0; iv<3; iv++) {  
  			    float  y = (float) tl.fitData.getItem(is,3*il+iv+10,100,getRunNumber()).mean;
  			    float ye = (float) tl.fitData.getItem(is,3*il+iv+10,100,getRunNumber()).meane;
  			    ((H2F)tl.Timeline.getItem(3*il+iv+1,0)).fill(runIndex,is,y);	
  			    ((H2F)tl.Timeline.getItem(3*il+iv+1,1)).fill(runIndex,is,ye);	
  			    float  y1 = (float) tl.fitData.getItem(is,3*il+iv+10,100,getRunNumber()).sigma;
  			    float y1e = (float) tl.fitData.getItem(is,3*il+iv+10,100,getRunNumber()).sigmae;
  			    ((H2F)tl.Timeline.getItem(3*il+iv+11,0)).fill(runIndex,is,Math.abs(y1));	
  			    ((H2F)tl.Timeline.getItem(3*il+iv+11,1)).fill(runIndex,is,y1e);	
  	        }
  		  }
        } 
        runIndex++;
    } 
    
    public String getTLtag() {
    	switch (TRpid) {
    	case  11: return "el";
    	case 211: return "pi"; 
    	case  22: return "ph";
    	}
    	return "el";
    }
    
    public void saveTimelines() {
    	System.out.println("ECt: Saving timelines");
    	String tag = getTLtag();
    	saveTimeLine(10,0,100,"StartTime","TIME");
    	saveTimeLine(20,1,100,"PCALU3","TIME");
    	saveTimeLine(1,10,100,"PCAL"+tag+"tU","TIME");
    	saveTimeLine(2,11,100,"PCAL"+tag+"tV","TIME");
    	saveTimeLine(3,12,100,"PCAL"+tag+"tW","TIME");
    	saveTimeLine(4,13,100,"ECIN"+tag+"tU","TIME");
    	saveTimeLine(5,14,100,"ECIN"+tag+"tV","TIME");
    	saveTimeLine(6,15,100,"ECIN"+tag+"tW","TIME");
    	saveTimeLine(7,16,100,"ECOU"+tag+"tU","TIME");
    	saveTimeLine(8,17,100,"ECOU"+tag+"tV","TIME");
    	saveTimeLine(9,18,100,"ECOU"+tag+"tW","TIME");
    	saveTimeLine(11,10,100,"PCAL"+tag+"sU","SIG");
    	saveTimeLine(12,11,100,"PCAL"+tag+"sV","SIG");
    	saveTimeLine(13,12,100,"PCAL"+tag+"sW","SIG");
    	saveTimeLine(14,13,100,"ECIN"+tag+"sU","SIG");
    	saveTimeLine(15,14,100,"ECIN"+tag+"sV","SIG");
    	saveTimeLine(16,15,100,"ECIN"+tag+"sW","SIG");
    	saveTimeLine(17,16,100,"ECOU"+tag+"sU","SIG");
    	saveTimeLine(18,17,100,"ECOU"+tag+"sV","SIG");
    	saveTimeLine(19,18,100,"ECOU"+tag+"sW","SIG");
    }
    
    public void plotTimeLines(int index) {
    	if(TLflag) {plotTimeLineSectors(index);} else {
    	plotPeakTimeLines(index);}
    } 
    
    public void plotSttTimeLines(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); 
        int           is = getActiveSector(); 
        
    	double[] tlmin = {85,210,0.02}, tlmax= {200,214,.06};
    	float[] tlmean = {89.06f,212f,0.04f};
    	String[]   tit = {" StartTime"," PCAL U3 Time"," #sigma(T) / T"};
    	int[]      ind = {0,2,1};
       
    	DataLine line1 = new DataLine(0,is,  runIndex+1,is);                   line1.setLineColor(5);
    	DataLine line2 = new DataLine(0,is+1,runIndex+1,is+1);                 line2.setLineColor(5);
    	DataLine line3 = new DataLine(runIndexSlider,1,  runIndexSlider,7);    line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,runIndexSlider+1,7);  line4.setLineColor(5);
    	
        c.clear(); c.divide(3, 2); 

        for (int i=0; i<2; i++) { int i3=i*3;
            double min=tlmin[i]; double max=tlmax[i]; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(i3); c.getPad(i3).setAxisRange(0,runIndex,1,7); c.getPad(i3).setTitleFontSize(18); c.getPad(i3).getAxisZ().setRange(min,max);
    		c.draw((H2F)tl.Timeline.getItem(10*(i+1),0));c.draw(line1);c.draw(line2);c.draw(line3);c.draw(line4);
             
    		c.cd(i3+1); c.getPad(i3+1).setAxisRange(-0.5,runIndex,min,max); c.getPad(i3+1).setTitleFontSize(18);
    		drawTimeLine(c,is,10*(i+1),tlmean[i],"Sector "+is+tit[i]);
    		
    		FitData       fd = tl.fitData.getItem(is,i==1?1:0,100,getRunNumber()); float tlm=(float)fd.mean;
//    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(tlmean[i>1?0:i]*0.8,tlmean[i>1?0:i]*1.2,0.,fd.getGraph().getMax()*1.1);
    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(tlm*0.8,tlm*1.2,0.,fd.getGraph().getMax()*1.1);
    		fd.getHist().getAttributes().setOptStat("1000100");
    		DataLine line6 = new DataLine(tlmean[i],-50,tlmean[i],fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);
    		c.draw(fd.getHist()); c.draw(fd.getGraph(),"same"); c.draw(line6);
        }
    }
    
    public void plotPeakTimeLines(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));                
        int           is = getActiveSector();
        int           iv = getActiveView();
        int 		  pc = getActivePC();
        
        FitData       fd = null;
        
    	String  v[] = {" U "," V "," W "};
       	String  t[] = {" Time Resol", "Time Resid"};
        
    	DataLine line1 = new DataLine(0,is,  runIndex+1,is);                     line1.setLineColor(5);
    	DataLine line2 = new DataLine(0,is+1,runIndex+1,is+1);                   line2.setLineColor(5);
    	DataLine line3 = new DataLine(runIndexSlider,  0,  runIndexSlider,  7);  line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,0,  runIndexSlider+1,7);  line4.setLineColor(5);
    	
        c.clear(); c.divide(3, 3); 

        for (int il=0; il<3; il++) {int i3=il*3;
            double min=-10f,max=3f; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(i3); c.getPad(i3).setAxisRange(0,runIndex,1,7); c.getPad(i3).getAxisZ().setLog(false); c.getPad(i3).setTitleFontSize(18); c.getPad(i3).getAxisZ().setRange(min,max);
    		c.draw((H2F)tl.Timeline.getItem(3*il+iv+(pc==0?11:1),0));c.draw(line1);c.draw(line2);c.draw(line3);c.draw(line4);
             
    		c.cd(i3+1); c.getPad(i3+1).setAxisRange(-0.5,runIndex,min,max); c.getPad(i3+1).setTitleFontSize(18);
    		drawTimeLine(c,is,3*il+iv+(pc==0?11:1),0f,"Sector "+is+v[iv]+t[pc] );
    		
    		fd = tl.fitData.getItem(is,3*il+iv+10,100,getRunNumber());
    		
    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(fd.getHist().getXaxis().min(),fd.getHist().getXaxis().max(),0.,fd.getGraph().getMax());  
            fd.getHist().getAttributes().setOptStat("1000100");
            DataLine line6 = new DataLine(0,-50,0,fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);
            c.draw(fd.getHist()); c.draw(fd.getGraph(),"same");  c.draw(line6);
        }    	
    } 
    
    public void plotTimeLineSectors(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int pc = getActivePC();
        int iv = getActiveView();
        int il = getActiveLayer();
       	String  v[] = {" U "," V "," W "};
       	String  l[] = {" PCAL "," ECIN "," ECOU "};
       	String  t[] = {" Time Resol", "Time Resid"};
        
        c.clear(); c.divide(3, 2);
    	for (int is=1; is<7; is++) {
    		double min=pc==0?0:-2 ; double max=pc==0?4:2; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(is-1); c.getPad(is-1).setAxisRange(-0.5,runIndex,min,max); c.getPad(is-1).setTitleFontSize(18);
    		drawTimeLine(c,is,3*il+iv+(pc==0?11:1),pc==0?0.5f:0f,"Sector "+is+l[il]+v[iv]+t[pc]);
    	}
    }
    
    @Override
    public void timerUpdate() {
    	
    }
    
    public static void main(String[] args) {
    	
    	ECt ect = new ECt("ECt",2385);
    	ect.setA0offset(2, -40);
    	ect.writeFile("A0offset",1,7,0,3,0,3);    
    	
    }
     
}
