package org.clas.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;

import org.clas.tools.FitData;
import org.clas.tools.ParallelSliceFitter;
import org.clas.viewer.DetectorMonitor;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;

import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;

import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import org.clas.service.ec.ECCluster;
import org.clas.service.ec.ECPeak;
import org.clas.service.ec.ECStrip;

import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECt extends DetectorMonitor {

    IndexedList<List<Float>>  tdcs = new IndexedList<List<Float>>(3);
    IndexedList<List<Float>> tlist = new IndexedList<List<Float>>(3);
    IndexedTable time=null, oftime=null, ftime=null, dtime=null, gtw=null, fo=null, fgo=null, tmf=null, tgo=null, gain=null, veff=null, rfT=null;
    
    int[]     npmt = {68,62,62,36,36,36,36,36,36};    
    int[]    npmts = new int[]{68,36,36};
    String[]   det = new String[]{"pcal","ecin","ecou"};
    String[]   pid = new String[]{"elec","pion","neut"};
    String[]     v = new String[]{"u","v","w"};  
    
    IndexedList<GraphErrors>  TDCSummary = new IndexedList<GraphErrors>(4);
    IndexedList<FitData>         TDCFits = new IndexedList<FitData>(4);
    IndexedList<Integer>        pathlist = new IndexedList<Integer>(3);    

    Boolean                isAnalyzeDone = false;
    Boolean               isResidualDone = false;
    Boolean                    isTMFDone = false;
    Boolean                   isGTMFDone = false;
    Boolean           isTimeLineFitsDone = false;
    
    int trigger_sect = 0;
    int        phase = 0;
    float        STT = 0;
    float         RF = 0;
    float    RF_TIME = 124.25f;
    
    static int tlnum;
    
    Boolean       isMC = false;    
    Boolean   isGoodTL = false;
    
//    static float TOFFSET = 600; 
   
    static float    BGOFFSET = 102;
    static float    FTOFFSET = 0;
    static float     TOFFSET = 0;
    static float   TOFFSETMC = 180;
    static float   TOFFSETER = 305;
    static float    RFPERIOD = 4.008f;
    static float BEAM_BUCKET = 2.004f;
    static float           c = 29.98f;
    static float    A0offset = 0f;  //RGM pass0 only
    static float    A0sector = 0;
    
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
                                 "Calib TMF",
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
        this.useECEnginePane(true);
        tlnum = getDetectorTabNames().indexOf("TL");
        this.init();
        this.localinit();
    }
    
    public ECt(String name, int runno) {
    	super(name);
    	initCCDB(runno);
    	
    }
    
    public void localinit() {
    	System.out.println(getDetectorName()+".localinit()");
    	tl.setFitData(Fits);
    	eng.engine.setGeomVariation("rga_spring2018");
    }  
    
    public void localclear() {
    	System.out.println(getDetectorName()+".localclear()");
    	isAnalyzeDone = false;
        isTimeLineFitsDone = false;
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
	    histosExist = true;
        setRunNumber(run);  runlist.add(run);    
        
        this.setNumberOfEvents(0);  
        
        int t1=-70, t2=150, t3=50, t4=250;
        createTDCHistos(10,-30.,30.,"T-TVERT-PATH/#beta*c (ns)"); 
        createBETAHistos(22);
        createTLHistos(tlnum,t1,t2,t3,t4);
        
        if(dropSummary) return;
        createTDCHistos(0,t3,t4,"TIME (ns)");
        createTDCHistos(1,t3,t4,"TIME (ns)");
        createTDCHistos(2,t3,t4,"TIME (ns)");
        createTDCHistos(3,t3,t4,"TIME (ns)");
        createTDCHistos(4,t3,t4,"TIME (ns)");    
        createTDCHistos(5,-15,15,"TIME (ns)");    
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
	    if(isResidualDone) plotResidualSummary(23);
        if(!dropSummary) {
        	updateFITS(26); updateFITS(27); plotEVSummary(30); plotT0Summary(31);
        	if(isTMFDone)      plotTMFSummary(24);
        	if(isGTMFDone)     plotGTMFSummary(25);
    	}
        if(!isTimeLineFitsDone) return;
        if(!isMC) plotTimeLines(29);
    }
    
    public void createBETAHistos(int k) {

        int run = getRunNumber();
        H1F h;
   
        DataGroup dg = new DataGroup(3,3);
        for (int id=1; id<4; id++) {
    	for (int il=1; il<4; il++) {
        for (int is=1; is<7 ; is++) {
           h = new H1F("beta-"+k+"-"+is+"-"+il+"-"+id+"-"+run,"beta-"+k+"-"+is+"-"+il+"-"+id+"-"+run,100,0.4,1.5);
           h.setTitleX(pid[id-1]+" "+det[il-1]+" beta "); h.setTitleY("Counts"); h.setLineColor(is==6?9:is);
           dg.addDataSet(h,id*3+il-4);
        }
    	}
        }
        this.getDataGroup().add(dg,0,0,k,run);  
                
    }
    
    public void createTLHistos(int k, int t1, int t2, int t3, int t4) {

        int run = getRunNumber(); 
        H1F h; 
        H2F h2;
        
        DataGroup dg = new DataGroup(6,3);

        for (int is=1; is<7; is++) {
        	h = new H1F("TL1_"+k+"-"+is+"-"+run,"TL1_"+k+"-"+is+"-"+run,200,t1,t2);
        	h.setTitleX("Sector "+is+" Start Time-fgo"); h.setTitleY("Counts"); h.setOptStat("1100");
            dg.addDataSet(h,is-1);    	
        	h = new H1F("TL2_"+k+"-"+is+"-"+run,"TL2_"+k+"-"+is+"-"+run,100,t3,t4);
        	h.setTitleX("Sector "+is+" PCAL U3 Cluster Time-fgo"); h.setTitleY("Counts"); h.setOptStat("1100");     
            dg.addDataSet(h,is+5);
            h2 = new H2F("TL3_"+k+"-"+is+"-"+run,"TL3_"+k+"-"+is+"-"+run,100,-2,2,62,1,63);
            h2.setTitleX("Sector "+is+" F1B Vertex Time"); h2.setTitleY("Paddle"); 
            dg.addDataSet(h2,is+11);
        }
        this.getDataGroup().add(dg,0,0,k,run);  
        createTDC1DHistos(k,-10.,10.,"T-TVERT-PATH/#beta*c (ns)"); 
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
                h = new H2F("uvw-pcal-u"+ip+"-s"+is+"-"+k+"-"+run,"uvw-pcal-u"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"U"+ip);       
                dg1.addDataSet(h,ip-1); dg1.addDataSet(f1,ip-1);
                uvw = (ip>15)?(30+(ip-15)):2*ip;
                ymx=(scaly)?uvw*4.5*1.23:ymax;
                if(scaly2) {ymx=0.82*ymax*(0.75+0.25*ip/npmt[1]) ;}
                h = new H2F("uvw-pcal-v"+ip+"-s"+is+"-"+k+"-"+run,"uvw-pcal-v"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"V"+ip);
                dg2.addDataSet(h,ip-1); dg2.addDataSet(f1,ip-1);
                h = new H2F("uvw-pcal-w"+ip+"-s"+is+"-"+k+"-"+run,"uvw-pcal-w"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"W"+ip); 
                dg3.addDataSet(h,ip-1); dg3.addDataSet(f1,ip-1);
     	    }
            this.getDataGroup().add(dg1,is,1,k,run); this.getDataGroup().add(dg2,is,2,k,run); this.getDataGroup().add(dg3,is,3,k,run);
            
            DataGroup dg4 = new DataGroup(6,6); DataGroup dg5 = new DataGroup(6,6); DataGroup dg6 = new DataGroup(6,6);        	         	   
            f1 = new F1D("p0"+is+2+k,"[a]",xmin*sca1-xoff,xmax*sca1-xoff); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
     	    for (int ip=1; ip<npmts[1]+1; ip++) {double ymx=(scaly)?ymax*ip/npmts[1]:ymax;int ybns=(scaly)?((ip>4)?ybins*ip/npmts[1]:5):ybins;
                if(scaly2) {ymx=0.75*ymax*(1-0.47*ip/npmts[1]);}
                h = new H2F("uvw-ecin-u"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecin-u"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin*sca1-xoff,xmax*sca1-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECIN "+xtxt);  h.setTitleY(ytxt+"U"+ip); 
                dg4.addDataSet(h,ip-1); dg4.addDataSet(f1,ip-1);
                if(scaly2) {ymx=0.55*ymax*(0.63+0.37*ip/npmts[1]) ;}
                h = new H2F("uvw-ecin-v"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecin-v"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin*sca1-xoff,xmax*sca1-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"V"+ip); 
                dg5.addDataSet(h,ip-1); dg5.addDataSet(f1,ip-1);
                if(scaly2) {ymx=0.75*ymax*(0.63+0.37*ip/npmts[1]) ;}
                h = new H2F("uvw-ecin-w"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecin-w"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin*sca1-xoff,xmax*sca1-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"W"+ip);
                dg6.addDataSet(h,ip-1); dg6.addDataSet(f1,ip-1);
                
     	    }
            this.getDataGroup().add(dg4,is,4,k,run); this.getDataGroup().add(dg5,is,5,k,run); this.getDataGroup().add(dg6,is,6,k,run);
     	   
            DataGroup dg7 = new DataGroup(6,6); DataGroup dg8 = new DataGroup(6,6); DataGroup dg9 = new DataGroup(6,6);        	         	   
            f1 = new F1D("p0"+is+3+k,"[a]",xmin*sca2-xoff,xmax*sca2-xoff); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
     	    for (int ip=1; ip<npmts[2]+1; ip++) {double ymx=(scaly)?ymax*ip/npmts[2]:ymax;int ybns=(scaly)?((ip>4)?ybins*ip/npmts[2]:5):ybins;
                if(scaly2) {ymx=0.75*ymax*(1-0.47*ip/npmts[2]);}
                h = new H2F("uvw-ecou-u"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecou-u"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin*sca2-xoff,xmax*sca2-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECOU "+xtxt); h.setTitleY(ytxt+"U"+ip);
                dg7.addDataSet(h,ip-1); dg7.addDataSet(f1,ip-1);
                if(scaly2) {ymx=0.5*ymax*(0.63+0.37*ip/npmts[2]) ;}
                h = new H2F("uvw-ecou-v"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecou-v"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin*sca2-xoff,xmax*sca2-xoff,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" ECOU "+xtxt);  h.setTitleY(ytxt+"V"+ip);
                dg8.addDataSet(h,ip-1); dg8.addDataSet(f1,ip-1);
                if(scaly2) {ymx=0.75*ymax*(0.63+0.37*ip/npmts[1]) ;}
                h = new H2F("uvw-ecou-w"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecou-w"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin*sca2-xoff,xmax*sca2-xoff,ybns,ymin,ymx);
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
            h = new H2F("tdc-pcal-u-"+is+"-"+k+"-"+run,"tdc-pcal-u-"+is+"-"+k+"-"+run,100, tmin, tmax, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("tdc-pcal-v-"+is+"-"+k+"-"+run,"tdc-pcal-v-"+is+"-"+k+"-"+run,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("tdc-pcal-w-"+is+"-"+k+"-"+run,"tdc-pcal-w-"+is+"-"+k+"-"+run,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
            
            h = new H2F("tdc-ecin-u_"+is+"-"+k+"-"+run,"tdc-ecin-u-"+is+"-"+k+"-"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("tdc-ecin-v-"+is+"-"+k+"-"+run,"tdc-ecin-v-"+is+"-"+k+"-"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("tdc-ecin-w-"+is+"-"+k+"-"+run,"tdc-ecin-w-"+is+"-"+k+"-"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
            
            h = new H2F("tdc-ecou-u-"+is+"-"+k+"-"+run,"tdc-ecou]-u-"+is+"-"+k+"-"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("tdc_ecou_v_"+is+"-"+k+"-"+run,"tdc-ecou-v-"+is+"-"+k+"-"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("tdc-ecou-w-"+is+"-"+k+"-"+run,"tdc-ecou-w-"+is+"-"+k+"-"+run,100, tmin, tmax, 36, 1., 37.);
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
            h = new H1F("tdc1-pcal-u-"+is+"-"+k+"-"+run,"tdc1-pcal-u-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" PCAL U "+txt);
            dg.addDataSet(h,0);  
            h = new H1F("tdc1-pcal-v-"+is+"-"+k+"-"+run,"tdc1-pcal-v-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" PCAL V "+txt);       
            dg.addDataSet(h,1);            
            h = new H1F("tdc1-pcal-w-"+is+"-"+k+"-"+run,"tdc1-pcal-w-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" PCAL W "+txt);  
            dg.addDataSet(h,2); 
            
            h = new H1F("tdc1-ecin-u-"+is+"-"+k+"-"+run,"tdc1-ecin-u-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECIN U "+txt);  
            dg.addDataSet(h,3);  
            h = new H1F("tdc1-ecin-v-"+is+"-"+k+"-"+run,"tdc1-ecin-v-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECIN V "+txt);        
            dg.addDataSet(h,4);            
            h = new H1F("tdc1-ecin-w-"+is+"-"+k+"-"+run,"tdc1-ecin-w-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECIN W "+txt);   
            dg.addDataSet(h,5); 
            
            h = new H1F("tdc1-ecou-u-"+is+"-"+k+"-"+run,"tdc1-ecou-u-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECOU U "+txt);   
            dg.addDataSet(h,6);  
            h = new H1F("tdc1-ecou-v-"+is+"-"+k+"-"+run,"tdc1-ecou-v-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECOU V "+txt);      
            dg.addDataSet(h,7);            
            h = new H1F("tdc1-ecou-w-"+is+"-"+k+"-"+run,"tdc1-ecou-w-"+is+"-"+k+"-"+run,100, tmin, tmax);
            h.setTitleX("Sector "+is+" ECOU W "+txt); 
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,0,k,run);
        }            

    } 
    
    public void initCCDB(int runno) {
    	System.out.println(getDetectorName()+".initCCDB("+runno+")");
        gain    = cm.getConstants(runno, "/calibration/ec/gain");
        time    = cm.getConstants(runno, "/calibration/ec/timing");
        oftime  = cm.getConstants(runno, "/calibration/ec/ftiming");
        ftime   = cm.getConstants(runno, "/calibration/ec/ftime");
        dtime   = cm.getConstants(runno, "/calibration/ec/dtime");
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
       
       BGOFFSET = isMC ? 102:0; 
       
       int elec_trigger_sect = isMC ? 5 : getElecTriggerSector();       
       int htcc_trigger_sect = getHTCCTriggerSector()-1;
       
       boolean goodECALSector = elec_trigger_sect>0 && elec_trigger_sect<7; 
       boolean goodHTCCSector = htcc_trigger_sect>0 && htcc_trigger_sect<7; 
       
       if(goodHTCCSector) trigger_sect = htcc_trigger_sect;
       if(goodECALSector) trigger_sect = elec_trigger_sect;
       
       boolean goodSector = goodECALSector || goodHTCCSector;
       
       if(!goodSector) return;
       
       phase  = getTriggerPhase(); 
       
       STT = event.hasBank("REC::Event") ? event.getBank("REC::Event").getFloat("startTime", 0):0;
       RF  = event.hasBank("REC::Event") ? event.getBank("REC::Event").getFloat("RFTime",    0):0;
   		 
       isGoodTL = event.hasBank("REC::Event")&&
    		      event.hasBank("REC::Particle")&&
    		      event.hasBank("REC::Calorimeter")&&
    		      STT>0;
    		     	   
       if(isGoodTL)     processTL(event); // TimeLine analysis for TL tab
       if(!dropSummary) processRaw(event);
       
       processRec(event);       
       
    }
    
    public void processTL(DataEvent event) {
    	
  	   int  run = getRunNumber();
  	   
  	   ((H1F) this.getDataGroup().getItem(0,0,tlnum,run).getData(trigger_sect-1).get(0)).fill(STT-FTOFFSET);
  	   
  	   if(event.hasBank("ECAL::clusters")){       
  		   DataBank  bank = event.getBank("ECAL::clusters");
  		   for(int loop = 0; loop < bank.rows(); loop++){
  			   int is = bank.getByte("sector", loop);
  			   int il = bank.getByte("layer", loop);
  			   if (is==trigger_sect && il==1){
  				   float    t = bank.getFloat("time",loop);
  				   int iU = (bank.getInt("coordU", loop)-4)/8+1;
  				   if(iU==3) ((H1F) this.getDataGroup().getItem(0,0,tlnum,run).getData(is+5).get(0)).fill(t+TOFFSET-FTOFFSET);
  			   }
  		   }
       }
  	   
  	   if(event.hasBank("REC::Scintillator")) {
  		   DataBank  bank = event.getBank("REC::Scintillator");
  		   for(int loop = 0; loop < bank.rows(); loop++){  		   
  			   int is = bank.getByte("sector", loop);
			   int il = bank.getByte("layer", loop);
			   int id = bank.getByte("detector", loop);
	           int     ip = bank.getShort("component",loop);
 			   if(is==trigger_sect && il==2 && id==12) {
 				   float  nrg = bank.getFloat("energy", loop);
  	               float path = bank.getFloat("path", loop);
  	               float time = bank.getFloat("time", loop);
	               ((H2F) this.getDataGroup().getItem(0,0,tlnum,run).getData(is+11).get(0)).fill(time-STT-path/29.97,ip);
  			   }
  		   }
  	   }
      
    }
    
    public void processRaw(DataEvent event) { //duplicates ECEngine process flow to check consistency
    	
 	   int  run = getRunNumber();
 	  
       tdcs.clear();
       
       if(event.hasBank("ECAL::tdc")==true){
           DataBank  bank = event.getBank("ECAL::tdc");
           for(int i = 0; i < bank.rows(); i++){
               int  is = bank.getByte("sector",i);
               int  il = bank.getByte("layer",i);
               int  ip = bank.getShort("component",i);   
               float tdcd  = bank.getInt("TDC",i)*tps-BGOFFSET; // raw time
               float tdcdc = tdcd-phase; // phase corrected time
               if(is>0 && is<7 && tdcd>0) {
                   if(!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Float>(),is,il,ip);    
              	       ((H2F) this.getDataGroup().getItem(is,0,0,run).getData(il-1).get(0)).fill(tdcd-FTOFFSET, ip); 
              	       ((H2F) this.getDataGroup().getItem(is,0,1,run).getData(il-1).get(0)).fill(tdcdc-FTOFFSET,ip); 
                       if (is==trigger_sect || isMC) {
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
                                                 + FTOFFSET-BGOFFSET;
               
               float tmax = 1000; float tdcm = 1000;
              
               if (tdcs.hasItem(is,il,ip)) { // sector,layer,component FADC/TDC match
            	   float radc = (float)Math.sqrt(adc);
                   for (float tdc : tdcs.getItem(is,il,ip)) {
                	   float tdif = tdc - (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - t;
                	   ((H2F) this.getDataGroup().getItem(is,0,8,run).getData(il-1).get(0)).fill(tdif,ip); // FADC t - TDC time
                	   if (Math.abs(tdif)<10 && tdif<tmax) {tmax = tdif; tdcm = tdc;}                	    
                   }
                   double a0 = time.getDoubleValue("a0", is, il, ip); 
                   double a2 = time.getDoubleValue("a2", is, il, ip);
                   double a3 = time.getDoubleValue("a3", is, il, ip);
                   double a4 = time.getDoubleValue("a4", is, il, ip);
                   
                   double tdcmc = tdcm - a0 -  (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - a2/radc - a3 - a4/Math.sqrt(radc);
          	       ((H2F) this.getDataGroup().getItem(is,0,3,run).getData(il-1).get(0)).fill(tdcm-FTOFFSET,ip);  // matched FADC/TDC
          	       ((H2F) this.getDataGroup().getItem(is,0,4,run).getData(il-1).get(0)).fill(tdcmc-FTOFFSET,ip); // calibrated time
          	       
          	       ((H2F) this.getDataGroup().getItem(is,0,5,run).getData(il-1).get(0)).fill(tdcmc-t+a0,ip); // calibrated TDC-FADC time
          	       if(isGoodTL && il==1 && ip==3) ((H1F) this.getDataGroup().getItem(0,0,tlnum,run).getData(is+5).get(0)).fill(tdcmc); 
               }
           }
       }    	
    }
    
    public void processRec(DataEvent event) { // process reconstructed timing
    	
       RF = event.hasBank("REC::Event") ? event.getBank("REC::Event").getFloat("RFTime",0):0;

       DataBank recRunRF  = null; float trf = 0;

       if(event.hasBank("RUN::rf"))  {
    	   recRunRF = event.getBank("RUN::rf");             
           for(int k = 0; k < recRunRF.rows(); k++){
        	   if(recRunRF.getInt("id", k)==1) trf = recRunRF.getFloat("time",k);
           }    
       }
       
       if(isMC && STT<0) STT=124.25f;
       
       if(!(STT>0)) return;
       
       int run = getRunNumber(); 
        
       if(dropBanks) dropBanks(event); // rerun ECEngine with updated CCDB constants
       
       List<ECStrip>     strips = eng.engine.getStrips();
       List<ECPeak>       peaks = eng.engine.getPeaks(); 
       List<ECCluster> clusters = eng.engine.getClusters();
       
       if (run>=2000) {phase=0; shiftTV[1]=0;} // Corrections needed until runs<4013 are recooked
       
       if(!dropSummary) {
/*    	   
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
*/       
       }
       
       if(!(event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("ECAL::clusters"))) return;
              
       DataBank bankc = event.getBank("REC::Calorimeter"); //entries ordered by pindex,layer
       DataBank bankp = event.getBank("REC::Particle");
       
//   	   HashMap<Integer,ArrayList<Integer>>  mapPndex = mapByIndex(bankc,"pindex"); //maps part bank index to REC::Calorimeter
//   	   HashMap<Integer,ArrayList<Integer>>  mapCndex = mapByIndex(bankc,"index");  //maps clus bank index to REC::Calorimeter 
       
/*       if(event.hasBank("ECAL::calibpass2")) {
    	   DataBank bank = event.getBank("ECAL::calibpass2");
    	   for (int loop =0; loop<bank.rows(); loop++) {
    		   int is = bank.getByte("sector",loop);
    		   int il = bank.getByte("layer", loop);
    		   float e2 = bank.getFloat("energy",loop);
    		   float t2 = bank.getFloat("time",loop);
    		   float ft2 = bank.getFloat("ftime",loop);
    		   short db1 = bank.getShort("dbstU",loop);
    		   short db2 = bank.getShort("dbstV",loop);
    		   short db3 = bank.getShort("dbstW",loop);
    		   float dt2u = bank.getFloat("recDTU",loop);
    		   float dt2v = bank.getFloat("recDTV",loop);
    		   float dt2w = bank.getFloat("recDTW",loop);
    		   float ft2u = bank.getFloat("recFTU",loop);
    		   float ft2v = bank.getFloat("recFTV",loop);
    		   float ft2w = bank.getFloat("recFTW",loop);
    		   System.out.println(is+" "+il+" "+e2+" "+t2+" "+ft2+" "+db1+" "+db2+" "+db3+" "+dt2u+" "+dt2v+" "+dt2w+" "+ft2u+" "+ft2v+" "+ft2w);
    	   }    
       }
*/
       pathlist.clear();
       tdcs.clear();
       
       for(int loop = 0; loop < bankc.rows(); loop++){ //loop over REC::Calorimeter
          int     is = bankc.getByte("sector",loop);
          int     il = bankc.getByte("layer",loop);
          int     in = bankc.getShort("index",loop); // index to cooked ECAL::clusters (now replaced by dropBanks re-cooking)
          int    det = bankc.getByte("detector",loop);
          if (det==7 && !pathlist.hasItem(is,il,in)) pathlist.add(loop,is,il,in); // associate ECAL::cluster index to REC::Calorimeter index               
       }
       
       int[] iip = new int[3], iid = new int[3], statp = new int[3];
       float[] tid = new float[3], lef = new float[3], add = new float[3];
             
       DataBank  bank1 = event.getBank("ECAL::clusters"); //entries ordered by sector,layer
       DataBank  bank3 = event.getBank("ECAL::peaks");
              
       for(int loop = 0; loop < bank1.rows(); loop++){ //loop over ECAL::clusters
           int is = bank1.getByte("sector", loop);
             if (true) {
               int     il =  bank1.getByte("layer", loop);
               float ener =  bank1.getFloat("energy",loop)*1000;
               float    t =  bank1.getFloat("time",loop);

               Point3D  pc = new Point3D(bank1.getFloat("x",loop),
            		                     bank1.getFloat("y",loop),
            		                     bank1.getFloat("z",loop));
               
               iip[0] = (bank1.getInt("coordU", loop)-4)/8+1; //strip number of peak U
               iip[1] = (bank1.getInt("coordV", loop)-4)/8+1; //strip number of peak V
               iip[2] = (bank1.getInt("coordW", loop)-4)/8+1; //strip number of peak W
               iid[0] =  bank1.getInt("idU",loop)-1; //peak index U
               iid[1] =  bank1.getInt("idV",loop)-1; //peak index V
               iid[2] =  bank1.getInt("idW",loop)-1; //peak index W	               
               lef[0] = getLeff(pc,getPeakline(iid[0],pc,bank3)); //readout distance U
               lef[1] = getLeff(pc,getPeakline(iid[1],pc,bank3)); //readout distance V
               lef[2] = getLeff(pc,getPeakline(iid[2],pc,bank3)); //readout distance W
               tid[0] = bank3.getFloat("time",iid[0]) - lef[0]/(float)veff.getDoubleValue("veff", is,il+0,iip[0]); //readout time subtracted U
               tid[1] = bank3.getFloat("time",iid[1]) - lef[1]/(float)veff.getDoubleValue("veff", is,il+1,iip[1]); //readout time subtracted V
               tid[2] = bank3.getFloat("time",iid[2]) - lef[2]/(float)veff.getDoubleValue("veff", is,il+2,iip[2]); //readout time subtracted W
               add[0] = bank3.getFloat("energy",iid[0]); //peak energy U
               add[1] = bank3.getFloat("energy",iid[1]); //peak energy V
               add[2] = bank3.getFloat("energy",iid[2]); //peak energy W
               
               int statc = bank1.getShort("status", loop);
               statp[0]  = bank3.getShort("status", iid[0]);
               statp[1]  = bank3.getShort("status", iid[1]);
               statp[2]  = bank3.getShort("status", iid[2]);
               
               if (pathlist.hasItem(is,il,loop)) {
                   int    pin = bankc.getShort("pindex", pathlist.getItem(is,il,loop));
                   float path = bankc.getFloat("path",   pathlist.getItem(is,il,loop)); 
                   int    pid = bankp.getInt("pid",pin);
                   float beta = bankp.getFloat("beta",pin); 
                   
//                   Point3D  vc = new Point3D(bankp.getFloat("vx",pin),
//                                             bankp.getFloat("vy",pin),
//                                             bankp.getFloat("vz",pin));
 
//                   if(pid==22 || pid==2112) path = (float) pc.distance(vc);
                   
                   int status = Math.abs(bankp.getInt("status",pin));
                   
                   for (int i=0; i<3; i++) { //loop over U,V,W
                	   float tu=0,tdc=0,tdcc=0,tdccc=0,leff=0,adc=0; int ip=0;
                	   if (clusters.size()>0) { // use ECEngine structure clusters
                         tu    = (float) clusters.get(loop).getTime(i);
                         ip    =         clusters.get(loop).getPeak(i).getMaxStrip();
                         adc   =         clusters.get(loop).getPeak(i).getMaxECStrip().getADC();
                         tdc   = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getRawTime(true)-TOFFSET;
                         tdcc  = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTWCTime();
//                         tdcc = tdc-tw[il+i-1]/(float)Math.sqrt(adc);
//                         tdccc = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTime(); 
                         leff  = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTdist();
//                         System.out.println("C Layer "+il+" "+clusters.get(loop).getPeak(i).getMaxECStrip().getLine().toString());
//                         System.out.println("C Layer "+il+" "+clusters.get(loop).getPeak(i).getLine().toString());
//                         System.out.println("P Layer "+il+" "+peaks.get(iid[i]).getLine().toString());
//                         System.out.println("P Layer "+il+" "+getPeakline(iid[i],pc,bank3).toString());
//                         System.out.println("Layer "+il+" "+tu+" "+tid[i]+" "+leff+ " "+lef[i]);
                	   } else { // use ECAL::clusters bank
                		 tu    = tid[i];
                		 ip    = iip[i];
                		 leff  = lef[i];
                         adc   = 10000*add[i]/(float)gain.getDoubleValue("gain", is, il+i, ip);
                	   }
                	   
                	   float vcorr = STT - phase + TVOffset;  // phase=0 unless STT is not phase corrected (early engineering runs)
                       float tvcor = tu  - vcorr;
                       float mybet = path/tvcor/c; // use ECAL beta for beta distribution plots and neutral residuals
                       
                       if(Math.abs(pid)==211)((H1F) this.getDataGroup().getItem(0,0,22,run).getData(getDet(il)+3).get(is-1)).fill(mybet);  
                       if(pid==22||pid==2112)((H1F) this.getDataGroup().getItem(0,0,22,run).getData(getDet(il)+6).get(is-1)).fill(mybet); 
                       if(pid==11)           ((H1F) this.getDataGroup().getItem(0,0,22,run).getData(getDet(il)  ).get(is-1)).fill(mybet); 
                      
                       float vel=c; 
                       if(Math.abs(pid)==211 || Math.abs(pid)==2212) vel=Math.abs(beta*c); //use EB beta for pion or proton calibration residuals                       
                       
                	   float pcorr = path/vel;
                       float resid = tvcor - pcorr;  
                       float residt = t-vcorr-pcorr;
                       
                       if(!tdcs.hasItem(is,il+i,ip)) tdcs.add(new ArrayList<Float>(),is,il+i,ip);
                           tdcs.getItem(is,il+i,ip).add(resid);
                           
                	   float  radc = (float) Math.sqrt(adc);                	   
                       float lcorr = leff/(float)veff.getDoubleValue("veff", is, il+i, ip);
                       
                       boolean goodSector = dropEsect ? is!=trigger_sect : is==trigger_sect;  
                       
                       boolean goodPID = TRpid==11211 ? Math.abs(pid)==211 && is!=trigger_sect || 
                    		                            Math.abs(pid)==11  && is==trigger_sect : //combine e- and pi+/pi-
                    		                            Math.abs(pid)==TRpid && goodSector;      //PID from menu selection
                                              
                       boolean goodStatus = status>=2000 && status<3000; 
                       boolean goodHisto  = goodPID && goodStatus;
                                              
                       if (goodHisto) {
                    	   ((H1F) this.getDataGroup().getItem(is,0,tlnum,run).getData(il+i-1).get(0)).fill(resid); //timelines
                           ((H2F) this.getDataGroup().getItem(is,0,   10,run).getData(il+i-1).get(0)).fill(resid, ip);
                       }
                       
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
                   } 
               }
           }
       }
       
       IndexGenerator ig = new IndexGenerator();
       
//       System.out.println("");
       for (Map.Entry<Long,List<Float>>  entry : tdcs.getMap().entrySet()){
           long hash = entry.getKey();
           int is = ig.getIndex(hash, 0);
           int il = ig.getIndex(hash, 1);
           int ip = ig.getIndex(hash, 2);
//           System.out.println(is+" "+il+" "+ip+" "+tdcs.getItem(is,il,ip));
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
                                 double max = g.getDataY(g.getDataSize(0)-1); double mn=Math.min(min, max); 
                                 double  mx = Math.max(min, max);                                                
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
        if (eng.useFADCTime &&  eng.usePass2Time) {writeFile("ftime",1,7,0,3,0,3);  writeFile("fveff",1,7,0,3,0,3);}
        if(!eng.useFADCTime &&  eng.usePass2Time) {writeFile("dtime",1,7,0,3,0,3);  writeFile("dveff",1,7,0,3,0,3);}
        if(!eng.useFADCTime && !eng.usePass2Time) {writeFile("timing",1,7,0,3,0,3); writeFile("effective_velocity",1,7,0,3,0,3);}    
        isAnalyzeDone = true;
    }
    
    public void analyzeResiduals() {
        getResidualSummary(1,7,0,3,0,3);
        System.out.println("analyzeResiduals Finished");
        if (eng.useFADCTime &&  eng.usePass2Time) writeFile("ftime_update",1,7,0,3,0,3);
        if(!eng.useFADCTime &&  eng.usePass2Time) writeFile("dtime_update",1,7,0,3,0,3);
        if(!eng.useFADCTime && !eng.usePass2Time) writeFile("timing_update",1,7,0,3,0,3);
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
                   g = graphShift(fitter2.getMeanSlices(),-tl.fitData.getItem(is,3*il+iv+1,ip,run).p0); //subtract t0 just fitted
                   g.getAttributes().setTitleX("Sector "+is+" "+det[il]+" "+v[iv]+(ip+1)); g.getAttributes().setTitleY(""); 
                   if(dumpGraphs) {
                	   String nam =eng.useFADCTime?"ftime/":"dtime/"+"s"+is+"l"+il+"v"+iv+"p"+ip;
                	   System.out.println("Writing "+filPath+"twplots/"+nam+".vec"); dumpGraph(filPath+"twplots/"+nam+".vec",g);
                   }
//            	   if(!eng.useFADCTime && il>2) tl.fitData.add(fitEngine(g,16,20),is,3*il+iv+1,ip+200,run); //TW fits to DISC time 
                   boolean doFit = true;
            	   if(is==6&&il==1&&iv==2&&ip==35) {System.out.println("Rejecting Sector "+is+" Layer "+il+" View "+iv+" PMT "+ip);doFit=false;}
            	   if(!eng.useFADCTime) tl.fitData.add(DTWFit(g,is,il,iv,doFit),is,3*il+iv+1,ip+200,run); //TW fits to DISC time           	   
            	   if( eng.useFADCTime) tl.fitData.add(FTWFit(g,is,il,iv),is,3*il+iv+1,ip+200,run);       //TW fits to FADC time            	   
                }
             }
          }
       } 
    }
    
    public FitData FTWFit(GraphErrors g, int is, int il, int iv) {
    	int fnum=13, fmin=10;
    	FitData fd = new FitData(g);
    	if(g.getDataSize(0)==0) return fd;
        fd.initFit(fnum,0,1,fmin,g.getDataX(g.getDataSize(0)-1)*1.05); 
        switch (il) {
        case 0: fd.initFunc(0,-0.5); fd.initFunc(1,20,18,22);  fd.initFunc(2,9,7,11);   fd.initFunc(3,170,160,180); fd.initFunc(4,20,15,25); break;
        case 1: fd.initFunc(0,-0.5); fd.initFunc(1,13.9,5,22); fd.initFunc(2,5.4,5,11); fd.initFunc(3,125,100,180); fd.initFunc(4,15,14,25); break;
        case 2: fd.initFunc(0,-0.5); fd.initFunc(1,17,10,18);  fd.initFunc(2,5,3,11);   fd.initFunc(3,145,110,180); fd.initFunc(4,15,14,20);
        }
        fd.doFit = true; 
        fd.fitGraph("",true,true); 
        return fd;
    }
   
    public FitData DTWFit(GraphErrors g, int is, int il, int iv, boolean doFit) {
    	int fitnum = 17;   
    	int fmin = il==0 ? 23:is==5 ? 12:15; if(il==2 && iv>0) fmin=10;
    	FitData fd = new FitData(g);
    	if(g.getDataSize(0)==0) return fd;
        fd.initFit(fitnum,0,1,fmin,g.getDataX(g.getDataSize(0)-1)*1.05); 
        switch (il) {
        case 0: fd.initFunc(0,-0.5);       fd.initFunc(1,25,25,27);  fd.initFunc(2,5,3,8);     fd.initFunc(3,160,140,180); fd.initFunc(4,30,15,35); 
                fd.initFunc(5, 15,10,18);  fd.initFunc(6,0.3,0,2); break;
        case 1: fd.initFunc(0,-0.5,-1,1);  fd.initFunc(1,13.9,5,22); fd.initFunc(2,5.4,2,11);  fd.initFunc(3,155,80,220);  fd.initFunc(4,15,20,145); 
                fd.initFunc(5, 0.5,0.3,1); fd.initFunc(6,0.7,0,1); break;
        case 2: fd.initFunc(0,-0.5);       fd.initFunc(1,20,19,21);  fd.initFunc(2,4,3.2,4.8); fd.initFunc(3,260,120,261); fd.initFunc(4,34,30,40); 
                fd.initFunc(5, 6,5.9,6.1); fd.initFunc(6,0.3,-1,1); 
        }
        if(il==2 && iv>0) {fd.initFunc(0,-0.5,-1,1);  fd.initFunc(1,13.9,5,22); fd.initFunc(2,5.4,2,11); fd.initFunc(3,155,80,220); fd.initFunc(4,15,20,145);         
                           fd.initFunc(5, 0.5,0.3,1); fd.initFunc(6,0.7,0,1);}
        fd.doFit = doFit; 
        fd.fitGraph("",true,true); 
        return fd;
    } 
    
    public void getResidualSummary(int is1, int is2, int il1, int il2, int iv1, int iv2) {
    	
        ParallelSliceFitter fitter;
        GraphErrors g1,g2;
        
        int run=getRunNumber();
        
        for (int is=is1; is<is2; is++) {            
           for (int il=il1; il<il2; il++) {
              for (int iv=iv1; iv<iv2; iv++) {
                  fitter = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(is,0,10,run).getData(3*il+iv).get(0));
                  fitter.setRange(-30,30); fitter.fitSlicesY(); 
                  g1 = sliceToGraph(fitter.getMeanSlices(),il,iv); 
                  g2 = sliceToGraph(fitter.getSigmaSlices(),il,iv);
                  GraphErrors mean = new GraphErrors("RESIDUAL_"+is+"-"+il+" "+iv,g1.getVectorX().getArray(),
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
                  GraphErrors mean = new GraphErrors("TMF-"+is+"-"+il+" "+iv,g1.getVectorX().getArray(),
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
			
			System.out.println(getDetectorName()+".writefile("+table+")");

			for (int is=is1; is<is2; is++) {
				for (int il=il1; il<il2; il++ ) {
					for (int iv=iv1; iv<iv2; iv++) {
						for (int ip=0; ip<npmt[3*il+iv]; ip++) {
							switch (table) {
							case "timing_update":      line =  getA0( is,il,iv,ip,il==il1 ? 0 : A0offset); break; //RGM pass0 EC Z-plane compensation
							case "dtime_update":       line =  getDA0(is,il,iv,ip,il==il1 ? 0 : A0offset); break; //RGM pass0 EC Z-plane compensation
							case "ftime_update":       line =  getFA0(is,il,iv,ip,il==il1 ? 0 : A0offset); break; //RGM pass0 EC Z-plane compensation
							case "timing":             line =  getTW(is,il,iv,ip);  break;
							case "dtime":              line =  getDTW(is,il,iv,ip,false); break; //true: update TWC but not offset
							case "ftime":              line =  getFTW(is,il,iv,ip); break;
							case "effective_velocity": line =  getEV(is,il,iv,ip);  break;
							case "fveff":              line =  getEV(is,il,iv,ip);  break;
							case "tmf_offset":         line =  getTMF(is,il,iv,ip); break;  
							case "fadc_offset":        line =  getGTMF(is,il,iv,ip);  
							case "setFTiming":         line =  getFTiming(is,il,iv,ip); break;
							case "A0offset":           line =  getNewA0(is,il,iv,ip); break;
							case "NewFTW":             line =  getNewFTW(is,il,iv,ip); break;
							case "NewDTW":             line =  getNewDTW(is,il,iv,ip); break;
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
	
	public String getTW(int is, int il, int iv, int ip) { //Absolute TDC residual calibration + pass1 time walk
		if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p0 //fits to corrected t vs leff
				+" 0.02345 "
				+(time.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p1)+" "
				+(time.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p0)+" "
				+(time.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p2)+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0";
		}
		
	}
	
	public String getDTW(int is, int il, int iv, int ip, boolean update) { //Absolute TDC residual calibration + pass2 time walk
		if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "+
		    	(update ? dtime.getDoubleValue("a0",is,3*il+iv+1,ip+1):tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p0)				
				+" 0.02345 "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p0+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p1+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p2+" "
			    +tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p3+" "
			    +tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p4+" "
			    +tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p5+" "
			    +tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p6+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0 0.0 0.0 0.0 0.0";
		}
		
	}
		
	public String getFTW(int is, int il, int iv, int ip) { //Absolute FADC residual calibration + pass 2 time walk
		if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p0 //fits to corrected t vs leff
				+" 0.02345 "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p0+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p1+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p2+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p3+" "
				+tl.fitData.getItem(is,3*il+iv+1,ip+200,getRunNumber()).p4+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0 0.0 0.0";
		}
		
	}

	public String getA0(int is, int il, int iv, int ip, float off) { //Relative residual calibration
		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(time.getDoubleValue("a0", is, 3*il+iv+1, ip+1) 
				+FitSummary.getItem(is,il,iv,getRunNumber()).getDataY(ip)+off)+" "  //residual correction to previous a0
				+" 0.02345 "
				+time.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+" " //replace with table
				+time.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+" " //replace with table
				+time.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+" ";//replace with table
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0";
		}
		
	}
	
	public String getDA0(int is, int il, int iv, int ip, float off) { //Relative residual calibration
		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(dtime.getDoubleValue("a0", is, 3*il+iv+1, ip+1) 
				+FitSummary.getItem(is,il,iv,getRunNumber()).getDataY(ip)+off)+" "  //residual correction to previous a0
				+" 0.02345 "
				+dtime.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+" "
				+dtime.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+" "
				+dtime.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+" "
				+dtime.getDoubleValue("a5", is, 3*il+iv+1, ip+1)+" "
				+dtime.getDoubleValue("a6", is, 3*il+iv+1, ip+1)+" "
			    +dtime.getDoubleValue("a7", is, 3*il+iv+1, ip+1)+" "
			    +dtime.getDoubleValue("a8", is, 3*il+iv+1, ip+1)+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0 0.0 0.0";
		}
		
	}
	
	public String getFA0(int is, int il, int iv, int ip, float off) { //Relative residual calibration
		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(ftime.getDoubleValue("a0", is, 3*il+iv+1, ip+1) 
				+FitSummary.getItem(is,il,iv,getRunNumber()).getDataY(ip)+off)+" "  //residual correction to previous a0
				+" 0.02345 "
				+ftime.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+" "
				+ftime.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+" "
				+ftime.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+" "
				+ftime.getDoubleValue("a5", is, 3*il+iv+1, ip+1)+" "
				+ftime.getDoubleValue("a6", is, 3*il+iv+1, ip+1)+" ";
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0 0.0 0.0";
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
	
	public String getNewFTW(int is, int il, int iv, int ip) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(oftime.getDoubleValue("a0", is, 3*il+iv+1, ip+1)+" "
				+" 0.02345 "
				+oftime.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+" "
				+oftime.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+" "
				+oftime.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+" "
				+"0.0"+" "
				+"0.0");		
	}
	
	public String getNewDTW(int is, int il, int iv, int ip) {
	    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
			+(oftime.getDoubleValue("a0", is, 3*il+iv+1, ip+1)+" "
			+" 0.02345 "
			+"0.0"+" "
			+"0.0"+" "
			+"0.0"+" "
			+"0.0"+" "
			+"0.0"+" "
			+"0.0"+" "
			+"0.0");		
    }	
	
	public String getFTiming(int is, int il, int iv, int ip) {
		return is+" "+(3*il+iv+1)+" "+(ip+1)+" 0.0 "+" 0.02345 "+" 0.0 "+" 0.0 "+" 0.0";
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
    	System.out.println(getDetectorName()+": Saving timelines");
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
    	
//    	ECt ect = new ECt("ECt",6684);
//    	ect.writeFile("NewFTW",1,7,0,3,0,3);    
//    	ECt ect = new ECt("ECt",2385);
//    	ect.writeFile("NewDTW",1,7,0,3,0,3);   	
    }
     
}
