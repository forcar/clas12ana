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

import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECt extends DetectorMonitor {

    IndexedList<List<Integer>>  tdcs = new IndexedList<List<Integer>>(3);
    IndexedList<List<Float>>  tresid = new IndexedList<List<Float>>(3);
    IndexedTable time=null, oftime=null, ftime=null, dtime=null, fveff=null, dveff=null, gtw=null;
    IndexedTable   fo=null, fgo=null, tmf=null, tgo=null, gain=null, veff=null, rfT=null, tmfcut=null;
    
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
    Boolean                   isTRESDone = false;
    Boolean           isTimeLineFitsDone = false;
    
    int trigger_sect = 0;
    int triggerPhase = 0;
    float        STT = 0;
    float         RF = 0;
    float    RF_TIME = 124.25f;
    
    static int tlnum;
    
    Boolean       isMC = false;    
    Boolean   isGoodTL = false;
    
//    static float TOFFSET = 600; 
   
    static float    BGOFFSET = 102;
    static float    FTOFFSET = 0;
    static float      TMFCUT = 0;
    static float     TOFFSET = 0;
    static float   TOFFSETMC = 180;
    static float   TOFFSETER = 305;
    static float    RFPERIOD = 4.008f;
    static float BEAM_BUCKET = 2.004f;
    static float           c = 29.98f;
    static float    A0offset = 0f;  //RGM pass0 only for EC Z-plane compensation
    static float    A0sector = 0;  
    static float         tps =  (float) 0.02345;
    float[] shiftTV = {0,40,0,0,0,0}; //Run 3050 t0 calibration  
    
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
                                 "T0",
                                 "STATUS",
                                 "TRES",
                                 "TRFITS");
        

        
        this.useCALUVWSECButtons(true);
        this.usePCCheckBox(true);
        this.useSliderPane(true);
        this.useECEnginePane(true);
        tlnum = getDetectorTabNames().indexOf("TL");
        this.init();
        this.localinit("rga_fall2018");
    }
    
    public ECt(String name, int runno) {
    	super(name);
    	initCCDB(runno);    	
    }
    
    @Override
    public void localinit(String variation) {
    	System.out.println(getDetectorName()+".localinit("+variation+")");
    	eng.engine.setGeomVariation(variation);
    	tl.setFitData(Fits);
    }  
    
    @Override
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
        setRunNumber(run);  runlist.add(run); isMC = run<100;       
        BGOFFSET = isMC ? 0:0; 
        
        this.setNumberOfEvents(0);  
        
        int t1=-70, t2=150, t3=50, t4=250;
        createTDCHistos(10,-30.,30.,"T-TVERT-PATH/#beta*c (ns)"); 
        createBETAHistos(22); createStatus(32);
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
        createTRESHistos(33,50,50,0,50000,-5,5,"FADC","RESID ");
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
        plotBetaHistos(22);
        plotStatusHistos(32);
        plotStatusHistos(33);
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
	    if(isTRESDone)     plotTRESSummary(34);
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
   
        for (int ic=0; ic<2; ic++) {
        for (int iv=0; iv<3; iv++) {
        DataGroup dg = new DataGroup(3,3);
        for (int id=1; id<4; id++) {
    	for (int il=1; il<4; il++) {
        for (int is=1; is<7 ; is++) {
           h = new H1F("beta-"+k+"-"+is+"-"+il+"-"+id+"-"+iv+"-"+ic+"-"+run,"beta-"+k+"-"+is+"-"+il+"-"+id+"-"+iv+"-"+ic+"-"+run,100,0.4,1.5);
           h.setTitleX(pid[id-1]+" "+det[il-1]+" "+v[iv]+" "+" beta "); h.setTitleY("Counts"); h.setLineColor(is==6?9:is);
           h.setOptStat("1000000");
           dg.addDataSet(h,id*3+il-4);
        }
    	}
        }
        this.getDataGroup().add(dg,ic,iv,k,run); 
        }
        }
                
    }    
    
    public void createStatus(int k) {
    	
        int run = getRunNumber(), n=0;
        H2F h;
    	DataGroup dg = new DataGroup(3,6);

    	for(int is=1; is<7; is++) {
    		for(int im=1; im<4; im++) {
    			int nx = im==1 ? 68:36;
    			h = new H2F("STATUS"+is+im, "STATUS"+is+im, nx, 1, nx+1, 3, 1, 4);
    			h.setTitleX("SECTOR "+is);h.setTitle(is==1 ? det[im-1]:" ");
    			dg.addDataSet(h, n);n++;
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
    
    public void createTRESHistos(int k, int xbins, int ybins, double xmin, double xmax, double ymin, double ymax, String xtxt, String ytxt) {
        H2F h ;  
        int run = getRunNumber(), n=0;
        float[] emx = {50000,30000,30000,30000,12000,25000,20000,5000,10000};
    	DataGroup dg = new DataGroup(9,6);
    	
        for (int is=1; is<7; is++) { 
        	for (int il=0; il<3; il++) {
            	for (int iv=0; iv<3; iv++) {
            		int i=3*il+iv;
            		h = new H2F("tres-"+k+"-"+det[il]+"-"+v[iv]+"s"+is+"-"+run,"tres-"+k+"-"+det[il]+"-"+v[iv]+"-s"+is+"-"+run,xbins,xmin,emx[i],ybins,ymin,ymax);
            		dg.addDataSet(h,n);n++;
            	}
        	}
        }
        this.getDataGroup().add(dg,0,0,k,run);
    }
    
    public void createUVWHistos(int k, int xbins, int ybins, double xmin, double xmax, double ymin, double ymax, String xtxt, String ytxt) {
    	
        H2F h;  F1D f1; 
        
        int uvw,ybns;
        double ymx,sca1=1.0,sca2=1.0; double xoff=0; boolean scaly=false; boolean scaly2=false;
        int run = getRunNumber();
        if (k==12) {sca1=0.7; sca2=0.6;}
        if (k==19||k==20) {scaly=true;}
        if (k==16||k==17) {scaly2=true;}
        
        for (int is=1; is<7; is++) {              	
            DataGroup dg1 = new DataGroup(9,8); DataGroup dg2 = new DataGroup(8,8); DataGroup dg3 = new DataGroup(8,8);        	    
            f1 = new F1D("p0"+is+1+k,"[a]",xmin,xmax); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);          
            for (int ip=1; ip<npmts[0]+1; ip++) {
            	uvw = (ip>52)?(52+(ip-52)*2):ip; ymx=scaly?2*uvw*4.5*0.51:ymax; ybns=scaly?((ip>6)?ybins*ip/npmts[0]:5):ybins;
                if(scaly2) {ymx=ymax*(1-0.4*ip/npmt[0]);}
                h = new H2F("uvw-pcal-u"+ip+"-s"+is+"-"+k+"-"+run,"uvw-pcal-u"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybns,ymin,ymx);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"U"+ip);       
                dg1.addDataSet(h,ip-1); dg1.addDataSet(f1,ip-1);
                uvw = (ip>15)?(30+(ip-15)):2*ip; ymx=scaly?uvw*4.5*1.23:ymax;
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
     	    for (int ip=1; ip<npmts[1]+1; ip++) {
     	    	ymx=scaly?ymax*ip/npmts[1]:ymax; ybns=scaly?((ip>4)?ybins*ip/npmts[1]:5):ybins;
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
     	    for (int ip=1; ip<npmts[2]+1; ip++) {
     	    	ymx=scaly?ymax*ip/npmts[2]:ymax; ybns=scaly?((ip>4)?ybins*ip/npmts[2]:5):ybins;
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
  
    public IndexedTable getVeff() {
      return (eng.useFADCTime ? fveff : (eng.usePass2Timing ? dveff:veff));
    }
    
    public void initCCDB(int runno) {
    	System.out.println(getDetectorName()+".initCCDB("+runno+")");
        gain    = cm.getConstants(runno, "/calibration/ec/gain");
        time    = cm.getConstants(runno, "/calibration/ec/timing");
        oftime  = cm.getConstants(runno, "/calibration/ec/ftiming");
        ftime   = cm.getConstants(runno, "/calibration/ec/ftime");
        dtime   = cm.getConstants(runno, "/calibration/ec/dtime");
        veff    = cm.getConstants(runno, "/calibration/ec/effective_velocity");
        fveff   = cm.getConstants(runno, "/calibration/ec/fveff");
        dveff   = cm.getConstants(runno, "/calibration/ec/dveff");
        fo      = cm.getConstants(runno, "/calibration/ec/fadc_offset");        // Crate/fiber FADC offsets 
        fgo     = cm.getConstants(runno, "/calibration/ec/fadc_global_offset"); // FADC capture window 
        tgo     = cm.getConstants(runno, "/calibration/ec/tdc_global_offset");  // TDC capture window
        gtw     = cm.getConstants(runno, "/calibration/ec/global_time_walk");   // Global time walk correction using raw ADC
        tmf     = cm.getConstants(runno, "/calibration/ec/tmf_offset");         // TDC-FADC offsets
        tmfcut  = cm.getConstants(runno, "/calibration/ec/tmf_window");         // TDC-FADC cut
        rfT     = ebcm.getConstants(runno, "/calibration/eb/rf/config");  
       
        FTOFFSET = (float)    fgo.getDoubleValue("global_offset",0,0,0);
        TOFFSET  = (float)    tgo.getDoubleValue("offset", 0,0,0);
        RFPERIOD = (float)    rfT.getDoubleValue("clock",1,1,1); 
        TMFCUT   = (float) tmfcut.getDoubleValue("window", 0,0,0); //acceptance window for TDC-FADC cut
    }
    
    @Override
    public void processEvent(DataEvent event) {  

       if( dropBanks) dropBanks(event); // rerun ECEngine with updated CCDB constants
       
       if(!eventFilter(event))     return;
       if(!triggerFilter())        return;   	
       if(!pidFilter(event,false)) return;
       
       if(!dropSummary) processRaw(event);

       STT = event.getBank("REC::Event").getFloat("startTime", 0);
             
       if(!isMC && STT<=0) return;
       if( isMC && STT<0) STT=124.25f;
  
       processTL(event); 
       processRec(event);       
       
    }
    
    public boolean eventFilter(DataEvent event) {
        
        boolean goodEvent = event.hasBank("REC::Event")&&
  		                    event.hasBank("REC::Particle")&&
  		                    event.hasBank("REC::Calorimeter"); 
        return goodEvent;
    }
    
    public boolean triggerFilter() {
        int elec_trigger_sect = isMC ? 5 : getElecTriggerSector(shiftTrigBits(getRunNumber()));       
        int htcc_trigger_sect = getHTCCTriggerSector()-1;
        
        boolean goodECALSector = elec_trigger_sect>0 && elec_trigger_sect<7; 
        boolean goodHTCCSector = htcc_trigger_sect>0 && htcc_trigger_sect<7; 
        
        if(goodHTCCSector) trigger_sect = htcc_trigger_sect;
        if(goodECALSector) trigger_sect = elec_trigger_sect;
        
        return goodECALSector;	
    }
        
    public boolean pidFilter(DataEvent event, Boolean flag) {
        
    	DataBank  RecPart = event.getBank("REC::Particle");
    	
        ArrayList<Integer>  eleCandi = new ArrayList<>();
        ArrayList<Integer> pionCandi = new ArrayList<>();

        for (int ipart=0; ipart<RecPart.rows(); ipart++) {
           
            final int    pid = RecPart.getInt("pid",ipart);
            final int status = RecPart.getInt("status",ipart);       

            final boolean isFD = (int)(Math.abs(status)/1000) == 2;

            if (isFD &&  pid==11)                 eleCandi.add(ipart);
            if (isFD && (pid==212 || pid==211))  pionCandi.add(ipart);

        }
        
        if (eleCandi.isEmpty()) return false;
        
        return (eleCandi.size()==1 && (flag ? pionCandi.size()>0:true)) ;
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
       triggerPhase  = getTriggerPhase(); 

       tdcs.clear();
       
       if(event.hasBank("ECAL::tdc")){
           DataBank  bank = event.getBank("ECAL::tdc");
           for(int i = 0; i < bank.rows(); i++){
               int  is = bank.getByte("sector",i);
               int  il = bank.getByte("layer",i);
               int  ip = bank.getShort("component",i);   
               int tdc = bank.getInt("TDC",i); 

               if(tdc>0 && is>0 && is<7) {
                   if(!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Integer>(),is,il,ip);
        	           tdcs.getItem(is,il,ip).add(tdc);
                   	   float tdcdc = tps*tdc-triggerPhase; // phase corrected time
//              	   ((H2F) this.getDataGroup().getItem(is,0,0,run).getData(il-1).get(0)).fill(tps*tdc-FTOFFSET-BGOFFSET, ip); 
              	       ((H2F) this.getDataGroup().getItem(is,0,1,run).getData(il-1).get(0)).fill(tdcdc-FTOFFSET,ip); 
              	       ((H2F) this.getDataGroup().getItem(is,0,2,run).getData(il-1).get(0)).fill(tdcdc-FTOFFSET,ip); // triggered time
               }
           }
       } 
        
       if(event.hasBank("ECAL::adc")){
           DataBank  bank = event.getBank("ECAL::adc");
           for(int i = 0; i < bank.rows(); i++){
               int  is = bank.getByte("sector",i);
               int  il = bank.getByte("layer",i);
               int  ip = bank.getShort("component",i);
               int adc = Math.abs(bank.getInt("ADC",i));
               if(adc==0) continue;
               float t = bank.getFloat("time",i) + (float) tmf.getDoubleValue("offset",is,il,ip)  // FADC-TDC offset (sector, layer, PMT) 
                                                 + (float)  fo.getDoubleValue("offset",is,il,0);  // FADC-TDC offset (sector, layer) 
               
               float ftdc_corr = t+FTOFFSET-BGOFFSET;
               
               float tmax = 1000; float tdc = 1000;
                            
               if (tdcs.hasItem(is,il,ip)) { // sector,layer,component FADC/TDC match
            	   float radc = (float)Math.sqrt(adc);
                   for (float tdcc : tdcs.getItem(is,il,ip)) {
                	   float tdif = tps*tdcc - triggerPhase - (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - ftdc_corr;
                	   ((H2F) this.getDataGroup().getItem(is,0,8,run).getData(il-1).get(0)).fill(tdif,ip); // TDC - FADC time
                	   if (Math.abs(tdif)<TMFCUT && tdif<tmax) {tmax = tdif; tdc = tps*tdcc-triggerPhase;}                	    
                   }
                   double da0 = dtime.getDoubleValue("a0", is, il, ip); 
                   double da2 = dtime.getDoubleValue("a2", is, il, ip);
                   double da3 = dtime.getDoubleValue("a3", is, il, ip);
                   double da4 = dtime.getDoubleValue("a4", is, il, ip);
                   double da5 = dtime.getDoubleValue("a5", is, il, ip);
                   double da6 = dtime.getDoubleValue("a6", is, il, ip);
                   double da7 = dtime.getDoubleValue("a7", is, il, ip);
                   double da8 = dtime.getDoubleValue("a8", is, il, ip);
                   
                   double fa0 = ftime.getDoubleValue("a0", is, il, ip); 
                   double fa2 = ftime.getDoubleValue("a2", is, il, ip);
                   double fa3 = ftime.getDoubleValue("a3", is, il, ip);
                   double fa4 = ftime.getDoubleValue("a4", is, il, ip);
                   double fa5 = ftime.getDoubleValue("a5", is, il, ip);
                   double fa6 = ftime.getDoubleValue("a6", is, il, ip);   
                   
                   double x = radc;
//                   double dcorr = isMC ? 0:0;
                   double dcorr = da2 + Math.exp(-(x-da3)/da4)+1-Math.exp(-(da5-x)/da6)-Math.exp(-(x-da3*0.95)/da7)*Math.pow(x,da8);
                   double fcorr = fa2 + Math.exp(-(x-fa3)/fa4)+1-Math.exp( (x-fa5)/fa6);;
                   double dtdcmc = tdc - (float)gtw.getDoubleValue("time_walk",is,il,0)/x - dcorr - da0;
                   double ftdcmc = ftdc_corr - fcorr - fa0;
                   double tdcmc  = dtdcmc;
          	       ((H2F) this.getDataGroup().getItem(is,0,3,run).getData(il-1).get(0)).fill(tdc-FTOFFSET,ip);  // matched FADC/TDC
          	       ((H2F) this.getDataGroup().getItem(is,0,4,run).getData(il-1).get(0)).fill(tdcmc-FTOFFSET,ip); // calibrated time
          	       ((H2F) this.getDataGroup().getItem(is,0,5,run).getData(il-1).get(0)).fill(dtdcmc-ftdcmc,ip); // calibrated TDC-FADC time
          	       if(isGoodTL && il==1 && ip==3) ((H1F) this.getDataGroup().getItem(0,0,tlnum,run).getData(is+5).get(0)).fill(tdcmc); 
               }
           }
       }    	
    }
  
    public void processRec(DataEvent event) { // process reconstructed timing

       int run = getRunNumber(); 
        
       if(dropBanks) dropBanks(event); // rerun ECEngine with updated CCDB constants       
       
       if(event.hasBank("ECAL::hits")){
           DataBank  bank = event.getBank("ECAL::hits");
           for(int loop = 0; loop < bank.rows(); loop++){
              int   is = bank.getByte("sector", loop);
              int   il = bank.getByte("layer", loop); 
              int   ip = bank.getByte("strip", loop);
              float  t = bank.getFloat("time", loop);
              float tt = t+TOFFSET-FTOFFSET;
              if (tt>0) { 
                ((H2F) this.getDataGroup().getItem(is,0, 0,run).getData(il-1               ).get(0)).fill(tt, ip); 
                ((H2F) this.getDataGroup().getItem( 0,0,32,run).getData(3*(is-1)+getDet(il)).get(0)).fill(ip,getLay(il));
              }
           }
      }
       
       if (run>=2000) {triggerPhase=0; shiftTV[1]=0;} // Corrections needed until runs<4013 are recooked
       
       if(!event.hasBank("ECAL::clusters")) return;
              
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
       tresid.clear();
       
       for(int loop = 0; loop < bankc.rows(); loop++){ //loop over REC::Calorimeter
          int     is = bankc.getByte("sector",loop);
          int     il = bankc.getByte("layer",loop);
          int     in = bankc.getShort("index",loop); // index to cooked ECAL::clusters (unless dropBanks re-cooking shuffles ordering)
          int    det = bankc.getByte("detector",loop);
          if (det==7 && !pathlist.hasItem(is,il,in)) pathlist.add(loop,is,il,in); // associate ECAL::cluster index to REC::Calorimeter index               
       }
       
       int[]   iip = new int[3],   iid = new int[3],   stp = new int[3];
       float[] tid = new float[3], lef = new float[3], add = new float[3];
             
       DataBank  bank1 = event.getBank("ECAL::clusters"); //entries ordered by sector,layer
       DataBank  bank3 = event.hasBank("ECAL::peaks") ? event.getBank("ECAL::peaks") : null;
       DataBank  bank4 = event.hasBank("ECAL::calib") ? event.getBank("ECAL::calib") : null;

       for(int loop = 0; loop < bank1.rows(); loop++){ //loop over ECAL::clusters
           
             if (true) {
               int     is =  bank1.getByte("sector", loop);
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
               
               if(bank4 != null) {
            	   tid[0] = getTime(bank4,"u",loop);
            	   tid[1] = getTime(bank4,"v",loop);
            	   tid[2] = getTime(bank4,"w",loop); 
               } else {
            	   tid[0] = bank3.getFloat("time",iid[0]) - lef[0]/(float)getVeff().getDoubleValue("veff", is,il+0,iip[0]); //readout time subtracted U
            	   tid[1] = bank3.getFloat("time",iid[1]) - lef[1]/(float)getVeff().getDoubleValue("veff", is,il+1,iip[1]); //readout time subtracted V
            	   tid[2] = bank3.getFloat("time",iid[2]) - lef[2]/(float)getVeff().getDoubleValue("veff", is,il+2,iip[2]); //readout time subtracted W
               }
               
               add[0] = bank3.getFloat("energy",iid[0]); //peak energy U
               add[1] = bank3.getFloat("energy",iid[1]); //peak energy V
               add[2] = bank3.getFloat("energy",iid[2]); //peak energy W 
               
               stp[0] = bank3.getShort("status", iid[0]);
               stp[1] = bank3.getShort("status", iid[1]);
               stp[2] = bank3.getShort("status", iid[2]);
               
               if (pathlist.hasItem(is,il,loop)) { //should always be true unless dropBanks=true            
                   int    pin = bankc.getShort("pindex", pathlist.getItem(is,il,loop));
                   float path = bankc.getFloat("path",   pathlist.getItem(is,il,loop)); 
                   int    pid = bankp.getInt("pid",pin);
                   float beta = bankp.getFloat("beta",pin);                    
                   int status = Math.abs(bankp.getInt("status",pin));
                   
                   for (int i=0; i<3; i++) { //loop over U,V,W
                	   float tu=0,tdc=0,tdcc=0,tdccc=0,leff=0,adc=0; int ip=0;
                	   if (eng.clusters.size()>0) { // use ECEngine clusters object
                         tu    = (float) eng.clusters.get(loop).getTime(i);
                         ip    =         eng.clusters.get(loop).getPeak(i).getMaxStrip();
                         adc   =         eng.clusters.get(loop).getPeak(i).getMaxECStrip().getADC();
                         tdc   = (float) eng.clusters.get(loop).getPeak(i).getMaxECStrip().getRawTime(true)-TOFFSET;
                         tdcc  = (float) eng.clusters.get(loop).getPeak(i).getMaxECStrip().getTWCTime(); 
                         leff  = (float) eng.clusters.get(loop).getPeak(i).getMaxECStrip().getTdist();
                	   } else { // use ECAL::clusters and ECAL::peaks  
                		 tu    = tid[i];
                		 ip    = iip[i];
                		 leff  = lef[i];
                         adc   = 10000*add[i]/(float)gain.getDoubleValue("gain", is, il+i, ip);
                	   }
                	   
                	   float vcorr = STT - triggerPhase + TVOffset;  // phase=0 unless STT is not phase corrected (early engineering runs)
                       float tvcor = tu  - vcorr;
                      
                       float mybet = path/tvcor/c; // use ECAL beta for beta distribution plots and neutral residual plots
                     
                       if(Math.abs(pid)==211)((H1F) this.getDataGroup().getItem(0,i,22,run).getData(getDet(il)+3).get(is-1)).fill(mybet);  
                       if(pid==22||pid==2112)((H1F) this.getDataGroup().getItem(0,i,22,run).getData(getDet(il)+6).get(is-1)).fill(mybet); 
                       if(pid==11)           ((H1F) this.getDataGroup().getItem(0,i,22,run).getData(getDet(il)  ).get(is-1)).fill(mybet); 
                       
                       if(Math.abs(t-tu)<0.001) { //Choose U,V,W time tu used for cluster time
                       if(Math.abs(pid)==211)((H1F) this.getDataGroup().getItem(1,i,22,run).getData(getDet(il)+3).get(is-1)).fill(mybet);  
                       if(pid==22||pid==2112)((H1F) this.getDataGroup().getItem(1,i,22,run).getData(getDet(il)+6).get(is-1)).fill(mybet); 
                       if(pid==11)           ((H1F) this.getDataGroup().getItem(1,i,22,run).getData(getDet(il)  ).get(is-1)).fill(mybet); 
                       }
                       
                       if (TRpid==2112 && pid==2112) { };
                      
                       float vel = (Math.abs(pid)==211 || Math.abs(pid)==2212) ? Math.abs(beta*c):c; //use EB beta for pion or proton calibration residuals                       
                       
                	   float pcorr = path/vel;
                       float resid = tvcor - pcorr;  
                       
                       if(!tresid.hasItem(is,il+i,ip)) tresid.add(new ArrayList<Float>(),is,il+i,ip);
                           tresid.getItem(is,il+i,ip).add(resid);
                           
                	   float  radc = (float) Math.sqrt(adc);  
                       float lcorr = leff/(float)getVeff().getDoubleValue("veff", is, il+i, ip);
                       
                       boolean goodSector = dropEsect ? is!=trigger_sect : is==trigger_sect;  
                       
                       boolean goodPID = TRpid==11211 ? Math.abs(pid)==211 && is!=trigger_sect || 
                    		                            Math.abs(pid)==11  && is==trigger_sect : //combine e- and pi+/pi-
                    		                            Math.abs(pid)==TRpid && goodSector;      //PID from menu selection
                                              
                       boolean goodStatus = status>=2000 && status<3000; 
                       boolean goodHisto  = goodPID && goodStatus;
                       
                       float offset = 0;
                       if (is==2 && run>3030 && run<3106) offset=30;
                                              
                       if (goodHisto) {
                    	   ((H1F) this.getDataGroup().getItem(is,0,tlnum,run).getData(il+i-1).get(0)).fill(resid); //timelines
                           ((H2F) this.getDataGroup().getItem(is,0,   10,run).getData(il+i-1).get(0)).fill(resid+offset, ip); //used for calibration
                       }
                       
                       if (!dropSummary && goodHisto) {  //Menu selection
                           ((H2F) this.getDataGroup().getItem(is,0,6,run).getData(il+i-1).get(0)).fill(tu+TOFFSET-FTOFFSET, ip); //peak times
                           ((H2F) this.getDataGroup().getItem(is,0,7,run).getData(il+i-1).get(0)).fill(t +TOFFSET-FTOFFSET, ip); //cluster times
                           ((H2F) this.getDataGroup().getItem(is,0,9,run).getData(il+i-1).get(0)).fill(tvcor, ip);      // TVertex corrected time
                           ((H2F) this.getDataGroup().getItem(is,il+i,11,run).getData(ip-1).get(0)).fill(path, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,12,run).getData(ip-1).get(0)).fill(ener, resid); //cluster energy
                           ((H2F) this.getDataGroup().getItem(is,il+i,13,run).getData(ip-1).get(0)).fill(leff, resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,14,run).getData(ip-1).get(0)).fill(tu+TOFFSET, resid);  
                           ((H2F) this.getDataGroup().getItem(is,il+i,15,run).getData(ip-1).get(0)).fill(adc,  resid);
                           ((H2F) this.getDataGroup().getItem(is,il+i,16,run).getData(ip-1).get(0)).fill(tdc- vcorr-pcorr-lcorr, radc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,17,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr-lcorr, radc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,18,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr-lcorr, adc);
                           ((H2F) this.getDataGroup().getItem(is,il+i,19,run).getData(ip-1).get(0)).fill(tdc- vcorr-pcorr, leff);
                           ((H2F) this.getDataGroup().getItem(is,il+i,20,run).getData(ip-1).get(0)).fill(tdcc-vcorr-pcorr, leff);
                           ((H2F) this.getDataGroup().getItem(is,il+i,21,run).getData(ip-1).get(0)).fill(vcorr+TOFFSET+(isMC?50:0),leff); 
                           if(TRESFILL(il+i-1)==ip) {
                              ((H2F) this.getDataGroup().getItem(0,0,33,run).getData(il+i-1+(is-1)*9).get(0)).fill(adc,resid);
                           }
                       }                        
                   } //loop over U,V,W
               } //pathList check
           } //dummy boolean
       } //loop over ECAL::clusters   	
    }
    
    public int TRESFILL(int tag) {
    	switch (tag) {
    	case 0: return 10;
    	case 1: return 52;
    	case 2: return 52;
    	case 3: return  4;
    	case 4: return 33;
    	case 5: return 33;
    	case 6: return  4;
    	case 7: return 33;
    	case 8: return 33;
    	}
    	return 0;
    }
    
    public float getTime(DataBank bank4, String uvw, int index) {
    	switch (uvw) {
    	case "u": return eng.useFADCTime ? bank4.getFloat("recFTU",index): bank4.getFloat("recDTU",index); 
    	case "v": return eng.useFADCTime ? bank4.getFloat("recFTV",index): bank4.getFloat("recDTV",index);
    	case "w": return eng.useFADCTime ? bank4.getFloat("recFTW",index): bank4.getFloat("recDTW",index);
    	}
    	return 0;    	
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
	   if(trfitEnable) analyzeTRES();
    
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
        System.out.println(getDetectorName()+".analyzeTimeLineFits Finished");
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
        if (eng.useFADCTime &&  eng.usePass2Timing) {writeFile("ftime",1,7,0,3,0,3);  writeFile("fveff",1,7,0,3,0,3);}
        if(!eng.useFADCTime &&  eng.usePass2Timing) {writeFile("dtime",1,7,0,3,0,3);  writeFile("dveff",1,7,0,3,0,3);}
        if(!eng.useFADCTime && !eng.usePass2Timing) {writeFile("timing",1,7,0,3,0,3); writeFile("effective_velocity",1,7,0,3,0,3);}    
        isAnalyzeDone = true;
    }
    
    public void analyzeResiduals() {
        getResidualSummary(1,7,0,3,0,3);
        if( eng.useFADCTime &&  eng.usePass2Timing) writeFile("ftime_update",1,7,0,3,0,3);
        if(!eng.useFADCTime &&  eng.usePass2Timing) writeFile("dtime_update",1,7,0,3,0,3);
        if(!eng.useFADCTime && !eng.usePass2Timing) writeFile("timing_update",1,7,0,3,0,3);
        System.out.println("analyzeResiduals Finished");
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
    
    public void analyzeTRES() {
    	analyzeTRESGraphs(1,7,0,3,0,3);
    	if (trfitEnable &&  eng.useFADCTime &&  eng.usePass2Timing) writeFile("ftres",1,7,0,3,0,3);    	
    	if (trfitEnable && !eng.useFADCTime &&  eng.usePass2Timing) writeFile("dtres",1,7,0,3,0,3);    	
        System.out.println("analyzeTRES Finished");
    	isTRESDone = true;
    }
    
    public GraphErrors graphShift(GraphErrors g, double shift) {
    	GraphErrors gnew = new GraphErrors();      	
    	for (int i=0; i<g.getDataSize(0); i++) {
           gnew.addPoint(g.getDataX(i),g.getDataY(i)+shift,g.getDataEX(i),g.getDataEY(i));
    	}
    	return gnew;
    }
    
    
    public void analyzeTRESGraphs(int is1, int is2, int il1, int il2, int iv1, int iv2) {
    	
        ParallelSliceFitter fitter;
        GraphErrors g;
        
        int run=getRunNumber(), n=0;
        
        for (int is=is1; is<is2; is++) {            
           for (int il=il1; il<il2; il++) {
              for (int iv=iv1; iv<iv2; iv++) {
                  System.out.println("Fitting Sector "+is+" Layer "+il+" View "+iv+" "+this.getDataGroup().hasItem(0,0,33,run)); 
                  fitter = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(0,0,33,run).getData(n).get(0));
                  fitter.setRange(-3,3); 
                  fitter.fitSlicesX(); g = fitter.getSigmaSlices();
                  g.setMarkerColor(1); g.setLineColor(1);
                  g.getAttributes().setTitleX("sector "+is+" "+det[il]+" "+v[iv]+" FADC"); g.getAttributes().setTitleY("sigma (ns)");
          	      tl.fitData.add(TRESFit(g,true),0,n,400,run); //Fits to timing residual sigma vs FADC
              	  n++;
              }
           }
        }   	
    }
       
    public FitData TRESFit(GraphErrors g, boolean doFit) {
    	
    	int fitnum = 16;
    	FitData fd = new FitData(g);
    	if(g.getDataSize(0)==0) return fd;
    	fd.initFit(fitnum,0,0,1,g.getDataX(g.getDataSize(0)-1)*1.05);
    	fd.doFit = doFit;
    	fd.fitGraph("", true, false);
    	return fd;
    	
    }
    
    
    public void plotTRESSummary(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));         
        int n=0;
        FitData fd;
        
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        	for (int il=0; il<3; il++) {
        		for (int iv=0; iv<3; iv++) { 
        			c.cd(n); c.getPad(n).getAxisY().setRange(0,2); 
        			c.draw(tl.fitData.getItem(0,n,400,getRunNumber()).getGraph()); n++; 
        		}
        	}
        }
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
               	   tl.fitData.add(fitEngine(g,6,0),is,3*il+iv+1,ip,run); //LEFF fits to obtain t0 and veff
              	   
                   fitter2 = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(is,3*il+iv+1,17,run).getData(ip).get(0));
                   fitter2.setRange(-10,50); fitter2.fitSlicesY();
                   g = graphShift(fitter2.getMeanSlices(),-tl.fitData.getItem(is,3*il+iv+1,ip,run).p0); //subtract t0 just fitted
                   g.getAttributes().setTitleX("Sector "+is+" "+det[il]+" "+v[iv]+(ip+1)); g.getAttributes().setTitleY(""); 
                   if(dumpGraphs) {
                	   String nam = (eng.useFADCTime ? "ftime/":"dtime/")+"s"+is+"l"+il+"v"+iv+"p"+ip;
                	   String out = filPath+"twplots/"+getRunGroup(run)+"/"+nam+".vec";
                	   System.out.println("Writing "+out); dumpGraph(out,g);
                   }
//            	   if(!eng.useFADCTime && il>2) tl.fitData.add(fitEngine(g,16,20),is,3*il+iv+1,ip+200,run); //TW fits to DISC time 
                   boolean doFit = true;
                   
            	   if(is==6&&il==1&&iv==2&&ip==35) {System.out.println("Rejecting Sector "+is+" Layer "+il+" View "+iv+" PMT "+ip);doFit=false;}
            	   
            	   if(!eng.useFADCTime) tl.fitData.add(DTWFit(g,is,il,iv,doFit),is,3*il+iv+1,ip+200,run); //TW fits to DISC time           	   
            	   if( eng.useFADCTime) tl.fitData.add(FTWFit(g,is,il,iv),      is,3*il+iv+1,ip+200,run); //TW fits to FADC time            	   
                }
             }
          }
       } 
    }
    
    public FitData FTWFit(GraphErrors g, int is, int il, int iv) { //time walk corrections for Ftime
    	int fitnum = 13;
    	int   fmin = 10;
    	FitData fd = new FitData(g);
    	if(g.getDataSize(0)==0) return fd;
        fd.initFit(fitnum,0,1,fmin,g.getDataX(g.getDataSize(0)-1)*1.05); 
        switch (il) {
        //rgb
        case 0: fd.initFunc(0,-0.5); fd.initFunc(1,20,15,22);  fd.initFunc(2,9,7,11);   fd.initFunc(3,170,160,180);fd.initFunc(4,20,15,25); break;
        case 1: fd.initFunc(0,-0.5); fd.initFunc(1,13.9,5,22); fd.initFunc(2,5.4,5,11); fd.initFunc(3,125,10,300); fd.initFunc(4,15,5,30); break;
        case 2: fd.initFunc(0,-0.5); fd.initFunc(1,13.9,5,22); fd.initFunc(2,5.4,2,11); fd.initFunc(3,125,10,300); fd.initFunc(4,15,5,30); break;
        //rga
//        case 0: fd.initFunc(0,-0.5); fd.initFunc(1,20,18,22);  fd.initFunc(2,9,7,11);   fd.initFunc(3,170,160,180); fd.initFunc(4,20,15,25); break;
//        case 1: fd.initFunc(0,-0.5); fd.initFunc(1,13.9,5,22); fd.initFunc(2,5.4,5,11); fd.initFunc(3,125,100,180); fd.initFunc(4,15,14,25); break;
//        case 2: fd.initFunc(0,-0.5); fd.initFunc(1,17,10,18);  fd.initFunc(2,5,3,11);   fd.initFunc(3,145,110,180); fd.initFunc(4,15,14,20);
        }
        fd.doFit = true; 
        fd.fitGraph("",true,true); 
        return fd;
    }
   
    public FitData DTWFit(GraphErrors g, int is, int il, int iv, boolean doFit) { //time walk corrections for Dtime
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
        GraphErrors g1,g2,g3;
        
        int run=getRunNumber();
        
        for (int is=is1; is<is2; is++) {            
           for (int il=il1; il<il2; il++) {
              for (int iv=iv1; iv<iv2; iv++) {
                  fitter = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(is,0,10,run).getData(3*il+iv).get(0));
                  fitter.setRange(-10,10); fitter.fitSlicesY(); 
                  g1 = sliceToGraph(fitter.getMeanSlices(),il,iv); 
                  g2 = sliceToGraph(fitter.getSigmaSlices(),il,iv);
                  g3 = sliceToGraph(fitter.getMeans(),il,iv);
                  GraphErrors mean = new GraphErrors("RESIDUAL_"+is+"-"+il+" "+iv,g1.getVectorX().getArray(),
                		                                                          g1.getVectorY().getArray(),
//	                                             getGoodMean(g1.getVectorY(),g3.getVectorY(),0.2).getArray(), 
                		                                                          new double[g1.getDataSize(0)],
                		                                                          g2.getVectorY().getArray()); 
                  mean.getAttributes().setTitleX("Sector "+is+" "+det[il]+" "+v[iv]); mean.getAttributes().setTitleY("");  
             	  FitSummary.add(mean, is,il,iv,run);
              }
           }
        } 
        
        System.out.println("getResidualSummary Finished");

    }
    
    public GraphErrors findGoodMean(GraphErrors g1, GraphErrors g2) {
    	GraphErrors g3 = null;
    	return g3;
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
    
    /* WRITE CALIBRATION FILES */    
    
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
							case "timing_update":      line =  getA0( is,il,iv,ip); break; 
							case "dtime_update":       line =  getDA0(is,il,iv,ip); break; 
							case "ftime_update":       line =  getFA0(is,il,iv,ip); break; 
							case "timing":             line =  getTW(is,il,iv,ip);  break; //pass1
							case "dtime":              line =  getDTW(is,il,iv,ip,false); break; //true: update TWC but not offset
							case "ftime":              line =  getFTW(is,il,iv,ip,false); break;
							case "effective_velocity": line =  getEV(is,il,iv,ip);  break;
							case "fveff":              line =  getEV(is,il,iv,ip);  break;
							case "tmf_offset":         line =  getTMF(is,il,iv,ip); break;  
							case "fadc_offset":        line =  getGTMF(is,il,iv,ip);  
							case "setFTiming":         line =  getFTiming(is,il,iv,ip); break;
							case "A0offset":           line =  getNewA0(is,il,iv,ip); break;
							case "DefFTW":             line =  getDefFTW(is,il,iv,ip); break;
							case "DefDTW":             line =  getDefDTW(is,il,iv,ip); break;
							case "dtres":              line =  getTRES(is,il,iv,ip); break; 
							case "ftres":              line =  getTRES(is,il,iv,ip); 
							}
							if (line!=null) {
//						      System.out.println(line);
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
	
	public String getTRES(int is, int il, int iv, int ip) { //Absolute TDC residual calibration + pass1 time walk
		if (ip>0) return null;
		int in = 3*il+iv+(is-1)*9;
		if(tl.fitData.hasItem(0,in,400,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+" "
					+tl.fitData.getItem(0,in,400,getRunNumber()).p0+" " 
					+tl.fitData.getItem(0,in,400,getRunNumber()).p1+" "  
					+tl.fitData.getItem(0,in,400,getRunNumber()).p2+" " 
					+tl.fitData.getItem(0,in,400,getRunNumber()).p3+" ";
		} else {
			return is+" "+(3*il+iv+1)+" 0.0 0.0 0.0 0.0";
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
		
	public String getFTW(int is, int il, int iv, int ip, boolean update) { //Absolute FADC residual calibration + pass 2 time walk
		if(tl.fitData.hasItem(is,3*il+iv+1,ip,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "+
			    (update ? ftime.getDoubleValue("a0",is,3*il+iv+1,ip+1):tl.fitData.getItem(is,3*il+iv+1,ip,getRunNumber()).p0)				
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

	public String getA0(int is, int il, int iv, int ip) { //Adjust residual offset only
		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(time.getDoubleValue("a0", is, 3*il+iv+1, ip+1) 
				+rejectLoVW(is,il,iv,ip))+" "  //residual correction to previous a0
				+" 0.02345 "
				+time.getDoubleValue("a2", is, 3*il+iv+1, ip+1)+" " //replace with table
				+time.getDoubleValue("a3", is, 3*il+iv+1, ip+1)+" " //replace with table
				+time.getDoubleValue("a4", is, 3*il+iv+1, ip+1)+" ";//replace with table
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" 10.0 0.02345"+" 0.0 0.0 0.0";
		}
		
	}
	
	public String getDA0(int is, int il, int iv, int ip) { //Adjust residual offset only
        float offset = 0; int run = getRunNumber();
        if (is==2 && run>3030 && run<3106) offset=30;
		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(dtime.getDoubleValue("a0", is, 3*il+iv+1, ip+1)-offset 
				+rejectLoVW(is,il,iv,ip))+" "  //residual correction to previous a0
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
	
	public String getFA0(int is, int il, int iv, int ip) { //Adjust residual offset only
        float offset = 0; int run = getRunNumber();
        if (is==2 && run>3030 && run<3106) offset=30;
		if(FitSummary.hasItem(is,il,iv,getRunNumber())) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+(ftime.getDoubleValue("a0", is, 3*il+iv+1, ip+1)-offset
				+rejectLoVW(is,il,iv,ip))+" "  //residual correction to previous a0
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
	
	public float rejectLoVW(int is, int il, int iv, int ip) { // low occupancy V,W strips
		if ((iv==1 || iv==2) && ip<0) return 0f;
		return (float) FitSummary.getItem(is,il,iv,getRunNumber()).getDataY(ip)+A0offset;
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
	
	public String getDefFTW(int is, int il, int iv, int ip) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+("0.0"+" "
				+" 0.02345 "
				+"0.0"+" "
				+"0.0"+" "
				+"0.0"+" "
				+"0.0"+" "
				+"0.0");		
	}
	
	public String getDefDTW(int is, int il, int iv, int ip) {
	    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
			+("0.0"+" "
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
    
    public void plotBetaHistos(int index) {
        int run = getRunNumber();
        drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActivePC(),getActiveView(),index,run));	    
    }  
    
    public void plotStatusHistos(int index) {
        int run = getRunNumber();
        drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,run));	        	
    }
    
    /* TIMELINES */
    
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
    
//      public static void main(String[] args) {
    	
//    	ECt ect = new ECt("ECt",6684);
//      	ect.writeFile("DefFTW",1,7,0,3,0,3);    
//      	ECt ect = new ECt("ECt",2385);
//     	ect.writeFile("DefDTW",1,7,0,3,0,3);   	
//     }
     
}
