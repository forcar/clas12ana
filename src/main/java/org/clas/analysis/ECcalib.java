package org.clas.analysis;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.clas.tools.FitData;
import org.clas.tools.ParallelSliceFitter;
import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
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
import org.clas.tools.Event;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECcalib extends DetectorMonitor {

	H2F        h = null;
	DataGroup dg = null;
	
    int is,la,ic,idet,nstr;
    
    float[][][][] ecmean = new float[6][3][3][68];
    float[][][][]  ecrms = new float[6][3][3][68];
    String[]         det = new String[]{"pcal","ecin","ecou"};
    String[]           v = new String[]{"u","v","w"};
    float[]         mipc = {30,30,48};  
    float[]         mipp = {10,10,16};  
    double[]         mxc = {60,60,96};  
    double[]         mxp = {20,20,32};  
    double[]     fitLimp = { 5, 3, 6,17,17,27};
    double[]     fitLimc = {20,17,35,40,48,75};
    int[]           npmt = {68,62,62,36,36,36,36,36,36};    
    int[]          npmts = new int[]{68,36,36};
    
    Boolean         isMC = false;  
        
    List<IndexedList<List<Float>>> RDIFmap = new ArrayList<IndexedList<List<Float>>>();
    IndexedList<Float>           PixLength = new IndexedList<Float>(3);    
    List<Float>                       pmap = new ArrayList<Float>();	
    List<Particle>                    part = new ArrayList<Particle>();
    
    IndexedTable time=null, offset=null, goffset=null, gain=null, shift=null, veff=null;
    
    IndexedList<List<Particle>>     ecpart = new IndexedList<List<Particle>>(2);    
    IndexedList<H1F>            VarSummary = new IndexedList<H1F>(4);
    
	public boolean goodELEC,goodPROT,goodPBAR,goodPIP,goodPIM,goodMIP,goodNEUT,goodPHOT,goodPHOTR,goodPHOT2,goodPIPP,goodPI0;
	
//	public List<Particle> elec_ecal  = new ArrayList<Particle>();
//	public List<Particle>  pim_ecal  = new ArrayList<Particle>();
//	public List<Particle>  pip_ecal  = new ArrayList<Particle>();	
//	public List<Particle>  mip_ecal  = new ArrayList<Particle>();	
//	public List<Particle>  prot_ecal = new ArrayList<Particle>();	

    Event     ev = new Event();
	H1F       h1 = null;
	
	//FTOF
	IndexedList<List<Float>>          tdcs = new IndexedList<List<Float>>(4);
	IndexedList<List<Float>>          adcs = new IndexedList<List<Float>>(4);
	IndexedList<List<Integer>>       lapmt = new IndexedList<List<Integer>>(3); 
	IndexedList<List<Integer>>       ltpmt = new IndexedList<List<Integer>>(3); 
	
	//ECAL
	IndexedList<List<Float>>         etdcs = new IndexedList<List<Float>>(4);
	IndexedList<List<Float>>         eadcs = new IndexedList<List<Float>>(4);
	
	
    public ECcalib(String name) {
        super(name);
        this.setDetectorTabNames("MIP",
                                 "UVW",
                                 "Fits",
                                 "Mean",
                                 "RMS",
                                 "Maps",
                                 "PID",
                                 "MOM",
                                 "PCAL/ECTOT",
                                 "PathIJ",
                                 "PIXEL",
                                 "Timeline",
                                 "ATT",
                                 "MIPvP",
                                 "VAR",
                                 "Align",
                                 "Align2");
        
        this.useRDIFButtons(true);
        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
    }
    
    public void localinit() {
    	System.out.println("ECcalib.localinit()");    	
        getPixLengthMap(outPath+"files/ECpixdepthtotal.dat");
    }  
    
    public void localclear() {
    	System.out.println("ECcalib:localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	FitSummary.clear();
    	Fits.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
        engine.setCCDBGain(!defaultGain);
    }
    
     @Override    
     public void createHistos(int run) {
	     histosExist = true;
	     System.out.println("ECcalib:createHistos("+run+")");
	     setRunNumber(run);
	     runlist.add(run);
	     createMIPHistos(0,1,50,0, 40," Peak Energy (MeV)");
	     createMIPHistos(0,2,50,0,100," Cluster Energy (MeV)");	     
	     if(dropSummary) return;
	     createXYHistos(5,80,-420,420,-420,420);    
	     createPIDHistos(6);
	     createMIPHistos(7,1,25,0,5.0," + Momentum (GeV)");
	     createMIPHistos(7,2,25,0,5.0," - Momentum (GeV)");
	     createMIPHistos(7,3,25,0,5.0," + Momentum (GeV)");
	     createMIPHistos(7,4,25,0,5.0," - Momentum (GeV)");
	     createMIPHistos(7,5,25,0,5.0," + Momentum (GeV)");
	     createMIPHistos(7,6,25,0,5.0," - Momentum (GeV)");
	     createPathHistos(9);
	     createPixHistos(10);	
	     createUVWHistos(12,25,0.,2.," MIP ");
	     createMIPvPHistos(13);
	     createBARHistos(15,1);
     }
     
     @Override       
     public void plotHistos(int run) {
    	 if(!histosExist) return;
    	 plotSummary(run);
    	 plotAnalysis(run);
     }
     
     public void plotSummary(int run) {
    	 if(dropSummary) return;
    	 setRunNumber(run);
    	 plotMIP(0);      	
    	 plotXY(5); 
    	 plotPIDSummary(6);
    	 plotMIPZ(7);
    	 plotPathSummary(9,getActivePC()==1?getActiveView()+1:0);
    	 plotPathSummary(10,0);
    	 plotUVW(12);    
    	 plotPIDSummary(13);
    	 plotMIP(15);
     }
     
     public void plotAnalysis(int run) {
    	 setRunNumber(run);
    	 if(!isAnalyzeDone) return;
    	 if(!dropSummary) {
    		 updateFITS(2); 
 //   		 updateFITS(15); 
    		 if(runlist.size()==3) dumpFiles("gain");
    		 if(TLname=="UVW") {
    			 /*plotMean()*/;plotVarSummary(14); plotMeanSummary(3); plotRmsSummary(4);
    			 } else {
    				 plotMeanHWSummary(3); plotRmsHWSummary(4);
             }	 
    	 }
//    	 if(!dropSummary) {updateFITS(2);plotMeanHWSummary(3); plotRmsHWSummary(4);}
    	 updateUVW(1);   
//    	 plotAlignSummary(16); 
    	 plotTimeLines(11);  	    
     }
     
     public void plotMean() {
    	 if(getActiveRDIF()==1 && runlist.size()==2) {getMeanRDIF(); dumpFiles("rdif");     plotMeanRDIF(3);                        return;} //plot RDIF if two runs present and dump RDIF
    	 if(getActiveRDIF()==0 && runlist.size()==1) {               dumpFiles("gain");     plotMeanSummary(3); plotVarSummary(14); return;} //dump PMT gains based on current analyzed run
    	 if(getActiveRDIF()==1 && runlist.size()==1) {               dumpFiles("rdifgain"); plotMeanSummary(3); plotVarSummary(14);}         //plot effect of RDIF correction and dump RDIF corrected CCDB gains
     }

     public void dumpFiles(String val) {
    	 if(dumpFiles) writeFile(val,1,7,0,3,0,3);
     }		      
     public void createXYHistos(int k, int nb, int bx1, int bx2, int by1, int by2) {
    	 
 	     int run = getRunNumber();
 	     
         String[] t = {"e","w","r"};
         
         for (int i=0; i<3; i++) { //i=0,1,2 event,weight,ratio
        	 dg = new DataGroup(3,2);
        	 for (int d=0; d<3; d++) { //id=0,1,2 pcal,ecin,ecou clusters
        		 h = new H2F("hi-"+det[d]+"-xyc-"+t[i]+"-"+k+"-"+run,"hi-"+det[d]+"-xyc-"+t[i]+"-"+k+"-"+run,nb,bx1,bx2,nb,by1,by2);
                 dg.addDataSet(h,d);  
	    	 }
             this.getDataGroup().add(dg,i,2,k,run);
         }

         for (int i=0; i<3; i++) { //i=0,1,2 event,weight,ratio
        	 dg = new DataGroup(3,3);
        	 for (int j=0; j<3; j++) { //j=0,1,2 U,V,W peaks
        		 for (int d=0; d<3; d++) { //id=0,1,2 pcal,ecin,ecou
        			 h = new H2F("hi-"+det[d]+"-xyp-"+v[j]+t[i]+"-"+k+"-"+run,"hi-"+det[d]+"-xyp-"+v[j]+t[i]+"-"+k+"-"+run,nb,bx1,bx2,nb,by1,by2);
                     dg.addDataSet(h,j+d*3);                      
        	     } 
        	 }
             this.getDataGroup().add(dg,i,1,k,run);
         }
          
     }
     
     public void createUVWHistos(int k, int ybins, double ymin, double ymax, String ytxt) {
     	
         F1D f1; 
         
         int run = getRunNumber();
         
         for (int is=1; is<7; is++) {      
             DataGroup dg1 = new DataGroup(8,8); DataGroup dg2 = new DataGroup(8,8); DataGroup dg3 = new DataGroup(8,8);        	    
             f1 = new F1D("p0"+is+1+k,"[a]",1,68); f1.setParameter(0,1); f1.setLineColor(1); f1.setLineStyle(1);
            
             for (int ip=1; ip<npmts[0]+1; ip++) {
                 h = new H2F("uvw-pcal-u"+ip+"-s"+is+"-"+k+"-"+run,"uvw-pcal-u"+ip+"-s"+is+"-"+k+"-"+run,68,1,69,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" PCAL W"); h.setTitleY(ytxt+"U"+ip);       
                 dg1.addDataSet(h,ip-1); dg1.addDataSet(f1,ip-1);
                 h = new H2F("uvw-pcal-v"+ip+"-s"+is+"-"+k+"-"+run,"uvw-pcal-v"+ip+"-s"+is+"-"+k+"-"+run,68,1,69,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" PCAL U"); h.setTitleY(ytxt+"V"+ip);
                 dg2.addDataSet(h,ip-1); dg2.addDataSet(f1,ip-1);
                 h = new H2F("uvw-pcal-w"+ip+"-s"+is+"-"+k+"-"+run,"uvw-pcal-w"+ip+"-s"+is+"-"+k+"-"+run,68,1,69,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" PCAL V"); h.setTitleY(ytxt+"W"+ip); 
                 dg3.addDataSet(h,ip-1); dg3.addDataSet(f1,ip-1);
      	     }
             this.getDataGroup().add(dg1,is,1,k,run); this.getDataGroup().add(dg2,is,2,k,run); this.getDataGroup().add(dg3,is,3,k,run);
             
             DataGroup dg4 = new DataGroup(6,6); DataGroup dg5 = new DataGroup(6,6); DataGroup dg6 = new DataGroup(6,6);        	         	   
             f1 = new F1D("p0"+is+2+k,"[a]",1,37); f1.setParameter(0,1); f1.setLineColor(1); f1.setLineStyle(1);
      	     for (int ip=1; ip<npmts[1]+1; ip++) {
                 h = new H2F("uvw-ecin-u"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecin-u"+ip+"-s"+is+"-"+k+"-"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECIN W");  h.setTitleY(ytxt+"U"+ip); 
                 dg4.addDataSet(h,ip-1); dg4.addDataSet(f1,ip-1);
                 h = new H2F("uvw-ecin-v"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecin-v"+ip+"-s"+is+"-"+k+"-"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECIN U"); h.setTitleY(ytxt+"V"+ip); 
                 dg5.addDataSet(h,ip-1); dg5.addDataSet(f1,ip-1);
                 h = new H2F("uvw-ecin-w"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecin-w"+ip+"-s"+is+"-"+k+"-"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECIN V"); h.setTitleY(ytxt+"W"+ip);
                 dg6.addDataSet(h,ip-1); dg6.addDataSet(f1,ip-1);                 
      	     }
             this.getDataGroup().add(dg4,is,4,k,run); this.getDataGroup().add(dg5,is,5,k,run); this.getDataGroup().add(dg6,is,6,k,run);
      	   
             DataGroup dg7 = new DataGroup(6,6); DataGroup dg8 = new DataGroup(6,6); DataGroup dg9 = new DataGroup(6,6);        	         	   
             f1 = new F1D("p0"+is+3+k,"[a]",1,37); f1.setParameter(0,1); f1.setLineColor(1); f1.setLineStyle(1);
      	     for (int ip=1; ip<npmts[2]+1; ip++) {
                 h = new H2F("uvw-ecou-u"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecou-u"+ip+"-s"+is+"-"+k+"-"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECOU W"); h.setTitleY(ytxt+"U"+ip);
                 dg7.addDataSet(h,ip-1); dg7.addDataSet(f1,ip-1);
                 h = new H2F("uvw-ecou-v"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecou-v"+ip+"-s"+is+"-"+k+"-"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECOU U");  h.setTitleY(ytxt+"V"+ip);
                 dg8.addDataSet(h,ip-1); dg8.addDataSet(f1,ip-1);
                 h = new H2F("uvw-ecou-w"+ip+"-s"+is+"-"+k+"-"+run,"uvw-ecou-w"+ip+"-s"+is+"-"+k+"-"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECOU V"); h.setTitleY(ytxt+"W"+ip);
                 dg9.addDataSet(h,ip-1); dg9.addDataSet(f1,ip-1);    
      	     }
             this.getDataGroup().add(dg7,is,7,k,run); this.getDataGroup().add(dg8,is,8,k,run); this.getDataGroup().add(dg9,is,9,k,run);     	    
         }        
     }     
    
     public void createMIPHistos(int k, int n, int nch, double x1, double x2, String txt) {
    	
	    int run = getRunNumber();
    
        for (int is=1; is<7; is++) {
            String tag = is+"-"+n+"-"+k+"-"+run;
            dg = new DataGroup(3,3);
            h = new H2F("mip-pcal-u-"+tag,"mip-pcal-u-"+tag, nch, x1, x2, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL U"+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("mip-pcal-v-"+tag,"mip-pcal-v-"+tag, nch, x1, x2, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL V"+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("mip-pcal-w-"+tag,"mip-pcal-w-"+tag, nch, x1, x2, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL W"+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
        
            h = new H2F("mip-ecin-u-"+tag,"mip-ecin-u-"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN U"+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("mip-ecin-v-"+tag,"mip-ecin-v-"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN V"+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("mip-ecin-w-"+tag,"mip-ecin-w-"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN W"+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
        
            h = new H2F("mip-ecou-u-"+tag,"mip-ecou-u-"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU U"+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("mip-ecou-v-"+tag,"mip-ecou-v-"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU V"+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("mip-ecou-w-"+tag,"mip-ecou-w-"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU W"+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,n,k,run);
        }            
    }
     
     public void createBARHistos(int k, int n) {
     	
	    int run = getRunNumber();
	    String lab1=" FTOF1B BAR";
    
        for (int is=1; is<7; is++) {
            String tag = is+"-"+n+"-"+k+"-"+run;
            dg = new DataGroup(3,5);
            int nch=62, x1=1, x2=63;
            		
            h = new H2F("ftof-pcal-u-"+tag,"ftof-pcal-u-"+tag, nch, x1, x2, 68, 1, 69);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("PCAL U Strip"); 
            dg.addDataSet(h,0);  
            h = new H2F("ftof-pcal-v-"+tag,"ftof-pcal-v-"+tag, nch, x1, x2, 62, 1, 63);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("PCAL V Strip");        
            dg.addDataSet(h,1);            
            h = new H2F("ftof-pcal-w-"+tag,"ftof-pcal-w-"+tag, nch, x1, x2, 62, 1, 63);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("PCAL W Strip");  
            dg.addDataSet(h,2); 
        
            h = new H2F("ftof-ecin-u-"+tag,"ftof-ecin-u-"+tag, nch, x1, x2, 36, 1, 37);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("ECIN U Strip");    
            dg.addDataSet(h,3);  
            h = new H2F("ftof-ecin-v-"+tag,"ftof-ecin-v-"+tag, nch, x1, x2, 36, 1, 37);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("ECIN V Strip");        
            dg.addDataSet(h,4);            
            h = new H2F("ftof-ecin-w-"+tag,"ftof-ecin-w-"+tag, nch, x1, x2, 36, 1, 37);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("ECIN W Strip");  
            dg.addDataSet(h,5); 
        
            h = new H2F("ftof-ecou-u-"+tag,"ftof-ecou-u-"+tag, nch, x1, x2, 36, 1, 37);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("ECOU U Strip");    
            dg.addDataSet(h,6);  
            h = new H2F("ftof-ecou-v-"+tag,"ftof-ecou-v-"+tag, nch, x1, x2, 36, 1, 37);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("ECOU V Strip");        
            dg.addDataSet(h,7);            
            h = new H2F("ftof-ecou-w-"+tag,"ftof-ecou-w-"+tag, nch, x1, x2, 36, 1, 37);
            h.setTitleX("Sector "+is+lab1); h.setTitleY("ECOU W Strip");  
            dg.addDataSet(h,8);   
            
            nch=36; x1=1; x2=37;
            h = new H2F("pcal-ecin-u-"+tag,"pcal-ecin-u-"+tag, 68, 1, 69, nch, x1, x2);
            h.setTitleY("Sector "+is+" ECIN U Strip"); h.setTitleX("PCAL U Strip"); 
            dg.addDataSet(h,9);  
            h = new H2F("pcal-ecin-v-"+tag,"pcal-ecin-v-"+tag, 62, 1, 63, nch, x1, x2);
            h.setTitleY("Sector "+is+" ECIN V Strip"); h.setTitleX("PCAL V Strip");        
            dg.addDataSet(h,10);            
            h = new H2F("pcal-ecin-w-"+tag,"pcal-ecin-w-"+tag, 62, 1, 63, nch, x1, x2);
            h.setTitleY("Sector "+is+" ECIN W Strip"); h.setTitleX("PCAL W Strip");  
            dg.addDataSet(h,11);
            
            h = new H2F("pcal-ecou-u-"+tag,"pcal-ecou-u-"+tag, 68, 1, 69, nch, x1, x2);
            h.setTitleY("Sector "+is+" ECOU U Strip"); h.setTitleX("PCAL U Strip"); 
            dg.addDataSet(h,12);  
            h = new H2F("pcal-ecou-v-"+tag,"pcal-ecou-v-"+tag, 62, 1, 63, nch, x1, x2);
            h.setTitleY("Sector "+is+" ECOU V Strip"); h.setTitleX("PCAL V Strip");        
            dg.addDataSet(h,13);            
            h = new H2F("pcal-ecou-w-"+tag,"pcal-ecou-w-"+tag, 62, 1, 63, nch, x1, x2);
            h.setTitleY("Sector "+is+" ECOU W Strip"); h.setTitleX("PCAL W Strip");  
            dg.addDataSet(h,14);
            
            this.getDataGroup().add(dg,is,n,k,run);          
        }            
    }  
     
     public void createVARHistos(int k, int n, int nch, double x1, double x2, String txt) {
     	
	    int run = getRunNumber();
    
        for (int is=1; is<7; is++) {
            String tag = is+"-"+n+"-"+k+"-"+run;
            dg = new DataGroup(3,3);
            h1 = new H1F("mip-pcal-u-"+tag,"mip-pcal-u-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" PCAL U"+txt);
            dg.addDataSet(h,0);  
            h1 = new H1F("mip-pcal-v-"+tag,"mip-pcal-v-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" PCAL V"+txt);      
            dg.addDataSet(h,1);            
            h1 = new H1F("mip-pcal-w-"+tag,"mip-pcal-w-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" PCAL W"+txt);
            dg.addDataSet(h,2); 
        
            h1 = new H1F("mip-ecin-u-"+tag,"mip-ecin-u-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" ECIN U"+txt);    
            dg.addDataSet(h,3);  
            h1 = new H1F("mip-ecin-v-"+tag,"mip-ecin-v-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" ECIN V"+txt);      
            dg.addDataSet(h,4);            
            h1 = new H1F("mip-ecin-w-"+tag,"mip-ecin-w-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" ECIN W"+txt);   
            dg.addDataSet(h,5); 
        
            h1 = new H1F("mip-ecou-u-"+tag,"mip-ecou-u-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" ECOU U"+txt); 
            dg.addDataSet(h,6);  
            h1 = new H1F("mip-ecou-v-"+tag,"mip-ecou-v-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" ECOU V"+txt);         
            dg.addDataSet(h,7);            
            h1 = new H1F("mip-ecou-w-"+tag,"mip-ecou-w-"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" ECOU W"+txt);  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,n,k,run);
        }            
    } 
     
    public void createPathHistos(int k) {
    	
	   int run = getRunNumber();
       
	   for (int id=0; id<3; id++) {
		   String d1 = det[id], d2=" "+d1.toUpperCase()+"/MIP";
       for (int il=0; il<4; il++) {
           dg = new DataGroup(6,3); 
       for (int is=1; is<7; is++) {
           String tag = is+"-"+il+"-"+k+"-"+run;
           h = new H2F("hi-"+d1+"-path12-"+tag,"hi-"+d1+"-path12-"+tag,40,0,2,70,32,50);
           h.setTitleX("Sector "+is+d2); h.setTitleY("Path12 (cm)");
           dg.addDataSet(h, is-1);  
           h = new H2F("hi-"+d1+"-path13-"+tag,"hi-"+d1+"-path13-"+tag,40,0,2,70,51,75);
           h.setTitleX("Sector "+is+d2); h.setTitleY("Path13 (cm)");
           dg.addDataSet(h, is+5);  
           h = new H2F("hi-"+d1+"-path23-"+tag,"hi-"+d1+"-path23-"+tag,40,0,2,70,18,35);
           h.setTitleX("Sector "+is+d2); h.setTitleY("Path23 (cm)");
           dg.addDataSet(h, is+11);  
       }
        this.getDataGroup().add(dg,il,id,k,run);
       }
	   }                
    }
    
    public void createPixHistos(int k) {
    	
 	    int run = getRunNumber();
        
 	    for (int id=0; id<3; id++) {
 	    	String d1 = det[id], d2=" "+d1.toUpperCase();
 	    	dg = new DataGroup(6,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"-"+k+"-"+run;
            h = new H2F("hi-"+d1+"-pix1-"+tag,"hi-"+d1+"-pix1-"+tag,50,1,200,9,3,12);
            h.setTitleX("Sector "+is+d2+" (MeV)"); h.setTitleY("No. Strips");
            dg.addDataSet(h, is-1);   
            h = new H2F("hi-"+d1+"-pix2-"+tag,"hi-"+d1+"-pix2-"+tag,70,32,50,9,3,12);
            h.setTitleX("Sector "+is+" Path12 (cm)"); h.setTitleY("No. Strips");
            dg.addDataSet(h, is+5);   
            h = new H2F("hi-"+d1+"-pix3-"+tag,"hi-"+d1+"-pix3-"+tag,70,51,75,9,3,12);
            h.setTitleX("Sector "+is+" Path13 (cm)"); h.setTitleY("No. Strips");
            dg.addDataSet(h, is+11);   
            h = new H2F("hi-"+d1+"-pix4-"+tag,"hi-"+d1+"-pix4-"+tag,70,18,35,9,3,12);
            h.setTitleX("Sector "+is+" Path23 (cm)"); h.setTitleY("No. Strips");
            dg.addDataSet(h, is+17);   
        }
        this.getDataGroup().add(dg,0,id,k,run);
 	    }    	
    }
    
    public void createPIDHistos(int k)  {
    
 	    int run = getRunNumber();
	    int is  = 0;
        String tag = is+"-"+run;
        
        dg = new DataGroup(4,2);
        
        F1D f1 = new F1D("f-1"+tag,"1/(1+[a]^2/x^2)^0.5", 0.41,3.5); f1.setParameter(0,0.13957); f1.setLineColor(1); f1.setLineStyle(1);    
        F1D f2 = new F1D("f-2"+tag,"1/(1+[a]^2/x^2)^0.5", 0.41,3.5); f2.setParameter(0,0.93827); f2.setLineColor(1); f2.setLineStyle(1);   
        h = new H2F("pid-pos-"+tag,"pid-pos-"+tag,100,0.,3.5,100,0.4,1.5);       h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 0); dg.addDataSet(f1,0); dg.addDataSet(f2,0); 
        h = new H2F("pid-neg-"+tag,"pid-neg-"+tag,100,0.,3.5,100,0.4,1.5);       h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 4); dg.addDataSet(f1,4); 
        h = new H2F("pid-fc-pos-"+tag,"pid-fc-pos-"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 1); dg.addDataSet(f1,1); dg.addDataSet(f2,1); 
        h = new H2F("pid-fc-neg-"+tag,"pid-fc-neg-"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 5); dg.addDataSet(f1,5); 
        h = new H2F("pid-fc-ppi-"+tag,"pid-fc-ppi-"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 2); dg.addDataSet(f1,2); dg.addDataSet(f2,2); 
        h = new H2F("pid-fc-npi-"+tag,"pid-fc-npi-"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 6); dg.addDataSet(f1,6); 
        h1 = new H1F("pid-fc1-ppi-"+tag,"pid-fc1-ppi-"+tag,100,-0.5,1.1); h1.setTitleX("Mass^2 (GeV)"); 
        dg.addDataSet(h1, 3);  
        h1 = new H1F("pid-fc1-npi-"+tag,"pid-fc1-npi-"+tag,100,-0.5,1.1); h1.setTitleX("Mass^2 (GeV)");  
        dg.addDataSet(h1, 7);  
       
        this.getDataGroup().add(dg,is,0,k,run);
    }
    
    public void createMIPvPHistos(int k) {
 	    int run = getRunNumber(); 
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"-"+k+"-"+run;
            h = new H2F("hi-pcal-dedx-pm"+tag,"hi-pcal-dedx-pm"+tag,50,0.3,6.,30,0.5,2);
            h.setTitleX("Sector "+is+" P (GeV)");
            h.setTitleY("MEAN/MIP");
            dg.addDataSet(h, is-1);   
            h = new H2F("hi-pcal-dedx-pp"+tag,"hi-pcal-dedx-pp"+tag,100,0.3,6.,30,0.5,2);
            h.setTitleX("Sector "+is+" P (GeV)");
            h.setTitleY("MEAN/MIP");
            dg.addDataSet(h, is+5);   
         }
        this.getDataGroup().add(dg,0,0,k,run); 
        
    }
        
    public void initCCDB(int runno) {
    	System.out.println("ECcalib.initCCDB("+runno+")");
        gain    = cm.getConstants(2,     "/calibration/ec/gain");       
        time    = cm.getConstants(runno, "/calibration/ec/timing");
        veff    = cm.getConstants(runno, "/calibration/ec/effective_velocity");
        offset  = cm.getConstants(runno, "/calibration/ec/fadc_offset");
        goffset = cm.getConstants(runno, "/calibration/ec/fadc_global_offset");    	
        shift   = cm.getConstants(runno, "/calibration/ec/torus_gain_shift");        
    }
    
	public void myinit(){
		goodPROT = false; goodPBAR = false; goodPIP  = false; goodNEUT = false; goodPHOT = false; goodPI0=false;  
    }
	
    public void processEvent(DataEvent event) {
    	
    	int run = getRunNumber();
    	
    	if(dropBanks) dropBanks(event);
    	 
    	ev.init(event);
        ev.isMC = (getRunNumber()<100) ? true:false;    	
    	ev.setHipoEvent(isHipo3Event);
    	ev.setEventNumber(getEventNumber());
    	ev.requireOneElectron(false);
    	ev.setElecTriggerSector(ev.getElecTriggerSector());
        
    	if(!ev.procEvent(event)) return;
    	
 	    this.myinit();
	    
        List<Particle>  part_ecal = new ArrayList<Particle>();
        
	    if(ev.isPhys) {part_ecal.addAll(makeELEC()); part_ecal.addAll(makePROT()); part_ecal.addAll(makeMIP());}
	
	    if(ev.isMuon)  part_ecal.addAll(makeMUON()); 
	    
	    fillHists(run, part_ecal, event);    	
    }
    
    public List<Particle> makeMUON() {
    	List<Particle> olist = new ArrayList<Particle>();        
        for (Particle p : ev.getParticle(13)) if(p.getProperty("status")>=2000 && p.getProperty("status")<3000) olist.add(p);         
        return olist;    	   	
    }
    
    public List<Particle> makeELEC() {
    	
    	List<Particle> olist = new ArrayList<Particle>();        
        for (Particle p : ev.getParticle(11)) if(p.getProperty("status")>=2000 && p.getProperty("status")<3000 && p.p()>0.5) olist.add(p);         
        return olist;    	
    }    
    
    public List<Particle> makeMIP() {   	
    	
    	List<Particle> olist = new ArrayList<Particle>();
    	olist.addAll(makePIP());
    	olist.addAll(makePIM());    	
    	return olist;   	
    }
    
    public List<Particle> makePIP() {
    	
    	List<Particle> olist = new ArrayList<Particle>();        
        for (Particle p : ev.getParticle(211)) if(p.getProperty("status")>=2000 && p.getProperty("status")<3000 && p.p()>0.1) olist.add(p);                 
        return olist;
    }    
    
    public List<Particle> makePIM() {
    	
    	List<Particle> olist = new ArrayList<Particle>();        
        for (Particle p : ev.getParticle(212)) if(p.getProperty("status")>=2000 && p.getProperty("status")<3000 && p.p()>0.1) olist.add(p);                
        return olist;    	
    } 
    
    public List<Particle> makePROT() {
    	
    	List<Particle> olist = new ArrayList<Particle>();
        for (Particle p : ev.getParticle(2212)) if(p.getProperty("status")>=2000 && p.getProperty("status")<3000 && p.p()>0.1) olist.add(p);        
        return olist;    	
    }
    
    public List<Particle>  makeNEUTRAL() {
    	
    	List<Particle> olist = new ArrayList<Particle>();
    	
        for (Particle p : ev.getParticle(2112)) {
            short status = (short) p.getProperty("status");         
        	if(status>=2000 && status<3000) olist.add(p);
        }   
        
        for (Particle p :ev.getParticle(22)) {
            short status = (short) p.getProperty("status");
            if(status>=2000 && status<3000 && p.e()>0.01) olist.add(p); 
        } 
        
        return olist;    	
    }
    
    public IndexedList<List<Particle>> getECALClusters(List<Particle> list) {
    	
        IndexedList<List<Particle>> olist = new IndexedList<List<Particle>>(2);       
    	for (Particle p : list) {
    		int ip = (int)p.getProperty("pindex");
    		for (Particle ec : ev.getECAL(ip)) {
    			int is = (int) ec.getProperty("sector");
    			int il = (int) ec.getProperty("layer");
    			if(ev.isPhys && il==1) fillPID(getRunNumber(),1,ip);
    			if (!olist.hasItem(is,il)) olist.add(new ArrayList<Particle>(), is,il); 
    			     olist.getItem(is,il).add(ec);
    		}
    	}
    	return olist;
    }
    
    public IndexedList<List<Particle>> filterECALClusters(int pid, IndexedList<List<Particle>> list) {
    	
        IndexedList<List<Particle>> olist = new IndexedList<List<Particle>>(1);       
		IndexGenerator ig = new IndexGenerator();
		
    	for (Map.Entry<Long,List<Particle>>  entry : list.getMap().entrySet()){
			int is = ig.getIndex(entry.getKey(), 0); 
			int ip = (int) entry.getValue().get(0).getProperty("pindex");
			if(entry.getValue().size()==1 && Math.abs(ev.part.get(ip).getProperty("ppid"))==pid) {
				if(!olist.hasItem(is)) olist.add(new ArrayList<Particle>(), is); 
				    olist.getItem(is).add(entry.getValue().get(0));
			}
    	}
    	return olist;    	
    }
    
    public void fillHists(int run, List<Particle> list, DataEvent event) {
    	
    	List<Vector3> rl = new ArrayList<Vector3>();  
		IndexGenerator ig = new IndexGenerator();
		       	
       	float[]   uvw = new float[3]; float[]  wuv = new float[3]; float[] ep = new float[3];
       	float[]  puvw = new float[3];
       	boolean[] fid = new boolean[3];
       	
		int is,il,ip,trigger=0,trig=TRpid;
		float ecl,x,y,z,d,wu,wv,ww,wsum,pmip=0,v12mag,v13mag,v23mag,beta,mass2=0;
		
		IndexedList<List<Particle>> ecpart = new IndexedList<List<Particle>>(1);
		
		if (ev.isPhys) {
	        trigger = (int) ev.part.get(0).getProperty("ppid");
			if (Math.abs(trigger)!=trig) return;
			fillPID(run,0,0);
		}
		
		ecpart = filterECALClusters(ev.isPhys ? 211:13,getECALClusters(list));
    	
    	for (Map.Entry<Long,List<Particle>>  entry : ecpart.getMap().entrySet()){ //loop over sectors
			is = ig.getIndex(entry.getKey(), 0);
            if(ecpart.getItem(is).size()==3) { //Require PCAL,ECIN,ECOU
            	
        		getFTOFADC(run,is,event); List<Integer> fbars = getFTOFBAR(100);
        		
            	rl.clear();
            	for (Particle p : entry.getValue()) rl.add(new Vector3(p.getProperty("x"),p.getProperty("y"),p.getProperty("z")));
            	
                Vector3  v1 = new Vector3(rl.get(0).x(),rl.get(0).y(),rl.get(0).z());
            	Vector3  v2 = new Vector3(rl.get(1).x(),rl.get(1).y(),rl.get(1).z());
            	Vector3  v3 = new Vector3(rl.get(2).x(),rl.get(2).y(),rl.get(2).z());
            	Vector3 v23 = new Vector3(v2.x(),v2.y(),v2.z());

            	v2.sub(v1); v23.sub(v3); v3.sub(v1);  
            
            	v12mag = (float) v2.mag();
            	v13mag = (float) v3.mag();
            	v23mag = (float) v23.mag();
            	
            	int[][] nucut = new int[][]{{60,60,60,60,60,60},{35,36,36,35,36,36},{35,36,35,35,35,36}}; //U near cuts
            	int[][] nvcut = new int[][]{{67,67,67,67,67,67},{36,36,35,36,36,36},{36,35,35,36,35,35}}; //V near cuts
            	int[][] nwcut = new int[][]{{60,60,61,60,60,61},{36,36,36,35,36,36},{34,35,35,35,35,35}}; //W near cuts
            	
            	int[][] fucut = new int[][]{{60,60,60,60,60,61},{35,36,35,35,36,36},{34,35,35,35,36,35}}; //U far cuts
            	int[][] fvcut = new int[][]{{60,60,60,60,60,60},{35,36,36,35,36,36},{35,36,36,35,35,36}}; //V far cuts
            	int[][] fwcut = new int[][]{{67,66,67,67,67,67},{36,35,35,36,36,36},{36,35,35,35,35,35}}; //W far cuts
            	
            	for (Particle ec : entry.getValue()) {	//Loop over PCAL, ECIN, ECOU			
            		ip     = (int)   ec.getProperty("pindex");
            		il     = (int)   ec.getProperty("layer"); int ild = getDet(il);
            		
            		x      = (float) ec.getProperty("x");
            		y      = (float) ec.getProperty("y");
            		z      = (float) ec.getProperty("z");
            		
            		uvw[0] = (float) ec.getProperty("iu");
            		uvw[1] = (float) ec.getProperty("iv");
            		uvw[2] = (float) ec.getProperty("iw");
            		
            		wuv[0] = uvw[2]; // use W strips to measure U readout distance
            		wuv[1] = uvw[0]; // use U strips to measure V readout distance
            		wuv[2] = uvw[1]; // use V strips to measure W readout distance 
            		
                	fid[0] = wuv[0]>nucut[ild][is-1] || uvw[1]>fucut[ild][is-1]; //flag near and far ends of U readout strip to reject corner clippers
                	fid[1] = wuv[1]>nvcut[ild][is-1] || uvw[2]>fvcut[ild][is-1]; //flag near and far ends of V readout strip to reject corner clippers
                	fid[2] = wuv[2]>nwcut[ild][is-1] || uvw[0]>fwcut[ild][is-1]; //flag near and far ends of W readout strip to reject corner clippers
                	
            		wu     = (float) ec.getProperty("du");
            		wv     = (float) ec.getProperty("dv");
            		ww     = (float) ec.getProperty("dw");
            		
                    wsum   = wu+wv+ww;
                    
            		float e_cz     = ec.hasProperty("cz")?(float) ec.getProperty("cz"):0;
            		d      = (il>1 && ev.isMuon && normPix && PixLength.hasItem((int)uvw[0],(int)uvw[1],(int)uvw[2]))? PixLength.getItem((int)uvw[0],(int)uvw[1],(int)uvw[2]):1f; 
            		ecl    = (float) ec.getProperty("energy")/d;	    	   
            		ep[0]  = (float) ec.getProperty("receu")/d;
            		ep[1]  = (float) ec.getProperty("recev")/d;
            		ep[2]  = (float) ec.getProperty("recew")/d;
            		
            		
            		if (ev.isPhys) {
            			pmip   = (float) ev.part.get(ip).p() * ev.part.get(ip).charge();
            			beta   = (float) ev.part.get(ip).getProperty("beta");
            			mass2  = (float) (pmip*pmip*(1f/(beta*beta)-1));
            		}
            		
            		if(ev.isPhys && il==1) fillPID(run,2,ip);
            		if(this.getDataGroup().hasItem(0,0,9,run)) fillPATH(is,il,run,ecl,ep,wsum,v12mag,v13mag,v23mag,e_cz,fid);
            		
//            		Boolean  pid = ev.isMuon ? true : pmip!=0 && Math.abs(mass2-0.0193)<0.03;
//            		Boolean  mip = il==1?(ev.isMuon ? v12mag<35&&wsum==3 : v13mag<56&&(wsum==3||wsum==4)) : v23mag<24&&wsum==3; 

            		Boolean mip = true, pid = true;
            		Boolean  pix = il==1 ? wu==1 && wv==1 && ww==1 : wsum==3 || wsum==4;
            		
            		if (mip && pid) fillMIP(is,il,run,uvw,wuv,fid,ecl,ep,pmip,x,y);
            		
            		if (pix) {
            		for (Integer bar : fbars) for (int i=0; i<3; i++) ((H2F) this.getDataGroup().getItem(is,1,15,run).getData(i+il-1).get(0)).fill(bar,uvw[i]);            		
            		if (il==1) {puvw[0]=uvw[0]; puvw[1]=uvw[1]; puvw[2]=uvw[2];}
            		if (il==4) for (int i=0; i<3; i++) ((H2F) this.getDataGroup().getItem(is,1,15,run).getData(i+il-1+6).get(0)).fill(puvw[i],uvw[i]);
            		if (il==7) for (int i=0; i<3; i++) ((H2F) this.getDataGroup().getItem(is,1,15,run).getData(i+il-1+6).get(0)).fill(puvw[i],uvw[i]);
            		}
            	}
			}  
    	}
    }
    
    
  //Target-PCAL: 6977.8 mm CCDB:/geometry/pcal/dist2tgt
  //Target-ECAL: 7303.3 mm CCDB:/geometry/ec/dist2tgt
  //PCAL-ECin: 325.5 mm
  //ECin-ECou: 14*2.2+15*10.0 = 180.8 mm
  //PCAL-ECou: 325.5+180.8= 506.3 mm
    
    public void fillPID(int run, int fill, int ip) {
    	
        if(ev.part.isEmpty()) return;
        
    	if(fill==0) {
    		for (Particle p : ev.part) {
    			if(p.getProperty("ppid")!=0 && p.getProperty("index")>0 && p.charge()!=0 && p.getProperty("status")>=2000) {
    				((H2F) this.getDataGroup().getItem(0,0,6,run).getData(p.charge()>0?0:4).get(0)).fill(p.p(),p.getProperty("beta"));
            		float beta = (float) p.getProperty("beta");
                    float mass2 =(float) (p.p()*p.p()*(1f/(beta*beta)-1));
                    ((H1F) this.getDataGroup().getItem(0,0,6,run).getData(p.charge()>0?3:7).get(0)).fill(mass2);
            	}
            }
        }

    	if(fill>0) {
    		Particle p = ev.part.get(ip);
    		float pm = (float) p.p(); float b = (float) p.getProperty("beta"); int pid = (int) p.getProperty("ppid"); int q = p.charge();
    		if(fill==1) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(q>0?     1:5).get(0)).fill(pm,b);   		
    		if(fill==2) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(pid==211?2:6).get(0)).fill(pm,b);
    	}
    }
    
    public void fillMIP(int is, int il, int run, float[] uvw, float[] wuv, boolean fid[], float ec, float[] ep, float p, float x, float y) {
    	
    	boolean fidc = !(fid[0] || fid[1] || fid[2]);
    	
    	int il3=il/3;
    	((H2F) this.getDataGroup().getItem(0,  2,5,run).getData(il3).get(0)).fill(-x,y,ec<mxc[il3]&&fidc?mipc[il3]:0);
    	((H2F) this.getDataGroup().getItem(1,  2,5,run).getData(il3).get(0)).fill(-x,y,ec<mxc[il3]&&fidc?ec:0); 
    	
    	if(il==1) ((H2F) this.getDataGroup().getItem(0,0,13,run).getData(is+(p<0?0:6)-1).get(0)).fill(Math.abs(p),fidc?ec:0/mipc[0]);

    	for (int i=0; i<3; i++) {
    		((H2F) this.getDataGroup().getItem(is,2,0,run).getData(i+il-1).get(0)).fill(fidc?ec:-10,uvw[i]);    	
        	if (!fid[i]) {
    			((H2F) this.getDataGroup().getItem(is,1,0,run).getData(i+il-1).get(0)).fill(ep[i],uvw[i]);	
    			((H2F) this.getDataGroup().getItem(is,i+il,12,run).getData((int)uvw[i]-1).get(0)).fill(wuv[i],ep[i]/mipp[il3]);	
//    	    	float z1 = ep[i]<mxp[il3]?mipp[il3]:0, z2 = ep[i]<mxp[il3]?ep[i]:0;
    			float z1 = ep[i]<mxp[il3]?1:0, z2 = ep[i]<mxp[il3]?ep[i]/mipp[il3]:0;
    			((H2F) this.getDataGroup().getItem(is,p<0?2:1,7,run).getData(i+il-1).get(0)).fill(Math.abs(p),uvw[i],z1);		  
    			((H2F) this.getDataGroup().getItem(is,p<0?4:3,7,run).getData(i+il-1).get(0)).fill(Math.abs(p),uvw[i],z2);		  
    			((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(i+il-1).get(0)).fill(-x,y,z1);
    			((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(i+il-1).get(0)).fill(-x,y,z2);
    		}
    	} 	
    }
    
    public void fillPATH(int is, int il, int run, float ec, float[] ep, float w, float v12mag, float v13mag, float v23mag, float cz, boolean fid[]) {
    	
    	boolean fidc = !(fid[0] || fid[1] || fid[2]);
    	
    	int il3=il/3;       
    	
    	for (int i=0; i<4; i++) {
    		boolean t = i==0 ? fidc:!fid[i-1];
    		((H2F) this.getDataGroup().getItem(i,il3,9,run).getData(is- 1).get(0)).fill(i==0?(t?ec/mipc[il3]:-10):(t?ep[i-1]/mipp[il3]:-10),v12mag);
    		((H2F) this.getDataGroup().getItem(i,il3,9,run).getData(is+ 5).get(0)).fill(i==0?(t?ec/mipc[il3]:-10):(t?ep[i-1]/mipp[il3]:-10),v13mag);
    		((H2F) this.getDataGroup().getItem(i,il3,9,run).getData(is+11).get(0)).fill(i==0?(t?ec/mipc[il3]:-10):(t?ep[i-1]/mipp[il3]:-10),v23mag);
    	}
     
        ((H2F) this.getDataGroup().getItem(0,il3,10,run).getData(is- 1).get(0)).fill(fidc?ec:0,w);
        ((H2F) this.getDataGroup().getItem(0,il3,10,run).getData(is+ 5).get(0)).fill(v12mag,w);    	
        ((H2F) this.getDataGroup().getItem(0,il3,10,run).getData(is+11).get(0)).fill(v13mag,w);    	
        ((H2F) this.getDataGroup().getItem(0,il3,10,run).getData(is+17).get(0)).fill(v23mag,w);    	
    }    
    
    public void getECALTDC(int run, int eis, DataEvent event) {
    	
    	etdcs.clear(); 
    	
        if(event.hasBank("ECAL::tdc")==true){
            DataBank  bank = event.getBank("ECAL::tdc");
            int rows = bank.rows();
            for(int i = 0; i < rows; i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  ip = bank.getShort("component",i);               
                float tdcd    = bank.getInt("TDC",i)*0.02345f;
                if(eis==is && tdcd>0) {
                    if(!etdcs.hasItem(is,il,ip)) etdcs.add(new ArrayList<Float>(),is,il,ip);
                        etdcs.getItem(is,il,ip).add(tdcd);       
                }
            }
        }    	
    }
    
    public void getECALADC(int run, int eis, DataEvent event) {
    	
   	 eadcs.clear(); 
    	
     if(event.hasBank("ECAL::adc")==true){
         DataBank  bank = event.getBank("ECAL::adc");
         int rows = bank.rows();
         for(int i = 0; i < rows; i++){
             int  is = bank.getByte("sector",i);            
             int  il = bank.getByte("layer",i);
             int  ip = bank.getShort("component",i);
             int adc = Math.abs(bank.getInt("ADC",i));
             float t = bank.getFloat("time",i);               
             int ped = bank.getShort("ped", i); 
             
             if(eis==is && adc>0) {
            	 if(!eadcs.hasItem(is,il,ip)) adcs.add(new ArrayList<Float>(),is,il,ip);
                 	 eadcs.getItem(is,il,ip).add((float) adc); 
             }
         }
     }             
    }
    
    public void getFTOFTDC(int run, int eis, DataEvent event) {
    	
    	tdcs.clear(); ltpmt.clear() ; 
        
        if(event.hasBank("FTOF::tdc")){
            DataBank  bank = event.getBank("FTOF::tdc");
            int rows = bank.rows();           
            for(int i = 0; i < rows; i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  lr = bank.getByte("order",i);                       
                int  ip = bank.getShort("component",i);
                float tdcd = bank.getInt("TDC",i)*0.02345f;  
                
                if(eis==is && tdcd>0) {
                	if(!tdcs.hasItem(is,il,lr-2,ip)) tdcs.add(new ArrayList<Float>(),is,il,lr-2,ip);
                    	tdcs.getItem(is,il,lr-2,ip).add(tdcd); 
                    	if (!ltpmt.hasItem(is,il,ip)) {
                    		ltpmt.add(new ArrayList<Integer>(),is,il,ip);
                    		ltpmt.getItem(is,il,ip).add(ip);
                    } 
                }
            }
        }        
    }
    
    public void getFTOFADC(int run, int eis, DataEvent event) { 
    	
    	 adcs.clear(); lapmt.clear();
    	 
    	 if(event.hasBank("FTOF::adc")){
            DataBank  bank = event.getBank("FTOF::adc");
            int rows = bank.rows();
            for(int i = 0; i < rows; i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  lr = bank.getByte("order",i);
                int  ip = bank.getShort("component",i);
                int adc = bank.getInt("ADC",i);
                float t = bank.getFloat("time",i);               
                int ped = bank.getShort("ped", i);
                
                if(eis==is && adc>100) {
                if(!adcs.hasItem(is,il,lr,ip)) adcs.add(new ArrayList<Float>(),is,il,lr,ip);
                    adcs.getItem(is,il,lr,ip).add((float) adc); 
                    if (!lapmt.hasItem(is,il,ip)) {
           	             lapmt.add(new ArrayList<Integer>(),is,il,ip);
                         lapmt.getItem(is,il,ip).add(ip);              
                    }
                }                   
            }
        }   	
    }
    
    public List<Integer> getFTOFBAR(float thr) {
 	   
        IndexGenerator ig = new IndexGenerator();
        List<Integer> bars = new ArrayList<Integer>(); 
        float gm=0;
        
        for (Map.Entry<Long,List<Integer>>  entry : lapmt.getMap().entrySet()){
            long hash = entry.getKey();
            int is = ig.getIndex(hash, 0);
            int il = ig.getIndex(hash, 1);
            int ip = ig.getIndex(hash, 2);
                   
         	if(adcs.hasItem(is,il,0,ip) && adcs.hasItem(is,il,1,ip)) {
         		gm = (float) Math.sqrt(adcs.getItem(is,il,0,ip).get(0)*adcs.getItem(is,il,1,ip).get(0));
         	}
         	if(gm>thr && il==2) bars.add(ip);
        }       
        return bars;       
    }    
    
    public void fillBAR(int run, DataEvent de) {
    	 if(!(de.hasBank("ECAL::adc") && de.hasBank("FTOF::adc"))) return;
    	 DataBank ecal = de.getBank("ECAL::adc");  
    	 DataBank ftof = de.getBank("FTOF::adc"); 
    	 for (int i=0; i<ftof.rows(); i++ ) {
             int  fis = ftof.getByte("sector", i);
             int  fil = ftof.getByte("layer", i); 
             int  fip = ftof.getShort("component", i);
             int  fad = ftof.getInt("ADC",i);
             if(fil==2) {
            	 for (int j=0; j<ecal.rows(); j++) {
            		 int  eis = ecal.getByte("sector", j);
            		 int  eil = ecal.getByte("layer", j); 
            		 int  eip = ecal.getShort("component", j); 
                     int  ead = ecal.getInt("ADC",j);
                     if(eis==fis) {
            			 ((H2F) this.getDataGroup().getItem(eis,1,15,run).getData(eil-1).get(0)).fill(fip,eip);  
           		 }
            	 }
             }
    	 }
    	
    }
    
    private void updateUVW(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int    pc = getActivePC();
        int    is = getActiveSector(); 
        
        c.clear();
        c.divide(3, 3);

        for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            c.cd(3*i+j); c.getPad(3*i+j).getAxisY().setLog(false); 
            c.draw(tl.fitData.getItem(is,i+10*pc*(j+1),0,getRunNumber()).getHist());
            c.draw(tl.fitData.getItem(is,i+10*pc*(j+1),0,getRunNumber()).getGraph(),"same");
        }
        }
        
    }
    
    public void updateFITS(int index) {
    	switch (index) {
    	case  2: updateFITS2(index); break;
    	case 15: updateFITS15(index);
    	}
    }
    
    public void updateFITS2(int index) {
       
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int    pc = getActivePC();
        int    is = getActiveSector(); 
        int     i = getActiveLayer();
        int     j = getActiveView();
        
        if (pc==2) pc=1;
        
        int    np = npmt[i*3+j];
        
        c.clear();
        if (np==36) c.divide(6,6); if (np>36) c.divide(8,9);
        
        for (int ip=0; ip<np ; ip++) {
            c.cd(ip); c.getPad(ip).getAxisY().setLog(false);
            c.draw(tl.fitData.getItem(is,i+10*(pc+1)*(pc+1)*(j+1),ip+1,getRunNumber()).getHist());
            c.draw(tl.fitData.getItem(is,i+10*(pc+1)*(pc+1)*(j+1),ip+1,getRunNumber()).getGraph(),"same");
       }
    }
    
    public void updateFITS15(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
    	GraphErrors g = null;
        
        int    is = getActiveSector();
        
        c.divide(3, 5);
        
        int[] np = {0,3,6,9,10,11,12,13,14};      
        
       	for (int n : np) {
       		if(Fits.hasItem(is,n,0,getRunNumber())) {
       			g=Fits.getItem(is,n,0,getRunNumber()).getGraph(); 
       			g.setMarkerColor(1); g.setLineColor(1);g.setMarkerSize(3); g.getFunction().setOptStat("110"); g.getFunction().setLineWidth(2);
       			H2F h = (H2F) this.getDataGroup().getItem(is,1,15,getRunNumber()).getData(n).get(0);
       			if(g.getDataSize(0)>0) {       				 
       				c.getPad(n).getAxisX().setRange(h.getXAxis().min(),h.getXAxis().max());
       				c.getPad(n).getAxisY().setRange(h.getYAxis().min(),h.getYAxis().max());
       				c.cd(n); c.draw(g,"same");
       			}
       		}       	
       	}            	
    }
    
    @Override
    public void plotEvent(DataEvent de) {
    	analyze();
    }

    public void analyze() {    
    	System.out.println(getDetectorName()+".analyze() ");
        fitGraphs(1,7,0,3,0,(dropSummary)?0:3);
        analyzeGraphs();
        if(!isAnalyzeDone) createTimeLineHistos();
        fillTimeLineHisto();
        System.out.println("Finished");
        isAnalyzeDone = true;
    }
    
    public void fitGraphs(int is1, int is2, int id1, int id2, int il1, int il2) {
    	
    	H2F h2=null, h2a=null, h2b=null; FitData fd=null;       
        int ipc=0,iipc=0, run=getRunNumber();
        double min=1,max=20,mip=10;
        System.out.println("Analyzing run "+run);
        for (int is=is1; is<is2; is++) {            
            for (int id=id1; id<id2; id++) {
            	for (int pc=0; pc<2; pc++) {
                    if(pc==0) {min = fitLimc[id]; max = fitLimc[id+3]; mip=mipc[id]; ipc=2;}
                    if(pc==1) {min = fitLimp[id]; max = fitLimp[id+3]; mip=mipp[id]; ipc=1;}  
                    h2 = CombineH2F((H2F) this.getDataGroup().getItem(is,ipc,0,run).getData(id*3+0).get(0),  //U
		                            (H2F) this.getDataGroup().getItem(is,ipc,0,run).getData(id*3+1).get(0),  //V        
		                            (H2F) this.getDataGroup().getItem(is,ipc,0,run).getData(id*3+2).get(0)); //W    
                    if(TLname=="UVW")       {h2a = rebinY(h2,npmt[id*3+0],npmt[id*3+1],npmt[id*3+2]);}  
                    if(TLname=="FADC Slot") {h2a = rebinY(h2,TimeSlice.get(TLname));}      
                    if(TLname=="HV Slot")   {h2a = rebinY(h2,TimeSlice.get(TLname));} 
                    int nb = h2a.getYAxis().getNBins();
            	    tl.setNYbins(id, nb);
                    for (int il=0; il<((pc==0)?1:nb); il++) {
                     	fd = fitEngine(h2a.sliceY(il),0,min,max,min,max); //Sector+UVW slices
                     	if(TLname=="UVW") fd.hist.getAttributes().setTitleX(((H2F) this.getDataGroup().getItem(is,ipc,0,run).getData(id*3+il).get(0)).getTitleX()); 
                	    tl.fitData.add(fd,is,id+10*pc*(il+1),0,run); 
                    }
                    for (int il=il1; il<il2; il++) { //PMT slices               	
                        h2b = (H2F) this.getDataGroup().getItem(is,ipc,0,run).getData(id*3+il).get(0);
                    	for (int i=0; i<npmt[id*3+il]; i++) tl.fitData.add(fitEngine(h2b.sliceY(i),0,min,max,min,max),is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); //PMT slices
            		    fitStore(is, id, il, pc, run, mip);
                    }   
               }
            }
        }
    }
    
    public void analyzeGraphs() {
    	
    	ParallelSliceFitter fitter1, fitter2;
        GraphErrors g = new GraphErrors();
        
        int run=getRunNumber(); 
        
        int[] np = {0,3,6,9,10,11,12,13,14};      
        int[] xy = {0,0,0,0, 0, 0, 0, 0, 0};
        int[] fl = {35,56,56,52,15,15,52,15,15};
        
        for (int is=1; is<7; is++) {
        	for (int i=0; i<np.length; i++ ) {
        		H2F h = (H2F) this.getDataGroup().getItem(is,1,15,run).getData(np[i]).get(0);
        		fitter1 = new ParallelSliceFitter(h);
        		fitter1.setRange(0,62); 
        		if(xy[i]==0) {fitter1.fitSlicesX(); g = fitter1.getMeanSlices(); g.setTitleX(h.getTitleX()); g.setTitleY(h.getTitleY());}
        		if(xy[i]==1) {fitter1.fitSlicesY(); g = fitter1.getMeanSlices(); g.setTitleX(h.getTitleY()); g.setTitleY(h.getTitleX());}
        		Fits.add(fitEngine(g,6,0,fl[i]), is,np[i],0,run);        		
        	}        	       	
        }
    }
    
    public void fitStore(int is, int id, int il, int pc, int run, double nrm) {
    	FitData fd=null;  
    	int np = npmt[id*3+il];
        double[]      x = new double[np]; double[]  ymean = new double[np]; double[] yrms = new double[np];
        double[]     xe = new double[np]; double[] ymeane = new double[np]; double[]   ye = new double[np]; 
        double[]  yMean = new double[np]; double[] yMeanc = new double[np]; double[]  ymeanc = new double[np];
        H1F h1 = new H1F("VAR1-"+is+"-"+id+" "+il,20,0.7,1.7);
        H1F h5 = new H1F("VAR5-"+is+"-"+id+" "+il,20,0.7,1.7);
        H1F h6 = new H1F("VAR6-"+is+"-"+id+" "+il,20,0.7,1.7);
        H1F h7 = new H1F("VAR7-"+is+"-"+id+" "+il,20,0.7,1.7);
        for (int i=0; i<np; i++) {
            x[i] = i+1; xe[i]=0; ye[i]=0; yrms[i]=0; 
            fd = tl.fitData.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); 
            fd.graph.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            fd.hist.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            double mean = fd.mean;       
            if(mean>0) yrms[i] = fd.sigma/mean; 
            	     yMeanc[i] = (fd.getMean()/nrm)*((pc==1)?shift.getDoubleValue("shift", is, 3*id+il+1, i+1):0); //torus correction
                     ymeanc[i] =         (mean/nrm)*((pc==1)?shift.getDoubleValue("shift", is, 3*id+il+1, i+1):0); //torus correction
                      yMean[i] = fd.getMean()/nrm;
                      ymean[i] =         mean/nrm;
                     ymeane[i] = fd.meane/nrm;
                     if(varcut(id,il,i)&&fd.integral>15) {h1.fill(ymean[i]); h5.fill(yMean[i]); h6.fill(ymeanc[i]); h7.fill(yMeanc[i]);}
        }
//       fd = fitEngine(h5,0.5,1.5,0.5,1.5); tl.fitData.add(fd,is,id+100*pc*(il+1),0,run); 
//       fd = fitEngine(h7,0.5,1.5,0.5,1.5); tl.fitData.add(fd,is,id+100*pc*(il+1),0,run); 
        GraphErrors  mean = new GraphErrors("MIP-"+is+"-"+id+" "+il,x,ymean,xe,ymeane);                   
        GraphErrors  Mean = new GraphErrors("MIP-"+is+"-"+id+" "+il,x,yMean,xe,ymeane);                   
        GraphErrors meanc = new GraphErrors("MIP-"+is+"-"+id+" "+il,x,ymeanc,xe,ymeane);                   
        GraphErrors Meanc = new GraphErrors("MIP-"+is+"-"+id+" "+il,x,yMeanc,xe,ymeane);                   
        GraphErrors   rms = new GraphErrors("MIP-"+is+"-"+id+" "+il,x,yrms,xe,ye);                  
        FitSummary.add(mean,  is,id+10*(pc+1)*(pc+1)*(il+1),1,run); VarSummary.add(h1,  is,id+10*(pc+1)*(pc+1)*(il+1),1,run);
        FitSummary.add(rms,   is,id+10*(pc+1)*(pc+1)*(il+1),2,run);                    
        FitSummary.add(Mean,  is,id+10*(pc+1)*(pc+1)*(il+1),5,run); VarSummary.add(h5,  is,id+10*(pc+1)*(pc+1)*(il+1),5,run);       	        
        FitSummary.add(meanc, is,id+10*(pc+1)*(pc+1)*(il+1),6,run); VarSummary.add(h6,  is,id+10*(pc+1)*(pc+1)*(il+1),6,run);       	        
        FitSummary.add(Meanc, is,id+10*(pc+1)*(pc+1)*(il+1),7,run); VarSummary.add(h7,  is,id+10*(pc+1)*(pc+1)*(il+1),7,run);
        
    }
    
    public boolean varcut(int id, int il, int ip) {
    	if (id==0 && il==0 && ip<5 && ip>67) return false;
    	if (id==0 && il>0  && ip<5 && ip>61) return false;
    	if (id>0  && ip<5  && ip>35)         return false;
    	return true;
    }
    
    public void plotVarSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = getActivePC();
        int            n = 0;
        
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        for (int il=0; il<3; il++) {           	
            H1F plot1 = VarSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),getActiveRDIF()==0?1:6,getRunNumber());
            H1F plot2 = VarSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),getActiveRDIF()==0?5:7,getRunNumber());

            plot1.setFillColor(2) ; plot2.setFillColor(38); plot1.getAttributes().setOptStat("1100");plot2.getAttributes().setOptStat("1100");           
           
            c.cd(n); c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16); c.getPad(n).setStatBoxFontSize(18);
            plot2.getAttributes().setTitleX("S"+is+" "+det[id]+" "+v[il].toUpperCase()+" MEAN/MIP");
            n++; c.draw(plot2);                   
        }
        }
        }        
    } 
    
    public void plotMeanSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = getActivePC();
        int            n = 0;
        
        double ymin=0.99f, ymax=1.01f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
         
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        for (int il=0; il<3; il++) {           	
            F1D f1 = new F1D("p0","[a]",0.,npmt[id*3+il]); f1.setParameter(0,1);
            GraphErrors plot1 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),getActiveRDIF()==0?1:6,getRunNumber());
            GraphErrors plot2 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),getActiveRDIF()==0?5:7,getRunNumber());
            plot1.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) plot1.getAttributes().setTitleY("MEAN / MIP");
            plot1.getAttributes().setTitleX("S"+is+" "+det[id]+" "+v[il].toUpperCase()+" PMT");
            n++; c.draw(plot1); c.draw(plot2,"same");
            f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        }
        }
        }        
    }        
    
    public void plotRmsSummary(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = getActivePC();
        int            n = 0;
        
        double ymin=0.3f, ymax=1.0f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
         
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        for (int il=0; il<3; il++) {           	
            GraphErrors plot1 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),2,getRunNumber());
            plot1.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) plot1.getAttributes().setTitleY("RMS / MEAN");
            plot1.getAttributes().setTitleX("S"+is+" "+det[id]+" "+v[il].toUpperCase()+" PMT");
            n++; c.draw(plot1);
        }
        }
        }          

    }
    
    public void plotRmsHWSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = getActivePC();
        int            n = 0;
        
        List<DataLine> lines = new ArrayList<DataLine>();
        
        Boolean t = TLname!="UVW";
        
        c.clear(); c.divide(3, 6);
                
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        	GraphErrors hwplot1 = new GraphErrors();
        	int m=0; lines.clear();
            for (int il=0; il<3; il++) {           	
                GraphErrors plot1 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),2,getRunNumber());
                 for (int ip=0; ip<npmt[id*3+il]; ip++) {m++;
        	        hwplot1.addPoint(m, plot1.getDataY(ip), plot1.getDataEX(ip), plot1.getDataEY(ip));
        	        if(Math.floorMod(t?m:ip, t?TimeSlice.get(TLname):npmt[id*3+il])==(t?1:0)) {
        	    	    DataLine line = new DataLine(m,0.0,m,0.5) ; line.setLineColor(1); line.setLineWidth(1); 
        	    	    lines.add(line);
        	        }
                }
        }
        c.cd(n); c.getPad(n).getAxisY().setRange((pc==0)?0.1:0.18, (pc==0)?0.3:0.5); 
        c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
        if(n==0||n==3||n==6||n==9||n==12||n==15) hwplot1.getAttributes().setTitleY("RMS / MEAN");
        hwplot1.getAttributes().setTitleX("S"+is+" "+det[id]+" PMT");
        n++; c.draw(hwplot1); for(DataLine line: lines) c.draw(line);
        F1D f1 = new F1D("p0","[a]",0.,m); f1.setParameter(0,1);
        f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        }
        }        
    }   
    
    public void plotMeanHWSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = getActivePC();
        int            n = 0;
        
        List<DataLine> lines = new ArrayList<DataLine>();
        
        Boolean t = TLname!="UVW";
        double ymin=0.99f, ymax=1.01f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;
        
        c.clear(); c.divide(3, 6);
                
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        	GraphErrors hwplot1 = new GraphErrors();
        	GraphErrors hwplot2 = new GraphErrors();
        	int m=0; lines.clear();
            for (int il=0; il<3; il++) {           	
                GraphErrors plot1 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),1,getRunNumber()); 
                GraphErrors plot2 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),5,getRunNumber());
                for (int ip=0; ip<npmt[id*3+il]; ip++) {m++;
        	        hwplot1.addPoint(m, plot1.getDataY(ip), plot1.getDataEX(ip), plot1.getDataEY(ip));
        	        hwplot2.addPoint(m, plot2.getDataY(ip), plot2.getDataEX(ip), plot2.getDataEY(ip));
        	        if(Math.floorMod(t?m:ip, t?TimeSlice.get(TLname):npmt[id*3+il])==(t?1:0)) {
        	        	int mm = (TLname=="HV Slot"&&id==2)?m+12:m;
        	    	    DataLine line = new DataLine(mm,ymin,mm,ymax) ; line.setLineColor(1); line.setLineWidth(1); 
        	    	    lines.add(line);
        	        }
                }
            }
            hwplot1.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==3||n==6||n==9||n==12||n==15) hwplot1.getAttributes().setTitleY("MEAN / MIP");
            hwplot1.getAttributes().setTitleX("S"+is+" "+det[id]+" PMT");
            n++; c.draw(hwplot1); c.draw(hwplot2,"same"); for(DataLine line: lines) c.draw(line);
            F1D f1 = new F1D("p0","[a]",0.,m); f1.setParameter(0,1);
            f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        }
        }        
    }   

    public void getMeanRDIF() {
        RDIFmap.clear();
        for (int pc=0; pc<2; pc++) {
        	RDIFmap.add(new IndexedList<List<Float>>(4));
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        for (int il=0; il<3; il++) {           	
            GraphErrors plot1 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),5,runlist.get(0));
            GraphErrors plot2 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),5,runlist.get(1));
            int m=0;
            for (int ip=0; ip<npmt[id*3+il]; ip++) {m++;
                double x1=plot1.getDataY(ip) ; double x2=plot2.getDataY(ip); double x1e=plot1.getDataEY(ip); double x2e=plot2.getDataEY(ip);
        		double  y1 = x1>0&&x2>0 ? x1/x2:1.0;
        		double y1e = x1>0&&x2>0 ? y1*Math.sqrt(Math.pow(x1e/x1,2)+Math.pow(x2e/x2,2)):ip>5?0.09:0;
                RDIFmap.get(pc).add(new ArrayList<Float>(), is,id,il,ip); 
                RDIFmap.get(pc).getItem(is,id,il,ip).add((float)m);
                RDIFmap.get(pc).getItem(is,id,il,ip).add((float)(y1e>0&&y1e<(ip>5?1:0.1)?y1:1));
                RDIFmap.get(pc).getItem(is,id,il,ip).add((float)plot1.getDataEX(ip));
                RDIFmap.get(pc).getItem(is,id,il,ip).add((float)(y1e>0&&y1e<(ip>5?1:0.1)?y1e:0));
            }            
        }
        }
        }    
        }
    	
    }
    
    public void plotMeanRDIF(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = getActivePC();
        int            n = 0;
        
        double ymin=0.99f, ymax=1.01f;
        ymin=ymin*lMin/250; ymax=ymax*lMax/250;   
        
        c.clear(); c.divide(9, 6);
        
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        for (int il=0; il<3; il++) {           	
        	GraphErrors hwplot1 = new GraphErrors();
            F1D f1 = new F1D("p0","[a]",0.,npmt[id*3+il]); f1.setParameter(0,1);
            for (int ip=0; ip<npmt[id*3+il]; ip++) {
            	hwplot1.addPoint(RDIFmap.get(pc).getItem(is,id,il,ip).get(0),
                                 RDIFmap.get(pc).getItem(is,id,il,ip).get(1),
                                 RDIFmap.get(pc).getItem(is,id,il,ip).get(2),
                                 RDIFmap.get(pc).getItem(is,id,il,ip).get(3));
            }            
            hwplot1.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin, ymax); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) hwplot1.getAttributes().setTitleY("#Delta MEAN / MEAN");
            hwplot1.getAttributes().setTitleX("S"+is+" "+det[id]+" "+v[il].toUpperCase()+" PMT");
            n++; c.draw(hwplot1); 
            f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        }
        }
        }        
    } 
    
    public void plotXY(int index) {        
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));        
        int          run = getRunNumber();
        
        c.clear(); c.divide(4,3);
        
    	for (int i=0; i<3; i++) {
            for (int j=0; j<4; j++) {int ind1=(j==3)?2:1;int ind2=(j==3)?i:3*i+j;
        	    H2F h1 = (H2F) this.getDataGroup().getItem(0,ind1,5,run).getData(ind2).get(0);   
                H2F h2 = (H2F) this.getDataGroup().getItem(1,ind1,5,run).getData(ind2).get(0); 
                H2F h3 = (H2F) this.getDataGroup().getItem(2,ind1,5,run).getData(ind2).get(0); 
                h3 = h2.divide(h2, h1); h3.setTitle(h1.getName());
                c.cd(4*i+j); c.getPad(4*i+j).getAxisZ().setLog(false); c.getPad(4*i+j).getAxisZ().setRange(0., 2.);
                c.draw(h3);            
        	}
    	}
   
    }
    
    public void plotMIPZ(int index) {        
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));        
        int          run = getRunNumber();
        
        c.clear(); c.divide(3,3);
        
    	for (int i=1; i<3; i++) {
    		for (int il=0; il<9; il++) {
    			for (int is=1; is<7; is++) {
    				H2F h1 = (H2F) this.getDataGroup().getItem(is,i,  index,run).getData(il).get(0);   
    				H2F h2 = (H2F) this.getDataGroup().getItem(is,i+2,index,run).getData(il).get(0); 
    				this.getDataGroup().getItem(is,i+4,index,run).getData(il).clear();
    				this.getDataGroup().getItem(is,i+4,index,run).getData(il).add(h2.divide(h2, h1));
    				((H2F)this.getDataGroup().getItem(is,i+4,index,run).getData(il).get(0)).setTitle(h1.getName());           
    			}
    		}
    	}
    	
    	zMin=0f; zMax=0.01;
    	int pc = getActivePC()==0?6:5;
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),pc,index,getRunNumber()));    	
   
    } 
    
    public void plotAlignSummary(int index) {
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
      	List<GraphErrors> gl = new ArrayList<GraphErrors>();
      	
      	int[]  np = {0,3,6,9,10,11,12,13,14};
//      int[] ncd = {0,3,1, 4,2, 5,6, 9,7,10,8,11,12,15,13,16,14,17}; //p0 alternate with p1  
      	int[] ncd = {0,9,1,10,2,11,3,12,4,13,5,14, 6,15, 7,16, 8,17}; //p0 first, then p1  
      	int[] col = {1,3,2}; int[] ms = {4,3,3};
      	float[] p1bw = {379.254f,379.413f,379.254f,379.095f,379.254f,379.096f};
      	float[]  pcw = {4.529f,4.517f,4.518f,4.526f,4.524f,4.528f};        
      	String opt,tit=" "; int nn = 0;
      	
      	GraphErrors rat = new GraphErrors(); rat.setMarkerColor(4); rat.setMarkerSize(4);
      	for (int iis=0; iis<6; iis++) rat.addPoint(iis+1, p1bw[iis]/pcw[iis]/62, 0, 0);
      	
      	c.clear(); c.divide(3, 6);
      	
  		for (int ip=0; ip<np.length; ip++) {int inp = np[ip];
  			for (int in=0; in<2; in++) {gl.clear(); float ysum=0; int nsum=0;float ymn=1000,ymx=-1000;
  			for (int ir=0; ir<runlist.size(); ir++) { int run = runlist.get(ir);
  				if(ir==0) tit =  ((H2F) this.getDataGroup().getItem(1,1,15,run).getData(inp).get(0)).getTitle();
 				gl.add(new GraphErrors()); gl.get(ir).setMarkerColor(col[ir]);gl.get(ir).setLineColor(col[ir]);gl.get(ir).setMarkerSize(ms[ir]);    
      	      	for (int is=1; is<7; is++)	{
      				float  y = (float) Fits.getItem(is,inp,0,run).graph.getFunction().parameter(in).value();  
      				float ye = (float) Fits.getItem(is,inp,0,run).graph.getFunction().parameter(in).error();
      				ysum+=y; ymx=y>ymx?y:ymx; ymn=y<ymn?y:ymn;    				
      				gl.get(ir).addPoint(is+ir*0.1, y, 0, ye); 
      				nsum++;
      			}
      		}  			
  			float ymid = ysum/nsum;
  	      	c.cd(ncd[nn]); c.getPad(ncd[nn]).getAxisY().setRange(ymid-(ymx-ymn)*1.2,ymid+(ymx-ymn)*1.2); 
  	      	gl.get(0).setTitleY("P"+in); gl.get(0).setTitle("BLACK: "+tit+" GREEN: rga-fall2018 RED: default");
  	      	opt=" "; for(GraphErrors g : gl) {c.draw(g,opt); opt="same";} if(ncd[nn]==9) c.draw(rat,opt); nn++;
  			}
      	}
    }
    
    public void plotMIP(int index) {
    	int pc = getActivePC()==0?2:1;
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),pc,index,getRunNumber()));
    }
        
    public void plotPIDSummary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));
    }
    
    public void plotPathSummary(int index, int flag) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(flag,getActiveLayer(),index,getRunNumber()));
    }
    
    public void plotUVW(int index) {
        int run = getRunNumber();
        drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index,run));	    
    } 

/*   TIMELINES */
    
    @Override
    public void createTimeLineHistos() {   
    	System.out.println("Initializing "+TLname+" timeline"); 
    	runIndex = 0;
    	tl.createTimeLineHisto(10,"PCAL Cluster Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto(20,"ECIN Cluster Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto(30,"ECOU Cluster Mean/MIP","Sector",TLmax,6,1,7);    	
    	tl.createTimeLineHisto( 1,"PCAL U Peak Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 2,"PCAL V Peak Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 3,"PCAL W Peak Mean/MIP","Sector",TLmax,6,1,7);    	
    	tl.createTimeLineHisto( 4,"ECIN U Peak Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 5,"ECIN V Peak Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 6,"ECIN W Peak Mean/MIP","Sector",TLmax,6,1,7);    	
    	tl.createTimeLineHisto( 7,"ECOU U Peak Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 8,"ECOU V Peak Mean/MIP","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto( 9,"ECOU W Peak Mean/MIP","Sector",TLmax,6,1,7);
    }
    
    public void fillTimeLineHisto() {    	
    	//clusters
		for (int is=1; is<7; is++) {
		  for (int il=0; il<3; il++) {
			  float  y = (float) tl.fitData.getItem(is,il,0,getRunNumber()).mean/mipc[il];
			  float ye = (float) tl.fitData.getItem(is,il,0,getRunNumber()).meane/mipc[il];
			  ((H2F)tl.Timeline.getItem((il+1)*10,0)).fill(runIndex,is,y);
			  ((H2F)tl.Timeline.getItem((il+1)*10,1)).fill(runIndex,is,ye);			  
		  }
		}
	
		//peaks
		for (int is=1; is<7; is++) {
		  for (int il=0; il<3; il++) {	
	        for (int iv=0; iv<3; iv++) {
			    float  y = (float) tl.fitData.getItem(is,il+10*(iv+1),0,getRunNumber()).mean/mipp[il];
			    float ye = (float) tl.fitData.getItem(is,il+10*(iv+1),0,getRunNumber()).meane/mipp[il];
			    ((H2F)tl.Timeline.getItem(3*il+iv+1,0)).fill(runIndex,is,y);	
			    ((H2F)tl.Timeline.getItem(3*il+iv+1,1)).fill(runIndex,is,ye);	
	        }
		  }
		}
		runIndex++;
    }
        
    public void saveTimelines() {
    	System.out.println("ECmip: Saving timelines");
    	saveTimeLine(10,0,0,"PCALmip","MIP");
    	saveTimeLine(20,1,0,"ECINmip","MIP");
    	saveTimeLine(30,2,0,"ECOUmip","MIP");
    	saveTimeLine(1,10,0,"PCALmipU","MIP");
    	saveTimeLine(2,20,0,"PCALmipV","MIP");
    	saveTimeLine(3,30,0,"PCALmipW","MIP");
    	saveTimeLine(4,11,0,"ECINmipU","MIP");
    	saveTimeLine(5,21,0,"ECINmipV","MIP");
    	saveTimeLine(6,31,0,"ECINmipW","MIP");
    	saveTimeLine(7,12,0,"ECOUmipU","MIP");
    	saveTimeLine(8,22,0,"ECOUmipV","MIP");
    	saveTimeLine(9,32,0,"ECOUmipW","MIP");

    }
    
    public void plotTimeLines(int index) {        
    	if(TLflag) {plotTimeLineSectors(index); } else {
    	if (getActivePC()==0) {plotClusterTimeLines(index);} else {plotPeakTimeLines(index);}}
    }
    
    public void plotClusterTimeLines(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); 
        int           is = getActiveSector(); 
        
        FitData       fd = null;
       
    	DataLine line1 = new DataLine(0,is,  runIndex+1,is);                   line1.setLineColor(5);
    	DataLine line2 = new DataLine(0,is+1,runIndex+1,is+1);                 line2.setLineColor(5);
    	DataLine line3 = new DataLine(runIndexSlider,1,  runIndexSlider,  7);  line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,runIndexSlider+1,7);  line4.setLineColor(5);

        c.clear(); c.divide(3, 3); 

        for (int il=0; il<3; il++) { int i3=il*3; 
            double min=0.99,max=1.01; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(i3); c.getPad(i3).setAxisRange(0,runIndex,1,7); c.getPad(i3).setTitleFontSize(18); c.getPad(i3).getAxisZ().setRange(min,max);
    		c.draw((H2F)tl.Timeline.getItem(10*(il+1),0));c.draw(line1);c.draw(line2);c.draw(line3);c.draw(line4);
             
    		c.cd(i3+1); c.getPad(i3+1).setAxisRange(-0.5,runIndex,min,max); c.getPad(i3+1).setTitleFontSize(18);
    		drawTimeLine(c,is,10*(il+1),1f,"Sector "+is+" Mean/MIP" );
    		
    		fd = tl.fitData.getItem(is,il,0,getRunNumber()); H1F h = fd.getHist().histClone("");
    		
    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(0.,h.getXaxis().max(),0.,fd.getGraph().getMax()*1.1);  
            h.getAttributes().setOptStat("1000100");
            DataLine line6 = new DataLine(mipc[il],-50,mipc[il],fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);            
            c.draw(h); c.draw(fd.getGraph(),"same");  c.draw(line6);            
        }
    }
     
    public void plotPeakTimeLines(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));                
        int           is = getActiveSector();
        int           iv = getActiveView();
        
        FitData       fd = null;
        
    	String  v[] = {" U "," V "," W "};
        
    	DataLine line1 = new DataLine(0,is,  runIndex+1,is);                     line1.setLineColor(5);
    	DataLine line2 = new DataLine(0,is+1,runIndex+1,is+1);                   line2.setLineColor(5);
    	DataLine line3 = new DataLine(runIndexSlider,  1,  runIndexSlider,  7);  line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,  runIndexSlider+1,7);  line4.setLineColor(5);
    	
        c.clear(); c.divide(3, 3); 

        for (int il=0; il<3; il++) {int i3=il*3;
            double min=0.99f,max=1.01f; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(i3); c.getPad(i3).setAxisRange(0,runIndex,1,7); c.getPad(i3).setTitleFontSize(18); c.getPad(i3).getAxisZ().setRange(min,max);
    		c.draw((H2F)tl.Timeline.getItem(3*il+iv+1,0));c.draw(line1);c.draw(line2);c.draw(line3);c.draw(line4);
             
    		c.cd(i3+1); c.getPad(i3+1).setAxisRange(-0.5,runIndex,min,max); c.getPad(i3+1).setTitleFontSize(18);
    		drawTimeLine(c,is,3*il+iv+1,1f,"Sector "+is+v[iv]+" Mean/MIP" );
    		
    		fd = tl.fitData.getItem(is,il+10*(iv+1),0,getRunNumber()); H1F h = fd.getHist().histClone("");
    		
    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(0.,h.getXaxis().max(),0.,fd.getGraph().getMax()*1.1);  
            h.getAttributes().setOptStat("1000100");
            DataLine line6 = new DataLine(mipp[il],-50,mipp[il],fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);
            c.draw(h); c.draw(fd.getGraph(),"same");  c.draw(line6);       
        }    	
    }
    
    public void plotTimeLineSectors(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int pc = getActivePC();
        int iv = getActiveView();
        int il = getActiveLayer();
       	String  v[] = {" U "," V "," W "};
       	String  l[] = {" PCAL "," ECIN "," ECOU "};
        
        c.clear(); c.divide(3, 2);
    	for (int is=1; is<7; is++) {
    		double min=0.99 ; double max=1.01; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(is-1); c.getPad(is-1).setAxisRange(-0.5,runIndex,min,max); c.getPad(is-1).setTitleFontSize(18);
    		drawTimeLine(c,is,(pc==0)?10*(il+1):3*il+iv+1,1f,"Sector "+is+((pc==0)?" Mean/MIP":l[il]+v[iv]+" Mean/MIP"));
    	}
    }
    
/* CALIBRATION FILES */    
    
    @Override
	public void writeFile(String table, int is1, int is2, int il1, int il2, int iv1, int iv2) {
		
    	if(!dumpFiles) return;
    	
		String path = "/Users/colesmith/CLAS12ANA/ECcalib/ccdb/";
		String line = new String();
		
		try { 
			File outputFile = new File(path+table+"_"+detcal[0]);
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("ECcalib.writefile("+table+")");

			for (int is=is1; is<is2; is++) {
				for (int il=il1; il<il2; il++ ) { //pcal,ecin,ecou
					for (int iv=iv1; iv<iv2; iv++) { //u,v,w
						for (int ip=0; ip<npmt[3*il+iv]; ip++) {
							switch (table) {
							case "gain":     line = getGAIN(is,il,iv,ip,runlist.get(il)); break;
							case "rdif":     line = getRDIF(is,il,iv,ip); break;
							case "rdifgain": line = getRDIFGAIN(is,il,iv,ip); 
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
    
	public String getGAIN(int is, int il, int iv, int ip, int run) {
		int pc = 1; 		
		if(tl.fitData.hasItem(is,il+10*(pc+1)*(pc+1)*(iv+1),ip+1,run)) {
			double     g = tl.fitData.getItem(is,il+10*(pc+1)*(pc+1)*(iv+1),ip+1,run).getMean()/mipp[il];
			double    ge = tl.fitData.getItem(is,il+10*(pc+1)*(pc+1)*(iv+1),ip+1,run).meane/mipp[il];
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				  +(gain.getDoubleValue("gain", is, 3*il+iv+1, ip+1)/(g<0.2?1.0:g))+" "
				  +ge;
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+"  "
			      + gain.getDoubleValue("gain", is, 3*il+iv+1, ip+1)
			      +" 0.0";			
		}
		
	}
	
	public String getRDIF(int is, int il, int iv, int ip) {	
		int pc = 1;
		if(RDIFmap.get(pc).hasItem(is,il,iv,ip) && il>0) {
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+RDIFmap.get(pc).getItem(is,il,iv,ip).get(1)+" "
				+RDIFmap.get(pc).getItem(is,il,iv,ip).get(3);
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+"  "
		        + "1.0"
		        +" 0.0";
		}	    
	}
	
	public String getRDIFGAIN(int is, int il, int iv, int ip) {	
		return is+" "+(3*il+iv+1)+" "+(ip+1)+" "
				+shift.getDoubleValue("shift", is, 3*il+iv+1, ip+1)*gain.getDoubleValue("gain",   is, 3*il+iv+1, ip+1)+" "
				+shift.getDoubleValue("shift", is, 3*il+iv+1, ip+1)*gain.getDoubleValue("gainErr",is, 3*il+iv+1, ip+1);    
	}

/*  
 * 
 *       
    public void createTimeLineGraph(int k, String tit, String xtit, String ytit) {    	
    	 String[] v = new String[]{"U","V","W"};
    	 for (int is=1; is<7; is++) {
    		for (int il=1; il<4; il++) {
    	       GraphErrors g = new GraphErrors(tit);
    	       g.setTitleX("Run Index"); g.setTitleY("Sector "+is+" "+v[il-1]+" "+ytit); g.setLineColor(1); g.setMarkerColor(1); g.setMarkerSize(3);
    	       Timeline.add(g,k,is);
    		}
    	}
    }
    
    public void createSectorGraphs() {
    	createSectorGraph(0,"PCAL MIP Cluster Energy","Sector","Mean/MIP");
    	createSectorGraph(1,"ECIN MIP Cluster Energy","Sector","Mean/MIP");
    	createSectorGraph(2,"ECOU MIP Cluster Energy","Sector","Mean/MIP");
    }
    
        
    public void createSectorGraph(int k, String tit, String xtit, String ytit) {  
    	
    	int ic=0;
    	for (int i=0; i<runlist.size(); i++) {
    		if (runlist.get(i)==getRunNumber()) {
            GraphErrors  g = new GraphErrors("MIP"); 
            ic=(i>8)?i-8:i+1;       
            g.setTitle(tit); g.setTitleX(xtit); g.setTitleY(ytit); g.setLineColor(ic); g.setMarkerColor(ic); g.setMarkerSize(5);  
            Timeline.add(g,k,i);
    		}
    	}
    }
    
    public void fillTimeLineGraph() {
    	
    	float mip[] = {30,30,48};
    	
    	for (int i=0; i<runlist.size(); i++) {
    		if (runlist.get(i)==getRunNumber()) {
    		for (int is=1; is<7; is++) {
    			double off = -runlist.size()*0.05+i*0.1;
    			for (int id=0; id<3; id++) {
    				float  y = (float) MipFits.getItem(is,3*id,0,getRunNumber()).mean/mip[id];
    				float ye = (float) MipFits.getItem(is,3*id,0,getRunNumber()).meane/mip[id];
    				((GraphErrors)Timeline.getItem(id,i)).addPoint(is+off, y, 0., ye);	
    			}	
    		}
    		}
    	}    	
    }    
    public void plotTimeLineGraph(int index) {
    	int[] xlab = {80, 140, 200, 260, 320, 380, 80, 140, 200, 260, 320, 380, 80, 140, 200, 260, 320, 380};
    	int[] ylab = {30, 30, 30, 30, 30, 30, 50, 50, 50, 50, 50, 50, 70, 70, 70, 70, 70, 70};

        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));                
        c.divide(3, 2); 
        
        F1D f1 = new F1D("p0","[a]",0.5,6.5); f1.setParameter(0,1);
               
        for (int id=0; id<3; id++) {        	
    		c.cd(id); c.getPad(id).setAxisRange(0.5, 6.5, 0.85, 1.15);c.getPad(id).setTitleFontSize(18);
            int ic=0;
          	for (int i=0; i<runlist.size(); i++) {
        		if (runlist.get(i)==getRunNumber()) {
            		c.draw((GraphErrors)Timeline.getItem(id,i),i==0?" ":"same");           		
        		    LatexText text = new LatexText("RUN "+runlist.get(i),xlab[i],ylab[i]);
        		    text.setFont("HanziPen TC");
        		    text.setFontSize(10);
        		    ic=(i>8)?i-8:i+1;
        		    text.setColor(ic);
        		    c.draw(text);
        		}
        	}
        	f1.setLineColor(3); f1.setLineWidth(3); c.draw(f1,"same");	
        }
        
    }
    
    */

/*    
    private void updateSummary() {
        
        DataGroup dg4 = this.getDataGroup().getItem(4,0,0);
        H2F h2;
        EmbeddedCanvas c = null;
        String id = det[getActiveLayer()];        
        c = this.summary.getCanvas("PCAL/ECTOT");
        c.divide(3,4);
        for (int is=1; is<7; is++) {
            h2 =dg4.getH2F("hi-"+id+"-path1-"+is);   
            c.cd(is-1); c.getPad(is-1).getAxisZ().setLog(true);       
            c.draw(h2);   
        }
        for (int is=1; is<7; is++) {
            h2 =dg4.getH2F("hi-"+id+"-path2-"+is);   
            c.cd(is-1+6); c.getPad(is-1+6).getAxisZ().setLog(true);       
            c.draw(h2);   
        }
        
        c.repaint();
    }
    
    private void updatePvsE() {
        
        DataGroup dg4 = this.getDataGroup().getItem(4,0,0);
        H2F h2;
        EmbeddedCanvas c = null; 
        c = this.summary.getCanvas("PvsE");
        c.divide(2,2);
        h2 =dg4.getH2F("hi-pcal-1");   
        c.cd(0); c.getPad(0).getAxisZ().setLog(true); c.draw(h2);
        h2 =dg4.getH2F("hi-ecali-1");   
        c.cd(1); c.getPad(1).getAxisZ().setLog(true); c.draw(h2);
        h2 =dg4.getH2F("hi-ecalo-1");   
        c.cd(2); c.getPad(2).getAxisZ().setLog(true); c.draw(h2);
        h2 =dg4.getH2F("hi-etot-1");   
        c.cd(3); c.getPad(3).getAxisZ().setLog(true); c.draw(h2);
        
        c.repaint();
        
    }
*/  
    public void getPixLengthMap(String filename) {   
        
    	System.out.println("ECmip.getPixLengthMap("+filename+")");
    	
        try{
            FileReader       file = new FileReader(filename);
            BufferedReader reader = new BufferedReader(file);
            int n = 0 ;
            while (n<1296) {
              String line = reader.readLine();
              String[] col = line.trim().split("\\s+"); 
              int i = Integer.parseInt(col[0]); 
              int j = Integer.parseInt(col[1]);
              int k = Integer.parseInt(col[2]);
              float d = Float.parseFloat(col[3]);
              PixLength.add(d,i,j,k);
              n++;
            }    
            reader.close();
            file.close();
         }  
         
         catch(FileNotFoundException ex) {
            ex.printStackTrace();            
         }     
         catch(IOException ex) {
             ex.printStackTrace();
         }

    }
    
    @Override
    public void timerUpdate() {
    	
    } 
  
}

