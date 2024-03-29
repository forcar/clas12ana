package org.clas.analysis;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.clas.tools.FitData;

import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;
import org.jlab.utils.groups.IndexedTable;

public class ECmip extends DetectorMonitor {
	
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
    
    IndexedList<List<Float>>          tdcs = new IndexedList<List<Float>>(4);
    IndexedList<List<Float>>          adcs = new IndexedList<List<Float>>(4);
    IndexedList<List<Integer>>       lapmt = new IndexedList<List<Integer>>(3); 
    IndexedList<List<Integer>>       ltpmt = new IndexedList<List<Integer>>(3);    
    
    public ECmip(String name) {
        super(name);
        this.setDetectorTabNames("MIP",
                                 "UVW",
                                 "Fits",
                                 "Mean",
                                 "RMS",
                                 "Maps",
                                 "PID",
                                 "MOM",
                                 "FTOF",
                                 "EFF",
                                 "PIXEL",
                                 "Timeline",
                                 "ATT",
                                 "MIPvP");
        
        this.useRDIFButtons(true);
        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
        this.useZSliderPane(true);
        this.init();
        this.localinit();
    }
    
    public void localinit() {
    	System.out.println("ECmip.localinit()");
    	tl.setFitData(Fits);    	
        getPixLengthMap(outPath+"files/ECpixdepthtotal.dat");
    }  
    
    public void localclear() {
    	System.out.println("ECmip:localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	FitSummary.clear();
    	Fits.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
     @Override    
     public void createHistos(int run) {
	     histosExist = true;
	     System.out.println("ECmip.createHistos("+run+")");
	     setRunNumber(run);
	     runlist.add(run);
	     createMIPHistos(0,1,25,0,40," Peak Energy (MeV)");
	     createMIPHistos(0,2,50,0,100," Cluster Energy (MeV)");	     
	     if(dropSummary) return;
	     createXYHistos(5,200,-420,420,-420,420);    
	     createPIDHistos(6);
	     createMIPHistos(7,1,25,0,5.0," + Momentum (GeV)");
	     createMIPHistos(7,2,25,0,5.0," - Momentum (GeV)");
	     createFTOFHistos(8);
//	     createEFFHistos(9,200,-420,420,-420,420);
//	     createPathHistos(9);
	     createPixHistos(10);	
	     createUVWHistos(12,25,0.,2.," MIP ");
	     createMIPvPHistos(13);
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
    	 plotXYSummary(5); 
    	 plotPIDSummary(6);
    	 plotMIP(7);
    	 plotPIDSummary(8);
//    	 plotPathSummary(9,getActivePC()==1?getActiveView()+1:0);
//    	 plotPathSummary(10,0);
    	 plotUVW(12);    
    	 plotPIDSummary(13);
     }
     
     public void plotAnalysis(int run) {
    	 setRunNumber(run);
    	 if(!isAnalyzeDone) return;
    	 if(!dropSummary) {updateFITS(2); if(TLname=="UVW") {plotMean(); plotRmsSummary(4);}else{plotMeanHWSummary(3); plotRmsHWSummary(4);}}
//    	 if(!dropSummary) {updateFITS(2);plotMeanHWSummary(3); plotRmsHWSummary(4);}
    	 updateUVW(1); plotTimeLines(11);    	    
     }
     
     public void plotMean() {
    	 if(getActiveRDIF()==1 && runlist.size()==2) {getMeanRDIF(); dumpGraphs("rdif");     plotMeanRDIF(3);    return;}    //plot RDIF if two runs present and dump RDIF
    	 if(getActiveRDIF()==0 && runlist.size()>=2) {                                       plotMeanSummary(3); return;} //
    	 if(getActiveRDIF()==0 && runlist.size()==1) {               dumpGraphs("gain");     plotMeanSummary(3); return;} //dump PMT gains based on current analyzed run
    	 if(getActiveRDIF()==1 && runlist.size()==1) {               dumpGraphs("rdifgain"); plotMeanSummary(3);}         //plot effect of RDIF correction and dump RDIF corrected CCDB gains
     }
     
     public void dumpGraphs(String val) {
    	 if(dumpGraphs) writeFile(val,1,7,0,3,0,3);
     }
     
     public void createFTOFHistos(int k) {
    	 
    	 int run = getRunNumber();
    	 
    	 String[] ids = {"P1A ", "P1B "};
    	 String[] ivs = {"GMEAN","TDIF"};
    	 int n = 0;
    	 dg = new DataGroup(6,4);
		 for (int id=0; id<2; id++) {
			 for (int iv=0; iv<2; iv++) {
				 for (int is=1; is<7; is++) {
					 h = new H2F(ids[id]+ivs[iv]+"_"+is+"_"+run,"",50,iv==0?0:-30,iv==0?3000:30,id==0?23:62,1,id==0?24:63);	
					 h.setTitleX("Sector "+is+" "+ids[id]+ivs[iv]); h.setTitleY("Counter Number");
					 dg.addDataSet(h,n); n++;
				 }
    		 }
    	 }
		 
         this.getDataGroup().add(dg,0,0,k,run);
     }
     
     public void createEFFHistos(int k, int nb, int bx1, int bx2, int by1, int by2) {
    	 
 	     int run = getRunNumber();  
 	     
 	     dg = new DataGroup(3,2);
 	     
 	     int n = 0;
 	     for (int i=0; i<1; i++) {
 	    	 for (int d=0; d<3; d++) {
        		 h = new H2F(det[d]+"_"+i+"_"+k+"_"+run,"",nb,bx1,bx2,nb,by1,by2);
                 dg.addDataSet(h,n);n++;  
	    	 }
         } 	     
         this.getDataGroup().add(dg,0,0,k,run);
     }
     
     public void createXYHistos(int k, int nb, int bx1, int bx2, int by1, int by2) {
    	 
 	     int run = getRunNumber();
 	     
         String[] t = {"e","w","r","reff1","reff2"};
         
         //clusters
         for (int i=0; i<5; i++) { //e,w,r,reff1,reff2
        	 dg = new DataGroup(3,2);
        	 for (int d=0; d<3; d++) { //pcal,ecin,ecou
        		 h = new H2F("hi_"+det[d]+"_xyc_"+t[i]+"_"+k+"_"+run,"hi_"+det[d]+"_xyc_"+t[i]+"_"+k+"_"+run,nb,bx1,bx2,nb,by1,by2);
                 dg.addDataSet(h,d);  
	    	 }
             this.getDataGroup().add(dg,i,2,k,run);
         }

         //peaks
         for (int i=0; i<3; i++) { //e,w,r
        	 dg = new DataGroup(3,3);
        	 for (int j=0; j<3; j++) { //u,v,w
        		 for (int d=0; d<3; d++) { //pcal,ecin,ecou
        			 h = new H2F("hi_"+det[d]+"_xyp_"+v[j]+t[i]+"_"+k+"_"+run,"hi_"+det[d]+"_xyp_"+v[j]+t[i]+"_"+k+"_"+run,nb,bx1,bx2,nb,by1,by2);
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
                 h = new H2F("uvw_pcal_u"+ip+"_s"+is+"_"+k+"_"+run,"uvw_pcal_u"+ip+"_s"+is+"_"+k+"_"+run,68,1,69,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" PCAL W"); h.setTitleY(ytxt+"U"+ip);       
                 dg1.addDataSet(h,ip-1); dg1.addDataSet(f1,ip-1);
                 h = new H2F("uvw_pcal_v"+ip+"_s"+is+"_"+k+"_"+run,"uvw_pcal_v"+ip+"_s"+is+"_"+k+"_"+run,68,1,69,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" PCAL U"); h.setTitleY(ytxt+"V"+ip);
                 dg2.addDataSet(h,ip-1); dg2.addDataSet(f1,ip-1);
                 h = new H2F("uvw_pcal_w"+ip+"_s"+is+"_"+k+"_"+run,"uvw_pcal_w"+ip+"_s"+is+"_"+k+"_"+run,68,1,69,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" PCAL V"); h.setTitleY(ytxt+"W"+ip); 
                 dg3.addDataSet(h,ip-1); dg3.addDataSet(f1,ip-1);
      	    }
             this.getDataGroup().add(dg1,is,1,k,run); this.getDataGroup().add(dg2,is,2,k,run); this.getDataGroup().add(dg3,is,3,k,run);
             
             DataGroup dg4 = new DataGroup(6,6); DataGroup dg5 = new DataGroup(6,6); DataGroup dg6 = new DataGroup(6,6);        	         	   
             f1 = new F1D("p0"+is+2+k,"[a]",1,37); f1.setParameter(0,1); f1.setLineColor(1); f1.setLineStyle(1);
      	    for (int ip=1; ip<npmts[1]+1; ip++) {
                 h = new H2F("uvw_ecin_u"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecin_u"+ip+"_s"+is+"_"+k+"_"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECIN W");  h.setTitleY(ytxt+"U"+ip); 
                 dg4.addDataSet(h,ip-1); dg4.addDataSet(f1,ip-1);
                 h = new H2F("uvw_ecin_v"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecin_v"+ip+"_s"+is+"_"+k+"_"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECIN U"); h.setTitleY(ytxt+"V"+ip); 
                 dg5.addDataSet(h,ip-1); dg5.addDataSet(f1,ip-1);
                 h = new H2F("uvw_ecin_w"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecin_w"+ip+"_s"+is+"_"+k+"_"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECIN V"); h.setTitleY(ytxt+"W"+ip);
                 dg6.addDataSet(h,ip-1); dg6.addDataSet(f1,ip-1);
                 
      	    }
             this.getDataGroup().add(dg4,is,4,k,run); this.getDataGroup().add(dg5,is,5,k,run); this.getDataGroup().add(dg6,is,6,k,run);
      	   
             DataGroup dg7 = new DataGroup(6,6); DataGroup dg8 = new DataGroup(6,6); DataGroup dg9 = new DataGroup(6,6);        	         	   
             f1 = new F1D("p0"+is+3+k,"[a]",1,37); f1.setParameter(0,1); f1.setLineColor(1); f1.setLineStyle(1);
      	    for (int ip=1; ip<npmts[2]+1; ip++) {
                 h = new H2F("uvw_ecou_u"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecou_u"+ip+"_s"+is+"_"+k+"_"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECOU W"); h.setTitleY(ytxt+"U"+ip);
                 dg7.addDataSet(h,ip-1); dg7.addDataSet(f1,ip-1);
                 h = new H2F("uvw_ecou_v"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecou_v"+ip+"_s"+is+"_"+k+"_"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECOU U");  h.setTitleY(ytxt+"V"+ip);
                 dg8.addDataSet(h,ip-1); dg8.addDataSet(f1,ip-1);
                 h = new H2F("uvw_ecou_w"+ip+"_s"+is+"_"+k+"_"+run,"uvw_ecou_w"+ip+"_s"+is+"_"+k+"_"+run,36,1,37,ybins,ymin,ymax);
                 h.setTitleX("Sector "+is+" ECOU V"); h.setTitleY(ytxt+"W"+ip);
                 dg9.addDataSet(h,ip-1); dg9.addDataSet(f1,ip-1);    
      	    }
             this.getDataGroup().add(dg7,is,7,k,run); this.getDataGroup().add(dg8,is,8,k,run); this.getDataGroup().add(dg9,is,9,k,run);     	    
         }        
     }     
    
     public void createMIPHistos(int k, int n, int nch, double x1, double x2, String txt) {
    	
	    int run = getRunNumber();
    
        for (int is=1; is<7; is++) {
            String tag = is+"_"+n+"_"+k+"_"+run;
            dg = new DataGroup(3,3);
            h = new H2F("mip_pcal_u_"+tag,"mip_pcal_u_"+tag, nch, x1, x2, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL U"+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("mip_pcal_v_"+tag,"mip_pcal_v_"+tag, nch, x1, x2, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL V"+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("mip_pcal_w_"+tag,"mip_pcal_w_"+tag, nch, x1, x2, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL W"+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
        
            h = new H2F("mip_ecin_u_"+tag,"mip_ecin_u_"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN U"+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("mip_ecin_v_"+tag,"mip_ecin_v_"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN V"+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("mip_ecin_w_"+tag,"mip_ecin_w_"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN W"+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
        
            h = new H2F("mip_ecou_u_"+tag,"mip_ecou_u_"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU U"+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("mip_ecou_v_"+tag,"mip_ecou_v_"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU V"+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("mip_ecou_w_"+tag,"mip_ecou_w_"+tag, nch, x1, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU W"+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,n,k,run);
        }            
    }  
     
    public void createPathHistos(int k) {
    	
	   int run = getRunNumber();
       
       for (int il=0; il<4; il++) {
           dg = new DataGroup(3,4); 
       for (int is=1; is<7; is++) {
           String tag = is+"_"+il+"_"+k+"_"+run;
           h = new H2F("hi_pcal_path1_"+tag,"hi_pcal_path1_"+tag,20,0.,2.,118,31.,50.);
           h.setTitleX("Sector "+is+" PCAL/MIP");
           h.setTitleY("Path12 (cm)");
           dg.addDataSet(h, is-1);  
           h = new H2F("hi_pcal_path2_"+tag,"hi_pcal_path2_"+tag,20,0.,2.,70,50.,70.);
           h.setTitleX("Sector "+is+" PCAL/MIP");
           h.setTitleY("Path13 (cm)");
           dg.addDataSet(h, is+5);  
        }
        this.getDataGroup().add(dg,il,0,k,run);
        }
       
        for (int il=0; il<4; il++) {
            dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"_"+il+"_"+k+"_"+run;
            h = new H2F("hi_ecin_path1_"+tag,"hi_ecin_path1_"+tag,20,0.,2.,70,50.,70.);
            h.setTitleX("Sector "+is+" ECIN/MIP");
            h.setTitleY("Path13 (cm)");
            dg.addDataSet(h, is-1);  
            h = new H2F("hi_ecin_path2_"+tag,"hi_ecin_path2_"+tag,20,0.,2.,66,17.,30.);
            h.setTitleX("Sector "+is+" ECIN/MIP");
            h.setTitleY("Path23 (cm)");
            dg.addDataSet(h, is+5);    
        }
        this.getDataGroup().add(dg,il,1,k,run);
        }
        
        for (int il=0; il<4; il++) {
            dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {        
            String tag = is+"_"+il+"_"+k+"_"+run;
            h = new H2F("hi_ecou_path1_"+tag,"hi_ecou_path1_"+tag,20,0.,2.,70,50.,70.);
            h.setTitleX("Sector "+is+" ECOU/MIP");
            h.setTitleY("Path13 (cm)");
            dg.addDataSet(h, is-1);  
            h = new H2F("hi_ecou_path2_"+tag,"hi_ecou_path2_"+tag,20,0.,2.,66,17.,30.);
            h.setTitleX("Sector "+is+" ECOU/MIP");
            h.setTitleY("Path23 (cm)");
            dg.addDataSet(h, is+5);      
        }
        this.getDataGroup().add(dg,il,2,k,run);
        }
                
    }
    
    public void createPixHistos(int k) {
    	
 	    int run = getRunNumber();
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_pcal_pix1_"+tag,"hi_pcal_pix1_"+tag,50,1.,200.,12,3.,15.);
            h.setTitleX("Sector "+is+" PCAL (MeV)");
            h.setTitleY("No. Strips");
            dg.addDataSet(h, is-1);   
            h = new H2F("hi_pcal_pix2_"+tag,"hi_pcal_pix2_"+tag,100,31.,40.,12,3.,15.);
            h.setTitleX("Sector "+is+" Path12 (cm)");
            h.setTitleY("No. Strips");
            dg.addDataSet(h, is+5);   
         }
        this.getDataGroup().add(dg,0,0,k,run);
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_ecin_pix_"+tag,"hi_ecin_pix_"+tag,50,1.,200.,12,3.,15.);
            h.setTitleX("Sector "+is+" ECin (MeV)");
            h.setTitleY("No. Strips");
            dg.addDataSet(h, is-1);   
            h = new H2F("hi_ecin_pix2_"+tag,"hi_ecin_pix2_"+tag,100,31.,40.,12,3.,15.);
            h.setTitleX("Sector "+is+" Path12 (cm)");
            h.setTitleY("No. Strips");
            dg.addDataSet(h, is+5);   
         }
        this.getDataGroup().add(dg,0,1,k,run);
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_ecou_pix_"+tag,"hi_ecou_pix_"+tag,50,1.,200.,12,3.,15.);
            h.setTitleX("Sector "+is+" ECou (MeV)");
            h.setTitleY("No. Strips");
            dg.addDataSet(h, is-1);   
            h = new H2F("hi_ecou_pix2_"+tag,"hi_ecou_pix2_"+tag,100,31.,40.,12,3.,15.);
            h.setTitleX("Sector "+is+" Path12 (cm)");
            h.setTitleY("No. Strips");
            dg.addDataSet(h, is+5);   
         }
        this.getDataGroup().add(dg,0,2,k,run);
    	
    }
    
    public void createPIDHistos(int k)  {
    
 	    int run = getRunNumber();
	    int is  = 0;
        String tag = is+"_"+run;
        
        dg = new DataGroup(2,3);
        
        F1D f1 = new F1D("f_1"+tag,"1/(1+[a]^2/x^2)^0.5", 0.41,3.5); f1.setParameter(0,0.13957); f1.setLineColor(1); f1.setLineStyle(1);    
        F1D f2 = new F1D("f_2"+tag,"1/(1+[a]^2/x^2)^0.5", 0.41,3.5); f2.setParameter(0,0.93827); f2.setLineColor(1); f2.setLineStyle(1);   
        h = new H2F("pid_pos_"+tag,"pid_pos_"+tag,100,0.,3.5,100,0.4,1.5);       h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 0); dg.addDataSet(f1,0); dg.addDataSet(f2,0); 
        h = new H2F("pid_neg_"+tag,"pid_neg_"+tag,100,0.,3.5,100,0.4,1.5);       h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 1); dg.addDataSet(f1,1); 
        h = new H2F("pid_fc_pos_"+tag,"pid_fc_pos_"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 2); dg.addDataSet(f1,2); dg.addDataSet(f2,2); 
        h = new H2F("pid_fc_neg_"+tag,"pid_fc_neg_"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 3); dg.addDataSet(f1,3); 
        h = new H2F("pid_fc_ppi_"+tag,"pid_fc_ppi_"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 4); dg.addDataSet(f1,4); dg.addDataSet(f2,4); 
        h = new H2F("pid_fc_npi_"+tag,"pid_fc_npi_"+tag,100,0.,3.5,100,0.4,1.5); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 5); dg.addDataSet(f1,5); 
        
        this.getDataGroup().add(dg,is,0,k,run);
    }
    
    public void createMIPvPHistos(int k) {
 	    int run = getRunNumber(); 
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_pcal_dedx_pm"+tag,"hi_pcal_dedx_pm"+tag,50,0.3,6.,30,0.5,2);
            h.setTitleX("Sector "+is+" P (GeV)");
            h.setTitleY("MEAN/MIP");
            dg.addDataSet(h, is-1);   
            h = new H2F("hi_pcal_dedx_pp"+tag,"hi_pcal_dedx_pp"+tag,100,0.3,6.,30,0.5,2);
            h.setTitleX("Sector "+is+" P (GeV)");
            h.setTitleY("MEAN/MIP");
            dg.addDataSet(h, is+5);   
         }
        this.getDataGroup().add(dg,0,0,k,run); 
        
    }
       
    public void initCCDB(int runno) {
    	System.out.println("ECmip.initCCDB("+runno+")");
        gain    = cm.getConstants(runno, "/calibration/ec/gain");
        time    = cm.getConstants(runno, "/calibration/ec/timing");
        veff    = cm.getConstants(runno, "/calibration/ec/effective_velocity");
        offset  = cm.getConstants(runno, "/calibration/ec/fadc_offset");
        goffset = cm.getConstants(runno, "/calibration/ec/fadc_global_offset");    	
        shift   = cm.getConstants(runno, "/calibration/ec/torus_gain_shift");        
    }
    
    public void processFTOFbank(DataEvent event) {
        
        float      tps =  (float) 0.02345;
        float     tdcd = 0;
        
        adcs.clear(); tdcs.clear(); ltpmt.clear() ; lapmt.clear();
        
        if(event.hasBank("FTOF::tdc")){
            DataBank  bank = event.getBank("FTOF::tdc");
            int rows = bank.rows();           
            for(int i = 0; i < rows; i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  lr = bank.getByte("order",i);                       
                int  ip = bank.getShort("component",i);
                tdcd = bank.getInt("TDC",i)*tps;  
                
                if(tdcd>0) {
                if(!tdcs.hasItem(is,il,lr-2,ip)) tdcs.add(new ArrayList<Float>(),is,il,lr-2,ip);
                    tdcs.getItem(is,il,lr-2,ip).add(tdcd); 
                    if (!ltpmt.hasItem(is,il,ip)) {
                 	    ltpmt.add(new ArrayList<Integer>(),is,il,ip);
                 	    ltpmt.getItem(is,il,ip).add(ip);
                    } 
                }
            }
        }
            
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
                
                if(true) {
                
                if(adc>0) {
                if(!adcs.hasItem(is,il,lr,ip)) adcs.add(new ArrayList<Float>(),is,il,lr,ip);
                    adcs.getItem(is,il,lr,ip).add((float) adc); 
                    if (!lapmt.hasItem(is,il,ip)) {
           	            lapmt.add(new ArrayList<Integer>(),is,il,ip);
                        lapmt.getItem(is,il,ip).add(ip);              
                    }
                } 
                
                } //isGoodSector?
            }
        }
    }
    
    public Boolean processFTOF() {
 	   
        IndexGenerator ig = new IndexGenerator();
        int run = getRunNumber();
        float gmtest=0;
        
        for (Map.Entry<Long,List<Integer>>  entry : lapmt.getMap().entrySet()){
            long hash = entry.getKey();
            int is = ig.getIndex(hash, 0);
            int il = ig.getIndex(hash, 1);
            int ip = ig.getIndex(hash, 2);
                   
         	if(adcs.hasItem(is,il,0,ip)&&adcs.hasItem(is,il,1,ip)) {
         	   float gm = (float) Math.sqrt(adcs.getItem(is,il,0,ip).get(0)*
                                      adcs.getItem(is,il,1,ip).get(0));
         	   ((H2F) this.getDataGroup().getItem(0,0,8,run).getData(is-1+(il-1)*12).get(0)).fill(gm, ip);         	
               if(il==2) gmtest=gm;	
         	}

        }
        
        for (Map.Entry<Long,List<Integer>>  entry : ltpmt.getMap().entrySet()){
            long hash = entry.getKey();
            int is = ig.getIndex(hash, 0);
            int il = ig.getIndex(hash, 1);
            int ip = ig.getIndex(hash, 2);
             	   
         	if(tdcs.hasItem(is,il,0,ip)&&tdcs.hasItem(is,il,1,ip)) {
         		float td = tdcs.getItem(is,il,0,ip).get(0)-
         				   tdcs.getItem(is,il,1,ip).get(0);
          	   ((H2F) this.getDataGroup().getItem(0,0,8,run).getData(is-1+(il-1)*12+6).get(0)).fill(td, ip);         	
         	}
        }       
        return gmtest>100;
    } 
    
    public void processEvent(DataEvent event) {
    	
       isMC = (getRunNumber()<100) ? true:false;
       
       IndexedList<List<Particle>> ecpart = new IndexedList<List<Particle>>(2);
       IndexedList<Vector3>            rl = new IndexedList<Vector3>(2);  
       
	   Boolean isMuon = true, goodEvt = false, goodRows = false, goodTrig = false, goodPC,goodECi,goodECo;       
       int nrows = 0, trigger = 0, trig = TRpid, run = getRunNumber();
       double startTime = -1000;

       DataBank bank1    = null;
       DataBank bank2    = null;
       DataBank bank3    = null;
       DataBank bank4    = null;
       DataBank runConf  = null;
       DataBank recPart  = null;
       DataBank recScin  = null;
       DataBank recCalo  = null;
       DataBank recTraj  = null;
       DataBank recEven  = null;
       DataBank ecalClus = null;
       DataBank ecalCali = null;
       DataBank ecalPeak = null;
       DataBank ftofRaw  = null;
       DataBank ftofTrue = null;
       
	   if(dropBanks) dropBanks(event);

       if(event.hasBank("RUN::config"))            runConf  = event.getBank("RUN::config");
       if(event.hasBank("REC::Particle"))          recPart  = event.getBank("REC::Particle");
       if(event.hasBank("REC::Scintillator"))      recScin  = event.getBank("REC::Scintillator");
       if(event.hasBank("REC::Calorimeter"))       recCalo  = event.getBank("REC::Calorimeter");
       if(event.hasBank("REC::Traj"))              recTraj  = event.getBank("REC::Traj");
       if(event.hasBank("REC::Event"))             recEven  = event.getBank("REC::Event");
       if(event.hasBank("ECAL::clusters"))         ecalClus = event.getBank("ECAL::clusters");
       if(event.hasBank("ECAL::calib"))            ecalCali = event.getBank("ECAL::calib");
       if(event.hasBank("ECAL::peaks"))            ecalPeak = event.getBank("ECAL::peaks");
       if(event.hasBank("MC::True"))               ftofTrue = event.getBank("MC::True");
       
       if(runConf == null) return;
       
       part.clear();
   
       processFTOFbank(event);  //For use as trigger for ECAL efficency studies
       goodTrig = processFTOF();
       
       float pcx=-1000,pcy=-1000,pcz,tmax=30;
       if(ftofTrue!=null) {
           for(int i=0; i < ftofTrue.rows(); i++) {
               if(goodTrig) {
                  float pcX = ftofTrue.getFloat("avgX",i);
                  float pcY = ftofTrue.getFloat("avgY",i);
                  float pcZ = ftofTrue.getFloat("avgZ",i);
                  float pcT = ftofTrue.getFloat("avgT",i);
                  if(pcT<tmax){pcx=pcX; pcy=pcY; pcz=pcZ ; tmax = pcT;}
               }
           }
       } 
       
//       ((H2F) this.getDataGroup().getItem(0,0,8,run).getData(24).get(0)).fill(-pcx*0.1,pcy*0.1);
       
       
       if (recPart!=null && recEven!=null) { //physics run data
    	    isMuon = false;  
    	    startTime = (isHipo3Event) ? recEven.getFloat("STTime", 0):
                                         recEven.getFloat("startTime", 0);           
            if(startTime > -100) {             	
                for(int i = 0; i < recPart.rows(); i++){              	
                    int      pid = recPart.getInt("pid", i);              
                    float     px = recPart.getFloat("px", i);
                    float     py = recPart.getFloat("py", i);
                    float     pz = recPart.getFloat("pz", i);
                    float     vx = recPart.getFloat("vx", i);
                    float     vy = recPart.getFloat("vy", i);
                    float     vz = recPart.getFloat("vz", i);
                    float   beta = recPart.getFloat("beta", i);
                    short status = (short) Math.abs(recPart.getShort("status", i));
                    if (i==0) trigger = Math.abs(pid);
                    Particle p = new Particle(); 
                    if (pid==0) {p.setProperty("index", i); p.setProperty("ppid", 0);}             
                    if (pid!=0) {
                    	p.initParticle(pid, px, py, pz, vx, vy, vz);                    	   
                    	p.setProperty("ppid", pid);
                        p.setProperty("status", status);
                        p.setProperty("beta", beta);
                        p.setProperty("index",i);                        
                    }
                    part.add(i,p);     //Lists do not support sparse indices !!!!!            
                }                
            }                       
       }
    		   
       if (!part.isEmpty() && trigger==trig) {
           for (Particle p: part) {
        	   if (p.getProperty("ppid")!=0 && p.getProperty("index")>0 && p.charge()!=0 && p.getProperty("status")>=2000 && p.getProperty("status")<3000) {
        		   ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(p.charge()>0?0:1).get(0)).fill(p.p(),p.getProperty("beta")); 
        	   }
           }
       }
 
       Boolean goodREC  = recCalo!=null, goodECAL = ecalClus!=null && ecalCali!=null;
       
       if( isMuon)                     goodEvt = goodECAL;
       if(!isMuon&&goodREC&&goodECAL) {goodEvt = true; bank1 = recCalo;}
       
       if (goodEvt) {
    	   
           bank2 = ecalClus; 
           bank3 = ecalCali;
           bank4 = ftofRaw;
           
           if ( isMuon) {goodRows = true;                                                   nrows = bank2.rows();}
           if (!isMuon) {goodRows = bank1.rows()==bank2.rows()&&bank1.rows()==bank3.rows(); nrows = bank1.rows();}
           
           int[] n1 = new int[6]; int[] n4 = new int[6]; int[] n7 = new int[6];
           float[][]   w1 = new float[6][20] ; 
           float[][]   w4 = new float[6][20] ; 
           float[][]   w7 = new float[6][20] ; 
           float[][]  e1c = new float[6][20]; float[][][]   e1p = new float[6][3][20]; 
           float[][]  e4c = new float[6][20]; float[][][]   e4p = new float[6][3][20]; 
           float[][]  e7c = new float[6][20]; float[][][]   e7p = new float[6][3][20]; 
           float[][]  p1c = new float[6][20]; float[][][]   p1p = new float[6][3][20]; 
           float[][]  p4c = new float[6][20]; float[][][]   p4p = new float[6][3][20]; 
           float[][]  p7c = new float[6][20]; float[][][]   p7p = new float[6][3][20];  
           float[][][] cU = new float[6][3][20];
           float[][][] cV = new float[6][3][20];
           float[][][] cW = new float[6][3][20];
           
           if(goodRows) { //When re-running ECEngine this may not always be true
          
           for(int loop = 0; loop < nrows; loop++){ //loop over REC::calorimeter (!isMuon) or ECAL::cluster (isMuon)
        	   int ic=0, is=0, il=0;
	           if (!isMuon) {ic = bank1.getShort("index",  loop); is = bank1.getByte("sector",loop); il = bank1.getByte("layer",loop);}
	           if ( isMuon) {ic = loop;                           is = bank2.getByte("sector",loop); il = bank2.getByte("layer",loop);}          
	           float   en = bank2.getFloat("energy",ic)*1000;	    	   
               float   ti = bank2.getFloat("time",ic)*1000;
               float    x = bank2.getFloat("x", ic);
               float    y = bank2.getFloat("y", ic);
               float    z = bank2.getFloat("z", ic);
               int     iU = (bank2.getInt("coordU", ic)-4)/8+1;
               int     iV = (bank2.getInt("coordV", ic)-4)/8+1;
               int     iW = (bank2.getInt("coordW", ic)-4)/8+1;
               float   wu = bank2.getFloat("widthU",ic);
               float   wv = bank2.getFloat("widthV",ic);
               float   ww = bank2.getFloat("widthW",ic);
               float  enu = bank3.getFloat("recEU", ic)*1000;
               float  env = bank3.getFloat("recEV", ic)*1000;
               float  enw = bank3.getFloat("recEW", ic)*1000;   
               
               float wsum = wu+wv+ww;
               float pm = -100, pp = -100;
               
               Particle p = new Particle(); 
               int ppid = 0;
               
               if(!isMuon) { //PID diagnostics
    	         int     ip = bank1.getShort("pindex", loop);
                 p.copy(part.get(ip));
                 p.setProperty("beta",part.get(ip).getProperty("beta"));
                 p.setProperty("status",part.get(ip).getProperty("status"));
                 p.setProperty("ppid",part.get(ip).getProperty("ppid"));
                 ppid = (int) part.get(ip).getProperty("ppid");
                 if (trigger==trig && ppid!=11 && ppid!=0 && p.charge()!=0 && p.getProperty("status")>=2000) {
                   p.setProperty("energy",en);
    	           p.setProperty("time",ti);
    	           p.setProperty("x",x);
    	           p.setProperty("y",y);
    	           p.setProperty("z",z);
                   p.setProperty("iU",iU);
                   p.setProperty("iV",iV);
                   p.setProperty("iW",iW);
                   p.setProperty("enu",enu);
                   p.setProperty("env",env);
                   p.setProperty("enw",enw);
                                
                   if (p.charge()>0) {pp = (float) p.p(); ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(2).get(0)).fill(pp,p.getProperty("beta"));}
                   if (p.charge()<0) {pm = (float) p.p(); ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(3).get(0)).fill(pm,p.getProperty("beta"));} 
               
                   if (ppid==+211) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(4).get(0)).fill(pp,p.getProperty("beta"));
                   if (ppid==-211) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(5).get(0)).fill(pm,p.getProperty("beta"));
                 
                   if (!ecpart.hasItem(is,il)) ecpart.add(new ArrayList<Particle>(), is,il);    	                
 	                    ecpart.getItem(is,il).add(p);  
                 }
               }
               
        
//               Boolean goodPID = isMuon ? true:(Math.abs(p.pid())==211&&p.getProperty("status")>=2000);
               Boolean goodPID = isMuon ? true:(Math.abs(ppid)==211 && p.getProperty("status")>=2000 && p.getProperty("status")<3000);
   	                         
               Vector3 r = new Vector3(x,y,z);
               
               goodPC = goodPID&&il==1&&n1[is-1]<20;  goodECi = goodPID&&il==4&&n4[is-1]<20;  goodECo = goodPID&&il==7&&n7[is-1]<20; 

               if (goodPC)  {e1c[is-1][n1[is-1]]=en; rl.add(r,is,0); cU[is-1][0][n1[is-1]]=iU; cV[is-1][0][n1[is-1]]=iV; cW[is-1][0][n1[is-1]]=iW; p1c[is-1][n1[is-1]]=pm;}
               if (goodECi) {e4c[is-1][n4[is-1]]=en; rl.add(r,is,1); cU[is-1][1][n4[is-1]]=iU; cV[is-1][1][n4[is-1]]=iV; cW[is-1][1][n4[is-1]]=iW; p4c[is-1][n4[is-1]]=pm;}
               if (goodECo) {e7c[is-1][n7[is-1]]=en; rl.add(r,is,2); cU[is-1][2][n7[is-1]]=iU; cV[is-1][2][n7[is-1]]=iV; cW[is-1][2][n7[is-1]]=iW; p7c[is-1][n7[is-1]]=pm;}
               if (goodPC)  {p1p[is-1][0][n1[is-1]]=pp;  p1p[is-1][1][n1[is-1]]=pp;  p1p[is-1][2][n1[is-1]]=pp;  w1[is-1][n1[is-1]]=wsum;}
               if (goodECi) {p4p[is-1][0][n4[is-1]]=pp;  p4p[is-1][1][n4[is-1]]=pp;  p4p[is-1][2][n4[is-1]]=pp;  w4[is-1][n4[is-1]]=wsum;}
               if (goodECo) {p7p[is-1][0][n7[is-1]]=pp;  p7p[is-1][1][n7[is-1]]=pp;  p7p[is-1][2][n7[is-1]]=pp;  w7[is-1][n7[is-1]]=wsum;}
               if (goodPC)  {e1p[is-1][0][n1[is-1]]=enu; e1p[is-1][1][n1[is-1]]=env; e1p[is-1][2][n1[is-1]]=enw; n1[is-1]++;}
               if (goodECi) {e4p[is-1][0][n4[is-1]]=enu; e4p[is-1][1][n4[is-1]]=env; e4p[is-1][2][n4[is-1]]=enw; n4[is-1]++;}
               if (goodECo) {e7p[is-1][0][n7[is-1]]=enu; e7p[is-1][1][n7[is-1]]=env; e7p[is-1][2][n7[is-1]]=enw; n7[is-1]++;}
               
           }
           }
             
            for (int is=0; is<6; is++) {
                int iis = is+1;
//                if (isGoodTrigger(iis)) {
//                if(n1[is]>=1&&n1[is]<=4&&n4[is]>=1&&n4[is]<=4) { //Cut out vertical cosmic rays
               boolean pcMuon = n1[is]==1 || n4[is]==1 || n7[is]==1;
               boolean   Muon = n1[is]==1 && n4[is]==1 && n7[is]==1;
                
               Boolean mtest = (trig==-1 || isMC) ? pcMuon : Muon; //7.29.2023
               
               if(trigger==(trig==-1?0:trig) && mtest  && iis==(isMC ? iis:iis)) { //Only one cluster in each layer to reject vertical cosmics
            	   
            	   
//                    Boolean goodU = Math.abs(cU[is][1][n4[is]]-cU[is][2][n7[is]])<=1;
//                    Boolean goodV = Math.abs(cV[is][1][n4[is]]-cV[is][2][n7[is]])<=1;
//                    Boolean goodW = Math.abs(cW[is][1][n4[is]]-cW[is][2][n7[is]])<=1;
//                    Boolean goodUVW = goodU&&goodV&&goodW;
//                if(is==1&&partRecEB==null) System.out.println("No particle found");
//                if(is==1&&partRecEB!=null) System.out.println("Found Sector= "+partRecEB.getProperty("sector")+" P= "+partRecEB.p());
//                if(is==1&&partRecEB!=null) System.out.println("Energy1,e1 "+partRecEB.getProperty("energy1")+" "+e1[is][0]);
//                Boolean goodPion = (is==1&&partRecEB!=null&&partRecEB.p()>0.7);
//                if(is==1&&!goodPion) {n1[is]=0;n4[is]=0;n7[is]=0;}
// Target-PCAL: 6977.8 mm CCDB:/geometry/pcal/dist2tgt
// Target-ECAL: 7303.3 mm CCDB:/geometry/ec/dist2tgt
// PCAL-ECin: 325.5 mm
// ECin-ECou: 14*2.2+15*10.0 = 180.8 mm
// PCAL-ECou: 325.5+180.8= 506.3 mm

                float v12mag = 0;
                float v13mag = 0;
                float v23mag = 0;
              	               	   
               	if(!isMC && trig>-1) {
                   Vector3  v1 = new Vector3(rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),rl.getItem(iis,0).z());
                   Vector3  v2 = new Vector3(rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),rl.getItem(iis,1).z());
                   Vector3  v3 = new Vector3(rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),rl.getItem(iis,2).z());
                   Vector3 v23 = new Vector3(v2.x(),v2.y(),v2.z());
           
                   v2.sub(v1); v23.sub(v3); v3.sub(v1);  
                   
                   v12mag = (float) v2.mag();
                   v13mag = (float) v3.mag();
                   v23mag = (float) v23.mag();
               	}            	   
            	
                H2F h;
                double ectot = e4c[is][0]+e7c[is][0] ; double etot = e1c[is][0]+ectot ;
/*                
                h = (H2F) this.getDataGroup().getItem(0,0,5,run).getData(3).get(0); h.fill(e1c[is][0],ectot);
                dg4.getH2F("hi_pcal_ectot_"+iis).fill(e1c[is][0],ectot);
                dg4.getH2F("hi_pcal_ectot_max_"+iis).fill(e1c[is][0],ectot);
                
                if (iis==2&&!pmap.isEmpty()) {                
                for (float p: pmap) {
                    dg4.getH2F("hi_pcal_1").fill(e1c[is][0],p);
                    dg4.getH2F("hi_ecali_1").fill(e4c[is][0],p);
                    dg4.getH2F("hi_ecalo_1").fill(e7c[is][0],p);
                    dg4.getH2F("hi_etot_1").fill(p,etot*1e-3/p);  
                }
                }
*/ 

                if(this.getDataGroup().hasItem(0,0,9,run)) {
/*            		
                for (int il=0; il<4; il++) {
                ((H2F) this.getDataGroup().getItem(il,0,9,run).getData(iis-1).get(0)).fill((il==0)?e1c[is][0]/mipc[0]:e1p[is][il-1][0]/mipp[0],v12mag);
                ((H2F) this.getDataGroup().getItem(il,0,9,run).getData(iis+5).get(0)).fill((il==0)?e1c[is][0]/mipc[0]:e1p[is][il-1][0]/mipp[0],v13mag);
                ((H2F) this.getDataGroup().getItem(il,1,9,run).getData(iis-1).get(0)).fill((il==0)?e4c[is][0]/mipc[1]:e4p[is][il-1][0]/mipp[1],v13mag);
                ((H2F) this.getDataGroup().getItem(il,1,9,run).getData(iis+5).get(0)).fill((il==0)?e4c[is][0]/mipc[1]:e4p[is][il-1][0]/mipp[1],v23mag);
                ((H2F) this.getDataGroup().getItem(il,2,9,run).getData(iis-1).get(0)).fill((il==0)?e7c[is][0]/mipc[2]:e7p[is][il-1][0]/mipp[2],v13mag);
                ((H2F) this.getDataGroup().getItem(il,2,9,run).getData(iis+5).get(0)).fill((il==0)?e7c[is][0]/mipc[2]:e7p[is][il-1][0]/mipp[2],v23mag);
                }
*/             
                ((H2F) this.getDataGroup().getItem(0,0,10,run).getData(iis-1).get(0)).fill(e1c[is][0],w1[is][0]);
                ((H2F) this.getDataGroup().getItem(0,1,10,run).getData(iis-1).get(0)).fill(e4c[is][0],w4[is][0]);
                ((H2F) this.getDataGroup().getItem(0,2,10,run).getData(iis-1).get(0)).fill(e7c[is][0],w7[is][0]);
//               ((H2F) this.getDataGroup().getItem(0,0,10,run).getData(iis+5).get(0)).fill(v12mag,e1c[is][0]-e4c[is][0]);
//               ((H2F) this.getDataGroup().getItem(0,1,10,run).getData(iis+5).get(0)).fill(v13mag,e1c[is][0]-e7c[is][0]+18);
//               ((H2F) this.getDataGroup().getItem(0,2,10,run).getData(iis+5).get(0)).fill(v23mag,e4c[is][0]-e7c[is][0]+18);
                ((H2F) this.getDataGroup().getItem(0,0,10,run).getData(iis+5).get(0)).fill(v12mag,w1[is][0]);
                ((H2F) this.getDataGroup().getItem(0,1,10,run).getData(iis+5).get(0)).fill(v12mag,w4[is][0]);
                ((H2F) this.getDataGroup().getItem(0,2,10,run).getData(iis+5).get(0)).fill(v12mag,w7[is][0]);
               
                }
                
                //Muon pixel cut too restrictive for physics pions in PCAL due to target/Bfield constrained tracks
                //For MC efficiency studies of EC do not constrain PCAL track
                Boolean pcaltest = (isMC) ? n1[is]==1 && w1[is][0]>=3 : (isMuon) ? v12mag<35 && w1[is][0]==3 : v13mag<56&&(w1[is][0]==3||w1[is][0]==4);
                Boolean ecintest = (isMC) ? n4[is]==1 && w4[is][0]>=3 : v23mag<21&&w4[is][0]==3; 
                Boolean ecoutest = (isMC) ? n7[is]==1 && w7[is][0]>=3 : v23mag<21&&w7[is][0]==3;
                
                if(pcaltest) { 
                for(int n=0; n<n1[is]; n++) {
                	int u = (int) cU[is][0][n]; int v = (int) cV[is][0][n]; int w = (int) cW[is][0][n]; 
                	float ec = e1c[is][n] ; float epu = e1p[is][0][n] ; float epv = e1p[is][1][n] ; float epw = e1p[is][2][n] ;
                	((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(0).get(0)).fill(ec,u);
                	((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(1).get(0)).fill(ec,v);
                	((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(2).get(0)).fill(ec,w);
                	
                	((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(0).get(0)).fill(epu,u);
                    ((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(1).get(0)).fill(epv,v);
                	((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(2).get(0)).fill(epw,w);    
                	
                	((H2F) this.getDataGroup().getItem(iis,1,12,run).getData(u-1).get(0)).fill(w,epu/mipp[0]);
                	((H2F) this.getDataGroup().getItem(iis,2,12,run).getData(v-1).get(0)).fill(u,epv/mipp[0]);
                	((H2F) this.getDataGroup().getItem(iis,3,12,run).getData(w-1).get(0)).fill(v,epw/mipp[0]);
               		
                	((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(0).get(0)).fill(p1c[is][n],   u);
                	((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(1).get(0)).fill(p1c[is][n],   v);
                	((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(2).get(0)).fill(p1c[is][n],   w);
                	((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(0).get(0)).fill(p1p[is][0][n],u);
                	((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(1).get(0)).fill(p1p[is][1][n],v);
                	((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(2).get(0)).fill(p1p[is][2][n],w); 
                	
                	((H2F) this.getDataGroup().getItem(0,0,13,run).getData(is).get(0)).fill(p1c[is][n],ec/mipc[0]); 
                	((H2F) this.getDataGroup().getItem(0,0,13,run).getData(is+6).get(0)).fill(p1p[is][0][n],ec/mipc[0]); 
               		
                	if(ecoutest) ((H2F) this.getDataGroup().getItem(4,  2,5,run).getData(0).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),ec<mxc[0]?1:0);
                	             ((H2F) this.getDataGroup().getItem(3,  2,5,run).getData(0).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),ec<mxc[0]?1:0);
//                	((H2F) this.getDataGroup().getItem(3,  2,5,run).getData(0).get(0)).fill(-pcx*0.1,pcy*0.1,ec<mxc[0]?1:0);
                	
                	((H2F) this.getDataGroup().getItem(0,  2,5,run).getData(0).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),ec<mxc[0]?mipc[0]:0);
                    ((H2F) this.getDataGroup().getItem(1,  2,5,run).getData(0).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),ec<mxc[0]?ec:0); 
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(0).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),epu<mxp[0]?mipp[0]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(0).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),epu<mxp[0]?epu:0);
                    ((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(1).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),epv<mxp[0]?mipp[0]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(1).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),epv<mxp[0]?epv:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(2).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),epw<mxp[0]?mipp[0]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(2).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),epw<mxp[0]?epw:0);
                }
                }
                
                if(ecintest) {
                int id = 1;
                for(int n=0; n<n4[is]; n++) {
                	int u = (int) cU[is][1][n]; int v = (int) cV[is][1][n]; int w = (int) cW[is][1][n]; 
                	float d = (PixLength.hasItem(u,v,w))? PixLength.getItem(u,v,w):1f; d=1;
                	float ec = e4c[is][n]/d ; float epu = e4p[is][0][n]/d ; float epv = e4p[is][1][n]/d ; float epw = e4p[is][2][n]/d ;
                	((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(3).get(0)).fill(ec,u);
                	((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(4).get(0)).fill(ec,v);
                	((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(5).get(0)).fill(ec,w);
                	((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(3).get(0)).fill(epu,u);
                	((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(4).get(0)).fill(epv,v);
                	((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(5).get(0)).fill(epw,w);
                	((H2F) this.getDataGroup().getItem(iis,4,12,run).getData(u-1).get(0)).fill(w,epu/mipp[1]);
                	((H2F) this.getDataGroup().getItem(iis,5,12,run).getData(v-1).get(0)).fill(u,epv/mipp[1]);
                	((H2F) this.getDataGroup().getItem(iis,6,12,run).getData(w-1).get(0)).fill(v,epw/mipp[1]);
                	
                	((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(3).get(0)).fill(p4c[is][n],   u);
                	((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(4).get(0)).fill(p4c[is][n],   v);
                	((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(5).get(0)).fill(p4c[is][n],   w);
                	((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(3).get(0)).fill(p4p[is][0][n],u);
                	((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(4).get(0)).fill(p4p[is][1][n],v);
                	((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(5).get(0)).fill(p4p[is][2][n],w);
                	
//                	((H2F) this.getDataGroup().getItem(iis,2,13,run).getData(1).get(0)).fill(p4c[is][n],ec/mipc[1]); 
//                	((H2F) this.getDataGroup().getItem(iis,1,13,run).getData(1).get(0)).fill(p4p[is][0][n],ec/mipc[1]); 
                   
//                	if(isMC) {id=0; d=1;  ec = e1c[is][n]/d; epu = e1p[is][0][n]/d; epv = e1p[is][1][n]/d; epw = e1p[is][2][n]/d; }
                	
                	if(pcaltest) ((H2F) this.getDataGroup().getItem(3,  2,5,run).getData(1).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),ec<mxc[id]?1:0);
//                	((H2F) this.getDataGroup().getItem(3,  2,5,run).getData(1).get(0)).fill(-pcx*0.1,pcy*0.1,ec<mxc[id]?1:0);
                	
                	((H2F) this.getDataGroup().getItem(0,  2,5,run).getData(1).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),ec<mxc[id]?mipc[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  2,5,run).getData(1).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),ec<mxc[id]?ec:0); 
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(3).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epu<mxp[id]?mipp[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(3).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epu<mxp[id]?epu:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(4).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epv<mxp[id]?mipp[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(4).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epv<mxp[id]?epv:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(5).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epw<mxp[id]?mipp[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(5).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epw<mxp[id]?epw:0);
                }
                }
                                
                if(ecoutest) {                
                int id = 2;
                for(int n=0; n<n7[is]; n++) {
                	int u = (int) cU[is][2][n]; int v = (int) cV[is][2][n]; int w = (int) cW[is][2][n]; 
                	float d = (PixLength.hasItem(u,v,w))? PixLength.getItem(u,v,w):1f; d=1;
                	float ec = e7c[is][n]/d ; float epu = e7p[is][0][n]/d ; float epv = e7p[is][1][n]/d ; float epw = e7p[is][2][n]/d ;
            		((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(6).get(0)).fill(ec,u);
            	    ((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(7).get(0)).fill(ec,v);
            		((H2F) this.getDataGroup().getItem(iis,2,0,run).getData(8).get(0)).fill(ec,w);
            		((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(6).get(0)).fill(epu,u);
            		((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(7).get(0)).fill(epv,v);
            		((H2F) this.getDataGroup().getItem(iis,1,0,run).getData(8).get(0)).fill(epw,w);
                	((H2F) this.getDataGroup().getItem(iis,7,12,run).getData(u-1).get(0)).fill(w,epu/mipp[2]);
                	((H2F) this.getDataGroup().getItem(iis,8,12,run).getData(v-1).get(0)).fill(u,epv/mipp[2]);
                	((H2F) this.getDataGroup().getItem(iis,9,12,run).getData(w-1).get(0)).fill(v,epw/mipp[2]);
            		    
            		((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(6).get(0)).fill(p7c[is][n],   u);
            		((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(7).get(0)).fill(p7c[is][n],   v);
            		((H2F) this.getDataGroup().getItem(iis,2,7,run).getData(8).get(0)).fill(p7c[is][n],   w);
            		((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(6).get(0)).fill(p7p[is][0][n],u);
            		((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(7).get(0)).fill(p7p[is][1][n],v);
            		((H2F) this.getDataGroup().getItem(iis,1,7,run).getData(8).get(0)).fill(p7p[is][2][n],w);
            		
//                	((H2F) this.getDataGroup().getItem(iis,2,13,run).getData(2).get(0)).fill(p7c[is][n],ec/mipc[2]); 
//                	((H2F) this.getDataGroup().getItem(iis,1,13,run).getData(2).get(0)).fill(p7p[is][0][n],ec/mipc[2]); 
                    
//                	if(isMC) {id=0; d=1; ec = e1c[is][n]/d; epu = e1p[is][0][n]/d; epv = e1p[is][1][n]/d; epw = e1p[is][2][n]/d; }
            		
                	if(pcaltest) ((H2F) this.getDataGroup().getItem(3,  2,5,run).getData(2).get(0)).fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),ec<mxc[id]?1:0);
                	             ((H2F) this.getDataGroup().getItem(4,  2,5,run).getData(2).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),ec<mxc[id]?1:0);
                	
//                	((H2F) this.getDataGroup().getItem(3,  2,5,run).getData(2).get(0)).fill(-pcx*0.1,pcy*0.1,ec<mxc[id]?1:0);
                	
                	((H2F) this.getDataGroup().getItem(0,  2,5,run).getData(2).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),ec<mxc[id]?mipc[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  2,5,run).getData(2).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),ec<mxc[id]?ec:0); 
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(6).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epu<mxp[id]?mipp[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(6).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epu<mxp[id]?epu:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(7).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epv<mxp[id]?mipp[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(7).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epv<mxp[id]?epv:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(8).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epw<mxp[id]?mipp[id]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(8).get(0)).fill(-rl.getItem(iis,id).x(),rl.getItem(iis,id).y(),epw<mxp[id]?epw:0);
                }
                }
                }
//            }
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
            c.cd(3*i+j); c.getPad(3*i+j).getAxisY().setLog(getLogY()); 
            c.draw(tl.fitData.getItem(is,i+10*pc*(j+1),0,getRunNumber()).getHist());
            c.draw(tl.fitData.getItem(is,i+10*pc*(j+1),0,getRunNumber()).getGraph(),"same");
        }
        }
        
    }   
    
    public void updateFITS(int index) {
       
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
    
    @Override
    public void plotEvent(DataEvent de) {
    	analyze();
    }

    public void analyze() {    
    	System.out.println(getDetectorName()+".analyze() ");
        fitGraphs(1,7,0,3,0,(dropSummary)?0:3);
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
    
    public void fitStore(int is, int id, int il, int pc, int run, double nrm) {
    	
    	int np = npmt[id*3+il];
        double[]      x = new double[np]; double[]  ymean = new double[np]; double[] yrms = new double[np];
        double[]     xe = new double[np]; double[] ymeane = new double[np]; double[]   ye = new double[np]; 
        double[]  yMean = new double[np]; double[] yMeanc = new double[np]; double[]  ymeanc = new double[np];
        for (int i=0; i<np; i++) {
            x[i] = i+1; xe[i]=0; ye[i]=0; yrms[i]=0; 
            FitData fd = tl.fitData.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); 
            fd.graph.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            fd.hist.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            double mean = fd.mean;                        
            if(mean>0) yrms[i] = fd.sigma/mean; 
            	     yMeanc[i] = (fd.getMean()/nrm)*((pc==1)?shift.getDoubleValue("shift", is, 3*id+il+1, i+1):0);
                     ymeanc[i] =         (mean/nrm)*((pc==1)?shift.getDoubleValue("shift", is, 3*id+il+1, i+1):0);
                      yMean[i] = fd.getMean()/nrm;
                      ymean[i] =         mean/nrm;
                     ymeane[i] = fd.meane/nrm;
        }
        GraphErrors mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,ymean,xe,ymeane);                   
        GraphErrors Mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yMean,xe,ymeane);                   
        GraphErrors meanc = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,ymeanc,xe,ymeane);                   
        GraphErrors Meanc = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yMeanc,xe,ymeane);                   
        GraphErrors  rms = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yrms,xe,ye);                  
        FitSummary.add(mean,  is,id+10*(pc+1)*(pc+1)*(il+1),1,run);
        FitSummary.add(rms,   is,id+10*(pc+1)*(pc+1)*(il+1),2,run);                    
        FitSummary.add(Mean,  is,id+10*(pc+1)*(pc+1)*(il+1),5,run);        	        
        FitSummary.add(meanc, is,id+10*(pc+1)*(pc+1)*(il+1),6,run);        	        
        FitSummary.add(Meanc, is,id+10*(pc+1)*(pc+1)*(il+1),7,run);        	        
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
    
    public void plotXYSummary(int index) {        
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));        
        int          run = getRunNumber();
        H2F h4 = null;
        c.clear(); c.divide(4,3);
        
    	for (int i=0; i<3; i++) {
            h4 = (H2F) this.getDataGroup().getItem(i==0?4:3,2,5,run).getData(i).get(0); 
            for (int j=0; j<4; j++) {int ind1=(j==3)?2:1;int ind2=(j==3)?i:3*i+j;
    	        H2F h1 = (H2F) this.getDataGroup().getItem(0,ind1,5,run).getData(ind2).get(0);   
                H2F h2 = (H2F) this.getDataGroup().getItem(1,ind1,5,run).getData(ind2).get(0); 
                H2F h3 = (H2F) this.getDataGroup().getItem(2,ind1,5,run).getData(ind2).get(0); 
                h3 = h3.divide(h2, h1); 
                if(j==3 && getActivePC()==1) h3 = h3.divide(h4, (H2F) this.getDataGroup().getItem(i==0?4:3,2,5,run).getData(i==0?2:0).get(0)); 
//                if(j==3 && getActivePC()==1) h3 = h3.divide(h4,  (H2F) this.getDataGroup().getItem(0,0,8,run).getData(24).get(0)); 
                h3.setTitle(h1.getName());
                c.cd(4*i+j); c.getPad(4*i+j).getAxisZ().setLog(false); c.getPad(4*i+j).getAxisZ().setRange(zMin*0.03, zMax*0.04);
                c.draw(h3);            
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
    	System.out.println("ECmip.createTimeLineHistos(): Initializing "+TLname+" timeline"); 
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
    	System.out.println("ECmip.saveTimelines()");
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
    		
    		fd = tl.fitData.getItem(is,il,0,getRunNumber());
    		
    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(0.,fd.getHist().getXaxis().max(),0.,fd.getGraph().getMax()*1.1);  
            fd.getHist().getAttributes().setOptStat("1000100");
            DataLine line6 = new DataLine(mipc[il],-50,mipc[il],fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);            
            c.draw(fd.getHist()); c.draw(fd.getGraph(),"same");  c.draw(line6);
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
    		
    		fd = tl.fitData.getItem(is,il+10*(iv+1),0,getRunNumber());
    		
    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(0.,fd.getHist().getXaxis().max(),0.,fd.getGraph().getMax()*1.1);  
            fd.getHist().getAttributes().setOptStat("1000100");
            DataLine line6 = new DataLine(mipp[il],-50,mipp[il],fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);
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
        
        c.clear(); c.divide(3, 2);
    	for (int is=1; is<7; is++) {
    		double min=0.99 ; double max=1.01; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(is-1); c.getPad(is-1).setAxisRange(-0.5,runIndex,min,max); c.getPad(is-1).setTitleFontSize(18);
    		drawTimeLine(c,is,(pc==0)?10*(il+1):3*il+iv+1,1f,"Sector "+is+((pc==0)?" Mean/MIP":l[il]+v[iv]+" Mean/MIP"));
    	}
    }
    
    @Override
	public void writeFile(String table, int is1, int is2, int il1, int il2, int iv1, int iv2) {
		
    	if(!dumpGraphs) return;
    	
		String path = "/Users/colesmith/CLAS12ANA/ECmip/ccdb/";
		String line = new String();
		
		try { 
			File outputFile = new File(path+table+"_"+detcal[0]);
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("ECmip.writeFile("+table+")");

			for (int is=is1; is<is2; is++) {
				for (int il=il1; il<il2; il++ ) {
					for (int iv=iv1; iv<iv2; iv++) {
						for (int ip=0; ip<npmt[3*il+iv]; ip++) {
							switch (table) {
							case "gain": line = getGAIN(is,il,iv,ip,getRunNumber()); break;
							case "rdif": line = getRDIF(is,il,iv,ip,getRunNumber()); break;
							case "rdifgain": line = getRDIFGAIN(is,il,iv,ip,getRunNumber()); 
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
		gain    = eng.engine.getConstantsManager().getConstants(run, "/calibration/ec/gain");
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
	
	public String getRDIF(int is, int il, int iv, int ip, int run) {	
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
	
	public String getRDIFGAIN(int is, int il, int iv, int ip, int run) {	
		gain    = eng.engine.getConstantsManager().getConstants(run, "/calibration/ec/gain");
		shift   = eng.engine.getConstantsManager().getConstants(run, "/calibration/ec/torus_gain_shift");
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
            h2 =dg4.getH2F("hi_"+id+"_path1_"+is);   
            c.cd(is-1); c.getPad(is-1).getAxisZ().setLog(true);       
            c.draw(h2);   
        }
        for (int is=1; is<7; is++) {
            h2 =dg4.getH2F("hi_"+id+"_path2_"+is);   
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
        h2 =dg4.getH2F("hi_pcal_1");   
        c.cd(0); c.getPad(0).getAxisZ().setLog(true); c.draw(h2);
        h2 =dg4.getH2F("hi_ecali_1");   
        c.cd(1); c.getPad(1).getAxisZ().setLog(true); c.draw(h2);
        h2 =dg4.getH2F("hi_ecalo_1");   
        c.cd(2); c.getPad(2).getAxisZ().setLog(true); c.draw(h2);
        h2 =dg4.getH2F("hi_etot_1");   
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
