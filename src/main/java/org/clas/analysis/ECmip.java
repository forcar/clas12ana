package org.clas.analysis;

import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.clas.tools.FitData;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorDescriptor;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.LatexText;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

public class ECmip extends DetectorMonitor {
	
    int is,la,ic,idet,nstr;
    
    float[][][][] ecmean = new float[6][3][3][68];
    float[][][][]  ecrms = new float[6][3][3][68];
    String[]         det = new String[]{"pcal","ecin","ecou"};
    String[]           v = new String[]{"u","v","w"};
    double[]        mipc = {30,30,48};  
    double[]        mipp = {10,10,16};  
    double[]         mxc = {60,60,96};  
    double[]         mxp = {20,20,32};  
    double[]     fitLimp = { 5, 3, 6,17,17,27};
    double[]     fitLimc = {20,17,35,40,48,75};
    int[]           npmt = {68,62,62,36,36,36,36,36,36};
    
    int[]    npmts = new int[]{68,36,36};
        
    IndexedList<GraphErrors>  MIPSummary = new IndexedList<GraphErrors>(4);
    IndexedList<FitData>         MipFits = new IndexedList<FitData>(4);
    IndexedList<GraphErrors>    Timeline = new IndexedList<GraphErrors>(2);
    IndexedList<Float>         PixLength = new IndexedList<Float>(3);
    Boolean                isAnalyzeDone = false;
    
    List<Float>   pmap = new ArrayList<Float>();	
    
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
                                 "PCAL/ECTOT",
                                 "PathIJ",
                                 "PIXEL",
                                 "Timeline",
                                 "ATT");
        
        this.usePCCheckBox(true);
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
    }
    
    public void localinit() {
        configEngine("muon"); 
        getPixLengthMap(outPath+"files/ECpixdepthtotal.dat");
    }  
    
     @Override    
     public void createHistos(int run) {
    	 System.out.println("ECMip:createHistos("+run+")");
	     setRunNumber(run);
	     runlist.add(run);
	     createMIPHistos(0,1,25,0,40," Peak Energy (MeV)");
	     createMIPHistos(0,2,50,0,100," Cluster Energy (MeV)");	     
	     createXYHistos(5,130,420);    
	     createPIDHistos(6);
	     createMIPHistos(7,1,25,0,5.0," Momentum (GeV)");
	     createMIPHistos(7,2,25,0,5.0," Momentum (GeV)");
	     createPathHistos(9);
	     createPixHistos(10);
	     createTimeLineGraphs();
	     createUVWHistos(12,25,0.,2.," MIP ");
     }
     
     @Override       
     public void plotHistos(int run) {
    	     setRunNumber(run);
    	     plotMIP(0);  
    	     if(isAnalyzeDone) {updateUVW(1); updateFITS(2); plotMeanSummary(3); plotRmsSummary(4);plotXYSummary(5);plotTimeline(11);}
    	     plotPIDSummary(6);
    	     plotMIP(7);
    	     plotPathSummary(9);
    	     plotPathSummary(10);
    	     plotUVW(12);
     }
     
     public void createXYHistos(int k, int nb, int bmx) {
    	 
 	     int run = getRunNumber();
         H2F h;  
         
         String[] t = {"e","w","r"};
         
         for (int i=0; i<3; i++) {
        	     DataGroup dg = new DataGroup(3,2);
	    	     for (int d=0; d<3; d++) {
                 h = new H2F("hi_"+det[d]+"_xyc_"+t[i]+"_"+k+"_"+run,"hi_"+det[d]+"_xyc_"+t[i]+"_"+k+"_"+run,nb,-bmx,bmx,nb,-bmx,bmx);
                 dg.addDataSet(h,d);  
	    	     }
             this.getDataGroup().add(dg,i,2,k,run);
         }

         for (int i=0; i<3; i++) {
             DataGroup dg = new DataGroup(3,3);
        	     for (int j=0; j<3; j++) {
        	    	     for (int d=0; d<3; d++) {
                      h = new H2F("hi_"+det[d]+"_xyp_"+v[j]+t[i]+"_"+k+"_"+run,"hi_"+det[d]+"_xyp_"+v[j]+t[i]+"_"+k+"_"+run,nb,-bmx,bmx,nb,-bmx,bmx);
                      dg.addDataSet(h,j+d*3);                      
        	         } 
        	     }
             this.getDataGroup().add(dg,i,1,k,run);
         }
          
     }
     
     public void createUVWHistos(int k, int ybins, double ymin, double ymax, String ytxt) {
     	
         H2F h;  F1D f1; 
         
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
        H2F h; 
    
        for (int is=1; is<7; is++) {
            String tag = is+"_"+n+"_"+k+"_"+run;
            DataGroup dg = new DataGroup(3,3);
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
    	
       DataGroup dg;
	   int run = getRunNumber();
       H2F h; 
       
       dg = new DataGroup(3,4); 
       for (int is=1; is<7; is++) {
           String tag = is+"_"+k+"_"+run;
           h = new H2F("hi_pcal_path1_"+tag,"hi_pcal_path1_"+tag,50,0.,100.,118,31.,50.);
           h.setTitleX("Sector "+is+" PCAL (MeV)");
           h.setTitleY("Path12 (cm)");
           dg.addDataSet(h, is-1);  
           h = new H2F("hi_pcal_path2_"+tag,"hi_pcal_path2_"+tag,50,0.,100.,70,50.,70.);
           h.setTitleX("Sector "+is+" PCAL (MeV)");
           h.setTitleY("Path13 (cm)");
           dg.addDataSet(h, is+5);  
        }
       this.getDataGroup().add(dg,0,0,k,run);
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_ecin_path1_"+tag,"hi_ecin_path1_"+tag,50,0.,100.,70,50.,70.);
            h.setTitleX("Sector "+is+" ECin (MeV)");
            h.setTitleY("Path13 (cm)");
            dg.addDataSet(h, is-1);  
            h = new H2F("hi_ecin_path2_"+tag,"hi_ecin_path2_"+tag,50,0.,100.,66,17.,30.);
            h.setTitleX("Sector "+is+" ECin (MeV)");
            h.setTitleY("Path23 (cm)");
            dg.addDataSet(h, is+5);    
        }
        this.getDataGroup().add(dg,0,1,k,run);
        
        dg = new DataGroup(3,4);         
        for (int is=1; is<7; is++) {        
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_ecou_path1_"+tag,"hi_ecou_path1_"+tag,50,0.,100.,70,50.,70.);
            h.setTitleX("Sector "+is+" ECou (MeV)");
            h.setTitleY("Path13 (cm)");
            dg.addDataSet(h, is-1);  
            h = new H2F("hi_ecou_path2_"+tag,"hi_ecou_path2_"+tag,50,0.,100.,66,17.,30.);
            h.setTitleX("Sector "+is+" ECou (MeV)");
            h.setTitleY("Path23 (cm)");
            dg.addDataSet(h, is+5);      
        }
        this.getDataGroup().add(dg,0,2,k,run);
                
    }
    
    public void createPixHistos(int k) {
    	
        DataGroup dg;
 	    int run = getRunNumber();
        H2F h; 
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            String tag = is+"_"+k+"_"+run;
            h = new H2F("hi_pcal_pix1_"+tag,"hi_pcal_pix1_"+tag,50,1.,200.,12,3.,15.);
            h.setTitleX("Sector "+is+" PCAL (MeV)");
            h.setTitleY("No. Strips");
            dg.addDataSet(h, is-1);   
            h = new H2F("hi_pcal_pix2_"+tag,"hi_pcal_pix2_"+tag,50,31.,50.,20,-25.,25.);
            h.setTitleX("Sector "+is+" Path12 (cm)");
            h.setTitleY("PCAL-ECin (MeV)");
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
            h = new H2F("hi_ecin_pix2_"+tag,"hi_ecin_pix2_"+tag,20,50.,60.,20,-25.,25.);
            h.setTitleX("Sector "+is+" Path13 (cm)");
            h.setTitleY("PCAL-ECou (MeV)");
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
            h = new H2F("hi_ecou_pix2_"+tag,"hi_ecou_pix2_"+tag,30,17.,30.,20,-25,25.);
            h.setTitleX("Sector "+is+" Path23 (cm)");
            h.setTitleY("ECin-ECou (MeV)");
            dg.addDataSet(h, is+5);   
         }
        this.getDataGroup().add(dg,0,2,k,run);
    	
    }
    
    public void createPIDHistos(int k) {
        DataGroup dg = new DataGroup(2,3);
	    int run = getRunNumber();
	    int is  = 0;
        H2F h; 
        String tag = is+"_"+run;
        
        F1D f1 = new F1D("f_1"+tag,"1/(1+[a]^2/x^2)^0.5", 0.41,3.5); f1.setParameter(0,0.13957); f1.setLineColor(1); f1.setLineStyle(1);   
        F1D f2 = new F1D("f_2"+tag,"1/(1+[a]^2/x^2)^0.5", 0.41,3.5); f2.setParameter(0,0.93827); f2.setLineColor(1); f2.setLineStyle(1);   
        h = new H2F("pid_pos_"+tag,"pid_pos_"+tag,100,0.,3.5,100,0.4,1.1);       h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 0); dg.addDataSet(f1,0); dg.addDataSet(f2,0); 
        h = new H2F("pid_neg_"+tag,"pid_neg_"+tag,100,0.,3.5,100,0.4,1.1);       h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 1); dg.addDataSet(f1,1); 
        h = new H2F("pid_fc_pos_"+tag,"pid_fc_pos_"+tag,100,0.,3.5,100,0.4,1.1); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 2); dg.addDataSet(f1,2); dg.addDataSet(f2,2); 
        h = new H2F("pid_fc_neg_"+tag,"pid_fc_neg_"+tag,100,0.,3.5,100,0.4,1.1); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 3); dg.addDataSet(f1,3); 
        h = new H2F("pid_fc_ppi_"+tag,"pid_fc_ppi_"+tag,100,0.,3.5,100,0.4,1.1); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 4); dg.addDataSet(f1,4); dg.addDataSet(f2,4); 
        h = new H2F("pid_fc_npi_"+tag,"pid_fc_npi_"+tag,100,0.,3.5,100,0.4,1.1); h.setTitleX("Momentum (GeV)"); h.setTitleY("BETA");
        dg.addDataSet(h, 5); dg.addDataSet(f1,5); 
        
        this.getDataGroup().add(dg,is,0,k,run);
    }
    
    public void createMiscHistos(int k) {
        DataGroup dg = new DataGroup(2,2);
	    int run = getRunNumber();
        H2F h; 
        String tag = is+"_"+k+"_"+run;
        
        h = new H2F("hi_pcal_1_"+tag,"hi_pcal_1_"+tag,60, 0., 100., 60, 0., 3.);
        dg.addDataSet(h, 0);  
        h = new H2F("hi_ecali_1_"+tag,"hi_ecali_1_"+tag,60, 0., 100., 60, 0., 3.);
        dg.addDataSet(h, 1);  
        h = new H2F("hi_ecalo_1_"+tag,"hi_ecalo_1_"+tag,60, 0., 100., 60, 0., 3.);
        dg.addDataSet(h, 2);  
        h = new H2F("hi_etot_1_"+tag,"hi_etot_1_"+tag,50, 0., 5., 70, 0.05, 0.45);
        dg.addDataSet(h, 3);            
        h = new H2F("hi_pcal_ectot_"+tag,"hi_pcal_ectot_"+tag,50,0.,100.,50,0.,200.);
        h.setTitleX("Sector "+is+" PCAL (MeV)");
        h.setTitleY("ECTOT (MeV)");
        dg.addDataSet(h, 4);  
        h = new H2F("hi_pcal_ectot_max_"+tag,"hi_pcal_ectot_max_"+tag,100,0.,200.,100,0.,300.);
        h.setTitleX("Sector "+is+" PCAL (MeV)");
        h.setTitleY("ECTOT (MeV)");        
        dg.addDataSet(h, 5);
        
        this.getDataGroup().add(dg,0,0,k,run);
    }
    
    public void processEvent(DataEvent event) {
    	
       IndexedList<List<Particle>> ecpart = new IndexedList<List<Particle>>(2);
       List<Particle>                part = new ArrayList<Particle>();
       IndexedList<Vector3>            rl = new IndexedList<Vector3>(2);  
	   Boolean                     isMuon = true;

       Boolean goodPC,goodECi,goodECo;       
       DataBank bank1 = null;
       int nrows = 0;
       Boolean goodEvt  = false, goodRows = false;
    	   
	   int run = getRunNumber();
	   
	   dropBanks(event);
	   
       if (event.hasBank("REC::Particle")) {
    	    isMuon = false;
            DataBank bank = event.getBank("REC::Particle");
            for(int loop = 0; loop < bank.rows(); loop++){
            	    if(bank.getInt("pid",loop)!=0) {
                    Particle p = new Particle(bank.getInt("pid", loop),
                                              bank.getFloat("px", loop),
                                              bank.getFloat("py", loop),
                                              bank.getFloat("pz", loop),
                                              bank.getFloat("vx", loop),
                                              bank.getFloat("vy", loop),
                                              bank.getFloat("vz", loop));
                    p.setProperty("beta",     bank.getFloat("beta", loop));
                    part.add(loop,p);
            	    }
            } 
        }
       
       for (Particle p: part) {
    	       if (p.charge()!=0) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(p.charge()>0?0:1).get(0)).fill(p.p(),p.getProperty("beta"));
       }
       
       Boolean goodREC  = event.hasBank("REC::Calorimeter");
       Boolean goodECAL = event.hasBank("ECAL::clusters")&&event.hasBank("ECAL::calib");
       
       if( isMuon)  goodEvt = goodECAL;
       if(!isMuon&&goodREC&&goodECAL) {goodEvt = true; bank1 = event.getBank("REC::Calorimeter");}
       
       if (goodEvt) {
           DataBank bank2 = event.getBank("ECAL::clusters"); 
           DataBank bank3 = event.getBank("ECAL::calib");
           
           if ( isMuon) {goodRows = true; nrows = bank2.rows();}
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
           
           if(goodRows) {
           
           for(int loop = 0; loop < nrows; loop++){
        	   int ic=0, is=0, il=0;
	           if (!isMuon) {ic = bank1.getShort("index",  loop); is = bank1.getByte("sector",loop); il = bank1.getByte("layer",loop);}
	           if ( isMuon) {ic =loop; is = bank2.getByte("sector",loop); il = bank2.getByte("layer",loop);}          
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
               Particle p = new Particle(); 
               float pm = 0;
               
               if(!isMuon) {
    	       int     ip = bank1.getShort("pindex", loop);
               p.copy(part.get(ip));
               p.setProperty("beta",part.get(ip).getProperty("beta"));
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
               
               pm = (float) p.p();
               
               if (p.charge()>0) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(2).get(0)).fill(pm,p.getProperty("beta"));
               if (p.charge()<0) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(3).get(0)).fill(pm,p.getProperty("beta"));
               
               if (p.pid()==211)  ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(4).get(0)).fill(pm,p.getProperty("beta"));
               if (p.pid()==-211) ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(5).get(0)).fill(pm,p.getProperty("beta"));
               }
               
               Boolean goodPID = isMuon ? true:Math.abs(p.pid())==211;
   	           
               if (!ecpart.hasItem(is,il)) ecpart.add(new ArrayList<Particle>(), is,il);    	                
    	            ecpart.getItem(is,il).add(p);    	   
               
               Vector3 r = new Vector3(x,y,z);
               
               goodPC = goodPID&&il==1&&n1[is-1]<20;  goodECi = goodPID&&il==4&&n4[is-1]<20;  goodECo = goodPID&&il==7&&n7[is-1]<20; 
               
               if (goodPC)  {e1c[is-1][n1[is-1]]=en; rl.add(r,is,0); cU[is-1][0][n1[is-1]]=iU; cV[is-1][0][n1[is-1]]=iV; cW[is-1][0][n1[is-1]]=iW; p1c[is-1][n1[is-1]]=pm;}
               if (goodECi) {e4c[is-1][n4[is-1]]=en; rl.add(r,is,1); cU[is-1][1][n4[is-1]]=iU; cV[is-1][1][n4[is-1]]=iV; cW[is-1][1][n4[is-1]]=iW; p4c[is-1][n4[is-1]]=pm;}
               if (goodECo) {e7c[is-1][n7[is-1]]=en; rl.add(r,is,2); cU[is-1][2][n7[is-1]]=iU; cV[is-1][2][n7[is-1]]=iV; cW[is-1][2][n7[is-1]]=iW; p7c[is-1][n7[is-1]]=pm;}
               if (goodPC)  {p1p[is-1][0][n1[is-1]]=pm; p1p[is-1][1][n1[is-1]]=pm; p1p[is-1][2][n1[is-1]]=pm; w1[is-1][n1[is-1]]=wsum;}
               if (goodECi) {p4p[is-1][0][n4[is-1]]=pm; p4p[is-1][1][n4[is-1]]=pm; p4p[is-1][2][n4[is-1]]=pm; w4[is-1][n4[is-1]]=wsum;}
               if (goodECo) {p7p[is-1][0][n7[is-1]]=pm; p7p[is-1][1][n7[is-1]]=pm; p7p[is-1][2][n7[is-1]]=pm; w7[is-1][n7[is-1]]=wsum;}
               if (goodPC)  {e1p[is-1][0][n1[is-1]]=enu; e1p[is-1][1][n1[is-1]]=env; e1p[is-1][2][n1[is-1]]=enw; n1[is-1]++;}
               if (goodECi) {e4p[is-1][0][n4[is-1]]=enu; e4p[is-1][1][n4[is-1]]=env; e4p[is-1][2][n4[is-1]]=enw; n4[is-1]++;}
               if (goodECo) {e7p[is-1][0][n7[is-1]]=enu; e7p[is-1][1][n7[is-1]]=env; e7p[is-1][2][n7[is-1]]=enw; n7[is-1]++;}
               
           }
           }
             
            for (int is=0; is<6; is++) {
                int iis = is+1;
                if (isGoodTrigger(iis)) {
//                if(n1[is]>=1&&n1[is]<=4&&n4[is]>=1&&n4[is]<=4) { //Cut out vertical cosmic rays
                if(n1[is]==1&n4[is]==1&&n7[is]==1) { //Cut out vertical cosmic rays
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
                     
                Vector3  v1 = new Vector3(rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),rl.getItem(iis,0).z());
                Vector3  v2 = new Vector3(rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),rl.getItem(iis,1).z());
                Vector3  v3 = new Vector3(rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),rl.getItem(iis,2).z());
                Vector3 v23 = new Vector3(v2.x(),v2.y(),v2.z());
        
                v2.sub(v1); v23.sub(v3); v3.sub(v1);  
                
                float v12mag = (float) v2.mag();
                float v13mag = (float) v3.mag();
                float v23mag = (float) v23.mag();
                    
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
               
                
               ((H2F) this.getDataGroup().getItem(0,0,9,run).getData(iis-1).get(0)).fill(e1c[is][0],v12mag);
               ((H2F) this.getDataGroup().getItem(0,0,9,run).getData(iis+5).get(0)).fill(e1c[is][0],v13mag);
               ((H2F) this.getDataGroup().getItem(0,1,9,run).getData(iis-1).get(0)).fill(e4c[is][0],v13mag);
               ((H2F) this.getDataGroup().getItem(0,1,9,run).getData(iis+5).get(0)).fill(e4c[is][0],v23mag);
               ((H2F) this.getDataGroup().getItem(0,2,9,run).getData(iis-1).get(0)).fill(e7c[is][0],v13mag);
               ((H2F) this.getDataGroup().getItem(0,2,9,run).getData(iis+5).get(0)).fill(e7c[is][0],v23mag);
             
               ((H2F) this.getDataGroup().getItem(0,0,10,run).getData(iis-1).get(0)).fill(e1c[is][0],w1[is][0]);
               ((H2F) this.getDataGroup().getItem(0,1,10,run).getData(iis-1).get(0)).fill(e4c[is][0],w4[is][0]);
               ((H2F) this.getDataGroup().getItem(0,2,10,run).getData(iis-1).get(0)).fill(e7c[is][0],w7[is][0]);
               ((H2F) this.getDataGroup().getItem(0,0,10,run).getData(iis+5).get(0)).fill(v12mag,e1c[is][0]-e4c[is][0]);
               ((H2F) this.getDataGroup().getItem(0,1,10,run).getData(iis+5).get(0)).fill(v13mag,e1c[is][0]-e7c[is][0]+18);
               ((H2F) this.getDataGroup().getItem(0,2,10,run).getData(iis+5).get(0)).fill(v23mag,e4c[is][0]-e7c[is][0]+18);

               
                if(v12mag<35&&w1[is][0]==3) {
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
                
                if(v23mag<21&&w4[is][0]==3) {
                for(int n=0; n<n4[is]; n++) {
                	int u = (int) cU[is][1][n]; int v = (int) cV[is][1][n]; int w = (int) cW[is][1][n]; 
                	float d = (PixLength.hasItem(u,v,w))? PixLength.getItem(u,v,w):1f;
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
                    
                	((H2F) this.getDataGroup().getItem(0,  2,5,run).getData(1).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),ec<mxc[1]?mipc[1]:0);
                	((H2F) this.getDataGroup().getItem(1,  2,5,run).getData(1).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),ec<mxc[1]?ec:0); 
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(3).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),epu<mxp[1]?mipp[1]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(3).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),epu<mxp[1]?epu:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(4).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),epv<mxp[1]?mipp[1]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(4).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),epv<mxp[1]?epv:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(5).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),epw<mxp[1]?mipp[1]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(5).get(0)).fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),epw<mxp[1]?epw:0);
                }
                }
                
                if(v23mag<21&&w7[is][0]==3) {                
                for(int n=0; n<n7[is]; n++) {
                	int u = (int) cU[is][2][n]; int v = (int) cV[is][2][n]; int w = (int) cW[is][2][n]; 
                	float d = (PixLength.hasItem(u,v,w))? PixLength.getItem(u,v,w):1f;
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
                    
                	((H2F) this.getDataGroup().getItem(0,  2,5,run).getData(2).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),ec<mxc[2]?mipc[2]:0);
                	((H2F) this.getDataGroup().getItem(1,  2,5,run).getData(2).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),ec<mxc[2]?ec:0); 
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(6).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),epu<mxp[2]?mipp[2]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(6).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),epu<mxp[2]?epu:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(7).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),epv<mxp[2]?mipp[2]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(7).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),epv<mxp[2]?epv:0);
                	((H2F) this.getDataGroup().getItem(0,  1,5,run).getData(8).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),epw<mxp[2]?mipp[2]:0);
                	((H2F) this.getDataGroup().getItem(1,  1,5,run).getData(8).get(0)).fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),epw<mxp[2]?epw:0);
                }
                }
                }
            }
            }
        }
   
    }

    private void updateUVW(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int   off = (getActivePC()==2) ? 0:2;
        int    is = getActiveSector(); 
        int   iis = is+10*off;         
        
        c.clear();
        c.divide(3, 3);

        for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            c.cd(3*i+j); c.getPad(3*i+j).getAxisY().setLog(false); 
            c.draw(MipFits.getItem(iis,i*3+j,0,getRunNumber()).getHist());
            c.draw(MipFits.getItem(iis,i*3+j,0,getRunNumber()).getGraph(),"same");
        }
        }
        
    }
            
    public void updateFITS(int index) {
       
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int   off = (getActivePC()==2) ? 0:2;
        int    is = getActiveSector(); 
        int   iis = is+10*off;         
        int    np = npmt[getActiveLayer()*3+getActiveView()];
        
        c.clear();
        if (np==36) c.divide(6,6); if (np>36) c.divide(8,9);
        
        for (int i=0; i<np ; i++) {
            c.cd(i); c.getPad(i).getAxisY().setLog(false);
            c.draw(MipFits.getItem(iis,getActiveLayer()*3+getActiveView(),i+1,getRunNumber()).getHist());
            c.draw(MipFits.getItem(iis,getActiveLayer()*3+getActiveView(),i+1,getRunNumber()).getGraph(),"same");
       }
       
//        if(isAnalyzeDone) plotMIPSummary(c);       	
    }
    
    @Override
    public void plotEvent(DataEvent de) {
    	    analyze();
    }

    public void analyze() {    
        System.out.println("I am in analyze()");
        analyzeGraphs(1,7,0,3,0,3,"c");
        analyzeGraphs(1,7,0,3,0,3,"p");
        fillTimeLine();
        System.out.println("Finished");
        isAnalyzeDone = true;
    }
    
    public void analyzeGraphs(int is1, int is2, int id1, int id2, int il1, int il2, String ro) {
        
        H2F h2=null;
        FitData fd = null;
        int off=0,ipc=0,run=getRunNumber();
        double min=1,max=20,mip=10;
        System.out.println("Analyzing run "+run);
        for (int is=is1; is<is2; is++) {            
            for (int id=id1; id<id2; id++) {
                if(ro.equals("c")) {min = fitLimc[id]; max = fitLimc[id+3]; off=0; mip=mipc[id]; ipc=2;}
                if(ro.equals("p")) {min = fitLimp[id]; max = fitLimp[id+3]; off=2; mip=mipp[id]; ipc=1;}  
                int iis = is+10*off;
                for (int il=0; il<3; il++) {
                    h2 = (H2F) this.getDataGroup().getItem(is,ipc,0,run).getData(3*id+il).get(0);
//                    h2 = dg4.getH2F("hi_"+det[id]+"_"+lay[il]+ro+"_"+is);
                    fd = new FitData(h2.projectionX().getGraph(),min,max); fd.setInt((int)h2.projectionX().getIntegral()); 
                    fd.setHist(h2.projectionX());
                    fd.graph.getAttributes().setTitleX(h2.getTitleX()); 
                    fd.hist.getAttributes().setTitleX(h2.getTitleX()); 
                    fd.initFit(min,max); fd.fitGraph(""); MipFits.add(fd,iis,id*3+il,0,run);                 
                }                    
                for (int il=il1; il<il2; il++) {
//                    System.out.println("ro:"+ro+" sector "+is+" det "+id+" lay "+il);
                    int np = npmt[id*3+il];
                    double[]  x = new double[np]; double[]  ymean = new double[np]; double[] yrms = new double[np];
                    double[] xe = new double[np]; double[] ymeane = new double[np]; double[]   ye = new double[np]; 
                    double[]  yMean = new double[np];
                    h2 = (H2F) this.getDataGroup().getItem(is,ipc,0,run).getData(id*3+il).get(0);
//                    h2 = dg4.getH2F("hi_"+det[id]+"_"+lay[il]+ro+"_"+is);
                    for (int i=0; i<np; i++) {                     
//                        System.out.println("sector "+is+" det "+id+" lay "+il+" pmt "+i);
                        fd = new FitData(h2.sliceY(i).getGraph(),min,max); fd.setInt((int)h2.sliceY(i).getIntegral()); 
                        fd.setHist(h2.sliceY(i));
                        fd.graph.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
                        fd.hist.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
                        fd.initFit(min,max); fd.fitGraph(""); MipFits.add(fd,iis,id*3+il,i+1,run);
                        x[i] = i+1; xe[i]=0; ye[i]=0; yrms[i]=0;
                        double mean = fd.mean;                        
                        if(mean>0) yrms[i] = fd.sigma/mean; 
                         yMean[i] = fd.getMean()/mip;
                         ymean[i] = mean/mip;
                        ymeane[i] = fd.meane/mip;
                    }
                    GraphErrors mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,ymean,xe,ymeane);                   
                    GraphErrors Mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yMean,xe,ymeane);                   
                    GraphErrors  rms = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yrms,xe,ye);                  
                    MIPSummary.add(mean, 1+off,is,id*3+il,run);
                    MIPSummary.add(rms,  2+off,is,id*3+il,run);                    
                    MIPSummary.add(Mean, 5+off,is,id*3+il,run);
                }
            }
        }
       
    }

    public void plotMeanSummary1(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.divide(3, 6);
        int   off = (getActivePC()==2) ? 0:2;
        int    id = getActiveLayer();
        
        for (int il=0; il<3; il++) {
            F1D f1 = new F1D("p0","[a]",0.,npmt[id*3+il]); f1.setParameter(0,1);
            int n = 3*il;
            for (int is=1; is<7; is++) {
               GraphErrors plot = MIPSummary.getItem(1+off,is,id,il);
               if (is==4) n=9+il*3;
               c.cd(n); c.getPad(n).getAxisY().setRange(0.5, 1.5); 
               c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
               if(n<3||(n>8 && n<12))  plot.getAttributes().setTitle("SECTOR "+is); 
               if(n==0||n==9) plot.getAttributes().setTitleY("MEAN / MIP");
               plot.getAttributes().setTitleX(det[id]+" "+v[il].toUpperCase()+" PMT");
               n++; c.draw(plot);
               f1.setLineColor(3); f1.setLineWidth(3); c.draw(f1,"same");
            }
        }
        
    }
    
    public void plotMeanSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.divide(9, 6);
        int   off = (getActivePC()==2) ? 0:2;        
        int n = 0;
        
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        for (int il=0; il<3; il++) {           	
            F1D f1 = new F1D("p0","[a]",0.,npmt[id*3+il]); f1.setParameter(0,1);
            GraphErrors plot = MIPSummary.getItem(1+off,is,id*3+il,getRunNumber());
            GraphErrors plot1 = MIPSummary.getItem(5+off,is,id*3+il,getRunNumber());
            plot1.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(0.5, 1.5); 
            c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==9||n==18||n==27||n==36||n==45) plot.getAttributes().setTitleY("MEAN / MIP");
            plot.getAttributes().setTitleX("SECTOR "+is+" "+det[id]+" "+v[il].toUpperCase()+" PMT");
            n++; c.draw(plot); c.draw(plot1,"same");
            f1.setLineColor(3); f1.setLineWidth(2); c.draw(f1,"same");
        }
        }
        }
        
    }
    
    public void plotRmsSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.divide(3, 6);
        int   off = (getActivePC()==2) ? 0:2;
        int    id = getActiveLayer();
        
        for (int il=0; il<3; il++) {
            int n = 3*il;
            for (int is=1; is<7; is++) {
               GraphErrors plot = MIPSummary.getItem(2+off,is,id*3+il,getRunNumber());
               if (is==4) n=9+il*3;
               c.cd(n); c.getPad(n).getAxisY().setRange(0.,1.0); 
               c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
               if(n<3||(n>8 && n<12))  plot.getAttributes().setTitle("SECTOR "+is); 
               if(n==0||n==9) plot.getAttributes().setTitleY("RMS / MEAN");
               plot.getAttributes().setTitleX(det[id]+" "+v[il].toUpperCase()+" PMT");
               n++; c.draw(plot);
            }
        }
    }
    
    public void plotXYSummary(int index) {        
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int   run = getRunNumber();
        if(getActivePC()==2) {
        	c.clear(); c.divide(3,2);
	    for (int i=0; i<3; i++) {
            H2F h1 = (H2F) this.getDataGroup().getItem(0,2,5,run).getData(i).get(0); 
            H2F h2 = (H2F) this.getDataGroup().getItem(1,2,5,run).getData(i).get(0);  
            H2F h3 = (H2F) this.getDataGroup().getItem(2,2,5,run).getData(i).get(0);  
            h3 = h2.divide(h2, h1); 
            c.cd(i); c.getPad(i).getAxisZ().setLog(false); c.getPad(i).getAxisZ().setRange(0., 2.);
            c.draw(h3);            
	    }
        }
        
        if(getActivePC()==1) {
        	c.clear(); c.divide(3,3);
    	    for (int i=0; i<3; i++) {
        	    for (int j=0; j<3; j++) {
                H2F h1 = (H2F) this.getDataGroup().getItem(0,1,5,run).getData(3*j+i).get(0);  
                H2F h2 = (H2F) this.getDataGroup().getItem(1,1,5,run).getData(3*j+i).get(0); 
                H2F h3 = (H2F) this.getDataGroup().getItem(2,1,5,run).getData(3*j+i).get(0); 
                h3 = h2.divide(h2, h1); 
                c.cd(3*j+i); c.getPad(3*j+i).getAxisZ().setLog(false); c.getPad(3*j+i).getAxisZ().setRange(0., 2.);
                c.draw(h3);            
        	    }
    	    }
        }      	
        
    }
    
    public void plotMIP(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),getActivePC(),index,getRunNumber()));
    }
        
    public void plotPIDSummary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));
    }
    
    public void plotPathSummary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,getActiveLayer(),index,getRunNumber()));
    }
    
    public void plotUVW(int index) {
        int run = getRunNumber();
        drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index,run));	    
    } 
    
    public void createTimeLineGraphs() {
    	createTimeLineGraph(0,"PCAL MIP Cluster Energy","Sector","Mean/MIP",24);
    	createTimeLineGraph(1,"ECIN MIP Cluster Energy","Sector","Mean/MIP",24);
    	createTimeLineGraph(2,"ECOU MIP Cluster Energy","Sector","Mean/MIP",24);
    }
    
    public void createTimeLineGraph(int k, String tit, String xtit, String ytit, int siz) {  
    	
    	int ic=0;
    	for (int i=0; i<runlist.size(); i++) {
    		if (runlist.get(i)==getRunNumber()) {
            GraphErrors  g = new GraphErrors("MIP"); 
            ic=i+1;if(i>8)ic=i-8;
            g.setTitle(tit); g.setTitleX(xtit); g.setTitleY(ytit); g.setLineColor(ic); g.setMarkerColor(ic); g.setMarkerSize(5);  
            Timeline.add(g,k,i);
    		}
    	}
    }
    
    public void fillTimeLine() {
    	
    	float mip[] = {30,30,48};
    	
    	for (int i=0; i<runlist.size(); i++) {
    		if (runlist.get(i)==getRunNumber()) {
    		for (int is=1; is<7; is++) {
    			double off = -runlist.size()*0.05+i*0.1;
    			for (int id=0; id<3; id++) {
    				float  y = (float) MipFits.getItem(is,3*id,0,getRunNumber()).mean/mip[id];
    				float ye = (float) MipFits.getItem(is,3*id,0,getRunNumber()).meane/mip[id];
    				Timeline.getItem(id,i).addPoint(is+off, y, 0., ye);	
    			}	
    		}
    		}
    	}
    	
    }
    
    public void plotTimeline(int index) {
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
        		c.draw(Timeline.getItem(id,i),i==0?" ":"same");
        		LatexText text = new LatexText("RUN "+runlist.get(i),xlab[i],ylab[i]);
        		text.setFont("HanziPen TC");
        		text.setFontSize(10);
        		ic=i+1;if(i>8)ic=i-8;
        		text.setColor(ic);
        		c.draw(text);
        		}
        	}
        	f1.setLineColor(3); f1.setLineWidth(3); c.draw(f1,"same");	
        }
        
    }
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

    @Override
    public void timerUpdate() {

    } 
  
}
