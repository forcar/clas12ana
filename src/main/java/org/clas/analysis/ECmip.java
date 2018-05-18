package org.clas.analysis;

import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.List;

import org.clas.fcmon.ec.ECGainsApp.FitData;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorDescriptor;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

public class ECmip extends DetectorMonitor {
	
    int is,la,ic,idet,nstr;
    
    float[][][][] ecmean = new float[6][3][3][68];
    float[][][][]  ecrms = new float[6][3][3][68];
    String[]         det = new String[]{"pcal","ecin","ecou"};
    String[]         lay = new String[]{"u","v","w"};
    double[]        mipc = {30,30,48};  
    double[]        mipp = {10,10,16};  
    double[]     fitLimp = { 5, 3, 6,17,17,27};
    double[]     fitLimc = {20,17,35,40,48,75};
    int[]           npmt = {68,62,62,36,36,36,36,36,36};
    
    IndexedList<GraphErrors>  MIPSummary = new IndexedList<GraphErrors>(4);
    IndexedList<FitData>         MipFits = new IndexedList<FitData>(4);
    Boolean                isAnalyzeDone = false;
    
    List<Float>   pmap = new ArrayList<Float>();	
    public ECmip(String name) {
        super(name);
        this.setDetectorTabNames("Summary",
        		                     "Fits",
                                 "Mean/Mip",
                                 "RMS/Mean",
                                 "PCAL/ECTOT",
                                 "PathIJ",
                                 "Maps");
        
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init();
    }
    
    public void OldcreateHistos(int run) {
        
	    setRunNumber(run);
        this.setNumberOfEvents(0); 
         
        DataGroup dg_mip = new DataGroup(1,1);
                
         H2F h2 = new H2F();
         
         int n = 0;
         
         for (int is=1; is<7; is++) {
             h2 = new H2F("hi_pcal_uc_"+is,"hi_pcal_uc_"+is,50, 0., 100., 68, 1., 69.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("coordU");    
             dg_mip.addDataSet(h2, 1); 
             h2 = new H2F("hi_pcal_vc_"+is,"hi_pcal_vc_"+is,50, 0., 100., 62, 1., 63.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("coordV");        
             dg_mip.addDataSet(h2, 1);         
             h2 = new H2F("hi_pcal_wc_"+is,"hi_pcal_wc_"+is,50, 0., 100., 62, 1., 63.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("coordW");  
             dg_mip.addDataSet(h2, 1); 
             
             h2 = new H2F("hi_ecin_uc_"+is,"hi_ecin_uc_"+is,50, 0., 100., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("coordU");        
             dg_mip.addDataSet(h2, 1); 
             h2 = new H2F("hi_ecin_vc_"+is,"hi_ecin_vc_"+is,50, 0., 100., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("coordV");        
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecin_wc_"+is,"hi_ecin_wc_"+is,50, 0., 100., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("coordW");  
             dg_mip.addDataSet(h2, 1);  
             
             h2 = new H2F("hi_ecou_uc_"+is,"hi_ecou_uc_"+is,50, 0., 100., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("coordU");        
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecou_vc_"+is,"hi_ecou_vc_"+is,50, 0., 100., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("coordV");        
             dg_mip.addDataSet(h2, 1); 
             h2 = new H2F("hi_ecou_wc_"+is,"hi_ecou_wc_"+is,50, 0., 100., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("coordW");
             dg_mip.addDataSet(h2, 1); 
             
             h2 = new H2F("hi_pcal_ectot_"+is,"hi_pcal_ectot_"+is,50,0.,100.,50,0.,200.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("ECTOT (MeV)");
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_pcal_ectot_max_"+is,"hi_pcal_ectot_max_"+is,100,0.,200.,100,0.,300.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("ECTOT (MeV)");
             dg_mip.addDataSet(h2, 1);  
             
             h2 = new H2F("hi_pcal_path1_"+is,"hi_pcal_path1_"+is,50,0.,100.,118,31.,90.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("Path12 (cm)");
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_pcal_path2_"+is,"hi_pcal_path2_"+is,50,0.,100.,70,50.,120.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("Path13 (cm)");
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecin_path1_"+is,"hi_ecin_path1_"+is,50,0.,100.,70,50.,120.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("Path13 (cm)");
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecin_path2_"+is,"hi_ecin_path2_"+is,50,0.,100.,66,17.,50.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("Path23 (cm)");
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecou_path1_"+is,"hi_ecou_path1_"+is,50,0.,100.,70,50.,120.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("Path13 (cm)");
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecou_path2_"+is,"hi_ecou_path2_"+is,50,0.,100.,66,17.,50.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("Path23 (cm)");
             dg_mip.addDataSet(h2, 1);  
         }   
         
         for (int is=1; is<7; is++) {
             h2 = new H2F("hi_pcal_up_"+is,"hi_pcal_up_"+is,25, 0., 40., 68, 1., 69.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("coordU");    
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_pcal_vp_"+is,"hi_pcal_vp_"+is,25, 0., 40., 62, 1., 63.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("coordV");        
             dg_mip.addDataSet(h2, 1);            
             h2 = new H2F("hi_pcal_wp_"+is,"hi_pcal_wp_"+is,25, 0., 40., 62, 1., 63.);
             h2.setTitleX("Sector "+is+" PCAL (MeV)");
             h2.setTitleY("coordW");  
             dg_mip.addDataSet(h2, 1);  
             
             h2 = new H2F("hi_ecin_up_"+is,"hi_ecin_up_"+is,25, 0., 40., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("coordU");        
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecin_vp_"+is,"hi_ecin_vp_"+is,25, 0., 40., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("coordV");        
             dg_mip.addDataSet(h2, 1); 
             h2 = new H2F("hi_ecin_wp_"+is,"hi_ecin_wp_"+is,25, 0., 40., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECin (MeV)");
             h2.setTitleY("coordW");  
             dg_mip.addDataSet(h2, 1);  
            
             h2 = new H2F("hi_ecou_up_"+is,"hi_ecou_up_"+is,25, 0., 40., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("coordU");        
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecou_vp_"+is,"hi_ecou_vp_"+is,25, 0., 40., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("coordV");        
             dg_mip.addDataSet(h2, 1);  
             h2 = new H2F("hi_ecou_wp_"+is,"hi_ecou_wp_"+is,25, 0., 40., 36, 1., 37.);
             h2.setTitleX("Sector "+is+" ECou (MeV)");
             h2.setTitleY("coordW");
             dg_mip.addDataSet(h2, 1); 
         } 
         
         DataGroup dg 
         h2 = new H2F("hi_pcal_xyc_e","hi_pcal_xyc_e",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_pcal_xyc_w","hi_pcal_xyc_w",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1); 
         h2 = new H2F("hi_ecin_xyc_e","hi_ecin_xyc_e",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecin_xyc_w","hi_ecin_xyc_w",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1); 
         h2 = new H2F("hi_ecou_xyc_e","hi_ecou_xyc_e",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecou_xyc_w","hi_ecou_xyc_w",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         
         h2 = new H2F("hi_pcal_xyp_ue","hi_pcal_xyp_ue",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_pcal_xyp_uw","hi_pcal_xyp_uw",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1); 
         h2 = new H2F("hi_pcal_xyp_ve","hi_pcal_xyp_ve",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_pcal_xyp_vw","hi_pcal_xyp_vw",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_pcal_xyp_we","hi_pcal_xyp_we",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_pcal_xyp_ww","hi_pcal_xyp_ww",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecin_xyp_ue","hi_ecin_xyp_ue",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecin_xyp_uw","hi_ecin_xyp_uw",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecin_xyp_ve","hi_ecin_xyp_ve",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecin_xyp_vw","hi_ecin_xyp_vw",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecin_xyp_we","hi_ecin_xyp_we",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecin_xyp_ww","hi_ecin_xyp_ww",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecou_xyp_ue","hi_ecou_xyp_ue",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecou_xyp_uw","hi_ecou_xyp_uw",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecou_xyp_ve","hi_ecou_xyp_ve",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecou_xyp_vw","hi_ecou_xyp_vw",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecou_xyp_we","hi_ecou_xyp_we",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecou_xyp_ww","hi_ecou_xyp_ww",100, -500., 500., 100, -500.,500);
         dg_mip.addDataSet(h2, 1);  
         
         h2 = new H2F("hi_pcal_1","hi_pcal_1",60, 0., 100., 60, 0., 3.);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecali_1","hi_ecali_1",60, 0., 100., 60, 0., 3.);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_ecalo_1","hi_ecalo_1",60, 0., 100., 60, 0., 3.);
         dg_mip.addDataSet(h2, 1);  
         h2 = new H2F("hi_etot_1","hi_etot_1",50, 0., 5., 70, 0.05, 0.45);
         dg_mip.addDataSet(h2, 1); 
         
         this.getDataGroup().clear();
         this.getDataGroup().add(dg_mip,4,0,0);        
     }  
    
     public void crateHistos() {
	     createMIPHistos( 1,25, 40,"(MeV)");
	     createMIPHistos(10,50,100,"(MeV)");
	     
	     createPathHistos(5);
    	 
     }
    
     public void createMIPHistos(int k, int nch, double x1, String txt) {
    	
	    int run = getRunNumber();
        H2F h;  
    
        for (int is=1; is<7; is++) {
            DataGroup dg = new DataGroup(3,3);
            h = new H2F("mip_pcal_u_"+is+"_"+k+"_"+run,"mip_pcal_u_"+is+"_"+k+"_"+run, nch, 0., x1, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("mip_pcal_v_"+is+"_"+k+"_"+run,"mip_pcal_v_"+is+"_"+k+"_"+run, nch, 0., x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("mip_pcal_w_"+is+"_"+k+"_"+run,"mip_pcal_w_"+is+"_"+k+"_"+run, nch, 0., x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
        
            h = new H2F("mip_ecin_u_"+is+"_"+k+"_"+run,"mip_ecin_u_"+is+"_"+k+"_"+run, nch, 0., x1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("mip_ecin_v_"+is+"_"+k+"_"+run,"mip_ecin_v_"+is+"_"+k+"_"+run, nch, 0., x1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("mip_ecin_w_"+is+"_"+k+"_"+run,"mip_ecin_w_"+is+"_"+k+"_"+run, nch, 0., x1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
        
            h = new H2F("mip_ecou_u_"+is+"_"+k+"_"+run,"mip_ecou_u_"+is+"_"+k+"_"+run, nch, 0., x1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("mip_ecou_v_"+is+"_"+k+"_"+run,"mip_ecou_v_"+is+"_"+k+"_"+run, nch, 0., x1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("mip_ecou_w_"+is+"_"+k+"_"+run,"mip_ecou_w_"+is+"_"+k+"_"+run, nch, 0., x1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,0,k,run);
        }            
    }  
     
    public void createPathHistos(int k) {
    	
    	   DataGroup dg;
	   int run = getRunNumber();
       H2F h; 
       
       dg = new DataGroup(3,4); 
       for (int is=1; is<7; is++) {
           h = new H2F("hi_pcal_path1_"+is+"_"+k+"_"+run,"hi_pcal_path1_"+is+"_"+k+"_"+run,50,0.,100.,118,31.,90.);
           h.setTitleX("Sector "+is+" PCAL (MeV)");
           h.setTitleY("Path12 (cm)");
           dg.addDataSet(h, is-1);  
           h = new H2F("hi_pcal_path2_"+is+"_"+k+"_"+run,"hi_pcal_path2_"+is+"_"+k+"_"+run,50,0.,100.,70,50.,120.);
           h.setTitleX("Sector "+is+" PCAL (MeV)");
           h.setTitleY("Path13 (cm)");
           dg.addDataSet(h, is+5);  
        }
        this.getDataGroup().add(dg,0,0,k,run);
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {
            h = new H2F("hi_ecin_path1_"+is+"_"+k+"_"+run,"hi_ecin_path1_"+is+"_"+k+"_"+run,50,0.,100.,70,50.,120.);
            h.setTitleX("Sector "+is+" ECin (MeV)");
            h.setTitleY("Path13 (cm)");
            dg.addDataSet(h, is-1);  
            h = new H2F("hi_ecin_path2_"+is+"_"+k+"_"+run,"hi_ecin_path2_"+is+"_"+k+"_"+run,50,0.,100.,66,17.,50.);
            h.setTitleX("Sector "+is+" ECin (MeV)");
            h.setTitleY("Path23 (cm)");
            dg.addDataSet(h, is+5);    
        }
        this.getDataGroup().add(dg,0,1,k,run);
        
        dg = new DataGroup(3,4); 
        for (int is=1; is<7; is++) {        
            h = new H2F("hi_ecou_path1_"+is+"_"+k+"_"+run,"hi_ecou_path1_"+is+"_"+k+"_"+run,50,0.,100.,70,50.,120.);
            h.setTitleX("Sector "+is+" ECou (MeV)");
            h.setTitleY("Path13 (cm)");
            dg.addDataSet(h, is-1);  
            h = new H2F("hi_ecou_path2_"+is+"_"+k+"_"+run,"hi_ecou_path2_"+is+"_"+k+"_"+run,50,0.,100.,66,17.,50.);
            h.setTitleX("Sector "+is+" ECou (MeV)");
            h.setTitleY("Path23 (cm)");
            dg.addDataSet(h, is+5);      
        }
        this.getDataGroup().add(dg,0,2,k,run);
        
    }
    
    public void addEvent(DataEvent event) {
        
        pmap.clear();

        if (event.hasBank("RECHB::Particle")) {
            DataBank bank = event.getBank("RECHB::Particle");
            for(int loop = 0; loop < bank.rows(); loop++){
                float px = bank.getFloat("px", loop);
                float py = bank.getFloat("py", loop);
                float pz = bank.getFloat("pz", loop);
                int q    = bank.getByte("charge", loop);                
                if (q!=0) pmap.add((float) Math.sqrt(px*px+py*py+pz*pz));
            } 
        }
        
        if (event.hasBank("MIP::event")) {
            DataBank bank = event.getBank("MIP::event");
            for (int loop = 0; loop < bank.rows(); loop++) {                
                if (bank.getInt("q",loop)<0) pmap.add((float)bank.getFloat("p", loop));
            }                   
        }
/*            
        if(recPartEB!=null && recDeteEB!=null) {
            int nrows = recPartEB.rows();
            for(int loop = 0; loop < nrows; loop++){
                int pidCode = 0;
                if(recPartEB.getByte("charge", loop)==-1) pidCode = -211;
                if(recPartEB.getByte("charge", loop)==+1) pidCode = +211;
                Boolean pidPion = pidCode==-211 || pidCode==+211;
                System.out.println(pidCode);
                if(pidPion) {
                recParticle = new Particle(pidCode,
                  recPartEB.getFloat("px", loop),
                  recPartEB.getFloat("py", loop),
                  recPartEB.getFloat("pz", loop),
                  recPartEB.getFloat("vx", loop),
                  recPartEB.getFloat("vy", loop),
                  recPartEB.getFloat("vz", loop));
                
                double energy1=0;
                double energy4=0;
                double energy7=0;

                for(int j=0; j<recDeteEB.rows(); j++) {
                    if(recDeteEB.getShort("pindex",j)==loop && recDeteEB.getShort("detector",j)==16) {
                        if(energy1 >= 0 && recDeteEB.getShort("layer",j) == 1) energy1 += recDeteEB.getFloat("energy",j);
                        if(energy4 >= 0 && recDeteEB.getShort("layer",j) == 4) energy4 += recDeteEB.getFloat("energy",j);
                        if(energy7 >= 0 && recDeteEB.getShort("layer",j) == 7) energy7 += recDeteEB.getFloat("energy",j);
                    }
                }
                
                recParticle.setProperty("energy1",energy1);
                recParticle.setProperty("energy4",energy4);
                recParticle.setProperty("energy7",energy7);
                
                if(partRecEB==null && pidPion) {
                    recParticle.setProperty("sector",recBankTB.getByte("sector", loop)*1.0);
                    partRecEB=recParticle;
                }
                }
            }
        }
*/        
//        event.removeBank("ECAL::clusters");
//        engine.processDataEvent(event);
          if (app.isSingleEvent()) {
              event.show();
              if(event.hasBank("ECAL::adc")) event.getBank("ECAL::adc").show();
              if(event.hasBank("ECAL::hits")) event.getBank("ECAL::hits").show();
          }
//        event.removeBank("ECAL::clusters");
//        event.show();
        
        // EC clusters
        Boolean goodPC,goodECi,goodECo;
        IndexedList<Vector3> rl = new IndexedList<Vector3>(2);
        if(event.hasBank("ECAL::clusters") && event.hasBank("ECAL::calib")){
            DataBank  bank1 = event.getBank("ECAL::clusters");
            DataBank  bank2 = event.getBank("ECAL::calib");
            int[] n1 = new int[6]; int[] n4 = new int[6]; int[] n7 = new int[6];
            float[][]  e1c = new float[6][20]; float[][][]   e1p = new float[6][3][20]; 
            float[][]  e4c = new float[6][20]; float[][][]   e4p = new float[6][3][20]; 
            float[][]  e7c = new float[6][20]; float[][][]   e7p = new float[6][3][20]; 
            float[][][] cU = new float[6][3][20];
            float[][][] cV = new float[6][3][20];
            float[][][] cW = new float[6][3][20];
            if(bank1.rows()==bank2.rows()) {
            int rows = bank1.rows();
            for(int loop = 0; loop < rows; loop++){
                int   is = bank1.getByte("sector", loop);
                int   il = bank1.getByte("layer", loop);
                float en = bank1.getFloat("energy",loop)*1000;
                float  x = bank1.getFloat("x", loop);
                float  y = bank1.getFloat("y", loop);
                float  z = bank1.getFloat("z", loop);
                int   iU = (bank1.getInt("coordU", loop)-4)/8;
                int   iV = (bank1.getInt("coordV", loop)-4)/8;
                int   iW = (bank1.getInt("coordW", loop)-4)/8;
                float enu = bank2.getFloat("recEU",loop)*1000;
                float env = bank2.getFloat("recEV",loop)*1000;
                float enw = bank2.getFloat("recEW",loop)*1000;
                Vector3 r = new Vector3(x,y,z);
               
                goodPC = il==1&&n1[is-1]<20;  goodECi = il==4&&n4[is-1]<20;  goodECo = il==7&&n7[is-1]<20; 
                if (goodPC)  {e1c[is-1][n1[is-1]]=en; rl.add(r,is,0); cU[is-1][0][n1[is-1]]=iU; cV[is-1][0][n1[is-1]]=iV; cW[is-1][0][n1[is-1]]=iW;}
                if (goodECi) {e4c[is-1][n4[is-1]]=en; rl.add(r,is,1); cU[is-1][1][n4[is-1]]=iU; cV[is-1][1][n4[is-1]]=iV; cW[is-1][1][n4[is-1]]=iW;}
                if (goodECo) {e7c[is-1][n7[is-1]]=en; rl.add(r,is,2); cU[is-1][2][n7[is-1]]=iU; cV[is-1][2][n7[is-1]]=iV; cW[is-1][2][n7[is-1]]=iW;}
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
                    
                DataGroup dg4 = this.getDataGroup().getItem(4,0,0);
                
                double ectot = e4c[is][0]+e7c[is][0] ; double etot = e1c[is][0]+ectot ;
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
                
                dg4.getH2F("hi_pcal_path1_"+iis).fill(e1c[is][0],v12mag);
                dg4.getH2F("hi_pcal_path2_"+iis).fill(e1c[is][0],v13mag);
                dg4.getH2F("hi_ecin_path1_"+iis).fill(e4c[is][0],v13mag);
                dg4.getH2F("hi_ecin_path2_"+iis).fill(e4c[is][0],v23mag);
                dg4.getH2F("hi_ecou_path1_"+iis).fill(e7c[is][0],v13mag);
                dg4.getH2F("hi_ecou_path2_"+iis).fill(e7c[is][0],v23mag);
                
                if(v12mag<34) {
                for(int n=0; n<n1[is]; n++) {
                		dg4.getH2F("hi_pcal_uc_"+iis).fill(e1c[is][n],cU[is][0][n]);
                		dg4.getH2F("hi_pcal_vc_"+iis).fill(e1c[is][n],cV[is][0][n]);
                		dg4.getH2F("hi_pcal_wc_"+iis).fill(e1c[is][n],cW[is][0][n]);
                		dg4.getH2F("hi_pcal_up_"+iis).fill(e1p[is][0][n],cU[is][0][n]);
                		dg4.getH2F("hi_pcal_vp_"+iis).fill(e1p[is][1][n],cV[is][0][n]);
                		dg4.getH2F("hi_pcal_wp_"+iis).fill(e1p[is][2][n],cW[is][0][n]);
                		dg4.getH2F("hi_pcal_xyc_e").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),mipc[0]);
                		dg4.getH2F("hi_pcal_xyc_w").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),e1c[is][n]);
                		dg4.getH2F("hi_pcal_xyp_ue").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),mipp[0]);
                		dg4.getH2F("hi_pcal_xyp_uw").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),e1p[is][0][n]);
                		dg4.getH2F("hi_pcal_xyp_ve").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),mipp[0]);
                		dg4.getH2F("hi_pcal_xyp_vw").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),e1p[is][1][n]);
                		dg4.getH2F("hi_pcal_xyp_we").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),mipp[0]);
                		dg4.getH2F("hi_pcal_xyp_ww").fill(-rl.getItem(iis,0).x(),rl.getItem(iis,0).y(),e1p[is][2][n]);
                }
                }
                
                if(v23mag<20) {
                for(int n=0; n<n4[is]; n++) {
                		dg4.getH2F("hi_ecin_uc_"+iis).fill(e4c[is][n],cU[is][1][n]);
                		dg4.getH2F("hi_ecin_vc_"+iis).fill(e4c[is][n],cV[is][1][n]);
                		dg4.getH2F("hi_ecin_wc_"+iis).fill(e4c[is][n],cW[is][1][n]);
                		dg4.getH2F("hi_ecin_up_"+iis).fill(e4p[is][0][n],cU[is][1][n]);
                		dg4.getH2F("hi_ecin_vp_"+iis).fill(e4p[is][1][n],cV[is][1][n]);
                		dg4.getH2F("hi_ecin_wp_"+iis).fill(e4p[is][2][n],cW[is][1][n]);
                		dg4.getH2F("hi_ecin_xyc_e").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),mipc[1]);
                		dg4.getH2F("hi_ecin_xyc_w").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),e4c[is][n]);
                		dg4.getH2F("hi_ecin_xyp_ue").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),mipp[1]);
                		dg4.getH2F("hi_ecin_xyp_uw").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),e4p[is][0][n]);
                		dg4.getH2F("hi_ecin_xyp_ve").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),mipp[1]);
                		dg4.getH2F("hi_ecin_xyp_vw").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),e4p[is][1][n]);
                		dg4.getH2F("hi_ecin_xyp_we").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),mipp[1]);
                		dg4.getH2F("hi_ecin_xyp_ww").fill(-rl.getItem(iis,1).x(),rl.getItem(iis,1).y(),e4p[is][2][n]);
                }
                for(int n=0; n<n7[is]; n++) {
                		dg4.getH2F("hi_ecou_uc_"+iis).fill(e7c[is][n],cU[is][2][n]);
                		dg4.getH2F("hi_ecou_vc_"+iis).fill(e7c[is][n],cV[is][2][n]);
                		dg4.getH2F("hi_ecou_wc_"+iis).fill(e7c[is][n],cW[is][2][n]);
                		dg4.getH2F("hi_ecou_up_"+iis).fill(e7p[is][0][n],cU[is][2][n]);
                		dg4.getH2F("hi_ecou_vp_"+iis).fill(e7p[is][1][n],cV[is][2][n]);
                		dg4.getH2F("hi_ecou_wp_"+iis).fill(e7p[is][2][n],cW[is][2][n]);
                		dg4.getH2F("hi_ecou_xyc_e").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),mipc[2]);
                		dg4.getH2F("hi_ecou_xyc_w").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),e7c[is][n]);
                		dg4.getH2F("hi_ecou_xyp_ue").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),mipp[2]);
                		dg4.getH2F("hi_ecou_xyp_uw").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),e7p[is][0][n]);
                		dg4.getH2F("hi_ecou_xyp_ve").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),mipp[2]);
                		dg4.getH2F("hi_ecou_xyp_vw").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),e7p[is][1][n]);
                		dg4.getH2F("hi_ecou_xyp_we").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),mipp[2]);
                		dg4.getH2F("hi_ecou_xyp_ww").fill(-rl.getItem(iis,2).x(),rl.getItem(iis,2).y(),e7p[is][2][n]);
                }
                }
                }
            }
            }
        }
   
    }
        
    private class FitData {
        
        GraphErrors graph = null;

        public double xmin;
        public double xmax;
        public double amp;
        public double mean;
        public double meane;
        public double sigma;
        public int integral;
        public int intmin=30;
        public int fitcol=4;
        public String optstat="1100";
        Boolean doFit = false;
        
        FitData(GraphErrors graph, double xmin, double xmax) {
            setGraph(graph);
            this.graph.setFunction(new F1D("f","[amp]*gaus(x,[mean],[sigma])", xmin, xmax));            
            this.graph.getFunction().setLineWidth(1);
        };
        
        public void setGraph(GraphErrors graph) {
            this.graph = graph;
        }
        
        public void setInt(int integral) {
            this.integral = integral;
            fitcol = 4;
            if(integral>intmin) {
                doFit = true;
                fitcol = 1;
            }
        }
        
        public void setIntMin(int intmin) {
            this.intmin = intmin;
        }
        
        public void initFit(double xmin, double xmax) {
            this.xmin = xmin;
            this.xmax = xmax;
            amp   = getMaxYIDataSet(graph,xmin,xmax);
            mean  = getMeanIDataSet(graph,xmin,xmax); 
            sigma = getRMSIDataSet(graph,xmin,xmax);
            graph.getFunction().setParameter(0, amp);
            graph.getFunction().setParameter(1, mean);
            graph.getFunction().setParameter(2, sigma);
            graph.getFunction().setRange(mean-2.5*sigma, mean+1.7*sigma);
         }
        
        public void fitGraph(String opt) {
            if (doFit) DataFitter.fit(graph.getFunction(), graph, opt);
            amp   = graph.getFunction().getParameter(0);
            mean  = graph.getFunction().parameter(1).value();
            meane = graph.getFunction().parameter(1).error();
            sigma = graph.getFunction().getParameter(2);   
            graph.getFunction().setLineColor(fitcol);
            graph.getFunction().setOptStat(opt=="Q"?"0":optstat);
        }
        
        public void plotGraph(EmbeddedCanvas c, int col) {
            graph.getFunction().setLineColor(col); 
            c.draw(graph);
        }
        
        public GraphErrors getGraph() {
            return graph;
        }
        
    }
    
    public void updateCanvas(DetectorDescriptor dd) {
        if (!app.doEng) return;
        this.is = dd.getSector();
        this.la = dd.getLayer();
        this.ic = dd.getComponent();   
        this.idet = ilmap;      
        
        if (la>3) return;
               
        updatePlots();
                
     }
    
    private void updatePlots() {        
        updateMIP(peaks);
        updateMIP(clusters);
        updateSummary();
        updatePvsE();
    }
    
    private void updateMIP(EmbeddedCanvasTabbed c) {
        
        DataGroup dg4 = this.getDataGroup().getItem(4,0,0);
        
        EmbeddedCanvas cp = c.getCanvas(c.getName());
        EmbeddedCanvas cf = c.getCanvas("Fits");
        String d = c.getName().substring(0,1).toLowerCase();
        int off = (d.equals("c")) ? 0:2;
        int    is = getActiveSector(); 
        int    iis=is+10*off;
        String il = lay[getActiveLayer()];
        String id = det[getActiveLayer()];        
        int    np = npmt[getActiveLayer()*3+getActiveView()];
        
        if(isAnalyzeDone==false||app.getInProcess()<2) return;        
//        if(app.getInProcess()==1) analyzeGraphs(is,is+1,activeDetector,activeDetector+1,activeLayer,activeLayer+1,d);
        
        cp.clear(); cp.divide(3, 2);
        
        for (int i=0; i<3; i++) {
            cp.cd(i); cp.getPad(i).getAxisZ().setLog(true);
            cp.draw(dg4.getH2F("hi_"+id+"_"+lay[i]+d+"_"+is));
            cp.cd(i+3); cp.getPad(i+3).getAxisY().setLog(false); 
            cp.draw(MipFits.getItem(iis,getActiveLayer(),i,0).getGraph());
        }
        
        cf.clear(); if (np==36) cf.divide(6,6); if (np>36)  cf.divide(8,9);
        
        for (int i=0; i<np ; i++) {
            cf.cd(i); cf.getPad(i).getAxisY().setLog(false);  
            cf.draw(MipFits.getItem(iis,getActiveLayer(),getActiveView(),i+1).getGraph());
        }
        
        if(isAnalyzeDone&&app.getInProcess()>1) plotMIPSummary(c);           
    }

    public void analyze() {    
        System.out.println("I am in analyze()");
        analyzeGraphs(1,7,0,3,0,3,"c");
        analyzeGraphs(1,7,0,3,0,3,"p");
        System.out.println("Finished");
        isAnalyzeDone = true;
    }
    
    public void analyzeGraphs(int is1, int is2, int id1, int id2, int il1, int il2, String ro) {
        
        DataGroup dg4 = this.getDataGroup().getItem(4,0,0);

        H2F h2=null;
        FitData fd = null;
        int off=0;
        double min=1,max=20,mip=10;
        System.out.println("I am in analyzeGraphs");
        for (int is=is1; is<is2; is++) {            
            for (int id=id1; id<id2; id++) {
                if(ro.equals("c")) {min = fitLimc[id]; max = fitLimc[id+3]; off=0; mip=mipc[id];}
                if(ro.equals("p")) {min = fitLimp[id]; max = fitLimp[id+3]; off=2; mip=mipp[id];}  
                int iis = is+10*off;
                for (int il=0; il<3; il++) {
                    h2 =dg4.getH2F("hi_"+det[id]+"_"+lay[il]+ro+"_"+is);
                    fd = new FitData(h2.projectionX().getGraph(),min,max); fd.setInt((int)h2.projectionX().getIntegral()); 
                    fd.graph.getAttributes().setTitleX(h2.getTitleX()); 
                    fd.initFit(min,max); fd.fitGraph("Q"); MipFits.add(fd,iis,id,il,0);                    
                }                    
                for (int il=il1; il<il2; il++) {
                    System.out.println("ro:"+ro+" sector "+is+" det "+id+" lay "+il);
                    int np = npmt[id*3+il];
                    double[]  x = new double[np]; double[]  ymean = new double[np]; double[] yrms = new double[np];
                    double[] xe = new double[np]; double[] ymeane = new double[np]; double[]   ye = new double[np]; 
                    h2 =dg4.getH2F("hi_"+det[id]+"_"+lay[il]+ro+"_"+is);
                    for (int i=0; i<np; i++) {                     
                        System.out.println("sector "+is+" det "+id+" lay "+il+" pmt "+i);
                        fd = new FitData(h2.sliceY(i).getGraph(),min,max); fd.setInt((int)h2.sliceY(i).getIntegral()); 
                        fd.graph.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+lay[il]+(i+1));
                        fd.initFit(min,max); fd.fitGraph("Q"); MipFits.add(fd,iis,id,il,i+1);
                        x[i] = i+1; xe[i]=0; ye[i]=0; yrms[i]=0;
                        double mean = fd.mean;                        
                        if(mean>0) yrms[i] = fd.sigma/mean; 
                         ymean[i] = mean/mip;
                        ymeane[i] = fd.meane/mip;
                    }
                    GraphErrors mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,ymean,xe,ymeane);                   
                    GraphErrors  rms = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yrms,xe,ye);                  
                    MIPSummary.add(mean, 1+off,is,id,il);
                    MIPSummary.add(rms,  2+off,is,id,il);                    
                }
            }
        }
        
    }
    
    public void plotMIPSummary(EmbeddedCanvasTabbed ct) {
        
        DataGroup dg4 = this.getDataGroup().getItem(4,0,0);
        EmbeddedCanvas c = ct.getCanvas("Mean/Mip"); c.divide(3, 6);
        String        ro = ct.getName().substring(0,1).toLowerCase();
        int          off = (ro.equals("c")) ? 0:2;        
        int           id = getActiveLayer();
        int          iil = getActiveView();
        String[]     det = {"PCAL","ECin","ECou"};

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
               plot.getAttributes().setTitleX(det[id]+" "+lay[il].toUpperCase()+" PMT");
               n++; c.draw(plot);
               f1.setLineColor(3); f1.setLineWidth(3); c.draw(f1,"same");
            }
        }
        
        c = ct.getCanvas("RMS/Mean"); c.divide(3, 6);

        for (int il=0; il<3; il++) {
            int n = 3*il;
            for (int is=1; is<7; is++) {
               GraphErrors plot = MIPSummary.getItem(2+off,is,id,il);
               if (is==4) n=9+il*3;
               c.cd(n); c.getPad(n).getAxisY().setRange(0.,1.0); 
               c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
               if(n<3||(n>8 && n<12))  plot.getAttributes().setTitle("SECTOR "+is); 
               if(n==0||n==9) plot.getAttributes().setTitleY("RMS / MEAN");
               plot.getAttributes().setTitleX(det[id]+" "+lay[il].toUpperCase()+" PMT");
               n++; c.draw(plot);
            }
        }
        
        c = ct.getCanvas("Maps"); c.cd(0) ; 
        c.getPad(0).getAxisZ().setLog(false); c.getPad(0).getAxisZ().setRange(0., 2.);  
        if(off==0) {       
            H2F h1 =dg4.getH2F("hi_"+det[id].toLowerCase()+"_"+"xyc_e");
            H2F h2 =dg4.getH2F("hi_"+det[id].toLowerCase()+"_"+"xyc_w");
            H2F h3 = h2.divide(h2, h1); c.draw(h3);
        } else {
            H2F h1 =dg4.getH2F("hi_"+det[id].toLowerCase()+"_"+"xyp_"+lay[iil]+"e");
            H2F h2 =dg4.getH2F("hi_"+det[id].toLowerCase()+"_"+"xyp_"+lay[iil]+"w");            
            H2F h3 = h2.divide(h2, h1); c.draw(h3);
        }
        
        c.repaint();
        
    }
    
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
    
    private double getMaxYIDataSet(IDataSet data, double min, double max) {
        double max1 = 0;
        double xMax = 0;
        for (int i = 0; i < data.getDataSize(0); i++) {
            double x = data.getDataX(i);
            double y = data.getDataY(i);
            if (x > min && x < max && y != 0) {
                if (y > max1) {
                    max1 = y;
                    xMax = x;
                }
            }
        }
        return max1;
    }
    
    private double getMeanIDataSet(IDataSet data, double min, double max) {
        int nsamples = 0;
        double sum = 0;
        double nEntries = 0;
        for (int i = 0; i < data.getDataSize(0); i++) {
            double x = data.getDataX(i);
            double y = data.getDataY(i);
            if (x > min && x < max && y != 0) {
                nsamples++;
                sum += x * y;
                nEntries += y;
            }
        }
        return sum / (double) nEntries;
    }
    
    private double getRMSIDataSet(IDataSet data, double min, double max) {
        int nsamples = 0;
        double mean = getMeanIDataSet(data, min, max);
        double sum = 0;
        double nEntries = 0;

        for (int i = 0; i < data.getDataSize(0); i++) {
            double x = data.getDataX(i);
            double y = data.getDataY(i);
            if (x > min && x < max && y != 0) {
                nsamples++;
                sum += Math.pow(x - mean, 2) * y;
                nEntries += y;
            }
        }
        return Math.sqrt(sum / (double) nEntries);
    }   
    
}
