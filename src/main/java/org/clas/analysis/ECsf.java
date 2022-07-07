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
import org.clas.tools.ParallelSliceFitter;
import org.clas.tools.SFFunction;
import org.clas.viewer.DetectorMonitor;
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
import org.jlab.utils.groups.IndexedTable;

public class ECsf extends DetectorMonitor {

	IndexedTable electron_sf = null;
	float[][] sfpar = new float[7][2];
	
    String[]  det = {"pcal","ecin","ecou"};
    int[]    npmt = {68,62,62,36,36,36,36,36,36};    
    String[]    v = new String[]{"u","v","w"};
    float     EB = 0;
    boolean isMC = false;
    int nelec=0, trigger_sect=0;
    
    public double[][] par = {{0.105,0.039},{0.099,0.040},{0.100,0.034},{0.093,0.044},{0.085,0.046},{0.113,0.028}};
    
    public ECsf(String name) {
        super(name);
        this.setDetectorTabNames("E/P",
                				 "SLC",
                				 "UVW",
                				 "Timing",
                                 "XY",
                                 "PMT Fits",
                                 "Summary",
                                 "PID Fits",
                                 "Timeline",
                                 "MC");

        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECsf.localinit()");
//        eng.config("muon");  
        tl.setFitData(Fits);
    }
    
    public void localclear() {
    	System.out.println("ECsf.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
    @Override
    public void createHistos(int run) {        
        histosExist = true;
        System.out.println("ECsf.createHistos("+run+")");
        readSF(outPath+"ECsf/electron_sf/electron_sf_"+run);
        EB = getBeamEnergy(run);
        System.out.println("ECsf: EB="+EB);
	    DataGroup dg = null;
	    
        setRunNumber(run);
        runlist.add(run);
        this.setNumberOfEvents(0);
        dg = new DataGroup(6,4);
        createSecHistos("E/P",0,0,50,0.0,EB*0.25,50,0.12,0.35,"ep_emf", " Measured Energy (GeV)", " E/P",dg);
        createSecHistos("E/P",0,1,50,0.2,EB,  50,0.1,0.35,"ep_pf",  " Momentum (GeV)",   " E/P",dg);
        createSecHistos("E/P",0,2,30,  3.,35.,50,0.1,0.35,"ep_thvf"," VertexTheta (deg)"," E/P",dg);
        createSecHistos("E/P",0,3,48,  3.,35.,50,0.1,0.35,"ep_thdf"," Detector Theta (deg)"," E/P",dg);
        dg = new DataGroup(6,4);
        createSecHistos("E/P",1,0,50,0.0,EB*0.25,50,0.12,0.35,"ep_em", " Measured Energy (GeV)", " E/P",dg);
        createSecHistos("E/P",1,1,50,0.2,EB,  50,0.1,0.35,"ep_p"," Momentum (GeV)",   " E/P",dg);
        createSecHistos("E/P",1,2,30,  3.,35.,50,0.1,0.35,"ep_thv"," VertexTheta (deg)"," E/P",dg);
        createSecHistos("E/P",1,3,48,  3.,35.,50,0.1,0.35,"ep_thd"," Detector Theta (deg)"," E/P",dg);
        createXYZHistos("XY","xyz"); createXYZHistos("XY","hxyz"); 
        createSLCHistos("SLC",2,0,15,0.,0.30/3,"PEAK PARTIAL SF",false);
        createSLCHistos("SLC",2,1,12,0.,0.15/3,"PEAK PARTIAL SF",false);
        createSLCHistos("SLC",2,2,12,0.,0.06/3,"PEAK PARTIAL SF",false);
        createSLCHistos("SLC",1,0,15,0.,0.30,"CLUSTER PARTIAL SF",false);
        createSLCHistos("SLC",1,1,12,0.,0.15,"CLUSTER PARTIAL SF",false);
        createSLCHistos("SLC",1,2,12,0.,0.05,"CLUSTER PARTIAL SF",false);
        createSLCHistos("SLC",0,0,25,0.1,0.35,"TOTAL SF",false);
        createSLCHistos("SLC",0,1,25,0.1,0.35,"TOTAL SF",false);
        createSLCHistos("SLC",0,2,25,0.1,0.35,"TOTAL SF",false);
        createSLCHistos("Timing",0,0,50,-5.,5.," T-TVERT-PATH/c (ns)",true);
        createSLCHistos("Timing",0,1,50,-5.,5.," T-TVERT-PATH/c (ns)",true);
        createSLCHistos("Timing",0,2,50,-5.,5.," T-TVERT-PATH/c (ns)",true);        
        createUVWHistos("UVW",0,25,25,0,EB,0.01,0.30/3,"P (GEV)","SF ");
        createUVWHistos("UVW",1,25,25,0,EB,0,0.15/3,"P (GEV)","SF ");
        createUVWHistos("UVW",2,25,25,0,EB,0,0.06/3,"P (GEV)","SF ");
        createPIDFitHistos("PID Fits",50,0.1,0.4," "," Corrected E/P");
        createMCHistos("MC");
    }

    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);
    	plotAnalysis(run);
    }
      
    public void plotSummary(int run) {
    	setRunNumber(run);
    	if(dropSummary) return;
    	plotSECHistos("E/P"); 
    	plotSLCHistos("SLC");
    	if(!dropSummary) plotXYZHistos("XY");
    	plotSLCHistos("Timing");
    	plotUVWHistos("UVW");
    	plotMCHistos("MC");    	
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;
    	plotPIDFits("PID Fits");
    	if(!dropSummary) {updateFits("PMT Fits");plotMeanHWSummary("Summary");}
    	plotTimeLines("Timeline");
    	if (dumpGraphs) dumpGraphs();
    }
    
    public void createPIDFitHistos(String tab,int xb, double x1, double x2, String txt1, String txt2) {
    	 int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
    	 H1F h;
    	 
    	 int n=0;
    	 DataGroup dg = new DataGroup(6,1);
    	 
    	 for (int is=1; is<7; is++) {
           h = new H1F(txt1+"_"+n+"_"+k+"_"+run, xb, x1, x2);
           h.setTitleX("Sector "+is+txt2);
           dg.addDataSet(h,n);n++;
         }     		 
         this.getDataGroup().add(dg,0,0,k,run);

    }
    
    public void createSecHistos(String tab, int pc, int st, int xb, double x1, double x2, int yb, double y1, double y2, String txt1, String txt2, String txt3, DataGroup dg) {
    	
        int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
        H2F h;  
        
        int n=st*6;
       
        for (int is=1; is<7; is++) {
            h = new H2F(txt1+"_"+n+"_"+k+"_"+run, xb, x1, x2, yb, y1, y2);
//            h.setTitle(txt1);
            h.setTitleX("Sector "+is+txt2);
            h.setTitleY(txt3);
            dg.addDataSet(h,n);n++;
        } 
        this.getDataGroup().add(dg,pc,0,k,run);

    } 
    
    public void createXYZHistos(String tab, String xyz) {
    	    
        int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);

        H2F h;
        DataGroup dg1 = new DataGroup(3,3);
        DataGroup dg2 = new DataGroup(3,3);
        DataGroup dg3 = new DataGroup(3,3);
        
        String tit1[] = {"PCAL EVENTS","ECIN EVENTS","ECOU EVENTS"};
        String tit2[] = {"PCAL PARTIAL SF","ECIN PARTIAL SF","ECOU PARTIAL SF"};
        String tit3[] = {"TOTAL SF PCAL","TOTAL SF ECIN","TOTAL SF ECOU"};
        
        
        int nb=200, x=400, y=400;
        
        for (int i=1; i<4; i++) {
            h = new H2F("ep_"+xyz+"_w"+i+"_"+run,  nb, -x, x, nb, -y,y);                        dg2.addDataSet(h,i-1); 
            h = new H2F("ep_"+xyz+"_ww"+i+"_"+run, nb, -x, x, nb, -y,y);                        dg3.addDataSet(h,i-1); 
            h = new H2F("ep_"+xyz+"_e"+i+"_"+run,  nb, -x, x, nb, -y,y); h.setTitle(tit1[i-1]); dg1.addDataSet(h,i-1); 
            h = new H2F("ep_"+xyz+"_sf"+i+"_"+run, nb, -x, x, nb, -y,y); h.setTitle(tit2[i-1]); dg1.addDataSet(h,i+2);  
            h = new H2F("ep_"+xyz+"_sff"+i+"_"+run,nb, -x, x, nb, -y,y); h.setTitle(tit3[i-1]); dg1.addDataSet(h,i+5);  
        }
        int num = xyz=="xyz"?0:1;
        this.getDataGroup().add(dg1, num,0,k,run);
        this.getDataGroup().add(dg2, num,1,k,run);
        this.getDataGroup().add(dg3, num,2,k,run);
    }
    
    public void createUVWHistos(String tab, int cal, int xbins, int ybins, double xmin, double xmax, double ymin, double ymax, String xtxt, String ytxt) {
    	
        H2F h;  

        int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
        
        for (int is=1; is<7; is++) {      
        	switch (cal) {
        	case 0:
            DataGroup dg1 = new DataGroup(9,8); DataGroup dg2 = new DataGroup(8,8); DataGroup dg3 = new DataGroup(8,8);        	          
            for (int ip=1; ip<npmt[0]+1; ip++) {            
                h = new H2F("uvw-pcal-u"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"U"+ip);       
                dg1.addDataSet(h,ip-1); 
                h = new H2F("uvw-pcal-v"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"V"+ip);
                dg2.addDataSet(h,ip-1); 
                h = new H2F("uvw-pcal-w"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"W"+ip); 
                dg3.addDataSet(h,ip-1); 
     	    }
            this.getDataGroup().add(dg1,is,1,k,run); this.getDataGroup().add(dg2,is,2,k,run); this.getDataGroup().add(dg3,is,3,k,run);
            break;
            
        	case 1:
            DataGroup dg4 = new DataGroup(6,6); DataGroup dg5 = new DataGroup(6,6); DataGroup dg6 = new DataGroup(6,6);        	         	               
     	    for (int ip=1; ip<npmt[3]+1; ip++) {               
                h = new H2F("uvw-ecin-u"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECIN "+xtxt);  h.setTitleY(ytxt+"U"+ip); 
                dg4.addDataSet(h,ip-1);              
                h = new H2F("uvw-ecin-v"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"V"+ip); 
                dg5.addDataSet(h,ip-1); 
                h = new H2F("uvw-ecin-w"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"W"+ip);
                dg6.addDataSet(h,ip-1);                 
     	    }
            this.getDataGroup().add(dg4,is,4,k,run); this.getDataGroup().add(dg5,is,5,k,run); this.getDataGroup().add(dg6,is,6,k,run);
     	    break;
     	    
        	case 2:
            DataGroup dg7 = new DataGroup(6,6); DataGroup dg8 = new DataGroup(6,6); DataGroup dg9 = new DataGroup(6,6);        	         	              
     	    for (int ip=1; ip<npmt[6]+1; ip++) {               
                h = new H2F("uvw-ecou-u"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECOU "+xtxt); h.setTitleY(ytxt+"U"+ip);
                dg7.addDataSet(h,ip-1);                
                h = new H2F("uvw-ecou-v"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECOU "+xtxt); h.setTitleY(ytxt+"V"+ip);
                dg8.addDataSet(h,ip-1); 
                h = new H2F("uvw-ecou-w"+ip+"-s"+is+"-"+k+"-"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECOU "+xtxt); h.setTitleY(ytxt+"W"+ip);
                dg9.addDataSet(h,ip-1);   
     	    }
            this.getDataGroup().add(dg7,is,7,k,run); this.getDataGroup().add(dg8,is,8,k,run); this.getDataGroup().add(dg9,is,9,k,run);   
        	}
        }        
    }
    
    public void createSLCHistos(String tab, int pc, int cal, int nch, double y0, double y1, String txt, Boolean fun) {
    	
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
        H2F h;  F1D f0,f1,f2;
        f0 = new F1D("f0","[a]",1,68); f0.setParameter(0,0); 
        f1 = new F1D("f1","[a]",1,62); f1.setParameter(0,0); 
        f2 = new F1D("f2","[a]",1,36); f2.setParameter(0,0); 
       
        DataGroup dg = new DataGroup(6,3);
   
        for (int is=1; is<7; is++) {
            switch (cal) {
            case 0:
            h = new H2F("adc_pcal_u_"+is+"_"+pc+"_"+k+"_"+run, 68, 1., 69., nch, y0, y1);
            h.setTitleX("Sector "+is+" PCAL "+"U Strip"); h.setTitleY(txt); 
            dg.addDataSet(h,is-1);    if(fun) dg.addDataSet(f0, is-1);
            h = new H2F("adc_pcal_v_"+is+"_"+pc+"_"+k+"_"+run, 62, 1., 63., nch, y0, y1);
            h.setTitleX("Sector "+is+" PCAL "+"V Strip"); h.setTitleY(txt);        
            dg.addDataSet(h,is-1+6);  if(fun) dg.addDataSet(f1, is-1+6);           
            h = new H2F("adc_pcal_w_"+is+"_"+pc+"_"+k+"_"+run, 62, 1., 63., nch, y0, y1);
            h.setTitleX("Sector "+is+" PCAL "+"W Strip"); h.setTitleY(txt);  
            dg.addDataSet(h,is-1+12); if(fun) dg.addDataSet(f1, is-1+12);
            break;
            case 1:
            h = new H2F("adc_ecin_u_"+is+"_"+pc+"_"+k+"_"+run, 36, 1., 37., nch, y0, y1);
            h.setTitleX("Sector "+is+" ECIN "+"U Strip"); h.setTitleY(txt);    
            dg.addDataSet(h,is-1);    if(fun) dg.addDataSet(f2, is-1); 
            h = new H2F("adc_ecin_v_"+is+"_"+pc+"_"+k+"_"+run, 36, 1., 37., nch, y0, y1);
            h.setTitleX("Sector "+is+" ECIN "+"V Strip"); h.setTitleY(txt);        
            dg.addDataSet(h,is-1+6);  if(fun) dg.addDataSet(f2, is-1+6);           
            h = new H2F("adc_ecin_w_"+is+"_"+pc+"_"+k+"_"+run, 36, 1., 37., nch, y0, y1);
            h.setTitleX("Sector "+is+" ECIN "+"W Strip"); h.setTitleY(txt);  
            dg.addDataSet(h,is-1+12); if(fun) dg.addDataSet(f2, is-1+12); 
            break;
            case 2:
            h = new H2F("adc_ecou_u_"+is+"_"+pc+"_"+k+"_"+run, 36, 1., 37., nch, y0, y1);
            h.setTitleX("Sector "+is+" ECOU "+"U Strip"); h.setTitleY(txt);    
            dg.addDataSet(h,is-1);    if(fun) dg.addDataSet(f2, is-1);
            h = new H2F("adc_ecou_v_"+is+"_"+pc+"_"+k+"_"+run, 36, 1., 37., nch, y0, y1);
            h.setTitleX("Sector "+is+" ECOU "+"V Strip"); h.setTitleY(txt);        
            dg.addDataSet(h,is-1+6);  if(fun) dg.addDataSet(f2, is-1+6);           
            h = new H2F("adc_ecou_w_"+is+"_"+pc+"_"+k+"_"+run, 36, 1., 37., nch, y0, y1);
            h.setTitleX("Sector "+is+" ECOU "+"W Strip"); h.setTitleY(txt);  
            dg.addDataSet(h,is-1+12); if(fun) dg.addDataSet(f2, is-1+12);  
            }
        } 

        this.getDataGroup().add(dg,cal,pc,k,run);

    } 
    
    public void createMCHistos(String tab) {
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
        H1F h;  
        
        DataGroup dg = new DataGroup(6,3);
        
        for (int is=1; is<7; is++) {
            h = new H1F("adc_pcal_u_"+is+"_"+k+"_"+run, 100,-3,3); h.setOptStat("1100");
            h.setTitleX("Sector "+is+" PCAL Vertex Time (ns) ");  
            dg.addDataSet(h,is-1);        	
            h = new H1F("adc_ecin_u_"+is+"_"+k+"_"+run, 100,-3,3); h.setOptStat("1100");
            h.setTitleX("Sector "+is+" ECIN Vertex Time (ns) ");  
            dg.addDataSet(h,is-1+6);        	
            h = new H1F("adc_ecou_u_"+is+"_"+k+"_"+run, 100,-3,3); h.setOptStat("1100");
            h.setTitleX("Sector "+is+" ECOU Vertex Time (ns) ");  
            dg.addDataSet(h,is-1+12);        	
        }
        
        this.getDataGroup().add(dg,0,0,k,run);
    }
    
    public void initCCDB(int runno) {
        electron_sf = ebcm.getConstants(runno, "/calibration/eb/electron_sf");    
    }
    
    public Point3D squeeze(Point3D xyz, int det, int e_sect) {        
        xyz.rotateZ(Math.toRadians(-60*(e_sect-1)));
        xyz.rotateY(Math.toRadians(-25.));
        xyz.translateXYZ(-(det==1?40:50),0,0);
        xyz.rotateY(Math.toRadians(25.));
        xyz.rotateZ(Math.toRadians(60*(e_sect-1)));
        return xyz;
    }
    
    @Override
    public void processEvent(DataEvent event) {
    	
        isMC = (getRunNumber()<100) ? true:false;
        trigger_sect = isMC ? (event.hasBank("ECAL::adc") ? event.getBank("ECAL::adc").getByte("sector",0):5) : getElecTriggerSector(); 
        
        boolean goodSector = trigger_sect>0 && trigger_sect<7; 
        
//        if(!goodSector) return;
    	
    	int run = getRunNumber();
    	DataGroup dg = this.getDataGroup().getItem(0,0,0,run);
    	
    	float Ebeam=EB, e_mom=0, e_theta=0, e_vz=0, e_ecal_E=0, lU=0, lV=0, lW=0; 
    	float[]     x_ecal = {-1000,-1000,-1000};
        float[]     y_ecal = {-1000,-1000,-1000};
    	float[]    hx_ecal = {-1000,-1000,-1000};
        float[]    hy_ecal = {-1000,-1000,-1000};
    	float[]  e_ecal_TH = new float[3];
    	float[]  e_ecal_EL = new float[4];
    	float[]   e_ecal_u = new float[3];
    	float[]   e_ecal_v = new float[3];
    	float[]   e_ecal_w = new float[3];
    	float[]     t_ecal = new float[4];
    	float[]    pa_ecal = new float[4];
    	int e_sect=0;
    	int[] iU = new int[3], idU = new int[3]; float[] tU = new float[3];
    	int[] iV = new int[3], idV = new int[3]; float[] tV = new float[3];
    	int[] iW = new int[3], idW = new int[3]; float[] tW = new float[3];
    	
        if(dropBanks) dropBanks(event);
        
//        if (event.hasBank("MC::Particle")) return;

        boolean goodEvent = event.hasBank("REC::Particle")&& event.hasBank("REC::Calorimeter");
        
        if (!goodEvent) return;
        
        float Tvertex = event.hasBank("REC::Event") ? (isHipo3Event ? event.getBank("REC::Event").getFloat("STTime", 0):
                                                                      event.getBank("REC::Event").getFloat("startTime", 0)):0;        
      	DataBank    reccal = event.getBank("REC::Calorimeter");
      	DataBank ecalclust = event.getBank("ECAL::clusters");
      	DataBank ecalpeaks = event.getBank("ECAL::peaks");
      	DataBank ecalcalib = event.getBank("ECAL::calib");
        DataBank    recpar = event.getBank("REC::Particle");
        
      	Map<Integer,List<Integer>> caloMap = loadMapByIndex(reccal,"pindex");
      	Map<Integer,List<Integer>> partMap = loadMapByIndex(recpar,"pid");
      	
        int tpid = 11;

        if (!partMap.containsKey(tpid) || partMap.get(tpid).size()!=1) return;      	
        
        for (int ipart : partMap.get(tpid)) {			
            float px = recpar.getFloat("px", ipart);
            float py = recpar.getFloat("py", ipart);
            float pz = recpar.getFloat("pz", ipart);
            float vz = recpar.getFloat("vz", ipart);				
            float ep = (float)Math.sqrt(px*px+py*py+pz*pz);			
            float th = (float)Math.toDegrees(Math.acos(pz/ep));
            short status = (short) Math.abs(recpar.getShort("status", ipart));
            boolean inDC = (status>=2000 && status<3000);
            e_mom    = ep;
            e_theta  = th;
            if(inDC && ep>0.01*EB && ep<EB && th>4 && Math.abs(vz)<200 ){
               for (int icalo : caloMap.get(ipart)) {
				    int  det = reccal.getInt("layer", icalo);
	                short ic = reccal.getShort("index",icalo);
                    float  x = reccal.getFloat("x",icalo);
                    float  y = reccal.getFloat("y",icalo);
                    float  z = reccal.getFloat("z",icalo);					         
                    float hx = reccal.getFloat("hx",icalo);
                    float hy = reccal.getFloat("hy",icalo);
                    float hz = reccal.getFloat("hz",icalo);					         
                    float pa = reccal.getFloat("path", icalo);
                    float  t = reccal.getFloat("time",icalo);
                    float  r = (float) Math.sqrt(x*x+y*y+z*z);
                    if(det==1) e_sect = reccal.getByte("sector",icalo);	                 
                    int ind = getDet(det);
                    Point3D  xyz = squeeze(new Point3D( x, y, z),det,e_sect);
                    Point3D hxyz = squeeze(new Point3D(hx,hy,hz),det,e_sect);
                    if (det==1) {
                    	lU = reccal.getFloat("lu",icalo);
                    	lV = reccal.getFloat("lv",icalo);
                    	lW = reccal.getFloat("lw",icalo);
                    }
                    x_ecal[ind]     = (float) xyz.x();
                    y_ecal[ind]     = (float) xyz.y();
                    hx_ecal[ind]    = (float) hxyz.x();
                    hy_ecal[ind]    = (float) hxyz.y();
                    t_ecal[ind]     = t-Tvertex-pa/29.98f;
                    e_ecal_TH[ind]  = (float) Math.toDegrees(Math.acos(z/r));	               
                    e_ecal_EL[ind] += ecalclust.getFloat("energy", ic);
                    e_ecal_EL[3]   += ecalclust.getFloat("energy", ic);
                    iU[ind]         = (ecalclust.getInt("coordU", ic)-4)/8+1;
                    iV[ind]         = (ecalclust.getInt("coordV", ic)-4)/8+1;
                    iW[ind]         = (ecalclust.getInt("coordW", ic)-4)/8+1; 
                    int ic1         =  ecalclust.getInt("idU",ic)-1; //peak index U
                    int ic2         =  ecalclust.getInt("idV",ic)-1; //peak index V
                    int ic3         =  ecalclust.getInt("idW",ic)-1; //peak index W
                    e_ecal_u[ind]   =  ecalcalib.getFloat("recEU",ic); //peak energy U
                    e_ecal_v[ind]   =  ecalcalib.getFloat("recEV",ic); //peak energy V
                    e_ecal_w[ind]   =  ecalcalib.getFloat("recEW",ic); //peak energy W
                    /*
                    idU[ind]        = ecalclust.getInt("idU",ic);
                    idV[ind]        = ecalclust.getInt("idV",ic);
                    idW[ind]        = ecalclust.getInt("idW",ic);
                    tU[ind]         = ecalpeaks.getFloat("time",idU[ind]);
                    tV[ind]         = ecalpeaks.getFloat("time",idV[ind]);
                    tW[ind]         = ecalpeaks.getFloat("time",idW[ind]);
                    Point3D  point1 = new Point3D(ecalpeaks.getFloat("xo",idU[ind]),
		                                          ecalpeaks.getFloat("yo",idV[ind]),
		                                          ecalpeaks.getFloat("zo",idW[ind]));
                    Point3D  point2 = new Point3D(ecalpeaks.getFloat("xe",idU[ind]),
		                                          ecalpeaks.getFloat("ye",idV[ind]),
		                                          ecalpeaks.getFloat("ze",idW[ind]));
                    Point3D   point = new Point3D(x,y,z);
                    Line3D     line = new Line3D(point1,point2); 
                    */
                    
               }
           }
        }		
		
//        if(e_sect!=trigger_sect) return;	
      
        float[] sff = new float[4];
        
        for (int i=0; i<4; i++) sff[i] = e_ecal_EL[i]/e_mom;
        
        float sf = sff[3];
        
        boolean     good_e = e_sect>0 && e_sect<7 && e_mom>EB*0.02 && sff[3] > 0.02;         
        boolean good_fiduc = lV>15 && lW>15;
        
//      boolean good_fiduc = iU[0]>2&&iV[0]<63&&iW[0]<63&&iU[1]>2&&iV[1]<36&&iW[1]<36&&iU[2]>2&&iV[2]<36&&iW[2]<36;
       
        if(!good_e)  return;
        
        if(goodSector && good_fiduc) {
        	((H2F) getDG(0,0,"E/P",run).getData(e_sect-1  ).get(0)).fill(e_ecal_EL[3], sf);
        	((H2F) getDG(0,0,"E/P",run).getData(e_sect-1+6).get(0)).fill(e_mom,sf);
            ((H2F) getDG(0,0,"E/P",run).getData(e_sect-1+12).get(0)).fill(e_theta,sf);
            ((H2F) getDG(0,0,"E/P",run).getData(e_sect-1+18).get(0)).fill(e_ecal_TH[0],sf); 
            ((H1F) getDG(0,0,"PID Fits",run).getData(e_sect-1).get(0)).fill(sf/getSFcorr(e_sect,e_ecal_EL[3]));
        }
    	((H2F) getDG(1,0,"E/P",run).getData(e_sect-1   ).get(0)).fill(e_ecal_EL[3], sf);
    	((H2F) getDG(1,0,"E/P",run).getData(e_sect-1+ 6).get(0)).fill(e_mom,sf);
        ((H2F) getDG(1,0,"E/P",run).getData(e_sect-1+12).get(0)).fill(e_theta,sf);
        ((H2F) getDG(1,0,"E/P",run).getData(e_sect-1+18).get(0)).fill(e_ecal_TH[0],sf);  
        
        int is = e_sect;
        
        for (int id=0; id<3; id++) {
            ((H2F) getDG(0,0,"XY",run).getData(id).get(0)).fill(-x_ecal[id], y_ecal[id],sff[id]<0.5?1f:0);
            ((H2F) getDG(0,1,"XY",run).getData(id).get(0)).fill(-x_ecal[id], y_ecal[id],sff[id]<0.5?sff[id]:0.);
            ((H2F) getDG(0,2,"XY",run).getData(id).get(0)).fill(-x_ecal[id], y_ecal[id],sf<0.5?sf:0.);
            ((H2F) getDG(id,2,"SLC",run).getData(is-1   ).get(0)).fill(iU[id], e_ecal_u[id]/e_mom);
            ((H2F) getDG(id,2,"SLC",run).getData(is-1+ 6).get(0)).fill(iV[id], e_ecal_v[id]/e_mom);
            ((H2F) getDG(id,2,"SLC",run).getData(is-1+12).get(0)).fill(iW[id], e_ecal_w[id]/e_mom);
            ((H2F) getDG(id,1,"SLC",run).getData(is-1   ).get(0)).fill(iU[id], sff[id]<0.5?sff[id]:0.);
            ((H2F) getDG(id,1,"SLC",run).getData(is-1+ 6).get(0)).fill(iV[id], sff[id]<0.5?sff[id]:0.);
            ((H2F) getDG(id,1,"SLC",run).getData(is-1+12).get(0)).fill(iW[id], sff[id]<0.5?sff[id]:0.);				  
            ((H2F) getDG(id,0,"SLC",run).getData(is-1   ).get(0)).fill(iU[id], sf<0.5?sf:0.);				  
            ((H2F) getDG(id,0,"SLC",run).getData(is-1+ 6).get(0)).fill(iV[id], sf<0.5?sf:0.);				  
            ((H2F) getDG(id,0,"SLC",run).getData(is-1+12).get(0)).fill(iW[id], sf<0.5?sf:0.);	
            ((H2F) getDG(id,0,"Timing",run).getData(is-1   ).get(0)).fill(iU[id], t_ecal[id]);				  
            ((H2F) getDG(id,0,"Timing",run).getData(is-1+ 6).get(0)).fill(iV[id], t_ecal[id]);				  
            ((H2F) getDG(id,0,"Timing",run).getData(is-1+12).get(0)).fill(iW[id], t_ecal[id]);  
            
            ((H2F) getDG(is,3*id+1,"UVW",run).getData(iU[id]-1).get(0)).fill(e_mom, e_ecal_u[id]/e_mom);
            ((H2F) getDG(is,3*id+2,"UVW",run).getData(iV[id]-1).get(0)).fill(e_mom, e_ecal_v[id]/e_mom);
            ((H2F) getDG(is,3*id+3,"UVW",run).getData(iW[id]-1).get(0)).fill(e_mom, e_ecal_w[id]/e_mom);           
            
            if(sff[id]>0) {
            if(id==0) ((H1F) getDG(0,0,"MC",run).getData(is-1).get(0)).fill(-t_ecal[id]);
            if(id==1) ((H1F) getDG(0,0,"MC",run).getData(is-1+6).get(0)).fill(-t_ecal[id]);
            if(id==2) ((H1F) getDG(0,0,"MC",run).getData(is-1+12).get(0)).fill(-t_ecal[id]);
            }
         }       

    }
    
    public DataGroup getDG(int i, int j, String tab, int run) {
    	return this.getDataGroup().getItem(i,j,getDetectorTabNames().indexOf(tab),run);
    }
    
    public float getSFcorr(int is, float p) {
    	return 1+sfpar[is][0]/p + sfpar[is-1][1]/p/p;
    }
    
    public void updateFits(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
    	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));    	
        int    pc = 0;
        int    is = getActiveSector(); 
        int     i = getActiveLayer();
        int     j = getActiveView();
        
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
        if(dumpGraphs) dumpGraphs();   
    }

    public void analyze() {    
    	System.out.println(getDetectorName()+".Analyze() ");
    	fitGraphs(1,7,0,(dropSummary)?0:3,0,(dropSummary)?0:3); 
        writeFile("electron_sf");
        if(!isAnalyzeDone) createTimeLineHistos();
    	fillTimeLineHisto();
        System.out.println("Finished");
        isAnalyzeDone = true;    	
    } 
       
    public void fitGraphs(int is1, int is2, int id1, int id2, int il1, int il2) {    	
        int run = getRunNumber();	 
        
        // SF v UVW Fits
        for (int is=is1; is<is2; is++) {
 //           tl.fitData.add(fitEngine(((H2F) getDG(0,0,"E/P",run).getData(is-1).get(0)).projectionY(),0,0.15,0.3,0.15,0.3,1.7,2.5),is,0,7,run); 
            tl.fitData.add(fitEngine(((H1F) getDG(0,0,"PID Fits",run).getData(is-1).get(0)),0,0.15,0.3,0.15,0.3,1.7,2.5),is,0,7,run); 
           for (int id=id1; id<id2; id++) {
               for (int il=il1; il<il2; il++) {
               	  for (int pc=0; pc<1; pc++) {
                     H2F h = (H2F) getDG(id,pc,"SLC",run).getData(is-1+6*il).get(0);
           	         for (int i=0; i<npmt[id*3+il]; i++) tl.fitData.add(fitEngine(h.sliceX(i),0,0.15,0.3,0.15,0.3,1.7,2.5),is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); 
        	         fitStore(is, id, il, pc, run, 1f);
               	  }
               } 
           }
        }
        
    	for (int is=1; is<7; is++) { 
    		
            ParallelSliceFitter fitter;
              
            if (!dropSummary) {
            	
            // E/P vs. measured energy        	
            fitter = new ParallelSliceFitter((H2F) getDG(0,0,"E/P",run).getData(is-1).get(0));
            fitter.setBackgroundOrder(0); fitter.setMin(0.18); fitter.setMax(0.30); fitter.fitSlicesX(); 
            FitSummary.add(fitter.getMeanSlices(),is, 0, 7, run); 
          
    	    GraphErrors MeanGraph = fitter.getMeanSlices();
    	    tl.fitData.add(fitEngine(MeanGraph,14,0.3,EB*0.23,0.3,EB*0.23),is,0,5,run); 
    	    
    	    // E/P vs. tracking momentum
    		fitter = new ParallelSliceFitter((H2F) getDG(0,0,"E/P",run).getData(is-1+6).get(0));
    	    fitter.setBackgroundOrder(0); fitter.setMin(0.18); fitter.setMax(0.30); fitter.fitSlicesX();
    	    FitSummary.add(fitter.getMeanSlices(), is, 0, 1, run); 
    	    
    	    GraphErrors meanGraph = fitter.getMeanSlices(); 
    	    tl.fitData.add(fitEngine(meanGraph,14,2,7,2,7),is,0,6,run); 
    	    
    	    //SF sigE/E Fits
    	    
    	    GraphErrors  sigGraph = fitter.getSigmaSlices(); 
   	    
            int npts = meanGraph.getDataSize(0)  ;
            double[] xm  = new double[npts];
            double[] ym  = new double[npts];
            double[] yme = new double[npts];
            double xs,ys,yse;
            
            for (int i=0; i<meanGraph.getDataSize(0); i++) {
          	   xm[i] = 1/Math.sqrt(meanGraph.getDataX(i));
          	   ym[i] = meanGraph.getDataY(i);
          	  yme[i] = meanGraph.getDataEY(i);      			  
            } 
                	    
    	    GraphErrors  res0Graph = new GraphErrors();
            res0Graph.setTitleX("Sector "+is+" Electron Energy (GeV) ");  res0Graph.setTitleY("#sigma(E)/E"); 
            res0Graph.setMarkerSize(4); meanGraph.setMarkerStyle(1);            
    	        	
    	    int n = 0;
    	    for (int i=0; i<sigGraph.getDataSize(0); i++) {
                double y = sigGraph.getDataY(i); double ye = sigGraph.getDataEY(i);
                if(ym[i]>0&&y>0) {           
                    xs = sigGraph.getDataX(i);  //sig(E)/E vs True Energy
            	    ys = y/ym[i]; //sigma(E)/E = sigma(E/P)*(P/E)
            	   yse = ys*Math.sqrt(Math.pow(ye/y,2)+Math.pow(yme[i]/ym[i],2));
            	   res0Graph.addPoint(xs, ys, 0., yse);
            	}
            }  
    	    
    	    FitSummary.add(res0Graph, is, 0, 2, run);   
    	    
            GraphErrors resGraph = new GraphErrors();
            resGraph.setTitleX("Sector "+is+" 1/sqrt(Electron Energy (GeV)) ");  resGraph.setTitleY("#sigma(E)/E"); 
            resGraph.setMarkerSize(4); meanGraph.setMarkerStyle(1);            
            
            n=0;
            for (int i=0; i<sigGraph.getDataSize(0); i++) {
                double y = sigGraph.getDataY(i); double ye = sigGraph.getDataEY(i);
                if(ym[i]>0&&y>0) {
                    xs = 1/Math.sqrt(sigGraph.getDataX(i));  //sig(E)/E vs True Energy
            	    ys = y/ym[i]; //sigma(E)/E = sigma(E/P)*(P/E)
            	   yse = ys*Math.sqrt(Math.pow(ye/y,2)+Math.pow(yme[i]/ym[i],2));
            	   resGraph.addPoint(xs, ys, 0., yse);
            	}
            }  
            
            FitSummary.add(resGraph, is, 0, 3, run);  
            
            GraphErrors res2Graph = new GraphErrors();
            res2Graph.setTitleX("Sector "+is+"  (1/GeV) ");  res2Graph.setTitleY("[#sigma(E)/E]^2");             
            res2Graph.setMarkerSize(4); meanGraph.setMarkerStyle(1);
            
            n=0;
            for (int i=0; i<sigGraph.getDataSize(0); i++) {
                double y = sigGraph.getDataY(i); double ye = sigGraph.getDataEY(i);
                if(ym[i]>0&&y>0) {
                    xs = 1/sigGraph.getDataX(i);  //sig(E)/E vs True Energy
            	    ys = Math.pow(y/ym[i],2); //sigma(E)/E = sigma(E/P)*(P/E)
            	   yse = 2*ys*Math.sqrt(Math.pow(ye/y,2)+Math.pow(yme[i]/ym[i],2));
            	   if(ys>0&&yse/ys<0.12) res2Graph.addPoint(xs, ys, 0., yse);
            	}
            }  
            System.out.println("Fit 4: "+is+" "+sigGraph.getDataSize(0));
            FitSummary.add(res2Graph, is, 0, 4, run);  
            tl.fitData.add(fitEngine(res2Graph,6,0.,0.6,0.,0.6),is,0,4,run); 
            }
    	}    	
    }
   
    public void fitStore(int is, int id, int il, int pc, int run, float nrm) {
    	int np = npmt[id*3+il];
        double[]      x = new double[np]; double[]  ymean = new double[np]; double[] yrms = new double[np];
        double[]     xe = new double[np]; double[] ymeane = new double[np]; double[]   ye = new double[np]; 
        double[]  yMean = new double[np]; 
        for (int i=0; i<np; i++) {
            x[i] = i+1; xe[i]=0; ye[i]=0; yrms[i]=0; 
            FitData fd = tl.fitData.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); 
            fd.graph.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            fd.hist.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            double mean = fd.mean;                        
            if(mean>0) yrms[i] = fd.sigma/mean; 
                      yMean[i] = fd.getMean()/nrm;
                      ymean[i] = mean/nrm;
                     ymeane[i] = fd.meane/nrm;
        }
        GraphErrors mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,ymean,xe,ymeane);                   
        GraphErrors Mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yMean,xe,ymeane);                   
        GraphErrors  rms = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yrms,xe,ye);                  
        FitSummary.add(mean, is,id+10*(pc+1)*(pc+1)*(il+1),1,run);
        FitSummary.add(rms,  is,id+10*(pc+1)*(pc+1)*(il+1),2,run);                    
        FitSummary.add(Mean, is,id+10*(pc+1)*(pc+1)*(il+1),5,run);        	        
    }    
    
    public void plotMeanHWSummary(String tab) {
        
    	int index=getDetectorTabNames().indexOf(tab);
    	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int           pc = 0;
        int            n = 0;
        
        List<DataLine> lines = new ArrayList<DataLine>();
        
        Boolean t = TLname!="UVW";
        double ymin=0.2f, ymax=0.3f;
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
        	        	DataLine line = new DataLine(m,ymin,m,ymax) ; line.setLineColor(1); line.setLineWidth(1); 
        	        	lines.add(line);
        	        }
                }
            }
            hwplot2.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin,ymax); c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==3||n==6||n==9||n==12||n==15) hwplot1.getAttributes().setTitleY("SF");
            hwplot1.getAttributes().setTitleX("SECTOR "+is+" "+det[id]+" PMT");
            F1D f1 = new F1D("p0","[a]",0.,m); f1.setParameter(0,0.25); f1.setLineColor(3); f1.setLineWidth(2);
            c.draw(hwplot1); /*c.draw(hwplot2,"same");*/ for(DataLine line: lines) c.draw(line); c.draw(f1,"same"); n++; 
        }
        }        
    }  
        
    public void plotPIDFits(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(6, 3);
        GraphErrors g1 = new GraphErrors(), g2 = new GraphErrors();
        String txt;
        
        int run = getRunNumber(); EB = getBeamEnergy(run);
  	   
    	for (int is=1; is<7; is++) {  
    	    SFFunction sf = new SFFunction("esf",11,is,ebccdb,0.1,2.5);  sf.setLineWidth(2) ; sf.setLineColor(1);
   		    txt = "Sector "+is+" Measured Energy (GeV)";
            if (FitSummary.hasItem(is,0,7,run)) {
            	c.cd(is-1);c.getPad().getAxisZ().setLog(true);  c.draw((H2F) getDG(0,0,"E/P",run).getData(is-1).get(0));
            	GraphPlot((GraphErrors)FitSummary.getItem(is,0,7,run),c,is-1,0.0f,EB*0.25f,0.15f,0.35f,1,6,1,txt," E/P","same"); //c.draw(sf,"same");
            	tl.fitData.getItem(is,0,5,run).graph.getFunction().setLineColor(20); 
            	tl.fitData.getItem(is,0,5,run).graph.getFunction().setLineWidth(6);
            	tl.fitData.getItem(is,0,5,run).graph.getFunction().setOptStat("1110");
            	tl.fitData.getItem(is,0,5,run).graph.getFunction().setRange(0.1,EB*0.25);
            	c.draw(tl.fitData.getItem(is,0,5,run).graph.getFunction(),"same");             	
            }
         
            c.cd(is-1+6);c.getPad(is-1+12).getAxisX().setRange(0.1, 0.4); 
            c.draw(Fits.getItem(is,0,7,getRunNumber()).getHist());
            c.draw(Fits.getItem(is,0,7,getRunNumber()).getGraph(),"same");   
                        
            GraphPlot((GraphErrors)tl.fitData.getItem(is,0,4,run).getGraph(),c,is-1+12,0.f,0.6f,0.001f,0.008f,1,4,1,"","",""); 
            g1.addPoint(is,Math.sqrt(tl.fitData.getItem(is,0,4,run).p1),0,Math.sqrt(tl.fitData.getItem(is,0,4,run).p1e));
            g2.addPoint(is,Math.sqrt(tl.fitData.getItem(is,0,4,run).p0),0,Math.sqrt(tl.fitData.getItem(is,0,4,run).p0e)); 
     	}
    }    
    
    public void plotFitSummary1(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
    	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(3, 2);
        int col[] = {1,2,3,4,5,7};
        int run = getRunNumber();
    	for (int is=1; is<7; is++) {  
    		String txt = "Sector "+is+" Electron Energy (GeV)";
            if (FitSummary.hasItem(1,is,run)) GraphPlot((GraphErrors)FitSummary.getItem(is,0,1,run),c,is-1,1.5f,10.f,0.18f,0.300f,col[is-1],4,1,txt," E/P","");
    	}
    }    
    
    public void plotFitSummary2(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(3, 2);        
        int col[] = {1,2,3,4,5,7};
        int run = getRunNumber();
    	for (int is=1; is<7; is++) {    		
            F1D f = new F1D("res","sqrt([a]*[a]/x+[b]*[b])",1.6,10.0); f.setLineColor(1); f.setLineWidth(3);
            f.setParameter(0, par[is-1][0]);f.setParameter(1, par[is-1][1]);
            if (FitSummary.hasItem(is,0,2,run)) GraphPlot((GraphErrors)FitSummary.getItem(is,0,2,run),c,is-1,1.5f,10.f,0.034f,0.088f,col[is-1],4,1,"","","");
            if (FitSummary.hasItem(is,0,2,run)) c.draw(f,"same");            
    	}    	
    }
    
    public void plotFitSummary3(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(3, 2);
        int col[] = {1,2,3,4,5,7};
        int run = getRunNumber();
    	for (int is=1; is<7; is++) {    		
            F1D f = new F1D("res","sqrt([a]*[a]*x*x+[b]*[b])",0.3,0.76); f.setLineColor(1); f.setLineWidth(3);
            f.setParameter(0, par[is-1][0]);f.setParameter(1, par[is-1][1]);
            if (FitSummary.hasItem(is,0,3,run)) GraphPlot((GraphErrors)FitSummary.getItem(is,0,3,run),c,is-1,0.3f,0.76f,0.04f,0.09f,col[is-1],4,1,"","","");
            if (FitSummary.hasItem(is,0,3,run)) c.draw(f,"same");            
    	}    	
    } 
    
    public void plotFitSummary4(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
    	GraphErrors g1 = new GraphErrors(), g2 = new GraphErrors();
        c.setGridX(false); c.setGridY(false); c.divide(6, 3);
        int col[] = {1,2,3,4,5,7};
        int run = getRunNumber();
     	for (int is=1; is<7; is++) {    		
//            F1D f = new F1D("res","[a]*[a]*x+[b]*[b]",0.,0.6); f.setLineColor(1); f.setLineWidth(3);
//            f.setParameter(0, par[is-1][0]);f.setParameter(1, par[is-1][1]);
//            if (FitSummary.hasItem(is,0,4,run)) GraphPlot((GraphErrors)FitSummary.getItem(is,0,4,run),c,is-1,0.f,0.6f,0.001f,0.008f,col[is-1],4,1,"","","");
              GraphPlot((GraphErrors)tl.fitData.getItem(is,0,4,run).getGraph(),c,is-1+12,0.f,0.6f,0.001f,0.008f,col[is-1],4,1,"","",""); 
              g1.addPoint(is,Math.sqrt(tl.fitData.getItem(is,0,4,run).p1),0,Math.sqrt(tl.fitData.getItem(is,0,4,run).p1e));
              g2.addPoint(is,Math.sqrt(tl.fitData.getItem(is,0,4,run).p0),0,Math.sqrt(tl.fitData.getItem(is,0,4,run).p0e));
//            if (FitSummary.hasItem(is,0,4,run)) c.draw(f,"same");            
    	}    	
//        GraphPlot(g1,c,7,0.5f,6.5f,0.0f,0.11f,1,6,1,"SECTOR","",""); GraphPlot(g2,c,7,0.5f,6.5f,0.0f,0.11f,1,7,2,"","","same"); 
    }  
    
    public void plotXYZHistos(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
 	    EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
 	    DataGroup dg = getDataGroup().getItem(0,0,index,getRunNumber());
 	    
 	    c.setGridX(false); c.setGridY(false);
 	    c.divide(3, 3);
        
	    H2F h2 ;
	    double[] zmin={0.1,0.025,0.002,0.22,0.22,0.22}, zmax={0.245,0.110,0.029,0.275,0.275,0.275};
	    
	    if(zMax<50) {
	    	zmin[0]=0.001*zMin; zmin[1]=0.001*zMin; zmin[2]=0.001*zMin; zmin[3]=0.001*zMin; zmin[4]=0.001*zMin; zmin[5]=0.001*zMin;
	    	zmax[0]=0.006*zMax; zmax[1]=0.004*zMax; zmax[2]=0.001*zMax; zmax[3]=0.010*zMax; zmax[4]=0.010*zMax; zmax[5]=0.010*zMax;	    	
	    }	    
	    
	    c.cd(0); c.getPad().getAxisZ().setLog(getLogZ()); h2 = (H2F) dg.getData(0).get(0); c.draw(h2);
	    c.cd(1); c.getPad().getAxisZ().setLog(getLogZ()); h2 = (H2F) dg.getData(1).get(0); c.draw(h2);
	    c.cd(2); c.getPad().getAxisZ().setLog(getLogZ()); h2 = (H2F) dg.getData(2).get(0); c.draw(h2);
	   
	    for (int i=0; i<6; i++) {c.cd(i+3); c.getPad().getAxisZ().setRange(zmin[i],zmax[i]); h2 = (H2F) dg.getData(i+3).get(0); c.draw(h2);}
   	
    }    
    
    public void plotSLCHistos(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
    	if(!getDataGroup().hasItem(getActiveSector(),getActivePC(),index,getRunNumber())) return;
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveLayer(),getActivePC(),index,getRunNumber()));	       	
    }
    
    public void plotUVWHistos(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
    	if(!getDataGroup().hasItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index,getRunNumber())) return;
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index,getRunNumber()));	       	
    } 
    
    public void plotSECHistos(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActivePC(),0,index,getRunNumber()));	       	
    } 
    
    public void plotMCHistos(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));	       	    	
    }
    
	public void writeFile(String table) {
		
		String line = new String();
		
		try { 
			File outputFile = new File(outPath+"ECsf/electron_sf/"+table+"_"+getRunNumber());
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("ECsf.writefile("+table+")");
			
			for (int is=0; is<7; is++) {
				switch (table) {							
				case "electron_sf": line =  getSF(is);  break;
				}
				if (line!=null) {
				      System.out.println(line);
				      outputBw.write(line);
				      outputBw.newLine();
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
	
	public String getSF(int is) {
		if(is==0) {
			return is+" 0  0 "
					+String.format("%.5f",electron_sf.getDoubleValue("sf1",is,0,0))+" "
					+String.format("%.5f",electron_sf.getDoubleValue("sf2",is,0,0))+" "
					+String.format("%.5f",electron_sf.getDoubleValue("sf3",is,0,0))+" "
					+String.format("%.5f",electron_sf.getDoubleValue("sf4",is,0,0))+" "
					+String.format("%.5f",electron_sf.getDoubleValue("sfs1",is,0,0))+" "
					+String.format("%.5f",electron_sf.getDoubleValue("sfs2",is,0,0))+" "
					+String.format("%.5f",electron_sf.getDoubleValue("sfs3",is,0,0))+" "
					+String.format("%.5f",electron_sf.getDoubleValue("sfs4",is,0,0));
		} else {
			return is+" 0  0 "
					+String.format("%.5f",tl.fitData.getItem(is,0,5,getRunNumber()).p0)+" "
					+" 1.0 "
					+String.format("%.5f",tl.fitData.getItem(is,0,5,getRunNumber()).p1)+" "
					+String.format("%.5f",tl.fitData.getItem(is,0,5,getRunNumber()).p2)+" "
					+String.format("%.5f",tl.fitData.getItem(is,0,7,getRunNumber()).sigma)+" "
					+"1.0 0.0 0.0";					
		}
	}
	
	public void readSF(String fname) {
		
        try{
            File f = new File(fname); 
            if (!f.exists()) return;
            BufferedReader reader = new BufferedReader(new FileReader(f));
			System.out.println("ECsf.readSF("+fname+")");
            int n = 0 ;
            while (n<7) {		
                String line = reader.readLine();
                String[] col = line.trim().split("\\s+"); 
                int is = Integer.parseInt(col[0]); 
                sfpar[is][0] = Float.parseFloat(col[5]);
                sfpar[is][1] = Float.parseFloat(col[6]);
                n++;
            }
            reader.close();   
        }
        catch(FileNotFoundException ex) {
            ex.printStackTrace();            
        }     
        catch(IOException ex) {
            ex.printStackTrace();
        }        
	}
        
    public void dumpGraphs() {
    	int run = getRunNumber();
    	for (int is=1; is<7; is++) {
    		dumpGraph(vecPath+"meanGraph10_"+is+"_"+run,(GraphErrors)FitSummary.getItem(7,is,run));
    		dumpGraph(vecPath+"meanGraph_"+is+"_"+run,  (GraphErrors)FitSummary.getItem(1,is,run));
    		dumpGraph(vecPath+"res0Graph10_"+is+"_"+run,(GraphErrors)FitSummary.getItem(2,is,run));
    		dumpGraph(vecPath+"resGraph_"+is+"_"+run,   (GraphErrors)FitSummary.getItem(3,is,run));
    		dumpGraph(vecPath+"res2Graph_"+is+"_"+run,  (GraphErrors)FitSummary.getItem(4,is,run));
    	}
    }
    
/*   TIMELINES */
    
    public void createTimeLineHistos() {   
    	System.out.println("Initializing "+TLname+" timeline"); 
    	runIndex = 0;
    	tl.createTimeLineHisto(10,"ECAL E/P","Sector",TLmax,6,1,7);
    	tl.createTimeLineHisto(20,"ECAL #sigma(E)/E","Sector",TLmax,6,1,7);
    }
    
    public void fillTimeLineHisto() {
        for (int is=1; is<7; is++) {
            float   y = (float) tl.fitData.getItem(is,0,7,getRunNumber()).mean; 
            float  ye = (float) tl.fitData.getItem(is,0,7,getRunNumber()).meane;			 
            float  ys = (float) tl.fitData.getItem(is,0,7,getRunNumber()).sigma;
            float yse = (float) tl.fitData.getItem(is,0,7,getRunNumber()).sigmae;			 
            ((H2F)tl.Timeline.getItem(10,0)).fill(runIndex,is,y);	
            ((H2F)tl.Timeline.getItem(10,1)).fill(runIndex,is,ye);   		
            ((H2F)tl.Timeline.getItem(20,0)).fill(runIndex,is,ys/y);	
            ((H2F)tl.Timeline.getItem(20,1)).fill(runIndex,is,(ys/y)*Math.sqrt(Math.pow(yse/ys,2)+Math.pow(ye/y,2)));   
            
        } 
        runIndex++;
    }
    
    public void saveTimelines() {
    	System.out.println("ECsf: Saving timelines");
    	saveTimeLine(10,0,7,"SampFrac","SF");
    	saveTimeLine(20,0,7,"Resolution","SF");
    }
    
    public void plotTimeLines(String tab) {
    	int index=getDetectorTabNames().indexOf(tab);
    	if(TLflag) {plotTimeLineSectors(index); } else {plotClusterTimeLines(index);}
    }   
    
    public void plotClusterTimeLines(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); 
        int           is = getActiveSector(); 
        FitData       fd = tl.fitData.getItem(is,0,7,getRunNumber());
        
    	double mean=0.25f;
       
    	DataLine line1 = new DataLine(0,is,  runIndex+1,is);                   line1.setLineColor(5);
    	DataLine line2 = new DataLine(0,is+1,runIndex+1,is+1);                 line2.setLineColor(5);
    	DataLine line3 = new DataLine(runIndexSlider,1,  runIndexSlider,7);    line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,runIndexSlider+1,7);  line4.setLineColor(5);
    	
        c.clear(); c.divide(3, 2); 

        for (int i=0; i<2; i++) { int i3=i*3;
            double min=(i==0)?0.24:0.04; double max=(i==0)?0.26:0.1; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(i3); c.getPad(i3).setAxisRange(0,runIndex,1,7); c.getPad(i3).setTitleFontSize(18); c.getPad(i3).getAxisZ().setRange(min,max);
    		c.draw((H2F)tl.Timeline.getItem(10*(i+1),0));c.draw(line1);c.draw(line2);c.draw(line3);c.draw(line4);
             
    		c.cd(i3+1); c.getPad(i3+1).setAxisRange(-0.5,runIndex,min,max); c.getPad(i3+1).setTitleFontSize(18);
    		drawTimeLine(c,is,10*(i+1),((i==0)?0.25f:0.06f),"Sector "+is+((i==0)?"  E / P":" #sigma(E) / E"));
    		    
    		c.cd(i3+2); c.getPad(i3+2).setAxisRange(0.1,0.4,0.,fd.getGraph().getMax()*1.1);
    		fd.getHist().getAttributes().setOptStat("1000100");
    		DataLine line6 = new DataLine(mean,-50,mean,fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);
    		c.draw(fd.getHist()); c.draw(fd.getGraph(),"same"); c.draw(line6);
        }
    }
     
    public void plotPeakTimeLines(int index) {
    	
    }
    
    public void plotTimeLineSectors(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int pc = getActivePC();
        c.clear(); c.divide(3, 2);
    	for (int is=1; is<7; is++) {
    		double min=(pc==0)?0.24:0.04; double max=(pc==0)?0.26:0.1; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(is-1); c.getPad(is-1).setAxisRange(-0.5,runIndex,min,max); c.getPad(is-1).setTitleFontSize(18);
    		drawTimeLine(c,is,10*(pc+1),((pc==0)?0.25f:0.06f),"Sector "+is+((pc==0)?"  E / P":"  #sigma(E) / E"));
    	}
    }

    @Override
    public void timerUpdate() {  
    	
        for(int i=0; i<3; i++) {
        H2F e  = (H2F) getDG(0,0,"XY",getRunNumber()).getData(i).get(0);
        H2F w  = (H2F) getDG(0,1,"XY",getRunNumber()).getData(i).get(0);
        H2F ww = (H2F) getDG(0,2,"XY",getRunNumber()).getData(i).get(0);
        for(int loop = 0; loop < e.getDataBufferSize(); loop++) {
        	    float ne = e.getDataBufferBin(loop);
            if (ne>0) {H2F h = (H2F) getDG(0,0,"XY",getRunNumber()).getData(i+3).get(0); h.setDataBufferBin(loop,w.getDataBufferBin(loop)/ne);}
            if (ne>0) {H2F h = (H2F) getDG(0,0,"XY",getRunNumber()).getData(i+6).get(0); h.setDataBufferBin(loop,ww.getDataBufferBin(loop)/ne);}
        }
        }
        
    }
    
}
