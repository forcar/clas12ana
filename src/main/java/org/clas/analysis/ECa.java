package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.analysis.ECPart.SFFunction;
import org.clas.tools.FitData;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.math.StatNumber;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.eb.EBCCDBConstants;
import org.jlab.rec.eb.SamplingFractions;
import org.jlab.utils.groups.IndexedList;


public class ECa extends DetectorMonitor {

    private final int[] npaddles = new int[]{68,62,62,36,36,36,36,36,36};
    
	int Nevts, Nelecs, Ntrigs, runNum;
	boolean[] trigger_bits;
	public float EB, Ebeam;
	public float RFtime1, RFtime2, startTime;
	public long  TriggerWord;
	public int   trig_part_ind, trig_sect, trig_track_ind;
	public int   e_part_ind, e_sect, e_track_ind, hasLTCC, ngammas;
	public int   pip_part_ind, pip_track_ind, pip_sect, pim_part_ind, pim_track_ind, pim_sect;
	public int   foundCVT, CVTcharge;
	public float e_mom, e_theta, e_phi, e_vx, e_vy, e_vz, e_Ivy, e_Ivz, e_ecal_X, e_ecal_Y, e_ecal_Z, e_ecal_E;
	public float pim_ecal_E;
	public float pip_ecal_E;
	public float e_ecal_TH[] = new float[3];
	public float e_ecal_EL[] = new float[3];
	public float pim_ecal_TH[] = new float[3];
	public float pim_ecal_EL[] = new float[3];
	public float pip_ecal_TH[] = new float[3];
	public float pip_ecal_EL[] = new float[3];
	public float iU[] = new float[3];
	public float iV[] = new float[3];
	public float iW[] = new float[3];	
	public float pim_iU[] = new float[3];
	public float pim_iV[] = new float[3];
	public float pim_iW[] = new float[3];	
	public float pip_iU[] = new float[3];
	public float pip_iV[] = new float[3];
	public float pip_iW[] = new float[3];	
	public float x_ecal[] = new float[3];
	public float y_ecal[] = new float[3];
	public float pim_x_ecal[] = new float[3];
	public float pim_y_ecal[] = new float[3];
	public float pip_x_ecal[] = new float[3];
	public float pip_y_ecal[] = new float[3];
	public float e_track_chi2, e_vert_time, e_vert_time_RF, e_Q2, e_xB, e_W;
	public float e_HTCC, e_LTCC, e_pcal_e, e_etot_e, e_TOF_X, e_TOF_Y, e_TOF_Z, e_HTCC_X, e_HTCC_Y, e_HTCC_Z;
	public float e_DCR1_X, e_DCR1_Y, e_DCR1_Z, e_DCR2_X, e_DCR2_Y, e_DCR2_Z, e_DCR3_X, e_DCR3_Y, e_DCR3_Z;
	public float g1_e, g1_theta, g1_phi, g2_e, g2_theta, g2_phi;
	public float pip_mom, pip_theta, pip_phi, pip_vx, pip_vy, pip_vz, pip_vert_time, pip_beta, pip_track_chi2;
	public float pim_mom, pim_theta, pim_phi, pim_vx, pim_vy, pim_vz, pim_vert_time, pim_beta, pim_track_chi2;
	public float CVT_mom, CVT_theta, CVT_phi, CVT_vz, CVT_chi2, CVT_pathlength;
	public int CVT_ndf;
	public LorentzVector VB, VT, Ve, VG1, VG2, VPI0, VPIP, VPIM;  
    
    public double[][] par = {{0.105,0.039},{0.099,0.040},{0.100,0.034},{0.093,0.044},{0.085,0.046},{0.113,0.028}};
    
    public ECa(String name) {
        super(name);
        this.setDetectorTabNames("E/P v P",
                                 "E/P v ThV",
                                 "E/P v ThD",
                                 "EPC/P v ThD",
                                 "EECi/P v ThD",
                                 "EECo/P v ThD",
                                 "E/P v XY",
                                 "E/P v UVW",
                                 "PIM v UVW",
                                 "PIP v UVW",
                                 "E/P v Em",
                                 "Fits E/P",
                                 "Fits Sig 1",
                                 "Fits Sig 2",
                                 "Fits Sig 3",
                                 "Timeline");
        this.usePCCheckBox(true);
        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
    	System.out.println("ECa.localinit()");
        EB=12; Ebeam = 10.6f;
      	VB = new LorentzVector(0,0,Ebeam,Ebeam);
	    VT = new LorentzVector(0,0,0,0.93827);
	    tl.setFitData(Fits);
	    configEngine("muon");  
        configEventBuilder();
    }
    
    public void localclear() {
    	System.out.println("ECa.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	slider.setValue(0);
    }
    
    @Override
    public void createHistos(int run) {        
	    System.out.println("ECa:createHistos("+run+")");
        setRunNumber(run);
        runlist.add(run);
        this.setNumberOfEvents(0);        
        createEOPHistos(10,50,0.,2.5,"ep_em", " Measured Energy (GeV)", " E/P");
        if(dropSummary) return;
        createEOPHistos(0,50,0.5,10.5,"ep_p",  " Momentum (GeV)",      " E/P");
        createEOPHistos(1,30,  6.,36.,"ep_thv"," VertexTheta (deg)",   " E/P");
        createEOPHistos(2,48,  3.,37.,"ep_thd"," Detector Theta (deg)"," E/P");
        createEOPHistos(3,48,  3.,27.,"ep_th0"," PC Theta (deg)",      "EPC / P");
        createEOPHistos(4,48,  3.,27.,"ep_th1"," ECIN Theta (deg)",    "EECi / P");
        createEOPHistos(5,48,  3.,27.,"ep_th2"," ECOU Theta (deg)",    "EECo / P");
        createXYZHistos(6);
        createADCHistos(7,25,0.5,0.3,0.1,"SF");
        createADCHistos(8,50,100.,100.,100.,"PIM (MeV)");
        createADCHistos(9,50,100.,100.,100.,"PIP (MeV)");       
    }

    @Override       
    public void plotHistos(int run) {
    	plotSummary(run);
    	plotAnalysis(run);
    }
      
    public void plotSummary(int run) {
    	setRunNumber(run);
    	plotEOPHistos(10);
   	    if(dropSummary) return;
    	plotEOPHistos(0);
    	plotEOPHistos(1);
    	plotEOPHistos(2);
    	plotEOPHistos(3);
    	plotEOPHistos(4);
    	plotEOPHistos(5);
    	plotXYZHistos(6);
    	plotADCHistos(7);
    	plotADCHistos(8);
    	plotADCHistos(9);
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;
    	plotFitSummary10(11);
    	if(!dropSummary) {plotFitSummary2(12); plotFitSummary3(13); plotFitSummary4(14);}
    	plotTimeLines(15);
    }
    
    public void plotADCHistos(int index) {
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),0,index,getRunNumber()));	       	
    }
    
    public void plotEOPHistos(int index) {
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));	       	
    }  
    
    public void plotXYZHistos(int index) {
 	    EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
 	    DataGroup dg = getDataGroup().getItem(0,0,index,getRunNumber());
 	    
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 3);
        
	    H2F h2 ;
	    
	    c.cd(0); c.getPad().getAxisZ().setLog(getLogZ());                h2 = (H2F) dg.getData(0).get(0); c.draw(h2);
	    c.cd(1); c.getPad().getAxisZ().setLog(getLogZ());                h2 = (H2F) dg.getData(1).get(0); c.draw(h2);
	    c.cd(2); c.getPad().getAxisZ().setLog(getLogZ());                h2 = (H2F) dg.getData(2).get(0); c.draw(h2);
	   
	    c.cd(3); c.getPad().getAxisZ().setRange(0.001*zMin, 0.006*zMax); h2 = (H2F) dg.getData(3).get(0); c.draw(h2);
	    c.cd(4); c.getPad().getAxisZ().setRange(0.001*zMin, 0.004*zMax); h2 = (H2F) dg.getData(4).get(0); c.draw(h2);
	    c.cd(5); c.getPad().getAxisZ().setRange(0.001*zMin, 0.001*zMax); h2 = (H2F) dg.getData(5).get(0); c.draw(h2);
	    
	    c.cd(6); c.getPad().getAxisZ().setRange(0.001*zMin, 0.010*zMax); h2 = (H2F) dg.getData(6).get(0); c.draw(h2);
	    c.cd(7); c.getPad().getAxisZ().setRange(0.001*zMin, 0.010*zMax); h2 = (H2F) dg.getData(7).get(0); c.draw(h2);
	    c.cd(8); c.getPad().getAxisZ().setRange(0.001*zMin, 0.010*zMax); h2 = (H2F) dg.getData(8).get(0); c.draw(h2);    	
    }
    
    public void plotFitSummary1(int index) {
 	    EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(3, 2);
        int col[] = {1,2,3,4,5,7};
	    int run = getRunNumber();
    	for (int is=1; is<7; is++) {  
    		String txt = "Sector "+is+" Electron Energy (GeV)";
            if (FitSummary.hasItem(1,is,run)) GraphPlot((GraphErrors)FitSummary.getItem(is,0,1,run),c,is-1,1.5f,10.f,0.18f,0.300f,col[is-1],4,1,txt," E/P","");
    	}
    }    
    
    public void plotFitSummary10(int index) {
 	    EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(3, 4);
        int col[] = {1,2,3,4,5,7};
	    int run = getRunNumber();    
	   
	    SFFunction sf = new SFFunction("esf",-11,eb.ccdb,0.1,2.5); 

    	for (int is=1; is<7; is++) {  
    		String txt = "Sector "+is+" Measured Energy (GeV)";
            if (FitSummary.hasItem(is,0,10,run)) {
            	GraphPlot((GraphErrors)FitSummary.getItem(is,0,10,run),c,is-1,0.0f,2.5f,0.18f,0.300f,col[is-1],4,1,txt," E/P",""); c.draw(sf,"same");
            }
            c.cd(is-1+6);c.getPad(is-1+6).getAxisX().setRange(0.1, 0.4); 
            c.draw(Fits.getItem(is,0,10,getRunNumber()).getHist());
            c.draw(Fits.getItem(is,0,10,getRunNumber()).getGraph(),"same");     	
     	}
    }
    
    public void plotFitSummary2(int index) {
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
    
    public void plotFitSummary3(int index) {
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
    
    public void plotFitSummary4(int index) {
 	    EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(3, 2);
        int col[] = {1,2,3,4,5,7};
	    int run = getRunNumber();
    	for (int is=1; is<7; is++) {    		
            F1D f = new F1D("res","[a]*[a]*x+[b]*[b]",0.,0.6); f.setLineColor(1); f.setLineWidth(3);
            f.setParameter(0, par[is-1][0]);f.setParameter(1, par[is-1][1]);
            if (FitSummary.hasItem(is,0,4,run)) GraphPlot((GraphErrors)FitSummary.getItem(is,0,4,run),c,is-1,0.f,0.6f,0.001f,0.008f,col[is-1],4,1,"","","");
            if (FitSummary.hasItem(is,0,4,run)) c.draw(f,"same");            
    	}    	
    }  
    
    public void GraphPlot(GraphErrors graph, EmbeddedCanvas c, int zone, float xmin, float xmax, float ymin, float ymax, int mcol, int msiz, int msty, String xtit, String ytit, String opt) {
    	c.cd(zone); c.getPad(zone).getAxisX().setRange(xmin, xmax); c.getPad(zone).getAxisY().setRange(ymin, ymax); 
    	graph.setMarkerColor(mcol); graph.setMarkerSize(msiz); graph.setMarkerStyle(msty); graph.setLineColor(mcol);
    	if(!(xtit=="")) graph.setTitleX(xtit); 
    	if(!(ytit=="")) graph.setTitleY(ytit); 
    	c.draw(graph,opt);
    }
    
    @Override
    public void plotEvent(DataEvent de) {
    	    analyze();
    }

    public void analyze() {    
    	System.out.println(getDetectorName()+".Analyze() ");
    	fitGraphs(); if(dumpGraphs) dumpGraphs();
        if(!isAnalyzeDone) createTimeLineHistos();
    	fillTimeLineHisto();
        System.out.println("Finished");
        isAnalyzeDone = true;
    } 
  
    public void dumpGraphs() {
    	int run = getRunNumber();
    	for (int is=1; is<7; is++) {
    		dumpGraph(vecPath+"meanGraph10_"+is+"_"+run,(GraphErrors)FitSummary.getItem(10,is,run));
    		dumpGraph(vecPath+"meanGraph_"+is+"_"+run,  (GraphErrors)FitSummary.getItem( 1,is,run));
    		dumpGraph(vecPath+"res0Graph10_"+is+"_"+run,(GraphErrors)FitSummary.getItem( 2,is,run));
    		dumpGraph(vecPath+"resGraph_"+is+"_"+run,   (GraphErrors)FitSummary.getItem( 3,is,run));
    		dumpGraph(vecPath+"res2Graph_"+is+"_"+run,  (GraphErrors)FitSummary.getItem( 4,is,run));
    	}
    }
    
    public void fitGraphs() {
    	
    	ParallelSliceFitter fitter;
	    int run = getRunNumber();
	    
    	for (int is=1; is<7; is++) {
            tl.fitData.add(fitEngine(((H2F)this.getDataGroup().getItem(0,0,10,run).getData(is-1).get(0)).projectionY(),3,0.2,0.3,0.2,0.3,1.7,2.5),is,0,10,run); 
                       	
            fitter = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(0,0,10,run).getData(is-1).get(0));
            fitter.setBackgroundOrder(1); fitter.setMin(0.18); fitter.setMax(0.32); fitter.fitSlicesX();
            FitSummary.add(fitter.getMeanSlices(),is, 0, 10, run); // E/P vs. measured energy
               
            if (!dropSummary) {
            	
    		fitter = new ParallelSliceFitter((H2F) this.getDataGroup().getItem(0,0,0,run).getData(is-1).get(0));
    	    fitter.setBackgroundOrder(1); fitter.setMin(0.18); fitter.setMax(0.32); fitter.fitSlicesX();
    	    FitSummary.add(fitter.getMeanSlices(), is, 0, 1, run);  // E/P vs. tracking momentum
    	    
    	    GraphErrors meanGraph = fitter.getMeanSlices();  
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
            	   res2Graph.addPoint(xs, ys, 0., yse);
            	}
            }  
            System.out.println("Fit 4: "+is+" "+sigGraph.getDataSize(0));
            FitSummary.add(res2Graph, is, 0, 4, run);  
            }
    	}    	
    }
    
    public void createEOPHistos(int k, int xb, double x1, double x2, String txt1, String txt2, String txt3) {
    	
	    int run = getRunNumber();
        H2F h;  
        DataGroup dg = new DataGroup(3,2);
    
        for (int is=1; is<7; is++) {
            h = new H2F(txt1+"_"+is+"_"+k+"_"+run, xb, x1, x2, 50, 0., 0.5);
            h.setTitleX("Sector "+is+txt2);
            h.setTitleY(txt3);
            dg.addDataSet(h,is-1);
            this.getDataGroup().add(dg,0,0,k,run);
        }            
    } 
    
    public void createXYZHistos(int k) {
    	    
        int run = getRunNumber();
        H2F h;
        String tit1[] = {"PCAL EVENTS","ECIN EVENTS","ECOUT EVENTS"};
        String tit2[] = {"PCAL SF","ECIN SF","ECOUT SF"};
        String tit3[] = {"SF v PCAL","SF v ECIN","SF v ECOUT"};
        
        DataGroup dg1 = new DataGroup(3,3);
        DataGroup dg2 = new DataGroup(3,3);
        DataGroup dg3 = new DataGroup(3,3);
        
        for (int i=1; i<4; i++) {
            h = new H2F("ep_xyc_w"+i+"_"+run,  200, -200., 200., 200, -200.,200);                        dg2.addDataSet(h,i-1); 
            h = new H2F("ep_xyc_ww"+i+"_"+run, 200, -200., 200., 200, -200.,200);                        dg3.addDataSet(h,i-1); 
            h = new H2F("ep_xyc_e"+i+"_"+run,  200, -200., 200., 200, -200.,200); h.setTitle(tit1[i-1]); dg1.addDataSet(h,i-1); 
            h = new H2F("ep_xyc_sf"+i+"_"+run, 200, -200., 200., 200, -200.,200); h.setTitle(tit2[i-1]); dg1.addDataSet(h,i+2);  
            h = new H2F("ep_xyc_sff"+i+"_"+run,200, -200., 200., 200, -200.,200); h.setTitle(tit3[i-1]); dg1.addDataSet(h,i+5);  
        }
        this.getDataGroup().add(dg1, 0,0,k,run);
        this.getDataGroup().add(dg2, 0,1,k,run);
        this.getDataGroup().add(dg3, 0,2,k,run);
    }
    
    public void createADCHistos(int k, int nch, double x1, double x2, double x3, String txt) {
    	
	    int run = getRunNumber();
        H2F h;  
    
        for (int is=1; is<7; is++) {
            DataGroup dg = new DataGroup(3,3);
            h = new H2F("adc_pcal_u_"+is+"_"+k+"_"+run,"adc_pcal_u_"+is+"_"+k+"_"+run, nch, 0., x1, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("adc_pcal_v_"+is+"_"+k+"_"+run,"adc_pcal_v_"+is+"_"+k+"_"+run, nch, 0., x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("adc_pcal_w_"+is+"_"+k+"_"+run,"adc_pcal_w_"+is+"_"+k+"_"+run, nch, 0., x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
        
            h = new H2F("adc_ecin_u_"+is+"_"+k+"_"+run,"adc_ecin_u_"+is+"_"+k+"_"+run, nch, 0., x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("adc_ecin_v_"+is+"_"+k+"_"+run,"adc_ecin_v_"+is+"_"+k+"_"+run, nch, 0., x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("adc_ecin_w_"+is+"_"+k+"_"+run,"adc_ecin_w_"+is+"_"+k+"_"+run, nch, 0., x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
        
            h = new H2F("adc_ecou_u_"+is+"_"+k+"_"+run,"adc_ecou_u_"+is+"_"+k+"_"+run, nch, 0., x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("adc_ecou_v_"+is+"_"+k+"_"+run,"adc_ecou_v_"+is+"_"+k+"_"+run, nch, 0., x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("adc_ecou_w_"+is+"_"+k+"_"+run,"adc_ecou_w_"+is+"_"+k+"_"+run, nch, 0., x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,0,k,run);
        }            
    } 
    
/*   TIMELINES */
    
    public void createTimeLineHistos() {   
    	runIndex = 0;
    	createTimeLineHisto(10,"PCAL E/P","Sector",6,1,7);
    	createTimeLineHisto(20,"ECIN E/P","Sector",6,1,7);
    	createTimeLineHisto(40,"ECAL #sigma(E)/E","Sector",6,1,7);
    	createTimeLineHisto(30,"ECAL E/P","Sector",6,1,7);
    }
   
    public void createTimeLineHisto(int k, String tit, String ytit, int ny, int ymin, int ymax) {
    	H2F h1 =  new H2F(tit,tit,451, 0., 451., ny, ymin, ymax);
    	h1.setTitleX("Run Index") ; h1.setTitleY(ytit); tl.Timeline.add(h1,k,0); //mean
    	H2F h2 =  new H2F(tit,tit,451, 0., 451., ny, ymin, ymax);
    	h2.setTitleX("Run Index") ; h2.setTitleY(ytit); tl.Timeline.add(h2,k,1); //error
    }
    
    public void fillTimeLineHisto() {
        for (int is=1; is<7; is++) {
            float   y = (float) Fits.getItem(is,0,10,getRunNumber()).mean; 
            float  ye = (float) Fits.getItem(is,0,10,getRunNumber()).meane;			 
            float  ys = (float) Fits.getItem(is,0,10,getRunNumber()).sigma;
            float yse = (float) Fits.getItem(is,0,10,getRunNumber()).sigmae;			 
            ((H2F)tl.Timeline.getItem(30,0)).fill(runIndex,is,y);	
            ((H2F)tl.Timeline.getItem(30,1)).fill(runIndex,is,ye);   		
            ((H2F)tl.Timeline.getItem(40,0)).fill(runIndex,is,ys/y);	
            ((H2F)tl.Timeline.getItem(40,1)).fill(runIndex,is,(ys/y)*Math.sqrt(Math.pow(yse/ys,2)+Math.pow(ye/y,2)));   		
        } 
        runIndex++;
    }
    
    public void plotTimeLines(int index) {
    	if (getActivePC()==0) {plotClusterTimeLines(index);} else {plotPeakTimeLines(index);}
    }
    
    public void plotClusterTimeLines(int index) {
    	
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); 
        GraphErrors   g2 = null;
        int           is = getActiveSector(); 
        FitData       fd = tl.fitData.getItem(is,0,10,getRunNumber());
        
    	double min=0,max=0,mean=0.25f;
       
    	DataLine line1 = new DataLine(0,is,  runIndex+1,is);                   line1.setLineColor(5);
    	DataLine line2 = new DataLine(0,is+1,runIndex+1,is+1);                 line2.setLineColor(5);
    	DataLine line3 = new DataLine(runIndexSlider,1,  runIndexSlider,7);    line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,runIndexSlider+1,7);  line4.setLineColor(5);
    	DataLine line5 = new DataLine(-0.5,mean,runIndex,mean);                line5.setLineColor(3); line3.setLineWidth(2);
    	
        c.clear(); c.divide(3, 2); 

        for (int i=2; i<4; i++) { int i3=i*3-6;
            min=(i==3)?0.04f:0.22f; max=(i==3)?0.1f:0.27f; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(i3); c.getPad(i3).setAxisRange(0,runIndex,1,7); c.getPad(i3).setTitleFontSize(18); c.getPad(i3).getAxisZ().setRange(min,max);
    		c.draw((H2F)tl.Timeline.getItem((i+1)*10,0));c.draw(line1);c.draw(line2);c.draw(line3);c.draw(line4);
             
    		c.cd(i3+1); c.getPad(i3+1).setAxisRange(-0.5,runIndex,min,max); c.getPad(i3+1).setTitleFontSize(18);
    		List<GraphErrors> gglist = getGraph(((H2F)tl.Timeline.getItem((i+1)*10,0)),((H2F)tl.Timeline.getItem((i+1)*10,1)),is-1); 
               
    		for (int ii=1; ii<gglist.size(); ii++) {    
        		gglist.get(ii).setTitleX("Run Index"); gglist.get(ii).setTitleY("Sector "+is+((i==3)?" #sigma(E)/E":" E/P"));
    			c.draw(gglist.get(ii),(ii==1)?" ":"same"); c.draw(line5);	
    		}
    		g2 = new GraphErrors(); g2.setMarkerSize(5); g2.setMarkerColor(4); g2.setLineColor(2);
    		g2.addPoint(runIndexSlider,gglist.get(0).getDataY(runIndexSlider),0,0); c.draw(g2,"same");
    		    
    		c.cd(i3+2); c.getPad(i3+2).getAxisY().setRange(0.,fd.getGraph().getMax()*1.1);
            fd.getHist().getAttributes().setOptStat("1000100");
            DataLine line6 = new DataLine(mean,-50,mean,fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);
            c.draw(fd.getHist()); c.draw(fd.getGraph(),"same"); c.draw(line6);
        }
    }
     
    public void plotPeakTimeLines(int index) {
    	
    }
        
    @Override
    public void processEvent(DataEvent event) {
    	
    	int run = getRunNumber();
    	DataGroup dg = this.getDataGroup().getItem(0,0,0,run);
        
//        if (this.getNumberOfEvents() >= super.eventResetTime_current[5] && super.eventResetTime_current[5] > 0){
//            resetEventListener();
//        }
        
		trig_part_ind=-1; e_part_ind=-1; trig_sect=0;
		e_sect=0; e_ecal_E = 0;  e_pcal_e=0; e_etot_e=0;
		trig_track_ind = -1; e_track_ind = -1; pim_track_ind = -1; pip_track_ind = -1;
		pim_sect=0 ; pip_sect=0; pip_ecal_E=0; pim_ecal_E=0; pip_part_ind = -1; pim_part_ind = -1;
		
		for (int i=0; i<3; i++) {
			e_ecal_TH[i]=0; e_ecal_EL[i]=0; x_ecal[i]=1000; y_ecal[i]=1000; iU[i]=0; iV[i]=0; iW[i]=0;
			pim_ecal_TH[i]=0; pim_ecal_EL[i]=0; pim_x_ecal[i]=1000; pim_y_ecal[i]=1000; pim_iU[i]=0; pim_iV[i]=0; pim_iW[i]=0;
			pip_ecal_TH[i]=0; pip_ecal_EL[i]=0; pip_x_ecal[i]=1000; pip_y_ecal[i]=1000; pip_iU[i]=0; pip_iV[i]=0; pip_iW[i]=0;
			
		}

		if(dropBanks) dropBanks(event);
        
        boolean isMC = event.hasBank("MC::Particle");

        int trigger_sect = 0;
        
        if (!isMC) {
		// TRIGGER BIT SECTOR
      	trigger_sect = getElecTriggerSector(); 
      	
      	// HTCC*PCAL Q<0
//      	if(event.hasBank("REC::Particle"))  trig_part_ind = makeTrigElectron(event.getBank("REC::Particle"),event); 
      	
        // GET TB TRACK SECTOR, TRIG_TRACK_IND OF TRIG_PART_IND
//		if(event.hasBank("REC::Track") && 
//		   event.hasBank("TimeBasedTrkg::TBTracks")) getTrigTBTrack(event.getBank("TimeBasedTrkg::TBTracks"),event.getBank("REC::Track")); 
		
		// IS TRIG_PART_IND ID=11?
		if(event.hasBank("REC::Particle")) e_part_ind = makeElectron(event.getBank("REC::Particle")); 
		
		// FIND PION FOR COMPARISON
//		if(event.hasBank("REC::Particle")) makePiPlusPimPID(event.getBank("REC::Particle"));

		if(e_part_ind==-1)return;	
	
		// GET TRACK INDEX of E_PART_IND
//		if(event.hasBank("REC::Track")) fillEBTrack(event.getBank("REC::Track"));	
		
//        LorentzVector VGS = new LorentzVector(0,0,0,0);
//        VGS.add(VB);
//        VGS.sub(Ve);
//        e_Q2 = (float) -VGS.mass2();
//        e_xB = e_Q2/(2f*0.93827f*(Ebeam-e_mom));
//        e_W  = (float) Math.sqrt(0.93827f*0.93827f + e_Q2*(1f/e_xB-1f) );
        
        // GET PCAL, ECin, ECou ENERGY ETC. OF E_PART_IND   
		if(event.hasBank("REC::Calorimeter")&&event.hasBank("ECAL::clusters")) {
			getElecEBECal(event.getBank("REC::Calorimeter"),event.getBank("ECAL::clusters"));
//			getPiEBECal(event.getBank("REC::Calorimeter"),event.getBank("ECAL::clusters"));
		}
		
        }
        
        
        if (isMC&&event.hasBank("ECAL::clusters")) {Ebeam=10.6f; trig_track_ind=0; trig_sect=2; e_mom = 7f; getElecEBECalMC(event.getBank("ECAL::clusters"));}
		
		// GET SECTOR OF TRACK INDEX
//		if(event.hasBank("TimeBasedTrkg::TBTracks")) getTBTrack(event.getBank("TimeBasedTrkg::TBTracks"));
		
//        System.out.println(e_Q2+" "+e_W+" "+e_mom+" "+e_ecal_E+" "+trig_track_ind+" "+e_sect);
		
        float    sf = e_ecal_E/e_mom;
        float[] sff = new float[3];
        
        sff[0] = e_ecal_EL[0]/e_mom;
        sff[1] = e_ecal_EL[1]/e_mom;
        sff[2] = e_ecal_EL[2]/e_mom;
        
        trig_track_ind = 0;
        boolean good_e = e_sect>0&&e_sect<7, good_pim = pim_sect>0&&pim_sect<7, good_pip = pip_sect>0&&pip_sect<7 ; 
        
        trig_sect = trigger_sect; //when bypassing makeTrigElectron
        Boolean good_fiduc = iU[0]>7&&iV[0]<61&&iW[0]<61&&iU[1]>2&&iV[1]<36&&iW[1]<36&&iU[2]>2&&iV[2]<36&&iW[2]<36;
        
		if(e_mom>Ebeam*0.02 && sf > 0.02 && trig_track_ind>-1 && e_sect==trig_sect){
			if(good_e){	
				if(good_fiduc) {
				((H2F) this.getDataGroup().getItem(0,0,0,run).getData(e_sect-1).get(0)).fill(e_mom,sf);
				((H2F) this.getDataGroup().getItem(0,0,10,run).getData(e_sect-1).get(0)).fill(e_ecal_E,sf);
			    }
				((H2F) this.getDataGroup().getItem(0,0,1,run).getData(e_sect-1).get(0)).fill(e_theta,sf);
				((H2F) this.getDataGroup().getItem(0,0,2,run).getData(e_sect-1).get(0)).fill(e_ecal_TH[0],sf);
			}
			for (int i=0; i<3; i++) {
				if(good_e){
					((H2F) this.getDataGroup().getItem(0,0,3+i,run).getData(e_sect-1).get(0)).fill(e_ecal_TH[i],sff[i]);
					((H2F) this.getDataGroup().getItem(0,0,6,run).getData(i).get(0)).fill(-x_ecal[i], y_ecal[i],sff[i]<0.5?1f:0);
					((H2F) this.getDataGroup().getItem(0,1,6,run).getData(i).get(0)).fill(-x_ecal[i], y_ecal[i],sff[i]<0.5?sff[i]:0.);
					((H2F) this.getDataGroup().getItem(0,2,6,run).getData(i).get(0)).fill(-x_ecal[i], y_ecal[i],sf<0.5?sf:0.);
					((H2F) this.getDataGroup().getItem(e_sect,0,7,run).getData(3*i+0).get(0)).fill(sff[i]<0.5?sff[i]:0., iU[i]);
					((H2F) this.getDataGroup().getItem(e_sect,0,7,run).getData(3*i+1).get(0)).fill(sff[i]<0.5?sff[i]:0., iV[i]);
					((H2F) this.getDataGroup().getItem(e_sect,0,7,run).getData(3*i+2).get(0)).fill(sff[i]<0.5?sff[i]:0., iW[i]);				  
				}
				if(good_pim) {
					((H2F) this.getDataGroup().getItem(pim_sect,0,8,run).getData(3*i+0).get(0)).fill(pim_ecal_EL[i], pim_iU[i]);
					((H2F) this.getDataGroup().getItem(pim_sect,0,8,run).getData(3*i+1).get(0)).fill(pim_ecal_EL[i], pim_iV[i]);
					((H2F) this.getDataGroup().getItem(pim_sect,0,8,run).getData(3*i+2).get(0)).fill(pim_ecal_EL[i], pim_iW[i]);				  
				}
				if(good_pip) {
					((H2F) this.getDataGroup().getItem(pip_sect,0,9,run).getData(3*i+0).get(0)).fill(pip_ecal_EL[i], pip_iU[i]);
					((H2F) this.getDataGroup().getItem(pip_sect,0,9,run).getData(3*i+1).get(0)).fill(pip_ecal_EL[i], pip_iV[i]);
					((H2F) this.getDataGroup().getItem(pip_sect,0,9,run).getData(3*i+2).get(0)).fill(pip_ecal_EL[i], pip_iW[i]);						  
				}
			}
		}
       
    }
    
	public int makePiPlusPimPID(DataBank bank){
		
		for(int k = 0; k < bank.rows(); k++){
			int pid = bank.getInt("pid", k);				
				if(pid==211){
					float px = bank.getFloat("px", k);
					float py = bank.getFloat("py", k);
					float pz = bank.getFloat("pz", k);
					pip_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
					pip_theta = (float)Math.toDegrees(Math.acos(pz/pip_mom));
					pip_phi = (float)Math.toDegrees(Math.atan2(py,px));
					pip_vx = bank.getFloat("vx", k);
					pip_vy = bank.getFloat("vy", k);
					pip_vz = bank.getFloat("vz", k);
					pip_beta = bank.getFloat("beta", k);
					if(pip_mom>0.4 && pip_theta<40 && pip_theta>6 && pip_beta>0){
						VPIP = new LorentzVector(px,py,pz,Math.sqrt(pip_mom*pip_mom+0.139*0.139));
						pip_part_ind = k;
					}
				}
				
				if(pid==-211){
					float px = bank.getFloat("px", k);
					float py = bank.getFloat("py", k);
					float pz = bank.getFloat("pz", k);
					pim_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
					pim_theta = (float)Math.toDegrees(Math.acos(pz/pim_mom));
					pim_phi = (float)Math.toDegrees(Math.atan2(py,px));
					pim_vx = bank.getFloat("vx", k);
					pim_vy = bank.getFloat("vy", k);
					pim_vz = bank.getFloat("vz", k);
					pim_beta = bank.getFloat("beta", k);
					//System.out.println(pim_mom+" , "+pim_theta+" , "+pim_beta);
					if(pim_mom>0.4 && pim_theta<40 && pim_theta>6 && pim_beta>0){
						VPIM = new LorentzVector(px,py,pz,Math.sqrt(pim_mom*pim_mom+0.139*0.139));
						pim_part_ind = k;
					}
				}				
		}

		return -1;
	}
    
	public int makeElectron(DataBank bank){
		for(int k = 0; k < bank.rows(); k++){
			int pid = bank.getInt("pid", k);
			byte q = bank.getByte("charge", k);
			float px = bank.getFloat("px", k);
			float py = bank.getFloat("py", k);
			float pz = bank.getFloat("pz", k);
			e_mom = (float)Math.sqrt(px*px+py*py+pz*pz);			
			e_theta = (float)Math.toDegrees(Math.acos(pz/e_mom));
			e_vz = bank.getFloat("vz", k);
			
			if( pid == 11 && e_mom>0.01*Ebeam && e_mom<EB && e_theta>6 && Math.abs(e_vz)<200 ){
				e_phi = (float)Math.toDegrees(Math.atan2(py,px));
				e_vx = bank.getFloat("vx", k);
				e_vy = bank.getFloat("vy", k);
				Ve = new LorentzVector(px,py,pz,e_mom);
				return k;
			}
		}
		return -1;
	}
	
	public int makeTrigElectron(DataBank bank, DataEvent event){
		
		for(int k = 0; k < bank.rows(); k++){
			int  pid = bank.getInt("pid", k);
			byte   q = bank.getByte("charge", k);
			float px = bank.getFloat("px", k);
			float py = bank.getFloat("py", k);
			float pz = bank.getFloat("pz", k);
			e_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
			e_theta = (float)Math.toDegrees(Math.acos(pz/e_mom));
			e_vz = bank.getFloat("vz", k);
			
			if( pid == 11 && e_mom>0.02*Ebeam && e_mom<EB && e_theta>6 && Math.abs(e_vz)<20 ){}
			
			if( q<0 && e_theta>6 ) {
				float e_ecal_E=0;
				if(event.hasBank("REC::Calorimeter")) {
					DataBank ECALbank = event.getBank("REC::Calorimeter");
					for(int l = 0; l < ECALbank.rows(); l++) {
//						System.out.println(l+" "+ECALbank.getInt("layer",l)+" "+k+" "+ECALbank.getShort("pindex",l)+" "+ECALbank.getByte("sector",l));
						if(ECALbank.getShort("pindex",l)==k) {
							if(ECALbank.getInt("layer",l)==1) {
								trig_sect = ECALbank.getByte("sector",l);
							}
//							System.out.println(x);
					         e_ecal_E += ECALbank.getFloat("energy",l);
						}
					}
				}
				float HTCCnphe = 0;
				if(event.hasBank("REC::Cherenkov")){
					DataBank HTCCbank = event.getBank("REC::Cherenkov");
					for(int l = 0; l < HTCCbank.rows(); l++){
						if(HTCCbank.getShort("pindex",l)==k && HTCCbank.getInt("detector",l)==15){
							HTCCnphe = HTCCbank.getFloat("nphe",l); //Change to FLOAT
//							HTCCnphe = HTCCbank.getShort("nphe",l); //Change to FLOAT
						}
					}
				}
//				System.out.println("TrigSect= "+trig_sect+" E= "+e_ecal_E+" "+HTCCnphe+" "+e_mom);
				if( HTCCnphe>1 && e_ecal_E/e_mom > 0.18){}
				if( HTCCnphe>1 && e_ecal_E/e_mom > 0.02){
					e_phi = (float)Math.toDegrees(Math.atan2(py,px));
					e_vx = bank.getFloat("vx", k);
					e_vy = bank.getFloat("vy", k);
					Ve = new LorentzVector(px,py,pz,e_mom);
					return k;
				}
			}
		}
		return -1;
	}	

    public int getDet(int layer) {
	    int[] il = {0,0,0,1,1,1,2,2,2}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
	    return il[layer-1];
	}	
    
    public void getElecEBECalMC(DataBank clust) {
    	
		for(int k = 0; k < clust.rows(); k++){
			int det = clust.getInt("layer", k);
            int sec = clust.getByte("sector", k);
            float x = clust.getFloat("x",k);
            float y = clust.getFloat("y",k);
            float z = clust.getFloat("z",k);					         
            float r = (float) Math.sqrt(x*x+y*y+z*z);
			if(sec==2){				
				 int ind = getDet(det);				
		         x_ecal[ind] = x; y_ecal[ind] = y;
		         e_ecal_TH[ind] = (float) Math.toDegrees(Math.acos(z/r));
		         float ee = clust.getFloat("energy", k);
		         e_ecal_EL[ind] += ee;
                 e_ecal_E       += ee;
                 if(ind==0) e_sect = sec;
                 iU[ind] = (clust.getInt("coordU", k)-4)/8+1;
                 iV[ind] = (clust.getInt("coordV", k)-4)/8+1;
                 iW[ind] = (clust.getInt("coordW", k)-4)/8+1;        
			}
		}    	
    }
	   
	public void getElecEBECal(DataBank bank, DataBank clust){
		
		for(int k = 0; k < bank.rows(); k++){
			   int det = bank.getInt("layer", k);
			short pind = bank.getShort("pindex",k);
			short bind = bank.getShort("index",k);
               float x = bank.getFloat("x",k);
               float y = bank.getFloat("y",k);
               float z = bank.getFloat("z",k);					         
               float r = (float) Math.sqrt(x*x+y*y+z*z);
               float e = bank.getFloat("energy",k);

			if(pind==e_part_ind){				
				 int ind = getDet(det);				
		         x_ecal[ind] = x; y_ecal[ind] = y;
		         e_ecal_TH[ind] = (float) Math.toDegrees(Math.acos(z/r));
		         float ee = clust.getFloat("energy", bind);
		         e_ecal_EL[ind] += ee;
                 e_ecal_E       += ee;
                 if(det==1) e_sect = bank.getByte("sector",k);
                 
                 iU[ind] = (clust.getInt("coordU", bind)-4)/8+1;
                 iV[ind] = (clust.getInt("coordV", bind)-4)/8+1;
                 iW[ind] = (clust.getInt("coordW", bind)-4)/8+1;        
			}

		}
	}
	
	public void getPiEBECal(DataBank bank, DataBank clust){
		
		for(int k = 0; k < bank.rows(); k++){
			   int det = bank.getInt("layer", k);
			short pind = bank.getShort("pindex",k);
			short bind = bank.getShort("index",k);
               float x = bank.getFloat("x",k);
               float y = bank.getFloat("y",k);
               float z = bank.getFloat("z",k);					         
               float r = (float) Math.sqrt(x*x+y*y+z*z);
               float e = bank.getFloat("energy",k)*1000;            
			if(pind==-211){				
				 int ind = getDet(det);				
		         pim_x_ecal[ind] = x; pim_y_ecal[ind] = y;
		         pim_ecal_TH[ind] = (float) Math.toDegrees(Math.acos(z/r));
		         pim_ecal_EL[ind] += e;
                 pim_ecal_E       += e;
                 if(det==0) pim_sect = bank.getByte("sector",k);
                 
                 pim_iU[ind] = (clust.getInt("coordU", bind)-4)/8+1;
                 pim_iV[ind] = (clust.getInt("coordV", bind)-4)/8+1;
                 pim_iW[ind] = (clust.getInt("coordW", bind)-4)/8+1;        
			}
			if(pind==211){				
				 int ind = getDet(det);				
		         pip_x_ecal[ind] = x; pim_y_ecal[ind] = y;
		         pip_ecal_TH[ind] = (float) Math.toDegrees(Math.acos(z/r));
		         pip_ecal_EL[ind] +=  e;
                 pip_ecal_E       +=  e;
                 if(det==0) pip_sect = bank.getByte("sector",k);
                
                pip_iU[ind] = (clust.getInt("coordU", bind)-4)/8+1;
                pip_iV[ind] = (clust.getInt("coordV", bind)-4)/8+1;
                pip_iW[ind] = (clust.getInt("coordW", bind)-4)/8+1;        
			}

		}
	}	
	public void fillEBTrack(DataBank bank){
		e_track_ind=-1;pip_track_ind=-1;pim_track_ind=-1;
		for(int k = 0; k < bank.rows(); k++){
			short pind = bank.getShort("pindex",k);
			if(pind==e_part_ind){
				e_track_chi2 = 	bank.getFloat("chi2",k);
				e_track_ind = bank.getShort("index",k);
			}
			if(pind==pip_part_ind){
				pip_track_chi2 = bank.getFloat("chi2",k);
				pip_track_ind = bank.getShort("index",k);
				//System.out.println("fillEBTrack found pip track "+pip_track_ind);
			}
			if(pind==pim_part_ind){
				pim_track_chi2 = bank.getFloat("chi2",k);
				pim_track_ind = bank.getShort("index",k);
				//System.out.println("fillEBTrack found pim track "+pim_track_ind);
			}
		}
		//System.out.println("fillEBTrack : "+pim_part_ind+" , "+pip_part_ind+" ; "+pim_track_ind+" , "+pip_track_ind);
	}
		
	public void getTBTrack(DataBank bank){ 
		 if(e_track_ind>-1){
			 e_track_chi2 = bank.getFloat("chi2" , e_track_ind);
			 e_sect = bank.getInt("sector", e_track_ind);
		 }
		 if(pip_track_ind>-1)pip_sect = bank.getInt("sector", pip_track_ind);
		 if(pim_track_ind>-1)pim_sect = bank.getInt("sector", pim_track_ind);
	}
	
    public void getTrigTBTrack(DataBank bank, DataBank recBank){
    	
        for(int k = 0; k < bank.rows(); k++){
        	    if(recBank.getShort("pindex",k)==trig_part_ind) trig_track_ind = recBank.getShort("index",k);
        }
        
        if(trig_track_ind>-1 && trig_sect == bank.getInt("sector", trig_track_ind)){
                e_track_chi2 = bank.getFloat("chi2" , trig_track_ind);
                e_sect = bank.getInt("sector", trig_track_ind);
                e_DCR1_X = bank.getFloat("c1_x" , trig_track_ind);
                e_DCR1_Y = bank.getFloat("c1_y" , trig_track_ind);
                e_DCR1_Z = bank.getFloat("c1_z" , trig_track_ind);
                e_DCR3_X = bank.getFloat("c3_x" , trig_track_ind);
                e_DCR3_Y = bank.getFloat("c3_y" , trig_track_ind);
                e_DCR3_Z = bank.getFloat("c3_z" , trig_track_ind);
                e_DCR2_X = bank.getFloat("t1_x" , trig_track_ind);
                e_DCR2_Y = bank.getFloat("t1_y" , trig_track_ind);
                e_DCR2_Z = bank.getFloat("t1_z" , trig_track_ind);
	            Vector3 DCR1POS = new Vector3(e_DCR2_X,e_DCR2_Y,e_DCR2_Z);
	            Vector3 DCR1DIR = new Vector3(bank.getFloat("t1_px" , trig_track_ind),bank.getFloat("t1_py" , trig_track_ind),bank.getFloat("t1_pz" , trig_track_ind));
	            DCR1POS.rotateZ( -3.141597f*(e_sect-1)/3f );
	            DCR1DIR.rotateZ( -3.141597f*(e_sect-1)/3f );
	            float er1X = (float)DCR1POS.x();
	            float er1Y = (float)DCR1POS.y();
	            float er1Z = (float)DCR1POS.z();
	            float er1dX = (float)DCR1DIR.x();
	            float er1dY = (float)DCR1DIR.y();
	            float er1dZ = (float)DCR1DIR.z();
	            e_Ivy = er1Y + er1dY * (0f-er1X) / er1dX;
	            e_Ivz = er1Z + er1dZ * (0f-er1X) / er1dX;
	            float checkPh1 = (float)Math.toDegrees(Math.atan2( er1dY , er1dX ));
	            float checkTh1 = (float)Math.toDegrees(Math.acos( er1dZ / Math.sqrt( er1dX*er1dX+er1dY*er1dY+er1dZ*er1dZ ) ));
	            float checkPh2 = (float)Math.toDegrees(Math.atan2( e_Ivy-er1Y , -er1X ));
	            float checkTh2 = (float)Math.toDegrees(Math.acos( (e_Ivz-er1Z) / Math.sqrt( er1X*er1X + (e_Ivy-er1Y)*(e_Ivy-er1Y) + (e_Ivz-er1Z)*(e_Ivz-er1Z) )  ));

        }
    }

    public static class SFFunction extends Func1D{
    	
        EBCCDBConstants ccdb = new EBCCDBConstants();
   	    DetectorParticle p = new DetectorParticle(); 
        int pid;
  	    
        public SFFunction(String name, int pid, EBCCDBConstants ccdb, double min, double max) {
            super(name, min, max);
            this.ccdb = ccdb;
            this.pid  = pid;
            
            p.addResponse(new CalorimeterResponse(1,1,0));
            p.getDetectorResponses().get(0).getDescriptor().setType(DetectorType.ECAL);
        }
        @Override
        public double evaluate(double x){        	 
        	 p.getDetectorResponses().get(0).setEnergy(x);
       	     return  SamplingFractions.getMean(pid, p, ccdb);
        }
    }

    @Override
    public void timerUpdate() {    	
        for(int i=0; i<3; i++) {
        H2F e  = (H2F) this.getDataGroup().getItem(0,0,6,getRunNumber()).getData(i).get(0);
        H2F w  = (H2F) this.getDataGroup().getItem(0,1,6,getRunNumber()).getData(i).get(0);
        H2F ww = (H2F) this.getDataGroup().getItem(0,2,6,getRunNumber()).getData(i).get(0);
        for(int loop = 0; loop < e.getDataBufferSize(); loop++) {
        	    float ne = e.getDataBufferBin(loop);
            if (ne>0) {H2F h = (H2F) this.getDataGroup().getItem(0,0,6,getRunNumber()).getData(i+3).get(0); h.setDataBufferBin(loop,w.getDataBufferBin(loop)/ne);}
            if (ne>0) {H2F h = (H2F) this.getDataGroup().getItem(0,0,6,getRunNumber()).getData(i+6).get(0); h.setDataBufferBin(loop,ww.getDataBufferBin(loop)/ne);}
        }
        }
    }
    
}
