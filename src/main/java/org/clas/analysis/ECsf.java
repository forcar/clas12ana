package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.clas.tools.FitData;
import org.clas.tools.ParallelSliceFitter;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H2F;
//import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.eb.EBCCDBConstants;
import org.jlab.rec.eb.SamplingFractions;

public class ECsf extends DetectorMonitor {

    String[]  det = {"pcal","ecin","ecou"};
    int[]    npmt = {68,62,62,36,36,36,36,36,36};    
    String[]    v = new String[]{"u","v","w"};    
    int nelec=0;
    
    public double[][] par = {{0.105,0.039},{0.099,0.040},{0.100,0.034},{0.093,0.044},{0.085,0.046},{0.113,0.028}};
    
    public ECsf(String name) {
        super(name);
        this.setDetectorTabNames("E/P v P",
                                 "E/P v ThV",
                                 "E/P v ThD",
                                 "EPC/P v ThD", 
                                 "EECi/P v ThD",
                                 "EECo/P v ThD",
                                 "E/P v XY",
                                 "E/P v Em",
                                 "Sector Fits",
                                 "UVW",
                                 "Fits",
                                 "Summary",
                                 "Timeline",
                                 "Fits Sig 1",
                                 "Fits Sig 2",
                                 "Fits Sig 3",
                                 "Timing");
        this.usePCCheckBox(true);
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();
    }
    
    public void localinit() {
        System.out.println("ECa.localinit()");
        configEngine("muon");  
        configEventBuilder();
        tl.setFitData(Fits);
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
        createEOPHistos(7,0,50,0.,2.5,"ep_em", " Measured Energy (GeV)", " E/P");
        if(dropSummary) return;
        createEOPHistos(0,0,50,0.5,10.5,"ep_p", " Momentum (GeV)",      " E/P");
        createEOPHistos(0,1,50,0.5,10.5,"ep_pnf",  " Momentum (GeV)",   " E/P");
        createEOPHistos(1,0,30,  6.,36.,"ep_thv"," VertexTheta (deg)",   " E/P");
        createEOPHistos(2,0,48,  3.,37.,"ep_thd"," Detector Theta (deg)"," E/P");
        createEOPHistos(3,0,48,  3.,27.,"ep_th0"," PC Theta (deg)",      "EPC / P");
        createEOPHistos(4,0,48,  3.,27.,"ep_th1"," ECIN Theta (deg)",    "EECi / P");
        createEOPHistos(5,0,48,  3.,27.,"ep_th2"," ECOU Theta (deg)",    "EECo / P");
        createXYZHistos(6);
        createADCHistos(9,0,25,0.,0.5,0.3,0.1,"PARTIAL SF");
        createADCHistos(9,1,25,0.05,0.4,0.4,0.4,"TOTAL SF");
        createADCHistos(16,0,100,-10.,10.,10.,10.," T-TVERT-PATH/c (ns)");
    }

    @Override       
    public void plotHistos(int run) {
    	plotSummary(run);
    	plotAnalysis(run);
    }
      
    public void plotSummary(int run) {
    	setRunNumber(run);
    	plotEOPHistos(7);
    	if(dropSummary) return;
    	plotEOPHistos(0);
    	plotEOPHistos(1);
    	plotEOPHistos(2);
    	plotEOPHistos(3);
    	plotEOPHistos(4);
    	plotEOPHistos(5);
    	plotXYZHistos(6);
    	plotADCHistos(9);
    	plotADCHistos(16);
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;
    	plotFitSummary8(8);
    	if(!dropSummary) {updateFits(10);plotMeanHWSummary(11); plotFitSummary2(13); plotFitSummary3(14); plotFitSummary4(15);}
    	plotTimeLines(12);
    }
    
    public void createEOPHistos(int k, int pc, int xb, double x1, double x2, String txt1, String txt2, String txt3) {
    	
        int run = getRunNumber();
        H2F h;  
        DataGroup dg = new DataGroup(3,2);
    
        for (int is=1; is<7; is++) {
            h = new H2F(txt1+"_"+is+"_"+k+"_"+run, xb, x1, x2, 50, 0., 0.5);
            h.setTitleX("Sector "+is+txt2);
            h.setTitleY(txt3);
            dg.addDataSet(h,is-1);
            this.getDataGroup().add(dg,pc,0,k,run);
        }            
    } 
    
    public void createXYZHistos(int k) {
    	    
        int run = getRunNumber();
        H2F h;
        String tit1[] = {"PCAL EVENTS","ECIN EVENTS","ECOU EVENTS"};
        String tit2[] = {"PCAL PARTIAL SF","ECIN PARTIAL SF","ECOU PARTIAL SF"};
        String tit3[] = {"TOTAL SF PCAL","TOTAL SF ECIN","TOTAL SF ECOU"};
        
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
    
    public void createADCHistos(int k, int pc, int nch, double x0, double x1, double x2, double x3, String txt) {
    	
        int run = getRunNumber();
        H2F h;  
    
        for (int is=1; is<7; is++) {
            DataGroup dg = new DataGroup(3,3);
            h = new H2F("adc_pcal_u_"+is+"_"+pc+"_"+k+"_"+run,"adc_pcal_u_"+is+"_"+k+"_"+run, nch, x0, x1, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("adc_pcal_v_"+is+"_"+pc+"_"+k+"_"+run,"adc_pcal_v_"+is+"_"+k+"_"+run, nch, x0, x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("adc_pcal_w_"+is+"_"+pc+"_"+k+"_"+run,"adc_pcal_w_"+is+"_"+k+"_"+run, nch, x0, x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
        
            h = new H2F("adc_ecin_u_"+is+"_"+pc+"_"+k+"_"+run,"adc_ecin_u_"+is+"_"+k+"_"+run, nch, x0, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("adc_ecin_v_"+is+"_"+pc+"_"+k+"_"+run,"adc_ecin_v_"+is+"_"+k+"_"+run, nch, x0, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("adc_ecin_w_"+is+"_"+pc+"_"+k+"_"+run,"adc_ecin_w_"+is+"_"+k+"_"+run, nch, x0, x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
        
            h = new H2F("adc_ecou_u_"+is+"_"+pc+"_"+k+"_"+run,"adc_ecou_u_"+is+"_"+k+"_"+run, nch, x0, x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("adc_ecou_v_"+is+"_"+pc+"_"+k+"_"+run,"adc_ecou_v_"+is+"_"+k+"_"+run, nch, x0, x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("adc_ecou_w_"+is+"_"+pc+"_"+k+"_"+run,"adc_ecou_w_"+is+"_"+k+"_"+run, nch, x0, x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,pc,k,run);
        }  

    } 
    
    @Override
    public void processEvent(DataEvent event) {
    	
    	int run = getRunNumber();
    	DataGroup dg = this.getDataGroup().getItem(0,0,0,run);
    	
    	float Ebeam=10.6f, EB=10.6f, e_mom=0, e_theta=0, e_vz=0, e_ecal_E=0; 
    	float[]    x_ecal = new float[3];
        float[]    y_ecal = new float[3];
    	float[] e_ecal_TH = new float[3];
    	float[] e_ecal_EL = new float[4];
    	float[]    t_ecal = new float[4];
    	float[]   pa_ecal = new float[4];
    	int e_sect=0;
    	int[] iU = new int[3];int[] idU = new int[3]; float[] tU = new float[3];
    	int[] iV = new int[3];int[] idV = new int[3]; float[] tV = new float[3];
    	int[] iW = new int[3];int[] idW = new int[3]; float[] tW = new float[3];
    	
        if(dropBanks) dropBanks(event);
        
        if (event.hasBank("MC::Particle")) return;

        int trigger_sect = 0;
        
        boolean goodEvent = event.hasBank("REC::Particle")&&event.hasBank("REC::Calorimeter");
        
        if (!goodEvent) return;
        
        float Tvertex = event.hasBank("REC::Event") ? (isHipo3Event ? event.getBank("REC::Event").getFloat("STTime", 0):
                                                                      event.getBank("REC::Event").getFloat("startTime", 0)):0;        
      	DataBank    reccal = event.getBank("REC::Calorimeter");
      	DataBank ecalclust = event.getBank("ECAL::clusters");
      	DataBank ecalpeaks = event.getBank("ECAL::peaks");
        DataBank    recpar = event.getBank("REC::Particle");
        
      	Map<Integer,List<Integer>> caloMap = loadMapByIndex(reccal,"pindex");
      	Map<Integer,List<Integer>> partMap = loadMapByIndex(recpar,"pid");
      	
      	trigger_sect = getElecTriggerSector(); 
        int tpid = 11;
        
      	if (!(trigger_sect>0)) return;
      	
      	if(!partMap.containsKey(tpid)) return;
      	if(partMap.get(tpid).size()!=1) return;
      	
//	    nelec++;System.out.println("Evnt "+getEventNumber()+" Nelec "+nelec);
      	    
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
           if(inDC && ep>0.01*Ebeam && ep<EB && th>6 && Math.abs(vz)<200 ){
               for (int icalo : caloMap.get(ipart)) {
				    int  det = reccal.getInt("layer", icalo);
	                short ic = reccal.getShort("index",icalo);
                    float  x = reccal.getFloat("x",icalo);
                    float  y = reccal.getFloat("y",icalo);
                    float  z = reccal.getFloat("z",icalo);					         
                    float pa = reccal.getFloat("path", icalo);
                    float  t = reccal.getFloat("time",icalo);
                    float  r = (float) Math.sqrt(x*x+y*y+z*z);
                    if(det==1) e_sect = reccal.getByte("sector",icalo);	                 
                    int ind = getDet(det);	
                    x_ecal[ind]     = x;
                    y_ecal[ind]     = y;
                    t_ecal[ind]     = t-Tvertex-pa/29.98f;
                    e_ecal_TH[ind]  = (float) Math.toDegrees(Math.acos(z/r));	               
                    e_ecal_EL[ind] += ecalclust.getFloat("energy", ic);
                    e_ecal_EL[3]   += ecalclust.getFloat("energy", ic);
                    iU[ind]         = (ecalclust.getInt("coordU", ic)-4)/8+1;
                    iV[ind]         = (ecalclust.getInt("coordV", ic)-4)/8+1;
                    iW[ind]         = (ecalclust.getInt("coordW", ic)-4)/8+1; 
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
		
        if(e_sect!=trigger_sect) return;	
      
        float[] sff = new float[4];
        
        for (int i=0; i<4; i++) sff[i] = e_ecal_EL[i]/e_mom;
        float sf = sff[3];
        
        boolean     good_e = e_sect>0&&e_sect<7;         
        boolean good_fiduc = iU[0]>2&&iV[0]<63&&iW[0]<63&&iU[1]>2&&iV[1]<36&&iW[1]<36&&iU[2]>2&&iV[2]<36&&iW[2]<36;
        
        if(e_mom>Ebeam*0.02 && sff[3] > 0.02){			
           if(good_e){	
              if(good_fiduc) {
                 ((H2F) this.getDataGroup().getItem(0,0,0,run).getData(e_sect-1).get(0)).fill(e_mom,sf);
                 ((H2F) this.getDataGroup().getItem(0,0,7,run).getData(e_sect-1).get(0)).fill(e_ecal_EL[3],sf);				
              }
              ((H2F) this.getDataGroup().getItem(1,0,0,run).getData(e_sect-1).get(0)).fill(e_mom,sf);
              ((H2F) this.getDataGroup().getItem(0,0,1,run).getData(e_sect-1).get(0)).fill(e_theta,sf);
              ((H2F) this.getDataGroup().getItem(0,0,2,run).getData(e_sect-1).get(0)).fill(e_ecal_TH[0],sf);
           }
           for (int i=0; i<3; i++) {
             if(good_e){
                ((H2F) this.getDataGroup().getItem(0,0,3+i,run).getData(e_sect-1).get(0)).fill(e_ecal_TH[i],sff[i]);
                ((H2F) this.getDataGroup().getItem(0,0,6,run).getData(i).get(0)).fill(-x_ecal[i], y_ecal[i],sff[i]<0.5?1f:0);
                ((H2F) this.getDataGroup().getItem(0,1,6,run).getData(i).get(0)).fill(-x_ecal[i], y_ecal[i],sff[i]<0.5?sff[i]:0.);
                ((H2F) this.getDataGroup().getItem(0,2,6,run).getData(i).get(0)).fill(-x_ecal[i], y_ecal[i],sf<0.5?sf:0.);
                ((H2F) this.getDataGroup().getItem(e_sect,0,9,run).getData(3*i+0).get(0)).fill(sff[i]<0.5?sff[i]:0., iU[i]);
                ((H2F) this.getDataGroup().getItem(e_sect,0,9,run).getData(3*i+1).get(0)).fill(sff[i]<0.5?sff[i]:0., iV[i]);
                ((H2F) this.getDataGroup().getItem(e_sect,0,9,run).getData(3*i+2).get(0)).fill(sff[i]<0.5?sff[i]:0., iW[i]);				  
                ((H2F) this.getDataGroup().getItem(e_sect,1,9,run).getData(3*i+0).get(0)).fill(sf<0.5?sf:0., iU[i]);				  
                ((H2F) this.getDataGroup().getItem(e_sect,1,9,run).getData(3*i+1).get(0)).fill(sf<0.5?sf:0., iV[i]);				  
                ((H2F) this.getDataGroup().getItem(e_sect,1,9,run).getData(3*i+2).get(0)).fill(sf<0.5?sf:0., iW[i]);				  
                ((H2F) this.getDataGroup().getItem(e_sect,0,16,run).getData(3*i+0).get(0)).fill(t_ecal[i], iU[i]);				  
                ((H2F) this.getDataGroup().getItem(e_sect,0,16,run).getData(3*i+1).get(0)).fill(t_ecal[i], iV[i]);				  
                ((H2F) this.getDataGroup().getItem(e_sect,0,16,run).getData(3*i+2).get(0)).fill(t_ecal[i], iW[i]);				  
             }
           }
        }

    }
    
    public void updateFits(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));    	
        int    pc = 1;
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
        if(!isAnalyzeDone) createTimeLineHistos();
    	fillTimeLineHisto();
        System.out.println("Finished");
        isAnalyzeDone = true;    	
    } 
       
    public void fitGraphs(int is1, int is2, int id1, int id2, int il1, int il2) {    	
        int run = getRunNumber();	 
        
        // SF v UVW Fits
        for (int is=is1; is<is2; is++) {
           tl.fitData.add(fitEngine(((H2F)this.getDataGroup().getItem(0,0,7,run).getData(is-1).get(0)).projectionY(),0,0.15,0.3,0.15,0.3,1.7,2.5),is,0,7,run); 
           for (int id=id1; id<id2; id++) {
               for (int il=il1; il<il2; il++) {
               	  for (int pc=1; pc<2; pc++) {
                     H2F h = (H2F) this.getDataGroup().getItem(is,pc,9,run).getData(id*3+il).get(0);
           	         for (int i=0; i<npmt[id*3+il]; i++) tl.fitData.add(fitEngine(h.sliceY(i),0,0.15,0.3,0.15,0.3,1.7,2.5),is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); 
        	         fitStore(is, id, il, pc, run, 1f);
               	  }
               } 
           }
        }
        
        //SF v P Fits
    	for (int is=1; is<7; is++) { 
    		
            ParallelSliceFitter fitter;
              
            if (!dropSummary) {
            	
            fitter = new ParallelSliceFitter((H2F)this.getDataGroup().getItem(0,0,7,run).getData(is-1).get(0));
            fitter.setBackgroundOrder(1); fitter.setMin(0.18); fitter.setMax(0.32); fitter.fitSlicesX(); 
            FitSummary.add(fitter.getMeanSlices(),is, 0, 7, run); // E/P vs. measured energy
            
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
    
    public void plotMeanHWSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = 1;
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
    
    public void plotFitSummary8(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(3, 4);
        int col[] = {1,2,3,4,5,7};
        int run = getRunNumber();    
	   
	    SFFunction sf = new SFFunction("esf",-11,eb.ccdb,0.1,2.5); 

    	for (int is=1; is<7; is++) {  
    		String txt = "Sector "+is+" Measured Energy (GeV)";
            if (FitSummary.hasItem(is,0,7,run)) {
            	GraphPlot((GraphErrors)FitSummary.getItem(is,0,7,run),c,is-1,0.0f,2.5f,0.18f,0.300f,col[is-1],4,1,txt," E/P",""); c.draw(sf,"same");
            }
            c.cd(is-1+6);c.getPad(is-1+6).getAxisX().setRange(0.1, 0.4); 
            c.draw(Fits.getItem(is,0,7,getRunNumber()).getHist());
            c.draw(Fits.getItem(is,0,7,getRunNumber()).getGraph(),"same");     	
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
    	GraphErrors g1 = new GraphErrors(), g2 = new GraphErrors();
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.setGridX(false); c.setGridY(false); c.divide(6, 2);
        int col[] = {1,2,3,4,5,7};
        int run = getRunNumber();
     	for (int is=1; is<7; is++) {    		
//            F1D f = new F1D("res","[a]*[a]*x+[b]*[b]",0.,0.6); f.setLineColor(1); f.setLineWidth(3);
//            f.setParameter(0, par[is-1][0]);f.setParameter(1, par[is-1][1]);
//            if (FitSummary.hasItem(is,0,4,run)) GraphPlot((GraphErrors)FitSummary.getItem(is,0,4,run),c,is-1,0.f,0.6f,0.001f,0.008f,col[is-1],4,1,"","","");
              GraphPlot((GraphErrors)tl.fitData.getItem(is,0,4,run).getGraph(),c,is-1,0.f,0.6f,0.001f,0.008f,col[is-1],4,1,"","",""); 
              g1.addPoint(is,Math.sqrt(tl.fitData.getItem(is,0,4,run).p1),0,Math.sqrt(tl.fitData.getItem(is,0,4,run).p1e));
              g2.addPoint(is,Math.sqrt(tl.fitData.getItem(is,0,4,run).p0),0,Math.sqrt(tl.fitData.getItem(is,0,4,run).p0e));
//            if (FitSummary.hasItem(is,0,4,run)) c.draw(f,"same");            
    	}    	
        GraphPlot(g1,c,7,0.5f,6.5f,0.0f,0.11f,1,6,1,"SECTOR","",""); GraphPlot(g2,c,7,0.5f,6.5f,0.0f,0.11f,1,7,2,"","","same"); 
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
    
    public void plotADCHistos(int index) {
    	if(!getDataGroup().hasItem(getActiveSector(),getActivePC(),index,getRunNumber())) return;
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),getActivePC(),index,getRunNumber()));	       	
    }
    
    public void plotEOPHistos(int index) {
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActivePC(),0,index,getRunNumber()));	       	
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
    
    public void plotTimeLines(int index) {
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
