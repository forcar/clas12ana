package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingUtilities;

import org.clas.viewer.DetectorMonitor;
import org.dom4j.CDATA;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
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

    ECEngine engine = new ECEngine();
    IndexedList<List<Float>> tdcs = new IndexedList<List<Float>>(3);
    IndexedTable time = null;
    String[] layer = new String[]{"pcal","ecin","ecou"};
    String[]  view = new String[]{"u","v","w"};   
    int[]    npmts = new int[]{68,36,36};
    
    int trigger_sect = 0;
    float phase = 0;
    
//  static float TOFFSET = 436; 
    static float TOFFSET = 125;
    static float veff = 18.1f;
    static float    c = 29.98f;
    float      tps =  (float) 0.02345;
    float[] shiftTV = {0,40,0,0,0,0}; //Run 3050 t0 calibration
    
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
        		                     "FADC-TIME",
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
        		                     "LEFF v TTW",
        		                     "LEFF v TVERT");       
        this.useSectorButtons(true);
        this.useSliderPane(true);
        engine.init();       
        engine.setVariation("default");
        engine.setVeff(veff);
        time = engine.getConstantsManager().getConstants(3050, "/calibration/ec/timing");
        this.init();
    }
    
    @Override
    public void createHistos(int run) {  
    	    setRunNumber(run);
        createTDCHistos(0,150,350,"TIME (ns)");
        createTDCHistos(1,150,350,"TIME (ns)");
        createTDCHistos(2,150,350,"TIME (ns)");
        createTDCHistos(3,150,350,"TIME (ns)");
        createTDCHistos(4,150,350,"TIME (ns)");    
        createTDCHistos(5,150,350,"TIME (ns)");    
        createTDCHistos(6,150,350,"TIME (ns)");    
        createTDCHistos(7,150,350,"TIME (ns)");    
        createTDCHistos(8,-50.,50.,"FADC-TIME (ns)");    
        createTDCHistos(9,  0.,50.,"T-TVERT (ns)");    
        createTDCHistos(10,-10.,10.,"T-TVERT-PATH/c (ns)"); 
        createUVWHistos(11,700,800,-5,5,"PATH","RESID ");
        createUVWHistos(12,0,3000,-5,5,"ENERGY","RESID ");
        createUVWHistos(13,0,430,-5,5,"LEFF","RESID ");
        createUVWHistos(14,190,220,-5,5,"T","RESID ");
        createUVWHistos(15,0,6000,-5,5,"ADC","RESID ");
        createUVWHistos(16,180,240,0,6000,"T","ADC ");
        createUVWHistos(17,180,240,0,6000,"TTW","ADC ");
        createUVWHistos(18,180,240,0,12000,"TTW","ADCHI ");
        createUVWHistos(19,200,270,0,430,"TTW","LEFF ");
        createUVWHistos(20,160,190,0,430,"TVERT","LEFF ");
    }

    @Override        
    public void plotHistos(int run) {  
    	    setRunNumber(run);
    	    plotTDCHistos(0);
    	    plotTDCHistos(1);    	    	    
    	    plotTDCHistos(2);    	    	    
    	    plotTDCHistos(3);    	    	    
    	    plotTDCHistos(4);    	    	    
    	    plotTDCHistos(5);    	    	    
    	    plotTDCHistos(6);    	    	    
    	    plotTDCHistos(7);    	    	    
    	    plotTDCHistos(8);    	    	    
    	    plotTDCHistos(9);    	    	    
    	    plotTDCHistos(10); 
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
    }

    @Override
    public void processEvent(DataEvent event) {  
    	
       trigger_sect = getElecTriggerSector(); 
       phase        = getTriggerPhase()*4;
       processRaw(event);
       processRec(event);
             
    }
    
    public void processRaw(DataEvent event) {
    	
 	   int run = getRunNumber();
       H2F h;
       float tdcd,tdcdc =  0;
        
       tdcs.clear();
        
       int trigger_sect = getElecTriggerSector(); 
       float      phase = getTriggerPhase()*4;
        
       if(event.hasBank("ECAL::tdc")==true){
           DataBank  bank = event.getBank("ECAL::tdc");
           int rows = bank.rows();
           for(int i = 0; i < rows; i++){
               int  is = bank.getByte("sector",i);
               int  il = bank.getByte("layer",i);
               int  ip = bank.getShort("component",i);               
               tdcd    = bank.getInt("TDC",i)*tps;
               tdcdc   = tdcd-phase;
               if(is>0&&is<7&&tdcd>0) {
                   if(!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Float>(),is,il,ip);    
              	       h = (H2F) this.getDataGroup().getItem(is,0,0,run).getData(il-1).get(0); h.fill(tdcd, ip);
              	       h = (H2F) this.getDataGroup().getItem(is,0,1,run).getData(il-1).get(0); h.fill(tdcdc,ip);
                       if (is==trigger_sect) {
                    	       tdcs.getItem(is,il,ip).add(tdcdc);
                  	       h = (H2F) this.getDataGroup().getItem(is,0,2,run).getData(il-1).get(0); h.fill(tdcdc,ip);
                       }
               }
           }
       } 
        
       if(event.hasBank("ECAL::adc")==true){
           DataBank  bank = event.getBank("ECAL::adc");
           int rows = bank.rows();
           for(int i = 0; i < rows; i++){
               int  is = bank.getByte("sector",i);
               int  il = bank.getByte("layer",i);
               int  ip = bank.getShort("component",i);
               int adc = Math.abs(bank.getInt("ADC",i));
               float t = bank.getFloat("time",i);               
               
               Float[] tdcc; float[] tdc; float tmax = 1000; float tdcm = 1000;
               
               if (tdcs.hasItem(is,il,ip)) {
        	           double a0 = time.getDoubleValue("a0", is, il, ip);
        	           double a2 = time.getDoubleValue("a2", is, il, ip);
                   List<Float> list = new ArrayList<Float>();
                   list = tdcs.getItem(is,il,ip); tdcc=new Float[list.size()]; list.toArray(tdcc);
                   tdc  = new float[list.size()];
                   for (int ii=0; ii<tdcc.length; ii++) {
              	      float tdif = (tdcc[ii]-TOFFSET)-t; 
             	      h = (H2F) this.getDataGroup().getItem(is,0,8,run).getData(il-1).get(0); h.fill(tdif,ip);
            	          if (Math.abs(tdif)<30&&tdif<tmax) {tmax = tdif; tdcm = tdcc[ii];}                	    
                   }
                   double tdcmc = tdcm - a0 - a2 / Math.sqrt(adc);
          	       h = (H2F) this.getDataGroup().getItem(is,0,3,run).getData(il-1).get(0); h.fill(tdcm,ip);
          	       h = (H2F) this.getDataGroup().getItem(is,0,4,run).getData(il-1).get(0); h.fill(tdcmc,ip);
               }
           }
       }    	
    }
    
    public void processRec(DataEvent event) {
  	   
    	   int run = getRunNumber();
    	   
       event.removeBank("ECAL::hits");        
       event.removeBank("ECAL::peaks");        
       event.removeBank("ECAL::clusters");        
       event.removeBank("ECAL::calib");
        
       engine.processDataEvent(event); 
        
       List<ECStrip>     strips = engine.getStrips();
       List<ECPeak>       peaks = engine.getPeaks();
       List<ECCluster> clusters = engine.getClusters();
       
       H2F h;
       
       if(event.hasBank("ECAL::hits")){
          	DataBank  bank = event.getBank("ECAL::hits");
            for(int loop = 0; loop < bank.rows(); loop++){
               int   is = bank.getByte("sector", loop);
               int   il = bank.getByte("layer", loop); 
               int   ip = bank.getByte("strip", loop);
               float  t = bank.getFloat("time", loop)-phase;
               if (is==trigger_sect) {
            	     h = (H2F) this.getDataGroup().getItem(is,0,5,run).getData(il-1).get(0); h.fill(t, ip);
               }
            }
       }
       
       if(!(event.hasBank("REC::Calorimeter") && event.hasBank("REC::Particle"))) return;

       float Tvertex = event.hasBank("REC::Event") ? event.getBank("REC::Event").getFloat("STTime", 0):0;
       
       IndexedList<Integer> pathlist = new IndexedList<Integer>(3);    
       
       DataBank  bank = event.getBank("REC::Calorimeter");
       
       for(int loop = 0; loop < bank.rows(); loop++){
    	       int   is = bank.getByte("sector", loop);
    	       int   il = bank.getByte("layer", loop);
    	       int   in = bank.getShort("index", loop);
    	       int  det = bank.getByte("detector", loop);
    	       if (det==7 && !pathlist.hasItem(is,il,in)) pathlist.add(loop,is,il,in);                 
       }
       
       if(event.hasBank("ECAL::clusters")){
            DataBank  bank1 = event.getBank("ECAL::clusters");
            DataBank  bank2 = event.getBank("ECAL::calib");
            for(int loop = 0; loop < bank1.rows(); loop++){
                int is = bank1.getByte("sector", loop);
                if (is==trigger_sect){
                    int     il = bank1.getByte("layer", loop);
                    float ener = bank1.getFloat("energy",loop)*1000;
                    float    t = bank1.getFloat("time",loop)-phase;
                    int iU = (bank1.getInt("coordU", loop)-4)/8;
                    int iV = (bank1.getInt("coordV", loop)-4)/8;
                    int iW = (bank1.getInt("coordW", loop)-4)/8;
                    if (pathlist.hasItem(is,il,loop)) {
                    	    DataBank  bankc = event.getBank("REC::Calorimeter");
                    	    DataBank  bankp = event.getBank("REC::Particle");
                    	    int    pin = bankc.getShort("pindex", pathlist.getItem(is,il,loop));
                    	    float path = bankc.getFloat("path",   pathlist.getItem(is,il,loop));
                    	    for (int i=0; i<3; i++) {                    	 
                             float tu    = (float) clusters.get(loop).getTime(i)-phase; 
                             int  ip     =         clusters.get(loop).getPeak(i).getMaxStrip();
                             int  adc    =         clusters.get(loop).getPeak(i).getMaxECStrip().getADC();
                             float tdc   = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getRawTime()-phase;
                             float tdcc  = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTWCTime()-phase;
                             float tdccc = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTime()-phase; 
                             float leff  = (float) clusters.get(loop).getPeak(i).getMaxECStrip().getTdist();
                             float tdif  = tu  - Tvertex + phase;
                             float tdifp = tdif - path/c;
                             float texp  = path/c + leff/veff + Tvertex - phase - shiftTV[is-1];
                             h = (H2F) this.getDataGroup().getItem(is,0,6,run).getData(il+i-1).get(0); h.fill(tu, ip);
                             h = (H2F) this.getDataGroup().getItem(is,0,7,run).getData(il+i-1).get(0); h.fill(t,  ip);
                             h = (H2F) this.getDataGroup().getItem(is,0,9,run).getData(il+i-1).get(0); h.fill(tdif, ip);
                             if (bankp.getInt("pid",pin)==11) {
                                h = (H2F) this.getDataGroup().getItem(is,   0,10,run).getData(il+i-1).get(0);  h.fill(tdifp, ip);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,11,run).getData(ip-1).get(0); h.fill(path, tdifp);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,12,run).getData(ip-1).get(0); h.fill(ener, tdifp);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,13,run).getData(ip-1).get(0); h.fill(leff, tdifp);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,14,run).getData(ip-1).get(0); h.fill(tu  -shiftTV[is-1], tdifp);  
                                h = (H2F) this.getDataGroup().getItem(is,il+i,15,run).getData(ip-1).get(0); h.fill(adc,  tdifp);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,16,run).getData(ip-1).get(0); h.fill(tdc -shiftTV[is-1]-leff/veff, adc);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,17,run).getData(ip-1).get(0); h.fill(tdcc-shiftTV[is-1]-leff/veff, adc);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,18,run).getData(ip-1).get(0); h.fill(tdcc-shiftTV[is-1]-leff/veff, adc);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,19,run).getData(ip-1).get(0); h.fill(tdcc, leff);
                                h = (H2F) this.getDataGroup().getItem(is,il+i,20,run).getData(ip-1).get(0); h.fill(Tvertex-shiftTV[is-1]-phase,leff); 	
                             } 	      
                    	    } 
                    }
                }
            }
       }    	
    }
    
    public void createUVWHistos(int k, int xmin, int xmax, int ymin, int ymax, String xtxt, String ytxt) {
    	
        H2F h;  F1D f1; 
        
        int xbins=50, ybins=50; double sca1=1.0,sca2=1.0;
        int run = getRunNumber();
        
        if (k==12) {sca1=0.5; sca2=0.3;}
        
        for (int is=1; is<7; is++) {      
            DataGroup dg1 = new DataGroup(8,8); DataGroup dg2 = new DataGroup(8,8); DataGroup dg3 = new DataGroup(8,8);        	    
            f1 = new F1D("p0"+is+1+k,"[a]",xmin,xmax); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
           
            for (int ip=1; ip<npmts[0]+1; ip++) {
                h = new H2F("uvw_pcal_u"+ip+"_s"+is+"_"+k,"uvw_pcal_u"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"U"+ip);       
                dg1.addDataSet(h,ip-1); dg1.addDataSet(f1,ip-1);
                h = new H2F("uvw_pcal_v"+ip+"_s"+is+"_"+k,"uvw_pcal_v"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"V"+ip);
                dg2.addDataSet(h,ip-1); dg2.addDataSet(f1,ip-1);
                h = new H2F("uvw_pcal_w"+ip+"_s"+is+"_"+k,"uvw_pcal_w"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin,xmax,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" PCAL "+xtxt); h.setTitleY(ytxt+"W"+ip); 
                dg3.addDataSet(h,ip-1); dg3.addDataSet(f1,ip-1);
     	    }
            this.getDataGroup().add(dg1,is,1,k,run); this.getDataGroup().add(dg2,is,2,k,run); this.getDataGroup().add(dg3,is,3,k,run);
            
            DataGroup dg4 = new DataGroup(6,6); DataGroup dg5 = new DataGroup(6,6); DataGroup dg6 = new DataGroup(6,6);        	         	   
            f1 = new F1D("p0"+is+2+k,"[a]",xmin*sca1,xmax*sca1); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
     	    for (int ip=1; ip<npmts[1]+1; ip++) {
                h = new H2F("uvw_ecin_u"+ip+"_s"+is+"_"+k,"uvw_ecin_u"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca1,xmax*sca1,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECIN "+xtxt);  h.setTitleY(ytxt+"U"+ip); 
                dg4.addDataSet(h,ip-1); dg4.addDataSet(f1,ip-1);
                h = new H2F("uvw_ecin_v"+ip+"_s"+is+"_"+k,"uvw_ecin_v"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca1,xmax*sca1,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"V"+ip); 
                dg5.addDataSet(h,ip-1); dg5.addDataSet(f1,ip-1);
                h = new H2F("uvw_ecin_w"+ip+"_s"+is+"_"+k,"uvw_ecin_w"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca1,xmax*sca1,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECIN "+xtxt); h.setTitleY(ytxt+"W"+ip);
                dg6.addDataSet(h,ip-1); dg6.addDataSet(f1,ip-1);
                
     	    }
            this.getDataGroup().add(dg4,is,4,k,run); this.getDataGroup().add(dg5,is,5,k,run); this.getDataGroup().add(dg6,is,6,k,run);
     	   
            DataGroup dg7 = new DataGroup(6,6); DataGroup dg8 = new DataGroup(6,6); DataGroup dg9 = new DataGroup(6,6);        	         	   
            f1 = new F1D("p0"+is+3+k,"[a]",xmin*sca2,xmax*sca2); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
     	    for (int ip=1; ip<npmts[2]+1; ip++) {
                h = new H2F("uvw_ecou_u"+ip+"_s"+is+"_"+k,"uvw_ecou_u"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca2,xmax*sca2,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECOU "+xtxt); h.setTitleY(ytxt+"U"+ip);
                dg7.addDataSet(h,ip-1); dg7.addDataSet(f1,ip-1);
                h = new H2F("uvw_ecou_v"+ip+"_s"+is+"_"+k,"uvw_ecou_v"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca2,xmax*sca2,ybins,ymin,ymax);
                h.setTitleX("Sector "+is+" ECOU "+xtxt);  h.setTitleY(ytxt+"V"+ip);
                dg8.addDataSet(h,ip-1); dg8.addDataSet(f1,ip-1);
                h = new H2F("uvw_ecou_w"+ip+"_s"+is+"_"+k,"uvw_ecou_w"+ip+"_s"+is+"_"+k+"_"+run,xbins,xmin*sca2,xmax*sca2,ybins,ymin,ymax);
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
            h = new H2F("tdc_pcal_u_"+is+"_"+k,"tdc_pcal_u_"+is+"_"+k+"_"+run,100, tmin, tmax, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("tdc_pcal_v_"+is+"_"+k,"tdc_pcal_v_"+is+"_"+k+"_"+run,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("tdc_pcal_w_"+is+"_"+k,"tdc_pcal_w_"+is+"_"+k+"_"+run,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
            
            h = new H2F("tdc_ecin_u_"+is+"_"+k,"tdc_ecin_u_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("tdc_ecin_v_"+is+"_"+k,"tdc_ecin_v_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("tdc_ecin_w_"+is+"_"+k,"tdc_ecin_w_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
            
            h = new H2F("tdc_ecou_u_"+is+"_"+k,"tdc_ecou_u_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("tdc_ecou_v_"+is+"_"+k,"tdc_ecou_v_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("tdc_ecou_w_"+is+"_"+k,"tdc_ecou_w_"+is+"_"+k+"_"+run,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,0,k,run);
        }            

    } 
  
    public void plotUVWHistos(int index) {
    	    int run = getRunNumber();
//        getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)).draw(getDataGroup().getItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index));    	
	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),3*getActiveLayer()+getActiveView()+1,index,run));	    
    }

    
    public void plotTDCHistos(int index) {
    	    int run = getRunNumber();
//        getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)).draw(getDataGroup().getItem(getActiveSector(),0,index));
    	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),0,index,run));	    
    }
    
    @Override
    public void timerUpdate() {
    	
    }
     
}
