package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.viewer.DetectorMonitor;
import org.dom4j.CDATA;
import org.jlab.groot.data.H2F;
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

public class ECrec extends DetectorMonitor {

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
	float      tps =  (float) 0.02345;
	float[] shiftTV = {0,20,0,0,0,0}; //Run 3050 t0 calibration
    
    public ECrec(String name) {
        super(name);
        this.setDetectorTabNames("Raw TDC",
        		                     "PhaseCorr TDC",
        		                     "Triggered TDC",
        		                     "Matched TDC",
        		                     "Calib TDC",
        		                     "Hit Time",
        		                     "Peak Time",
        		                     "Cluster Time",
        		                     "fADC-TDC",
        		                     "TIME-STT",
        		                     "RESID v STRIP",
        		                     "RESID v PATH",
        		                     "RESID v ENERGY",
        		                     "LEFF v TDC ",
        		                     "DTEXP v TDC ",
        		                     "LEFF v TDCTWC ",
        		                     "DTEXP v TDCTWC ");

        
        this.useSectorButtons(true);
        this.useSliderPane(true);
        engine.init();       
        engine.setVariation("default");
        time = engine.getConstantsManager().getConstants(3050, "/calibration/ec/timing");
        this.init(false);
    }
    
    @Override
    public void createHistos() {  
        createTDCHistos(0,0,0,150,350,"TDC (ns)");
        createTDCHistos(0,0,1,150,350,"TDC (ns)");
        createTDCHistos(0,0,2,150,350,"TDC (ns)");
        createTDCHistos(0,0,3,150,350,"TDC (ns)");
        createTDCHistos(0,0,4,150,350,"TDC (ns)");    
        createTDCHistos(0,0,5,150,350,"TDC (ns)");    
        createTDCHistos(0,0,6,150,350,"TDC (ns)");    
        createTDCHistos(0,0,7,150,350,"TDC (ns)");    
        createTDCHistos(0,0,8,-50.,50.,"FADC-TDC (ns)");    
        createTDCHistos(0,0,9,0.,50.,"TIME-Tvert (ns)");    
        createTDCHistos(0,0,10,-10.,10.,"TIME-Tvert-PATH/c (ns)"); 
        createUVWHistos(0,0,11,700,800,-5,5,"PATH","RESID ");
        createUVWHistos(0,0,12,0,3000,-5,5,"ENERGY","RESID ");
        createUVWHistos(0,0,13,200,270,0,450,"TDC","LEFF ");
        createUVWHistos(0,0,14,190,300,-8,30,"TDC","TEXP ");
        createUVWHistos(0,0,15,200,270,0,450,"TDCC","LEFF ");
        createUVWHistos(0,0,16,190,300,-8,30,"TDCC","TEXP ");
    }

    @Override        
    public void plotHistos() {  
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
    }
    
    public int getDet(int layer) {
        int[] il = {0,0,0,1,1,1,2,2,2}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
     }
    
    public int getLay(int layer) {
        int[] il = {1,2,3,1,2,3,1,2,3}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
    }
    
    @Override
    public void processEvent(DataEvent event) {  
    	
       trigger_sect = getElecTriggerSector(); 
       phase = getTriggerPhase()*4;
       processRaw(event);
       processRec(event);
//       processEB(event);
              
    }
    
    public void processRaw(DataEvent event) {
    	
  	   DataGroup dg0 = this.getDataGroup().getItem(0,0,0);
 	   DataGroup dg1 = this.getDataGroup().getItem(0,0,1);
 	   DataGroup dg2 = this.getDataGroup().getItem(0,0,2);
 	   DataGroup dg3 = this.getDataGroup().getItem(0,0,3);
 	   DataGroup dg4 = this.getDataGroup().getItem(0,0,4);
 	   DataGroup dg8 = this.getDataGroup().getItem(0,0,8);
 	   
       float     tdcd,tdcdc =  0;
        
       tdcs.clear();
        
       int trigger_sect = getElecTriggerSector(); 
       float phase = getTriggerPhase()*4;
        
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
                       dg0.getH2F("tdc_"+layer[getDet(il)]+"_"+view[getLay(il)-1]+"_"+is).fill(tdcd,  ip+0.5);               
                       dg1.getH2F("tdc_"+layer[getDet(il)]+"_"+view[getLay(il)-1]+"_"+is).fill(tdcdc, ip+0.5);                      
                       if (is==trigger_sect) {
                    	       tdcs.getItem(is,il,ip).add(tdcdc);
                           dg2.getH2F("tdc_"+layer[getDet(il)]+"_"+view[getLay(il)-1]+"_"+is).fill(tdcdc, ip+0.5);  
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
               int ped = bank.getShort("ped", i); 
               
               Float[] tdcc; float[] tdc; float tmax = 1000; float tdcm = 1000;
               
               if (tdcs.hasItem(is,il,ip)) {
        	           double a0 = time.getDoubleValue("a0", is, il, ip);
        	           double a2 = time.getDoubleValue("a2", is, il, ip);
                   List<Float> list = new ArrayList<Float>();
                   list = tdcs.getItem(is,il,ip); tdcc=new Float[list.size()]; list.toArray(tdcc);
                   tdc  = new float[list.size()];
                   for (int ii=0; ii<tdcc.length; ii++) {
              	      float tdif = (tdcc[ii]-TOFFSET)-t; 
                      dg8.getH2F("tdc_"+layer[getDet(il)]+"_"+view[getLay(il)-1]+"_"+is).fill(tdif, ip+0.5); 
            	          if (Math.abs(tdif)<30&&tdif<tmax) {tmax = tdif; tdcm = tdcc[ii];}                	    
                   }
                   double tdcmc = tdcm - a0 - a2 / Math.sqrt(adc);
                   dg3.getH2F("tdc_"+layer[getDet(il)]+"_"+view[getLay(il)-1]+"_"+is).fill(tdcm, ip+0.5); 
                   dg4.getH2F("tdc_"+layer[getDet(il)]+"_"+view[getLay(il)-1]+"_"+is).fill(tdcmc,ip+0.5); 
               }
           }
       }    	
    }
    
    public void processRec(DataEvent event) {
    	
  	   DataGroup  dg5 = this.getDataGroup().getItem(0,0,5);
  	   DataGroup  dg6 = this.getDataGroup().getItem(0,0,6);
  	   DataGroup  dg7 = this.getDataGroup().getItem(0,0,7);
  	   DataGroup  dg9 = this.getDataGroup().getItem(0,0,9);
  	   DataGroup dg10 = this.getDataGroup().getItem(0,0,10);
  	   DataGroup dg11 = this.getDataGroup().getItem(0,0,11);
  	   DataGroup dg12 = this.getDataGroup().getItem(0,0,12);
  	   DataGroup dg13 = this.getDataGroup().getItem(0,0,13);
  	   DataGroup dg14 = this.getDataGroup().getItem(0,0,14);
  	   DataGroup dg15 = this.getDataGroup().getItem(0,0,15);
  	   DataGroup dg16 = this.getDataGroup().getItem(0,0,16);
  	   
       event.removeBank("ECAL::hits");        
       event.removeBank("ECAL::peaks");        
       event.removeBank("ECAL::clusters");        
       event.removeBank("ECAL::calib");
        
       engine.processDataEvent(event); 
        
       List<ECStrip>     strips = engine.getStrips();
       List<ECPeak>       peaks = engine.getPeaks();
       List<ECCluster> clusters = engine.getClusters();
       
       if(event.hasBank("ECAL::hits")){
          	DataBank  bank = event.getBank("ECAL::hits");
            for(int loop = 0; loop < bank.rows(); loop++){
               int   is = bank.getByte("sector", loop);
               int   il = bank.getByte("layer", loop); 
               int   ip = bank.getByte("strip", loop);
               float  t = bank.getFloat("time", loop)-phase;
               if (is==trigger_sect) dg5.getH2F("tdc_"+layer[getDet(il)]+"_"+view[getLay(il)-1]+"_"+is).fill(t, ip+0.5); 
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
                int   is = bank1.getByte("sector", loop);
                int   il = bank1.getByte("layer", loop);
                float ener = bank1.getFloat("energy",loop)*1000;
                float  t = bank1.getFloat("time",loop)-phase;
                int   iU = (bank1.getInt("coordU", loop)-4)/8;
                int   iV = (bank1.getInt("coordV", loop)-4)/8;
                int   iW = (bank1.getInt("coordW", loop)-4)/8;
                
                float    tu = (float) clusters.get(loop).getTime(0)-phase;
                float    tv = (float) clusters.get(loop).getTime(1)-phase;
                float    tw = (float) clusters.get(loop).getTime(2)-phase;
                int      iu =         clusters.get(loop).getPeak(0).getMaxStrip();
                int      iv =         clusters.get(loop).getPeak(1).getMaxStrip();
                int      iw =         clusters.get(loop).getPeak(2).getMaxStrip();
                float  tdcu = (float) clusters.get(loop).getPeak(0).getMaxECStrip().getRawTime()-phase;
                float  tdcv = (float) clusters.get(loop).getPeak(1).getMaxECStrip().getRawTime()-phase;
                float  tdcw = (float) clusters.get(loop).getPeak(2).getMaxECStrip().getRawTime()-phase;
                float tdcuc = (float) clusters.get(loop).getPeak(0).getMaxECStrip().getTWCTime()-phase;
                float tdcvc = (float) clusters.get(loop).getPeak(1).getMaxECStrip().getTWCTime()-phase;
                float tdcwc = (float) clusters.get(loop).getPeak(2).getMaxECStrip().getTWCTime()-phase;
                
                float leffu = (float) clusters.get(loop).getPeak(0).getMaxECStrip().getTdist();
                float leffv = (float) clusters.get(loop).getPeak(1).getMaxECStrip().getTdist();
                float leffw = (float) clusters.get(loop).getPeak(2).getMaxECStrip().getTdist();                 
                
                float tudif = tu-Tvertex+phase;
                float tvdif = tv-Tvertex+phase;
                float twdif = tw-Tvertex+phase;
                
                if (is==trigger_sect){
                	    dg6.getH2F("tdc_"+layer[getDet(il)]+"_u_"+is).fill(tu,iu); //Peak times
                	    dg6.getH2F("tdc_"+layer[getDet(il)]+"_v_"+is).fill(tv,iv);
                	    dg6.getH2F("tdc_"+layer[getDet(il)]+"_w_"+is).fill(tw,iw);
                	    dg7.getH2F("tdc_"+layer[getDet(il)]+"_u_"+is).fill(t,iu); //Cluster times
                	    dg7.getH2F("tdc_"+layer[getDet(il)]+"_v_"+is).fill(t,iv);
                	    dg7.getH2F("tdc_"+layer[getDet(il)]+"_w_"+is).fill(t,iw);
                	    dg9.getH2F("tdc_"+layer[getDet(il)]+"_u_"+is).fill(tudif,iu); //TIME-STT
                	    dg9.getH2F("tdc_"+layer[getDet(il)]+"_v_"+is).fill(tvdif,iv);
                	    dg9.getH2F("tdc_"+layer[getDet(il)]+"_w_"+is).fill(twdif,iw);
                  
                    if (pathlist.hasItem(is,il,loop)) {
                    	    DataBank  bankc = event.getBank("REC::Calorimeter");
                    	    DataBank  bankp = event.getBank("REC::Particle");
                    	    int     ip = bankc.getShort("pindex", pathlist.getItem(is,il,loop));
                    	    if (bankp.getInt("pid",ip)==11) {
                        	    float path = bankc.getFloat("path", pathlist.getItem(is,il,loop));
                        	    float tudifp = tudif-path/29.97f;
                        	    float tvdifp = tvdif-path/29.97f;
                        	    float twdifp = twdif-path/29.97f;
                        	    float texpu = path/29.97f + leffu/18.1f + Tvertex - phase - shiftTV[is-1];
                        	    float texpv = path/29.97f + leffv/18.1f + Tvertex - phase - shiftTV[is-1];
                        	    float texpw = path/29.97f + leffw/18.1f + Tvertex - phase - shiftTV[is-1];   
 
                    	    	    dg10.getH2F("tdc_"+layer[getDet(il)]+"_u_"+is).fill(tudifp,iu); //TIME-STT-PATH
                    	    	    dg10.getH2F("tdc_"+layer[getDet(il)]+"_v_"+is).fill(tvdifp,iv);
                    	    	    dg10.getH2F("tdc_"+layer[getDet(il)]+"_w_"+is).fill(twdifp,iw);
                    	    	    dg11.getH2F("tdc_"+layer[getDet(il)]+"_u"+iu+"_s"+is).fill(path,tudifp); 
                    	    	    dg11.getH2F("tdc_"+layer[getDet(il)]+"_v"+iv+"_s"+is).fill(path,tvdifp);
                    	    	    dg11.getH2F("tdc_"+layer[getDet(il)]+"_w"+iw+"_s"+is).fill(path,twdifp);                    	    	    
                    	    	    dg12.getH2F("tdc_"+layer[getDet(il)]+"_u"+iu+"_s"+is).fill(ener,tudifp);
                    	    	    dg12.getH2F("tdc_"+layer[getDet(il)]+"_v"+iv+"_s"+is).fill(ener,tvdifp);
                    	    	    dg12.getH2F("tdc_"+layer[getDet(il)]+"_w"+iw+"_s"+is).fill(ener,twdifp);
                    	    	    
                    	    	    dg13.getH2F("tdc_"+layer[getDet(il)]+"_u"+iu+"_s"+is).fill(tdcu,leffu);
                    	    	    dg13.getH2F("tdc_"+layer[getDet(il)]+"_v"+iv+"_s"+is).fill(tdcv,leffv);
                    	    	    dg13.getH2F("tdc_"+layer[getDet(il)]+"_w"+iw+"_s"+is).fill(tdcw,leffw);
                    	    	    dg14.getH2F("tdc_"+layer[getDet(il)]+"_u"+iu+"_s"+is).fill(tdcu,tdcu-texpu);
                    	    	    dg14.getH2F("tdc_"+layer[getDet(il)]+"_v"+iv+"_s"+is).fill(tdcv,tdcv-texpv);
                    	    	    dg14.getH2F("tdc_"+layer[getDet(il)]+"_w"+iw+"_s"+is).fill(tdcw,tdcw-texpw);                                	                    	    	    
                    	    	    dg15.getH2F("tdc_"+layer[getDet(il)]+"_u"+iu+"_s"+is).fill(tdcuc,leffu);
                    	    	    dg15.getH2F("tdc_"+layer[getDet(il)]+"_v"+iv+"_s"+is).fill(tdcvc,leffv);
                    	    	    dg15.getH2F("tdc_"+layer[getDet(il)]+"_w"+iw+"_s"+is).fill(tdcwc,leffw);
                    	    	    dg16.getH2F("tdc_"+layer[getDet(il)]+"_u"+iu+"_s"+is).fill(tdcuc,tdcuc-texpu);
                    	    	    dg16.getH2F("tdc_"+layer[getDet(il)]+"_v"+iv+"_s"+is).fill(tdcvc,tdcvc-texpv);
                    	    	    dg16.getH2F("tdc_"+layer[getDet(il)]+"_w"+iw+"_s"+is).fill(tdcwc,tdcwc-texpw);                                	                    	    	    
                    	    } 
                    }
                }
            }
       }    	
    }
    
    public void createUVWHistos(int i, int j, int k, int xmin, int xmax, int ymin, int ymax, String xtxt, String ytxt ) {
    	
        DataGroup dg = new DataGroup(3,2);
        H2F h;   
        
        int xbins=50, ybins=50; double sca1=1.0,sca2=1.0;
        
        if (k==12) {sca1=0.5; sca2=0.3;}
        
        for (int is=1; is<7; is++) {
        	   
        	   for (int ip=1; ip<npmts[0]+1; ip++) {
               h = new H2F("tdc_pcal_u"+ip+"_s"+is,"tdc_pcal_u"+ip+"_s"+is,xbins,xmin,xmax,ybins,ymin,ymax);
               h.setTitleX("Sector "+is+" PCAL "+xtxt);
               h.setTitleY(ytxt+"U"+ip);
               dg.addDataSet(h, 1);  
               h = new H2F("tdc_pcal_v"+ip+"_s"+is,"tdc_pcal_v"+ip+"_s"+is,xbins,xmin,xmax,ybins,ymin,ymax);
               h.setTitleX("Sector "+is+" PCAL "+xtxt);
               h.setTitleY(ytxt+"V"+ip);
               dg.addDataSet(h, 1);  
               h = new H2F("tdc_pcal_w"+ip+"_s"+is,"tdc_pcal_w"+ip+"_s"+is,xbins,xmin,xmax,ybins,ymin,ymax);
               h.setTitleX("Sector "+is+" PCAL "+xtxt);
               h.setTitleY(ytxt+"W"+ip); 
               dg.addDataSet(h, 1);  
        	   }
        	   
        	   for (int ip=1; ip<npmts[1]+1; ip++) {
                   h = new H2F("tdc_ecin_u"+ip+"_s"+is,"tdc_ecin_u"+ip+"_s"+is,xbins,xmin*sca1,xmax*sca1,ybins,ymin,ymax);
                   h.setTitleX("Sector "+is+" ECIN "+xtxt);
                   h.setTitleY(ytxt+"U"+ip);
                   dg.addDataSet(h, 1);
                   h = new H2F("tdc_ecin_v"+ip+"_s"+is,"tdc_ecin_v"+ip+"_s"+is,xbins,xmin*sca1,xmax*sca1,ybins,ymin,ymax);
                   h.setTitleX("Sector "+is+" ECIN "+xtxt);
                   h.setTitleY(ytxt+"V"+ip);
                   dg.addDataSet(h, 1);
                   h = new H2F("tdc_ecin_w"+ip+"_s"+is,"tdc_ecin_w"+ip+"_s"+is,xbins,xmin*sca1,xmax*sca1,ybins,ymin,ymax);
                   h.setTitleX("Sector "+is+" ECIN "+xtxt);
                   h.setTitleY(ytxt+"W"+ip);
                   dg.addDataSet(h, 1);
           }
        	   
        	   for (int ip=1; ip<npmts[2]+1; ip++) {
                   h = new H2F("tdc_ecou_u"+ip+"_s"+is,"tdc_ecou_u"+ip+"_s"+is,xbins,xmin*sca2,xmax*sca2,ybins,ymin,ymax);
                   h.setTitleX("Sector "+is+" ECOU "+xtxt);
                   h.setTitleY(ytxt+"U"+ip);
                   dg.addDataSet(h, 1);
                   h = new H2F("tdc_ecou_v"+ip+"_s"+is,"tdc_ecou_v"+ip+"_s"+is,xbins,xmin*sca2,xmax*sca2,ybins,ymin,ymax);
                   h.setTitleX("Sector "+is+" ECOU "+xtxt);
                   h.setTitleY(ytxt+"V"+ip);
                   dg.addDataSet(h, 1);
                   h = new H2F("tdc_ecou_w"+ip+"_s"+is,"tdc_ecou_w"+ip+"_s"+is,xbins,xmin*sca2,xmax*sca2,ybins,ymin,ymax);
                   h.setTitleX("Sector "+is+" ECOU "+xtxt);
                   h.setTitleY(ytxt+"W"+ip);
                   dg.addDataSet(h, 1);
           }
        }
        
        this.getDataGroup().add(dg,i,j,k);       
        
    }
  
    public void createTDCHistos(int i, int j, int k, double tmin, double tmax, String txt) {
    	
        DataGroup dg = new DataGroup(3,2);
        H2F h;  
        
        for (int is=1; is<7; is++) {
            h = new H2F("tdc_pcal_u_"+is,"tdc_pcal_u_"+is,100, tmin, tmax, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt);
            h.setTitleY("U");    
            dg.addDataSet(h, 1);  
            h = new H2F("tdc_pcal_v_"+is,"tdc_pcal_v_"+is,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt);
            h.setTitleY("V");        
            dg.addDataSet(h, 1);            
            h = new H2F("tdc_pcal_w_"+is,"tdc_pcal_w_"+is,100, tmin, tmax, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt);
            h.setTitleY("W");  
            dg.addDataSet(h, 1); 
            
            h = new H2F("tdc_ecin_u_"+is,"tdc_ecin_u_"+is,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt);
            h.setTitleY("U");    
            dg.addDataSet(h, 1);  
            h = new H2F("tdc_ecin_v_"+is,"tdc_ecin_v_"+is,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt);
            h.setTitleY("V");        
            dg.addDataSet(h, 1);            
            h = new H2F("tdc_ecin_w_"+is,"tdc_ecin_w_"+is,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt);
            h.setTitleY("W");  
            dg.addDataSet(h, 1); 
            
            h = new H2F("tdc_ecou_u_"+is,"tdc_ecou_u_"+is,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt);
            h.setTitleY("U");    
            dg.addDataSet(h, 1);  
            h = new H2F("tdc_ecou_v_"+is,"tdc_ecou_v_"+is,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt);
            h.setTitleY("V");        
            dg.addDataSet(h, 1);            
            h = new H2F("tdc_ecou_w_"+is,"tdc_ecou_w_"+is,100, tmin, tmax, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt);
            h.setTitleY("W");  
            dg.addDataSet(h, 1);              
        }    
        
        this.getDataGroup().add(dg,i,j,k);

    } 
    
    public void plotUVWHistos(int index) {
    	
 	   H2F h2; 
	   DataGroup     dg = this.getDataGroup().getItem(0,0,index);	   
       EmbeddedCanvas c = this.getDetectorCanvas().getCanvas(this.getDetectorTabNames().get(index)); 
       
       int s = getActiveSector();       
       int l = getActiveLayer();
       int v = getActiveView();
             
       c.setGridX(false); c.setGridY(false);
       
	   int ipmax=0;
       switch (l) {
       case 0: c.divide(8, 8); ipmax=65; break;
       case 1: c.divide(6, 6); ipmax=37; break;
       case 2: c.divide(6, 6); ipmax=37;
       }
       
       double x1 = dg.getH2F("tdc_"+layer[l]+"_"+view[v]+1+"_s"+s).getXAxis().min();
       double x2 = dg.getH2F("tdc_"+layer[l]+"_"+view[v]+1+"_s"+s).getXAxis().max();
       double y1 = dg.getH2F("tdc_"+layer[l]+"_"+view[v]+1+"_s"+s).getYAxis().min();
       double y2 = dg.getH2F("tdc_"+layer[l]+"_"+view[v]+1+"_s"+s).getYAxis().max();
       
       F1D f1 = new F1D("p0","[a]",x1,x2); f1.setParameter(0,0); f1.setLineColor(1); f1.setLineStyle(1);
       double a0 = 100;
      
       for (int ip=1; ip<ipmax; ip++ ) {
    	      if(index==16) a0 = time.getDoubleValue("a0", s, 3*l+(v+1), ip)+shiftTV[s-1];
    	      c.cd(ip-1); c.getPad().getAxisZ().setLog(getLogZ()); h2 = dg.getH2F("tdc_"+layer[l]+"_"+view[v]+ip+"_s"+s); c.draw(h2);    	      
    	      c.draw(f1,"same"); 
    	      if(index==16&&a0<y2&&a0>y1) {F1D f2 = new F1D("p0","[a]",x1,x2); f2.setLineColor(2); f2.setLineStyle(1);f2.setParameter(0,a0);c.draw(f2,"same");}
       }
       
    }
    
    public void plotTDCHistos(int index) {
    	
    	   H2F h2;
 	   DataGroup     dg = this.getDataGroup().getItem(0,0,index);	   
       EmbeddedCanvas c = this.getDetectorCanvas().getCanvas(this.getDetectorTabNames().get(index));
       
       c.setGridX(false); c.setGridY(false);
       c.divide(3, 3);
       
       int s = getActiveSector();
       
       c.cd(0); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_pcal_u_"+s); c.draw(h2);
       c.cd(1); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_pcal_v_"+s); c.draw(h2);
       c.cd(2); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_pcal_w_"+s); c.draw(h2);
       c.cd(3); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_ecin_u_"+s); c.draw(h2);
       c.cd(4); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_ecin_v_"+s); c.draw(h2);
       c.cd(5); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_ecin_w_"+s); c.draw(h2);
       c.cd(6); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_ecou_u_"+s); c.draw(h2);
       c.cd(7); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_ecou_v_"+s); c.draw(h2);
       c.cd(8); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("tdc_ecou_w_"+s); c.draw(h2);	
    	
    }
    
    @Override
    public void timerUpdate() {
    	
    }
     
}
