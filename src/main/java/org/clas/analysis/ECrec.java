package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.service.ec.ECEngine;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;

public class ECrec extends DetectorMonitor {

    ECEngine engine = new ECEngine();
    IndexedList<List<Float>> tdcs = new IndexedList<List<Float>>(3);
    IndexedTable time;
    String[] layer = new String[]{"pcal","ecin","ecou"};
    String[] view  = new String[]{"u","v","w"};   
    int trigger_sect = 0;
    float phase = 0;
    
//  static float TOFFSET = 436; 
    static float TOFFSET = 125; 
    
    public ECrec(String name) {
        super(name);
        this.setDetectorTabNames("Raw TDC","PhaseCorr TDC","Triggered TDC","Matched TDC","Calib TDC","Hit Time","Peak Time","Cluster Time","TDIF","STTDIF");
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init(false);
        engine.init();       
        engine.setVariation("default");
        time = engine.getConstantsManager().getConstants(3050, "/calibration/ec/timing");
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
        createTDCHistos(0,0,9,0.,50.,"TDC-START (ns)");    
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
 	   
 	   float      tps =  (float) 0.02345;
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
    	
  	   DataGroup dg5 = this.getDataGroup().getItem(0,0,5);
  	   DataGroup dg6 = this.getDataGroup().getItem(0,0,6);
  	   DataGroup dg7 = this.getDataGroup().getItem(0,0,7);
  	   DataGroup dg9 = this.getDataGroup().getItem(0,0,9);
  	   
       event.removeBank("ECAL::hits");        
       event.removeBank("ECAL::peaks");        
       event.removeBank("ECAL::clusters");        
       event.removeBank("ECAL::calib");
        
       engine.processDataEvent(event); 
        
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
      
       float stt = 0;
       
       if(event.hasBank("REC::Event")){
    	      DataBank bank = event.getBank("REC::Event");
    	      stt = bank.getFloat("STTime", 0);
       }
       		        
       if(event.hasBank("ECAL::clusters")){
            DataBank  bank = event.getBank("ECAL::clusters");
        	    DataBank  bank2 = event.getBank("ECAL::peaks");
            for(int loop = 0; loop < bank.rows(); loop++){
                int   is = bank.getByte("sector", loop);
                int   il = bank.getByte("layer", loop);
                float  e = bank.getFloat("energy",loop)*1000;
                float  t = bank.getFloat("time",loop)-phase;
                int   iU = (bank.getInt("coordU", loop)-4)/8;
                int   iV = (bank.getInt("coordV", loop)-4)/8;
                int   iW = (bank.getInt("coordW", loop)-4)/8;
                int  idU = bank.getByte("idU", loop)-1;
                int  idV = bank.getByte("idV", loop)-1;
                int  idW = bank.getByte("idW", loop)-1;
                float tu = bank2.getFloat("time", idU)-phase;
                float tv = bank2.getFloat("time", idV)-phase;
                float tw = bank2.getFloat("time", idW)-phase;
                float tdif = t-stt+phase;
                if (is==trigger_sect) {
                  dg6.getH2F("tdc_"+layer[getDet(il)]+"_"+"u_"+is).fill(tu,iU+1.5); 
                  dg6.getH2F("tdc_"+layer[getDet(il)]+"_"+"v_"+is).fill(tv,iV+1.5); 
                  dg6.getH2F("tdc_"+layer[getDet(il)]+"_"+"w_"+is).fill(tw,iW+1.5); 
                  dg7.getH2F("tdc_"+layer[getDet(il)]+"_"+"u_"+is).fill(t,iU+1.5); 
                  dg7.getH2F("tdc_"+layer[getDet(il)]+"_"+"v_"+is).fill(t,iV+1.5); 
                  dg7.getH2F("tdc_"+layer[getDet(il)]+"_"+"w_"+is).fill(t,iW+1.5); 
                  dg9.getH2F("tdc_"+layer[getDet(il)]+"_"+"u_"+is).fill(tdif,iU+1.5); 
                  dg9.getH2F("tdc_"+layer[getDet(il)]+"_"+"v_"+is).fill(tdif,iV+1.5); 
                  dg9.getH2F("tdc_"+layer[getDet(il)]+"_"+"w_"+is).fill(tdif,iW+1.5); 
                }
            }
       }    	
    }
    
    void processEB(DataEvent event) {
    	
    	    int setdetsector = 1;
    	
		if(event.hasBank("REC::Particle") 
		&& event.hasBank("REC::Calorimeter")  
		&& event.hasBank("REC::Scintillator")  
		&& event.hasBank("REC::Event")  
		&& event.hasBank("ECAL::clusters") 
		&& event.hasBank("ECAL::adc")  
		&& event.hasBank("ECAL::tdc")  
		&& event.hasBank("RUN::rf")) {
			HipoDataBank recpart = (HipoDataBank) event.getBank("REC::Particle");
			float p = 0;
			int npart  = 0;
			int nelec  = 0;
			for(int parti = 0; parti < recpart.rows(); parti++){
				if(recpart.getInt("pid", parti) == 11){
					nelec++;
					HipoDataBank reccal = (HipoDataBank) event.getBank("REC::Calorimeter");
					for(int i = 0; i < reccal.rows(); i++){
						short pindex  = reccal.getShort("pindex", i);
						byte detector = reccal.getByte("detector", i);
						byte sector   = reccal.getByte("sector", i);
						if(pindex != parti &&  detector == 7) npart++;
	
					}					   
                }
            }
            System.out.println("elec,part = "+nelec +" "+npart );
        }
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
    
}
