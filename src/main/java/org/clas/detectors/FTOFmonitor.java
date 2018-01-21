package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;


public class FTOFmonitor  extends DetectorMonitor {

    private final int[] npaddles = new int[]{23,62,5};

    FTOFHits          ftofHits[] = new FTOFHits[3];    
    
    public FTOFmonitor(String name) {
        super(name);
        
        this.setDetectorTabNames("ADC Occupancies", "TDC Occupancies","ADC Histograms", "FADC timing", "TDC Histograms","GMEAN");
        this.useSectorButtons(true);
        this.init(false);   // set to true for picture on left side
        ftofHits[0] = new FTOFHits("PANEL1A");
        ftofHits[1] = new FTOFHits("PANEL1B");
        ftofHits[2] = new FTOFHits("PANEL2");
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").divide(2, 3);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").setGridX(false);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").setGridY(false);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").divide(2, 3);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").setGridX(false);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").setGridY(false);
        this.getDetectorCanvas().getCanvas("ADC Histograms").divide(2, 3);
        this.getDetectorCanvas().getCanvas("ADC Histograms").setGridX(false);
        this.getDetectorCanvas().getCanvas("ADC Histograms").setGridY(false);
        this.getDetectorCanvas().getCanvas("FADC timing").divide(2, 3);
        this.getDetectorCanvas().getCanvas("FADC timing").setGridX(false);
        this.getDetectorCanvas().getCanvas("FADC timing").setGridY(false);
        this.getDetectorCanvas().getCanvas("TDC Histograms").divide(2, 3);
        this.getDetectorCanvas().getCanvas("TDC Histograms").setGridX(false);
        this.getDetectorCanvas().getCanvas("TDC Histograms").setGridY(false);
        this.getDetectorCanvas().getCanvas("GMEAN").divide(2, 3);
        this.getDetectorCanvas().getCanvas("GMEAN").setGridX(false);
        this.getDetectorCanvas().getCanvas("GMEAN").setGridY(false);
        
        String[] stacks = new String[]{"P1A","P1B","P2"};
        String[] views = new String[]{"Left","Right"};   
        
        H2F sumStackp1 = new H2F("sum_p1","sum_p1",12,0.5,6.5, 2,-0.5,1.5);
        sumStackp1.setTitleX("sector P1A and P1B");
        sumStackp1.setTitleY("left / right PMT");
        sumStackp1.setTitle("FTOF panel 1a and 1b");
        H2F sumStackp2 = new H2F("sum_p2","sum_p2",6,0.5,6.5, 2,-0.5,1.5);
        sumStackp2.setTitleX("sector");
        sumStackp2.setTitleY("left / right PMT");
        sumStackp2.setTitle("FTOF panel 2");
            
        DataGroup sum = new DataGroup(2,1); 
        sum.addDataSet(sumStackp1, 0);
        sum.addDataSet(sumStackp2, 1);
        this.setDetectorSummary(sum);
        
        for(int lay=0; lay < 3; lay++) {
            DataGroup dg = new DataGroup(2,2);
        	for(int ord=0; ord < 2; ord++) {
            H2F occADC = new H2F("occADC"+lay+ord, "lay/ord " + lay + ord + " Occupancy", 6, 0.5, 6.5, npaddles[lay], 1, npaddles[lay]+1);
            occADC.setTitle(stacks[lay]+" "+views[ord]+" PMTS");
            occADC.setTitleY("paddle");
            occADC.setTitleX("sector");
            H2F occTDC = new H2F("occTDC"+lay+ord, "lay/ord " + lay + ord + " Occupancy", 6, 0.5, 6.5, npaddles[lay], 1, npaddles[lay]+1);
            occTDC.setTitle(stacks[lay]+" "+views[ord]+" PMTS");
            occTDC.setTitleY("paddle");
            occTDC.setTitleX("sector");           
            dg.addDataSet(occADC, 0);
            dg.addDataSet(occTDC, 0);
            this.getDataGroup().add(dg,0,lay,0);
        }
        }
        
        for(int sec=1; sec < 7; sec++) {
        for(int lay=0; lay < 3; lay++) {
            DataGroup dg = new DataGroup(2,3);
        	for(int ord=0; ord < 2; ord++) {
            H2F datADC = new H2F("datADC"+sec+lay+ord, "sec/lay/ord "+sec+lay+ord+" ADC", 100, 0., 4000., npaddles[lay], 1, npaddles[lay]+1);
            datADC.setTitleY(stacks[lay] + " " + views[ord] + " PMTS");
            datADC.setTitleX("ADC Channel");
            datADC.setTitle("Sector "+sec);
            H2F timeFADC = new H2F("timeFADC"+sec+lay+ord, "sec/lay/ord "+sec+lay+ord+" FADC", 80, 0., 400., npaddles[lay], 1, npaddles[lay]+1);
            timeFADC.setTitleY(stacks[lay] + " " + views[ord] + " PMTS");
            timeFADC.setTitleX("FADC timing");
            timeFADC.setTitle("Sector "+sec);
            H2F datTDC = new H2F("datTDC"+sec+lay+ord, "sec/lay/ord "+sec+lay+ord+" TDC", 100, 450., 850., npaddles[lay], 1, npaddles[lay]+1);
            datTDC.setTitleY(stacks[lay] + " " + views[ord] + " PMTS");
            datTDC.setTitleX("TDC Channel");
            datTDC.setTitle("Sector "+sec);
            dg.addDataSet(datADC, 0);
            dg.addDataSet(timeFADC, 1);
            dg.addDataSet(datTDC, 2);
        }
            H2F GMEAN = new H2F("GMEAN"+sec+lay, "sec/lay"+sec+lay+" GMEAN", 100, 0., 6000.,this.npaddles[lay], 1, npaddles[lay]+1 );
            GMEAN.setTitleY(stacks[lay] +" PMTS");
            GMEAN.setTitleX("GMEAN");
            GMEAN.setTitle("Sector "+sec);
            dg.addDataSet(GMEAN, 3);
            H2F TDIF = new H2F("TDIF"+sec+lay, "sec/lay"+sec+lay+" TDIF", 100, -40., 40.,this.npaddles[lay], 1, npaddles[lay]+1 );
            TDIF.setTitleY(stacks[lay] +" PMTS");
            TDIF.setTitleX("TLeft-TRight");
            TDIF.setTitle("Sector "+sec);
            dg.addDataSet(TDIF, 4);
            this.getDataGroup().add(dg,sec,lay,0);
        }
        }
    }

    @Override
    
    public void drawDetector() {
        
        double FTOFSize = 500.0;
        int[]     widths   = new int[]{6,15,25};
        int[]     lengths  = new int[]{6,15,25};

        String[]  names    = new String[]{"FTOF 1A","FTOF 1B","FTOF 2"};
        for(int sector = 1; sector <= 6; sector++){
            double rotation = Math.toRadians((sector-1)*(360.0/6)+90.0);
            for(int layer = 1; layer <=3; layer++){
                int width  = widths[layer-1];
                int length = lengths[layer-1];
                for(int paddle = 1; paddle <= npaddles[layer-1]; paddle++){
                    DetectorShape2D shape = new DetectorShape2D();
                    shape.getDescriptor().setType(DetectorType.FTOF);
                    shape.getDescriptor().setSectorLayerComponent(sector, layer, paddle);
                    shape.createBarXY(20 + length*paddle, width);
                    shape.getShapePath().translateXYZ(0.0, 40 + width*paddle , 0.0);
                    shape.getShapePath().rotateZ(rotation);
                    this.getDetectorView().getView().addShape(names[layer-1], shape);
                }
            }
        }
        this.getDetectorView().setName("FTOF");
        this.getDetectorView().updateBox();
                
    }
            
    @Override
    public void plotHistos() {
    	
        for(int lay=0; lay<3; lay++) {
        	for(int ord=0; ord<2; ord++) {
            this.getDetectorCanvas().getCanvas("ADC Occupancies").cd(lay*2+ord);
            this.getDetectorCanvas().getCanvas("ADC Occupancies").getPad(lay*2+ord).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("ADC Occupancies").draw(this.getDataGroup().getItem(0,lay,0).getH2F("occADC"+lay+ord));
            this.getDetectorCanvas().getCanvas("TDC Occupancies").cd(lay*2+ord);
            this.getDetectorCanvas().getCanvas("TDC Occupancies").getPad(lay*2+ord).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("TDC Occupancies").draw(this.getDataGroup().getItem(0,lay,0).getH2F("occTDC"+lay+ord));
        }
        }
        
        for(int sec=1; sec<7; sec++) {
            if(getActiveSector()==sec) {
               for(int lay=0; lay < 3; lay++) {
            	   for(int ord=0; ord < 2; ord++) {
                   this.getDetectorCanvas().getCanvas("ADC Histograms").cd(lay*2+ord);
                   this.getDetectorCanvas().getCanvas("ADC Histograms").getPad(lay*2+ord).getAxisZ().setLog(getLogZ());
                   this.getDetectorCanvas().getCanvas("ADC Histograms").draw(this.getDataGroup().getItem(sec,lay,0).getH2F("datADC"+sec+lay+ord));
                   this.getDetectorCanvas().getCanvas("FADC timing").cd(lay*2+ord);
                   this.getDetectorCanvas().getCanvas("FADC timing").getPad(lay*2+ord).getAxisZ().setLog(getLogZ());
                   this.getDetectorCanvas().getCanvas("FADC timing").draw(this.getDataGroup().getItem(sec,lay,0).getH2F("timeFADC"+sec+lay+ord));
                   this.getDetectorCanvas().getCanvas("TDC Histograms").cd(lay*2+ord);
                   this.getDetectorCanvas().getCanvas("TDC Histograms").getPad(lay*2+ord).getAxisZ().setLog(getLogZ());
                   this.getDetectorCanvas().getCanvas("TDC Histograms").draw(this.getDataGroup().getItem(sec,lay,0).getH2F("datTDC"+sec+lay+ord));
               }
      	           this.getDetectorCanvas().getCanvas("GMEAN").cd(lay*2+0);
                   this.getDetectorCanvas().getCanvas("GMEAN").getPad(lay*2+0).getAxisZ().setLog(getLogZ());
                   this.getDetectorCanvas().getCanvas("GMEAN").draw(this.getDataGroup().getItem(sec,lay,0).getH2F("GMEAN"+sec+lay));
      	           this.getDetectorCanvas().getCanvas("GMEAN").cd(lay*2+1);
                   this.getDetectorCanvas().getCanvas("GMEAN").getPad(lay*2+1).getAxisZ().setLog(getLogZ());
                   this.getDetectorCanvas().getCanvas("GMEAN").draw(this.getDataGroup().getItem(sec,lay,0).getH2F("TDIF"+sec+lay));
               }
    	       }
        }
        
        this.getDetectorCanvas().getCanvas("ADC Occupancies").update();
        this.getDetectorCanvas().getCanvas("TDC Occupancies").update();
        this.getDetectorCanvas().getCanvas("ADC Histograms").update();
        this.getDetectorCanvas().getCanvas("TDC Histograms").update();
        this.getDetectorCanvas().getCanvas("GMEAN").update();
//      this.getDetectorView().getView().repaint();
//      this.getDetectorView().update();

    }           

    @Override
    public void processEvent(DataEvent event) {
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[9] && super.eventResetTime_current[9] > 0){
            resetEventListener();
        }
        
		//if (!testTriggerMask()) return;
	    
	    clear(0); clear(1); clear(2);
    	
        if(event.hasBank("FTOF::adc")==true){
            DataBank  bank = event.getBank("FTOF::adc");
            int rows = bank.rows();
            for(int i = 0; i < rows; i++){
                int    sector = bank.getByte("sector",i);
                int     layer = bank.getByte("layer",i);
                int    paddle = bank.getShort("component",i);
                int       ADC = bank.getInt("ADC",i);
                float    time = bank.getFloat("time",i);
                int     order = bank.getByte("order",i);  
                
                int lay=layer-1; int ord=order-0;
                if(ADC>0 && isGoodECALTrigger(sector)) {
                	  this.getDataGroup().getItem(0,lay,0).getH2F("occADC"+lay+ord).fill(sector*1.0,paddle*1.0);
                	  this.getDataGroup().getItem(sector,lay,0).getH2F("datADC"+sector+lay+ord).fill(ADC,paddle*1.0);
                	  if(time > 1) this.getDataGroup().getItem(sector,lay,0).getH2F("timeFADC"+sector+lay+ord).fill(time,paddle*1.0);
                	  if(layer == 1) this.getDetectorSummary().getH2F("sum_p1").fill(sector-0.25, order*1.0);
                          if(layer == 2) this.getDetectorSummary().getH2F("sum_p1").fill(sector+0.25, order*1.0); 
                          this.getDetectorSummary().getH2F("sum_p2").fill(sector*1.0, order*1.0); 
                          storeADCHits(lay,sector-1,ord,paddle,ADC,time);
                }
            }
        }
        if(event.hasBank("FTOF::tdc")==true){
            DataBank  bank = event.getBank("FTOF::tdc");
            int rows = bank.rows();
            for(int i = 0; i < rows; i++){
                int    sector = bank.getByte("sector",i);
                int     layer = bank.getByte("layer",i);
                int    paddle = bank.getShort("component",i);
                int       TDC = bank.getInt("TDC",i);
                int     order = bank.getByte("order",i); 
                
                int lay=layer-1; int ord=order-2;
                if(TDC>0 && isGoodECALTrigger(sector)) {
                	   this.getDataGroup().getItem(0,lay,0).getH2F("occTDC"+lay+ord).fill(sector*1.0,paddle*1.0);
                	   this.getDataGroup().getItem(sector,lay,0).getH2F("datTDC"+sector+lay+ord).fill(TDC*0.02345-triggerPhase*4,paddle*1.0);
                   storeTDCHits(lay,sector-1,ord,paddle,(float)(TDC*0.02345-triggerPhase*4));
                }
            }
        }
        
        for(int sec=1; sec<7; sec++) {
          	if (isGoodECALTrigger(sec)) {
          	    for(int il=0; il<3; il++) {
        		        getGM(il,sec-1,this.getDataGroup().getItem(sec,il,0).getH2F("GMEAN"+sec+il));
        		        getTD(il,sec-1,this.getDataGroup().getItem(sec,il,0).getH2F("TDIF"+sec+il));
                 }
          	}
        }
    }
    
    public class FTOFHits {
    	    public String  detName = null;
        int        nha[][] = new    int[6][2];
        int        nht[][] = new    int[6][2];
        int    strra[][][] = new    int[6][2][62]; 
        int    strrt[][][] = new    int[6][2][62]; 
        int     adcr[][][] = new    int[6][2][62];      
        float   tdcr[][][] = new  float[6][2][62]; 
        float     tf[][][] = new  float[6][2][62]; 
        float     ph[][][] = new  float[6][2][62]; 
        
        public FTOFHits(String det) {
        	   detName = det;
        }
    }  
    
    public void clear(int idet) {
        
        for (int is=0 ; is<6 ; is++) {
            for (int il=0 ; il<2 ; il++) {
                ftofHits[idet].nha[is][il] = 0;
                ftofHits[idet].nht[is][il] = 0;
                for (int ip=0 ; ip<62 ; ip++) {
                    ftofHits[idet].strra[is][il][ip] = 0;
                    ftofHits[idet].strrt[is][il][ip] = 0;
                    ftofHits[idet].adcr[is][il][ip]  = 0;
                    ftofHits[idet].tdcr[is][il][ip]  = 0;
                    ftofHits[idet].tf[is][il][ip]    = 0;
                    ftofHits[idet].ph[is][il][ip]    = 0;
                }
            }
        } 
    }
        
    public void storeTDCHits(int idet, int is, int il, int ip, float tdc) {

        if(tdc>450&&tdc<850){
              ftofHits[idet].nht[is][il]++; int inh = ftofHits[idet].nht[is][il];
              if (inh>npaddles[idet]) inh=npaddles[idet];            
              ftofHits[idet].tdcr[is][il][inh-1]  = tdc;
              ftofHits[idet].strrt[is][il][inh-1] = ip;                                 
        }
    }
    
    public void storeADCHits(int idet, int is, int il, int ip, int adc, float time) {
        
        if(adc>50){
              ftofHits[idet].nha[is][il]++; int inh = ftofHits[idet].nha[is][il];
              if (inh>npaddles[idet]) inh=npaddles[idet];
              ftofHits[idet].adcr[is][il][inh-1]  = adc;
              ftofHits[idet].tf[is][il][inh-1]    = time;
              ftofHits[idet].strra[is][il][inh-1] = ip;
        } 
    }
    
    public void getGM(int idet, int is, H2F h2) {
        int     iL = ftofHits[idet].nha[is][0];
        int     iR = ftofHits[idet].nha[is][1];
        int    ipL = ftofHits[idet].strra[is][0][0];
        int    ipR = ftofHits[idet].strra[is][1][0];
        double adL = ftofHits[idet].adcr[is][0][0];
        double adR = ftofHits[idet].adcr[is][1][0];
    	    if ((iL==1&&iR==1)&&(ipL==ipR)) h2.fill(Math.sqrt(adL*adR), ipL*1.0);
    }
    
    public void getTD(int idet, int is, H2F h2) {
        int     iL = ftofHits[idet].nht[is][0];
        int     iR = ftofHits[idet].nht[is][1];
        int    ipL = ftofHits[idet].strrt[is][0][0];
        int    ipR = ftofHits[idet].strrt[is][1][0];
        double tdL = ftofHits[idet].tdcr[is][0][0];
        double tdR = ftofHits[idet].tdcr[is][1][0];
    	   if ((iL==1&&iR==1)&&(ipL==ipR)) h2.fill(tdL-tdR, ipL*1.0);
    }    
    
    @Override
    public void timerUpdate() {

    }


}
