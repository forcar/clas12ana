package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */

public class CTOFmonitor  extends DetectorMonitor {        
    
    public CTOFmonitor(String name) {
        super(name);
        
        this.setDetectorTabNames("ADC Occupancies", "TDC Occupancies");
        this.init(false);
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",96,0.5,96.5);
        summary.setTitleX("PMT");
        summary.setTitleY("CTOF hits");
        summary.setTitle("CTOF");
        summary.setFillColor(38);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        
        H1F occADCL = new H1F("occADCL", "occADCL", 48, 0.5, 48.5);
        occADCL.setTitleX("PMT Upstream");
        occADCL.setTitleY("Counts");
        occADCL.setFillColor(38);
        H1F occADCR = new H1F("occADCR", "occADCR", 48, 0.5, 48.5);
        occADCR.setTitleX("PMT Downstream");
        occADCR.setTitleY("Counts");
        occADCR.setFillColor(38);
        H2F adcL = new H2F("adcL", "adcL", 50, 0, 5000, 48, 0.5, 48.5);
        adcL.setTitleX("ADC Upstream - amplitude");
        adcL.setTitleY("PMT Upstream");
        H2F adcR = new H2F("adcR", "adcR", 50, 0, 5000, 48, 0.5, 48.5);
        adcR.setTitleX("ADC Downstream - amplitude");
        adcR.setTitleY("PMT Downstream");   
        H2F fadcL_time = new H2F("fadcL_time", "fadcL_time", 80, 0, 400, 48, 0.5, 48.5);
        fadcL_time.setTitleX("FADC Upstream - timing");
        fadcL_time.setTitleY("PMT Upstream");
        H2F fadcR_time = new H2F("fadcR_time", "fadcR_time", 80, 0, 400, 48, 0.5, 48.5);
        fadcR_time.setTitleX("FADC Downstream - timing");
        fadcR_time.setTitleY("PMT Downstream");  
        H1F occTDCL = new H1F("occTDCL", "occTDCL", 48, 0.5, 48.5);
        occTDCL.setTitleX("PMT Upstream");
        occTDCL.setTitleY("Counts");
        occTDCL.setFillColor(38);
        H1F occTDCR = new H1F("occTDCR", "occTDCR", 48, 0.5, 48.5);
        occTDCR.setTitleX("PMT Downstream");
        occTDCR.setTitleY("Counts");
        occTDCR.setFillColor(38);
        H2F tdcL = new H2F("tdcL", "tdcL", 50, 0, 50000, 48, 0.5, 48.5);
        tdcL.setTitleX("TDC Upstream - amplitude");
        tdcL.setTitleY("PMT Upstream");
        H2F tdcR = new H2F("tdcR", "tdcR", 50, 0, 50000, 48, 0.5, 48.5);
        tdcR.setTitleX("TDC Downstream - amplitude");
        tdcR.setTitleY("PMT Downstream"); 
        
        DataGroup dg = new DataGroup(2,5);
        dg.addDataSet(occADCL, 0);
        dg.addDataSet(occADCR, 1);
        dg.addDataSet(adcL, 2);
        dg.addDataSet(adcR, 3);
        dg.addDataSet(fadcL_time, 4);
        dg.addDataSet(fadcR_time, 5);
        dg.addDataSet(occTDCL, 6);
        dg.addDataSet(occTDCR, 7);
        dg.addDataSet(tdcL, 8);
        dg.addDataSet(tdcR, 9);
        this.getDataGroup().add(dg,0,0,0);
    }
        
    @Override
    public void plotHistos() {        
        // plotting histos
        this.getDetectorCanvas().getCanvas("ADC Occupancies").divide(2, 3);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").setGridX(false);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").setGridY(false);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").divide(2, 2);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").setGridX(false);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").setGridY(false);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").cd(0);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH1F("occADCL"));
        this.getDetectorCanvas().getCanvas("ADC Occupancies").cd(1);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH1F("occADCR"));
        this.getDetectorCanvas().getCanvas("ADC Occupancies").cd(2);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").getPad(2).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("adcL"));
        this.getDetectorCanvas().getCanvas("ADC Occupancies").cd(3);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").getPad(3).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("adcR"));
        this.getDetectorCanvas().getCanvas("ADC Occupancies").cd(4);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").getPad(4).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadcL_time"));
        this.getDetectorCanvas().getCanvas("ADC Occupancies").cd(5);
        this.getDetectorCanvas().getCanvas("ADC Occupancies").getPad(5).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadcR_time"));
        this.getDetectorCanvas().getCanvas("ADC Occupancies").update();
        this.getDetectorCanvas().getCanvas("TDC Occupancies").cd(0);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH1F("occTDCL"));
        this.getDetectorCanvas().getCanvas("TDC Occupancies").cd(1);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH1F("occTDCR"));
        this.getDetectorCanvas().getCanvas("TDC Occupancies").cd(2);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").getPad(2).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdcL"));
        this.getDetectorCanvas().getCanvas("TDC Occupancies").cd(3);
        this.getDetectorCanvas().getCanvas("TDC Occupancies").getPad(3).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdcR"));
        this.getDetectorCanvas().getCanvas("TDC Occupancies").update();
    }

    @Override
    public void processEvent(DataEvent event) {
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[3] && super.eventResetTime_current[3] > 0){
            resetEventListener();
        }
        
		//if (!testTriggerMask()) return;
                
        if(event.hasBank("CTOF::adc")==true){
	    DataBank bank = event.getBank("CTOF::adc");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector  = bank.getByte("sector", loop);
                int layer   = bank.getByte("layer", loop);
                int comp    = bank.getShort("component", loop);
                int order   = bank.getByte("order", loop);
                int adc     = bank.getInt("ADC", loop);
                float time  = bank.getFloat("time", loop);
//                System.out.println("ROW " + loop + " SECTOR = " + sector + " LAYER = " + layer + " COMPONENT = " + comp + " ORDER + " + order +
//                      " ADC = " + adc + " TIME = " + time); 
                if(adc>0) {
                    if(order==0) {
                        this.getDataGroup().getItem(0,0,0).getH1F("occADCL").fill(comp*1.0);
                        this.getDataGroup().getItem(0,0,0).getH2F("adcL").fill(adc*1.0,comp*1.0);
                        if(time > 1) this.getDataGroup().getItem(0,0,0).getH2F("fadcL_time").fill(time*1.0,comp*1.0);
                    }
                    else if(order==1) {
                        this.getDataGroup().getItem(0,0,0).getH1F("occADCR").fill(comp*1.0);
                        this.getDataGroup().getItem(0,0,0).getH2F("adcR").fill(adc*1.0,comp*1.0);
                        if(time > 1) this.getDataGroup().getItem(0,0,0).getH2F("fadcR_time").fill(time*1.0,comp*1.0);
                    }
                    
                    this.getDetectorSummary().getH1F("summary").fill((order*48+comp)*1.0); 
                }
                
	    }
    	}
        if(event.hasBank("CTOF::tdc")==true){
            DataBank  bank = event.getBank("CTOF::tdc");
            int rows = bank.rows();
            for(int i = 0; i < rows; i++){
                int    sector = bank.getByte("sector",i);
                int     layer = bank.getByte("layer",i);
                int      comp = bank.getShort("component",i);
                int       tdc = bank.getInt("TDC",i);
                int     order = bank.getByte("order",i); // order specifies left-right for ADC
//                           System.out.println("ROW " + i + " SECTOR = " + sector
//                                 + " LAYER = " + layer + " PADDLE = "
//                                 + paddle + " TDC = " + TDC);    
                if(tdc>0) {
                    if(order==2) {
                        this.getDataGroup().getItem(0,0,0).getH1F("occTDCL").fill(comp*1.0);
                        this.getDataGroup().getItem(0,0,0).getH2F("tdcL").fill(tdc*1.0,comp*1.0);
                    }
                    else if(order==3) {
                        this.getDataGroup().getItem(0,0,0).getH1F("occTDCR").fill(comp*1.0);
                        this.getDataGroup().getItem(0,0,0).getH2F("tdcR").fill(tdc*1.0,comp*1.0);
                    }
                }
            }
        }        
    }

    @Override
    public void timerUpdate() {

    }


}
