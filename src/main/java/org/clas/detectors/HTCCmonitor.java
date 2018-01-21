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

public class HTCCmonitor  extends DetectorMonitor {
        
    
    public HTCCmonitor(String name) {
        super(name);
        
        this.setDetectorTabNames("Occupancies", "ADC spectra", "FADC timing spectra", "TDC spectra");
        this.init(false);
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        this.getDetectorCanvas().getCanvas("Occupancies").divide(1, 2);
        this.getDetectorCanvas().getCanvas("Occupancies").setGridX(false);
        this.getDetectorCanvas().getCanvas("Occupancies").setGridY(false);
        this.getDetectorCanvas().getCanvas("ADC spectra").divide(3, 2);
        this.getDetectorCanvas().getCanvas("ADC spectra").setGridX(false);
        this.getDetectorCanvas().getCanvas("ADC spectra").setGridY(false);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").divide(3, 2);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").setGridX(false);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").setGridY(false);
        this.getDetectorCanvas().getCanvas("TDC spectra").divide(3, 2);
        this.getDetectorCanvas().getCanvas("TDC spectra").setGridX(false);
        this.getDetectorCanvas().getCanvas("TDC spectra").setGridY(false);
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("HTCC hits");
        summary.setTitle("HTTC");
        summary.setFillColor(36);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        H2F occADC = new H2F("occADC", "occADC", 6, 0.5, 6.5, 8, 0.5, 8.5);
        occADC.setTitleY("ring-PMT");
        occADC.setTitleX("sector");
        occADC.setTitle("ADC Occupancy");
        H2F occTDC = new H2F("occTDC", "occTDC", 6, 0.5, 6.5, 8, 0.5, 8.5);
        occTDC.setTitleY("ring-PMT");
        occTDC.setTitleX("sector");
        occTDC.setTitle("TDC Occupancy");
        
        H2F adc1 = new H2F("adc_s1", "adc_s1", 150, 0, 15000, 8, 0.5, 8.5);
        adc1.setTitleX("ADC - value");
        adc1.setTitleY("PMT");
        adc1.setTitle("sector 1");
        H2F adc2 = new H2F("adc_s2", "adc_s2", 150, 0, 15000, 8, 0.5, 8.5);
        adc2.setTitleX("ADC - value");
        adc2.setTitleY("PMT");
        adc2.setTitle("sector 2");
        H2F adc3 = new H2F("adc_s3", "adc_s3", 150, 0, 15000, 8, 0.5, 8.5);
        adc3.setTitleX("ADC - value");
        adc3.setTitleY("PMT");
        adc3.setTitle("sector 3");
        H2F adc4 = new H2F("adc_s4", "adc_s4", 150, 0, 15000, 8, 0.5, 8.5);
        adc4.setTitleX("ADC - value");
        adc4.setTitleY("PMT");
        adc4.setTitle("sector 4");
        H2F adc5 = new H2F("adc_s5", "adc_s5", 150, 0, 15000, 8, 0.5, 8.5);
        adc5.setTitleX("ADC - value");
        adc5.setTitleY("PMT");
        adc5.setTitle("sector 5");
        H2F adc6 = new H2F("adc_s6", "adc_s6", 150, 0, 15000, 8, 0.5, 8.5);
        adc6.setTitleX("ADC - value");
        adc6.setTitleY("PMT");
        adc6.setTitle("sector 6");
        
        H2F fadc_time1 = new H2F("fadc_time_s1", "fadc_time_s1", 80, 0, 400, 8, 0.5, 8.5);
        fadc_time1.setTitleX("FADC timing");
        fadc_time1.setTitleY("PMT");
        fadc_time1.setTitle("sector 1");
        H2F fadc_time2 = new H2F("fadc_time_s2", "fadc_time_s2", 80, 0, 400, 8, 0.5, 8.5);
        fadc_time2.setTitleX("FADC timing");
        fadc_time2.setTitleY("PMT");
        fadc_time2.setTitle("sector 2");
        H2F fadc_time3 = new H2F("fadc_time_s3", "fadc_time_s3", 80, 0, 400, 8, 0.5, 8.5);
        fadc_time3.setTitleX("FADC timing");
        fadc_time3.setTitleY("PMT");
        fadc_time3.setTitle("sector 3");
        H2F fadc_time4 = new H2F("fadc_time_s4", "fadc_time_s4", 80, 0, 400, 8, 0.5, 8.5);
        fadc_time4.setTitleX("FADC timing");
        fadc_time4.setTitleY("PMT");
        fadc_time4.setTitle("sector 4");
        H2F fadc_time5 = new H2F("fadc_time_s5", "fadc_time_s5", 80, 0, 400, 8, 0.5, 8.5);
        fadc_time5.setTitleX("FADC timing");
        fadc_time5.setTitleY("PMT");
        fadc_time5.setTitle("sector 5");
        H2F fadc_time6 = new H2F("fadc_time_s6", "fadc_time_s6", 80, 0, 400, 8, 0.5, 8.5);
        fadc_time6.setTitleX("FADC timing");
        fadc_time6.setTitleY("PMT");
        fadc_time6.setTitle("sector 6");
        
        H2F tdc1 = new H2F("tdc_s1", "tdc_s1", 400, 0, 40000, 8, 0.5, 8.5);
        tdc1.setTitleX("TDC time");
        tdc1.setTitleY("PMT");
        tdc1.setTitle("sector 1");
        H2F tdc2 = new H2F("tdc_s2", "tdc_s2", 400, 0, 40000, 8, 0.5, 8.5);
        tdc2.setTitleX("TDC time");
        tdc2.setTitleY("PMT");
        tdc2.setTitle("sector 2");
        H2F tdc3 = new H2F("tdc_s3", "tdc_s3", 400, 0, 40000, 8, 0.5, 8.5);
        tdc3.setTitleX("TDC time");
        tdc3.setTitleY("PMT");
        tdc3.setTitle("sector 3");
        H2F tdc4 = new H2F("tdc_s4", "tdc_s4", 400, 0, 40000, 8, 0.5, 8.5);
        tdc4.setTitleX("TDC time");
        tdc4.setTitleY("PMT");
        tdc4.setTitle("sector 4");
        H2F tdc5 = new H2F("tdc_s5", "tdc_s5", 400, 0, 40000, 8, 0.5, 8.5);
        tdc5.setTitleX("TDC time");
        tdc5.setTitleY("PMT");
        tdc5.setTitle("sector 5");
        H2F tdc6 = new H2F("tdc_s6", "tdc_s6", 400, 0, 40000, 8, 0.5, 8.5);
        tdc6.setTitleX("TDC time");
        tdc6.setTitleY("PMT");
        tdc6.setTitle("sector 6");
        
        DataGroup dg = new DataGroup(5,4);
        dg.addDataSet(occADC, 0);
        dg.addDataSet(occTDC, 1);
        dg.addDataSet(adc1, 2);
        dg.addDataSet(adc2, 3);
        dg.addDataSet(adc3, 4);
        dg.addDataSet(adc4, 5);
        dg.addDataSet(adc5, 6);
        dg.addDataSet(adc6, 7);
        dg.addDataSet(fadc_time1, 8);
        dg.addDataSet(fadc_time2, 9);
        dg.addDataSet(fadc_time3, 10);
        dg.addDataSet(fadc_time4, 11);
        dg.addDataSet(fadc_time5, 12);
        dg.addDataSet(fadc_time6, 13);
        dg.addDataSet(tdc1, 14);
        dg.addDataSet(tdc2, 15);
        dg.addDataSet(tdc3, 16);
        dg.addDataSet(tdc4, 17);
        dg.addDataSet(tdc5, 18);
        dg.addDataSet(tdc6, 19);
        this.getDataGroup().add(dg,0,0,0);
    }
        
    @Override
    public void plotHistos() {
        // plotting histos
        this.getDetectorCanvas().getCanvas("Occupancies").cd(0);
        this.getDetectorCanvas().getCanvas("Occupancies").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("occADC"));
        this.getDetectorCanvas().getCanvas("Occupancies").cd(1);
        this.getDetectorCanvas().getCanvas("Occupancies").getPad(1).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("occTDC"));
        this.getDetectorCanvas().getCanvas("Occupancies").update();
        
        this.getDetectorCanvas().getCanvas("ADC spectra").cd(0);
        this.getDetectorCanvas().getCanvas("ADC spectra").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("adc_s1"));
        this.getDetectorCanvas().getCanvas("ADC spectra").cd(1);
        this.getDetectorCanvas().getCanvas("ADC spectra").getPad(1).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("adc_s2"));
        this.getDetectorCanvas().getCanvas("ADC spectra").cd(2);
        this.getDetectorCanvas().getCanvas("ADC spectra").getPad(2).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("adc_s3"));
        this.getDetectorCanvas().getCanvas("ADC spectra").cd(3);
        this.getDetectorCanvas().getCanvas("ADC spectra").getPad(3).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("adc_s4"));
        this.getDetectorCanvas().getCanvas("ADC spectra").cd(4);
        this.getDetectorCanvas().getCanvas("ADC spectra").getPad(4).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("adc_s5"));
        this.getDetectorCanvas().getCanvas("ADC spectra").cd(5);
        this.getDetectorCanvas().getCanvas("ADC spectra").getPad(5).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("adc_s6"));
        this.getDetectorCanvas().getCanvas("ADC spectra").update();
        
        this.getDetectorCanvas().getCanvas("FADC timing spectra").cd(0);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC timing spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s1"));
        this.getDetectorCanvas().getCanvas("FADC timing spectra").cd(1);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").getPad(1).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC timing spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s2"));
        this.getDetectorCanvas().getCanvas("FADC timing spectra").cd(2);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").getPad(2).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC timing spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s3"));
        this.getDetectorCanvas().getCanvas("FADC timing spectra").cd(3);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").getPad(3).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC timing spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s4"));
        this.getDetectorCanvas().getCanvas("FADC timing spectra").cd(4);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").getPad(4).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC timing spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s5"));
        this.getDetectorCanvas().getCanvas("FADC timing spectra").cd(5);
        this.getDetectorCanvas().getCanvas("FADC timing spectra").getPad(5).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC timing spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s6"));
        this.getDetectorCanvas().getCanvas("FADC timing spectra").update();

        this.getDetectorCanvas().getCanvas("TDC spectra").cd(0);
        this.getDetectorCanvas().getCanvas("TDC spectra").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdc_s1"));
        this.getDetectorCanvas().getCanvas("TDC spectra").cd(1);
        this.getDetectorCanvas().getCanvas("TDC spectra").getPad(1).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdc_s2"));
        this.getDetectorCanvas().getCanvas("TDC spectra").cd(2);
        this.getDetectorCanvas().getCanvas("TDC spectra").getPad(2).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdc_s3"));
        this.getDetectorCanvas().getCanvas("TDC spectra").cd(3);
        this.getDetectorCanvas().getCanvas("TDC spectra").getPad(3).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdc_s4"));
        this.getDetectorCanvas().getCanvas("TDC spectra").cd(4);
        this.getDetectorCanvas().getCanvas("TDC spectra").getPad(4).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdc_s5"));
        this.getDetectorCanvas().getCanvas("TDC spectra").cd(5);
        this.getDetectorCanvas().getCanvas("TDC spectra").getPad(5).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("TDC spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdc_s6"));
        this.getDetectorCanvas().getCanvas("TDC spectra").update();
        
        
        
        this.getDetectorView().getView().repaint();
        this.getDetectorView().update();
    }

    @Override
    public void processEvent(DataEvent event) {
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[11] && super.eventResetTime_current[11] > 0){
            resetEventListener();
        }
        
		//if (!testTriggerMask()) return;
        
        // process event info and save into data group
        if(event.hasBank("HTCC::adc")==true){
	    DataBank bank = event.getBank("HTCC::adc");
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
                if(adc>0 && time>0) {
                    this.getDataGroup().getItem(0,0,0).getH2F("occADC").fill(sector*1.0,((comp-1)*2+layer)*1.0);
                    if(sector == 1) this.getDataGroup().getItem(0,0,0).getH2F("adc_s1").fill(adc*1.0,((comp-1)*2+layer)*1.0);
                    if(sector == 1) this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s1").fill(time,((comp-1)*2+layer)*1.0);
                    if(sector == 2) this.getDataGroup().getItem(0,0,0).getH2F("adc_s2").fill(adc*1.0,((comp-1)*2+layer)*1.0);
                    if(sector == 2) this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s2").fill(time,((comp-1)*2+layer)*1.0);
                    if(sector == 3) this.getDataGroup().getItem(0,0,0).getH2F("adc_s3").fill(adc*1.0,((comp-1)*2+layer)*1.0);
                    if(sector == 3) this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s3").fill(time,((comp-1)*2+layer)*1.0);
                    if(sector == 4) this.getDataGroup().getItem(0,0,0).getH2F("adc_s4").fill(adc*1.0,((comp-1)*2+layer)*1.0);
                    if(sector == 4) this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s4").fill(time,((comp-1)*2+layer)*1.0);
                    if(sector == 5) this.getDataGroup().getItem(0,0,0).getH2F("adc_s5").fill(adc*1.0,((comp-1)*2+layer)*1.0);
                    if(sector == 5) this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s5").fill(time,((comp-1)*2+layer)*1.0);
                    if(sector == 6) this.getDataGroup().getItem(0,0,0).getH2F("adc_s6").fill(adc*1.0,((comp-1)*2+layer)*1.0);
                    if(sector == 6) this.getDataGroup().getItem(0,0,0).getH2F("fadc_time_s6").fill(time,((comp-1)*2+layer)*1.0);
                    this.getDetectorSummary().getH1F("summary").fill(sector*1.0);
                }
	    }
    	}
        if(event.hasBank("HTCC::tdc")==true){
            DataBank  bank = event.getBank("HTCC::tdc");
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
                if(tdc>0) this.getDataGroup().getItem(0,0,0).getH2F("occTDC").fill(sector*1.0,((comp-1)*2+layer)*1.0);
                
                if(sector == 1) this.getDataGroup().getItem(0,0,0).getH2F("tdc_s1").fill(tdc,((comp-1)*2+layer)*1.0);
                if(sector == 2) this.getDataGroup().getItem(0,0,0).getH2F("tdc_s2").fill(tdc,((comp-1)*2+layer)*1.0);
                if(sector == 3) this.getDataGroup().getItem(0,0,0).getH2F("tdc_s3").fill(tdc,((comp-1)*2+layer)*1.0);
                if(sector == 4) this.getDataGroup().getItem(0,0,0).getH2F("tdc_s4").fill(tdc,((comp-1)*2+layer)*1.0);
                if(sector == 5) this.getDataGroup().getItem(0,0,0).getH2F("tdc_s5").fill(tdc,((comp-1)*2+layer)*1.0);
                if(sector == 6) this.getDataGroup().getItem(0,0,0).getH2F("tdc_s6").fill(tdc,((comp-1)*2+layer)*1.0);
                
//              this.getDetectorSummary().getH1F("summary").fill(sector*1.0);                
            }
        }        
    }

    @Override
    public void timerUpdate() {

    }


}
