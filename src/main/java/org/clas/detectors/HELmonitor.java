package org.clas.detectors;

import java.util.ArrayList;
import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

public class HELmonitor extends DetectorMonitor {

    
    public HELmonitor(String name) {
        super(name);
        this.setDetectorTabNames("Helicity");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",30,-1.5,1.5);
        summary.setTitleX("helicity");
        summary.setTitleY("counts");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        H1F hel = new H1F("helicity","helicity", 30,-1.5,1.5);
        hel.setTitleX("helicity");
        hel.setTitleY("Counts");
        hel.setFillColor(33);
        
        DataGroup dg = new DataGroup(1,1);
        dg.addDataSet(hel, 0);

        this.getDataGroup().add(dg, 0,0,0);
    }
        
    @Override
    public void plotHistos() {
        // initialize canvas and plot histograms
        this.getDetectorCanvas().getCanvas("Helicity").divide(1, 1);
        this.getDetectorCanvas().getCanvas("Helicity").setGridX(false);
        this.getDetectorCanvas().getCanvas("Helicity").setGridY(false);
        this.getDetectorCanvas().getCanvas("Helicity").cd(0);
        this.getDetectorCanvas().getCanvas("Helicity").draw(this.getDataGroup().getItem(0,0,0).getH1F("helicity"));
        this.getDetectorCanvas().getCanvas("Helicity").update();
        
    }

    @Override
    public void processEvent(DataEvent event) {
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[16] && super.eventResetTime_current[16] > 0){
            resetEventListener();
        }
   
    // process event info and save into data group
        
        if(event.hasBank("HEADER::info")==true){
	    DataBank bank = event.getBank("HEADER::info");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                
                int helicity  = bank.getByte("helicity", loop);

                this.getDataGroup().getItem(0,0,0).getH1F("hel").fill(helicity);
                this.getDetectorSummary().getH1F("summary").fill(helicity);
            }
	}
               
        
    }

    @Override
    public void timerUpdate() {

    }
    
}
