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

public class FCUPmonitor extends DetectorMonitor {
    
    
    public FCUPmonitor(String name) {
        super(name);
        this.setDetectorTabNames("Faraday Cup");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",1000,0.0,1000);
        summary.setTitleX("Faraday Cup");
        summary.setTitleY("counts");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        H1F fcup = new H1F("fcup","fcup", 1000,0.0,1000);
        fcup.setTitleX("faraday cup");
        fcup.setTitleY("Counts");
        H1F fcup_gated = new H1F("fcup_gated","fcup_gated", 1000,0.0,1000);
        fcup_gated.setTitleX("faraday cup gated");
        fcup_gated.setTitleY("Counts");
        
        DataGroup dg = new DataGroup(1,2);
        dg.addDataSet(fcup, 0);
        dg.addDataSet(fcup_gated, 0);

        this.getDataGroup().add(dg, 0,0,0);
    }
        
    @Override
    public void plotHistos() {
        // initialize canvas and plot histograms
        this.getDetectorCanvas().getCanvas("Faraday Cup").divide(1, 2);
        this.getDetectorCanvas().getCanvas("Faraday Cup").setGridX(false);
        this.getDetectorCanvas().getCanvas("Faraday Cup").setGridY(false);
        this.getDetectorCanvas().getCanvas("Faraday Cup").cd(0);
        this.getDetectorCanvas().getCanvas("Faraday Cup").draw(this.getDataGroup().getItem(0,0,0).getH1F("fcup"));
        this.getDetectorCanvas().getCanvas("Faraday Cup").cd(1);
        this.getDetectorCanvas().getCanvas("Faraday Cup").draw(this.getDataGroup().getItem(0,0,0).getH1F("fcup_gated"));
        this.getDetectorCanvas().getCanvas("Faraday Cup").update();
        
    }

    @Override
    public void processEvent(DataEvent event) {

        if (this.getNumberOfEvents() >= super.eventResetTime_current[17] && super.eventResetTime_current[17] > 0){
            resetEventListener();
        }
        
        if(event.hasBank("HEADER::info")==true){
	    DataBank bank = event.getBank("HEADER::info");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                
                int fc   = bank.getByte("fc", loop);
                int fcg  = bank.getByte("fcg", loop);

                this.getDataGroup().getItem(0,0,0).getH1F("fcup").fill(fc);
                this.getDataGroup().getItem(0,0,0).getH1F("fcup_gated").fill(fcg);
                
                this.getDetectorSummary().getH1F("summary").fill(fc);
            }
	}
        
        
    }

    @Override
    public void timerUpdate() {

    }
    
}
