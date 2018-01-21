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


public class RFmonitor extends DetectorMonitor {
    
    private double tdc2Time = 0.023436;
    
    public RFmonitor(String name) {
        super(name);
        this.setDetectorTabNames("RF TDCs","RF Time");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",6,0.5,6.5);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        H1F rf1 = new H1F("rf1","rf1", 100,0.,120000);
        rf1.setTitleX("RF1 tdc");
        rf1.setTitleY("Counts");
        rf1.setFillColor(33);
        H1F rf2 = new H1F("rf2","rf2", 100,0.,120000);
        rf2.setTitleX("RF2 tdc");
        rf2.setTitleY("Counts");
        rf2.setFillColor(36);
        H1F rfdiff = new H1F("rfdiff","rfdiff", 160, 2.,5.);
        rfdiff.setTitleX("RF diff");
        rfdiff.setTitleY("Counts");
        F1D fdiff = new F1D("fdiff","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        fdiff.setParameter(0, 0);
        fdiff.setParameter(1, 0);
        fdiff.setParameter(2, 1.0);
        fdiff.setLineWidth(2);
        fdiff.setLineColor(2);
        fdiff.setOptStat("1111");
        H1F rfdiffAve = new H1F("rfdiffAve","rfdiffAve", 160, 3.,4.);
        rfdiffAve.setTitleX("RF diff");
        rfdiffAve.setTitleY("Counts");
        F1D fdiffAve = new F1D("fdiffAve","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        fdiffAve.setParameter(0, 0);
        fdiffAve.setParameter(1, 0);
        fdiffAve.setParameter(2, 1.0);
        fdiffAve.setLineWidth(2);
        fdiffAve.setLineColor(2);
        fdiffAve.setOptStat("1111");
        H1F rf1rawdiff = new H1F("rf1rawdiff","rf1rawdiff", 100, 6800.,6900.);
        rf1rawdiff.setTitleX("RF1 diff");
        rf1rawdiff.setTitleY("Counts");
        F1D f1rawdiff = new F1D("f1rawdiff","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f1rawdiff.setParameter(0, 0);
        f1rawdiff.setParameter(1, 0);
        f1rawdiff.setParameter(2, 1.0);
        f1rawdiff.setLineWidth(2);
        f1rawdiff.setLineColor(2);
        f1rawdiff.setOptStat("1111");
        H1F rf2rawdiff = new H1F("rf2rawdiff","rf2rawdiff", 100, 6800.,6900.);
        rf2rawdiff.setTitleX("RF2 diff");
        rf2rawdiff.setTitleY("Counts");
        F1D f2rawdiff = new F1D("f2rawdiff","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f2rawdiff.setParameter(0, 0);
        f2rawdiff.setParameter(1, 0);
        f2rawdiff.setParameter(2, 1.0);
        f2rawdiff.setLineWidth(2);
        f2rawdiff.setLineColor(2);
        f2rawdiff.setOptStat("1111");
        H2F rf1rawdiffrf1 = new H2F("rf1rawdiffrf1","rf1rawdiffrf1", 100,0.,120000, 25, 6800.,6900.);
        rf1rawdiffrf1.setTitleX("RF1 tdc");
        rf1rawdiffrf1.setTitleY("RF1 diff");
        H2F rf2rawdiffrf2 = new H2F("rf2rawdiffrf2","rf2rawdiffrf2", 100,0.,120000, 25, 6800.,6900.);
        rf2rawdiffrf2.setTitleX("RF2 tdc");
        rf2rawdiffrf2.setTitleY("RF2 diff");
        H1F rf1diff = new H1F("rf1diff","rf1diff", 160, 158.,162.);
        rf1diff.setTitleX("RF1 diff (ns)");
        rf1diff.setTitleY("Counts");
        F1D f1diff = new F1D("f1diff","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f1diff.setParameter(0, 0);
        f1diff.setParameter(1, 0);
        f1diff.setParameter(2, 1.0);
        f1diff.setLineWidth(2);
        f1diff.setLineColor(2);
        f1diff.setOptStat("1111");
        H1F rf2diff = new H1F("rf2diff","rf2diff", 160, 158.,162.);
        rf2diff.setTitleX("RF2 diff (ns)");
        rf2diff.setTitleY("Counts");
        F1D f2diff = new F1D("f2diff","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f2diff.setParameter(0, 0);
        f2diff.setParameter(1, 0);
        f2diff.setParameter(2, 1.0);
        f2diff.setLineWidth(2);
        f2diff.setLineColor(2);
        f2diff.setOptStat("1111");
        H2F timeRF1 = new H2F("timeRF1","timeRF1",100,0.,240, 100, 3., 4.);
        timeRF1.setTitleX("RF1 (ns)");
        timeRF1.setTitleY("RF diff (ns)");
        H2F timeRF2 = new H2F("timeRF2","timeRF2",100,0.,240, 100, 3., 4.);
        timeRF2.setTitleX("RF2 (ns)");
        timeRF2.setTitleY("RF diff (ns)");
 
        DataGroup dg = new DataGroup(1,12);
        dg.addDataSet(rf1, 0);
        dg.addDataSet(rf1rawdiff,1);
        dg.addDataSet(f1rawdiff, 1);
        dg.addDataSet(rf1rawdiffrf1,2);
        dg.addDataSet(rf2, 3);
        dg.addDataSet(rf2rawdiff,4);
        dg.addDataSet(f2rawdiff, 4);
        dg.addDataSet(rf2rawdiffrf2,5);        
        dg.addDataSet(rfdiff, 6);
        dg.addDataSet(fdiff,  6);
        dg.addDataSet(rf1diff,7);
        dg.addDataSet(f1diff, 7);
        dg.addDataSet(rf2diff,8);
        dg.addDataSet(f2diff, 8);
        dg.addDataSet(rfdiffAve, 9);
        dg.addDataSet(fdiffAve,  9);
        dg.addDataSet(timeRF1, 10);
        dg.addDataSet(timeRF2, 11);
        this.getDataGroup().add(dg, 0,0,0);
    }
        
    @Override
    public void plotHistos() {
        // initialize canvas and plot histograms
        this.getDetectorCanvas().getCanvas("RF TDCs").divide(3, 2);
        this.getDetectorCanvas().getCanvas("RF TDCs").setGridX(false);
        this.getDetectorCanvas().getCanvas("RF TDCs").setGridY(false);
        this.getDetectorCanvas().getCanvas("RF TDCs").cd(0);
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getH1F("rf1"));
        this.getDetectorCanvas().getCanvas("RF TDCs").cd(1);
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getH1F("rf1rawdiff"));
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getF1D("f1rawdiff"),"same");
        this.getDetectorCanvas().getCanvas("RF TDCs").cd(2);
        this.getDetectorCanvas().getCanvas("RF TDCs").getPad(2).getAxisZ().setLog(true);
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getH2F("rf1rawdiffrf1"));
        this.getDetectorCanvas().getCanvas("RF TDCs").cd(3);
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getH1F("rf2"));
        this.getDetectorCanvas().getCanvas("RF TDCs").cd(4);
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getH1F("rf2rawdiff"));
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getF1D("f2rawdiff"),"same");
        this.getDetectorCanvas().getCanvas("RF TDCs").cd(5);
        this.getDetectorCanvas().getCanvas("RF TDCs").getPad(5).getAxisZ().setLog(true);
        this.getDetectorCanvas().getCanvas("RF TDCs").draw(this.getDataGroup().getItem(0,0,0).getH2F("rf2rawdiffrf2"));
        this.getDetectorCanvas().getCanvas("RF TDCs").update();
        this.getDetectorCanvas().getCanvas("RF Time").divide(3, 2);
        this.getDetectorCanvas().getCanvas("RF Time").setGridX(false);
        this.getDetectorCanvas().getCanvas("RF Time").setGridY(false);
        this.getDetectorCanvas().getCanvas("RF Time").cd(0);
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getH1F("rfdiff"));
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getF1D("fdiff"),"same");
        this.getDetectorCanvas().getCanvas("RF Time").cd(1);
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getH1F("rf1diff"));
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getF1D("f1diff"),"same");
        this.getDetectorCanvas().getCanvas("RF Time").cd(2);
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getH1F("rf2diff"));
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getF1D("f2diff"),"same");
        this.getDetectorCanvas().getCanvas("RF Time").cd(3);
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getH1F("rfdiffAve"));
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getF1D("fdiffAve"),"same");
        this.getDetectorCanvas().getCanvas("RF Time").cd(4);
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getH2F("timeRF1"));
        this.getDetectorCanvas().getCanvas("RF Time").cd(5);
        this.getDetectorCanvas().getCanvas("RF Time").draw(this.getDataGroup().getItem(0,0,0).getH2F("timeRF2"));
        this.getDetectorCanvas().getCanvas("RF Time").update();
    }

    @Override
    public void processEvent(DataEvent event) {
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[15] && super.eventResetTime_current[15] > 0){
            resetEventListener();
        }
        
        // process event info and save into data group
        ArrayList<Integer> rf1 = new ArrayList();
        ArrayList<Integer> rf2 = new ArrayList();
        if(event.hasBank("RF::tdc")==true){
            DataBank  bank = event.getBank("RF::tdc");
            this.getDetectorOccupancy().addTDCBank(bank);
            int rows = bank.rows();
            for(int i = 0; i < rows; i++){
                int    sector = bank.getByte("sector",i);
                int     layer = bank.getByte("layer",i);
                int      comp = bank.getShort("component",i);
                int       TDC = bank.getInt("TDC",i);
                int     order = bank.getByte("order",i); 
                if(order==2) {
                    if(comp==1) {
                        this.getDataGroup().getItem(0,0,0).getH1F("rf1").fill(TDC*1.0);
                        rf1.add(TDC);
                    }
                    else {
                        this.getDataGroup().getItem(0,0,0).getH1F("rf2").fill(TDC*1.0);
                        rf2.add(TDC);
                    }
                }
                this.getDetectorSummary().getH1F("summary").fill(sector*1.0);
            }
        }
//        System.out.println(rf1.size() + " " +rf2.size());
        for(int i=0; i<rf1.size()-1; i++) {
            this.getDataGroup().getItem(0,0,0).getH1F("rf1rawdiff").fill((rf1.get(i+1)-rf1.get(i))*1.0);
            this.getDataGroup().getItem(0,0,0).getH2F("rf1rawdiffrf1").fill(rf1.get(i),(rf1.get(i+1)-rf1.get(i))*1.0);
            this.getDataGroup().getItem(0,0,0).getH1F("rf1diff").fill((rf1.get(i+1)-rf1.get(i))*tdc2Time);
        }
        for(int i=0; i<rf2.size()-1; i++) {
            this.getDataGroup().getItem(0,0,0).getH1F("rf2rawdiff").fill((rf2.get(i+1)-rf2.get(i))*1.0);
            this.getDataGroup().getItem(0,0,0).getH2F("rf2rawdiffrf2").fill(rf2.get(i),(rf2.get(i+1)-rf2.get(i))*1.0);
            this.getDataGroup().getItem(0,0,0).getH1F("rf2diff").fill((rf2.get(i+1)-rf2.get(i))*tdc2Time);
        }

        if(rf1.size()==rf2.size()) {
            double rfTime1 = 0;
            double rfTime2 = 0;
            for(int i=0; i<rf1.size(); i++) {
                this.getDataGroup().getItem(0,0,0).getH1F("rfdiff").fill((rf1.get(i)-rf2.get(i))*tdc2Time);
                rfTime1 += rf1.get(i)*tdc2Time - i*80*2.004;
                rfTime2 += rf2.get(i)*tdc2Time - i*80*2.004;
            }
            rfTime1 /=rf1.size();
            rfTime2 /=rf2.size();            
            this.getDataGroup().getItem(0,0,0).getH1F("rfdiffAve").fill(rfTime1-rfTime2);
            this.getDataGroup().getItem(0,0,0).getH2F("timeRF1").fill(rfTime1,rfTime1-rfTime2);
            this.getDataGroup().getItem(0,0,0).getH2F("timeRF2").fill(rfTime2,rfTime1-rfTime2);            
        }
    }

    @Override
    public void timerUpdate() {
//        System.out.println("Updating RF");
        this.fitRF(this.getDataGroup().getItem(0,0,0).getH1F("rf1rawdiff"),this.getDataGroup().getItem(0,0,0).getF1D("f1rawdiff"));
        this.fitRF(this.getDataGroup().getItem(0,0,0).getH1F("rf2rawdiff"),this.getDataGroup().getItem(0,0,0).getF1D("f2rawdiff"));
        this.fitRF(this.getDataGroup().getItem(0,0,0).getH1F("rf1diff"),   this.getDataGroup().getItem(0,0,0).getF1D("f1diff"));
        this.fitRF(this.getDataGroup().getItem(0,0,0).getH1F("rf2diff"),   this.getDataGroup().getItem(0,0,0).getF1D("f2diff"));
        this.fitRF(this.getDataGroup().getItem(0,0,0).getH1F("rfdiff"),    this.getDataGroup().getItem(0,0,0).getF1D("fdiff"));
        this.fitRF(this.getDataGroup().getItem(0,0,0).getH1F("rfdiffAve"), this.getDataGroup().getItem(0,0,0).getF1D("fdiffAve"));
    }

    public void fitRF(H1F hirf, F1D f1rf) {
        double mean  = hirf.getDataX(hirf.getMaximumBin());
        double amp   = hirf.getBinContent(hirf.getMaximumBin());
        double sigma = hirf.getRMS();
        f1rf.setParameter(0, amp);
        f1rf.setParameter(1, mean);
        f1rf.setParameter(2, sigma);
        f1rf.setRange(mean-3.*sigma,mean+3.*sigma);
        DataFitter.fit(f1rf, hirf, "Q"); //No options uses error for sigma        
    }
    
}
