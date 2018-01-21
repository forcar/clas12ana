package org.clas.detectors;

import java.util.ArrayList;

import org.clas.viewer.DetectorMonitor;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

public class TRIGGERmonitor extends DetectorMonitor {
	
	String tbit = "Trigger Bits: ECAL.PCAL.HTCC(0)    ECAL.PCAL.HTCC(1-6)    HTCC(7-12)    PCAL(13-18)    ECAL(19-24)   HT.PC(25)   HT.EC(26)   PC.EC(27)   FTOF.PC(28)  FT100(29)  FT500(30)  1K Pulser(31)";
	
    public TRIGGERmonitor(String name) {
        super(name);
        this.setDetectorTabNames("Trigger Bits", "EC peak trigger", "EC cluster trigger", "HTCC cluster trigger", "FTOF Cluster trigger");
        this.useSectorButtons(true);
        this.init(false);
        this.testTrigger = true;
    }

    
    @Override
    public void createHistos() {
        // create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",100,0.5,10000);
        summary.setFillColor(4);
        summary.setTitleX("trigger");
        summary.setTitleY("counts");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        
        H1F trig = new H1F(tbit,tbit, 32,-0.5,31.5);
        trig.setFillColor(4);      
        trig.setTitleX("Trigger Bits");
        trig.setTitleY("Counts");
        
        DataGroup dg = new DataGroup(4,3);
               
        for(int sec=1; sec <= 6; sec++) {
            
            /// ECAL peak
            
            H1F hist_ecpeak_energy = new H1F("ecpeak_energy"+sec, "ecpeak_energy"+sec, 120, 0, 12);
            hist_ecpeak_energy.setTitleX("EC peak energy");
            hist_ecpeak_energy.setTitleY("counts");
            hist_ecpeak_energy.setTitle("Sector "+sec);
            
            H1F hist_ecpeak_time = new H1F("ecpeak_time"+sec, "ecpeak_time"+sec, 1000, 0, 1000);
            hist_ecpeak_time.setTitleX("EC peak time");
            hist_ecpeak_time.setTitleY("counts");
            hist_ecpeak_time.setTitle("Sector "+sec);
            
            H1F hist_ecpeak_coord = new H1F("ecpeak_coord"+sec, "ecpeak_coord"+sec, 200, 0., 200.);
            hist_ecpeak_coord.setTitleX("EC peak coordinate");
            hist_ecpeak_coord.setTitleY("counts");
            hist_ecpeak_coord.setTitle("Sector "+sec);

            // ECAL cluster
            
            H1F hist_eccluster_energy = new H1F("eccluster_energy"+sec, "eccluster_energy"+sec, 120, 0, 12);
            hist_eccluster_energy.setTitleX("EC cluster energy");
            hist_eccluster_energy.setTitleY("counts");
            hist_eccluster_energy.setTitle("Sector "+sec);
            
            H1F hist_eccluster_time = new H1F("eccluster_time"+sec, "eccluster_time"+sec, 1000, 0, 1000);
            hist_eccluster_time.setTitleX("EC cluster time");
            hist_eccluster_time.setTitleY("counts");
            hist_eccluster_time.setTitle("Sector "+sec);
            
            H2F hist_eccluster_coord = new H2F("eccluster_coord"+sec, "eccluster_coord"+sec, 200, 0, 200, 200, 0, 200);
            hist_eccluster_coord.setTitleX("EC cluster coordinate");
            hist_eccluster_coord.setTitleY("counts");
            hist_eccluster_coord.setTitle("Sector "+sec);
            
            
             /// HTTC
        
            H1F hist_httc_cluster_mask_high = new H1F("httc_cluster_mask_high"+sec,"httc_cluster_mask_high"+sec, 200,0,200);
            hist_httc_cluster_mask_high.setTitleX("httc_cluster_mask_high");
            hist_httc_cluster_mask_high.setTitleY("Counts");
            hist_httc_cluster_mask_high.setTitle("Sector "+sec);
            H1F hist_httc_cluster_mask_low = new H1F("httc_cluster_mask_low"+sec,"httc_cluster_mask_low"+sec, 200,0,200);
            hist_httc_cluster_mask_low.setTitleX("httc_cluster_mask_low");
            hist_httc_cluster_mask_low.setTitleY("Counts");
            hist_httc_cluster_mask_low.setTitle("Sector "+sec);
        
            /// FTOF
        
            H1F hist_ftof_cluster_mask_high = new H1F("ftof_cluster_mask_high"+sec,"ftof_cluster_mask_high"+sec, 200,0,200);
            hist_ftof_cluster_mask_high.setTitleX("ftof_cluster_mask_high");
            hist_ftof_cluster_mask_high.setTitleY("Counts");
            hist_ftof_cluster_mask_high.setTitle("Sector "+sec);
            H1F hist_ftof_cluster_mask_low = new H1F("ftof_cluster_mask_low"+sec,"ftof_cluster_mask_low"+sec, 200,0,200);
            hist_ftof_cluster_mask_low.setTitleX("ftof_cluster_mask_low");
            hist_ftof_cluster_mask_low.setTitleY("Counts");
            hist_ftof_cluster_mask_low.setTitle("Sector "+sec);
            
            
            dg.addDataSet(hist_ecpeak_energy, 0);
            dg.addDataSet(hist_ecpeak_time, 1);
            dg.addDataSet(hist_ecpeak_coord, 2);
            dg.addDataSet(hist_eccluster_energy, 3);
            dg.addDataSet(hist_eccluster_time, 4);
            dg.addDataSet(hist_eccluster_coord, 5);
            dg.addDataSet(hist_httc_cluster_mask_high, 6);
            dg.addDataSet(hist_httc_cluster_mask_low, 7);
            dg.addDataSet(hist_ftof_cluster_mask_high, 8);
            dg.addDataSet(hist_ftof_cluster_mask_low, 9);

            this.getDataGroup().add(dg,sec,0,0);

        }
        
        dg.addDataSet(trig, 10);

        this.getDataGroup().add(dg, 0,0,0);
    }
        
    
    @Override
        public void drawDetector() {

    }
     
        
    @Override
    public void plotHistos() {
        // initialize canvas and plot histograms
        this.getDetectorCanvas().getCanvas("Trigger Bits").divide(1, 1);
        this.getDetectorCanvas().getCanvas("Trigger Bits").setGridX(false);
        this.getDetectorCanvas().getCanvas("Trigger Bits").setGridY(false);
        this.getDetectorCanvas().getCanvas("Trigger Bits").cd(0);
        this.getDetectorCanvas().getCanvas("Trigger Bits").getPad().getAxisY().setLog(true);
        this.getDetectorCanvas().getCanvas("Trigger Bits").draw(this.getDataGroup().getItem(0,0,0).getH1F(tbit));
        this.getDetectorCanvas().getCanvas("Trigger Bits").update();
        
        this.getDetectorCanvas().getCanvas("EC peak trigger").divide(3, 1);
        this.getDetectorCanvas().getCanvas("EC peak trigger").setGridX(false);
        this.getDetectorCanvas().getCanvas("EC peak trigger").setGridY(false);
        this.getDetectorCanvas().getCanvas("EC cluster trigger").divide(3, 1);
        this.getDetectorCanvas().getCanvas("EC cluster trigger").setGridX(false);
        this.getDetectorCanvas().getCanvas("EC cluster trigger").setGridY(false);
        this.getDetectorCanvas().getCanvas("HTCC cluster trigger").divide(1, 2);
        this.getDetectorCanvas().getCanvas("HTCC cluster trigger").setGridX(false);
        this.getDetectorCanvas().getCanvas("HTCC cluster trigger").setGridY(false);
        this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").divide(1, 2);
        this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").setGridX(false);
        this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").setGridY(false);
        
        
        for(int sec=1; sec<7; sec++) {
            if(getActiveSector()==sec) {
                this.getDetectorCanvas().getCanvas("EC peak trigger").cd(0);
                this.getDetectorCanvas().getCanvas("EC peak trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("ecpeak_energy"+sec));
                this.getDetectorCanvas().getCanvas("EC peak trigger").cd(1);
                this.getDetectorCanvas().getCanvas("EC peak trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("ecpeak_time"+sec));
                this.getDetectorCanvas().getCanvas("EC peak trigger").cd(2);
                this.getDetectorCanvas().getCanvas("EC peak trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("ecpeak_coord"+sec));
                this.getDetectorCanvas().getCanvas("EC peak trigger").update();
                
                this.getDetectorCanvas().getCanvas("EC cluster trigger").cd(0);
                this.getDetectorCanvas().getCanvas("EC cluster trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("eccluster_energy"+sec));
                this.getDetectorCanvas().getCanvas("EC cluster trigger").cd(1);
                this.getDetectorCanvas().getCanvas("EC cluster trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("eccluster_time"+sec));
                this.getDetectorCanvas().getCanvas("EC cluster trigger").cd(2);
                this.getDetectorCanvas().getCanvas("EC cluster trigger").draw(this.getDataGroup().getItem(0,0,0).getH2F("eccluster_coord"+sec));
                this.getDetectorCanvas().getCanvas("EC cluster trigger").update();
        
                this.getDetectorCanvas().getCanvas("HTCC cluster trigger").cd(0);
                this.getDetectorCanvas().getCanvas("HTCC cluster trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("httc_cluster_mask_high"+sec));
                this.getDetectorCanvas().getCanvas("HTCC cluster trigger").cd(1);
                this.getDetectorCanvas().getCanvas("HTCC cluster trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("httc_cluster_mask_low"+sec));
                this.getDetectorCanvas().getCanvas("HTCC cluster trigger").update();

                this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").cd(0);
                this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("ftof_cluster_mask_high"+sec));
                this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").cd(1);
                this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").draw(this.getDataGroup().getItem(0,0,0).getH1F("ftof_cluster_mask_low"+sec));
                this.getDetectorCanvas().getCanvas("FTOF Cluster trigger").update();
            }
        }
        
    }

    @Override
    public void processEvent(DataEvent event) {

        if (this.getNumberOfEvents() >= super.eventResetTime_current[18] && super.eventResetTime_current[18] > 0){
            resetEventListener();
        }  
        
		//if (!testTriggerMask()) return;
        
        for (int i=0; i<32; i++) if(isTrigBitSet(i)) this.getDataGroup().getItem(0,0,0).getH1F(tbit).fill(i);
        
        for(int sec=1; sec<=6; sec++) {
            if (isGoodECALTrigger(sec)) {
                this.getDataGroup().getItem(sec,0,0).getH1F("ecpeak_energy"+sec).fill(0);
                this.getDataGroup().getItem(sec,0,0).getH1F("ecpeak_time"+sec).fill(0);
                this.getDataGroup().getItem(sec,0,0).getH1F("ecpeak_coord"+sec).fill(0);
                
                this.getDataGroup().getItem(sec,0,0).getH1F("eccluster_energy"+sec).fill(0);
                this.getDataGroup().getItem(sec,0,0).getH1F("eccluster_time"+sec).fill(0);
                this.getDataGroup().getItem(sec,0,0).getH2F("eccluster_coord"+sec).fill(0,0);
                
                this.getDataGroup().getItem(sec,0,0).getH1F("httc_cluster_mask_high"+sec).fill(0);
                this.getDataGroup().getItem(sec,0,0).getH1F("httc_cluster_mask_low"+sec).fill(0);
                
                this.getDataGroup().getItem(sec,0,0).getH1F("ftof_cluster_mask_high"+sec).fill(0);
                this.getDataGroup().getItem(sec,0,0).getH1F("ftof_cluster_mask_low"+sec).fill(0);
            }
        }
                       
    }

    
    @Override
    public void timerUpdate() {

    }

    
}
