package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;


public class BSTmonitor extends DetectorMonitor {
    

    public BSTmonitor(String name) {
        super(name);
        this.setDetectorTabNames("Occupancies 2D", "Occupancies 1D", "Multiplicity");
        this.init(false);
    }
    
    @Override
    public void createHistos() {
        // create histograms
        this.setNumberOfEvents(0);
        H2F summary = new H2F("summary","summary",256,0.5,256.5, 84,0.5,84.5);
        summary.setTitleX("strip");
        summary.setTitleY("module");
        summary.setTitle("BST");
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);

        H2F occ_reg1l1 = new H2F("occ_reg1_l1", "occ_reg1_l1", 256, 0.5, 256.5, 10, 0.5, 10.5);
        occ_reg1l1.setTitleX("strip");
        occ_reg1l1.setTitleY("sector");
        occ_reg1l1.setTitle("region 1 - bottom layer");

        H2F occ_reg2l1 = new H2F("occ_reg2_l1", "occ_reg2_l1", 256, 0.5, 256.5, 14, 0.5, 14.5);
        occ_reg2l1.setTitleX("strip");
        occ_reg2l1.setTitleY("sector");
        occ_reg2l1.setTitle("region 2 - bottom layer");
        
        H2F occ_reg3l1 = new H2F("occ_reg3_l1", "occ_reg3_l1", 256, 0.5, 256.5, 18, 0.5, 18.5);
        occ_reg3l1.setTitleX("strip");
        occ_reg3l1.setTitleY("sector");
        occ_reg3l1.setTitle("region 3 - bottom layer");

        H2F occ_reg1l2 = new H2F("occ_reg1_l2", "occ_reg1_l2", 256, 0.5, 256.5, 10, 0.5, 10.5);
        occ_reg1l2.setTitleX("strip");
        occ_reg1l2.setTitleY("sector");
        occ_reg1l2.setTitle("region 1 - top layer");

        H2F occ_reg2l2 = new H2F("occ_reg2_l2", "occ_reg2_l2", 256, 0.5, 256.5, 14, 0.5, 14.5);
        occ_reg2l2.setTitleX("strip");
        occ_reg2l2.setTitleY("sector");
        occ_reg2l2.setTitle("region 2 - top layer");
        
        H2F occ_reg3l2 = new H2F("occ_reg3_l2", "occ_reg3_l2", 256, 0.5, 256.5, 18, 0.5, 18.5);
        occ_reg3l2.setTitleX("strip");
        occ_reg3l2.setTitleY("sector");
        occ_reg3l2.setTitle("region 3 - top layer");
        
        
        H1F occ_reg1l11d = new H1F("occ_reg1_l1_1d", "occ_reg1_l1_1d", 2560, 0.5, 2560.5);
        occ_reg1l11d.setTitleX("strip");
        occ_reg1l11d.setTitleY("hits");
        occ_reg1l11d.setTitle("region 1 - bottom layer");

        H1F occ_reg2l11d = new H1F("occ_reg2_l1_1d", "occ_reg2_l1_1d", 3584, 0.5, 3584.5);
        occ_reg2l11d.setTitleX("strip");
        occ_reg2l11d.setTitleY("hits");
        occ_reg2l11d.setTitle("region 2 - bottom layer");
        
        H1F occ_reg3l11d = new H1F("occ_reg3_l1_1d", "occ_reg3_l1_1d", 4608, 0.5, 4608.5);
        occ_reg3l11d.setTitleX("strip");
        occ_reg3l11d.setTitleY("hits");
        occ_reg3l11d.setTitle("region 3 - bottom layer");

        H1F occ_reg1l21d = new H1F("occ_reg1_l2_1d", "occ_reg1_l2_1d", 2560, 0.5, 2560.5);
        occ_reg1l21d.setTitleX("strip");
        occ_reg1l21d.setTitleY("hits");
        occ_reg1l21d.setTitle("region 1 - top layer");

        H1F occ_reg2l21d = new H1F("occ_reg2_l2_1d", "occ_reg2_l2_1d", 3584, 0.5, 3584.5);
        occ_reg2l21d.setTitleX("strip");
        occ_reg2l21d.setTitleY("hits");
        occ_reg2l21d.setTitle("region 2 - top layer");
        
        H1F occ_reg3l21d = new H1F("occ_reg3_l2_1d", "occ_reg3_l2_1d", 4608, 0.5, 4608.5);
        occ_reg3l21d.setTitleX("strip");
        occ_reg3l21d.setTitleY("hits");
        occ_reg3l21d.setTitle("region 3 - top layer");
        
        H1F multi = new H1F("multi", "multi", 350, -0.5, 349.5);
        multi.setTitleX("hit multiplicity");
        multi.setTitleY("counts");
        multi.setTitle("Multiplicity of BST channels");
        
        DataGroup dg = new DataGroup(4,4);
        dg.addDataSet(occ_reg1l1, 0);
        dg.addDataSet(occ_reg2l1, 1);
        dg.addDataSet(occ_reg3l1, 2);
        dg.addDataSet(occ_reg1l2, 3);
        dg.addDataSet(occ_reg2l2, 4);
        dg.addDataSet(occ_reg3l2, 5);
        dg.addDataSet(occ_reg1l11d, 6);
        dg.addDataSet(occ_reg2l11d, 7);
        dg.addDataSet(occ_reg3l11d, 8);
        dg.addDataSet(occ_reg1l21d, 9);
        dg.addDataSet(occ_reg2l21d, 10);
        dg.addDataSet(occ_reg3l21d, 11);
        dg.addDataSet(multi, 12);
        this.getDataGroup().add(dg,0,0,0);
        
    }
        
    @Override
    public void plotHistos() {
        // initialize canvas and plot histograms
        
        this.getDetectorCanvas().getCanvas("Occupancies 2D").divide(3, 2);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").setGridX(false);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").setGridY(false);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").cd(0);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies 2D").draw(this.getDataGroup().getItem(0,0,0).getH2F("occ_reg1_l1"));
        this.getDetectorCanvas().getCanvas("Occupancies 2D").cd(1);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").getPad(1).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies 2D").draw(this.getDataGroup().getItem(0,0,0).getH2F("occ_reg2_l1"));
        this.getDetectorCanvas().getCanvas("Occupancies 2D").cd(2);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").getPad(2).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies 2D").draw(this.getDataGroup().getItem(0,0,0).getH2F("occ_reg3_l1"));
        this.getDetectorCanvas().getCanvas("Occupancies 2D").cd(3);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").getPad(3).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies 2D").draw(this.getDataGroup().getItem(0,0,0).getH2F("occ_reg1_l2"));
        this.getDetectorCanvas().getCanvas("Occupancies 2D").cd(4);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").getPad(4).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies 2D").draw(this.getDataGroup().getItem(0,0,0).getH2F("occ_reg2_l2"));
        this.getDetectorCanvas().getCanvas("Occupancies 2D").cd(5);
        this.getDetectorCanvas().getCanvas("Occupancies 2D").getPad(5).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies 2D").draw(this.getDataGroup().getItem(0,0,0).getH2F("occ_reg3_l2"));
        this.getDetectorCanvas().getCanvas("Occupancies 2D").update();
        
        this.getDetectorCanvas().getCanvas("Occupancies 1D").divide(3, 2);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").setGridX(false);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").setGridY(false);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").cd(0);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").draw(this.getDataGroup().getItem(0,0,0).getH1F("occ_reg1_l1_1d"));
        this.getDetectorCanvas().getCanvas("Occupancies 1D").cd(1);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").draw(this.getDataGroup().getItem(0,0,0).getH1F("occ_reg2_l1_1d"));
        this.getDetectorCanvas().getCanvas("Occupancies 1D").cd(2);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").draw(this.getDataGroup().getItem(0,0,0).getH1F("occ_reg3_l1_1d"));
        this.getDetectorCanvas().getCanvas("Occupancies 1D").cd(3);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").draw(this.getDataGroup().getItem(0,0,0).getH1F("occ_reg1_l2_1d"));
        this.getDetectorCanvas().getCanvas("Occupancies 1D").cd(4);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").draw(this.getDataGroup().getItem(0,0,0).getH1F("occ_reg2_l2_1d"));
        this.getDetectorCanvas().getCanvas("Occupancies 1D").cd(5);
        this.getDetectorCanvas().getCanvas("Occupancies 1D").draw(this.getDataGroup().getItem(0,0,0).getH1F("occ_reg3_l2_1d"));
        this.getDetectorCanvas().getCanvas("Occupancies 1D").update();
        
        this.getDetectorCanvas().getCanvas("Multiplicity").divide(1, 1);
        this.getDetectorCanvas().getCanvas("Multiplicity").setGridX(false);
        this.getDetectorCanvas().getCanvas("Multiplicity").setGridY(false);
        this.getDetectorCanvas().getCanvas("Multiplicity").cd(0);
        this.getDetectorCanvas().getCanvas("Multiplicity").draw(this.getDataGroup().getItem(0,0,0).getH1F("multi"));
        this.getDetectorCanvas().getCanvas("Multiplicity").update();
        
        
    }

    @Override
    public void processEvent(DataEvent event) {
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[1] && super.eventResetTime_current[1] > 0){
	    resetEventListener();
        }
        
		//if (!testTriggerMask()) return;
        
        // process event info and save into data group
        if(event.hasBank("BST::adc")==true){
            DataBank  bank = event.getBank("BST::adc");
            this.getDetectorOccupancy().addTDCBank(bank);
            int rows = bank.rows();
            
            this.getDataGroup().getItem(0,0,0).getH1F("multi").fill(rows);
            
            for(int i = 0; i < rows; i++){
                int    sector = bank.getByte("sector",i);
                int     layer = bank.getByte("layer",i);
                int      comp = bank.getShort("component",i);
                int       ADC = bank.getInt("ADC",i);
                float    time = bank.getFloat("time",i);
                int     order = bank.getByte("order",i);
             
                if(ADC > 0){
                    
                    if(layer == 1){ 
                        this.getDataGroup().getItem(0,0,0).getH2F("occ_reg1_l1").fill(comp, sector);
                        this.getDataGroup().getItem(0,0,0).getH1F("occ_reg1_l1_1d").fill((sector-1)*256+comp);
                        this.getDetectorSummary().getH2F("summary").fill(comp, sector);
                    }
                    if(layer == 3){ 
                        this.getDataGroup().getItem(0,0,0).getH2F("occ_reg2_l1").fill(comp, sector);
                        this.getDataGroup().getItem(0,0,0).getH1F("occ_reg2_l1_1d").fill((sector-1)*256+comp);
                        this.getDetectorSummary().getH2F("summary").fill(comp, sector+20);
                    }
                    if(layer == 5){ 
                        this.getDataGroup().getItem(0,0,0).getH2F("occ_reg3_l1").fill(comp, sector);
                        this.getDataGroup().getItem(0,0,0).getH1F("occ_reg3_l1_1d").fill((sector-1)*256+comp);
                        this.getDetectorSummary().getH2F("summary").fill(comp, sector+48);
                    }
                    if(layer == 2){ 
                        this.getDataGroup().getItem(0,0,0).getH2F("occ_reg1_l2").fill(comp, sector);
                        this.getDataGroup().getItem(0,0,0).getH1F("occ_reg1_l2_1d").fill((sector-1)*256+comp);
                        this.getDetectorSummary().getH2F("summary").fill(comp, sector+10);
                    }
                    if(layer == 4){ 
                        this.getDataGroup().getItem(0,0,0).getH2F("occ_reg2_l2").fill(comp, sector);
                        this.getDataGroup().getItem(0,0,0).getH1F("occ_reg2_l2_1d").fill((sector-1)*256+comp);
                        this.getDetectorSummary().getH2F("summary").fill(comp, sector+34);
                    }
                    if(layer == 6){ 
                        this.getDataGroup().getItem(0,0,0).getH2F("occ_reg3_l2").fill(comp, sector);
                        this.getDataGroup().getItem(0,0,0).getH1F("occ_reg3_l2_1d").fill((sector-1)*256+comp);
                        this.getDetectorSummary().getH2F("summary").fill(comp, sector+66);
                    }
                    
                }
            }
       }
                
    }
    

    @Override
    public void timerUpdate() {
        
    }

}
