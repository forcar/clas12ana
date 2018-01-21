package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.detector.base.DetectorCollection;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;


public class FTCALmonitor  extends DetectorMonitor {        

    private final int nCrystal = 22;
    private final double crystalWidth = 15;
    private DetectorCollection<Integer> crystals = new DetectorCollection<Integer>();
    

    public FTCALmonitor(String name) {
        super(name);
        
        this.setDetectorTabNames("FADC Occupancies");
        this.init(false);
        this.initCrystalArray();
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        
        H1F summary = new H1F("summary","summary",332,0.5,332.5);
        summary.setTitleX("Crystal");
        summary.setTitleY("FTCAL hits");
        summary.setTitle("FTCAL");
        summary.setFillColor(38);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        
        /*
        H2F statistics = new H2F("statistics","statistics", 100, 0, 100, 350,0.5,350.5);
        statistics.setTitleX("time");
        statistics.setTitleY("number of active channels");
        statistics.setTitle("FTCAL");
        DataGroup stat = new DataGroup(1,1);
        stat.addDataSet(statistics, 0);
        this.setDetectorSummary(stat);
        int active_channels = 0;
        */
        
        H2F occFADC2D = new H2F("occFADC_2D", "occFADC_2D", 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2, 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2);
        occFADC2D.setTitleX("Crystal -X");
        occFADC2D.setTitleY("Crystal Y");
        
        H2F pedFADC2D = new H2F("pedFADC_2D", "pedFADC2D", 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2, 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2);
        pedFADC2D.setTitleX("Crystal -X");
        pedFADC2D.setTitleY("Crystal Y");
        
        H2F pedFADC2Dtmp1 = new H2F("pedFADC_2Dtmp1", "pedFADC_2Dtmp1", 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2, 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2);
        pedFADC2Dtmp1.setTitleX("Crystal -X");
        pedFADC2Dtmp1.setTitleY("Crystal Y");
        
        H2F pedFADC2Dtmp2 = new H2F("pedFADC_2Dtmp2", "pedFADC_2Dtmp2", 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2, 22, -crystalWidth*nCrystal/2, crystalWidth*nCrystal/2);
        pedFADC2Dtmp2.setTitleX("Crystal -X");
        pedFADC2Dtmp2.setTitleY("Crystal Y");
        
        H1F occFADC = new H1F("occFADC", "occFADC", 332, 0.5, 332.5);
        occFADC.setTitleX("Crystal");
        occFADC.setTitleY("Counts");
        occFADC.setFillColor(38);
        
        H2F fadc = new H2F("fadc", "fadc", 50, 0, 5000, 332, 0.5, 332.5);
        fadc.setTitleX("FADC - amplitude");
        fadc.setTitleY("Crystal");
        H2F fadc_time = new H2F("fadc_time", "fadc_time", 50, 0, 500, 332, 0.5, 332.5);
        fadc_time.setTitleX("FADC - time");
        fadc_time.setTitleY("Crystal");
        
        DataGroup dg = new DataGroup(3,3);
        dg.addDataSet(occFADC2D, 0);
        dg.addDataSet(pedFADC2D, 1);
        dg.addDataSet(pedFADC2Dtmp1, 2);
        dg.addDataSet(pedFADC2Dtmp2, 3);
        dg.addDataSet(fadc, 4);
        dg.addDataSet(fadc_time, 5);
        dg.addDataSet(occFADC, 6);
        
        this.getDataGroup().add(dg,0,0,0);
    }
        
    @Override
    public void plotHistos() {        
        // plotting histos
        this.getDetectorCanvas().getCanvas("FADC Occupancies").divide(2, 2);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").setGridX(false);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").setGridY(false);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").cd(0);
//        this.getDetectorCanvas().getCanvas("FADC Occupancies").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("occFADC_2D"));
        this.getDetectorCanvas().getCanvas("FADC Occupancies").cd(1);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").getPad(1).getAxisZ().setRange(0.01, 2.0);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("pedFADC_2D"));
        this.getDetectorCanvas().getCanvas("FADC Occupancies").cd(2);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").getPad(2).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc"));
        this.getDetectorCanvas().getCanvas("FADC Occupancies").cd(3);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").getPad(3).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("FADC Occupancies").draw(this.getDataGroup().getItem(0,0,0).getH2F("fadc_time"));
        this.getDetectorCanvas().getCanvas("FADC Occupancies").update();
    }

    @Override
    public void processEvent(DataEvent event) {

       if (this.getNumberOfEvents() >= super.eventResetTime_current[7] && super.eventResetTime_current[7] > 0){
            resetEventListener();
        }

        // process event info and save into data group
        if(event.hasBank("FTCAL::adc")==true){
	    DataBank bank = event.getBank("FTCAL::adc");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector  = bank.getByte("sector", loop);
                int layer   = bank.getByte("layer", loop);
                int comp    = bank.getShort("component", loop);
                int order   = bank.getByte("order", loop);
                int adc     = bank.getInt("ADC", loop);
                float time  = bank.getFloat("time", loop);
                int ped     = bank.getShort("ped", loop);
                
                int IDY = ((int) comp/22) + 1;
                int IDX = comp + 1 - (IDY -1)*22;    
                double x = (IDX-0.5-nCrystal/2)*crystalWidth;
                double y = (IDY-0.5-nCrystal/2)*crystalWidth;
               
                
             //   System.out.println("ROW " + loop + " SECTOR = " + sector + " LAYER = " + layer + " COMPONENT = " + comp + " ORDER + " + order +
              //        " ADC = " + adc + " TIME = " + time); 
                if(adc>0) {
                        this.getDataGroup().getItem(0,0,0).getH2F("occFADC_2D").fill(-x,y);
                        this.getDataGroup().getItem(0,0,0).getH2F("pedFADC_2Dtmp1").fill(-x,y,ped);
                        this.getDataGroup().getItem(0,0,0).getH2F("pedFADC_2Dtmp2").fill(-x,y,ped*ped);
                        this.getDataGroup().getItem(0,0,0).getH1F("occFADC").fill(crystals.get(1, 1, comp));
                        this.getDataGroup().getItem(0,0,0).getH2F("fadc").fill(adc*1.0,crystals.get(1, 1, comp));
                        this.getDataGroup().getItem(0,0,0).getH2F("fadc_time").fill(time*1.0,crystals.get(1, 1, comp));
                        this.getDetectorSummary().getH1F("summary").fill(crystals.get(1, 1, comp));                          
                }
                
	    }
    	}
            
    }

    @Override
    public void timerUpdate() {
        if(this.getNumberOfEvents()>0) {
            H2F rawN = this.getDataGroup().getItem(0,0,0).getH2F("occFADC_2D");
            H2F raw1 = this.getDataGroup().getItem(0,0,0).getH2F("pedFADC_2Dtmp1");
            H2F raw2 = this.getDataGroup().getItem(0,0,0).getH2F("pedFADC_2Dtmp2");
            H2F ped = this.getDataGroup().getItem(0,0,0).getH2F("pedFADC_2D");
            for(int loop = 0; loop < rawN.getDataBufferSize(); loop++){
                double nev   = rawN.getDataBufferBin(loop);
                double noise = 0;
                if(nev>0) {
                    double ped1 = raw1.getDataBufferBin(loop)/nev;
                    double ped2 = raw2.getDataBufferBin(loop)/nev;
                    noise = Math.sqrt(ped2-ped1*ped1);
                }
                ped.setDataBufferBin(loop,(float) noise);
            }
        }

    }

    private void initCrystalArray(){
        int icrystal=1;
        for (int component = 0; component < nCrystal*nCrystal; component++) {
            if(doesThisCrystalExist(component)) {
                icrystal++;
                crystals.add(1, 1, component, icrystal);
            }
        }
    }
    
    private boolean doesThisCrystalExist(int id) {

        boolean crystalExist=false;
        int iy = id / nCrystal;
        int ix = id - iy * nCrystal;

        double xcrystal = (nCrystal - ix - 0.5);
        double ycrystal = (nCrystal - iy - 0.5);
        double rcrystal = Math.sqrt(Math.pow(xcrystal - 11, 2.0) + Math.pow(ycrystal - 11, 2.0));
        if (rcrystal > 4 && rcrystal < 11) {
            crystalExist=true;
        }
        return crystalExist;
    }


}