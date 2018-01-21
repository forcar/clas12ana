package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;


public class FTTRKmonitor  extends DetectorMonitor {  
    
    private int nstrip = 768;
    private int nlayer = 4;
    private int numberOfStripsPerChip = 64 ;
    private int numberOfChips = 12;
    
    public FTTRKmonitor(String name) {
        super(name);
        
        this.setDetectorTabNames("Occupancies_2D", "Occupancies_1D", "Average Time Maximum", "ADC and time spectra");
        this.init(false);
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        
        H1F summary = new H1F("summary","summary",nstrip*nlayer,0.5,nstrip*nlayer+0.5);
        summary.setTitleX("channel");
        summary.setTitleY("FTTRK hits");
        summary.setTitle("FTTRK");
        summary.setFillColor(38);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        
        H2F occADC2D = new H2F("occADC_2D", "occADC_2D", nstrip, 0.5, nstrip+0.5, nlayer, 0.5, nlayer+0.5);
        occADC2D.setTitleX("channel");
        occADC2D.setTitleY("layer");
        
        for(int ilayer=1; ilayer<=nlayer; ilayer++) {
            H1F occADCl = new H1F("occADC_layer" + ilayer, "occADC_layer"  + ilayer, nstrip, 0.5, nstrip+0.5);
            occADCl.setTitleX("channel");
            occADCl.setTitleY("Counts");
            occADCl.setFillColor(38);
            occADCl.setTitle("layer " + ilayer);
            H1F timeMax = new H1F("timeMax_layer" + ilayer, "timeMax_layer"  + ilayer, numberOfChips, 0.5, numberOfChips+0.5);
            timeMax.setTitleX("chip");
            timeMax.setTitleY("Counts");
            timeMax.setFillColor(4);
            timeMax.setTitle("layer " + ilayer);            
            H1F timeMaxTmp1 = new H1F("timeMaxTmp1_layer" + ilayer, "timeMaxTmp1_layer"  + ilayer, numberOfChips, 0.5, numberOfChips+0.5);
            timeMaxTmp1.setTitleX("chip");
            timeMaxTmp1.setTitleY("Counts");
            timeMaxTmp1.setFillColor(4);
            timeMaxTmp1.setTitle("layer " + ilayer); 
            H1F timeMaxTmp2 = new H1F("timeMaxTmp2_layer" + ilayer, "timeMaxTmp2_layer"  + ilayer, numberOfChips, 0.5, numberOfChips+0.5);
            timeMaxTmp2.setTitleX("chip");
            timeMaxTmp2.setTitleY("Counts");
            timeMaxTmp2.setFillColor(4);
            timeMaxTmp2.setTitle("layer " + ilayer);             
            DataGroup dg = new DataGroup(1,2);
            dg.addDataSet(occADCl, 0);
            dg.addDataSet(timeMax, 1);
            dg.addDataSet(timeMaxTmp1, 1);
            dg.addDataSet(timeMaxTmp2, 1);
            this.getDataGroup().add(dg,0,ilayer,0);
        }
       
        H2F adc = new H2F("adc", "adc", 50, 0, 1000, nstrip*nlayer,0.5,nstrip*nlayer+0.5);
        adc.setTitleX("ADC - amplitude");
        adc.setTitleY("channel");
        H2F tdc = new H2F("tdc", "tdc", 50, 0, 1000, nstrip*nlayer,0.5,nstrip*nlayer+0.5);
        tdc.setTitleX("time");
        tdc.setTitleY("channel");
        
        DataGroup dg = new DataGroup(1,3);
        dg.addDataSet(occADC2D, 0);
        dg.addDataSet(adc, 1);
        dg.addDataSet(tdc, 1);
        this.getDataGroup().add(dg,0,0,0);
    }
        
    @Override
    public void plotHistos() {        
        // plotting histos
        this.getDetectorCanvas().getCanvas("Occupancies_2D").divide(1, 1);
        this.getDetectorCanvas().getCanvas("Occupancies_2D").setGridX(false);
        this.getDetectorCanvas().getCanvas("Occupancies_2D").setGridY(false);
        this.getDetectorCanvas().getCanvas("Occupancies_2D").cd(0);
        this.getDetectorCanvas().getCanvas("Occupancies_2D").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("Occupancies_2D").draw(this.getDataGroup().getItem(0,0,0).getH2F("occADC_2D"));
        this.getDetectorCanvas().getCanvas("Occupancies_1D").divide(2, 2);
        this.getDetectorCanvas().getCanvas("Occupancies_1D").setGridX(false);
        this.getDetectorCanvas().getCanvas("Occupancies_1D").setGridY(false);
        for(int ilayer=1; ilayer<=nlayer; ilayer++) {        
            this.getDetectorCanvas().getCanvas("Occupancies_1D").cd(0 + ilayer -1);
            this.getDetectorCanvas().getCanvas("Occupancies_1D").draw(this.getDataGroup().getItem(0,ilayer,0).getH1F("occADC_layer" + ilayer));
        }
        this.getDetectorCanvas().getCanvas("Average Time Maximum").divide(2, 2);
        this.getDetectorCanvas().getCanvas("Average Time Maximum").setGridX(false);
        this.getDetectorCanvas().getCanvas("Average Time Maximum").setGridY(false);
        for(int ilayer=1; ilayer<=nlayer; ilayer++) {        
            this.getDetectorCanvas().getCanvas("Average Time Maximum").cd(0 + ilayer -1);
            this.getDetectorCanvas().getCanvas("Average Time Maximum").draw(this.getDataGroup().getItem(0,ilayer,0).getH1F("timeMax_layer" + ilayer));
        }
        this.getDetectorCanvas().getCanvas("ADC and time spectra").divide(1, 2);
        this.getDetectorCanvas().getCanvas("ADC and time spectra").setGridX(false);
        this.getDetectorCanvas().getCanvas("ADC and time spectra").setGridY(false);
        this.getDetectorCanvas().getCanvas("ADC and time spectra").cd(0);
        this.getDetectorCanvas().getCanvas("ADC and time spectra").getPad(0).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC and time spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("adc"));
        this.getDetectorCanvas().getCanvas("ADC and time spectra").cd(1);
        this.getDetectorCanvas().getCanvas("ADC and time spectra").getPad(1).getAxisZ().setLog(getLogZ());
        this.getDetectorCanvas().getCanvas("ADC and time spectra").draw(this.getDataGroup().getItem(0,0,0).getH2F("tdc"));
        this.getDetectorCanvas().getCanvas("ADC and time spectra").update();
    }

    @Override
    public void processEvent(DataEvent event) {

        if (this.getNumberOfEvents() >= super.eventResetTime_current[10] && super.eventResetTime_current[10] > 0){
            resetEventListener();
        }
        
        // process event info and save into data group
        
        if(event.hasBank("FTTRK::adc")==true){
	    DataBank bank = event.getBank("FTTRK::adc");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector  = bank.getByte("sector", loop);
                int layer   = bank.getByte("layer", loop);
                int comp    = bank.getInt("component", loop);
                int order   = bank.getByte("order", loop);
                int adc     = bank.getInt("ADC", loop);
                float time  = bank.getFloat("time", loop);
                int channel = comp + (layer-1)*nstrip;
                int dream   = ((int) comp/this.numberOfStripsPerChip) + 1;
                        
//                System.out.println("ROW " + loop + " SECTOR = " + sector + " LAYER = " + layer + " COMPONENT = " + comp +
//                      " ADC = " + adc); 
                if(adc>0) {
                    
                        this.getDataGroup().getItem(0,0,0).getH2F("occADC_2D").fill(comp*1.0, layer*1.0);
                        
                        this.getDataGroup().getItem(0,layer,0).getH1F("occADC_layer" + layer).fill(comp*1.0);
                        this.getDataGroup().getItem(0,layer,0).getH1F("timeMaxTmp1_layer" + layer).fill(dream*1.0,1.0);
                        this.getDataGroup().getItem(0,layer,0).getH1F("timeMaxTmp2_layer" + layer).fill(dream*1.0,time);
                        
                        this.getDataGroup().getItem(0,0,0).getH2F("adc").fill(adc*1.0,channel*1.0);
                        this.getDataGroup().getItem(0,0,0).getH2F("tdc").fill(time*1.0,channel*1.0);
                }
                this.getDetectorSummary().getH1F("summary").fill(channel*1.0);
	    }
    	}
                
                
    }

    @Override
    public void timerUpdate() {
        if(this.getNumberOfEvents()>0) {
            for(int layer=1; layer <=nlayer; layer++) {
                H1F raw1 = this.getDataGroup().getItem(0,layer,0).getH1F("timeMaxTmp1_layer"+layer);
                H1F raw2 = this.getDataGroup().getItem(0,layer,0).getH1F("timeMaxTmp2_layer"+layer);
                H1F ave = this.getDataGroup().getItem(0,layer,0).getH1F("timeMax_layer"+layer);
                for(int loop = 0; loop < raw1.getDataSize(0); loop++){
                    ave.setBinContent(loop, raw2.getBinContent(loop)/raw1.getBinContent(loop));
                }
            }
        }

    }


}