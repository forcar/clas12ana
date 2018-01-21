package org.clas.detectors;

import java.util.Arrays;
import org.clas.viewer.DetectorMonitor;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;

public class FTHODOmonitor extends DetectorMonitor {

    ConstantsManager ccdb = new ConstantsManager();
    IndexedTable geometry = null;
    double tileWidth = 15;

    public FTHODOmonitor(String name) {
        super(name);

        ccdb.setVariation("default");
        ccdb.init(Arrays.asList(new String[]{
            "/geometry/ft/fthodo"}));
        geometry = ccdb.getConstants(11, "/geometry/ft/fthodo");

        this.setDetectorTabNames("FADC Occupancies", "FADC spectra");
        this.init(false);
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.geometry.show();
        this.setNumberOfEvents(0);

        H1F summary = new H1F("summary", "summary", 232, 0.5, 232.5);
        summary.setTitleX("PMT");
        summary.setTitleY("FTHODO hits");
        summary.setTitle("FTHODO");
        summary.setFillColor(38);
        DataGroup sum = new DataGroup(1, 1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);

        for (int ilayer = 1; ilayer <= 2; ilayer++) {
            H2F occFADC2D = new H2F("occFADC_2D_l" + ilayer, "occFADC_2D_l" + ilayer, 24, -12 * tileWidth, 12 * tileWidth, 24, -12 * tileWidth, 12 * tileWidth);
            occFADC2D.setTitleX("Tile -X");
            occFADC2D.setTitleY("Tile Y");
            occFADC2D.setTitle("layer " + ilayer);

            H1F occFADC = new H1F("occFADC_l" + ilayer, "occFADC_l" + ilayer, 116, 0.5, 116.5);
            occFADC.setTitleX("SiPM");
            occFADC.setTitleY("Counts");
            occFADC.setFillColor(38);
            occFADC.setTitle("layer " + ilayer);

            H2F fadc = new H2F("fadc_l" + ilayer, "fadc_l" + ilayer, 65, 0, 65000, 116, 0.5, 116.5);
            fadc.setTitleX("FADC - amplitude");
            fadc.setTitleY("SiPM");
            fadc.setTitle("layer " + ilayer);

            H2F fadc_time = new H2F("fadc_time_l" + ilayer, "fadc_time_l" + ilayer, 50, 0, 500, 116, 0.5, 116.5);
            fadc_time.setTitleX("FADC - time");
            fadc_time.setTitleY("SiPM");
            fadc_time.setTitle("layer " + ilayer);

            DataGroup dg = new DataGroup(2, 4);
            dg.addDataSet(occFADC2D, 0);
            dg.addDataSet(occFADC, 1);
            dg.addDataSet(fadc, 2);
            dg.addDataSet(fadc_time, 3);

            this.getDataGroup().add(dg, 0, ilayer, 0);
        }
    }

    @Override
    public void plotHistos() {
        // plotting histos
        this.getDetectorCanvas().getCanvas("FADC Occupancies").divide(2, 2);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").setGridX(false);
        this.getDetectorCanvas().getCanvas("FADC Occupancies").setGridY(false);
        for (int ilayer = 1; ilayer <= 2; ilayer++) {
            this.getDetectorCanvas().getCanvas("FADC Occupancies").cd(0 + (ilayer - 1) * 2);
//            this.getDetectorCanvas().getCanvas("FADC Occupancies").getPad(0 + (ilayer - 1) * 2).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("FADC Occupancies").draw(this.getDataGroup().getItem(0, ilayer, 0).getH2F("occFADC_2D_l" + ilayer));
            this.getDetectorCanvas().getCanvas("FADC Occupancies").cd(1 + (ilayer - 1) * 2);
            this.getDetectorCanvas().getCanvas("FADC Occupancies").draw(this.getDataGroup().getItem(0, ilayer, 0).getH1F("occFADC_l" + ilayer));
        }
        this.getDetectorCanvas().getCanvas("FADC spectra").divide(2, 2);
        this.getDetectorCanvas().getCanvas("FADC spectra").setGridX(false);
        this.getDetectorCanvas().getCanvas("FADC spectra").setGridY(false);
        for (int ilayer = 1; ilayer <= 2; ilayer++) {
            this.getDetectorCanvas().getCanvas("FADC spectra").cd(0 + (ilayer - 1) * 2);
            this.getDetectorCanvas().getCanvas("FADC spectra").getPad(0 + (ilayer - 1) * 2).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("FADC spectra").draw(this.getDataGroup().getItem(0, ilayer, 0).getH2F("fadc_l" + ilayer));
            this.getDetectorCanvas().getCanvas("FADC spectra").cd(1 + (ilayer - 1) * 2);
            this.getDetectorCanvas().getCanvas("FADC spectra").getPad(1 + (ilayer - 1) * 2).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("FADC spectra").draw(this.getDataGroup().getItem(0, ilayer, 0).getH2F("fadc_time_l" + ilayer));
        }
    }

    @Override
    public void processEvent(DataEvent event) {

        if (this.getNumberOfEvents() >= super.eventResetTime_current[8] && super.eventResetTime_current[8] > 0){
            resetEventListener();
        }

        // process event info and save into data group
        if (event.hasBank("FTHODO::adc") == true) {
            DataBank bank = event.getBank("FTHODO::adc");
            int rows = bank.rows();
            for (int loop = 0; loop < rows; loop++) {
                int sector = bank.getByte("sector", loop);
                int layer = bank.getByte("layer", loop);
                int comp = bank.getShort("component", loop);
                int order = bank.getByte("order", loop);
                int adc = bank.getInt("ADC", loop);
                double time = bank.getFloat("time", loop);
                double x = geometry.getDoubleValue("x", sector, layer, comp);
                double y = geometry.getDoubleValue("y", sector, layer, comp);
                double w = geometry.getDoubleValue("width", sector, layer, comp);
                int tile = this.getTileId(sector, comp);

//             System.out.println("ROW " + loop + " SECTOR = " + sector + " LAYER = " + layer + " COMPONENT = " + comp + " ORDER + " + order +
//             " ADC = " + adc + " TIME = " + time); 
                if (adc > 0) {

                    this.fillTile2D(this.getDataGroup().getItem(0, layer, 0).getH2F("occFADC_2D_l" + layer), -x, y, w);

                    this.getDataGroup().getItem(0, layer, 0).getH1F("occFADC_l" + layer).fill(tile);
                    this.getDataGroup().getItem(0, layer, 0).getH2F("fadc_l" + layer).fill(adc * 1.0, tile);
                    this.getDataGroup().getItem(0, layer, 0).getH2F("fadc_time_l" + layer).fill(time * 1.0, tile);
                    this.getDetectorSummary().getH1F("summary").fill(tile + 116 * (layer - 1));

                }

            }
        }

    }

    @Override
    public void timerUpdate() {

    }

    private void fillTile2D(H2F h2, double x, double y, double width) {
        if (width == tileWidth) {
            h2.fill(x, y);
        } else {
            double delta = tileWidth / 2;
            h2.fill(x - delta, y - delta);
            h2.fill(x - delta, y + delta);
            h2.fill(x + delta, y - delta);
            h2.fill(x + delta, y + delta);
        }
    }

    private int getTileId(int sector, int component) {
        int tile = -1;
        switch (sector) {
            case 1:
                tile = component + 0;
                break;
            case 2:
                tile = component + 9;
                break;
            case 3:
                tile = component + 29;
                break;
            case 4:
                tile = component + 38;
                break;
            case 5:
                tile = component + 58;
                break;
            case 6:
                tile = component + 67;
                break;
            case 7:
                tile = component + 87;
                break;
            case 8:
                tile = component + 96;
                break;
            default:
                tile = -1;
                break;
        }
        return tile;
    }
}