package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

public class BMTmonitor extends DetectorMonitor {
	int numberOfSamples;
	int timeSampling;

	int maxNumberLayer;
	int maxNumberSector;
	int maxNumberStrips;
	int numberOfStripsPerChip;
	int numberOfChips;
	int binDivision;
	int numberStrips[];
        int numberDream[];

	double hitNumber[][][];
	int hitNumberDream [];
	int timeMax[][][][];
	int layerHit[][];
	boolean isZ[];
	boolean mask[][][];
	
	public BMTmonitor(String name) {
		super(name);
		
		numberOfSamples = 16;
		timeSampling = 40;
		maxNumberLayer = 6;
		maxNumberSector = 3;
		maxNumberStrips = 1152;
		numberOfStripsPerChip = 64 ;
                
		numberOfChips = 240;
		isZ = new boolean[maxNumberLayer +1];
		mask = new boolean[maxNumberSector + 1][maxNumberLayer + 1][maxNumberStrips + 1];
		
                hitNumberDream = new int [numberOfChips+1];
		numberStrips = new int[maxNumberLayer + 1];
                numberDream = new int[maxNumberLayer + 1];
		
		numberStrips[1] = 896;
		numberStrips[2] = 640;
		numberStrips[3] = 640;
		numberStrips[4] = 1024;
		numberStrips[5] = 768;
		numberStrips[6] = 1152;
		
                numberDream[1] = 14;
                numberDream[2] = 10;
                numberDream[3] = 10;
                numberDream[4] = 16;
                numberDream[5] = 12;
                numberDream[6] = 18; 
                
		isZ[1] = false;
		isZ[2] = true;
		isZ[3] = true;
		isZ[4] = false;
		isZ[5] = true;
		isZ[6] = false;
		
		for (int sector = 1; sector <= maxNumberSector; sector++) {
			for (int layer = 1; layer <= maxNumberLayer; layer++) {
				for (int component = 1 ; component <= numberStrips[layer]; component++){
					mask[sector][layer][component]=true;
				}
			}
		}
		mask[1][1][80] = false;
		mask[1][1][82] = false;
		
		mask[2][1][179] = false;
		mask[2][1][192] = false;
		mask[2][1][544] = false;
		mask[2][1][545] = false;
		
		mask[3][4][757] = false;
		
		mask[3][5][102] = false; 
		mask[3][5][735] = false;
		
		this.setDetectorTabNames("Occupancies", "Occupancy C", "Occupancy Z", "TimeOfMax", "Multiplicity");
		this.init(false);
	}

	@Override
	public void createHistos() {
	
		// create histograms
		this.setNumberOfEvents(0);

		H2F summary = new H2F("summary","summary",maxNumberStrips, 0, maxNumberStrips, maxNumberLayer*maxNumberSector,0,maxNumberLayer*maxNumberSector);
		summary.setTitleX("strips");
		summary.setTitleY("detector");
		summary.setTitle("BMT");
		DataGroup sum = new DataGroup(1,1);
		sum.addDataSet(summary, 0);
		this.setDetectorSummary(sum);
	
		H2F occupancyHisto = new H2F("Occupancies","Occupancies",maxNumberStrips, 0, maxNumberStrips, maxNumberLayer*maxNumberSector,0,maxNumberLayer*maxNumberSector);
		occupancyHisto.setTitleX("Strips");
		occupancyHisto.setTitleY("Detector");
                H1F histmulti = new H1F("multi", "multi", 200, -0.5, 199.5);
                histmulti.setTitleX("hit multiplicity");
                histmulti.setTitleY("counts");
                histmulti.setTitle("Multiplicity of BMT channels"); 
		DataGroup occupancyGroup = new DataGroup("");
		occupancyGroup.addDataSet(occupancyHisto, 0);
                occupancyGroup.addDataSet(histmulti, 0);
		this.getDataGroup().add(occupancyGroup, 0, 0, 0);
		
		H1F timeOfMaxHisto = new H1F("TimeOfMax","TimeOfMax",numberOfChips,0,numberOfChips);
		timeOfMaxHisto.setTitleX("Electronic chip");
		timeOfMaxHisto.setTitleY("Time of max adc");
		timeOfMaxHisto.setFillColor(4);
		DataGroup timeOfMaxGroup = new DataGroup("");
		timeOfMaxGroup.addDataSet(timeOfMaxHisto, 0);
		this.getDataGroup().add(timeOfMaxGroup, 0, 0, 1);
		
		for (int sector = 1; sector <= maxNumberSector; sector++) {
			for (int layer = 1; layer <= maxNumberLayer; layer++) {
				H1F hitmapHisto = new H1F("Occupancy Layer " + layer + " Sector " + sector, "Occupancy Layer " + layer + " Sector " + sector,
						(numberStrips[layer])+1, 0., (double) (numberStrips[layer])+1);
				hitmapHisto.setTitleX("Strips (Layer " + layer  + " Sector " + sector+")");
				hitmapHisto.setTitleY("Nb of hits");
				if (isZ[layer]){
					hitmapHisto.setFillColor(4);
				}else{
					hitmapHisto.setFillColor(8);
				}
				DataGroup hitmapGroup = new DataGroup("");
				hitmapGroup.addDataSet(hitmapHisto, 0);
				this.getDataGroup().add(hitmapGroup, sector, layer,2);
			}
		}
                        
	}

	@Override
	public void plotHistos() {
		
		this.getDetectorCanvas().getCanvas("Occupancies").setGridX(false);
		this.getDetectorCanvas().getCanvas("Occupancies").setGridY(false);
		this.getDetectorCanvas().getCanvas("Occupancies").setAxisTitleSize(12);
		this.getDetectorCanvas().getCanvas("Occupancies").setAxisLabelSize(12);
                this.getDetectorCanvas().getCanvas("Occupancies").getPad(0).getAxisZ().setLog(getLogZ());
		this.getDetectorCanvas().getCanvas("Occupancies").draw(this.getDataGroup().getItem(0, 0, 0).getH2F("Occupancies"));
		this.getDetectorCanvas().getCanvas("Occupancies").update();

		this.getDetectorCanvas().getCanvas("TimeOfMax").setGridX(false);
		this.getDetectorCanvas().getCanvas("TimeOfMax").setGridY(false);
		this.getDetectorCanvas().getCanvas("TimeOfMax").setAxisTitleSize(12);
		this.getDetectorCanvas().getCanvas("TimeOfMax").setAxisLabelSize(12);
		this.getDetectorCanvas().getCanvas("TimeOfMax").draw(this.getDataGroup().getItem(0, 0, 1).getH1F("TimeOfMax"));
		this.getDetectorCanvas().getCanvas("TimeOfMax").update();

		this.getDetectorCanvas().getCanvas("Occupancy C").divide(maxNumberSector, maxNumberLayer/2);
		this.getDetectorCanvas().getCanvas("Occupancy C").setGridX(false);
		this.getDetectorCanvas().getCanvas("Occupancy C").setGridY(false);
		this.getDetectorCanvas().getCanvas("Occupancy C").setAxisTitleSize(12);
		this.getDetectorCanvas().getCanvas("Occupancy C").setAxisLabelSize(12);
		
		this.getDetectorCanvas().getCanvas("Occupancy Z").divide(maxNumberSector, maxNumberLayer/2);
		this.getDetectorCanvas().getCanvas("Occupancy Z").setGridX(false);
		this.getDetectorCanvas().getCanvas("Occupancy Z").setGridY(false);
		this.getDetectorCanvas().getCanvas("Occupancy Z").setAxisTitleSize(12);
		this.getDetectorCanvas().getCanvas("Occupancy Z").setAxisLabelSize(12);
		
		for (int sector = 1; sector <= maxNumberSector; sector++) {
			for (int layer = 1; layer <= maxNumberLayer; layer++) {
				int column=maxNumberSector - sector;
				int row;
				int numberOfColumns=maxNumberSector;
				switch (layer) {
				case 1: row=2; break;
				case 2: row=2; break;
				case 3: row=1; break;
				case 4: row=1; break;
				case 5: row=0; break;
				case 6: row=0; break;
				default:row=-1;break;
				}
				if (isZ[layer]){
					this.getDetectorCanvas().getCanvas("Occupancy Z").cd(column + numberOfColumns * row);
					this.getDetectorCanvas().getCanvas("Occupancy Z").draw(
							this.getDataGroup().getItem(sector, layer, 2).getH1F("Occupancy Layer " + layer + " Sector " + sector));
				}else{
					this.getDetectorCanvas().getCanvas("Occupancy C").cd(column + numberOfColumns * row);
					this.getDetectorCanvas().getCanvas("Occupancy C").draw(
							this.getDataGroup().getItem(sector, layer, 2).getH1F("Occupancy Layer " + layer + " Sector " + sector));
				}
			}
		}
		this.getDetectorCanvas().getCanvas("Occupancy Z").update();
		this.getDetectorCanvas().getCanvas("Occupancy C").update();
                
                this.getDetectorCanvas().getCanvas("Multiplicity").divide(1, 1);
                this.getDetectorCanvas().getCanvas("Multiplicity").setGridX(false);
                this.getDetectorCanvas().getCanvas("Multiplicity").setGridY(false);
                this.getDetectorCanvas().getCanvas("Multiplicity").cd(0);
                this.getDetectorCanvas().getCanvas("Multiplicity").draw(this.getDataGroup().getItem(0,0,0).getH1F("multi"));
                this.getDetectorCanvas().getCanvas("Multiplicity").update();
                
	}

	public void processEvent(DataEvent event) {
            
		if (this.getNumberOfEvents() >= super.eventResetTime_current[0] && super.eventResetTime_current[0] > 0){
		    resetEventListener();
		}
        
		//if (!testTriggerMask()) return;
            
		if (event.hasBank("BMT::adc") == true) {
			DataBank bank = event.getBank("BMT::adc");
                        
                        this.getDataGroup().getItem(0,0,0).getH1F("multi").fill(bank.rows());
                        
			for (int i = 0; i < bank.rows(); i++) {
				int sector = bank.getByte("sector", i);
				int layer = bank.getByte("layer", i);
				int strip = bank.getShort("component", i);
				float timeOfMax = bank.getFloat("time", i);
                                
				if (strip < 0 || !mask[sector][layer][strip]){
					continue;
				}
                                int dream=0;
                                int dreamLayer=0;
                                for (int layerNb=1; layerNb<layer; layerNb++){
                                    dreamLayer = dreamLayer + maxNumberSector * numberDream[layerNb];
                                }
                                int dreamSector = (sector-1) * numberDream[layer];
                                int dreamTile = (strip-1) / numberOfStripsPerChip+1;
                                dream = dreamLayer + dreamSector + dreamTile;
                                hitNumberDream[dream]++;		
				this.getDataGroup().getItem(sector, layer, 2).getH1F("Occupancy Layer " + layer + " Sector " + sector)
				.fill(strip);
				this.getDataGroup().getItem(0, 0, 0).getH2F("Occupancies").fill(strip,3*(layer-1)+(sector-1),1);
                                this.getDetectorSummary().getH2F("summary").fill(strip,3*(layer-1)+(sector-1),1);
		
                                double timeMaxDreamAvg = this.getDataGroup().getItem(0, 0, 1).getH1F("TimeOfMax").getBinContent(dream);
                                this.getDataGroup().getItem(0, 0, 1).getH1F("TimeOfMax").setBinContent(dream,timeMaxDreamAvg + (timeOfMax - timeMaxDreamAvg)/hitNumberDream[dream]);
			}
		}
	}
}
