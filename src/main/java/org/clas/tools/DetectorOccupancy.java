package org.clas.tools;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jlab.detector.base.DetectorCollection;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

/**
 *
 * @author gavalian
 */
public class DetectorOccupancy {
    
    DetectorCollection<DetectorMeasurement>  occupancyCollection = 
            new DetectorCollection<DetectorMeasurement>();
    
    DetectorCollection<DetectorMeasurement>  valueCollection = 
            new DetectorCollection<DetectorMeasurement>();
    
    private int maxLayers     = 9;
    private int maxComponents = 72;
    public  int[] ADCWindow = {20,300};
    public  int[] TDCWindow = {200,300};
    ECConstants ecc = new ECConstants();
    
    public DetectorOccupancy(){
        
    }
    
    public DetectorOccupancy(int max_layers, int max_components){
        this.maxLayers     = max_layers;
        this.maxComponents = max_components;
    }
    
    public void addADCBank(DataBank bank){
        int nrows = bank.rows();
        for(int row = 0; row < nrows; row++){
            int    sector = bank.getByte(  "sector",    row);
            int     layer = bank.getByte(  "layer",     row);
            int component = bank.getShort( "component", row);
            float sca = (float) ((sector==5)?ecc.SCALE5[layer-1]:ecc.SCALE[layer-1]);
            int       adc = (int)(bank.getInt("ADC", row)/sca);
            

            if(adc>this.ADCWindow[0] && adc<this.ADCWindow[1]){
                if(occupancyCollection.hasEntry(sector, layer, component)==true){
                    this.occupancyCollection.get(sector, layer, component).incrementADC(adc);                
                } else {
                    DetectorMeasurement measure = new DetectorMeasurement();
                    measure.incrementADC(adc);
                    this.occupancyCollection.add(sector, layer, component, measure);
                }
            }
            if(valueCollection.hasEntry(sector, layer, component)==true){
               this.valueCollection.get(sector, layer, component).storeADCV(adc);
            } else {
               DetectorMeasurement measure = new DetectorMeasurement();
               measure.storeADCV(adc);
               this.valueCollection.add(sector, layer, component, measure);
            }                            
        }
    }
    
    public void addTDCBank(DataBank bank){
        int nrows = bank.rows();
        for(int row = 0; row < nrows; row++){
            int    sector = bank.getByte(  "sector",    row);
            int     layer = bank.getByte(  "layer",     row);
            int component = bank.getShort( "component", row);
            int       tdc = (int)(bank.getInt("TDC", row)*0.02345f);
            if(tdc>this.TDCWindow[0] && tdc<this.TDCWindow[1]){
            	if(occupancyCollection.hasEntry(sector, layer, component)==true){
            		this.occupancyCollection.get(sector, layer, component).incrementTDC(tdc);                
            	} else {
            		DetectorMeasurement measure = new DetectorMeasurement();
            		measure.incrementTDC(tdc);
            		this.occupancyCollection.add(sector, layer, component, measure);
            	}
            }
            if(valueCollection.hasEntry(sector, layer, component)==true){
                this.valueCollection.get(sector, layer, component).storeTDCV(tdc);
            } else {
                DetectorMeasurement measure = new DetectorMeasurement();
                measure.storeTDCV(tdc);
                this.valueCollection.add(sector, layer, component, measure);
            }
        }
    }
    
    public int getTDC(int sector, int layer, int component){
        if(this.occupancyCollection.hasEntry(sector, layer, component)==false)
            return 0;
        return this.occupancyCollection.get(sector, layer, component).TDCCount;
    }
    
    public int getADC(int sector, int layer, int component){
        if(this.occupancyCollection.hasEntry(sector, layer, component)==false)
            return 0;
        return this.occupancyCollection.get(sector, layer, component).ADCCount;
    }
    
    public int getTDCV(int sector, int layer, int component){
        if(this.valueCollection.hasEntry(sector, layer, component)==false)
            return 0;
        return this.valueCollection.get(sector, layer, component).TDCValue;
    }
    
    public int getADCV(int sector, int layer, int component){
        if(this.valueCollection.hasEntry(sector, layer, component)==false)
            return 0;
        return this.valueCollection.get(sector, layer, component).ADCValue;
    }
    
    public int getTDCVC(int sector, int layer, int component){
        if(this.occupancyCollection.hasEntry(sector, layer, component)==false)
            return 0;
        return this.occupancyCollection.get(sector, layer, component).TDCVCount;
    }
    
    public int getADCVC(int sector, int layer, int component){
        if(this.occupancyCollection.hasEntry(sector, layer, component)==false)
            return 0;
        return this.occupancyCollection.get(sector, layer, component).ADCVCount;
    }   
    
    public int getMaxADC(){
        int max = 0;
        List<DetectorMeasurement> measures = this.occupancyCollection.getList();
         for(DetectorMeasurement m : measures){
             if(m.ADCCount>max){
                 max = m.ADCCount;
             }
        }
        return max;
    }
    
    public int getMaxTDC(){
        int max = 0;
        List<DetectorMeasurement> measures = this.occupancyCollection.getList();
         for(DetectorMeasurement m : measures){
             if(m.TDCCount>max){
                 max = m.TDCCount;
             }
        }
        return max;
    }
    
    public void reset(){
        List<DetectorMeasurement> measures = this.occupancyCollection.getList();
        for(DetectorMeasurement m : measures){
            m.reset();
        }
    }
    
    public void resetValue(){
        List<DetectorMeasurement> measures = this.valueCollection.getList();
        for(DetectorMeasurement m : measures){
            m.reset();
        }
    }        
    
    public GraphErrors  getOccupancyGraph(){
        int maxADC = this.getMaxADC();
        Set<Long> keySet = this.occupancyCollection.getKeys();
        GraphErrors graph = new GraphErrors();
        graph.setMarkerSize(0); 
        graph.setLineColor(4);
        
        Set<Integer>  sectors = this.getCollection().getSectors();
        
        for(Integer sector : sectors){
            Set<Integer> layers = this.getCollection().getLayers(sector);
            for(Integer layer : layers){
                Set<Integer> components = this.getCollection().getComponents(sector, layer);
                for(Integer component : components){
                    double x = (sector) * 1.0 + (layer-1)*1.0/maxLayers + (component-1)*1.0/maxLayers/maxComponents;
                    DetectorMeasurement measure = this.occupancyCollection.get(sector, layer, component);
                    double intencity = ((double) measure.ADCCount)/maxADC;
                    graph.addPoint((double) x, 0.0, 0.0, intencity);
                }
            }
        }        
        return graph;
    }
    
    public GraphErrors  getOccupancyGraphTDC(){
        int maxTDC = this.getMaxTDC();
        System.out.println(" TDC max = " + maxTDC);
        Set<Long> keySet = this.occupancyCollection.getKeys();
        GraphErrors graph = new GraphErrors();
        graph.setMarkerSize(0); 
        graph.setLineColor(4);
        
        Set<Integer>  sectors = this.getCollection().getSectors();
        
        for(Integer sector : sectors){
            Set<Integer> layers = this.getCollection().getLayers(sector);
            for(Integer layer : layers){
                Set<Integer> components = this.getCollection().getComponents(sector, layer);
                for(Integer component : components){
                    double x = (sector) * 1.0 + (layer-1)*1.0/maxLayers + (component-1)*1.0/maxLayers/maxComponents;
                    DetectorMeasurement measure = this.occupancyCollection.get(sector, layer, component);
                    double intencity = ((double) measure.TDCCount)/maxTDC;
                    graph.addPoint((double) x, 0.0, 0.0, intencity);
                    System.out.println(" MAX = " + maxTDC + " " + sector + " " + layer + " " + component + " " + measure.TDCCount);
                }
            }
        }        
        return graph;
    }
    
    public DetectorCollection getCollection(){ return this.occupancyCollection;}
    
    public DetectorCollection getValueCollection(){ return this.valueCollection;}
    
    public static class DetectorMeasurement {
        
        int  ADCCount = 0;
        int  TDCCount = 0;
        int  ADCValue = 0;
        int  TDCValue = 0;
        int ADCVCount = 0;
        int TDCVCount = 0;
        
        public DetectorMeasurement(){
            
        }
        
        public void reset(){
            this.ADCCount = 0;
            this.TDCCount = 0;
            this.ADCValue = 0;
            this.TDCValue = 0;
            this.ADCVCount = 0;
            this.TDCVCount = 0;
        }
        
        public void incrementADC(int val){
            this.ADCCount++;
            this.ADCVCount+=val;
        }
        
        public void incrementTDC(int val){
            this.TDCCount++;
            this.TDCVCount+=val;
        }
        
        public void storeADCV(int val){
            this.ADCValue = val;
        }
        
        public void storeTDCV(int val){
            this.TDCValue = val;
        }
        
    }
    
    public static void main(String[] args){
        HipoDataSource reader = new HipoDataSource();
        reader.open("/Users/gavalian/Work/Software/Release-9.0/COATJAVA/ecreconstruction/gemc_dis.hipo");
        DetectorOccupancy occupancyECAL = new DetectorOccupancy();
        DetectorOccupancy occupancyFTOF = new DetectorOccupancy(3,63);
        DetectorOccupancy occupancyDC   = new DetectorOccupancy(36,112);
        
        while(reader.hasEvent()==true){
            DataEvent event = reader.getNextEvent();
            if(event.hasBank("ECAL::adc")==true){
                DataBank bank = event.getBank("ECAL::adc");
                occupancyECAL.addADCBank(bank);
            }
            if(event.hasBank("FTOF::adc")==true){
                DataBank bank = event.getBank("FTOF::adc");
                occupancyFTOF.addADCBank(bank);
            }
            if(event.hasBank("DC::tdc")==true){
                DataBank bank = event.getBank("DC::tdc");
                occupancyDC.addTDCBank(bank);
            }
        }
        
        GraphErrors graphECAL = occupancyECAL.getOccupancyGraph();
        GraphErrors graphFTOF = occupancyFTOF.getOccupancyGraph();
        GraphErrors graphDC   = occupancyDC.getOccupancyGraphTDC();
        TCanvas c1 = new TCanvas("c1",800,300);
        c1.divide(1, 3);
        c1.cd(0);
        graphECAL.setTitle("ECAL");
        c1.draw(graphECAL);
        c1.getCanvas().getPad(0).setAxisRange(1.0, 7.0, -1.2, 1.2);
        c1.cd(1);
        graphFTOF.setTitle("FTOF");
        c1.draw(graphFTOF);
        c1.getCanvas().getPad(1).setAxisRange(1.0, 7.0, -1.2, 1.2);
        c1.cd(2);
        graphDC.setTitle("DC");
        c1.draw(graphDC);
        c1.getCanvas().getPad(2).setAxisRange(1.0, 7.0, -1.2, 1.2);
    }
}
