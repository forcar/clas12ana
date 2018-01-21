package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */

public class TRKmonitor extends DetectorMonitor {
    
    // reconstruction engines
//    DCHBEngine enHB = new DCHBEngine();
//	
//    DCTBEngine enTB = new DCTBEngine();

 
    public TRKmonitor(String name) {
        super(name);
        
        // initializa monitoring panel
        this.setDetectorTabNames("Tracking Statistics");
        this.init(true);
        
        // initialize reconstruction engines
//        enHB.init();
//        enTB.init();

    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").divide(2, 4);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").setGridX(false);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").setGridY(false);
        
        H1F sumHBT = new H1F("sumHBT","sumHBT",6,0.5,6.5);
        sumHBT.setTitleX("sector");
        sumHBT.setTitleY("HB tracks");
        sumHBT.setFillColor(34);
        H1F sumTBT = new H1F("sumTBT","sumTBT",6,0.5,6.5);
        sumTBT.setTitleX("sector");
        sumTBT.setTitleY("TB tracks");
        sumTBT.setFillColor(32);
        DataGroup summary = new DataGroup(2,1);
        summary.addDataSet(sumHBT, 0);
        summary.addDataSet(sumTBT, 1);
        this.setDetectorSummary(summary);

        H2F clusterSizeHBT      = new H2F("clusterSizeHBT", "clusterSizeHBT", 15, 3, 18, 6, 0.5, 6.5);
        clusterSizeHBT.setTitleX("HBT cluster size");
        clusterSizeHBT.setTitleY("sector");
        H2F numberOfHitsHBT     = new H2F("numberOfHitsHBT", "numberOfHitsHBT", 100, 0, 200, 6, 0.5, 6.5);
        numberOfHitsHBT.setTitleX("HBT hits");
        numberOfHitsHBT.setTitleY("sector");
        H2F numberOfClustersHBT = new H2F("numberOfClustersHBT", "numberOfClustersHBT", 60, 0, 60, 6, 0.5, 6.5);
        numberOfClustersHBT.setTitleX("HBT clusters");
        numberOfClustersHBT.setTitleY("sector");
        H2F numberOfCrossesHBT  = new H2F("numberOfCrossesHBT", "numberOfCrossesHBT", 20, 0, 20, 6, 0.5, 6.5);
        numberOfCrossesHBT.setTitleX("HBT crosses");
        numberOfCrossesHBT.setTitleY("sector");
        H2F clusterSizeTBT      = new H2F("clusterSizeTBT", "clusterSizeTBT", 15, 3, 18, 6, 0.5, 6.5);
        clusterSizeTBT.setTitleX("TBT cluster size");
        clusterSizeTBT.setTitleY("sector");
        H2F numberOfHitsTBT     = new H2F("numberOfHitsTBT", "numberOfHitsTBT", 100, 0, 200, 6, 0.5, 6.5);
        numberOfHitsTBT.setTitleX("TBT hits");
        numberOfHitsTBT.setTitleY("sector");
        H2F numberOfClustersTBT = new H2F("numberOfClustersTBT", "numberOfClustersTBT", 60, 0, 60, 6, 0.5, 6.5);
        numberOfClustersTBT.setTitleX("TBT clusters");
        numberOfClustersTBT.setTitleY("sector");
        H2F numberOfCrossesTBT  = new H2F("numberOfCrossesTBT", "numberOfCrossesTBT", 20, 0, 20, 6, 0.5, 6.5);
        numberOfCrossesTBT.setTitleX("TBT crosses");
        numberOfCrossesTBT.setTitleY("sector");

        DataGroup dg = new DataGroup(2,4);
        dg.addDataSet(clusterSizeHBT, 0);
        dg.addDataSet(numberOfHitsHBT, 2);
        dg.addDataSet(numberOfClustersHBT, 4);
        dg.addDataSet(numberOfCrossesHBT, 6);
        dg.addDataSet(clusterSizeTBT, 1);
        dg.addDataSet(numberOfHitsTBT, 3);
        dg.addDataSet(numberOfClustersTBT, 5);
        dg.addDataSet(numberOfCrossesTBT, 7);
        this.getDataGroup().add(dg, 0,0,0);
    }

    @Override
    public void drawDetector() {
        // Load the Constants
//        Constants.Load(true, true, 0);
//        GeometryLoader.Load(10, "default");

        for(int s =0; s< 6; s++)
            for(int slrnum = 5; slrnum > -1; slrnum--) {
                //DetectorShape2D module = new DetectorShape2D(DetectorType.DC,s,slrnum,0);

                 DetectorShape2D module = new DetectorShape2D();
                 module.getDescriptor().setType(DetectorType.DC);
                 module.getDescriptor().setSectorLayerComponent((s+1), (slrnum+1), 1);

//                module.getShapePath().addPoint(GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(0).getLine().origin().x(),  
//                        -GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(0).getLine().origin().y(),  0.0);
//                module.getShapePath().addPoint(GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(0).getLine().end().x(),  
//                        -GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(0).getLine().end().y(),  0.0);
//                module.getShapePath().addPoint(GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(111).getLine().end().x(),  
//                        -GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(111).getLine().end().y(),  0.0);
//                module.getShapePath().addPoint(GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(111).getLine().origin().x(),  
//                        -GeometryLoader.dcDetector.getSector(0).getSuperlayer(slrnum).getLayer(0).getComponent(111).getLine().origin().y(),  0.0);


                if(slrnum%2==1)
                        module.setColor(180-slrnum*15,180,255);
                if(slrnum%2==0)
                        module.setColor(255-slrnum*15,182,229, 200);

                module.getShapePath().translateXYZ(110.0+((int)(slrnum/2))*50, 0, 0);
                module.getShapePath().rotateZ(s*Math.toRadians(60.));

                this.getDetectorView().getView().addShape("DC",module);			

            }
        this.getDetectorView().setName("DC"); 
        //detectorViewDC.updateBox();
    }

    public void PlotCrosses(DataEvent event){
    		
        if(event.hasBank("HitBasedTrkg::HBCrosses")==true){
            DataBank bank = event.getBank("HitBasedTrkg::HBCrosses");
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
            	int sector = bank.getByte("sector", loop);
            	int region = bank.getByte("region", loop);
            	
            	float x = bank.getFloat("x", loop);
            	float y = bank.getFloat("y", loop);
            	
            	Point3D hit = new Point3D(x, -y, 0);
            	hit.translateXYZ(110.0+((int)((region-1)/2))*50, 0, 0);
                hit.rotateZ((sector-1)*Math.toRadians(60.));

                DetectorShape2D module = new DetectorShape2D(DetectorType.DC,sector,region,0);
                System.out.println(" Point "+hit.toString());
                module.getShapePath().addPoint( hit.x(), hit.y(),0 );
                this.getDetectorView().getView().addShape("region", module);	            	
            }
        } 
       
        if(event.hasBank("TimeBasedTrkg::TBCrosses")==true){
            DataBank bank = event.getBank("TimeBasedTrkg::TBCrosses");
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
            	int sector = bank.getByte("sector", loop);
            	int region = bank.getByte("region", loop);
            	
            	float x = bank.getFloat("x", loop);
            	float y = bank.getFloat("y", loop);
            	
            	Point3D hit = new Point3D(x, -y, 0);
            	hit.translateXYZ(110.0+((int)((region-1)/2))*50, 0, 0);
                hit.rotateZ((sector-1)*Math.toRadians(60.));
            	 
            	//detectorViewDC2.addHit(hit.x(), hit.y(), 0, 150,50,50);
            	//detectorViewDC2.repaint();	
            	
            }
        } 
        
    }    
    
    @Override
    public void plotHistos() {
        // initialize canvas and plot histograms
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(0);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("clusterSizeHBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(2);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("numberOfHitsHBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(4);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("numberOfClustersHBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(6);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("numberOfCrossesHBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(1);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("clusterSizeTBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(3);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("numberOfHitsTBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(5);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("numberOfClustersTBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").cd(7);
        this.getDetectorCanvas().getCanvas("Tracking Statistics").draw(this.getDataGroup().getItem(0,0,0).getH2F("numberOfCrossesTBT"));
        this.getDetectorCanvas().getCanvas("Tracking Statistics").update();        
    }
    
    @Override
    public void plotEvent(DataEvent event) {
        this.PlotCrosses(event);
        this.getDetectorView().update();

    }
    
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        // create Header bank
//        enHead.processDataEvent(eventd);
        // run HBT
//        enHB.processDataEvent(eventd);		
        // Processing TBT
//        enTB.processDataEvent(eventd);
        
        if(event.hasBank("HitBasedTrkg::HBTracks")==true){
	    DataBank bank = event.getBank("HitBasedTrkg::HBTracks");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector =bank.getByte("sector",loop);
                this.getDetectorSummary().getH1F("sumHBT").fill(sector*1.0);
	    }
        }
        if(event.hasBank("TimeBasedTrkg::TBTracks")==true){
	    DataBank bank = event.getBank("TimeBasedTrkg::TBTracks");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector =bank.getByte("sector",loop);
                this.getDetectorSummary().getH1F("sumTBT").fill(sector*1.0);
	    }
        }
        
        if(event.hasBank("HitBasedTrkg::HBCrosses")==true){
            int[] nHits     = new int[]{0,0,0,0,0,0};
            int[] nClusters = new int[]{0,0,0,0,0,0};
            int[] nCrosses  = new int[]{0,0,0,0,0,0};
	    DataBank bank = event.getBank("HitBasedTrkg::HBHits");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector =bank.getByte("sector",loop);
                nHits[sector-1]++;
	    }
            for(int sector=1; sector<=6; sector++) if(nHits[sector-1]>0) this.getDataGroup().getItem(0,0,0).getH2F("numberOfHitsHBT").fill(nHits[sector-1],sector);
            bank = event.getBank("HitBasedTrkg::HBClusters");
	    rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector = bank.getByte("sector",loop);
                int size   = bank.getByte("size",loop);
                nClusters[sector-1]++;
                this.getDataGroup().getItem(0,0,0).getH2F("clusterSizeHBT").fill(size,sector);
	    }
            for(int sector=1; sector<=6; sector++) if(nClusters[sector-1]>0) this.getDataGroup().getItem(0,0,0).getH2F("numberOfClustersHBT").fill(nClusters[sector-1],sector);
	    bank = event.getBank("HitBasedTrkg::HBCrosses");
	    rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector =bank.getByte("sector",loop);
                nCrosses[sector-1]++;
	    }
            for(int sector=1; sector<=6; sector++) if(nCrosses[sector-1]>0) this.getDataGroup().getItem(0,0,0).getH2F("numberOfCrossesHBT").fill(nCrosses[sector-1],sector);
	}
        if(event.hasBank("TimeBasedTrkg::TBCrosses")==true){
            int[] nHits     = new int[]{0,0,0,0,0,0};
            int[] nClusters = new int[]{0,0,0,0,0,0};
            int[] nCrosses  = new int[]{0,0,0,0,0,0};
	    DataBank bank = event.getBank("TimeBasedTrkg::TBHits");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector =bank.getByte("sector",loop);
                nHits[sector-1]++;
	    }
            for(int sector=1; sector<=6; sector++) if(nHits[sector-1]>0) this.getDataGroup().getItem(0,0,0).getH2F("numberOfHitsTBT").fill(nHits[sector-1],sector);
            bank = event.getBank("TimeBasedTrkg::TBClusters");
	    rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector = bank.getByte("sector",loop);
                int size   = bank.getByte("size",loop);
                nClusters[sector-1]++;
                this.getDataGroup().getItem(0,0,0).getH2F("clusterSizeTBT").fill(size,sector);
	    }
            for(int sector=1; sector<=6; sector++) if(nClusters[sector-1]>0)  this.getDataGroup().getItem(0,0,0).getH2F("numberOfClustersTBT").fill(nClusters[sector-1],sector);
	    bank = event.getBank("TimeBasedTrkg::TBCrosses");
	    rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int sector =bank.getByte("sector",loop);
                nCrosses[sector-1]++;
	    }
            for(int sector=1; sector<=6; sector++) if(nCrosses[sector-1]>0) this.getDataGroup().getItem(0,0,0).getH2F("numberOfCrossesTBT").fill(nCrosses[sector-1],sector);
	}       
    }

    @Override
    public void timerUpdate() {

    }

 

}
