package org.clas.service.ec;

import org.jlab.detector.base.DetectorDescriptor;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;

/**
 *
 * @author gavalian
 */

public class ECStrip implements Comparable {
    
	/*This is the time to consider two nearby strips (same detector view!) in time. 
	 * There's no correction (yet) for hit position along the bar
	 * However, two bars are near-by and, if the hit is really a matching hit, it is basically in the same position
	 * in the two bars. Hence, the time difference should be very small
	 */
	
	private DetectorDescriptor  desc = new DetectorDescriptor(DetectorType.ECAL);
    
    private int                iADC = 0;
    private int                iTDC = 0;
    private float             iTADC = 0;
    private double            iGain = 1.0;
    private double     iADC_to_MEV  = 1.0/10000.0;
    private double    iAttenLengthA = 1.0;
    private double    iAttenLengthB = 50000.0;
    private double    iAttenLengthC = 0.0;
    private double    fTim00,iTim00 = 0; // Global TDC offset
    private double    fTimA0,iTimA0 = 0; // Offset in ns (before applying a1)
    private double    fTimA1,iTimA1 = 1; // ns -> TDC conv. factor (TDC = ns/a1)
    private double    fTimA2,iTimA2 = 0; // time-walk factor (time_ns = time_ns + a2/sqrt(adc))
    private double    fTimA3,iTimA3 = 0; // 0
    private double    fTimA4,iTimA4 = 0; // 0
    private double    fTimA5,iTimA5 = 0; // 0
    private double    fTimA6,iTimA6 = 0; // 0
    private int        triggerPhase = 0;
    private double             veff = 18.1; // Effective velocity of scintillator light (cm/ns)
    private int                  id = -1;   // ID (row number) of the corresponding hit in the ADC bank for truth matching
    private int              peakID = -1;   // ID of peak this hit belongs to (all strips belonging to a peak have this same number)
    private int           clusterId = -1;   // Id (row number) of the cluster that this hit belongs to
    private int             stripId = -1;   // Id (row number) of the peak striplist that this hit belongs to
	
    private double            tdist = 0;
    private double            edist = 0;
    
    private Line3D         stripLine = new Line3D();
    private double stripDistanceEdge = 0.0;
    
    private static final double coincTIME = 25.; //ns. 	
    private double                   time = 0;
    private double              fgtw,dgtw = 0; //global time walk correction
    
    private TimeCorrection             tc = null; 
    
    public ECStrip(int sector, int layer, int component){
        desc.setSectorLayerComponent(sector, layer, component);       
        tc = ECCommon.useFADCTime ? new ExtendedTWCFTime() : new ExtendedTWCTime();
    }
	
    public DetectorDescriptor getDescriptor(){
    	return desc;
    }
    
    public ECStrip setADC(int adc){
        iADC = adc;
        return this;
    }
    
    public ECStrip setTDC(int tdc){
        iTDC = tdc;
        return this;
    }
    
    public ECStrip setTADC(float tdc) {
    	iTADC = tdc;
    	return this;
    }
	
    public int getADC(){
        return iADC;
    }

    public int getTDC(){
        return ECCommon.useFADCTime||iTDC==0 ? (int) (iTADC/iTimA1) : iTDC;
    }
    
    public double getRawTime(){
       	return tc.getRawTime();
    }
    
    public double getPhaseCorrectedTime() { 
        return tc.getPhaseCorrectedTime();
    }
    
    public double getRawTime(boolean phaseCorrection) {
 	    return phaseCorrection ? tc.getPhaseCorrectedTime():tc.getRawTime();
    }
    
    public double getTWCTime() {
    	return tc.getTWCTime();    	
    }
    
    public double getTime() {
    	return tc.getTime();
    }
    
    abstract class TimeCorrection {
    	public abstract double getRawTime();
    	public abstract double getPhaseCorrectedTime();
    	public abstract double getTWCTime();    	
    	public abstract double getTime();
    }
    
    public class SimpleTWCTime extends TimeCorrection {    	
        public double getRawTime(){
           	return iTDC * iTimA1;
        }
        
        public double getPhaseCorrectedTime() { 
            return iTDC * iTimA1 - triggerPhase;
        } 
         	
        public double getTWCTime() {
        	double radc = Math.sqrt(iADC);
          	return  getPhaseCorrectedTime() - iTimA2/radc;
        }  
        
    	public double getTime() {
        	double radc = Math.sqrt(iADC);
    		return getPhaseCorrectedTime() - iTimA2/radc - iTimA0;
    	}

    }
    
    public class ExtendedTWCTime extends TimeCorrection {    	
        public double getRawTime(){
           	return iTDC * iTimA1;
        }
        
        public double getPhaseCorrectedTime() { 
            return iTDC * iTimA1 - triggerPhase;
        } 
        
        public double getTWCTime() {
        	double radc = Math.sqrt(iADC);
          	return getPhaseCorrectedTime() - dgtw/radc - iTimA2/radc - iTimA3 - iTimA4/Math.sqrt(radc)          - iTim00;          	
        } 
        
    	public double getTime() {
        	double radc = Math.sqrt(iADC);
          	return getPhaseCorrectedTime() - dgtw/radc - iTimA2/radc - iTimA3 - iTimA4/Math.sqrt(radc) - iTimA0 - iTim00;          	
        }           	
    } 
    
    public class ExtendedTWCFTime extends TimeCorrection {    	
        public double getRawTime(){
           	return iTADC;
        }
        
        public double getPhaseCorrectedTime() {         	
            return iTADC;  
        }        
        
    	public double getTWCTime() {
        	double x = Math.sqrt(iADC);
        	double corr = 0;
        	if(fTimA2!=0 && fTimA4!=0) corr = fTimA2 + Math.exp(-(x-fTimA3)/fTimA4)+1-Math.exp( (x-fTimA5)/fTimA6);
          	return getRawTime() - fgtw/x - corr          - fTim00;          	
    	} 
    	
    	public double getTime() {
        	double radc = Math.sqrt(iADC);
        	double corr = 0;
        	if(fTimA2!=0 && fTimA4!=0) corr = fTimA2 + Math.exp(-(radc-fTimA3)/fTimA4)+1-Math.exp( (radc-fTimA5)/fTimA6);
          	return getRawTime() - fgtw/radc - corr - fTimA0 - fTim00;          	
    	}
    	
    }     
    
    
    public double getEnergy(){
        return iADC*iGain*iADC_to_MEV;
    }
    
    public void setDistanceEdge(double val){
        stripDistanceEdge = val;
    }
    
    public double getDistanceEdge(){
        return stripDistanceEdge;
    }
  
    public void setPeakId(int val){ 
    	peakID = val;
    }
    
    public int getPeakId(){
      	return peakID;      	
    }

    public void setID(int val){
        id = val;    
    }

    public int getID(){
        return id;
    }

    public void setStripID(int val){
        stripId = val;    
    }

    public int getStripID(){
        return stripId;
    }
    
    public void setClusterId(int val){
        clusterId = val;
    }

    public int getClusterId(){
        return clusterId;
    }
    
    public void setAttenuation(double a, double b, double c){
        iAttenLengthA = a;
        iAttenLengthB = b;
        iAttenLengthC = c;
    }
    
    public void setTriggerPhase(int val) {
        triggerPhase = val;
    } 
    
    public void setVeff(double val) {
        veff = val;
    }
    
    public double getVeff() {
 	    return veff;
    }  
    
    public void setGain(double val){
        iGain = val;
    }
    
    public void setDTiming(double a0, double a1, double a2, double a3, double a4) {
        iTimA0 = a0;
        iTimA1 = a1;
        iTimA2 = a2;
        iTimA3 = a3;
        iTimA4 = a4;
    } 
    
    public void setFTiming(double a0, double a1, double a2, double a3, double a4, double a5, double a6) {
    	fTimA0 = a0;
    	fTimA1 = a1;
    	fTimA2 = a2;
    	fTimA3 = a3;
    	fTimA4 = a4;
    	fTimA5 = a5;
    	fTimA6 = a6;
    }    
    
    public void setDtimeGlobalTimeWalk(double val) {
    	dgtw = val;
    }
    
    public void setFtimeGlobalFTimeWalk(double val) {
    	fgtw = val;
    }  
	
    public void setDtimeGlobalTimingOffset(double val) {
        iTim00 = val;
    }
    
    public void setFtimeGlobalTimingOffset(double val) {
        fTim00 = val;
    }
    
    public double[] getTiming() {
        double[] array = new double[5];
        array[0] = iTimA0;
        array[1] = iTimA1;
        array[2] = iTimA2;
        array[3] = iTimA3;
        array[4] = iTimA4;
        return array;    
    }
    
    public double[] getFTiming() {
        double[] array = new double[7];
        array[0] = fTimA0;
        array[1] = fTimA1;
        array[2] = fTimA2;
        array[3] = fTimA3;
        array[4] = fTimA4;
        array[5] = fTimA5;
        array[6] = fTimA6;
        return array;    
    }	
    
    public double getEnergy(Point3D point){
        edist = point.distance(stripLine.end());
        return iADC*iGain*iADC_to_MEV/getEcorr(-edist);
    }
    
    public double getEcorr(double dist) {
        return iAttenLengthA*Math.exp(dist/iAttenLengthB) + iAttenLengthC;    	
    }
    
    public double getTime(Point3D point) {		
        tdist = point.distance(stripLine.end());
        time =  getTime() - tdist/veff;
        return time;
    } 
	
    public double getPointTime() {
        return time;
    }
	
    public double getEdist() {return edist;}
	
    public double getTdist() {return tdist;}
    
    public Line3D getLine()  {return stripLine;}
    
    public boolean isNeighbour(ECStrip strip){
        if(strip.getDescriptor().getSector() == desc.getSector() &&
           strip.getDescriptor().getLayer()  == desc.getLayer()){
           if(Math.abs(strip.getDescriptor().getComponent()-desc.getComponent())<=ECCommon.touchID) return true;
        }
        return false;
    }
    
    public boolean isInTime(ECStrip strip) {
        if (Math.abs(getTime() - strip.getTime()) < ECStrip.coincTIME) return true;
        return false;
    } 
    
    public int compareTo(Object o) {
        ECStrip ob = (ECStrip) o;
        if (ECCommon.stripSortMethod==0) {
            if(ob.getDescriptor().getSector()     < desc.getSector())    return  1;
            if(ob.getDescriptor().getSector()     > desc.getSector())    return -1;
            if(ob.getDescriptor().getLayer()      < desc.getLayer())     return  1;
            if(ob.getDescriptor().getLayer()      > desc.getLayer())     return -1;
            if(ob.getDescriptor().getComponent() <  desc.getComponent()) return  1;
            if(ob.getDescriptor().getComponent() == desc.getComponent()) return  0;
        } else {
        	if(ob.getEnergy()                       > getEnergy())            return  1;
        	if(ob.getEnergy()                       < getEnergy())            return -1;
        }
        return -1;
    }
    
    @Override
    public String toString(){
        StringBuilder str = new StringBuilder();
        str.append(String.format("----> strip (%3d %3d %3d) ADC/TDC  %5d %5d  ENERGY=%6.4f TIME=%6.2f DIST=%6.2f PEAKID=%2d", 
                desc.getSector(),desc.getLayer(),desc.getComponent(),
                iADC,iTDC,getEnergy(),getTime(),getTdist(),getStripID()));
        str.append(String.format("  GAIN (%5.3f) ATT (%5.1f)", 
                iGain,iAttenLengthB));
        return str.toString();
    }
}
