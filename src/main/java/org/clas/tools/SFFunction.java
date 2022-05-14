package org.clas.tools;

import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.math.Func1D;
import org.jlab.rec.eb.EBCCDBConstants;
import org.jlab.rec.eb.SamplingFractions;

public class SFFunction extends Func1D{
	
    EBCCDBConstants ccdb = new EBCCDBConstants();
    DetectorParticle p = new DetectorParticle(); 
    int pid;
	    
    public SFFunction(String name, int pid, int is, EBCCDBConstants ccdb, double min, double max) {
        super(name, min, max);
        this.ccdb = ccdb;
        this.pid  = pid;
        
        p.addResponse(new CalorimeterResponse(1,1,0));
        p.getDetectorResponses().get(0).getDescriptor().setSector(is);
        p.getDetectorResponses().get(0).getDescriptor().setType(DetectorType.ECAL);
    }
    @Override
    public double evaluate(double x){        	 
    	 p.getDetectorResponses().get(0).setEnergy(x);
   	     return  SamplingFractions.getMean(pid, p, ccdb);
    }
}
