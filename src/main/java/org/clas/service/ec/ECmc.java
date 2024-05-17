package org.clas.service.ec;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.utils.groups.IndexedTable;

public class ECmc {
	
	int sec,lay,pmt;
	Random r = new Random();
		
	double ecc_ADC_GeV_to_evio     = 1e-4; // MIP calibration: 10(15) ch/MeV for ECAL(PCAL) 50 for ECAL Sector 5
	double ecc_pmtQE               = 0.27    ;
	double ecc_pmtDynodeGain       = 4.0     ;
	double ecc_pmtFactor           = Math.sqrt(1 + 1/(ecc_pmtDynodeGain-1));
    double     pmtPEYld            = 0;
    double ADC                     = 0;
    double A,B,C,D,E,G,fthr;
    double da0,da1,da2,da3,da4,da5,da6,da7,da8,gtw;
    double fa0,fa1,fa2,fa3,fa4,fa5,fa6,tgo,FTOFFSET,tmf,fo;
    double ftres0,ftres1,ftres2,ftres3;
    double dtres0,dtres1,dtres2,dtres3;
    double def0,def1,def2;
    int fadc=0; float DTIME_raw=0, FTIME_raw=0;
    int currentrun = 0; boolean skipInit = false;
    
    ConstantsManager manager;
    Map<String,IndexedTable> ccdb = new HashMap<>();

	public ECmc(){
	}

    public void initCCDB(int run, ConstantsManager manager) {
    	if(run==currentrun) return;
    	System.out.println("ECmc.initCCDB("+run+") + Current run: "+currentrun);
    	ccdb.clear();
    	this.manager = manager;        
        ccdb.put("gain",  manager.getConstants(run, "/calibration/ec/gain"));
        ccdb.put("ftime", manager.getConstants(run, "/calibration/ec/ftime"));     
        ccdb.put("dtime", manager.getConstants(run, "/calibration/ec/dtime"));     
        ccdb.put("gtw",   manager.getConstants(run, "/calibration/ec/global_time_walk"));
        ccdb.put("tgo",   manager.getConstants(run, "/calibration/ec/tdc_global_offset"));		
        ccdb.put("fo",    manager.getConstants(run, "/calibration/ec/fadc_offset"));         
        ccdb.put("tmf",   manager.getConstants(run, "/calibration/ec/tmf_offset"));         
        ccdb.put("fgo",   manager.getConstants(run, "/calibration/ec/fadc_global_offset"));  
        ccdb.put("fthr",  manager.getConstants(run, "/calibration/ec/fthr"));
        ccdb.put("deff",  manager.getConstants(run, "/calibration/ec/deff"));
        ccdb.put("ftres", manager.getConstants(run, "/calibration/ec/ftres"));
        ccdb.put("dtres", manager.getConstants(run, "/calibration/ec/dtres"));
        initCCDB();
        currentrun=run;
    }
	
	public void initCCDB() {		
		   G = ccdb.get("gain").getDoubleValue("gain",sec,lay,pmt);		
	     fa0 = ccdb.get("ftime").getDoubleValue("a0", sec,lay,pmt);		
	     fa1 = ccdb.get("ftime").getDoubleValue("a1", sec,lay,pmt);		
	     fa2 = ccdb.get("ftime").getDoubleValue("a2", sec,lay,pmt);		
	     fa3 = ccdb.get("ftime").getDoubleValue("a3", sec,lay,pmt);		
	     fa4 = ccdb.get("ftime").getDoubleValue("a4", sec,lay,pmt);		
	     fa5 = ccdb.get("ftime").getDoubleValue("a5", sec,lay,pmt);		
	     fa6 = ccdb.get("ftime").getDoubleValue("a6", sec,lay,pmt);
	     da0 = ccdb.get("dtime").getDoubleValue("a0", sec,lay,pmt);		
	     da1 = ccdb.get("dtime").getDoubleValue("a1", sec,lay,pmt);		
	     da2 = ccdb.get("dtime").getDoubleValue("a2", sec,lay,pmt);		
	     da3 = ccdb.get("dtime").getDoubleValue("a3", sec,lay,pmt);		
	     da4 = ccdb.get("dtime").getDoubleValue("a4", sec,lay,pmt);		
	     da5 = ccdb.get("dtime").getDoubleValue("a5", sec,lay,pmt);		
	     da6 = ccdb.get("dtime").getDoubleValue("a6", sec,lay,pmt);		
	     da7 = ccdb.get("dtime").getDoubleValue("a7", sec,lay,pmt);		
	     da8 = ccdb.get("dtime").getDoubleValue("a8", sec,lay,pmt);	

	     tgo = ccdb.get("tgo").getDoubleValue("offset", 0,0,0);	 
	     tmf = ccdb.get("tmf").getDoubleValue("offset", sec,lay,pmt);		
	      fo = ccdb.get("fo").getDoubleValue("offset",  sec,lay,0);	
	FTOFFSET = ccdb.get("fgo").getDoubleValue("global_offset", 0,0,0); 		      
         gtw = ccdb.get("gtw").getDoubleValue("time_walk",sec,lay,0);
	      
	    def0 = ccdb.get("deff").getDoubleValue("a0",  sec,lay,pmt);		
		def1 = ccdb.get("deff").getDoubleValue("a1",  sec,lay,pmt);		
		def2 = ccdb.get("deff").getDoubleValue("a2",  sec,lay,pmt);		
	    fthr = ccdb.get("fthr").getDoubleValue("thr", sec,lay,pmt);		

	    ftres0 = ccdb.get("ftres").getDoubleValue("p0",  sec,lay,0);		
	    ftres1 = ccdb.get("ftres").getDoubleValue("p1",  sec,lay,0);		
	    ftres2 = ccdb.get("ftres").getDoubleValue("p2",  sec,lay,0);		
	    ftres3 = ccdb.get("ftres").getDoubleValue("p3",  sec,lay,0);
	    
	    dtres0 = ccdb.get("dtres").getDoubleValue("p0",  sec,lay,0);		
	    dtres1 = ccdb.get("dtres").getDoubleValue("p1",  sec,lay,0);		
	    dtres2 = ccdb.get("dtres").getDoubleValue("p2",  sec,lay,0);		
	    dtres3 = ccdb.get("dtres").getDoubleValue("p3",  sec,lay,0);			      
	}
		
	public void setSLP(int sector, int layer, int component) {
		sec = sector;
		lay = layer;
		pmt = component;
		initCCDB();
	}
	
	public void addADC(int val) {
		fadc = val;
	}
	
	public void addDTDC(float val) {
		DTIME_raw = val;
	}	
	
	public void addFTDC(float val) {
		FTIME_raw = val;
	}
	
	public double dgtzFADC() {
		ADC=0;
		pmtPEYld = lay < 4 ? 11.5 : 3.5; //pcal:ecal
		double Etota = fadc * 1000 * ecc_ADC_GeV_to_evio * G;
		double EC_npe = shootP(Etota*pmtPEYld); //number of photoelectrons
		if (EC_npe > 0) {
			double sigma  = Math.sqrt(EC_npe)*ecc_pmtFactor;
			ADC = shootG(EC_npe,sigma)/1000./ecc_ADC_GeV_to_evio/G/pmtPEYld;
		}
		return ADC/10<fthr ? 0:ADC;
	}
	
	public double dgtzFTDC() {
		double radc = Math.sqrt(ADC);
		double ftim =  (fa4==0||fa6==0) ? 0 : fa0 + fa2 + Math.exp(-(radc-fa3)/fa4)+1-Math.exp( (radc-fa5)/fa6);
		double ftime_in_ns = FTIME_raw + ftim + tgo - FTOFFSET - tmf - fo;
		ftime_in_ns = shootG(ftime_in_ns,getTRES(ADC,ftres0,ftres2,ftres3,ftres1));
		return ftime_in_ns;
	}
	
	public double dgtzDTDC() {
		double radc = Math.sqrt(ADC);
		double dtim = (da4==0||da6==0||da7==0) ? 0 : da0 + gtw/radc +da2+Math.exp(-(radc-da3)/da4)+1-Math.exp(-(da5-radc)/da6)-Math.exp(-(radc-da3*0.95)/da7)*Math.pow(radc,da8);
		double dtime_in_ns = DTIME_raw + dtim + tgo;
		dtime_in_ns = shootG(dtime_in_ns,getTRES(ADC,dtres0,dtres2,dtres3,dtres1));
		if (def0>0 && dtime_in_ns > 0 && Math.random()> 1/Math.pow(1+Math.exp(-def0*(ADC/10-def1)),def2)) dtime_in_ns = 0; // DSC/TDC thresholds
		return  dtime_in_ns;
	}
	
	public double shootP(double mean) {		
		PoissonDistribution p = new PoissonDistribution(mean);
		return p.sample();
	}
	
	public double shootG(double mean, double sigma) {		
		return r.nextGaussian()*sigma + mean;
	}
		
    public double getTRES(double x, double p0, double p1, double p2, double p3) {
        return (p0*Math.exp(Math.pow(x,p3)/120000) + p1/x + p2/Math.pow(x,0.5));
    }
}
