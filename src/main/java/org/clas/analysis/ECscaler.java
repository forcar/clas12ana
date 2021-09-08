package org.clas.analysis;

import java.util.LinkedList;
import java.util.Set;

import org.clas.viewer.DetectorMonitor;
import org.jlab.detector.base.DetectorCollection;
import org.jlab.detector.base.DetectorOccupancy;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;


public class ECscaler extends DetectorMonitor {

    public static DetectorOccupancy  occupancyECAL = new DetectorOccupancy();
    public static int              occCounts = 0;
    public static int                 occMax = 10000;
    public static int                 occLo  =  2;
    public static int                 occHi  = 20;
    public static int                 occHL  = occHi-occLo+1;
    public static int                    nev = 1;
    
    public DetectorCollection<LinkedList<Integer>> fifoa = new DetectorCollection<LinkedList<Integer>>();
    public DetectorCollection<LinkedList<Integer>> fifot = new DetectorCollection<LinkedList<Integer>>();
    public DetectorCollection<Integer> anorm = new DetectorCollection<Integer>();
    public DetectorCollection<Integer> tnorm = new DetectorCollection<Integer>();
    
    HipoDataSync  writer = null;		
    
    int[]           npmt = {68,62,62,36,36,36,36,36,36};    
    String[]         det = new String[]{"pcal","ecin","ecou"};
    String[]           v = new String[]{"u","v","w"};
    
    public ECscaler(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("ECAL");

        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);

        this.init();
        this.initEPICS();
        this.localinit();
        this.localclear();
        this.initFIFO("ECAL",1,7);
    }
    
    public void localinit() {
        System.out.println("ECscaler.localinit()");
        writer = new HipoDataSync();
        writer.open("/Users/colesmith/CLAS12ANA/ECscaler/ECscaler.hipo");
    }
    
    public void localclear() {
    	System.out.println("ECscaler.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	slider.setValue(0);
    }
    
    public void openOutput() {
        writer = new HipoDataSync();
        writer.open("/Users/colesmith/CLAS12ANA/ECscaler/ECscaler.hipo");    	
    }
    
    public void initFIFO(String detName, int is1, int is2) {
        System.out.println("ECscaler:initFIFO():");
        for (int is=is1; is<is2 ; is++) {
            for (int il=1; il<layMap.get(detName).length+1 ; il++) {
                for (int ic=1; ic<nlayMap.get(detName)[il-1]+1; ic++) {
                	fifoa.add(is, il, ic,new LinkedList<Integer>());
                    fifot.add(is, il, ic,new LinkedList<Integer>());
                    anorm.add(is, il, ic,0);
                    tnorm.add(is, il, ic,0);
                }
            }
        }    	
    }
        
    @Override
    public void createHistos(int run) {
    	System.out.println("ECscaler:createHistos("+run+")");
    	if(dumpFiles) openOutput();
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;  
    	createECAL(2000);    	
    }
    
    public void createECAL(int nx) {    	
    	for(int is=1; is<7; is++) {
    		for(int im=1; im<4; im++) {
    			for(int iv=1; iv<4; iv++) {
    				int il = iv+3*(im-1);
    	    		dgm.add("ECAL",1,2, is, il, getRunNumber());
    	    		int ny=npmt[il-1]; String tity=det[im-1]+" "+v[iv-1]+" pmt"; String s = "SECTOR "+is; 
    	    		dgm.makeH2("ADC"+is+il, nx,1,nx+1,ny,1,ny+1, -1, s+" ADC", "sample", tity);
    	    		if(isNorm) dgm.cc("ADC"+is+il, false, false, 0, 0, -3, 3);
    	    		dgm.makeH2("TDC"+is+il, nx,1,nx+1,ny,1,ny+1, -1, s+" TDC", "sample", tity);
    	    		if(isNorm) dgm.cc("TDC"+is+il, false, false, 0, 0, -3, 3);
    			}
    		}
    	}
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plot("ECAL");
    }
    
	public void plot(String tabname) { 
		int is = getActiveSector(); int il = getActiveView()+3*getActiveLayer()+1;
    	dgm.drawGroup(tabname,is, il, getRunNumber());
    }
	
    @Override
    public void plotEvent(DataEvent de) {
    	System.out.println(getDetectorName()+".plotEvent");  
        analyze();         
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
    	fillHistos("ECAL",1,7);
    	if(dumpFiles) writer.close();
    	isAnalyzeDone = true;
    }	
    
    @Override
    public void processEvent(DataEvent de) {	    	
        if(de.hasBank("ECAL::adc")) occupancyECAL.addADCBank(de.getBank("ECAL::adc"));
        if(de.hasBank("ECAL::tdc")) occupancyECAL.addTDCBank(de.getBank("ECAL::tdc"));
        if(occCounts>=occMax) writeScalerBank(de); 
        fillFifo(de);       
        occCounts++;
    }
    
    public boolean inNormWindow() {
    	return occCounts>=occLo && occCounts<=occHi;
    }
    
    public void fillFifo(DataEvent de) {
        if(de.hasBank("ECAL::scaler")) {
        	DataBank  bank = de.getBank("ECAL::scaler");
        	for(int i=0; i<bank.rows(); i++) {
        		int is = bank.getByte("sector", i);
        		int il = bank.getByte("layer", i);
        		int ic = bank.getShort("component", i);
        		int ia = bank.getInt("acount", i);
        		int it = bank.getInt("tcount", i);
        		if(il>0) {
        			fifoa.get(is, il, ic).add(ia); fifot.get(is, il, ic).add(it);
        		    if(inNormWindow()) {
        		    	anorm.add(is,il,ic,anorm.get(is,il,ic)+ia);tnorm.add(is,il,ic,tnorm.get(is,il,ic)+it);
        		    }
        		}
        	}
       }    	
    }
    
    public void fillHist(DataEvent de) {
        if(de.hasBank("ECAL::scaler")) {
        	DataBank  bank = de.getBank("ECAL::scaler");
        	for(int i=0; i<bank.rows(); i++) {
        		int is = bank.getByte("sector", i);
        		int il = bank.getByte("layer", i);
        		int ic = bank.getShort("component", i);
        		int ia = bank.getInt("acount", i);
        		int it = bank.getInt("tcount", i);
        		if(il>0) {dgm.fill("ADC"+is+il,nev,ic,ia); dgm.fill("TDC"+is+il,nev,ic,it);}
        	}
    		nev++;
       }    	
    } 
    
    public void fillHistos(String detName, int is1, int is2) {   	
    	for (int is=is1; is<is2 ; is++) {
    		for (int il=1; il<layMap.get(detName).length+1 ; il++) {
    			for (int ic=1; ic<nlayMap.get(detName)[il-1]+1; ic++) {    	
            		Integer fa[] = new Integer[fifoa.get(is, il, ic).size()];
     				Integer ft[] = new Integer[fifot.get(is, il, ic).size()];
    				fifoa.get(is, il, ic).toArray(fa);
    				fifot.get(is, il, ic).toArray(ft);    				
    				for (int it=0; it<fa.length; it++) {
    					float y = (float)((float)(fa[it]-anorm.get(is, il, ic)/occHL)/Math.sqrt(fa[it]));
    					dgm.fill("ADC"+is+il,it+1,ic,isNorm ? y:fa[it]);
    				}
    				for (int it=0; it<ft.length; it++) {
    					float y = (float)((float)(ft[it]-tnorm.get(is, il, ic)/occHL)/Math.sqrt(ft[it]));
    					dgm.fill("TDC"+is+il,it+1,ic,isNorm ? y:ft[it]);
    				}
    			}
    		}
    	}
    }
    
    public int getLayer(int il) {
    	int[] lay = {0,0,0,1,1,1,2,2,2};
    	return lay[il-1];
    }
    
    public int getIndex(int is, int il, int ip) {  
        int off[] = {0,68,130,0,36,72,0,36,72};
        int sca[] = {192,192,192,216,216,216,216,216,216};
        return is>0 ? (is-1)*sca[il-1]+off[il-1]+ip:0;
    }  
	
    private void writeScalerBank(DataEvent de) {
    	
    	DetectorCollection dc = occupancyECAL.getCollection();
    	
    	Set<Integer>  sectors = dc.getSectors(); 
    	DataBank bank = de.createBank("ECAL::scaler", 2448); 

    	int n=0;
    	for(Integer sector : sectors){
    		Set<Integer> layers = dc.getLayers(sector);
    		for(Integer layer : layers){   		
    			Set<Integer> components = dc.getComponents(sector, layer);
    			for(Integer component : components){  
    				bank.setByte("sector",     n, (byte) (int) sector);
    				bank.setByte("layer",      n, (byte) (int) layer);
    				bank.setShort("component", n, (short)(int) component);
    				bank.setInt("event",       n, getEventNumber());
    				bank.setInt("acount",      n, occupancyECAL.getADC(sector, layer, component));
    				bank.setInt("tcount",      n, occupancyECAL.getTDC(sector, layer, component));
    				n++;
    			}
    		}
    	}
    	
    	occupancyECAL.reset();
    	occCounts = 0;
    	de.appendBanks(bank);
        if(dumpFiles) writer.writeEvent(de);
    }	
	
}