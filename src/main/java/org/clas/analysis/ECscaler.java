package org.clas.analysis;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.clas.tools.DetectorOccupancy;
import org.clas.viewer.DetectorMonitor;
import org.jlab.detector.base.DetectorCollection;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;


public class ECscaler extends DetectorMonitor {

    public static DetectorOccupancy  occupancyECAL = new DetectorOccupancy();
    public static int              occCounts = 0;
    public static int                 occMax = 10002;
    public static int                 occLo  =  1;
    public static int                 occHi  =  100;
    public static int                 occHL  = occHi-occLo+1;
    public static int                    nev = 1;
    public static int               evn_last = 0;
    public static long              tim_last = 0;
    
    public DetectorCollection<LinkedList<Integer>> fifoa = new DetectorCollection<LinkedList<Integer>>();
    public DetectorCollection<LinkedList<Integer>> fifot = new DetectorCollection<LinkedList<Integer>>();
    public DetectorCollection<Integer> anorm = new DetectorCollection<Integer>();
    public DetectorCollection<Integer> tnorm = new DetectorCollection<Integer>();
    
    public  List<Integer>                    evnlist = new ArrayList<Integer>();
    public  List<Integer>                    evrlist = new ArrayList<Integer>();
    
    HipoDataSync  writer = null;		
    
    int[]           npmt = {68,62,62,36,36,36,36,36,36};    
    String[]         det = new String[]{" PCAL"," ECIN"," ECOU"};
    String[]           v = new String[]{"U","V","W"};
    
    public ECscaler(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("RAW","NORM","TimeLine");

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
        occupancyECAL.ADCThreshold = 500;
    }
    
    public void localclear() {
    	System.out.println("ECscaler.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	evnlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(1);
    }
    
    public void openOutput(String file) {
    	System.out.println("ECscaler:openOutput("+file+")");
    	writer = new HipoDataSync();
        writer.open(file);    	
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
    	if(dumpFiles) openOutput(outPath+"ECscaler/ECscaler-"+run+".hipo");
    	setRunNumber(run);
    	runlist.add(run);    	
    	histosExist = true;  
    	TLmax = getTotalEvents();
    	createRAW(TLmax);
    	createNORM(TLmax);    	
    }
    
    public void createRAW(int nx) {    	
    	for(int is=1; is<7; is++) { //sectors 1-6
    		for(int im=1; im<4; im++) { //modules pcal,ecin,ecou
    			for(int iv=1; iv<4; iv++) { //views u,v,w
    				int il = iv+3*(im-1); //superlayer 1-9
    	    		dgm.add("RAW",1,2, is, il, getRunNumber());
    	    		int ny=npmt[il-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";  
    	    		dgm.makeH2("ADC"+is+il, nx,1,nx+1,ny,1,ny+1, -1," ADC COUNTS", "RUN INDEX", tity);    	    		
    	    		dgm.makeH2("TDC"+is+il, nx,1,nx+1,ny,1,ny+1, -1, "TDC COUNTS", "RUN INDEX", tity);    	    		
    			}
    		}
    	}
    }
    
    public void createNORM(int nx) {    	
    	for(int is=1; is<7; is++) {
    		for(int im=1; im<4; im++) {
    			for(int iv=1; iv<4; iv++) {
    				int il = iv+3*(im-1);
    	    		dgm.add("NORM",1,2, is, il, getRunNumber());
    	    		int ny=npmt[il-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";
    	    		dgm.makeH2("NADC"+is+il, nx,1,nx+1,ny,1,ny+1, -1,"NORM ADC COUNTS", "RUN INDEX", tity);
    	    		dgm.cc("NADC"+is+il, false, false, 0, 0, -3, 3);
    	    		dgm.makeH2("NTDC"+is+il, nx,1,nx+1,ny,1,ny+1, -1,"NORM TDC COUNTS", "RUN INDEX", tity);
    	    		dgm.cc("NTDC"+is+il, false, false, 0, 0, -3, 3);
    			}
    		}
    	}
    }    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);
    }
    
    public void plotSummary(int run) {
   	    if(dropSummary) return;
        setRunNumber(run);
    	plot("RAW");
    	plot("NORM");
    }
    
    @Override
    public void plotScalers() {
    	plotTimeLine("TimeLine");
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
    	fillHistFromFifo("ECAL",1,7);
    	if(dumpFiles) writer.close();
    	isAnalyzeDone = true;
    }	
    
    @Override
    public void processEvent(DataEvent de) {	    	
        if(de.hasBank("ECAL::scaler")) {fillFifoFromBank(de);return;}
        if(de.hasBank("ECAL::adc")) occupancyECAL.addADCBank(de.getBank("ECAL::adc"));
        if(de.hasBank("ECAL::tdc")) occupancyECAL.addTDCBank(de.getBank("ECAL::tdc"));
        if(occCounts>=occMax) {
        	fillFifoFromData(de);
        	if(dumpFiles) {fillBankFromData(de); writer.writeEvent(de);}
        	occupancyECAL.reset(); occCounts = 0;
        }      
       occCounts++;
    }
    
    public boolean inNormWindow() {
    	return occCounts>=occLo && occCounts<=occHi;
    }
    
    public void fillFifoFromBank(DataEvent de) {
    	int  run = de.getBank("RUN::config").getInt("run", 0);
    	int  evn = de.getBank("RUN::config").getInt("event", 0);
    	long tim = de.getBank("RUN::config").getLong("timestamp",0);
    	if (run==0) return;
    	int ev_rate = (int) (0.25e9*(evn-evn_last)/(tim-tim_last));
        runlist.add(occCounts,run); evnlist.add(occCounts,evn);evrlist.add(occCounts,evn_last>0?ev_rate:0); 
        occCounts++;
        evn_last = evn; tim_last = tim;
        DataBank  bank = de.getBank("ECAL::scaler");
        for(int i=0; i<bank.rows(); i++) {
        	int is = bank.getByte("sector", i);
        	int il = bank.getByte("layer", i);
        	int ic = bank.getShort("component", i);
        	int ia = bank.getInt("acount", i);
        	int it = bank.getInt("tcount", i);
        	if(il>0) {
        		fifoa.get(is, il, ic).add(ia); 
        		fifot.get(is, il, ic).add(it);
        		   if(inNormWindow()) {
        			   anorm.add(is,il,ic,anorm.get(is,il,ic)+ia);
        			   tnorm.add(is,il,ic,tnorm.get(is,il,ic)+it);
        		   }
        	}
        }
    	
    }
    
    public void fillFifoFromData(DataEvent de) {
    	
    	DetectorCollection dc = occupancyECAL.getCollection();
    	
    	Set<Integer>  sectors = dc.getSectors(); 
    	for(Integer is : sectors){
    		Set<Integer> layers = dc.getLayers(is);
    		for(Integer il : layers){   		
    			Set<Integer> components = dc.getComponents(is,il);
    			for(Integer ic : components){
    				fifoa.get(is, il, ic).add(occupancyECAL.getADC(is,il,ic)); 
    				fifot.get(is, il, ic).add(occupancyECAL.getTDC(is,il,ic));
        		    if(inNormWindow()) {
        		    	anorm.add(is,il,ic,anorm.get(is,il,ic)+occupancyECAL.getADC(is,il,ic));
        		    	tnorm.add(is,il,ic,tnorm.get(is,il,ic)+occupancyECAL.getTDC(is,il,ic));
        		    }        		    
    			}
    		}
    	}   				    	
    	
    } 
    
    public void fillHistFromBank(DataEvent de) {
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
   
    public void fillHistFromFifo(String detName, int is1, int is2) {   	
    	for (int is=is1; is<is2 ; is++) {
    		for (int il=1; il<layMap.get(detName).length+1 ; il++) {
    			for (int ic=1; ic<nlayMap.get(detName)[il-1]+1; ic++) {    	
            		Integer fa[] = new Integer[fifoa.get(is, il, ic).size()];
     				Integer ft[] = new Integer[fifot.get(is, il, ic).size()];
    				fifoa.get(is, il, ic).toArray(fa);
    				fifot.get(is, il, ic).toArray(ft);    				
    				for (int it=0; it<fa.length; it++) {
    					float y = (float)((float)(fa[it]-anorm.get(is, il, ic)/occHL)/Math.sqrt(fa[it]));
    					dgm.fill( "ADC"+is+il,it+1,ic,fa[it]);
    					dgm.fill("NADC"+is+il,it+1,ic,y);
    					    				}
    				for (int it=0; it<ft.length; it++) {
    					float y = (float)((float)(ft[it]-tnorm.get(is, il, ic)/occHL)/Math.sqrt(ft[it]));
    					dgm.fill( "TDC"+is+il,it+1,ic,ft[it]);
    					dgm.fill("NTDC"+is+il,it+1,ic,y);
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
    
    private void fillBankFromFifo(DataEvent de) {
	
    }
	
    private void fillBankFromData(DataEvent de) {
    	
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
    				bank.setInt("acount",      n, occupancyECAL.getADC(sector, layer, component));
    				bank.setInt("tcount",      n, occupancyECAL.getTDC(sector, layer, component));
    				n++;
    			}
    		}
    	}
    	if(de.hasBank("ECAL::adc")) de.removeBanks("ECAL::adc","ECAL::tdc");
    	de.appendBanks(bank);
    }
    
/*   TIMELINES */
    
    public void plotTimeLine(String tab) {
    	
        EmbeddedCanvas c = getCanvas(tab); c.clear(); c.divide(2, 2);
        
        int is = getActiveSector(), il = getActiveView()+3*getActiveLayer()+1;  
        
        H2F h2a = dgm.getH2F("ADC"+is+il); 
        H2F h2t = dgm.getH2F("TDC"+is+il); 
        
        ArrayList<H1F> aslice = h2a.getSlicesX();  
        ArrayList<H1F> tslice = h2t.getSlicesX();
        
    	DataLine line3 = new DataLine(runIndexSlider,  1,  runIndexSlider,  npmt[il-1]+1);  line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,  runIndexSlider+1,npmt[il-1]+1);  line4.setLineColor(5);  
    	
    	String tit1 = "RUN "+runlist.get(runIndexSlider-1)+"   EVENT "+evnlist.get(runIndexSlider-1);
    	String tit2 = "   EV/SEC "+evrlist.get(runIndexSlider-1);
    	H1F h1a = aslice.get(runIndexSlider-1); h1a.setTitle(tit1+(runlist.size()>0?tit2:" "));
    	H1F h1t = tslice.get(runIndexSlider-1); h1t.setTitle(" ");
    	    	
    	H1F h1ar = sumSlices(aslice,200,300); h1ar.setLineColor(3); float amax = (float) h1ar.getMax()*1.1f;
    	H1F h1tr = sumSlices(tslice,200,300); h1tr.setLineColor(3); float tmax = (float) h1tr.getMax()*1.1f;
    	
    	
    	c.cd(0); c.getPad().setTitleFontSize(24); c.draw(h2a); c.draw(line3);c.draw(line4);
    	c.cd(1); c.getPad().setTitleFontSize(24); c.getPad().getAxisY().setRange(0, amax); c.draw(h1a); c.draw(h1ar,"same");
    	c.cd(2); c.getPad().setTitleFontSize(24); c.draw(h2t); c.draw(line3);c.draw(line4);
    	c.cd(3); c.getPad().setTitleFontSize(24); c.getPad().getAxisY().setRange(0, tmax); c.draw(h1t); c.draw(h1tr,"same");
    	    	
    }
    
    public H1F sumSlices(ArrayList<H1F> list, int i1, int i2) {
    	H1F h = list.get(i1-1).histClone("dum");
    	for (int i=i1; i<=i2; i++) h.add(list.get(i));
    	h.normalize(i2-i1+1);
    	return h;
    }
	
}