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
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.utils.groups.IndexedList;

public class ECscaler extends DetectorMonitor {

    static int              occCounts = 0;
    static int                 occMax = 10002, counter=0;
    static int                 occLo  =  1;
    static int                 occHi  =  100;
    static int                 occHL  = occHi-occLo+1;
    static int                    nev = 1;
    static int               evn_last = 0;
    static long              tim_last = 0;
    static int                prevRun = 0;
    static DataEvent        prevEvent = null;
    static boolean          singleRun = false;
    static boolean        scalerEvent = false;
    
    static DetectorOccupancy               occupancyECAL = new DetectorOccupancy(); 
    
    DetectorCollection<LinkedList<Integer>> fifoac = new DetectorCollection<LinkedList<Integer>>();
    DetectorCollection<LinkedList<Integer>> fifotc = new DetectorCollection<LinkedList<Integer>>();
    DetectorCollection<LinkedList<Integer>> fifoav = new DetectorCollection<LinkedList<Integer>>();
    DetectorCollection<LinkedList<Integer>> fifotv = new DetectorCollection<LinkedList<Integer>>();
    DetectorCollection<LinkedList<Integer>> fifotr = new DetectorCollection<LinkedList<Integer>>();
    DetectorCollection<Integer>              anorm = new DetectorCollection<Integer>();
    DetectorCollection<Integer>              tnorm = new DetectorCollection<Integer>();
    DetectorCollection<Float>                 asum = new DetectorCollection<Float>();
    DetectorCollection<Float>                 tsum = new DetectorCollection<Float>(); 
    List<Integer>                          evnlist = new ArrayList<Integer>();
    List<Integer>                          evrlist = new ArrayList<Integer>();   
    IndexedList<ArrayList<H1F>>             ATData = new IndexedList<>(3);  
    IndexedList<H1F>                       NATData = new IndexedList<>(3);  
    
    HipoDataSync  writer = null;		
    
    int[]           npmt = {68,62,62,36,36,36,36,36,36};    
    String[]         det = new String[]{" PCAL"," ECIN"," ECOU"};
    String[]           v = new String[]{"U","V","W"};
    
    double aYL,aYS,aYR,tYL,tYS,tYR;
    
    public ECscaler(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("ATDATA","TIMELINE","STATUS","TRIGGER","SUMMARY");

        this.useSCALERButtons(true);
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
        occupancyECAL.ADCWindow[0]=10;  occupancyECAL.ADCWindow[1]=100;
        occupancyECAL.TDCWindow[0]=200; occupancyECAL.TDCWindow[1]=300;
    }
    
    public void localclear() {
    	System.out.println("ECscaler.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
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
                    fifoac.add(is, il, ic, new LinkedList<Integer>());
                    fifotc.add(is, il, ic, new LinkedList<Integer>());
                    fifoav.add(is, il, ic, new LinkedList<Integer>());
                    fifotv.add(is, il, ic, new LinkedList<Integer>());
                    anorm.add(is, il, ic,0);
                    tnorm.add(is, il, ic,0);
                }
            }
        }
        for (int ib=0; ib<32; ib++) fifotr.add(0, 0, ib, new LinkedList<Integer>());
    }
        
    @Override
    public void createHistos(int run) {
    	System.out.println("ECstatus:createHistos("+run+")");
    	if(dumpFiles) openOutput(outPath+"ECscaler/ECscaler-"+run+".hipo");
    	setRunNumber(run);    	   	
    	histosExist = true;  
    	TLmax = getTotalEvents()>occMax ? getTotalEvents()/occMax : getTotalEvents();    
    	if(useATDATA) createATDATA();
    	createTIMELINE(0,TLmax);
    	createTIMELINE(1,TLmax); 
    	createTIMELINE(2,TLmax);
    	createTIMELINE(3,TLmax); 
    	createTIMELINE(4,TLmax);
    	createTIMELINE(5,TLmax);
    	createTRIGGER(0);
    	createSTATUS(run);
    }
    
    public void createATDATA() {
    	for (int is=1; is<7; is++) {
    		for (int im=1; im<4; im++) {
    			int hl = 10*is+im;
    			String tity="SEC"+is+det[im-1]+" PMT";
    			dgm.add("ATDATA", 2, 1, is, im, 0);
    			dgm.makeH2("ADD"+hl, 100, 1, 600, im==1?192:108, 1, im==1?193:109, -1, "","ADC", tity);
    			dgm.makeH2("TDD"+hl, 100, 100, 400, im==1?192:108, 1, im==1?193:109, -1, "","TDC", tity);
    		}
    	}
    	
    }
    
    public void createTIMELINE(int st, int nx) {    	
    	switch (st) {
    	case 0:    	
        	for(int is=1; is<7; is++) { //sectors 1-6
        		for(int im=1; im<4; im++) { //modules pcal,ecin,ecou
        			for(int iv=1; iv<4; iv++) { //views u,v,w
        				int sl = iv+3*(im-1); //superlayer 1-9
        				int hl = 10*is+sl;    //hyperlayer 11-69
        				dgm.add("TIMELINE",1,2, hl, st, getRunNumber());
        				int ny=npmt[sl-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";  
        				dgm.makeH2("ADC"+hl, nx,0,nx,ny,1,ny+1, -1,"ADC COUNTS", "RUN INDEX", tity);    	    		
        				dgm.makeH2("TDC"+hl, nx,0,nx,ny,1,ny+1, -1,"TDC COUNTS", "RUN INDEX", tity);	
        			}
        		}
        	}
        	break;
    	case 1:
        	for(int is=1; is<7; is++) {
        		for(int im=1; im<4; im++) {
        			for(int iv=1; iv<4; iv++) {
        				int sl = iv+3*(im-1); //superlayer 1-9
        				int hl = 10*is+sl;    //hyperlayer 11-69
        	    		dgm.add("TIMELINE",1,2, hl, st, getRunNumber());
        	    		int ny=npmt[sl-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";
        	    		dgm.makeH2("NADC"+hl, nx,0,nx,ny,1,ny+1, -1,"NORM ADC COUNTS", "RUN INDEX", tity);
        	    		dgm.cc("NADC"+hl, false, false, 0, 0, -3, 3);
        	    		dgm.makeH2("NTDC"+hl, nx,0,nx,ny,1,ny+1, -1,"NORM TDC COUNTS", "RUN INDEX", tity);
        	    		dgm.cc("NTDC"+hl, false, false, 0, 0, -3, 3);
        			}
        		}
        	}
        	break;
    	case 2:    	
        	for(int is=1; is<7; is++) { //sectors 1-6
        		for(int im=1; im<4; im++) { //modules pcal,ecin,ecou
        			for(int iv=1; iv<4; iv++) { //views u,v,w
        				int sl = iv+3*(im-1); //superlayer 1-9
        				int hl = 10*is+sl;    //hyperlayer 11-69
        				dgm.add("TIMELINE",1,2, hl, st, getRunNumber());
        				int ny=npmt[sl-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";  
        				dgm.makeH2("VADC"+hl, nx,0,nx,ny,1,ny+1, -1,"ADC VALUES", "RUN INDEX", tity);    	    		
        				dgm.cc("VADC"+hl, false, false, 0, 0, 0, 100);
        				dgm.makeH2("VTDC"+hl, nx,0,nx,ny,1,ny+1, -1,"TDC VALUES", "RUN INDEX", tity);	
        				dgm.cc("VTDC"+hl, false, false, 0, 0, 200, 300);
        			}
        		}
        	}
        	break;
    	case 3:
        	for(int is=1; is<7; is++) {
        		for(int im=1; im<4; im++) {
        			for(int iv=1; iv<4; iv++) {
        				int sl = iv+3*(im-1); //superlayer 1-9
        				int hl = 10*is+sl;    //hyperlayer 11-69
        	    		dgm.add("TIMELINE",1,2, hl, st, getRunNumber());
        	    		int ny=npmt[sl-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";
        	    		dgm.makeH2("NVADC"+hl, nx,0,nx,ny,1,ny+1, -1,"NORM ADC VALUES", "RUN INDEX", tity);
        	    		dgm.cc("NVADC"+hl, false, false, 0, 0, 0.5f,1.5f);
        	    		dgm.makeH2("NVTDC"+hl, nx,0,nx,ny,1,ny+1, -1,"NORM TDC VALUES", "RUN INDEX", tity);
        	    		dgm.cc("NVTDC"+hl, false, false, 0, 0, -3, 3);
        			}
        		}
        	}
        	break;
    	case 4:
    		dgm.add("TIMELINE", 1, 2, 0, st, getRunNumber());
    		dgm.makeH2("TRIG",  nx, 0, nx, 32,-0.5,31.5, -1, "TRIGGER COUNTS", "RUN INDEX", "TRIGGER BITS");        	
    		break;
    	case 5:
    		dgm.add("TIMELINE", 1, 2, 0, st, getRunNumber());
    		dgm.makeH2("NTRIG", nx, 0, nx, 32,-0.5,31.5, -1, "NORM TRIGGER COUNTS", "RUN INDEX", "TRIGGER BITS");    	        	        	
    		dgm.cc("NTRIG", false, false, 0, 0, -3, 3);	
    	}
    	
    }
    
    public void createTRIGGER(int run) {
    	dgm.add("TRIGGER", 1, 2, 0, 0, run);
    	dgm.makeH1("trigmon",32,-0.5,31.5,-1,"TRIGGER","TRIGGER BITS","",1,1);
    }
    
    public void createSTATUS(int run) {    	
    	dgm.add("STATUS",3,6,0,0,run);
    	for(int is=1; is<7; is++) {
    		for(int im=1; im<4; im++) {
    			int nx = im==1 ? 68:36;
    			dgm.makeH2("STATUS"+is+im,nx,1,nx+1,3,1,4, -1,is==1 ? det[im-1]:" ", "SECTOR "+is, " "); 
    			dgm.cc("STATUS"+is+im, false, false, 0, 0, 0, 1);
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
        plot("TRIGGER");
        if(useATDATA) plot("ATDATA");
    }
    
    @Override
    public void plotScalers(int run) {
    	if(!isAnalyzeDone) return;
    	setRunNumber(run);
    	plotTimeLine("TIMELINE");
    	plotStatus("STATUS");
    	plotTLSummary("SUMMARY");
    }
    
	public void plot(String tabname) { 		
//		if(tabname=="STATUS")  {dgm.drawGroup(tabname,0, 0, getRunNumber()); return;}
		if(tabname=="TRIGGER") dgm.drawGroup(tabname, 0,0,0);
		int is = getActiveSector(); int im = getActiveLayer()+1;
    	dgm.drawGroup(tabname, is, im, 0);
    }
	
    @Override
    public void plotEvent(DataEvent de) {
    	System.out.println(getDetectorName()+".plotEvent");  
        analyze();         
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() ");  
    	fillHistFromFifo("ECAL",1,7);
    	analyzeSTATUS("ECAL",1,7);
    	analyzeNORM("ECAL",1,7);
    	if(dumpFiles) writer.close();
    	isAnalyzeDone = true;
    }	
    
    @Override
    public void processEvent(DataEvent de) {	 
    	
        if (!testTriggerMask()) {System.out.println("TriggerMask Fail Run "+getRunNumber()); return;}
        
        int run = getNewRunNumber(de);
    	
    	if(run<1) return;
    	
    	setRunNumber(run);
    	singleRun = prevRun==run;   
    	
        if (de.hasBank("ECAL::scaler")) {doScalerEvent(de); return;}
        
        if(occCounts>=occMax || (prevRun>0 && run!=prevRun)) {       	
        	fillFifoFromData();
        	doWriteEvent();
        	occupancyECAL.reset(); dgm.getH1F("trigmon").reset(); occCounts = 0;
        } 
        
        if(de.hasBank("ECAL::adc")) occupancyECAL.addADCBank(de.getBank("ECAL::adc"));
        if(de.hasBank("ECAL::tdc")) occupancyECAL.addTDCBank(de.getBank("ECAL::tdc"));
        
        fillTRIGGER(); 
                
        occCounts++;
        prevEvent = de;
        
        if(useATDATA) {fillATDATA(de); occupancyECAL.resetValue();}
    }
    
    public int getNewRunNumber(DataEvent de) {
    	prevRun = getRunNumber(); 
    	return de.hasBank("RUN::config") ? de.getBank("RUN::config").getInt("run", 0):-1;  
    }
    
    public void doScalerEvent(DataEvent de) {
    	processRUNCONFIG(de); fillFifoFromBank(de);     	
    }
    
    public void doWriteEvent() {
    	if(dumpFiles) {fillBankFromData(prevEvent); writer.writeEvent(prevEvent);}
    }
    
    @Override
    public void doSTOPEvent(DataEvent de) {
    	doWriteEvent();
    	if (de.hasBank("ECAL::scaler")) {
    		setRunNumber(getNewRunNumber(de)); 
    		doScalerEvent(de);
    	}
    }
    
    public boolean inNormWindow() {
    	return occCounts>=occLo && occCounts<=occHi;
    }
    
    public void processRUNCONFIG(DataEvent de) {
    	int  evn = de.getBank("RUN::config").getInt("event", 0);
    	long tim = de.getBank("RUN::config").getLong("timestamp",0);
    	int ev_rate = (int) (0.25e9*(evn-evn_last)/(tim-tim_last));
        runlist.add(occCounts,getRunNumber()); 
        evnlist.add(occCounts,evn);evrlist.add(occCounts,evn_last>0?ev_rate:0); 
        occCounts++;
        evn_last = evn; tim_last = tim;   	
    }
        
	public void fillTRIGGER() {
		for (int j=0; j<32; j++) if(isTrigBitSet(j)) dgm.fill("trigmon",j); 
	}
    
    public void fillATDATA(DataEvent de) {
    	
    	DetectorCollection dc = occupancyECAL.getValueCollection();
    	
    	Set<Integer>  sectors = dc.getSectors(); 
    	for(Integer is : sectors){
    		Set<Integer> layers = dc.getLayers(is);
    		for(Integer il : layers){   		
    			Set<Integer> components = dc.getComponents(is,il);
    			for(Integer ic : components){  
    				int hl = 10*is+getDET(il);
    				dgm.fill("ADD"+hl, occupancyECAL.getADCV(is, il, ic),getPMT(il,ic));
    				dgm.fill("TDD"+hl, occupancyECAL.getTDCV(is, il, ic),getPMT(il,ic));
    			}
    		}
    	}
    }
    
    public void fillFifoFromBank(DataEvent de) {

    	DataBank bank = null;
    	
        bank = de.getBank("ECAL::scaler");
        for (int i=0; i<bank.rows(); i++) {
        	int  is = bank.getByte("sector", i);
        	int  il = bank.getByte("layer", i);
        	int  ic = bank.getShort("component", i);
        	int  ia = bank.getInt("acount", i);
        	int  it = bank.getInt("tcount", i);
        	int iav = bank.getInt("avalue", i);
        	int itv = bank.getInt("tvalue", i);
        	if(il>0) {
        		fifoac.get(is, il, ic).add(ia); 
        		fifotc.get(is, il, ic).add(it);
        		fifoav.get(is, il, ic).add(ia>0?iav/ia:0); 
        		fifotv.get(is, il, ic).add(it>0?itv/it:0);
        		   if(inNormWindow()) {
        			   anorm.add(is,il,ic,anorm.get(is,il,ic)+ia);
        			   tnorm.add(is,il,ic,tnorm.get(is,il,ic)+it);
        		   }
        	}
        }
        
        bank = de.getBank("ECAL::trigger");
        for (int i=0; i<bank.rows(); i++) {
        	int ib = bank.getShort("bit", i);
        	int ic = bank.getInt("counts", i);
        	fifotr.get(0,0,ib).add(ic);
        }

    }
    
    public int getPMT(int il, int ic) {
    	int off[] = {0,68,130,0,36,72,0,36,72};
    	return ic+off[il-1];
    }
    
    public int getDET(int il) {
    	int det[] = {1,1,1,2,2,2,3,3,3};
    	return det[il-1];
    }

    public void fillFifoFromData() {
    	
    	DetectorCollection dc = occupancyECAL.getCollection();
    	
    	Set<Integer>  sectors = dc.getSectors(); 
    	for(Integer is : sectors){
    		Set<Integer> layers = dc.getLayers(is);
    		for(Integer il : layers){   		
    			Set<Integer> components = dc.getComponents(is,il);
    			for(Integer ic : components){
    				fifoac.get(is,il,ic).add(occupancyECAL.getADC(is,il,ic)); 
    				fifotc.get(is,il,ic).add(occupancyECAL.getTDC(is,il,ic));
    				fifoav.get(is,il,ic).add(occupancyECAL.getADCVC(is,il,ic)); 
    				fifotv.get(is,il,ic).add(occupancyECAL.getTDCVC(is,il,ic));
        		    if(inNormWindow()) {
        		    	anorm.add(is,il,ic,anorm.get(is,il,ic)+occupancyECAL.getADC(is,il,ic));
        		    	tnorm.add(is,il,ic,tnorm.get(is,il,ic)+occupancyECAL.getTDC(is,il,ic));
        		    }        		    
    			}
    		}
    	}   				    	
    	for (int ib=0; ib<32; ib++) fifotr.get(0,0,ib).add((int)dgm.getH1F("trigmon").getBinContent(ib));
    } 
    
    public void fillHistFromBank(DataEvent de) { //this bypasses fifo creation and precludes using NORM feature
        if(de.hasBank("ECAL::scaler")) {
        	DataBank  bank = de.getBank("ECAL::scaler");
        	for(int i=0; i<bank.rows(); i++) {
        		int is = bank.getByte("sector", i);
        		int il = bank.getByte("layer", i);
        		int ic = bank.getShort("component", i);
        		int ia = bank.getInt("acount", i);
        		int it = bank.getInt("tcount", i);
        		int iav = bank.getInt("avalue", i);
        		int itv = bank.getInt("tvalue", i);
        		int hl = 10*is+il;
        		if(il>0) {dgm.fill( "ADC"+hl,nev,ic,ia);            dgm.fill( "TDC"+hl,nev,ic,it);
        		          dgm.fill("VADC"+hl,nev,ic,ia>0?iav/ia:0); dgm.fill("VTDC"+hl,nev,ic,it>0?itv/it:0);}
        	}
    		nev++;
       }    	
    } 
   
    public void fillHistFromFifo(String detName, int is1, int is2) {   	
        System.out.println("ECscaler:fillHistFromFifo("+detName+","+is1+","+is2+")");
    	for (int is=is1; is<is2 ; is++) {
    		for (int il=1; il<layMap.get(detName).length+1 ; il++) {
    			int hl = 10*is+il;
				dgm.getH2F( "ADC"+hl).reset();  dgm.getH2F( "TDC"+hl).reset(); 
				dgm.getH2F("NADC"+hl).reset();  dgm.getH2F("NTDC"+hl).reset(); 
				dgm.getH2F("VADC"+hl).reset();  dgm.getH2F("VTDC"+hl).reset(); 
				dgm.getH2F("NVADC"+hl).reset(); dgm.getH2F("NVTDC"+hl).reset(); 
    			for (int ic=1; ic<nlayMap.get(detName)[il-1]+1; ic++) {    	
            		Integer fa[]  = new Integer[fifoac.get(is, il, ic).size()];
     				Integer ft[]  = new Integer[fifotc.get(is, il, ic).size()];
            		Integer fav[] = new Integer[fifoav.get(is, il, ic).size()];
     				Integer ftv[] = new Integer[fifotv.get(is, il, ic).size()];
    				fifoac.get(is, il, ic).toArray(fa);
    				fifotc.get(is, il, ic).toArray(ft);
    				fifoav.get(is, il, ic).toArray(fav);
    				fifotv.get(is, il, ic).toArray(ftv);
    				for (int it=0; it<fa.length; it++) {
    					float y = (float)((float)(fa[it]-getNorm(0,is,il,ic))/Math.sqrt(fa[it]));
    					dgm.fill( "ADC"+hl,it,ic,fa[it]);     					
    					dgm.fill("NADC"+hl,it,ic,y);
    				}   				
    				for (int it=0; it<ft.length; it++) {
    					float y = (float)((float)(ft[it]-getNorm(1,is,il,ic))/Math.sqrt(ft[it]));
    					dgm.fill( "TDC"+hl,it,ic,ft[it]); 
    					dgm.fill("NTDC"+hl,it,ic,y);
    				}
    				for (int it=0; it<fav.length; it++) {
    					float y = (float)((float)(fav[it]/getNorm(20,is,il,ic)));    					
    					dgm.fill( "VADC"+hl,it,ic,fav[it]);
    					dgm.fill("NVADC"+hl,it,ic,y);
    				}
    				for (int it=0; it<ftv.length; it++) {
    					float y = (float)((float)(ftv[it]-getNorm(21,is,il,ic)));    					
    					dgm.fill( "VTDC"+hl,it,ic,ftv[it]);
    					dgm.fill("NVTDC"+hl,it,ic,y);
    				} 	
    			}
    		}
    	}
    	
    	dgm.getH2F("TRIG").reset(); dgm.getH2F("NTRIG").reset();
    	for (int ib=0; ib<32; ib++) {
    		Integer ftr[] = new Integer[fifotr.get(0, 0, ib).size()]; 
    		fifotr.get(0, 0, ib).toArray(ftr);
    		for (int it=0; it<ftr.length; it++) {
    			float y = (float)((float)(ftr[it]-getNorm(40,0,0,ib+1))/Math.sqrt(ftr[it]));
    			dgm.fill("TRIG", it, ib, ftr[it]);
    			dgm.fill("NTRIG",it, ib, y);
    		}
    	}
    }
    
    private void fillBankFromFifo(DataEvent de) {
	
    }
	
    private void fillBankFromData(DataEvent de) {
    	
    	DataBank bank = null;
    	DetectorCollection dc = occupancyECAL.getCollection();
    	
    	Set<Integer>  sectors = dc.getSectors(); 
    	bank = de.createBank("ECAL::scaler", 2448); 

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
    				bank.setInt("avalue",      n, occupancyECAL.getADCVC(sector, layer, component));
    				bank.setInt("tvalue",      n, occupancyECAL.getTDCVC(sector, layer, component));
    				n++;
    			}
    		}
    	}
    	
    	if(de.hasBank("ECAL::adc")) de.removeBanks("ECAL::adc","ECAL::tdc");
    	
    	de.appendBanks(bank);
    	
    	bank = de.createBank("ECAL::trigger", 32); 

    	for (int i=0; i<32; i++) {
    		bank.setShort("bit", 	i, (short) i);    		
    		bank.setInt("counts",	i, (int) dgm.getH1F("trigmon").getBinContent(i));    		
    	}
    	
    	de.appendBanks(bank);
    	
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

    
    public void plotStatus(String tab) {           
    	EmbeddedCanvas c = getCanvas(tab); c.clear(); c.divide(3, 6);
        int n=0;
        for (int is=1; is<7; is++) {
            for (int im=1; im<4; im++) {      
                c.cd(n); c.getPad(n).getAxisZ().setRange(0.0, 1.0); c.draw(dgm.getH2F("STATUS"+is+im)); n++;                   
            }
        }            
    }
    
    public void analyzeSTATUS(String detName, int is1, int is2) {    
        System.out.println("ECscaler:AnalyzeSTATUS("+detName+","+is1+","+is2+")");
        asum.clear(); tsum.clear();
    	for (int sl=1; sl<layMap.get(detName).length+1 ; sl++) {
    		for (int ip=1; ip<nlayMap.get(detName)[sl-1]+1; ip++) {
    			float aint = 0, tint = 0; int acnt=0, tcnt=0;
                for (int is=is1; is<is2; is++) {
                	int hl = 10*is+sl;
                    asum.add(is,sl,ip,(float)dgm.getH2F("ADC"+hl).sliceY(ip-1).integral());
                    tsum.add(is,sl,ip,(float)dgm.getH2F("TDC"+hl).sliceY(ip-1).integral());
                    acnt+=(asum.get(is, sl, ip)>0?1:0);
                    tcnt+=(tsum.get(is, sl, ip)>0?1:0);
            		aint+=(asum.get(is, sl, ip)>0?asum.get(is, sl, ip):0);
            		tint+=(tsum.get(is, sl, ip)>0?tsum.get(is, sl, ip):0); 
                }
                asum.add(7,sl,ip,aint/acnt);
                tsum.add(7,sl,ip,tint/tcnt);
            }
        }
  	
        for(int is=is1; is<is2; is++) {
            for (int im=1; im<4 ; im++) {            	
            	for(int il=1; il<4; il++) {
                	int sl = il+(im-1)*3;  
//                	dgm.getH2F("STATUS"+is+im).reset();
                	for (int ip=1; ip<nlayMap.get(detName)[sl-1]+1; ip++) {                 		
                		Integer status = getStatus(is,sl,ip);
//                		calib.setIntValue(status,"status", is, sl, ip);  
                		dgm.fill("STATUS"+is+im,(float)ip, (float)il, getPlotStatus(status));
                	}
                }   
            }
        }      
    }
    
    public Integer getStatus(int is, int sl, int ip) {
    	float Asum = asum.get(7, sl, ip), A = asum.get(is, sl, ip);
    	float Tsum = tsum.get(7, sl, ip), T = tsum.get(is, sl, ip);
    	float aAsym = (A-Asum)/(A+Asum), tAsym = (T-Tsum)/(T+Tsum);
    	Boolean   badA = A==0 && T>0;   //dead ADC good TDC
        Boolean   badT = T==0 && A>0;   //dead TDC good ADC
        Boolean  badAT = A==0 && T==0;  //dead PMT or HV
    	Boolean nnbadA = aAsym < -0.85; //dead but noisy
    	Boolean nnbadT = tAsym < -0.85; //dead but noisy
    	Boolean  nbadA = aAsym < -0.30; //low gain, bad cable, high threshold
    	Boolean  nbadT = tAsym < -0.30; //low gain, bad cable, high threshold
    	Boolean  pbadA = aAsym >  0.30; //noisy, light leak
    	Boolean  pbadT = tAsym >  0.30; //noisy, light leak
    	if (badA && !badT) return 1;
    	if (badT && !badA) return 2;
    	if (badAT)         return 3;
    	if (nnbadA)        return 1;
    	if (nnbadT)        return 2;
    	if (nbadA)         return 4;
    	if (nbadT)         return 5;
    	if (pbadA)         return 6;
    	if (pbadT)         return 7;
        return 0;
    }
      
    public double getPlotStatus(int stat) {            
        switch (stat) 
        {
        case 0: return 0.0;  
        case 1: return 0.60; 
        case 2: return 0.48;  
        case 3: return 0.02;  
        case 4: return 0.75;
        case 5: return 0.85; 
        case 6: return 0.99;
        case 7: return 1.10; 
        }
    return 0.48;
        
    }
    
    public Boolean goodA() {return aYS>10000;}
    public Boolean goodT() {return tYS>200;}
    public Boolean badA()  {return (aYS<10000)&&(tYS>200);}
    public Boolean badT()  {return (tYS<200)&&(aYS>10000);}
    public Boolean badAT() {return (tYS<200)&&(aYS<10000);}
    public double  getLL() {return (goodT()) ? tYL/tYS:0.0;}  
    
/*   TIMELINES */
    
    public void plotTimeLine(String tab) {
    	
    	H1F h1a=new H1F(),h1ar=new H1F(),h1t=new H1F(),h1tr=new H1F();
    	
        EmbeddedCanvas c = getCanvas(tab); c.clear(); c.divide(2, 2);
        
        int as = getActiveSCAL();
        int it = as*10;
        int is = as==4 ? 0:getActiveSector();
        int sl = as==4 ? 0:getActiveView()+3*getActiveLayer()+1; 
    	int hl = 10*is + sl;
    	int st = as + (dNorm ? 1:0);
        String opstat = as==4 ? "":"1000000";
        
    	DataLine line3 = new DataLine(runIndexSlider,  as==4?-0.5:1,  runIndexSlider,  (as==4?31.5:npmt[sl-1])+1);  line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,as==4?-0.5:1,  runIndexSlider+1,(as==4?31.5:npmt[sl-1])+1);  line4.setLineColor(5);  
    	
    	String tit1 = "RUN "+runlist.get(runIndexSlider)+"   EVENT "+evnlist.get(runIndexSlider);
    	String tit2 = "   EV/SEC "+evrlist.get(runIndexSlider);
    	
    	h1a = ATData.getItem(is,sl,it+0).get(runIndexSlider); h1a.setTitle(tit1+(singleRun?tit2:" ")); h1a.setFillColor(21); h1a.setOptStat(opstat);
    	if(as<4) h1t = ATData.getItem(is,sl,it+1).get(runIndexSlider); h1t.setTitle(" ");              h1t.setFillColor(21); h1t.setOptStat(opstat); 
    		
    	float amax = (float) h1a.getMax()*1.3f, tmax = (float) h1t.getMax()*1.3f, amin=0, tmin=0;
    	
    	if (normrun>0) {    		
    		         h1ar = NATData.getItem(is,sl,it+0); h1ar.setLineColor(3); amax = (float) h1ar.getMax()*1.3f; h1ar.setLineWidth(5); h1ar.setOptStat(opstat);
    		if(as<4) h1tr = NATData.getItem(is,sl,it+1); h1tr.setLineColor(3); tmax = (float) h1tr.getMax()*1.3f; h1tr.setLineWidth(5); h1tr.setOptStat(opstat);
    	}
    	
    	if(as==2) {amin=30; amax=80; tmin=220; tmax=280;}
    	
    	c.cd(0); c.getPad().setTitleFontSize(24); dgm.draw("TIMELINE", c, hl, st, 0); c.draw(line3); c.draw(line4);
    	c.cd(1); c.getPad().setTitleFontSize(24); c.getPad().getAxisY().setRange(amin, amax); 
    	if(h1a.getIntegral()>0) c.draw(h1a); if (normrun>0) c.draw(h1ar,"same");
    	if(as==4) return;
    	
    	c.cd(2); c.getPad().setTitleFontSize(24); dgm.draw("TIMELINE", c, hl, st, 1); c.draw(line3); c.draw(line4);
    	c.cd(3); c.getPad().setTitleFontSize(24); c.getPad().getAxisY().setRange(tmin, tmax); 
    	if(h1t.getIntegral()>0) c.draw(h1t); if (normrun>0) c.draw(h1tr,"same");
    	    	
    }
    
    public void plotTLSummary(String tab) {
        EmbeddedCanvas c = getCanvas(tab); c.clear(); c.divide(3, 6);
        
        int as = getActiveSCAL();
        int is = getActiveSector();  
    	int st = as +(dNorm ? 1:0); 
    	int n = 0;
    	for (int il = 0; il<3 ; il++) {
    		for (int in = 0; in<2; in++) {
    			for (int iv = 0; iv<3; iv++) {
    				int sl = iv+3*il+1;
    				int hl = 10*is + sl;
    				c.cd(n); c.getPad().setTitleFontSize(24); dgm.draw("TIMELINE", c, hl, st, in); n++;
    			}
    		}
    	}
    }  
    	
    public void analyzeNORM(String detName, int is1, int is2) {
        System.out.println("ECscaler:analyzeNORM("+detName+","+is1+","+is2+")");
    	analyzeATData(detName,is1,is2);
    }
    
    public void analyzeATData(String detName, int is1, int is2) {
        System.out.println("ECscaler:analyzeATData("+detName+","+is1+","+is2+")");
        ATData.clear();
        for (int is=is1; is<is2; is++) {
        	for (int sl=1; sl<layMap.get(detName).length+1 ; sl++) { int hl=10*is+sl;
    		ATData.add(new ArrayList<H1F>(),is,sl, 0); ATData.getItem(is,sl, 0).addAll(dgm.getH2F("ADC"+hl).getSlicesX());
    		ATData.add(new ArrayList<H1F>(),is,sl, 1); ATData.getItem(is,sl, 1).addAll(dgm.getH2F("TDC"+hl).getSlicesX());
    		ATData.add(new ArrayList<H1F>(),is,sl,20); ATData.getItem(is,sl,20).addAll(dgm.getH2F("VADC"+hl).getSlicesX());
    		ATData.add(new ArrayList<H1F>(),is,sl,21); ATData.getItem(is,sl,21).addAll(dgm.getH2F("VTDC"+hl).getSlicesX());
        	}	
        }
		ATData.add(new ArrayList<H1F>(),0,0,40);  ATData.getItem(0,0,40).addAll(dgm.getH2F("TRIG").getSlicesX());
        
    }
    
    public void getATNData(String detName, int is1, int is2) {
        System.out.println("ECscaler:getATNData("+detName+","+is1+","+is2+")");
    	int i1=normrun, i2=i1+normrng-1;
    	NATData.clear();
    	for (int is=is1; is<is2; is++) {
    		for (int sl=1; sl<layMap.get(detName).length+1 ; sl++) {
    			NATData.add(new H1F(),is,sl, 0); NATData.add(sumSlices(ATData.getItem(is,sl, 0),i1,i2),is,sl, 0); 
    			NATData.add(new H1F(),is,sl, 1); NATData.add(sumSlices(ATData.getItem(is,sl, 1),i1,i2),is,sl, 1);     	
    			NATData.add(new H1F(),is,sl,20); NATData.add(sumSlices(ATData.getItem(is,sl,20),i1,i2),is,sl,20); 
    			NATData.add(new H1F(),is,sl,21); NATData.add(sumSlices(ATData.getItem(is,sl,21),i1,i2),is,sl,21);     	
    		}	
    	}
		NATData.add(new H1F(),0,0,40);  NATData.add(sumSlices(ATData.getItem(0,0,40),i1,i2),0,0,40);
    } 
    
    public float getNorm(int at, int is, int sl, int ic) {
    	if(at==0 && normrun==0) return anorm.get(is, sl, ic)/occHL;
    	if(at==1 && normrun==0) return tnorm.get(is, sl, ic)/occHL;
    	if(normrun>0)  return (float) NATData.getItem(is,sl,at).getBinContent(ic-1);
    	return 0f;
    }
    
    @Override
    public void NormRunFunction() {
    	getATNData("ECAL",1,7);
    	fillHistFromFifo("ECAL",1,7);    	
    }
    
    public H1F sumSlices(ArrayList<H1F> list, int i1, int i2) {
    	H1F h = list.get(i1-1).histClone("dum");
    	for (int i=i1; i<i2; i++) h.add(list.get(i));
    	h.normalize(i2-i1+1);
    	return h;
    }
  
	
}