package org.clas.analysis;

import java.awt.Dimension;
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
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.utils.groups.IndexedList;


public class ECscaler extends DetectorMonitor {

    static int              occCounts = 0;
    static int                 occMax = 10002;
    static int                 occLo  =  1;
    static int                 occHi  =  100;
    static int                 occHL  = occHi-occLo+1;
    static int                    nev = 1;
    static int               evn_last = 0;
    static long              tim_last = 0;
    static int                prevRun = 0;
    static DataEvent        prevEvent = null;
    static boolean          singleRun = false;
    
    static DetectorOccupancy               occupancyECAL = new DetectorOccupancy();    
    DetectorCollection<LinkedList<Integer>> fifoa = new DetectorCollection<LinkedList<Integer>>();
    DetectorCollection<LinkedList<Integer>> fifot = new DetectorCollection<LinkedList<Integer>>();
    DetectorCollection<Integer>             anorm = new DetectorCollection<Integer>();
    DetectorCollection<Integer>             tnorm = new DetectorCollection<Integer>();
    DetectorCollection<Float>                asum = new DetectorCollection<Float>();
    DetectorCollection<Float>                tsum = new DetectorCollection<Float>(); 
    List<Integer>                         evnlist = new ArrayList<Integer>();
    List<Integer>                         evrlist = new ArrayList<Integer>();
    
    IndexedList<ArrayList<H1F>>            AData = new IndexedList<ArrayList<H1F>>(2);  
    IndexedList<ArrayList<H1F>>            TData = new IndexedList<>(2);  
    IndexedList<H1F>                      ANData = new IndexedList<>(2);  
    IndexedList<H1F>                      TNData = new IndexedList<>(2);  
    
    HipoDataSync  writer = null;		
    
    int[]           npmt = {68,62,62,36,36,36,36,36,36};    
    String[]         det = new String[]{" PCAL"," ECIN"," ECOU"};
    String[]           v = new String[]{"U","V","W"};
    
    double aYL,aYS,aYR,tYL,tYS,tYR;
    
    public ECscaler(String name) {
        super(name);
        
        dgmActive=true; 
        setDetectorTabNames("RAW","NORM","TimeLine","STATUS");

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
    	System.out.println("ECstatus:createHistos("+run+")");
    	if(dumpFiles) openOutput(outPath+"ECscaler/ECscaler-"+run+".hipo");
    	setRunNumber(run);    	   	
    	histosExist = true;  
    	TLmax = getTotalEvents()>occMax?getTotalEvents()/occMax:getTotalEvents();     	
    	createRAW(TLmax);
    	createNORM(TLmax); 
    	createSTATUS(run);
    }
    
    public void createRAW(int nx) {    	
    	for(int is=1; is<7; is++) { //sectors 1-6
    		for(int im=1; im<4; im++) { //modules pcal,ecin,ecou
    			for(int iv=1; iv<4; iv++) { //views u,v,w
    				int sl = iv+3*(im-1); //superlayer 1-9
    	    		dgm.add("RAW",1,2, is, sl, getRunNumber());
    	    		int ny=npmt[sl-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";  
    	    		dgm.makeH2("ADC"+is+sl, nx,1,nx+1,ny,1,ny+1, -1," ADC COUNTS", "RUN INDEX", tity);    	    		
    	    		dgm.makeH2("TDC"+is+sl, nx,1,nx+1,ny,1,ny+1, -1, "TDC COUNTS", "RUN INDEX", tity);    	    		
    			}
    		}
    	}
    }
    
    public void createNORM(int nx) {    	
    	for(int is=1; is<7; is++) {
    		for(int im=1; im<4; im++) {
    			for(int iv=1; iv<4; iv++) {
    				int sl = iv+3*(im-1);
    	    		dgm.add("NORM",1,2, is, sl, getRunNumber());
    	    		int ny=npmt[sl-1]; String tity="SEC"+is+det[im-1]+" "+v[iv-1]+" PMT";
    	    		dgm.makeH2("NADC"+is+sl, nx,1,nx+1,ny,1,ny+1, -1,"NORM ADC COUNTS", "RUN INDEX", tity);
    	    		dgm.cc("NADC"+is+sl, false, false, 0, 0, -3, 3);
    	    		dgm.makeH2("NTDC"+is+sl, nx,1,nx+1,ny,1,ny+1, -1,"NORM TDC COUNTS", "RUN INDEX", tity);
    	    		dgm.cc("NTDC"+is+sl, false, false, 0, 0, -3, 3);
    			}
    		}
    	}
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
    	plot("RAW");
    	plot("NORM");
    	plotStatus("STATUS");
    }
    
    @Override
    public void plotScalers(int run) {
    	setRunNumber(run);
    	plotTimeLine("TimeLine");
    }
    
	public void plot(String tabname) { 
		if(tabname=="Status") {dgm.drawGroup(tabname,0, 0, getRunNumber()); return;}
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
    	analyzeSTATUS("ECAL",1,7);
    	analyzeNORM("ECAL",1,7);
    	if(dumpFiles) writer.close();
    	isAnalyzeDone = true;
    }	
    
    @Override
    public void processEvent(DataEvent de) {	    	
    	prevRun = getRunNumber(); 
    	int run = de.hasBank("RUN::config") ? de.getBank("RUN::config").getInt("run", 0):-1;
    	if(run<1) return;
    	
    	setRunNumber(run);
    	singleRun = prevRun==run;
    	
        if(de.hasBank("ECAL::scaler")) {fillFifoFromBank(de);return;}
        
//        if(getElecTrigger()!=1) return;
        
        if(occCounts>=occMax || run!=prevRun) {
        	fillFifoFromData();
        	if(dumpFiles) {fillBankFromData(prevEvent); writer.writeEvent(prevEvent);}
        	occupancyECAL.reset(); occCounts = 0;
        }      
        if(de.hasBank("ECAL::adc")) occupancyECAL.addADCBank(de.getBank("ECAL::adc"));
        if(de.hasBank("ECAL::tdc")) occupancyECAL.addTDCBank(de.getBank("ECAL::tdc"));
        
        occCounts++;
        prevEvent = de;
    }
    
    public boolean inNormWindow() {
    	return occCounts>=occLo && occCounts<=occHi;
    }
    
    public void fillFifoFromBank(DataEvent de) {
    	
    	int  evn = de.getBank("RUN::config").getInt("event", 0);
    	long tim = de.getBank("RUN::config").getLong("timestamp",0);
    	
    	int ev_rate = (int) (0.25e9*(evn-evn_last)/(tim-tim_last));
        runlist.add(occCounts,getRunNumber()); 
        evnlist.add(occCounts,evn);evrlist.add(occCounts,evn_last>0?ev_rate:0); 
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
    
    public void fillFifoFromData() {
    	
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
        System.out.println("ECscaler:fillHistFromFifo("+detName+","+is1+","+is2+")");
    	for (int is=is1; is<is2 ; is++) {
    		for (int il=1; il<layMap.get(detName).length+1 ; il++) {
				dgm.getH2F("ADC"+is+il).reset();  dgm.getH2F("TDC"+is+il).reset(); 
				dgm.getH2F("NADC"+is+il).reset(); dgm.getH2F("NTDC"+is+il).reset(); 
    			for (int ic=1; ic<nlayMap.get(detName)[il-1]+1; ic++) {    	
            		Integer fa[] = new Integer[fifoa.get(is, il, ic).size()];
     				Integer ft[] = new Integer[fifot.get(is, il, ic).size()];
    				fifoa.get(is, il, ic).toArray(fa);
    				fifot.get(is, il, ic).toArray(ft);
    				for (int it=0; it<fa.length; it++) {
    					float y = (float)((float)(fa[it]-getNorm(0,is,il,ic))/Math.sqrt(fa[it]));
    					dgm.fill( "ADC"+is+il,it+1,ic,fa[it]);
    					dgm.fill("NADC"+is+il,it+1,ic,y);
    				}
    				for (int it=0; it<ft.length; it++) {
    					float y = (float)((float)(ft[it]-getNorm(1,is,il,ic))/Math.sqrt(ft[it]));
    					dgm.fill( "TDC"+is+il,it+1,ic,ft[it]);
    					dgm.fill("NTDC"+is+il,it+1,ic,y);
    				}
    			}
    		}
    	}
    	System.out.println("return");
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
                    asum.add(is,sl,ip,(float)dgm.getH2F("ADC"+is+sl).sliceY(ip-1).integral());
                    tsum.add(is,sl,ip,(float)dgm.getH2F("TDC"+is+sl).sliceY(ip-1).integral());
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
                	int    sl = il+(im-1)*3;  
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
    	
        EmbeddedCanvas c = getCanvas(tab); c.clear(); c.divide(2, 2);
                
        int is = getActiveSector(), sl = getActiveView()+3*getActiveLayer()+1;  
        
    	DataLine line3 = new DataLine(runIndexSlider,  1,  runIndexSlider,  npmt[sl-1]+1);  line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,  runIndexSlider+1,npmt[sl-1]+1);  line4.setLineColor(5);  
    	
    	String tit1 = "RUN "+runlist.get(runIndexSlider-1)+"   EVENT "+evnlist.get(runIndexSlider-1);
    	String tit2 = "   EV/SEC "+evrlist.get(runIndexSlider-1);
    	
    	H1F h1a = AData.getItem(is,sl).get(runIndexSlider-1); h1a.setTitle(tit1+(singleRun?tit2:" ")); h1a.setFillColor(21); h1a.setOptStat("1000000");
    	H1F h1t = TData.getItem(is,sl).get(runIndexSlider-1); h1t.setTitle(" ");                       h1t.setFillColor(21); h1t.setOptStat("1000000");
    	
    	float amax = (float) h1a.getMax()*1.3f, tmax = (float) h1t.getMax()*1.3f;
    	
    	H1F h1ar=new H1F(),h1tr=new H1F();   
    	
    	if (normrun>0) {    		
    		h1ar = ANData.getItem(is,sl); h1ar.setLineColor(3); h1ar.setLineWidth(5); amax = (float) h1ar.getMax()*1.3f; h1ar.setOptStat("1000000");
    		h1tr = TNData.getItem(is,sl); h1tr.setLineColor(3); h1tr.setLineWidth(5); tmax = (float) h1tr.getMax()*1.3f; h1tr.setOptStat("1000000");
    	}
      	
    	c.cd(0); c.getPad().setTitleFontSize(24); c.draw(dgm.getH2F("ADC"+is+sl)); c.draw(line3);c.draw(line4);
    	c.cd(1); c.getPad().setTitleFontSize(24); c.getPad().getAxisY().setRange(0, amax); 
    	if(h1a.getIntegral()>0) c.draw(h1a); if (normrun>0) c.draw(h1ar,"same");
    	c.cd(2); c.getPad().setTitleFontSize(24); c.draw(dgm.getH2F("TDC"+is+sl)); c.draw(line3);c.draw(line4);
    	c.cd(3); c.getPad().setTitleFontSize(24); c.getPad().getAxisY().setRange(0, tmax); 
    	if(h1t.getIntegral()>0) c.draw(h1t); if (normrun>0) c.draw(h1tr,"same");
    	    	
    }
    
    public void analyzeNORM(String detName, int is1, int is2) {
        System.out.println("ECscaler:analyzeNORM("+detName+","+is1+","+is2+")");
    	analyzeATData(detName,is1,is2);
    }
    
    public void analyzeATData(String detName, int is1, int is2) {
        System.out.println("ECscaler:analyzeATData("+detName+","+is1+","+is2+")");
        AData.clear(); TData.clear();
        for (int is=is1; is<is2; is++) {  
			for (int sl=1; sl<layMap.get(detName).length+1 ; sl++) {
    		    AData.add(new ArrayList<H1F>(),is,sl);  AData.getItem(is,sl).addAll(dgm.getH2F("ADC"+is+sl).getSlicesX());
    		    TData.add(new ArrayList<H1F>(),is,sl);  TData.getItem(is,sl).addAll(dgm.getH2F("TDC"+is+sl).getSlicesX());    	
    		}
		}
    }
    public void getATNData(String detName, int is1, int is2) {
        System.out.println("ECscaler:getATNData("+detName+","+is1+","+is2+")");
    	int i1=normrun, i2=i1+normrng-1;
    	ANData.clear(); TNData.clear();
		for (int is=is1; is<is2; is++) {  
			for (int sl=1; sl<layMap.get(detName).length+1 ; sl++) {
    			ANData.add(new H1F(),is,sl); ANData.add(sumSlices(AData.getItem(is,sl),i1,i2),is,sl); 
    			TNData.add(new H1F(),is,sl); TNData.add(sumSlices(TData.getItem(is,sl),i1,i2),is,sl);     	
    		}
		}
    } 
    
    public float getNorm(int at, int is, int sl, int ic) {
    	if(at==0 && normrun==0) return anorm.get(is, sl, ic)/occHL;
    	if(at==1 && normrun==0) return tnorm.get(is, sl, ic)/occHL;
    	if(at==0 && normrun>0)  return (float) ANData.getItem(is,sl).getBinContent(ic-1);
    	if(at==1 && normrun>0)  return (float) TNData.getItem(is,sl).getBinContent(ic-1);
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