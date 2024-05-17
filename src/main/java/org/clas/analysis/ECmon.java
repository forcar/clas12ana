package org.clas.analysis;

import org.clas.service.ec.ECmc;
import org.clas.tools.FitData;
import org.clas.viewer.DetectorMonitor;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ECmon extends DetectorMonitor {
	
	IndexedList<List<Float>> tdcs = new IndexedList<List<Float>>(3);
    public static double[] SCALE  = {10,10,10,10,10,10,10,10,10}; // Fitter.ADC/SCALE is plotted and fitted in ECMon
    public static double[] SCALE5 = {10,10,10,5,5,5,5,5,5};       // Sector 5 ECAL uses EMI PMTs near max voltage (not used in ECmon!)
    IndexedTable time=null, ftime=null, gtw=null, fo=null, fgo=null, tmf=null, tmfcut=null, tgo=null, gain=null, veff=null, tt=null, fadc=null, rfT=null;
	float FTOFFSET, TOFFSET, TMFCUT;
    static float tps = (float) 0.02345;
    String cut1 = "ABS(TDC-FTDC)<10 ", cut2 = " 65<FTDC<120";
	String[] det = {"PCAL","ECIN","ECOU"}; 
	String[] uvw = {"U","V","W"};
	int bcut[] = {60,125}, bcutMC[] = {100,200};
    int[]   npmt = {68,62,62,36,36,36,36,36,36}; 
    int nc = 0;
    
    ECmc mc = new ECmc(); 
    
    public ECmon(String name) {
        super(name);
    	
        dgmActive=true; 
        this.setDetectorTabNames("AT",
        		                 "DEFF",
        		                 "PED",
        		                 "RMS",
        		                 "DEFF2",
        		                 "AT4");

        this.useCALUVWSECButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();       
    }

    public void localinit() {
        System.out.println(getDetectorName()+".localinit()"); 
        tl.setFitData(Fits);
    } 
    
    public void localclear() {
    	System.out.println(getDetectorName()+".localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	FitSummary.clear();
    	Fits.clear();
    	tl.fitData.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
    public void init(int run) {
	    setRunNumber(run);
	    runlist.add(run);     	
    }
    
    @Override  
    public void createHistos(int run) {
    	histosExist = true;
    	System.out.println(getDetectorName()+".createHistos("+run+")");
	    init(run);
	    createAT(0);
	    createDEFF(0);
	    createPED(0);
	    createRMS(0);
    }
    
    public void createAT(int st) {
	
    	for(int is=1; is<7; is++) {
    	for(int id=0; id<3; id++) {String tit = "Sector "+is+" "+det[id]+" ";
    	
    	int np = id==0?68:36; 
    	
    	dgm.add("AT",3,6,10*is+id,st,getRunNumber());    	
        for(int iv=0; iv<3; iv++) dgm.makeH2("at0"+is+id+iv,100, 0,150,100,0,240,-1,tit+uvw[iv]+" STRIPS","FADC","FTDC");
        for(int iv=0; iv<3; iv++) dgm.makeH2("at1"+is+id+iv,100, 0,150,100,0,240,-1,cut1,"FADC","TDC");
        for(int iv=0; iv<3; iv++) dgm.makeH2("at3"+is+id+iv,100, 0,150,100,-11,11,-1,cut1+cut2,"FADC","TDC-FTDC");
        for(int iv=0; iv<3; iv++) dgm.makeH2("at2"+is+id+iv,100, 0,150,100,0,240,-1," ","FADC","TDC");
        for(int iv=0; iv<3; iv++) dgm.makeH2("at4"+is+id+iv,100, 0,150,np,1,np+1,-1," ","FADC","PMT");
        for(int iv=0; iv<3; iv++) dgm.makeH2("at5"+is+id+iv,100, 0,150,np,1,np+1,-1," ","FADC","PMT");
        
    	}
    	}   	
    }
    
    public void createDEFF(int st) {

    	for(int is=1; is<7; is++) {
    	for(int id=0; id<3; id++) {String tit = "Sector "+is+" "+det[id]+" ";   
    		
        int np = id==0?68:36;
        
        dgm.add("DEFF",3,2,10*is+id,st,getRunNumber());               
        for(int iv=0; iv<3; iv++) {dgm.makeGE("ef1"+is+id+iv,0,cut1+cut2,"FADC","TDC EFFICIENCY");
                                 //dgm.makeF1D("fe1"+is+id+iv, "p0",0,150,1.0);
                                   dgm.cc("ef1"+is+id+iv,false,false,0,1.05f,0,0);} 
        
        for(int iv=0; iv<3; iv++) {dgm.makeH2("ef2"+is+id+iv,100,0,150,np,1,np+1,-1," ","FADC","PMT");                                  
                                   dgm.cc("ef2"+is+id+iv,false,false,0,0,0.1f,1.0f);}            
    	}
    	}   	
    }
    
    public void createPED(int st) {
    	
        for (int is=1; is<7 ; is++) {
        for (int id=0; id<3 ; id++) {String tit = "Sector "+is+" "+det[id]+" ";
        
        dgm.add("PED",3,3,10*is+id,st,getRunNumber());
        for (int iv=0; iv<3 ; iv++) { int sl = 3*id+iv; int np=npmt[sl];       
        dgm.makeH2("ped1"+is+id+iv,110,-10,100,np,1,np+1,-1,tit+uvw[iv]+" STRIPS ","PED (Measured-REF)","PMT");}

        for (int iv=0; iv<3 ; iv++) { int sl = 3*id+iv; int np=npmt[sl];
        dgm.makeH1("ped20"+is+id+iv,np,1,np+1,-1," ","PMT");
        dgm.makeH1("ped21"+is+id+iv,np,1,np+1,-2," ","PMT");
        dgm.makeH1("ped22"+is+id+iv,np,1,np+1,-2," ","PMT");
        dgm.makeH1("ped23"+is+id+iv,np,1,np+1,-2," ","PMT");       
        dgm.makeH1("ped24"+is+id+iv,np,1,np+1,-2," ","PMT");} 
        
        for (int iv=0; iv<3 ; iv++) { int sl = 3*id+iv;        
        dgm.makeGE("gped21"+is+id+iv,-1,"BLK:6 RED:15 GRN:30 BLU:50","PMT","BAD PEDESTAL FRACTION");
        dgm.cc("gped21"+is+id+iv,false,false,0f,0.3f,0,0);
        dgm.makeGE("gped22"+is+id+iv,-2,"","PMT","BAD PEDESTAL FRACTION");
        dgm.makeGE("gped23"+is+id+iv,-2,"","PMT","BAD PEDESTAL FRACTION");
        dgm.makeGE("gped24"+is+id+iv,-2,"","PMT","BAD PEDESTAL FRACTION");}          
        
        } 
        }
    }   
    
    public void createRMS(int st) {

    	for (int is=1; is<7 ; is++) {
        for (int id=0; id<3 ; id++) {String tit = "Sector "+is+" "+det[id]+" ";
        
        dgm.add("RMS",3,2,10*is+id,st,getRunNumber());
        for (int iv=0; iv<3 ; iv++) { int sl = 3*id+iv; int np=npmt[sl];
        dgm.makeH2("rms1"+is+id+iv,50,0,10,np,1,np+1,-1," ","RMS","PMT");} 
        
        }
        }
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plotSummary(run);
       	plotAnalysis(run);   	
    }
    
    public void plotSummary(int run) {
        setRunNumber(run);
        plot("AT");
        plot("PED");
        plot("RMS");
    }
    
    public void plotAnalysis(int run) {
        setRunNumber(run);
    	if(!isAnalyzeDone) return;  
     	plot("PED");
    	plot("DEFF");
    	plotUVW("DEFF2",6);
    	plotUVW("AT4",12);

       	if(dumpFiles) dumpFiles();
   }
        
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,10*getActiveSector()+getActiveLayer(),0,getRunNumber());
    }
    
    public void plotUVW(String name, int isoff) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(name);
        
        int    is = getActiveSector(); 
        int     i = getActiveLayer();
        int     j = getActiveView();
        
        H1F h; FitData fd; 
        int    np = npmt[i*3+j];
        String det[] = {" PCAL"," ECIN"," ECOU"}, uvw[] = {" U"," V"," W"};
        
        c.clear();
        if (np==36) c.divide(6,6); if (np==68) c.divide(9,8); if (np==62) c.divide(8,8); 
        
        for (int ip=0; ip<np ; ip++) {
            c.cd(ip); c.getPad(ip).getAxisY().setLog(false); c.setAxisLabelSize(14); c.setAxisTitleSize(14);
            fd = tl.fitData.getItem(is+isoff,i*3+j,ip+1,getRunNumber()); 
            h=fd.getHist();  h.setOptStat("0"); h.setTitleX("Sector "+is+det[i]+(j<3?uvw[j]+(ip+1):" ")+" FADC"); h.setTitle(" ");
            fd.getGraph().getFunction().setLineColor(2);
            c.draw(h); c.draw(fd.getGraph().getFunction(),"same");
       }
    }
    
    public void dumpFiles() {
//    	writeFile("ped",0,0,0,0,0,0);
    	writeFile("deff",1,7,0,3,0,3);
    	writeFile("at4",1,7,0,3,0,3);
    }
    
    @Override         
    public void plotEvent(DataEvent event) {
    	analyze();
    }
    
    @Override   
    public void initCCDB(int runno) {
    	System.out.println(getDetectorName()+".initCCDB("+runno+") "); 
    	fo      = cm.getConstants(runno, "/calibration/ec/fadc_offset");        //Crate/fiber FADC offsets 
        fgo     = cm.getConstants(runno, "/calibration/ec/fadc_global_offset"); //FADC capture window 
        gtw     = cm.getConstants(runno, "/calibration/ec/global_time_walk");   //Global time walk correction using raw ADC
        tmfcut  = cm.getConstants(runno, "/calibration/ec/tmf_window");
        tmf     = cm.getConstants(runno, "/calibration/ec/tmf_offset");         //TDC-FADC offsets
        time    = cm.getConstants(runno, "/calibration/ec/timing");
        ftime   = cm.getConstants(runno, "/calibration/ec/ftiming"); 
        fadc    = cm.getConstants(runno, "/daq/fadc/ec"); 
        tt      = cm.getConstants(runno, "/daq/tt/ec");
        
        FTOFFSET = (float) fgo.getDoubleValue("global_offset",0,0,0);
        TMFCUT   = (float) tmfcut.getDoubleValue("window", 0,0,0); //acceptance window for TDC-FADC cut
        
        getReverseTT(cm,runno,"/daq/tt/ec");
        
	    if (useGPP) mc.initCCDB(runno,cm);
    }
    
    public void getDeff(int is1, int is2, int id1, int id2) {
       
       for (int is=is1; is<is2; is++) {
       for (int id=id1; id<id2; id++) {
       for (int il=0; il<3; il++) {
    	   String tag = ""+is+id+il;
    	   H2F h2at2 = dgm.getH2F("at2"+tag);
    	   H2F h2at4 = dgm.getH2F("at4"+tag);
    	   H2F h2at5 = dgm.getH2F("at5"+tag);
    	   H1F h1at2 = projectionX(h2at2,50,240);
    	   H1F h1at4 = h2at4.projectionX();

    	   dgm.getGE("ef1"+tag).reset();
    	   dgm.getGE("ef1"+tag).copy(H1F.divide(h1at2, h1at4).getGraph());
    	   dgm.getGE("ef1"+tag).setMarkerColor(1); dgm.getGE("ef1"+tag).setLineColor(2);
    	    
    	   dgm.getH2F("ef2"+tag).reset();
    	   dgm.getH2F("ef2"+tag).add(H2F.divide(h2at5, h2at4));
    	   
    	   H1F h1p20 = dgm.getH1F("ped20"+tag);
    	   H1F h1p21 = dgm.getH1F("ped21"+tag);
    	   H1F h1p22 = dgm.getH1F("ped22"+tag);
    	   H1F h1p23 = dgm.getH1F("ped23"+tag);
    	   H1F h1p24 = dgm.getH1F("ped24"+tag);
    	   
    	   dgm.getGE("gped21"+tag).reset(); dgm.getGE("gped21"+tag).copy(H1F.divide(h1p21, h1p20).getGraph());
    	   dgm.getGE("gped22"+tag).reset(); dgm.getGE("gped22"+tag).copy(H1F.divide(h1p22, h1p20).getGraph());
    	   dgm.getGE("gped23"+tag).reset(); dgm.getGE("gped23"+tag).copy(H1F.divide(h1p23, h1p20).getGraph());
    	   dgm.getGE("gped24"+tag).reset(); dgm.getGE("gped24"+tag).copy(H1F.divide(h1p24, h1p20).getGraph());
    	   
    	   dgm.getGE("gped21"+tag).setMarkerColor(1); dgm.getGE("gped21"+tag).setLineColor(1);
    	   dgm.getGE("gped22"+tag).setMarkerColor(2); dgm.getGE("gped22"+tag).setLineColor(2);
    	   dgm.getGE("gped23"+tag).setMarkerColor(3); dgm.getGE("gped23"+tag).setLineColor(3);
    	   dgm.getGE("gped24"+tag).setMarkerColor(4); dgm.getGE("gped24"+tag).setLineColor(4);    	   
       }
       }
       }
    } 
        
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() "); 
    	getDeff(1,7,0,3);    	
    	FitUVW ef1=new FitUVW("ef1",19,0,true); 
    	FitUVW ef2=new FitUVW("ef2",19,6,false); 
    	FitUVW at4=new FitUVW("at4",20,12,false);
    	ef1.setx(10,140,10,140); ef1.doFit("INIT",1,7,0,3,0,3); ef1.doFit("FIT",1,7,0,3,0,3);
    	ef2.setx(10,140,10,140); ef2.doFit("INIT",1,7,0,3,0,3); ef2.doFit("FIT",1,7,0,3,0,3);
    	at4.setx(8,60,8,60);     at4.doFit("INIT",1,7,0,1,0,3); at4.doFit("FIT",1,7,0,1,0,3); at4.ff=21;
    	at4.setx(8,30,8,30);     at4.doFit("INIT",1,7,1,3,0,3); at4.doFit("FIT",1,7,1,3,0,3);   	
    	
    	isAnalyzeDone = true;
    }
    
    public class FitUVW {
    	
    	String DataSet, tag;
    	int ff=0, off=0, run=getRunNumber(); float[] x; 
   	    boolean SLV=false, SLVP=false, INIT=false, FIT=false; 	    
    	
    	public FitUVW(String name, int val1, int val2, boolean slv) {
    		DataSet=name;  ff=val1; off=val2;  SLV=slv; SLVP=!slv; 		
    	}
    	
    	public void setx(float...lim) {
    		x=lim;
    	}
    	
    	public void doFit(String cmd, int...lim) {
    		switch (cmd) {
    		case "INIT": INIT=true;   FIT=false; break;
    		case  "FIT": INIT=false;  FIT=true;
    		}
    		processFit(lim);
    	}

        public void processFit(int[] lim) {  
            for (int is=lim[0]; is<lim[1]; is++) {
            	System.out.println("FitUVW.processFit ff="+ff+" Sector "+is);
            	for (int il=lim[2]; il<lim[3]; il++) {
            		for (int iv=lim[4]; iv<lim[5]; iv++) {
            			if(INIT) init(is,il,iv);
            			if(FIT)   fit(is,il,iv);
                    }        		
            	}        	
            }                  
        }
            
        public void init(int is, int il, int iv) {
        	tag = DataSet+""+is+il+iv; 
        	if(SLV)                                      tl.fitData.add(initFit(dgm.getGE(tag),           ff,x[0],x[1],x[2],x[3]),is+off,il,     iv, run);
        	if(SLVP) for (int i=0; i<npmt[il*3+iv]; i++) tl.fitData.add(initFit(dgm.getH2F(tag).sliceY(i),ff,x[0],x[1],x[2],x[3]),is+off,il*3+iv,i+1,run);
        }       
        
        public void fit(int is, int il, int iv) {
        	if(SLV &&                                       goodFit(is+off,il,iv))       fitGraph(is+off,il,iv);
        	if(SLVP) for (int i=0; i<npmt[il*3+iv]; i++) if(goodFit(is+off,il*3+iv,i+1)) fitGraph(is+off,il*3+iv,i+1);
        }
        
        public boolean goodFit(int a, int b, int c) {return tl.fitData.hasItem(a,b,c,run);}
        public void   fitGraph(int a, int b, int c) {       tl.fitData.getItem(a,b,c,run).fitGraph("",cfitEnable,fitVerbose);}
    
        public FitData initFit(GraphErrors g, int ff, double pmin, double pmax, double fmin, double fmax) {
            FitData fd = new FitData(g);        
            fd.initFit(ff,pmin,pmax,fmin,fmax); 
            fd.doFit = g.getDataSize(0)==0 ? false:true;
            return fd;
        }  

        public FitData initFit(H1F h, int ff, double pmin, double pmax, double fmin, double fmax) {
            FitData fd = new FitData(h.getGraph());        
            fd.setInt((int)h.getIntegral());
            String tit = h.getTitle();
            fd.setHist(h);
            fd.graph.getAttributes().setTitleX(h.getTitleX()); 
            fd.hist.getAttributes().setTitleX(h.getTitleX());
            fd.hist.setTitle(tit);
            fd.graph.setTitle(tit);
            fd.initFit(ff,pmin,pmax,fmin,fmax);
            return fd;
        }
        
    }                 

    public int getDet(int layer) {
        int[] il = {0,0,0,1,1,1,2,2,2}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
     }
    
    public int getLay(int layer) {
        int[] il = {1,2,3,1,2,3,1,2,3}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
    }   
    
    @Override     
    public void processEvent(DataEvent event) {

        tdcs.clear();
        
        boolean isCosmic = getRunNumber()==11581;
        boolean     isMC = event.hasBank("MC::Particle");
        
        if (event.hasBank("ECAL::tdc")==true){
            DataBank  bank = event.getBank("ECAL::tdc");
            for(int i = 0; i < bank.rows(); i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  ip = bank.getShort("component",i);               
                float tdc = bank.getInt("TDC",i);

                if (tdc>0) {
                    if (!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Float>(),is,il,ip);
                         tdcs.getItem(is,il,ip).add(tdc);       
                }
            }
        }
        
        if (event.hasBank("ECAL::adc")==true) {
            DataBank  bank = event.getBank("ECAL::adc");
            for(int i = 0; i < bank.rows(); i++){
                int     is = bank.getByte("sector",i);
                int     il = bank.getByte("layer",i);
                int     ip = bank.getShort("component",i);
                int    adc = Math.abs(bank.getInt("ADC",i));
                int    ped = bank.getShort("ped",i);
                float tadc = bank.getFloat("time", i); 
               
                if (useGPP) {                	
                    mc.setSLP(is, il, ip); 
                	mc.addADC(adc); mc.addFTDC(tadc);
                	adc = (int) mc.dgtzFADC(); tadc = (float) mc.dgtzFTDC();
                }                
                
                if (adc==0) continue;
                
                float t = tadc + (float) tmf.getDoubleValue("offset",is,il,ip) // TDC-FADC offset (sector, layer, PMT)
                               + (float)  fo.getDoubleValue("offset",is,il,0); // TDC-FADC offset (sector, layer) 
               
                boolean bcut2 = t > (isMC ? bcutMC[0]:bcut[0]) && t < (isMC ? bcutMC[1]:bcut[1]);
                
                String tag = ""+is+getDet(il)+(getLay(il)-1);
                
                float  sca = (float) SCALE[il-1];
                float sadc = (float)(adc/sca); 
                
                float  tmax = 1000; float tdc = 0;
                
                float ftdc_corr = t+FTOFFSET;
                
                if (tdcs.hasItem(is,il,ip)) {
                    float radc = (float)Math.sqrt(adc);
                    for (float tdcc : tdcs.getItem(is,il,ip)) {
                   	    if (useGPP) {mc.addDTDC(tps*tdcc); tdcc = (float)mc.dgtzDTDC()/tps;}
                    	float tcor = tps*tdcc -  getTriggerPhase() - (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - FTOFFSET;
                        float tdif = tcor - t; 
                        if (Math.abs(tdif)<TMFCUT && tdif<tmax && bcut2) {tmax = tdif; tdc = tcor;} //here we only need tcor, not tdcc as in ECCommon
                        dgm.fill("at1"+tag,sadc,tcor); 
                    }
                } 
                
                if (sadc>0) {
                	dgm.fill("at0"+tag,sadc,t);
                	dgm.fill("at2"+tag,sadc,tdc);
                	dgm.fill("at3"+tag,sadc,tmax);
               	    if (bcut2)        dgm.fill("at4"+tag,sadc,ip);
                    if (bcut2&&tdc>0) dgm.fill("at5"+tag,sadc,ip);
                } 
                
                if (adc>(il>1?30:20) && rtt.hasItem(is,il,ip,0)) {
                    int[] dum = (int[]) rtt.getItem(is,il,ip,0);
                    int pdif = ped-(int) fadc.getIntValue("pedestal",dum[0],dum[1],dum[2]);
                    dgm.fill("ped1"+tag,pdif,ip);   
                    dgm.fill("ped20"+tag,ip);
                    if(pdif>6 ) dgm.fill("ped21"+tag,ip);
                    if(pdif>15) dgm.fill("ped22"+tag,ip);
                    if(pdif>30) dgm.fill("ped23"+tag,ip);
                    if(pdif>50) dgm.fill("ped24"+tag,ip);
                    dgm.fill("rms1"+tag,dgm.getH2F("ped1"+tag).getSlicesY().get(ip-1).getRMS(),ip);                   
                }
            }
        }
    }   

/* CALIBRATION FILES */ 

    @Override
	public void writeFile(String table, int is1, int is2, int il1, int il2, int iv1, int iv2) {
   	
		String path = filPath+"ccdb/";
		String line = new String();
		
		try { 
			File outputFile = new File(path+table+"_"+getRunNumber());
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println(getDetectorName()+".writefile("+table+")"+" RUN:"+getRunNumber());
			
			if (is1==0 && is2==0) { //daq tables
				
			int[] crate = {1,3,7,9,13,15,19,21,25,27,31,33};
			int[]  slot = {3,4,5,6,7,8,9,10,13,14,15,16,17,18};

	        for (int cr : crate) {
	        for (int sl : slot) {
	        for (int ch=0; ch<16; ch++) {							
	        	switch (table) {
			    case "ped": line = getPED(cr,sl,ch); break;
	        	}
	        	if(!skip(cr,sl,ch)) {
	        		outputBw.write(line);
	        		outputBw.newLine();
	        	}
	        }
	        }
	        }
	        
			} else { //calibration tables
				
			for (int is=is1; is<is2; is++) {
		    for (int il=il1; il<il2; il++ ) { //pcal,ecin,ecou
			for (int iv=iv1; iv<iv2; iv++) { //u,v,w
			for (int ip=0; ip<npmt[3*il+iv]; ip++) {
				switch (table) {
				case "deff": line = getDEFF(is,il,iv,ip,getRunNumber()); break; 
				case  "at4": line =  getAT4(is,il,iv,ip,getRunNumber()); 
				}
				outputBw.write(line);
				outputBw.newLine();
			}
			}
			}
			}
			
			}			
			outputBw.close();
			outputFw.close();
		}
		catch(IOException ex) {
			System.out.println("Error writing file '" );                   
			ex.printStackTrace();
		}

	}
	
	public String getDEFF(int is, int il, int iv, int ip, int run) {
		if(tl.fitData.hasItem(is+6,3*il+iv,ip+1,run)) {
			double deff0 = tl.fitData.getItem(is+6,3*il+iv,ip+1,run).p0;
			double deff1 = tl.fitData.getItem(is+6,3*il+iv,ip+1,run).p1;
			double deff2 = tl.fitData.getItem(is+6,3*il+iv,ip+1,run).p2;
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "+deff0+" "+deff1+" "+deff2;
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+"  0.0 0.0 0.0";		
		}
	}
	
	public String getAT4(int is, int il, int iv, int ip, int run) {
		
		if(tl.fitData.hasItem(is+12,3*il+iv,ip+1,run)) {
			double at4 = tl.fitData.getItem(is+12,3*il+iv,ip+1,run).p1;
		    return is+" "+(3*il+iv+1)+" "+(ip+1)+" "+at4;
		} else {
			return is+" "+(3*il+iv+1)+" "+(ip+1)+"  0.0";		
		}
	} 
  
    public boolean skip(int cr, int sl, int ch) {
    	boolean pc1 = cr==3||cr==9||cr==15||cr==21||cr==27||cr==33;
    	boolean pc2 = sl==17||sl==18;
    	return pc1 && pc2;
    }
    
    public double getDRM(int cr, int sl, int ch) {
    	
    	if(tt.hasEntry(cr,sl,ch)) {
    		int is = tt.getIntValue("sector",    cr,sl,ch);
            int il = tt.getIntValue("layer",     cr,sl,ch);
            int ip = tt.getIntValue("component", cr,sl,ch);
            String tag = ""+is+getDet(il)+(getLay(il)-1);
            return dgm.getH2F("ped1"+tag).getSlicesY().get(ip-1).getMean();
    	}
    	return 0;
    	
    }
    
	public String getPED(int cr, int sl, int ch) { 
		
		double D = fadc.getDoubleValue("pedestal",cr,sl,ch);
		int E    = fadc.getIntValue("nsb",cr,sl,ch);
		int F    = fadc.getIntValue("nsa",cr,sl,ch);
		int G    = fadc.getIntValue("tet",cr,sl,ch);
		int H    = fadc.getIntValue("window_offset",cr,sl,ch);
		int I    = fadc.getIntValue("window_size",cr,sl,ch);

		double DRM = getDRM(cr,sl,ch);

		return cr+" "+sl+" "+ch+" "+(D-DRM)+" "+E+" "+F+" "+G+" "+H+" "+I;				 
	}
	
    @Override
    public void timerUpdate() {   	
    	getDeff(getActiveSector(),getActiveSector()+1,getActiveLayer(),getActiveLayer()+1);
    } 
    
    public static double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
                                            t * ( 1.00002368 +
                                            t * ( 0.37409196 + 
                                            t * ( 0.09678418 + 
                                            t * (-0.18628806 + 
                                            t * ( 0.27886807 + 
                                            t * (-1.13520398 + 
                                            t * ( 1.48851587 + 
                                            t * (-0.82215223 + 
                                            t * ( 0.17087277))))))))));
        if (z >= 0) return  ans;
        else        return -ans;
    }
    
    public static double erfc(double z){
        return 1.0-erf(z);
    }
    
    public static double erf(double p0, double p1, double p2, double p3,double x){
        return p0 + p1*erfc((x-p2)/(p3*Math.sqrt(2.0)));
    }
    // fractional error less than x.xx * 10 ^ -4.
    // Algorithm 26.2.17 in Abromowitz and Stegun, Handbook of Mathematical.
    public static double erf2(double z) {
        double t = 1.0 / (1.0 + 0.47047 * Math.abs(z));
        double poly = t * (0.3480242 + t * (-0.0958798 + t * (0.7478556)));
        double ans = 1.0 - poly * Math.exp(-z*z);
        if (z >= 0) return  ans;
        else        return -ans;
    }

    public static double Phi(double z) {
        return 0.5 * (1.0 + erf(z / (Math.sqrt(2.0))));
    }
    
 //   public static void main(String[] args){ 
 //   	ECmon mon = new ECmon("ECmon");
 //   	mon.dumpFiles();0.70721
 //   }

/* 
    public static void main(String[] args){
        F1D func3 = new F1D("f1d","erf(x,0,0.70721)",-10,10);
        for(double x = func3.getMin(); x < func3.getMax(); x+= 0.05){
            double value = 2-func3.evaluate(x);
            double value2 = erf(x);
            double value3 = erf2(x);
            System.out.println("x = " + x + "  value = " + (value-1) +" "+value2+" "+value3);
            //func.show();
        }        
    }
    */
    
}
