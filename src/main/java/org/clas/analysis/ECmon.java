package org.clas.analysis;

import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ECmon extends DetectorMonitor {
	
	IndexedList<List<Float>> tdcs = new IndexedList<List<Float>>(3);
    public static double[] SCALE  = {10,10,10,10,10,10,10,10,10}; // Fitter.ADC/SCALE is plotted and fitted in ECMon
    public static double[] SCALE5 = {10,10,10,5,5,5,5,5,5};       // Sector 5 ECAL uses EMI PMTs near max voltage
    IndexedTable time=null, ftime=null, gtw=null, fo=null, fgo=null, tmf=null, tmfcut=null, tgo=null, gain=null, veff=null, tt=null, fadc=null, rfT=null;
	float FTOFFSET, TOFFSET, TMFCUT;
    static float tps = (float) 0.02345;
    String cut1 = "ABS(TDC-FTDC)<10 ", cut2 = " 65<FTDC<120";
	String[] det = {"PCAL","ECIN","ECOU"}; 
	String[] uvw = {"U","V","W"};
    int[]   npmt = {68,62,62,36,36,36,36,36,36}; 
    int nc = 0;
    
    public ECmon(String name) {
        super(name);
    	
        dgmActive=true; 
        this.setDetectorTabNames("AT",
        		                 "EFF",
        		                 "PED",
        		                 "RMS");
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
	    createEFF(0);
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
    
    public void createEFF(int st) {

    	for(int is=1; is<7; is++) {
    	for(int id=0; id<3; id++) {String tit = "Sector "+is+" "+det[id]+" ";   
    		
        int np = id==0?68:36;
        
        dgm.add("EFF",3,2,10*is+id,st,getRunNumber());               
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
    	if(dumpFiles) dumpFiles("ped");
    	plot("PED");
    	plot("EFF");
    }
    
    public void dumpFiles(String val) {
    	writeFile(val,0,0,0,0,0,0);
    }
    
    @Override         
    public void plotEvent(DataEvent event) {
    	analyze();
    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() "); 
    	if(dumpFiles) dumpFiles("ped");  
    	getDeff(1,7,0,3);
    	isAnalyzeDone = true;
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
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,10*getActiveSector()+getActiveLayer(),0,getRunNumber());
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
        
        if(event.hasBank("ECAL::tdc")==true){
            DataBank  bank = event.getBank("ECAL::tdc");
            for(int i = 0; i < bank.rows(); i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  ip = bank.getShort("component",i);               
                float tdc = bank.getInt("TDC",i)*tps-(isCosmic ? 293:0);
                if(tdc>0) {
                    if(!tdcs.hasItem(is,il,ip)) tdcs.add(new ArrayList<Float>(),is,il,ip);
                        tdcs.getItem(is,il,ip).add(tdc);       
                }
            }
        }
        
        if(event.hasBank("ECAL::adc")==true){
            DataBank  bank = event.getBank("ECAL::adc");
            for(int i = 0; i < bank.rows(); i++){
                int  is = bank.getByte("sector",i);
                int  il = bank.getByte("layer",i);
                int  ip = bank.getShort("component",i);
                int adc = Math.abs(bank.getInt("ADC",i));
                int ped = bank.getShort("ped",i);
                float t = bank.getFloat("time",i) + (float) tmf.getDoubleValue("offset",is,il,ip) 
                                                  + (float)  fo.getDoubleValue("offset",is,il,0); 
                
                boolean bcut2 = t > (isCosmic?110:65) && t < (isCosmic?170:120);
                
                String tag = ""+is+getDet(il)+(getLay(il)-1);
                
                float  sca = (float) ((is==5)?SCALE5[il-1]:SCALE[il-1]);
                float sadc = (float)(adc/sca); 
                
                float  tmax = 1000; float tdc = 0;
                
                double a0 = time.getDoubleValue("a0", is, il, ip); 
                double a2 = time.getDoubleValue("a2", is, il, ip);
                double a3 = time.getDoubleValue("a3", is, il, ip);
                double a4 = time.getDoubleValue("a4", is, il, ip);
                
                float ftdc_corr = t+FTOFFSET;
                
                if (tdcs.hasItem(is,il,ip)) {
                    float radc = (float)Math.sqrt(adc);
                    for (float tdcc : tdcs.getItem(is,il,ip)) {
                    	float tcor = tdcc -  getTriggerPhase() - (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - FTOFFSET;
//                        float tdcmc = tdcm - a0 -  (float)gtw.getDoubleValue("time_walk",is,il,0)/radc - a2/radc - a3 - a4/Math.sqrt(radc);
                        float tdif = tcor - t; 
                        if (Math.abs(tdif)<TMFCUT && tdif<tmax && bcut2) {tmax = tdif; tdc = tcor;} 
                        dgm.fill("at1"+tag,sadc,tcor); 
                    }
                } 
                
                if(sadc>0) {
                	dgm.fill("at0"+tag,sadc,t);
                	dgm.fill("at3"+tag,sadc,tmax);
                	dgm.fill("at2"+tag,sadc,tdc);
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
			
			int[] crate = {1,3,7,9,13,15,19,21,25,27,31,33};
			int[]  slot = {3,4,5,6,7,8,9,10,13,14,15,16,17,18};

	        for (int cr : crate) {
	        for (int sl : slot) {
	        for (int ch=0; ch<16; ch++) {							
	        	switch (table) {
			    case "ped": line = getPED(cr,sl,ch); break;
	        	}
	        	if(!skip(cr,sl,ch)) {
	        		System.out.println(line);
	        		outputBw.write(line);
	        		outputBw.newLine();
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
}
