package org.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.clas.tools.DataProvider;
import org.clas.tools.FitData;
import org.clas.tools.TOFPaddle;
import org.clas.tools.TimeLine;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.tree.TreeFile;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

import Jampack.Inv;

public class ECpi0 extends DetectorMonitor{
	
	H1F                              h1 = null;      
	H2F                              h2 = null;
	DataGroup                        dg = null;
	
    String[]                        det = {"pcal","ecin","ecou"};
    String[]                      dname = {"PCAL ","ECIN ","ECOU "};
    String[]                          v = new String[]{"u","v","w"};    
    int[]                          npmt = {68,62,62,36,36,36,36,36,36};    
    int[]                         iidet = {1,4,7};
	
	ECPart                         part = new ECPart();
    List<TOFPaddle>          paddleList = null;
    List<List<DetectorResponse>>    res = new ArrayList<List<DetectorResponse>>();   
    List<DetectorResponse>  ecClusters  = null;
    Map<String,Integer>            smap = new HashMap<String,Integer>();  
    
    int                        runIndex = 0;
       
    TreeFile                       tree = null; 
    float[]                        rowf = null; 
    int                         mcevent = 0;
    
    float                          mpi0 = 134.976f;
    float                          pmin = 100f;
    float                          pmax = 170f;
    float                          fmin = 30f;
    float                          fmax = 240f;
    
    public ECpi0(String name) {
        super(name);
        this.setDetectorTabNames("PI0",
        		                 "UVW",
                                 "SIJ",
        		                 "OPAE",
                                 "OPAX",                                 
                                 "IMvOPA",
                                 "IMvE1E2",
                                 "IMvEPI0",
                                 "TIJ", 
                                 "EPI0vTh",
                                 "XY",
                                 "FTOF",
                                 "MCPHOT",
                                 "Fits",
                                 "Summary",
                                 "Timeline");
        
        this.usePCCheckBox(true);
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
    }
    
    public void localinit() {
    	System.out.println("ECpi0.localinit()");
        configEngine("pi0");
        part.setGeom("2.5");  
        part.setConfig("pi0");  
        part.setGoodPhotons(1212);   
    	tl.setFitData(Fits);
//        openTree();
    }
    
    public void localclear() {
    	System.out.println("ECpi0.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	tl.Timeline.clear();
    	slider.setValue(0);
    }  
    
    public void openTree() {
    	tree  = new TreeFile(jawPath+"mctrue.hipo","T","ev:avgT:avgZ:vz:mvz:aR:vR:mvR:tedep");
    	rowf = new float[9];    	
    }
    
    @Override
    public void createHistos(int run) { 
	    System.out.println("ECpi0:createHistos("+run+")");
        setRunNumber(run);
        runlist.add(run);
        createUVWHistos(0, 1, 30,50.,250.," Inv. Mass (MeV)"); // invm vs. photon 1 strips
    	createUVWHistos(0, 2, 30,50.,250.," Inv. Mass (MeV)"); // invm vs. photon 2 strips
    	create1DHistos(1,75,5.,400.,"uvw","Two Photon Inv. Mass (MeV)"); // Sector ij photons i==j
    	createSIJHistos(2,130,5.,700.,"sij","Two Photon Inv. Mass (MeV)"); // Sector ij photons i!=j
    	create2DHistos(3,50,0.,20.,50,0.,5.0,"opae","Two Photon Opening Angle (deg)","E1*E2 (GeV^2)");
    	addFunctions(3,"IM1","0.13495*0.13495/2/(1-cos(x*3.14159/180.))",3.65,20.,1,2);
    	addFunctions(3,"IM2","0.12495*0.12495/2/(1-cos(x*3.14159/180.))",3.4,20.,5,2);
    	addFunctions(3,"IM3","0.14495*0.14495/2/(1-cos(x*3.14159/180.))",4.0,20.,5,2);
    	create2DHistos(4,40,0.,20.,50,-1.,1.,      "opax",   "Two Photon Opening Angle (deg)","X:(E1-E2)/(E1+E2)");
    	create2DHistos(5,60,-150.,150.,60,1.,20.,  "imopa",  "Inv. Mass Error (MeV)",         "Two Photon Opening Angle (deg)");
    	create2DHistos(6,60,-150.,150.,60,0.,8.,   "ime1e2", "Inv. Mass Error (MeV)",         "E1*E2 (GeV^2)");
    	create2DHistos(7,60,-150.,150.,60,0.,10.,  "imepi",  "Inv. Mass Error (MeV)",         "Pizero Energy (GeV)");
    	create2DHistos(8,60,-150.,150.,60,-15.,15.,"tij",    "Inv. Mass Error (MeV)",         "Time Difference (Phot1-Phot2) (ns)");
    	create2DHistos(9,60,3.,32.,60,0.,10.,      "pite",   "Pizero Theta (deg)",            "Pizero Energy (GeV)");
    	createXYHistos(10,1,60,410);
    	createXYHistos(10,2,60,410);
    	create1DHistos(11,100,0.,50.,"ftof","Energy (MeV)");
    	createMCHistos(12);
    }
    
    @Override       
    public void plotHistos(int run) {
   	   plotSummary(run);
   	   plotAnalysis(run);
    } 
          
    public void plotSummary(int run) {
   	    if(dropSummary) return;
        setRunNumber(run);
        plotUVW(0);
        if(isAnalyzeDone) {updateUVW(1);}else{plotPI0Summary(1);}
        plotPI0Summary(2);    	
        plotPI0Summary(3);    	
        plotPI0Summary(4);   
        plotPI0Summary(5);
        plotPI0Summary(6);
        plotPI0Summary(7);
        plotPI0Summary(8);
        plotPI0Summary(9);
        plotXYSummary(10);
        plotPI0Summary(11);
//        plotMCPHOT(12);
    }
    
    public void plotAnalysis(int run) {
   	   setRunNumber(run);
   	   if(!isAnalyzeDone) return;
   	   if(!dropSummary) {updateFits(13);plotMeanHWSummary(14);} 
       updateUVW(1); plotTimeLines(15);
    }  
    
    public void createPHOTHistos(int k) {
    	
 	   int run = getRunNumber();

       int is=0, n=0;
       String tag = is+"_"+n+"_"+k+"_"+run;
       h1 = new H1F("pi0_pcal_u_"+tag,"pi0_pcal_u_"+tag, 100, 0., 2.);
       dg = new DataGroup(1,1);
       dg.addDataSet(h1, 0);
       this.getDataGroup().add(dg,is,n,k,run);        	
    }
    
    public void createUVWHistos(int k, int n, int nch, double x1, double x2, String txt) {
    	
	   int run = getRunNumber();
   
       for (int is=1; is<7; is++) {
           String tag = is+"_"+n+"_"+k+"_"+run;
           dg = new DataGroup(3,3);
           h2 = new H2F("pi0_pcal_u_"+tag,"pi0_pcal_u_"+tag, nch, x1, x2, 68, 1., 69.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("PCAL U Strips"); 
           dg.addDataSet(h2,0);  
           h2 = new H2F("pi0_pcal_v_"+tag,"pi0_pcal_v_"+tag, nch, x1, x2, 62, 1., 63.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("PCAL V Strips");        
           dg.addDataSet(h2,1);            
           h2 = new H2F("pi0_pcal_w_"+tag,"pi0_pcal_w_"+tag, nch, x1, x2, 62, 1., 63.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("PCAL W Strips");  
           dg.addDataSet(h2,2); 
       
           h2 = new H2F("pi0_ecin_u_"+tag,"pi0_ecin_u_"+tag, nch, x1, x2, 36, 1., 37.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("ECIN U Strips");    
           dg.addDataSet(h2,3);  
           h2 = new H2F("pi0_ecin_v_"+tag,"pi0_ecin_v_"+tag, nch, x1, x2, 36, 1., 37.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("ECIN V Strips");        
           dg.addDataSet(h2,4);            
           h2 = new H2F("pi0_ecin_w_"+tag,"pi0_ecin_w_"+tag, nch, x1, x2, 36, 1., 37.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("ECIN W Strips");  
           dg.addDataSet(h2,5); 
       
           h2 = new H2F("pi0_ecou_u_"+tag,"pi0_ecou_u_"+tag, nch, x1, x2, 36, 1., 37.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("ECOU U Strips");    
           dg.addDataSet(h2,6);  
           h2 = new H2F("pi0_ecou_v_"+tag,"pi0_ecou_v_"+tag, nch, x1, x2, 36, 1., 37.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("ECOU V Strips");        
           dg.addDataSet(h2,7);            
           h2 = new H2F("pi0_ecou_w_"+tag,"pi0_ecou_w_"+tag, nch, x1, x2, 36, 1., 37.);
           h2.setTitleX("Sector "+is+txt); h2.setTitleY("ECOU W Strips");  
           dg.addDataSet(h2,8);   
           this.getDataGroup().add(dg,is,n,k,run);
       }            
   }
    
    public void create1DHistos(int k, int nch, double x1, double x2, String var, String txt) {
    	
	    int run = getRunNumber();
        dg = new DataGroup(3,2);
        
        for (int is=1; is<7; is++) {
            String tag = var+"_"+is+"_"+k+"_"+run;
            h1 = new H1F("pi0_"+tag,"pi0_"+tag, nch, x1, x2);
            h1.setTitleX("Sector "+is+" "+txt);  
            dg.addDataSet(h1,is-1);   
        }
        this.getDataGroup().add(dg,0,0,k,run);
    }
    
    public void create2DHistos(int k, int nchx, double x1, double x2, int nchy, double y1, double y2, String var, String txtx, String txty) {
    	
	    int run = getRunNumber();
        DataGroup dg = new DataGroup(3,2);
        
        for (int is=1; is<7; is++) {
            String tag = var+"_"+is+"_"+k+"_"+run;
            h2 = new H2F("pi0_"+tag,"pi0_"+tag, nchx, x1, x2, nchy, y1, y2);
            h2.setTitleX("Sector "+is+" "+txtx); h2.setTitleY(txty); 
            dg.addDataSet(h2,is-1);   
        }
        this.getDataGroup().add(dg,0,0,k,run);    	
    }
    
    public void createXYHistos(int k, int n, int nb, int bmx) {
   	 
	    int run = getRunNumber();
        
        String[] t = {"e","w","r"};
        
        for (int i=0; i<3; i++) {
           DataGroup dg = new DataGroup(3,2);
           for (int d=0; d<3; d++) {
              h2 = new H2F("pi0_"+det[d]+"_xy_"+t[i]+"_"+n+"_"+k+"_"+run,"hi_"+det[d]+"_xy_"+t[i]+"_"+k+"_"+run,nb,-bmx,bmx,nb,-bmx,bmx);
              h2.setTitleX("Photon "+n+" X(cm)"); h2.setTitleY("Photon "+n+" Y(cm)");
              dg.addDataSet(h2,d);  
	       }
           this.getDataGroup().add(dg,i,n,k,run);
        }
    }    
    
    public void addFunctions(int k, String fnam, String f, double x1, double x2, int lcol, int lwid) {
    	
	    int run = getRunNumber();
	    String tag = "_"+run;
	    this.getDataGroup().getItem(0,0,k,run);
        F1D f1 = new F1D(fnam+tag,f,x1,x2); f1.setLineColor(lcol); f1.setLineWidth(lwid);
    	for (int is=1; is<7; is++) {
    		this.getDataGroup().getItem(0,0,k,run).addDataSet(f1, is-1);
    	}
    	
    }
       
    public void createSIJHistos(int k, int nch, double x1, double x2, String var, String txt) {
    	
	    int run = getRunNumber();
	    
        int n = 0;
        DataGroup dg = new DataGroup(3,5);
        
        for (int is=1; is<7; is++) {   
            for (int j=is+1; j<7; j++) {
                String tag = var+"_"+is+"_"+j+"_"+k+"_"+run;
                h1 = new H1F("pi0_"+tag,"pi0_"+tag, nch, x1, x2);
                h1.setTitleX("Sector Pair "+is+j+" "+txt);  
                dg.addDataSet(h1,n);  n++; 
                smap.put(is+"_"+j,n);
            }
        }
        this.getDataGroup().add(dg,0,0,k,run);    	
    }
    
    public void createMCHistos(int k) {
 	   int run = getRunNumber();
       double xmin=-20, xmax=20.;
       double ymin= 20.,ymax=30.;
       String tag = k+"_"+run;
       for (int idet=0; idet<3; idet++) {
           DataGroup dg = new DataGroup(4,2);
    	   if (idet==0) {
              h2 = new H2F("pi0_mc80_"+det[idet]+"_"+tag, 50, 2*xmin,2*xmax, 3, 1., 4.);  
              dg.addDataSet(h2, 0);
              h2 = new H2F("pi0_mc81_"+det[idet]+"_"+tag, 50, xmin,xmax,100, ymin,ymax);  
              dg.addDataSet(h2, 1);
              h2 = new H2F("pi0_mc82_"+det[idet]+"_"+tag, 50, xmin,xmax,100, ymin,ymax);  
              dg.addDataSet(h2, 2);
              h2 = new H2F("pi0_mc83_"+det[idet]+"_"+tag, 50, 2*xmin,2*xmax,100, ymin,ymax);  
              dg.addDataSet(h2, 3);
              h2 = new H2F("pi0_mc84_"+det[idet]+"_"+tag, 50, xmin,xmax,100, ymin,ymax);  
              dg.addDataSet(h2, 4);
    	   }
           h2 = new H2F("pi0_mc85_"+det[idet]+"_"+tag, 50, -1.5, 1.5,3, 1., 4.);  
           dg.addDataSet(h2, 5);
           h2 = new H2F("pi0_mc86_"+det[idet]+"_"+tag, 50, -1, 1, 100, ymin,ymax);  
           dg.addDataSet(h2, 6);
           h2 = new H2F("pi0_mc87_"+det[idet]+"_"+tag, 50, 0., 3.5, 40, 0.15, 0.35);  
           dg.addDataSet(h2, 7);
           h2 = new H2F("pi0_mc88_"+det[idet]+"_"+tag, 50, 0.98, 1.004, 50, 699.,720.);  
           dg.addDataSet(h2, 8);
           h2 = new H2F("pi0_mc89_"+det[idet]+"_"+tag, 50, 0.96, 1.01, 50, 699.,704.);  
           dg.addDataSet(h2, 9);
           h1 = new H1F("pi0_mc810_"+det[idet]+"_"+tag, 200,580.,680.);  
           dg.addDataSet(h1, 10);
           h1 = new H1F("pi0_mc811_"+det[idet]+"_"+tag, 200,580.,680.);  
           dg.addDataSet(h1, 11);
           h2 = new H2F("pi0_mc812_"+det[idet]+"_"+tag, 50,700,730,50,700,703);  
           dg.addDataSet(h2, 12);
           h2 = new H2F("pi0_mc813_"+det[idet]+"_"+tag, 50,23,24,60,23,24.5);  
           dg.addDataSet(h2, 13);
           h2 = new H2F("pi0_mc814_"+det[idet]+"_"+tag, 50,695,703,60,23,24.5);  
           dg.addDataSet(h2, 14);
           h2 = new H2F("pi0_mc815_"+det[idet]+"_"+tag, 50, 0, 1, 50, 699.,720.);  
           dg.addDataSet(h2, 15);
           this.getDataGroup().add(dg,0,idet,k,run);
       }
           
    }
    
    public void processMC(DataEvent event, DataBank ecBank) {
    	
        int run = getRunNumber();
        
        double refE=0,refP=0,refTH=25,tmax=30;
        double pcx=0,pcy=0,pcz=0;
    	
        if(event.hasBank("MC::Particle")==true) {
            DataBank bank = event.getBank("MC::Particle");
            int   pid = bank.getInt("pid", 0);
            float ppx = bank.getFloat("px",0);
            float ppy = bank.getFloat("py",0);
            float ppz = bank.getFloat("pz",0);
            double  rm = (pid==111) ? 0.1349764:0.;
            refP  = Math.sqrt(ppx*ppx+ppy*ppy+ppz*ppz);  
            refE  = Math.sqrt(refP*refP+rm*rm);    
            refTH = Math.acos(ppz/refP)*180/Math.PI;
        }
            
        if(event.hasBank("MC::True")==true) {
            DataBank bank = event.getBank("MC::True");
            rowf[0] = mcevent++;
            for(int i=0; i < bank.rows(); i++) {
                float pcT = bank.getFloat("avgT",i); rowf[1]=pcT;
                float pcX = bank.getFloat("avgX",i);
                float pcY = bank.getFloat("avgY",i);
                float pcZ = bank.getFloat("avgZ",i); rowf[2]=pcZ;
                float  vx = bank.getFloat("vx",i);  
                float  vy = bank.getFloat("vy",i);  
                float  vz = bank.getFloat("vz",i);   rowf[3]=vz; 
                float mvx = bank.getFloat("mvx",i);  
                float mvy = bank.getFloat("mvy",i);  
                float mvz = bank.getFloat("mvz",i);  rowf[4]=mvz; 
                float totE= bank.getFloat("totEdep",i); rowf[8]=totE;
                rowf[5]=(float) Math.sqrt(pcX*pcX+pcY*pcY+pcZ*pcZ);
                rowf[6]=(float) Math.sqrt(vx*vx+vy*vy+vz*vz);
                rowf[7]=(float) Math.sqrt(mvx*mvx+mvy*mvy+mvz*mvz);
                tree.addRow(rowf);
                if(pcT<tmax){pcx=pcX; pcy=pcY; pcz=pcZ ; tmax = pcT;}
            }
        } 
        
        for(int i=0; i< ecBank.rows(); i++) {
            int sector = ecBank.getByte("sector", i);
            int  layer = ecBank.getByte("layer",  i);
            double x = ecBank.getFloat("x", i);
            double y = ecBank.getFloat("y", i);
            double z = ecBank.getFloat("z", i);        	
        	if (layer==1&&sector==2) {
                ((H1F) this.getDataGroup().getItem(0,0,12,run).getData(10).get(0)).fill(z);        		
                ((H1F) this.getDataGroup().getItem(0,0,12,run).getData(11).get(0)).fill(0.1*pcz);        		
        	}
        }

        double    mcR = Math.sqrt(pcx*pcx+pcy*pcy+pcz*pcz);
        double    mcB = tmax<30 ? 0.1*mcR/tmax/29.98:0;
        double mcThet = Math.asin(Math.sqrt(pcx*pcx+pcy*pcy)/mcR)*180/Math.PI;
        
        ((H2F) this.getDataGroup().getItem(0,0,12,run).getData(5).get(0)).fill(refTH-mcThet,1.);
        
        res.clear();
        
        for (int idet=0; idet<3; idet++) {
            res.add(part.eb.getUnmatchedResponses(ecClusters, DetectorType.ECAL,iidet[idet]));
            for(int i = 0; i < res.get(idet).size(); i++){
                int        is = res.get(idet).get(i).getDescriptor().getSector();
                double energy = res.get(idet).get(i).getEnergy();
                double      X = res.get(idet).get(i).getPosition().x();
                double      Y = res.get(idet).get(i).getPosition().y();
                double      Z = res.get(idet).get(i).getPosition().z();
                double      T = res.get(idet).get(i).getTime();
                double    pcR = Math.sqrt(X*X+Y*Y+Z*Z);
                double pcThet = Math.asin(Math.sqrt(X*X+Y*Y)/pcR)*180/Math.PI;
                double   beta = T>0 ? pcR/T/29.98:0.;
                if (false&&idet==0) {
                    System.out.println("GEMC: "+pcx+" "+pcy+" "+pcz+" idet="+idet);
                    System.out.println("Cluster: "+X+" "+Y+" "+Z+" idet="+idet);
                    System.out.println("GEMC R: "+0.1*mcR);
                    System.out.println("Cluster R: "+pcR);
                    System.out.println("GEMC thet: "+mcThet);
                    System.out.println("Cluster thet: "+pcThet);
                    System.out.println("GEMC T: "+tmax);
                    System.out.println("Cluster T: "+T);
                    System.out.println("GEMC beta: "+mcB);
                    System.out.println("Cluster beta: "+beta);
                    System.out.println(" ");
                }
                if(idet==0) {
                    ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(0).get(0)).fill(0.1*pcx-X,1.);
                    ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(0).get(0)).fill(0.1*pcy-Y,2.);
                    ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(0).get(0)).fill(0.1*pcz-Z,3.);

                    ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(1).get(0)).fill(0.1*pcx-X,refTH);
                    ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(2).get(0)).fill(0.1*pcy-Y,refTH);
                    ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(3).get(0)).fill(0.1*pcz-Z,refTH);
                    if((0.1*pcz-Z)>0) {
                        ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(8).get(0)).fill(mcB,0.1*mcR);
                    	((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(12).get(0)).fill(0.1*mcR,pcR);
                        ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(13).get(0)).fill(tmax,T);
                        ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(4).get(0)).fill(0.1*mcR-pcR,refTH);
                        ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(14).get(0)).fill(pcR,T);
                        ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(15).get(0)).fill(energy,0.1*mcR);
                    }
                }
                ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(5).get(0)).fill(refTH-pcThet,2.);
                ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(6).get(0)).fill(refTH-pcThet,refTH); //pcThet-refTH
                ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(7).get(0)).fill(refE*1e-3,energy/refE,1.); // Layer Cluster Normalized Energy
                ((H2F) this.getDataGroup().getItem(0,idet,12,run).getData(9).get(0)).fill(beta,pcR);
            }
        }
    
    
    }

    @Override
    public void processEvent(DataEvent event) {  
    	
        int run = getRunNumber();
        
        if (dropBanks) dropBanks(event);

        ecClusters = part.readEC(event);  
        
        if (ecClusters.size()==0) return;
        
        DataBank ecBank = event.getBank("ECAL::clusters");
               
        if(event.hasBank("MC::Particle")==true) {processMC(event, ecBank); return;}
        
        double[]  esum = {0,0,0,0,0,0};
        int[][]  nesum = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};    
  	    double pcalE[] = new double[6];
                
        res.clear();
        
        for (int idet=0; idet<3; idet++) {
            res.add(part.eb.getUnmatchedResponses(ecClusters, DetectorType.ECAL,iidet[idet]));
            for(int i = 0; i < res.get(idet).size(); i++){
                int        is = res.get(idet).get(i).getDescriptor().getSector();
                double energy = res.get(idet).get(i).getEnergy();
                double      X = res.get(idet).get(i).getPosition().x();
                double      Y = res.get(idet).get(i).getPosition().y();
                double      Z = res.get(idet).get(i).getPosition().z();
                double    pcR = Math.sqrt(X*X+Y*Y+Z*Z);
                double pcThet = Math.asin(Math.sqrt(X*X+Y*Y)/pcR)*180/Math.PI;
//                ecPix[idet].strips.hmap2.get("H2_a_Hist").get(is,4,0).fill(energy*1e3,4,1.);          // Layer Cluster Energy
                if(idet==0) pcalE[is-1] += energy*1e3;
                if(energy*1e3>10) {esum[is-1]+=energy*1e3; nesum[idet][is-1]++;}
            }
        }
        
        if(event.hasBank("MIP::event")){          
            DataBank bank = event.getBank("MIP::event");
            for(int i=0; i < bank.rows(); i++) part.mip[i]=bank.getByte("mip", i);
            
        } else if (event.hasBank("FTOF::adc")) {
            paddleList = DataProvider.getPaddleList(event);          
            double[] thresh = {500,1000,1000}; 
            for (int i=0; i<6; i++) part.mip[i]=0;       
            if (paddleList!=null) {
                for (TOFPaddle paddle : paddleList){           
                    int toflay = paddle.getDescriptor().getLayer();            
                    int   isec = paddle.getDescriptor().getSector();   
                    double gmean = paddle.geometricMean();
                    if(toflay==2)((H1F) this.getDataGroup().getItem(0,0,11,run).getData(isec-1).get(0)).fill(gmean); 
                    part.mip[isec-1] = (gmean>thresh[toflay-1]) ? 1:0;
                }
            }
            
        } else if (event.hasBank("REC::Scintillator")) {
            double[] thresh = {7,8,8}; 
            for (int i=0; i<6; i++) part.mip[i]=0;       
            DataBank bank = event.getBank("REC::Scintillator");
            for (int i=0; i<bank.rows(); i++) {
      		   if (bank.getByte("detector", i)==12) {
                  int   toflay = bank.getByte("layer", i);
                  int     isec = bank.getByte("sector", i);
                  float energy = bank.getFloat("energy",i);
                  if(toflay==2)((H1F) this.getDataGroup().getItem(0,0,11,run).getData(isec-1).get(0)).fill(energy); 
                  part.mip[isec-1] = (toflay<3&&energy>thresh[toflay-1]) ? 1:0;
               }
            }      
        } 
        
        part.getNeutralResponses(ecClusters);
        
        for (int is=1; is<7; is++) {
           
            if (part.mip[is-1]!=1) {  // No FTOF MIP in sector
                double invmass = Math.sqrt(part.getTwoPhotonInvMass(is));
                double    inv3 = invmass*1e3;
                double    invd = (invmass-part.mpi0)*1e3;
                double     opa = Math.acos(part.cth)*180/3.14159;
                
                boolean    ivmcut = inv3>pmin && inv3<pmax;
                boolean badPizero = part.X>1 && opa<0;
        
                if(invmass>0&&part.iis[0]>0&&part.iis[1]>0&&!badPizero) {                                                    
                    if(part.iis[0]< part.iis[1]) ((H1F) this.getDataGroup().getItem(0,0,2,run).getData(smap.get(part.iis[0]+"_"+part.iis[1])-1).get(0)).fill(invmass*1e3);   
                    
                    if(part.iis[0]==part.iis[1]) { //Both photons in same sector
                    	
                    	((H1F) this.getDataGroup().getItem(0,0,1,run).getData(part.iis[0]-1).get(0)).fill(inv3);   
                    	((H2F) this.getDataGroup().getItem(0,0,5,run).getData(part.iis[0]-1).get(0)).fill(invd,opa);                    
                    	((H2F) this.getDataGroup().getItem(0,0,6,run).getData(part.iis[0]-1).get(0)).fill(invd,part.e1c*part.e2c);                    
                    	((H2F) this.getDataGroup().getItem(0,0,7,run).getData(part.iis[0]-1).get(0)).fill(invd,Math.sqrt(part.tpi2));                    
                    	((H2F) this.getDataGroup().getItem(0,0,8,run).getData(part.iis[0]-1).get(0)).fill(invd,part.t[0][0]-part.t[1][0]);  
                    	
                        if(nesum[0][is-1]>1 && nesum[1][is-1]>0) {
                    	   ((H2F) this.getDataGroup().getItem(0,0,3,run).getData(part.iis[0]-1).get(0)).fill(opa,part.e1c*part.e2c);
                       	   ((H2F) this.getDataGroup().getItem(0,0,4,run).getData(part.iis[0]-1).get(0)).fill(opa,part.X);
                        }                                   
                        if (ivmcut) {
                        	((H2F) this.getDataGroup().getItem(0,0,9,run).getData(part.iis[0]-1).get(0)).fill(Math.acos(part.cpi0)*180/Math.PI,Math.sqrt(part.tpi2));                    
                        }  
                        for (int id=0; id<3; id++) {
                            for (int im=0; im<2; im++) {
                                if(ecBank!=null) {
                                    float ipU = (ecBank.getInt("coordU", part.iip[im][id])-4)/8+1;
                                    float ipV = (ecBank.getInt("coordV", part.iip[im][id])-4)/8+1;
                                    float ipW = (ecBank.getInt("coordW", part.iip[im][id])-4)/8+1;
                                    ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(id*3+0).get(0)).fill(inv3,ipU);
                                    ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(id*3+1).get(0)).fill(inv3,ipV);
                                    ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(id*3+2).get(0)).fill(inv3,ipW);    
                                }
                                if (ivmcut) {
                                	((H2F) this.getDataGroup().getItem(id,im+1,10,run).getData(0).get(0)).fill(-part.x[im][id], part.y[im][id],1.);
                                    ((H2F) this.getDataGroup().getItem(id,im+1,10,run).getData(1).get(0)).fill(-part.x[im][id], part.y[im][id],invmass/part.mpi0);
                                }
                            }
                       }    
                   }
                }
            }        
        }        
    }
    
    private void updateUVW(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        GraphErrors g;
        
        c.clear();
        c.divide(3, 2);
        
        for (int is=1; is<7; is++) {
        	g=tl.fitData.getItem(is,0,0,getRunNumber()).getGraph();
            c.cd(is-1); c.getPad(is-1).getAxisY().setRange(0.,g.getMax()*1.1);
            c.draw(tl.fitData.getItem(is,0,0,getRunNumber()).getHist());
            c.draw(g,"same");
            DataLine line = new DataLine(mpi0,0.,mpi0,g.getMax()*1.1); line.setLineColor(3); c.draw(line); 
        }      
        
    } 
    
    private void updateFits(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int    pc = getActivePC();
        int    is = getActiveSector(); 
        int     i = getActiveLayer();
        int     j = getActiveView();
        
        int    np = npmt[i*3+j];        
        
        c.clear();
        if (np==36) c.divide(6,6); if (np>36) c.divide(8,9);
        
        for (int ip=0; ip<np ; ip++) {
            c.cd(ip); c.getPad(ip).getAxisY().setLog(false);
            c.draw(tl.fitData.getItem(is,i+10*(pc+1)*(pc+1)*(j+1),ip+1,getRunNumber()).getHist());
            c.draw(tl.fitData.getItem(is,i+10*(pc+1)*(pc+1)*(j+1),ip+1,getRunNumber()).getGraph(),"same");
       }
    }  
    
    @Override
    public void plotEvent(DataEvent de) {
    	analyze();
    }

    public void analyze() {    
        System.out.println("ECpi0.analyze()");
        fitGraphs(1,7,0,(dropSummary)?0:3,0,(dropSummary)?0:3);
        if(!isAnalyzeDone) createTimeLineHistos();
        fillTimeLineHisto();
        System.out.println("Finished");
        isAnalyzeDone = true;       
    }
    
    public void fitGraphs(int is1, int is2, int id1, int id2, int il1, int il2) {
        int run = getRunNumber();
        for (int is=is1; is<is2; is++) {
           tl.fitData.add(fitEngine((H1F)this.getDataGroup().getItem(0,0,1,run).getData(is-1).get(0),3,pmin,pmax,fmin,fmax),is,0,0,run); 
           for (int id=id1; id<id2; id++) {
               for (int il=il1; il<il2; il++) {
               	  for (int pc=0; pc<2; pc++) {
                     H2F h = (H2F) this.getDataGroup().getItem(is,pc+1,0,run).getData(id*3+il).get(0);
           	         for (int i=0; i<npmt[id*3+il]; i++) tl.fitData.add(fitEngine(h.sliceY(i),3,pmin,pmax,55,fmax),is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); //PMT slices
       		         fitStore(is, id, il, pc, run, mpi0);
               	  }
               } 
           }
        }
    }
    
    public void fitStore(int is, int id, int il, int pc, int run, float nrm) {
    	int np = npmt[id*3+il];
        double[]      x = new double[np]; double[]  ymean = new double[np]; double[] yrms = new double[np];
        double[]     xe = new double[np]; double[] ymeane = new double[np]; double[]   ye = new double[np]; 
        double[]  yMean = new double[np]; 
        for (int i=0; i<np; i++) {
            x[i] = i+1; xe[i]=0; ye[i]=0; yrms[i]=0; 
            FitData fd = tl.fitData.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),i+1,run); 
            fd.graph.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            fd.hist.getAttributes().setTitleX("Sector "+is+" "+det[id]+" "+v[il]+(i+1));
            double mean = fd.mean;                        
            if(mean>0) yrms[i] = fd.sigma/mean; 
                      yMean[i] = fd.getMean()/nrm;
                      ymean[i] = mean/nrm;
                     ymeane[i] = fd.meane/nrm;
        }
        GraphErrors mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,ymean,xe,ymeane);                   
        GraphErrors Mean = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yMean,xe,ymeane);                   
        GraphErrors  rms = new GraphErrors("MIP_"+is+"_"+id+" "+il,x,yrms,xe,ye);                  
        FitSummary.add(mean, is,id+10*(pc+1)*(pc+1)*(il+1),1,run);
        FitSummary.add(rms,  is,id+10*(pc+1)*(pc+1)*(il+1),2,run);                    
        FitSummary.add(Mean, is,id+10*(pc+1)*(pc+1)*(il+1),5,run);        	        
    }   
    
    public void plotMeanHWSummary(int index) {
        
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int           pc = getActivePC();
        int            n = 0;
        
        List<DataLine> lines = new ArrayList<DataLine>();
        
        Boolean t = TLname!="UVW";
        float ymin=0.8f, ymax=1.2f;
        
        c.clear(); c.divide(3, 6);
        
        for (int is=1; is<7; is++) {
        for (int id=0; id<3; id++) {
        	GraphErrors hwplot1 = new GraphErrors();
        	GraphErrors hwplot2 = new GraphErrors();
        	int m=0; lines.clear();
            for (int il=0; il<3; il++) {           	
                F1D f1 = new F1D("p0","[a]",0.,npmt[id*3+il]); f1.setParameter(0,1);
                GraphErrors plot1 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),1,getRunNumber());
                GraphErrors plot2 = FitSummary.getItem(is,id+10*(pc+1)*(pc+1)*(il+1),5,getRunNumber());
                for (int ip=0; ip<npmt[id*3+il]; ip++) {m++;
        	        hwplot1.addPoint(m, plot1.getDataY(ip), plot1.getDataEX(ip), plot1.getDataEY(ip));
        	        hwplot2.addPoint(m, plot2.getDataY(ip), plot2.getDataEX(ip), plot2.getDataEY(ip));
        	        if(Math.floorMod(t?m:ip, t?TimeSlice.get(TLname):npmt[id*3+il])==(t?1:0)) {
        	        	DataLine line = new DataLine(m,ymin,m,ymax) ; line.setLineColor(1); line.setLineWidth(1); 
        	        	lines.add(line);
        	        }
                }
            }
            hwplot2.setMarkerColor(1);
            c.cd(n); c.getPad(n).getAxisY().setRange(ymin,ymax); c.getPad(n).setAxisTitleFontSize(14); c.getPad(n).setTitleFontSize(16);
            if(n==0||n==3||n==6||n==9||n==12||n==15) hwplot1.getAttributes().setTitleY("MEAN / MASS");
            hwplot1.getAttributes().setTitleX("SECTOR "+is+" "+det[id]+" PMT");
            F1D f1 = new F1D("p0","[a]",0.,m); f1.setParameter(0,1); f1.setLineColor(3); f1.setLineWidth(2);
            c.draw(hwplot1); /*c.draw(hwplot2,"same");*/ for(DataLine line: lines) c.draw(line); c.draw(f1,"same"); n++; 
        }
        }        
    }    
    
    public void plotXYSummary(int index) {        
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int   run = getRunNumber();
        
        c.clear(); c.divide(3,2);
        
        for (int idet=0; idet<3; idet++) {
           H2F h1 = (H2F) this.getDataGroup().getItem(idet,getActivePC()==0?2:1,index,run).getData(0).get(0); 
           H2F h2 = (H2F) this.getDataGroup().getItem(idet,getActivePC()==0?2:1,index,run).getData(1).get(0);  
           H2F h3 = (H2F) this.getDataGroup().getItem(idet,getActivePC()==0?2:1,index,run).getData(2).get(0);  
            
           h3 = h2.divide(h2, h1); h3.setTitle("MEAN/MASS "+dname[idet]);
           c.cd(idet);   c.getPad(idet).getAxisZ().setLog(true); c.draw(h1);
           c.cd(idet+3); c.getPad(idet+3).getAxisZ().setLog(false); c.getPad(idet+3).getAxisZ().setRange(0.7, 1.3); c.draw(h3);
        }	
   
    }   
    
    public void plotMCPHOT(int index) {
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int   run = getRunNumber();
        
        c.clear(); c.divide(3,5);
        
        H1F h1 ; H2F h2;
	    int ii=0;
	  
        h1 = ((H2F) this.getDataGroup().getItem(0,0,12,run).getData(0).get(0)).sliceY(0) ;  
        h1.setTitleX("PCAL X: GEMC-Cluster (cm)"); h1.setFillColor(2); 
        h1.setOptStat(Integer.parseInt("1001100")); h1.setTitle(" "); c.cd(ii); ii++; c.draw(h1);

        h1 = ((H2F) this.getDataGroup().getItem(0,0,12,run).getData(0).get(0)).sliceY(1) ;  
        h1.setTitleX("PCAL Y: GEMC-Cluster (cm)"); h1.setFillColor(2);
        h1.setOptStat(Integer.parseInt("1001100")); h1.setTitle(" "); c.cd(ii); ii++; c.draw(h1);
        
        h1 = ((H2F) this.getDataGroup().getItem(0,0,12,run).getData(0).get(0)).sliceY(2) ;  
        h1.setTitleX("PCAL Z: GEMC-Cluster (cm)"); h1.setFillColor(2);
        h1.setOptStat(Integer.parseInt("11001100")); h1.setTitle(" "); c.cd(ii); ii++; c.draw(h1);
        
        h2 = (H2F) this.getDataGroup().getItem(0,0,12,run).getData(1).get(0) ;  
        h2.setTitleX("PCAL X: GEMC-Cluster (cm)");  h2.setTitleY("True Theta (deg)");
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
        
        h2 = (H2F) this.getDataGroup().getItem(0,0,12,run).getData(2).get(0) ;  
        h2.setTitleX("PCAL Y: GEMC-Cluster (cm)");  h2.setTitleY("True Theta (deg)");
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
        
        h2 = (H2F) this.getDataGroup().getItem(0,0,12,run).getData(3).get(0) ;  
        h2.setTitleX("PCAL Z: GEMC-Cluster (cm)");  h2.setTitleY("True Theta (deg)");
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);

        h1 = ((H2F) this.getDataGroup().getItem(0,0,12,run).getData(5).get(0)).sliceY(0) ;  
        h1.setTitleX(dname[0]+"Theta: True-GEMC (deg)");  h1.setFillColor(2); 
        h1.setOptStat(Integer.parseInt("1001100")); h1.setTitle(" "); c.cd(ii); ii++; c.draw(h1);
        
        h1 = ((H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(5).get(0)).sliceY(1) ;  
        h1.setTitleX(dname[getActiveLayer()]+"Theta: True-Cluster (deg) ");  h1.setFillColor(2); 
        h1.setOptStat(Integer.parseInt("1001100")); h1.setTitle(" "); c.cd(ii); ii++; c.draw(h1);
          
        h2 =  (H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(6).get(0);  
        h2.setTitleX(dname[getActiveLayer()]+"Theta: True-Cluster (deg)");  h2.setTitleY("True Theta (deg)");  
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
 	    	
        h2 =  (H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(8).get(0);  
        h2.setTitleX(dname[0]+"GEMC BETA");  h2.setTitleY("GEMC R (cm)");  
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
	
        h2 =  (H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(9).get(0);  
        h2.setTitleX(dname[getActiveLayer()]+" BETA");  h2.setTitleY("Cluster R (cm)");  
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
/*        
        h1 = (H1F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(10).get(0) ;  
        h1.setTitleX(dname[getActiveLayer()]+"PCAL Z (cm) ");  h1.setFillColor(2); 
        h1.setOptStat(Integer.parseInt("1001100")); h1.setTitle(" "); c.cd(ii); ii++; c.draw(h1);
        
        h1 = (H1F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(11).get(0) ;  
        h1.setTitleX(dname[getActiveLayer()]+"GEMC Z (cm) ");  h1.setFillColor(2); 
        h1.setOptStat(Integer.parseInt("1001100")); h1.setTitle(" "); c.cd(ii); ii++; c.draw(h1);
*/        
        h2 =  (H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(12).get(0);  
        h2.setTitleX(dname[getActiveLayer()]+" GEMC R (cm)");  h2.setTitleY("PCAL R (cm)");  
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
        
        h2 =  (H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(13).get(0);  
        h2.setTitleX(dname[getActiveLayer()]+" GEMC T (ns)");  h2.setTitleY("PCAL T (ns)");  
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
        
        h2 =  (H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(14).get(0);  
        h2.setTitleX(dname[getActiveLayer()]+" PCAL R (cm)");  h2.setTitleY("PCAL T (ns)");  
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
        
        h2 =  (H2F) this.getDataGroup().getItem(0,getActiveLayer(),12,run).getData(15).get(0);  
        h2.setTitleX(dname[getActiveLayer()]+" PCAL E (GeV)");  h2.setTitleY("GEMC R (cm)");  
        h2.setTitle(" "); c.cd(ii); c.getPad(ii).getAxisZ().setLog(true); ii++; c.draw(h2);
    }
    
    public void plotUVW(int index) {    	
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),2,index,getRunNumber()));
    }   
    
    public void plotPI0Summary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));
    }
    
    public void plotMCSummary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,getActiveLayer(),index,getRunNumber()));
    }
    
/*   TIMELINES */
    
    @Override
    public void createTimeLineHistos() {   
    	System.out.println("Initializing "+TLname+" timeline"); 
    	runIndex = 0;   	
    	tl.createTimeLineHisto(10,"Pi0 Mean/Mass","Sector",451,6,1,7);
    	tl.createTimeLineHisto(20,"Pi0 #sigma(E)/E","Sector",451,6,1,7);
    }  
    
    public void fillTimeLineHisto() {    	
		for (int is=1; is<7; is++) {
		   float   y = (float) tl.fitData.getItem(is,0,0,getRunNumber()).mean;
		   float  ye = (float) tl.fitData.getItem(is,0,0,getRunNumber()).meane;			 
   		   float  ys = (float) tl.fitData.getItem(is,0,0,getRunNumber()).sigma;
   		   float yse = (float) tl.fitData.getItem(is,0,0,getRunNumber()).sigmae;			 
           ((H2F)tl.Timeline.getItem(10,0)).fill(runIndex,is,y/mpi0);	
           ((H2F)tl.Timeline.getItem(10,1)).fill(runIndex,is,ye/mpi0);	
           ((H2F)tl.Timeline.getItem(20,0)).fill(runIndex,is,ys/y);	
           ((H2F)tl.Timeline.getItem(20,1)).fill(runIndex,is,(ys/y)*Math.sqrt(Math.pow(yse/ys,2)+Math.pow(ye/y,2)));   		
		}			
		runIndex++;
    }
    
    public void plotTimeLines(int index) {
          
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)); 
        GraphErrors   g2 = null;
        int           is = getActiveSector(); 
        FitData       fd = tl.fitData.getItem(is,0,0,getRunNumber());
       
        c.clear(); c.divide(3, 2); 
    	DataLine line1 = new DataLine(0,is,  runIndex+1,is);                   line1.setLineColor(5);
    	DataLine line2 = new DataLine(0,is+1,runIndex+1,is+1);                 line2.setLineColor(5);
    	DataLine line3 = new DataLine(runIndexSlider,1,  runIndexSlider,7);    line3.setLineColor(5);
    	DataLine line4 = new DataLine(runIndexSlider+1,1,runIndexSlider+1,7);  line4.setLineColor(5);
    	DataLine line5 = new DataLine(-0.5,1,runIndex,1);                      line5.setLineColor(3); line3.setLineWidth(2);

        for (int i=0; i<2; i++) { int i3=i*3; //int nb = tl.getNYbins(i); 
            double min=(i==1)?0.07f:0.96f; double max=(i==1)?0.11f:1.04f; if(doAutoRange){min=min*lMin/250; max=max*lMax/250;}
    		c.cd(i3); c.getPad(i3).setAxisRange(0,runIndex,1,7); c.getPad(i3).setTitleFontSize(18); c.getPad(i3).getAxisZ().setRange(min,max);
    		c.draw((H2F)tl.Timeline.getItem((i+1)*10,0));c.draw(line1);c.draw(line2);c.draw(line3);c.draw(line4);
             
    		c.cd(i3+1); c.getPad(i3+1).setAxisRange(-0.5,runIndex,min,max); c.getPad(i3+1).setTitleFontSize(18);
    		List<GraphErrors> gglist = getGraph(((H2F)tl.Timeline.getItem((i+1)*10,0)),((H2F)tl.Timeline.getItem((i+1)*10,1)),is-1); 
               
    		for (int ii=1; ii<gglist.size(); ii++) {    
        		gglist.get(ii).setTitleX("Run Index"); gglist.get(ii).setTitleY("Sector "+is+((i==1)?" Pi0 #sigma(E)/E":" Pi0 Mean/Mass"));
    			c.draw(gglist.get(ii),(ii==1)?" ":"same"); c.draw(line5);	
    		}
    		g2 = new GraphErrors(); g2.setMarkerSize(5); g2.setMarkerColor(4); g2.setLineColor(2);
    		g2.addPoint(runIndexSlider,gglist.get(0).getDataY(runIndexSlider),0,0); c.draw(g2,"same");
    		
    		c.cd(i3+2); c.getPad(i3+2).getAxisY().setRange(0.,fd.getGraph().getMax()*1.1);
            fd.getHist().getAttributes().setOptStat("1000100");
            DataLine line6 = new DataLine(mpi0,-50,mpi0,fd.getGraph().getMax()*1.5); line6.setLineColor(3); line6.setLineWidth(2);
            c.draw(fd.getHist()); c.draw(fd.getGraph(),"same"); c.draw(line6);
        }
    } 
    
    @Override 
    public void timerUpdate() {
    	
    }
}
