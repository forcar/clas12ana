package org.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.clas.tools.DataProvider;
import org.clas.tools.FitData;
import org.clas.tools.TOFPaddle;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.CalorimeterResponse;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
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
	
	ECPart part = new ECPart();
    List<TOFPaddle>     paddleList = null;
    List<List<DetectorResponse>>   res = new ArrayList<List<DetectorResponse>>();   
    List<DetectorResponse> ecClusters  = null;
    Map<String,Integer> smap = new HashMap<String,Integer>();
    String[]         det = new String[]{"pcal","ecin","ecou"};
    String[]           v = new String[]{"u","v","w"};
    
    IndexedList<FitData>  Fits = new IndexedList<FitData>(4);    
    Boolean                isAnalyzeDone = false;
    TreeFile tree = null; float[] rowf = null; int nevent=0;
    
    double pcx,pcy,pcz,refE=0,refP=0,refTH=25;    
    int[]   iidet = {1,4,7};
    float ipU,ipV,ipW;   

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
                                 "MCPHOT");
        
        this.usePCCheckBox(true);
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init();
        localinit();
    }
    
    public void localinit() {
        configEngine("pi0");
        part.setGeom("2.5");  
        part.setConfig("pi0");  
        part.setGoodPhotons(1212);   
        openTree();
    }
    
    public void openTree() {
    	tree  = new TreeFile(jawPath+"mctrue.hipo","T","ev:avgT:avgZ:vz:mvz:aR:vR:mvR:tedep");
    	rowf = new float[9];    	
    }
    
    @Override
    public void createHistos(int run) { 
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
    	create1DHistos(11,100,0.,2500,"ftof","GMEAN (ADC ch.)");
    	createMCHistos(12);
    }
    
    @Override        
    public void plotHistos(int run) {
        setRunNumber(run);
        plotUVW(0);
        plotPI0Summary(1);
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
        plotMCPHOT(12);
    }
    
    public void createPHOTHistos(int k) {
    	
 	   int run = getRunNumber();
       H1F h; 
       int is=0, n=0;
       String tag = is+"_"+n+"_"+k+"_"+run;
       h = new H1F("pi0_pcal_u_"+tag,"pi0_pcal_u_"+tag, 100, 0., 2.);
       DataGroup dg = new DataGroup(1,1);
       dg.addDataSet(h, 0);
       this.getDataGroup().add(dg,is,n,k,run);        	
    }
    
    public void createUVWHistos(int k, int n, int nch, double x1, double x2, String txt) {
    	
	   int run = getRunNumber();
       H2F h; 
   
       for (int is=1; is<7; is++) {
           String tag = is+"_"+n+"_"+k+"_"+run;
           DataGroup dg = new DataGroup(3,3);
           h = new H2F("pi0_pcal_u_"+tag,"pi0_pcal_u_"+tag, nch, x1, x2, 68, 1., 69.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("PCAL U Strips"); 
           dg.addDataSet(h,0);  
           h = new H2F("pi0_pcal_v_"+tag,"pi0_pcal_v_"+tag, nch, x1, x2, 62, 1., 63.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("PCAL V Strips");        
           dg.addDataSet(h,1);            
           h = new H2F("pi0_pcal_w_"+tag,"pi0_pcal_w_"+tag, nch, x1, x2, 62, 1., 63.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("PCAL W Strips");  
           dg.addDataSet(h,2); 
       
           h = new H2F("pi0_ecin_u_"+tag,"pi0_ecin_u_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("ECIN U Strips");    
           dg.addDataSet(h,3);  
           h = new H2F("pi0_ecin_v_"+tag,"pi0_ecin_v_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("ECIN V Strips");        
           dg.addDataSet(h,4);            
           h = new H2F("pi0_ecin_w_"+tag,"pi0_ecin_w_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("ECIN W Strips");  
           dg.addDataSet(h,5); 
       
           h = new H2F("pi0_ecou_u_"+tag,"pi0_ecou_u_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("ECOU U Strips");    
           dg.addDataSet(h,6);  
           h = new H2F("pi0_ecou_v_"+tag,"pi0_ecou_v_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("ECOU V Strips");        
           dg.addDataSet(h,7);            
           h = new H2F("pi0_ecou_w_"+tag,"pi0_ecou_w_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+txt); h.setTitleY("ECOU W Strips");  
           dg.addDataSet(h,8);   
           this.getDataGroup().add(dg,is,n,k,run);
       }            
   }
    
    public void create1DHistos(int k, int nch, double x1, double x2, String var, String txt) {
    	
	    int run = getRunNumber();
        H1F h; 
        DataGroup dg = new DataGroup(3,2);
        
        for (int is=1; is<7; is++) {
            String tag = var+"_"+is+"_"+k+"_"+run;
            h = new H1F("pi0_"+tag,"pi0_"+tag, nch, x1, x2);
            h.setTitleX("Sector "+is+" "+txt);  
            dg.addDataSet(h,is-1);   
        }
        this.getDataGroup().add(dg,0,0,k,run);
    }
    
    public void create2DHistos(int k, int nchx, double x1, double x2, int nchy, double y1, double y2, String var, String txtx, String txty) {
    	
	    int run = getRunNumber();
        H2F h;     	
        DataGroup dg = new DataGroup(3,2);
        
        for (int is=1; is<7; is++) {
            String tag = var+"_"+is+"_"+k+"_"+run;
            h = new H2F("pi0_"+tag,"pi0_"+tag, nchx, x1, x2, nchy, y1, y2);
            h.setTitleX("Sector "+is+" "+txtx); h.setTitleY(txty); 
            dg.addDataSet(h,is-1);   
        }
        this.getDataGroup().add(dg,0,0,k,run);    	
    }
    
    public void createXYHistos(int k, int n, int nb, int bmx) {
   	 
	    int run = getRunNumber();
        H2F h;  
        
        String[] t = {"e","w","r"};
        
        for (int i=0; i<3; i++) {
           DataGroup dg = new DataGroup(3,2);
           for (int d=0; d<3; d++) {
              h = new H2F("pi0_"+det[d]+"_xy_"+t[i]+"_"+n+"_"+k+"_"+run,"hi_"+det[d]+"_xy_"+t[i]+"_"+k+"_"+run,nb,-bmx,bmx,nb,-bmx,bmx);
              h.setTitleX("Photon "+n+" X(cm)"); h.setTitleY("Photon "+n+" Y(cm)");
              dg.addDataSet(h,d);  
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
        H1F h; 
        int n = 0;
        DataGroup dg = new DataGroup(3,5);
        
        for (int is=1; is<7; is++) {   
            for (int j=is+1; j<7; j++) {
                String tag = var+"_"+is+"_"+j+"_"+k+"_"+run;
                h = new H1F("pi0_"+tag,"pi0_"+tag, nch, x1, x2);
                h.setTitleX("Sector Pair "+is+j+" "+txt);  
                dg.addDataSet(h,n);  n++; 
                smap.put(is+"_"+j,n);
            }
        }
        this.getDataGroup().add(dg,0,0,k,run);    	
    }
    
    public void createMCHistos(int k) {
 	   int run = getRunNumber();
       H1F h1;
       H2F h2; 
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
    	
    	double tmax=30;
    	
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
            rowf[0] = nevent++;
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
        
        dropBanks(event);

        ecClusters = part.readEC(event);  
        
        if (ecClusters.size()==0) return;
        
        DataBank ecBank = event.getBank("ECAL::clusters");
               
        if(event.hasBank("MC::Particle")==true) processMC(event, ecBank);
        
        if (true) return;
        
        double[] esum = {0,0,0,0,0,0};
        int[][] nesum = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};    
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
                  int toflay = bank.getByte("layer", i);
                  int   isec = bank.getByte("sector", i);
                  part.mip[isec-1] = (toflay<3&&bank.getFloat("energy",i)>thresh[toflay-1]) ? 1:0;
               }
            }      
        } 
        
        part.getNeutralResponses(ecClusters);
        
        for (int is=1; is<7; is++) {
           
            if (part.mip[is-1]!=1) {  // No FTOF MIP in sector
//            	((H1F) this.getDataGroup().getItem(0,0,12,run).getData(0).get(0)).fill(part.e1);
//          	    System.out.println(part.e1+" "+Math.acos(part.cth1)*180/3.14159);
                double invmass = Math.sqrt(part.getTwoPhotonInvMass(is));
                double     opa = Math.acos(part.cth)*180/3.14159;
                
//              boolean badPizero = part.X>0.5 && opa<8;
                boolean badPizero = part.X>1 && opa<0;
        
                if(invmass>0&&part.iis[0]>0&&part.iis[1]>0&&!badPizero) {                                                    
                    if(part.iis[0]< part.iis[1]) ((H1F) this.getDataGroup().getItem(0,0,2,run).getData(smap.get(part.iis[0]+"_"+part.iis[1])-1).get(0)).fill(invmass*1e3);   
                    
                    if(part.iis[0]==part.iis[1]) {
                    	
                    	((H1F) this.getDataGroup().getItem(0,0,1,run).getData(part.iis[0]-1).get(0)).fill(invmass*1e3);   
                    	double tdif = part.t[0][0]-part.t[1][0];
                    	((H2F) this.getDataGroup().getItem(0,0,5,run).getData(part.iis[0]-1).get(0)).fill((invmass-part.mpi0)*1e3,opa);                    
                    	((H2F) this.getDataGroup().getItem(0,0,6,run).getData(part.iis[0]-1).get(0)).fill((invmass-part.mpi0)*1e3,part.e1c*part.e2c);                    
                    	((H2F) this.getDataGroup().getItem(0,0,7,run).getData(part.iis[0]-1).get(0)).fill((invmass-part.mpi0)*1e3,Math.sqrt(part.tpi2));                    
                    	((H2F) this.getDataGroup().getItem(0,0,8,run).getData(part.iis[0]-1).get(0)).fill((invmass-part.mpi0)*1e3,tdif);                    
                        if(nesum[0][is-1]>1 && nesum[1][is-1]>0) {
//                            ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,4,0).fill(esum[is-1],7,1.);          // Total Cluster Energy            
//                            ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,4,2).fill(part.e1,part.SF1,1.);      // S.F. vs. meas. photon energy            
                    	   ((H2F) this.getDataGroup().getItem(0,0,3,run).getData(part.iis[0]-1).get(0)).fill(opa,part.e1c*part.e2c);
                       	   ((H2F) this.getDataGroup().getItem(0,0,4,run).getData(part.iis[0]-1).get(0)).fill(opa,part.X);
                        }            
                          
                        if (invmass>0.100 && invmass<0.180) {
                        	((H2F) this.getDataGroup().getItem(0,0,9,run).getData(part.iis[0]-1).get(0)).fill(Math.acos(part.cpi0)*180/Math.PI,Math.sqrt(part.tpi2));                    

/*                            ecPix[0].strips.hmap1.get("H1_a_Hist").get(is,4,0).fill((float)(1e3*(Math.sqrt(part.tpi2)-refE))); // Pizero total energy error
                            ecPix[0].strips.hmap1.get("H1_a_Hist").get(is,4,1).fill(Math.acos(part.cpi0)*180/3.14159-refTH);     // Pizero theta angle error
                            ecPix[0].strips.hmap1.get("H1_a_Hist").get(is,4,2).fill((float)part.X);                              // Pizero energy asymmetry
                            ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,4,4).fill(opa,(float)part.X);      
                            ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance11,1,1.); // Pizero photon 1 PCAL-ECinner cluster error
                            ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance12,2,1.); // Pizero photon 2 PCAL-ECinner cluster error
                            ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance21,3,1.); // Pizero photon 1 PCAL-ECouter cluster error
                            ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance22,4,1.); // Pizero photon 2 PCAL-ECouter cluster error   
                            */   
                            for (int il=0; il<3; il++) {
                            	for (int im=0; im<2; im++) {
                            		if(ecBank!=null) {
                                      ipU = (ecBank.getInt("coordU", part.iip[im][il])-4)/8+1;
                                      ipV = (ecBank.getInt("coordV", part.iip[im][il])-4)/8+1;
                                      ipW = (ecBank.getInt("coordW", part.iip[im][il])-4)/8+1;
                                      ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(0+il*3).get(0)).fill(invmass*1e3,ipU);
                                      ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(1+il*3).get(0)).fill(invmass*1e3,ipV);
                                      ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(2+il*3).get(0)).fill(invmass*1e3,ipW);
                            		}
                                    ((H2F) this.getDataGroup().getItem(il,im+1,10,run).getData(0).get(0)).fill(-part.x[im][il], part.y[im][il],1.);
                                    ((H2F) this.getDataGroup().getItem(il,im+1,10,run).getData(1).get(0)).fill(-part.x[im][il], part.y[im][il],invmass/part.mpi0);
                            	}
                            }
                        }
                    }
                
            }
            }
        
        }
                
    }
    
    @Override
    public void plotEvent(DataEvent de) {
    	    analyze();
    }

    public void analyze() {    
        System.out.println("I am in analyze()");
        analyzeGraphs(1,7,0,0,0,0,"c");
        fillTimeLineGraph();
        tree.close();
        System.out.println("Finished");
        isAnalyzeDone = true;
    }
    
    public void analyzeGraphs(int is1, int is2, int id1, int id2, int il1, int il2, String ro) {
    	FitData fd = null ; H2F h2 = null;
        int run=getRunNumber();
        for (int is=is1; is<is2; is++) {
        	H1F h1 = (H1F) this.getDataGroup().getItem(0,0,1,run).getData(is-1).get(0);
        	fd = new FitData(h1.getGraph()); fd.setInt((int)h1.getIntegral());
        	fd.setHist(h1);
            fd.graph.getAttributes().setTitleX(h1.getTitleX()); 
            fd.hist.getAttributes().setTitleX(h1.getTitleX()); 
            fd.initFit(0,40.,200.); fd.fitGraph(""); Fits.add(fd,is,0,0,run); 
        }
    }
    
    public void plotXYSummary(int index) {        
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        
        int   run = getRunNumber();
        
        c.clear(); c.divide(3,2);
        
        for (int idet=0; idet<3; idet++) {
           H2F h1 = (H2F) this.getDataGroup().getItem(idet,getActivePC(),index,run).getData(0).get(0); 
           H2F h2 = (H2F) this.getDataGroup().getItem(idet,getActivePC(),index,run).getData(1).get(0);  
           H2F h3 = (H2F) this.getDataGroup().getItem(idet,getActivePC(),index,run).getData(2).get(0);  
            
           h3 = h2.divide(h2, h1); 
           c.cd(idet);   c.getPad(idet).getAxisZ().setLog(true); c.draw(h1);
           c.cd(idet+3); c.getPad(idet+3).getAxisZ().setLog(false); c.getPad(idet+3).getAxisZ().setRange(0.7, 1.3); c.draw(h3);
        }	
   
    }   
    
    public void plotMCPHOT(int index) {
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        int   run = getRunNumber();
        
        c.clear(); c.divide(3,5);
        
        String dname[] = {"PCAL ","ECIN ","ECOU "};
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
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),getActivePC(),index,getRunNumber()));
    }   
    
    public void plotPI0Summary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));
    }
    
    public void plotMCSummary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,getActiveLayer(),index,getRunNumber()));
    }
    
    @Override 
    public void timerUpdate() {
    	
    }
}
