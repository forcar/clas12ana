package org.clas.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.clas.tools.DataProvider;
import org.clas.tools.TOFPaddle;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import Jampack.Inv;

public class ECpi0 extends DetectorMonitor{
	
	ECPart part = new ECPart();
    List<TOFPaddle>     paddleList = null;
    List<List<DetectorResponse>>   res = new ArrayList<List<DetectorResponse>>();   
    List<DetectorResponse> ecClusters  = null;
    Map<String,Integer> smap = new HashMap<String,Integer>();
    String[]         det = new String[]{"pcal","ecin","ecou"};
    String[]           v = new String[]{"u","v","w"};
    
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
                                 "FTOF");
        
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
    }
    
    @Override
    public void createHistos(int run) { 
    	createPI0Histos(0, 1, 30,50.,250.); // invm vs. photon 1 strips
    	createPI0Histos(0, 2, 30,50.,250.); // invm vs. photon 2 strips
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
    }
    
    @Override        
    public void plotHistos(int run) {
        setRunNumber(run);
        plotPI0(0);
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
    }
    
    public void createPI0Histos(int k, int n, int nch, double x1, double x2) {
    	
	   int run = getRunNumber();
       H2F h; 
   
       for (int is=1; is<7; is++) {
           String tag = is+"_"+n+"_"+k+"_"+run;
           DataGroup dg = new DataGroup(3,3);
           h = new H2F("pi0_pcal_u_"+tag,"pi0_pcal_u_"+tag, nch, x1, x2, 68, 1., 69.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("PCAL U Strips"); 
           dg.addDataSet(h,0);  
           h = new H2F("pi0_pcal_v_"+tag,"pi0_pcal_v_"+tag, nch, x1, x2, 62, 1., 63.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("PCAL V Strips");        
           dg.addDataSet(h,1);            
           h = new H2F("pi0_pcal_w_"+tag,"pi0_pcal_w_"+tag, nch, x1, x2, 62, 1., 63.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("PCAL W Strips");  
           dg.addDataSet(h,2); 
       
           h = new H2F("pi0_ecin_u_"+tag,"pi0_ecin_u_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("ECIN U Strips");    
           dg.addDataSet(h,3);  
           h = new H2F("pi0_ecin_v_"+tag,"pi0_ecin_v_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("ECIN V Strips");        
           dg.addDataSet(h,4);            
           h = new H2F("pi0_ecin_w_"+tag,"pi0_ecin_w_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("ECIN W Strips");  
           dg.addDataSet(h,5); 
       
           h = new H2F("pi0_ecou_u_"+tag,"pi0_ecou_u_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("ECOU U Strips");    
           dg.addDataSet(h,6);  
           h = new H2F("pi0_ecou_v_"+tag,"pi0_ecou_v_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("ECOU V Strips");        
           dg.addDataSet(h,7);            
           h = new H2F("pi0_ecou_w_"+tag,"pi0_ecou_w_"+tag, nch, x1, x2, 36, 1., 37.);
           h.setTitleX("Sector "+is+" Inv. Mass (MeV)"); h.setTitleY("ECOU W Strips");  
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

    @Override
    public void processEvent(DataEvent event) {  
    	
        int run = getRunNumber();
     
        dropBanks(event);

        if(!event.hasBank("ECAL::clusters")) return;
        
        DataBank ecBank = event.getBank("ECAL::clusters");
        
        ecClusters = part.readEC(event);  
        
        double[] esum = {0,0,0,0,0,0};
        int[][] nesum = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};    
  	    double pcalE[] = new double[6];
        
        if (ecClusters.size()==0) return;
        
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
                                    ipU = (ecBank.getInt("coordU", part.iip[im][il])-4)/8+1;
                                    ipV = (ecBank.getInt("coordV", part.iip[im][il])-4)/8+1;
                                    ipW = (ecBank.getInt("coordW", part.iip[im][il])-4)/8+1;
                                    ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(0+il*3).get(0)).fill(invmass*1e3,ipU);
                                    ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(1+il*3).get(0)).fill(invmass*1e3,ipV);
                                    ((H2F) this.getDataGroup().getItem(part.iis[im],im+1,0,run).getData(2+il*3).get(0)).fill(invmass*1e3,ipW);
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
    
    public void plotPI0(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),getActivePC(),index,getRunNumber()));
    }   
    
    public void plotPI0Summary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));
    }
    
    @Override 
    public void timerUpdate() {
    	
    }
}
