package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;

import org.clas.tools.DataProvider;
import org.clas.tools.TOFPaddle;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

public class ECpi0 extends DetectorMonitor{
	
	ECPart part = new ECPart();
    List<TOFPaddle>     paddleList = null;
    List<List<DetectorResponse>>   res = new ArrayList<List<DetectorResponse>>();    
    int[]   iidet = {1,4,7};
	
    public ECpi0(String name) {
        super(name);
        this.setDetectorTabNames("IVM",
                                 "SIJ");        
        this.usePCCheckBox(true);
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init();
    }
    
    @Override
    public void init() {
        configEngine("pi0");
        part.setGeom("2.5");  
        part.setConfig("pi0");  
        part.setGoodPhotons(1212);    	
    }
    
    @Override
    public void createHistos(int run) { 
    	createIVMHistos(0, 75,5.,400.,"invm","Two Photon Invariant Mass (MeV)"); // Sector ij photons i==j
    	createSIJHistos(1,130,5.,700.,"invm","Two Photon Invariant Mass (MeV)"); // Sector ij photons i!=j
    }
    
    @Override        
    public void plotHistos(int run) {
        setRunNumber(run);
        plotPI0Summary(0);
        plotPI0Summary(1);    	
    }
    
    public void createIVMHistos(int k, int nch, double x1, double x2, String var, String txt) {
    	
	    int run = getRunNumber();
        H1F h; 
        int n = 0;
        DataGroup dg = new DataGroup(3,2);
        
        for (int is=1; is<7; is++) {
            String tag = var+"_"+is+"_"+k+"_"+run;
            h = new H1F("pi0_"+tag,"pi0_"+tag, nch, x1, x2);
            h.setTitleX("Sector "+is+" "+txt);  
            dg.addDataSet(h,is-1);   
        }
        this.getDataGroup().add(dg,0,0,k,run);
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
            }
        }
        this.getDataGroup().add(dg,0,0,k,run);    	
    }
    
    @Override
    public void processEvent(DataEvent event) {   
        int run = getRunNumber();
 	   
        if(event.hasBank("ECAL::hits")) {
           event.removeBank("ECAL::hits");        
           event.removeBank("ECAL::peaks");        
           event.removeBank("ECAL::clusters");        
           event.removeBank("ECAL::calib");
           event.removeBank("ECAL::moments");
        } 
        
        engine.processDataEvent(event); 
        
        DataBank ecBank = event.getBank("ECAL::clusters");
        
        List<DetectorResponse> ecClusters = part.readEC(event);  
        
        double[] esum = {0,0,0,0,0,0};
        int[][] nesum = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};    
  	    double pcalE[] = new double[6];
        
        if (ecClusters.size()==0) return;
        
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
                    part.mip[isec-1] = (paddle.geometricMean()>thresh[toflay-1]) ? 1:0;
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
                
                boolean badPizero = part.X>0.5 && opa<8;
                if(part.iis[0]>0&&part.iis[1]>0&&!badPizero) {
                                  
                    if(part.iis[0]<=part.iis[1]) ((H1F) this.getDataGroup().getItem(part.iis[0],0,0,run).getData(part.iis[1]).get(0)).fill((float)invmass*1e3);
                    if(part.iis[0]> part.iis[1]) ((H1F) this.getDataGroup().getItem(part.iis[1],0,0,run).getData(part.iis[0]).get(0)).fill((float)invmass*1e3);
/*                
                if(part.iis[0]==part.iis[1]) {
                    
                if(nesum[0][is-1]>1 && nesum[1][is-1]>0) {
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,4,0).fill(esum[is-1],7,1.);          // Total Cluster Energy            
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,4,2).fill(part.e1,part.SF1,1.);      // S.F. vs. meas. photon energy            
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,7,2).fill(opa,part.e1c*part.e2c,1.); // E1*E2 vs opening angle            
                }            
                          
                if (invmass>0.100 && invmass<0.180) {
                    ecPix[0].strips.hmap1.get("H1_a_Hist").get(is,4,0).fill((float)(1e3*(Math.sqrt(part.tpi2)-refE))); // Pizero total energy error
                    ecPix[0].strips.hmap1.get("H1_a_Hist").get(is,4,1).fill(Math.acos(part.cpi0)*180/3.14159-refTH);   // Pizero theta angle error
                    ecPix[0].strips.hmap1.get("H1_a_Hist").get(is,4,2).fill((float)part.X);                            // Pizero energy asymmetry
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,4,4).fill(opa,(float)part.X);      
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance11,1,1.); // Pizero photon 1 PCAL-ECinner cluster error
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance12,2,1.); // Pizero photon 2 PCAL-ECinner cluster error
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance21,3,1.); // Pizero photon 1 PCAL-ECouter cluster error
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,9,4).fill(part.distance22,4,1.); // Pizero photon 2 PCAL-ECouter cluster erro    
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(1,9,5).fill(-part.x1, part.y1,1.);
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(1,9,6).fill(-part.x1, part.y1,invmass/part.mpi0);
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(1,9,7).fill(-part.x2, part.y2,1.);
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(1,9,8).fill(-part.x2, part.y2,invmass/part.mpi0);
                    float ipU,ipV,ipW;
                    ipU = (ecBank.getInt("coordU", part.iip[0])-4)/8;
                    ipV = (ecBank.getInt("coordV", part.iip[0])-4)/8;
                    ipW = (ecBank.getInt("coordW", part.iip[0])-4)/8;
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,10,1).fill(invmass*1e3,ipU);
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,10,2).fill(invmass*1e3,ipV);
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,10,3).fill(invmass*1e3,ipW);
                    ipU = (ecBank.getInt("coordU", part.iip[1])-4)/8;
                    ipV = (ecBank.getInt("coordV", part.iip[1])-4)/8;
                    ipW = (ecBank.getInt("coordW", part.iip[1])-4)/8;
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,10,4).fill(invmass*1e3,ipU);
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,10,5).fill(invmass*1e3,ipV);
                    ecPix[0].strips.hmap2.get("H2_a_Hist").get(is,10,6).fill(invmass*1e3,ipW);
                }
            }
*/                
            }
            }
        
        }
                
    }
    
    public void plotPI0Summary(int index) {
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));
    }
}
