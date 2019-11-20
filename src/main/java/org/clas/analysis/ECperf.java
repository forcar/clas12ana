package org.clas.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.clas.tools.Event;
import org.clas.tools.ParallelSliceFitter;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECperf extends DetectorMonitor {
	
	Event    ev = new Event();
    ECPart part = new ECPart();  
	
	DataGroup dg = null;
	int nc=0;
	
	public boolean goodPROT,goodPBAR,goodPIP,goodPIM,goodNEUT,goodPHOT,goodPHOTR,goodPHOT2,goodPIPP,goodPI0;
	String[]   det = new String[]{"PCAL ","ECIN ","ECOU "};
	String[] scdet = new String[]{"P1A ","P1B ","P2 "};
	String[]   xyz = new String[]{"X","Y","Z"}; 
	
	public int Nevts, Nelecs, Ntrigs, runNum;
	public int[] Ntrigs_sect = new int[6];
    public int[] Nelecs_sect = new int[6];
	float n_neut_ecal = 0, n_neut_mm = 0, n_neut_mm_save=0;

	public float EB, Eb, Mp;
	public float RFT, STT;
	public long TriggerWord;
	public float rfPeriod;
	public int rf_large_integer;
	public float lU,lV,lW,cZ;
	public int iU,iV,iW;
	
	public LorentzVector VB, VT, Ve, VGS, Vprot, Vpbar, Vpip, Vpim;
	public int   e_part_ind, e_sect, e_FTOF_pad1a, e_FTOF_pad1b, e_HTCC_bin_phi, e_HTCC_bin_theta;
	public float e_mom, e_the, e_phi, e_vx, e_vy, e_vz, e_cz, e_x, e_y;
	public float e_xB, e_Q2, e_W;

    public float e_ecal_esum,e_ecal_pcsum,e_ecal_ecsum;
        
	public int   prot_part_ind;
	public float prot_mom, prot_the, prot_phi, prot_vx, prot_vy, prot_vz, prot_beta;
	public float pbar_mom, pbar_the, pbar_phi, pbar_vx, pbar_vy, pbar_vz, pbar_beta;

	public int   pim_part_ind, pip_part_ind, pip_FTOF_pad1b;
	public float pip_mom, pip_the, pip_phi, pip_vx, pip_vy, pip_vz, pip_beta;
	public float pip_ecal_esum;
	
	public float pim_mom, pim_the, pim_phi, pim_vx, pim_vy, pim_vz, pim_beta;

	public int   G1_part_ind, G2_part_ind, G1_pcal_ind, G2_pcal_ind, G1_cal_layers, G2_cal_layers;
	public int   G1_sec,G2_sec,G1_lay,G2_lay;
	public float G1_mom, G1_e, G1_the, G1_phi, G2_mom, G2_e, G2_the, G2_phi;

	public float elast_dPhi, elast_EB;
	public float epip_dPhi, epip_MM, ep_dPhi, epbar_dPhi, ep_MM, epbar_MM;
	public float pi0_mass, pi0_mom, pi0_e, pi0_the, pi0_phi, pi0_open, pi0_cx, pi0_cy;
	
	public float ecal_pi0_mass, ecal_pi0_mom, ecal_pi0_e, ecal_pi0_the, ecal_pi0_phi;
	public float ecal_pi0_opa, ecal_pi0_cx, ecal_pi0_cy, ecal_pi0_X;	
	public float neut_mom,neut_the,neut_phi,neut_cx,neut_cy;
	public float ecal_neut_the,ecal_neut_phi,ecal_neut_beta,ecal_neut_cx,ecal_neut_cy;
	public float ecal_phot_the,ecal_phot_phi,ecal_phot_beta,ecal_phot_nrg;
	public int   ecal_neut_sec, ecal_phot_sec, ecal_pi0_sec;
	public int[] ecal_neut_esum = new int[6];
	public int[] ecal_phot_esum = new int[6];
	
	public IndexedList<Float> elec_ecal_resid = new IndexedList<Float>(3);
	public List<Particle> pim_ecal  = new ArrayList<Particle>();
	public List<Particle> pip_ecal  = new ArrayList<Particle>();
	public List<Particle> prot_ecal = new ArrayList<Particle>();
	public List<Particle> pbar_ecal = new ArrayList<Particle>();
	public List<Particle> phot_ecal = new ArrayList<Particle>();
	public List<Particle> neut_ecal = new ArrayList<Particle>();
	public List<NeutralMeson> pi0_ecal = new ArrayList<NeutralMeson>();
	public IndexedList<Float> ecal_rad   = new IndexedList<Float>(3);
	public IndexedList<Float> elec_ftof_resid = new IndexedList<Float>(3);
	
    public float[][]  counter = new float[6][3];
   
	public IndexedTable rfTable;	
	GraphErrors neuteff = null;
	
    public ECperf(String name) {
        super(name);
        this.setDetectorTabNames("ECkin",
        						 "ECelec",
                				 "SCelec",
                				 "ECprot",
                				 "ECpbar",
        		                 "ECpip",
              				     "ECpim",
        		                 "ECpi0",
        		                 "ECneut",
        		                 "ECphot",
        		                 "ECtime");

        this.useCALButtons(true);
        this.use123Buttons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();      
    }
    
    public void localinit() {
    	System.out.println("ECperf.localinit()");
        configEngine("muon"); 
    	tl.setFitData(Fits);  
        part.setGeom("2.5");  
        part.setConfig("pi0");  
        part.setGoodPhotons(1212);    	
        neuteff = getGraph(outPath+"files/neuteff.vec",50);
    }  
    
    public void localclear() {
    	System.out.println("ECperf:localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	FitSummary.clear();
    	Fits.clear();
    	tl.Timeline.clear();
    	slider.setValue(0);
    } 
    
    public void dstinit(int run) {
	    System.out.println("ECperf:dstinit("+run+")");
		runNum = run;
		Nevts=0;Nelecs=0;Ntrigs=0;
		Ntrigs_sect = new int[6];
		Nelecs_sect = new int[6];
		for(int s=0;s<6;s++){Ntrigs_sect[s]=0;Nelecs_sect[s]=0;}
//		trigger_bits = new boolean[32];
		Mp = 0.93827f;
        Eb = EB = getBeamEnergy(run);
        System.out.println("Eb="+Eb+" run="+runNum);
		rfPeriod = 4.008f;
        engine.getConstantsManager().init(Arrays.asList(new String[]{"/daq/tt/fthodo","/calibration/eb/rf/config"}));
        rfTable = engine.getConstantsManager().getConstants(runNum,"/calibration/eb/rf/config");
        if (rfTable.hasEntry(1, 1, 1)){
            System.out.println(String.format("RF period from ccdb for run %d: %f",runNum,rfTable.getDoubleValue("clock",1,1,1)));
            rfPeriod = (float)rfTable.getDoubleValue("clock",1,1,1);
        }
		rf_large_integer = 1000;    	
        VB = new LorentzVector(0,0,Eb,Eb);
        VT = new LorentzVector(0,0,0,Mp);    
    }
    
    public void createECkin(int st) {
    	
    	String tab = "ECkin", tag = null;
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
		
    	switch (st) {
        
        case 0: 
        dg = new DataGroup(6,4);
		for(int is=1;is<7;is++){    	
	        tag = is+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_1_",tag,60,5,12,60,EB-1,EB,   "","#theta (^o)","p (GeV)"),is-1);
	        dg.addDataSet(makeH2(tab+"_2_",tag,60,0.6,1.2,60,EB-1,EB,"","W (GeV)","p (GeV)"),is-1+6);
	        dg.addDataSet(makeH2(tab+"_3_",tag,60,0.6,1.2,60,5,12,   "","W (GeV)","#theta (^o)"),is-1+12);
	        dg.addDataSet(makeH2(tab+"_4_",tag,60,0.6,1.2,60,-20,30, "","W (GeV)","#phi (^o"),is-1+18);
		}
		break;
		
        case 1:
        dg = new DataGroup(6,3);
        for(int is=1; is<7; is++) {
	        tag = is+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH1(tab+"_1_",tag,100,0,450,"Sector "+is,"LU (cm)"),is-1);
	        dg.addDataSet(makeH1(tab+"_1_",tag,100,0,450," ","LV (cm)"),is-1+6);
	        dg.addDataSet(makeH1(tab+"_1_",tag,100,0,450," ","LW (cm)"),is-1+12);
        	
        }
        break;
        
        case 2:
        dg = new DataGroup(3,2);
        for(int iv=0; iv<1; iv++) {
	        tag = iv+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_1_",tag,200,-400,400,200,-400,400,det[iv],"X (CM)","Y(CM)"),iv);        	
        }
        
        
    	}
    	
    	this.getDataGroup().add(dg,0,st,k,run);      
    	
    }
    
    public void createECelec(int st) {
    	
    	String tab = "ECelec", tag = null;
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	 
		F1D     f1 = null;
		
    	switch (st) {
        
        case 0: 
        dg = new DataGroup(6,3);
		for(int is=1;is<7;is++){
	        tag = is+"_"+st+"_"+k+"_"+run;
			dg.addDataSet(makeH2(tab+"_1_",tag,60,0,EB,60,0.15,0.35,"","p (GeV)","E/P"),is-1);
			dg.addDataSet(makeH2(tab+"_2_",tag,68,1,69,50,0.15,0.35,"","PCAL U STRIP","E/P"),is-1+6);
			dg.addDataSet(makeH2(tab+"_3_",tag,70,0,EB/6,70,0,EB/5, "","PC (GeV)","EC (GeV)"),is-1+12);	
		}
		
		break;
		
        case 1:        
//		f1 = new F1D("H_e_EC_resid_f+"+run,"[a]",5,35); //angle
		f1 = new F1D("H_e_EC_resid_f+"+run,"[a]",0,EB);  //momentum
//		f1 = new F1D("H_e_EC_resid_f+"+run,"[a]",0.8,1.0);  //cz
		f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);		    	
		for(int i=0;i<3;i++) { //pcal,ecin,ecou
//			float ylim = (i==0)?5:10;
			int inn=0; dg = new DataGroup(6,3);        			 
			for(int n=0; n<3; n++) { //x,y,z
				float ylim1 = (n<2)?5:1; float ylim2=(n<2)?5:1;
				for(int is=1;is<7;is++){  //sector 
					tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
//					dg.addDataSet(makeH2(tag,"H_e_EC_resid_",60,5,35,40,-ylim,ylim,"","S"+is+" #thete_e","DC"+xyz[n]+"-"+det[i]),inn);
					dg.addDataSet(makeH2(tab+"_1_",tag,60,0,EB,40,-ylim1,ylim2,"","S"+is+" p_e","DC"+xyz[n]+"-"+det[i]),inn);
//					dg.addDataSet(makeH2(tag,"H_e_EC_resid_",60,0,EB,40,-ylim1,ylim2,"","S"+is+" cz_e","DC"+xyz[n]+"-"+det[i]),inn);
					dg.addDataSet(f1, inn); inn++; 
				}
			}
		 	this.getDataGroup().add(dg,i,st,k,run);
		}
	 	return;		
		
        case 2:
        dg = new DataGroup(6,6);        
        String lab2[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};
		f1 = new F1D("H_e_EC_rad2_f+"+run,"[a]",5,35); 
        f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
		for(int i=0;i<3;i++) { //pcal,ecin,ecou
			for(int n=0; n<2; n++) { //thet,phi
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,5,35,60,-5,20,"","S"+is+" #theta_e",lab2[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
			}
		}	
		}        	
        
    	break;
    	
    	case 3:
    	dg = new DataGroup(6,4);    	
        String lab3[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};    
    	f1 = new F1D("H_e_EC_rad3_f+"+run,"[a]",0.0,0.5); 
    	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
    	for(int i=0;i<3;i++) { //pcal,ecin,ecou
    		for(int n=0; n<2; n++) { //thet,phi
    			for(int is=1;is<7;is++){  //sector  	
    				tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,0,0.5,60,-5,20,"","S"+is+" "+det[i]+" E",lab3[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
    			}
    			}	
    		}
    	
    	break;
    	
    	case 4:
    	dg = new DataGroup(6,4);    	
        String lab4[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};    
    	f1 = new F1D("H_e_EC_rad4_f+"+run,"[a]",0.0,EB); 
    	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
    	for(int i=0;i<3;i++) { //pcal,ecin,ecou
    		for(int n=0; n<2; n++) { //thet,phi
    			for(int is=1;is<7;is++){  //sector  	
    				tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,0,EB,60,-5,20,"","S"+is+" p_e",lab4[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
    			}
    			}	
    		}
    	
    	break;
    	
    	case 5:
        dg = new DataGroup(3,2); GraphErrors g = null;
    	f1 = new F1D("H_e_EC_zero_f+"+run,"[a]",0,7); 
    	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
        for (int i=0; i<6; i++) { //pcal,ecin,ecou (resid) pcal,ecin,ecou (rad) 
			tag = i+"_"+"1"+"_"+st+"_"+k+"_"+run;
        	g = new GraphErrors(tag); g.setMarkerSize(5); g.setMarkerColor(1); g.setLineColor(1);
        	dg.addDataSet(g, i); 
			tag = i+"_"+"2"+"_"+st+"_"+k+"_"+run;
        	g = new GraphErrors(tag); g.setMarkerSize(5); g.setMarkerColor(2); g.setLineColor(2);
        	dg.addDataSet(g, i); 
        	dg.addDataSet(f1, i);
        	
        }
        
    	}
    	
    	this.getDataGroup().add(dg,0,st,k,run);      
	
    } 
    
    public void createSCelec(int st) {
    	
    	String tab = "SCelec", tag = null;
    	int    run=getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	
		F1D     f1 = null;    
		
    	switch (st) {
    	
        case 1:        
		f1 = new F1D("H_e_SC_resid_f+"+run,"[a]",5,35); 
		f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);		    	
		for(int i=0;i<2;i++) { //p1a,p1b
	        int inn = 0; dg = new DataGroup(6,3);        
			for(int n=0; n<3; n++) { //x,y,z
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;					
					dg.addDataSet(makeH2(tab+"_1_",tag,60,5,35,40,-5,5,"","S"+is+" #theta_e","DC"+xyz[n]+"-"+scdet[i]), inn);
					dg.addDataSet(f1, inn); inn++;  
				}
			}	
		 	this.getDataGroup().add(dg,i,st,k,run);
		}
		
		return;
		
    	case 5:
            dg = new DataGroup(2,2); GraphErrors g = null;
        	f1 = new F1D("H_e_SC_zero_f+"+run,"[a]",0,7); 
        	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
            for (int i=0; i<2; i++) { //p1a,p1b
    			tag = i+"_"+"1"+"_"+st+"_"+k+"_"+run;
            	g = new GraphErrors(tag); g.setMarkerSize(5); g.setMarkerColor(1); g.setLineColor(1);
            	dg.addDataSet(g, i); 
    			tag = i+"_"+"2"+"_"+st+"_"+k+"_"+run;
            	g = new GraphErrors(tag); g.setMarkerSize(5); g.setMarkerColor(2); g.setLineColor(2);
            	dg.addDataSet(g, i); 
            	dg.addDataSet(f1, i);           	
            }
    	}
    	
    	this.getDataGroup().add(dg,0,st,k,run);      
    	
    }  
    
    public void createECprot(int st) {
    	
    	String tab = "ECprot", tag = null;   	
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	
    	switch (st) {
    	   
        case 0:
        dg = new DataGroup(6,2);
        for(int is=1;is<7;is++){ 
        	tag = is+"_"+st+"_"+k+"_"+run;        	
        	dg.addDataSet(makeH2(tab+"_0_",tag,50,-0.2,0.8,50,0.9,3,"S"+is,"MM^2 (GeV^2)","W (GeV)"),is-1);
        	dg.addDataSet(makeH1(tab+"_1_",tag,100,-0.2,0.8," ","MM^2 (GeV^2)"),is-1+6);
        	dg.addDataSet(makeH1(tab+"_2_",tag,100,-0.2,0.8," ","MM^2 (GeV^2)"),is-1+6);
        }
        break;
        case 1:
        dg = new DataGroup(4,3); int n=0;
    	tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,EB,100,0,40,        "","e- p (GeV)","e- #theta (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,4,100,0,50,         "","prot p (GeV)","prot #theta (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-15,15,    "","e vz (cm)","prot vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-180,180,100,-15,15,  "","#phi (^o)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,   0, 80,100,-15,15,  "","#theta (^o)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-15,15,    "","vz (cm)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-180,180,100,-180,180,"","#phi (^o)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,   0, 80,100,-180,180,"","#theta (^o)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-180,180,  "","vz (cm)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,7,100,0.2,1.05,     "","p (GeV)","#beta"),n);n++;
    	}    	
    	this.getDataGroup().add(dg,0,st,k,run);      
   	
    }
    
    public void createECpbar(int st) {
    	
    	String tab = "ECpbar", tag = null;   	
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	
    	switch (st) {
    	   
        case 0:
        dg = new DataGroup(6,2);
        for(int is=1;is<7;is++){ 
        	tag = is+"_"+st+"_"+k+"_"+run;        	
        	dg.addDataSet(makeH2(tab+"_0_",tag,50,-0.2,0.8,50,0.9,3,"S"+is,"MM^2 (GeV^2)","W (GeV)"),is-1);
        	dg.addDataSet(makeH1(tab+"_1_",tag,100,-0.2,0.8," ","MM^2 (GeV^2)"),is-1+6);
        	dg.addDataSet(makeH1(tab+"_2_",tag,100,-0.2,0.8," ","MM^2 (GeV^2)"),is-1+6);
        }
        break;
        case 1:
        dg = new DataGroup(4,3); int n=0;
    	tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,EB,100,0,40,        "","e- p (GeV)","e- #theta (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,4,100,0,50,         "","pbar p (GeV)","pbar #theta (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-15,15,    "","e vz (cm)","pbar vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-180,180,100,-15,15,  "","#phi (^o)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,   0, 80,100,-15,15,  "","#theta (^o)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-15,15,    "","vz (cm)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-180,180,100,-180,180,"","#phi (^o)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,   0, 80,100,-180,180,"","#theta (^o)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-180,180,  "","vz (cm)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,7,100,0.2,1.05,     "","p (GeV)","#beta"),n);n++;
    	}    	
    	this.getDataGroup().add(dg,0,st,k,run);      
   	
    } 
    
    public void createECpip(int st) {
    	
    	String tab = "ECpip", tag = null;   	
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	
    	switch (st) {
    	   
        case 0:
        dg = new DataGroup(6,2);
        for(int is=1;is<7;is++){ 
        	tag = is+"_"+st+"_"+k+"_"+run;        	
        	dg.addDataSet(makeH2(tab+"_0_",tag,100,0,4,50,0.9,3,"S"+is,"MM^2 (GeV^2)","W (GeV)"),is-1);
        	dg.addDataSet(makeH1(tab+"_1_",tag,100,0,4," ","MM^2 (GeV^2)"),is-1+6);
           	dg.addDataSet(makeH1(tab+"_2_",tag,100,0,4," ","MM^2 (GeV^2)"),is-1+6);
        }
        break;
        case 1:
        dg = new DataGroup(4,3); int n=0;
    	tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,EB,100,0,40,        "","e- p (GeV)","e- #theta (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,4,100,0,50,         "","#pi+ p (GeV)","#pi+ #theta (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-15,15,    "","e vz (cm)","pip vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-180,180,100,-15,15,  "","#phi (^o)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,   0, 80,100,-15,15,  "","#theta (^o)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-15,15,    "","vz (cm)","#Delta vz (cm)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-180,180,100,-180,180,"","#phi (^o)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,   0, 80,100,-180,180,"","#theta (^o)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,-15,15,100,-180,180,  "","vz (cm)","#Delta#phi (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,EB,100,0.8,1.1,     "","p (GeV)","#beta"),n);n++;
    	}    	
    	this.getDataGroup().add(dg,0,st,k,run);      
   	
    }
    
    public void createECpim(int st) {

    	String tab = "ECpim", tag = null;   	
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);

		F1D     f1 = null; 
		
    	switch (st) {
    	
        case 1:        
		f1 = new F1D("H_pim_EC_resid_f+"+run,"[a]",5,35); 
		f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);		    	
		for(int i=0;i<3;i++) { //pcal,ecin,ecou
	        int inn = 0 ; dg = new DataGroup(6,3); 
			float ylim = (i==0)?5:10;
			for(int n=0; n<3; n++) { //x,y,z
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,5,35,40,-ylim,ylim,"","S"+is+" #theta_pim","DC"+xyz[n]+"-"+det[i]), inn);
					dg.addDataSet(f1, inn); inn++;  
				}
			}	
		 	this.getDataGroup().add(dg,i,st,k,run);
		}
		
		return;
		
    	case 5:
            dg = new DataGroup(3,2); GraphErrors g = null;
        	f1 = new F1D("H_pim_EC_zero_f+"+run,"[a]",0,7); 
        	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
            for (int i=0; i<3; i++) { //pcal,ecin,ecou
    			tag = i+"_"+"1"+"_"+st+"_"+k+"_"+run;
            	g = new GraphErrors(tag); g.setMarkerSize(5); g.setMarkerColor(1); g.setLineColor(1);
            	dg.addDataSet(g, i); 
    			tag = i+"_"+"2"+"_"+st+"_"+k+"_"+run;
            	g = new GraphErrors(tag); g.setMarkerSize(5); g.setMarkerColor(2); g.setLineColor(2);
            	dg.addDataSet(g, i); 
            	dg.addDataSet(f1, i);
            	
            }		
		
    	}
    	this.getDataGroup().add(dg,0,st,k,run);      
        
    }
    
    public void createECpi0(int st) {
    	
    	String tab = "ECpi0", tag = null;   	
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab); 
    	
    	switch (st) {
    	   
        case 1:		
        dg = new DataGroup(4,3); int n=0;
        tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0.5,5,50,0,40,      "pizero","p_mm (GeV)","#theta_mm (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.5,0.5,50,-0.5,0.5,    " ","cx_mm - cx_ecal","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                " ","cx_mm - cx_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                " ","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,15,50,0,8,             " ","Opening Angle (^o)","E1*E2 (GeV^2)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,5.5,50,0,5.5,          " ","p_mm (GeV)","p_ecal (GeV)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,60,0,0.3,                   " ","IVM (GeV)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,60,0,0.3,                   " ","IVM (GeV)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.6,0.0,50,-0.3,0.3,    " ","cx_mm","cy_mm"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.6,0.0,50,-0.3,0.3,    " ","cx_ecal","cy_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,0,5.5,                   " ","p_mm (GeV)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,0,5.5,                   " ","p_mm (GeV)"),n-1);n++;  
    	((H1F)dg.getData(10).get(1)).setFillColor(4);
    	}
    	this.getDataGroup().add(dg,0,st,k,run);      
    	
    } 
    
    public void createECneut(int st) {
    	
    	String tab = "ECneut", tag = null;   	
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
        F1D f1 = new F1D("neut","1/(1+[a]^2/x^2)^0.5", 0.21,2.5); 
        f1.setParameter(0,0.93957); f1.setLineColor(0); f1.setLineStyle(1); f1.setLineWidth(2);  
        
    	switch (st) {
    	   
        case 1:		
        dg = new DataGroup(4,3); int n=0;
        tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,2.5,50,0,40,     "neutron","p_mm (GeV)","#theta_mm (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.5,0.5,50,-0.5,0.5,    " ","cx_mm - cx_ecal","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                "cycut.M2cut ","cx_mm - cx_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                "cxcut.M2cut ","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,2.5,50,0.1,1.1,  "no electron","p_mm (GeV)","#beta_ecal"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,2.5,50,0.1,1.1,            " ","p_mm (GeV)","#beta_ecal"),n);n++;
    	dg.addDataSet(f1, n-1);
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.2,2.0,          "no electron","Mass^2 (GeV^2)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.2,2.0,    "cx,cy cut","Mass^2 (GeV^2)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.6,0.0,50,-0.3,0.3,    " ","cx_mm","cy_mm"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.6,0.0,50,-0.3,0.3,    " ","cx_ecal","cy_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,0,2.5,                   " ","p_mm (GeV)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,0,2.5,                   " ","p_mm (GeV)"),n-1);n++;  
    	((H1F)dg.getData(10).get(1)).setFillColor(4);
    	}
    	this.getDataGroup().add(dg,0,st,k,run);      
    	
    }
   	
    public void createECphot(int st) {
    	
    	String tab = "ECphot", tag = null;   	
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	
    	switch (st) {
    	   
        case 0:    	
        dg = new DataGroup(5,2); int n=0;
        tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,500,60,5,25,   "electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,300,50,0,35,   "electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
       	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,300,50,0,35,"no electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,25,50,0,25,    "electron","#theta_ecal (deg)","#theta_elec (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,25,50,0,25, "no electron","#theta_ecal (deg)","#theta_elec (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,80,-180,180,80,5,22,"electron","#phi_ecal (deg)",  "#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,80,-180,180,80,5,22,"electron","#phi_e (deg)"   ,  "#theta_ecal (deg)"),n);n++;
    	}
    	this.getDataGroup().add(dg,0,st,k,run);  
    	
    }
    
    @Override    
    public void createHistos(int run) {  
	    System.out.println("ECperf:createHistos("+run+")");
    	setRunNumber(run);
    	runlist.add(run);
    	dstinit(run);
    	
    	createECkin(0);
    	createECkin(1);
    	createECkin(2);
    	createECelec(0);
    	createECelec(1);
    	createECelec(2);
    	createECelec(3);
    	createECelec(4);
    	createECelec(5);
    	createSCelec(1);
    	createSCelec(5);
    	createECprot(0);
    	createECprot(1);
    	createECpbar(0);
    	createECpbar(1);
    	createECpim(1);
    	createECpim(5);
    	createECpip(0);
    	createECpip(1);
    	createECneut(1);
    	createECphot(0);
    	createECpi0(1);
    }
    
	public void myinit(){
		goodPROT = false; goodPBAR = false; goodPIP  = false; goodNEUT = false; goodPHOT = false; goodPI0=false;  
		G1_part_ind=G2_part_ind = -1; G1_mom=G2_mom=G1_sec=G2_sec=G1_lay=G2_lay=0;
    }

    public void processEvent(DataEvent event) {
    	
    	ev.setHipoEvent(isHipo3Event);
    	ev.setEventNumber(getEventNumber());
    	ev.setMC(event.hasBank("MC::Event"));
        if(getRunNumber()==5700) ev.setTimeShift(2f);
        
    	if(!ev.procEvent(event)) return;
    	
 	    this.myinit();
	    
	    if(!makeELEC()) return;

	    goodPROT  = makePROT();
	    goodPBAR  = makePBAR();
	    goodPIP   = makePIP();
	    goodPIM   = makePIM();
	    goodNEUT  = makeNEUTRAL();
	    goodPHOT  = makePHOT();
	    goodPHOTR = makePHOTR();
	    
	    fillHists();
    } 	

    // MAKE
    
    public boolean makeELEC(){
    	
        List<Particle> ec = ev.getParticle(11);
        
        if(ec.size()==0 || ec.size()>3) return false;
		         
    	boolean good_fiduc1 = false, good_fiduc2 = false, good_fiduc3 = false; 
        e_ecal_esum = 0f;e_ecal_pcsum=0; e_ecal_ecsum=0;
        elec_ecal_resid.clear();
        elec_ftof_resid.clear();  
        
        Particle epart = ec.get(0);
                
        short   status = (short) epart.getProperty("status");
        boolean   inDC = (status>=2000 && status<3000);
        
        if(!inDC) return false;
        
        e_mom = (float) epart.p();      
        e_vz  = (float) epart.vz();
        
        List<Particle> elecECAL = ev.getECAL((int)epart.getProperty("pindex"));
    	List<Particle> elecFTOF = ev.getFTOF((int)epart.getProperty("pindex"));
    	
    	for (Particle p : elecECAL) {    		
    		e_sect   = (int)   p.getProperty("sector");
    		float en = (float) p.getProperty("energy");   		
    		int  ind = getDet((int) p.getProperty("layer"));
    		      iU = p.hasProperty("iu")?(int)p.getProperty("iu"):0;
    		      iV = p.hasProperty("iv")?(int)p.getProperty("iv"):0;
    		      iW = p.hasProperty("iw")?(int)p.getProperty("iw"):0; 
    		int   iS = (int)p.getProperty("sector");
    		if(ind==0) {
    			lU = p.hasProperty("lu")?(int)p.getProperty("lu"):0;
  		      	lV = p.hasProperty("lv")?(int)p.getProperty("lv"):0;
  		      	lW = p.hasProperty("lw")?(int)p.getProperty("lw"):0; 
    			for (Particle psc : elecFTOF) {
    	    		int scind = (int) psc.getProperty("layer");
    		        Point3D xyz = getResidual(psc);
    	            elec_ftof_resid.add((float)xyz.x(),iS,0,scind-1);
    	    		elec_ftof_resid.add((float)xyz.y(),iS,1,scind-1);
    	    		elec_ftof_resid.add((float)xyz.z(),iS,2,scind-1); 
    			}
    		}
    		Point3D xyz = getResidual(p);	        
    		elec_ecal_resid.add((float)xyz.x(),iS,0,ind);
    		elec_ecal_resid.add((float)xyz.y(),iS,1,ind);
    		elec_ecal_resid.add((float)xyz.z(),iS,2,ind);
    		elec_ecal_resid.add(p.hasProperty("iu")?(float)p.getProperty("iu"):0,iS,3,ind);
    		if(ind>-1) e_ecal_esum  += en;
    		if(ind==0) e_ecal_pcsum  = en;
    		if(ind>0)  e_ecal_ecsum += en;
    		if (ind==0) good_fiduc1 = iU>2&&iV<62&&iW<62;
    	    if (ind==1) good_fiduc2 = iU>2&&iV<36&&iW<36;
    	    if (ind==2) good_fiduc3 = iU>2&&iV<36&&iW<36;   		
   	    }
          
    	if (fiduCuts && !((good_fiduc1)||(good_fiduc1&&good_fiduc2)||(good_fiduc1&&good_fiduc2&&good_fiduc3))) return false;
   	
//    	if (fiduCuts && !(good_fiduc1&&good_fiduc2&&good_fiduc3)) return false;
    	
        if(Math.abs(e_vz+3)<12 && e_mom>0.5){
    		e_sect = (int)   elecECAL.get(0).getProperty("sector");
    		e_x    = (float) elecECAL.get(0).getProperty("x");
    		e_y    = (float) elecECAL.get(0).getProperty("y");
            e_the  = (float) Math.toDegrees(epart.theta());
            e_phi  = (float) Math.toDegrees(epart.phi());
            e_vx   = (float) epart.vx(); 
            e_vy   = (float) epart.vy();
            e_cz   = elecECAL.get(0).hasProperty("cz")?(float) elecECAL.get(0).getProperty("cz"):0;
            Ve     = epart.vector();
            VGS = new LorentzVector(0,0,0,0);                	     
            VGS.add(VB);               	         
            VGS.sub(Ve);                	                   	    
    	    e_Q2 = (float) -VGS.mass2();
            e_xB = e_Q2/(2f*Mp*(Eb-e_mom));
            e_W  = (float) Math.sqrt(Mp*Mp + e_Q2*(1f/e_xB-1f) );
            return true;
         }  
        
        return false;

    }
    
    public boolean makePIM() {
    	
        List<Particle> nlist = ev.getParticle(212);
        if(nlist.size()==0) return false;
        
        pim_ecal.clear();
        
        for (Particle p : nlist) {
        	short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);
            if(inDC && p.p()>0.5) pim_ecal.add(p);
        }        
        return pim_ecal.size()>0;
    }

    
    public boolean makePIPP() {
    	
        List<Particle> nlist = ev.getParticle(211);       
        if (nlist.size()==0) return false;
        
        Particle pipart = nlist.get(0);

        pip_mom  = (float) pipart.p();
        short status = (short) pipart.getProperty("status");
        boolean inDC = (status>=2000 && status<3000);
        
        pip_ecal_esum = 0f;
    	for (Particle p : nlist) pip_ecal_esum += p.getProperty("energy");        	
        if(inDC && pip_mom>0.5){ 
            pip_the  = (float) Math.toDegrees(pipart.theta());
            pip_phi  = (float) Math.toDegrees(pipart.phi());
            pip_vz   = (float) pipart.vz();
        	pip_beta = (float) pipart.getProperty("beta");
        	Vpip     =         pipart.vector();
            return true;       	
        }
        return false;
    } 
    
    public boolean makePIP() {
    	
        List<Particle> nlist = ev.getParticle(211);
        if(nlist.size()==0) return false;
        
        pip_ecal.clear();   
        
        for (Particle p : nlist) {
            short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);
        	if(inDC && p.p()>0.5) pip_ecal.add(p); 
        }        
        return pip_ecal.size()>0;
    }
    
    public boolean makePROT() {
    	
        List<Particle> nlist = ev.getParticle(2212);
        if(nlist.size()==0) return false;
        
        prot_ecal.clear();        
        
        for (Particle p : nlist) {
            short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);
        	if(inDC && p.p()>0.5) prot_ecal.add(p);        	
        }        
        return prot_ecal.size()>0;
    }
    
    public boolean makePBAR() {
    	
        List<Particle> nlist = ev.getParticle(2213);
        if(nlist.size()==0) return false;
        
        pbar_ecal.clear();        
        
        for (Particle p : nlist) {
            short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);
        	if(inDC && p.p()>0.5) pbar_ecal.add(p);        	
        }        
        return pbar_ecal.size()>0;
    }
    
    public boolean makePHOT() {
    	
        List<Particle> nlist = ev.getParticle(22);
        if(nlist.size()==0) return false;
        
        phot_ecal.clear();
        
        for (Particle p : nlist) {
            short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);
            if(inDC && p.e()>0.05) phot_ecal.add(p); 

        }        
        return phot_ecal.size()>0;
    } 
    
    public boolean makeNEUT() {
    	
        List<Particle> nlist = ev.getParticle(2112);
        if(nlist.size()==0) return false;
        
        neut_ecal.clear();
        
        for (Particle p : nlist) {
            short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);            
        	if(inDC) neut_ecal.add(p);
        }        
        return neut_ecal.size()>0;
    } 
    
    public boolean makeNEUTRAL() {
    	
    	List<Particle> nlist = null;
    	
        neut_ecal.clear();
        
        nlist = ev.getParticle(22);
        
        for (Particle p : nlist) {
            short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);
            if(inDC && p.e()>0.05) neut_ecal.add(p); 

        } 
        
        nlist = ev.getParticle(2112);
        
        for (Particle p : nlist) {
            short status = (short) p.getProperty("status");
            boolean inDC = (status>=2000 && status<3000);            
        	if(inDC) neut_ecal.add(p);
        }   
        
        return neut_ecal.size()>0;    	
    }
    
    public boolean makeNEUTT() {
    	
        List<Particle> nlist = ev.getParticle(2112);
        if (nlist.size()==0) return false;
        
        for (int is=0; is<6; is++) ecal_neut_esum[is]=0;
        
		List<Particle> neutECAL = ev.getECAL((int)nlist.get(0).getProperty("pindex"));
        
        ecal_neut_sec  = (int)   neutECAL.get(0).getProperty("sector");
        ecal_neut_the  = (float) Math.toDegrees(nlist.get(0).theta());
        ecal_neut_phi  = (float) Math.toDegrees(nlist.get(0).phi());
        ecal_neut_cx   = (float) (nlist.get(0).px()/nlist.get(0).p());
        ecal_neut_cy   = (float) (nlist.get(0).py()/nlist.get(0).p());
        ecal_neut_beta = (float)  nlist.get(0).getProperty("beta");
        for (Particle p : neutECAL) ecal_neut_esum[0] += p.getProperty("energy"); 
        	
        return true;

    }
/*    
    public boolean makePI0part(DataEvent event) {
    	int n = 0;
    	
        ecClusters = part.readEC(event,"REC::Calorimeter");  
        if (ecClusters.size()==0) return false;
        
        part.getNeutralResponses(ecClusters);
        for (int is=1; is<7; is++) {
            ecal_pi0_mass = (float)Math.sqrt(part.getTwoPhotonInvMass(is));
            if(ecal_pi0_mass>0&&part.iis[0]==part.iis[1]) {
            	VG1 = part.VG1;   	
            	VG2 = part.VG2;   	
            	VPI0 = new LorentzVector(0,0,0,0);
            	VPI0.add(VG1);
				VPI0.add(VG2);
				ecal_pi0_sec  = is;
				ecal_pi0_mass = (float)VPI0.mass();
				ecal_pi0_e    = (float)VPI0.e();
				ecal_pi0_mom  = (float)VPI0.p();
				ecal_pi0_the  = (float)Math.toDegrees(VPI0.theta());
				ecal_pi0_phi  = (float)Math.toDegrees(VPI0.phi());
				ecal_pi0_opa = (float)Vangle(VG1.vect(),VG2.vect());
				n++;	
            }
        }
        
        return n==1&&ecal_pi0_mass>0.05;
	
    }
*/


    public boolean makeOldPI0() {
    	
    	if(!goodPROT) return false;
    	if(!goodPHOT) return false;
        int n = 0;        
        for (Particle p : phot_ecal) {
	        List<Particle> photECAL = ev.getECAL((int)p.getProperty("pindex"));
//			System.out.println("PHOT: "+photECAL.get(0).getProperty("x")+" "+photECAL.get(0).getProperty("y"));
			if( p.p()>0.2 && Math.toDegrees(p.theta())>6 && G1_mom < p.p()){
				G1_part_ind = n;
				G1_sec = (int) photECAL.get(0).getProperty("sector");
				G1_lay = (int) photECAL.get(0).getProperty("layer");
				G1_mom = (float) p.p();
				G1_e   = (float) p.e();
				G1_the = (float) Math.toDegrees(p.theta());
				G1_phi = (float) Math.toDegrees(p.phi());
//				System.out.println("PHT1: "+G1_sec+" "+G1_lay+" "+G1_the+" "+G1_phi);
//				VG1 = new LorentzVector(p.px(),p.py(),p.pz(),p.p());
			}
			if( G1_part_ind>-1 && n!=G1_part_ind && p.p()>0.2 && Math.toDegrees(p.theta())>6 && G2_mom < p.p() && p.p() < G1_mom){
				G2_part_ind = n;
				G2_sec = (int) photECAL.get(0).getProperty("sector");
				G2_lay = (int) photECAL.get(0).getProperty("layer");
				G2_mom = (float)p.p();
				G2_e   = (float)p.e();
				G2_the = (float) Math.toDegrees(p.theta());
				G2_phi = (float) Math.toDegrees(p.phi());
//				System.out.println("PHT2: "+G2_sec+" "+G2_lay+" "+G2_the+" "+G2_phi);
//				VG2 = new LorentzVector(p.px(),p.py(),p.pz(),p.p());
			}  
			n++;
        }
        if(G1_part_ind>-1 && G2_part_ind>-1){
//				VPI0 = new LorentzVector(0,0,0,0);
//				VPI0.add(VG1);
//				VPI0.add(VG2);
				ecal_pi0_sec  = G1_sec;
//				ecal_pi0_mass = (float)VPI0.mass();
//				ecal_pi0_e    = (float)VPI0.e();
//				ecal_pi0_mom  = (float)VPI0.p();
//				ecal_pi0_the  = (float)Math.toDegrees(VPI0.theta());
//				ecal_pi0_phi  = (float)Math.toDegrees(VPI0.phi());
//				System.out.println("PIZ0: "+G1_sec+" "+G2_sec+" "+ecal_pi0_the+" "+ecal_pi0_phi);
//				ecal_pi0_opa  = (float)Vangle(VG1.vect(),VG2.vect());
				ecal_pi0_X    = (float)((G1_e-G2_e)/(G1_e+G2_e));
//		        System.out.println(G1_part_ind+" "+G2_part_ind+" "+ecal_pi0_mass+" "+ecal_pi0_X);
				return ecal_pi0_mass>0.05;
//				return (ecal_pi0_opa>3 && ecal_pi0_opa>9 * (1 - ecal_pi0_e/4) && 
//						ecal_pi0_the>8 && ecal_pi0_mass>0.05 && ecal_pi0_mass<0.5); 		
		}
		
		return false;
        
    }


    public boolean makePI0() {
    	
    	if(!goodPHOT) return false;
    	pi0_ecal.clear();
        int n = 0;        
        for (Particle p : phot_ecal) {
			if( p.p()>0.2 && Math.toDegrees(p.theta())>6 && G1_mom < p.p()){
				G1_mom = (float) p.p();
				G1_part_ind = n;
			}
			if( G1_part_ind>-1 && n!=G1_part_ind && p.p()>0.2 && Math.toDegrees(p.theta())>6 && G2_mom < p.p() && p.p() < G1_mom){
				G2_mom = (float) p.p();
				G2_part_ind = n;
			}
			n++;			
        }
        if(!(G1_part_ind>-1 && G2_part_ind>-1)) return false;
        
		NeutralMeson nm = new NeutralMeson();
		nm.addPhoton(phot_ecal.get(G1_part_ind));			
		nm.addPhoton(phot_ecal.get(G2_part_ind));			
		pi0_ecal.add(nm);
		
		return true;
       
    }

    /*
    public boolean makePI0() {
    	
    	if(!findPhotons()) return false; 
    	
    	pi0_ecal.clear();
  	
    	for (int i=0; i<phot_ecal.size()-1; i++) {
    		for (int j=i+1; j<phot_ecal.size(); j++) {
        		NeutralMeson nm = new NeutralMeson();
        		nm.addPhoton(phot_ecal.get(i));			
       		    nm.addPhoton(phot_ecal.get(j));			
				pi0_ecal.add(nm);
    		}    		
    	}
    	return pi0_ecal.size()>0;
    }
*/    
    public class NeutralMeson {
    	
    	public float ecal_pi0_mass;
    	public float ecal_pi0_e;
    	public float ecal_pi0_mom;
    	public float ecal_pi0_the;
    	public float ecal_pi0_phi;
    	public float ecal_pi0_opa;
    	public float ecal_pi0_X;
    	
    	public  LorentzVector VPI0;
    	public  LorentzVector VG1;
    	public  LorentzVector VG2;
    	
    	public List<Particle>  plist = new ArrayList<Particle>();
    	
    	public void addPhoton(Particle p) {
    		plist.add(p);
    	}
    	
    	public Particle getPhoton(int n) {
    		return plist.get(n);    		
    	}
    	
    	public int getPhotonSector(int n) {
    		return (int) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("sector");
    	}
    	
    	public int getPhotonLayer(int n) {
    		return (int) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("layer");
    	}
    	
    	public void filter() {
    		
    	}
    	
    	public boolean getMeson() {
    		Particle p1 = plist.get(0);
    		Particle p2 = plist.get(1);
			VG1 = new LorentzVector(p1.px(),p1.py(),p1.pz(),p1.p());
			VG2 = new LorentzVector(p2.px(),p2.py(),p2.pz(),p2.p());
			VPI0 = new LorentzVector(0,0,0,0);
			VPI0.add(VG1);
			VPI0.add(VG2);
			return true;
    	}
    	
    	public boolean getMesonKin() {
    		if(!getMeson()) return false;
			this.ecal_pi0_mass = (float)VPI0.mass();
			this.ecal_pi0_e    = (float)VPI0.e();
			this.ecal_pi0_mom  = (float)VPI0.p();
			this.ecal_pi0_the  = (float)Math.toDegrees(VPI0.theta());
			this.ecal_pi0_phi  = (float)Math.toDegrees(VPI0.phi());
			this.ecal_pi0_opa  = (float)Vangle(VG1.vect(),VG2.vect());
			this.ecal_pi0_X    = (float)((VG1.e()-VG2.e())/(VG1.e()+VG2.e()));
			return this.ecal_pi0_mass > 0.08 && getPhotonSector(0)==getPhotonSector(1);			
    	}
    	    	    	
    }
    
    public boolean makePHOTT() {
    	
        List<Particle> nlist = ev.getParticle(22);
        if (nlist.size()==0) return false;
                
        for (int is=0; is<6; is++) ecal_phot_esum[is]=0;
        if(nlist.size()==1) {
        	ecal_phot_sec  = (int)   nlist.get(0).getProperty("sector");
        	ecal_phot_the  = (float) Math.toDegrees(nlist.get(0).theta());
        	ecal_phot_phi  = (float) Math.toDegrees(nlist.get(0).phi());
        	ecal_phot_beta = (float) nlist.get(0).getProperty("beta");
        	for (Particle p : nlist) ecal_phot_esum[ecal_phot_sec-1] += p.getProperty("energy");        	    
        	return true;
        }
    	return false;
    } 
    
    public boolean makePHOTR() {
    	
        List<Particle> nlist = ev.getParticle(22);
        if (nlist.size()==0) return false;
        ecal_rad.clear();
        for (Particle p : nlist) {
			float the = (float) Math.toDegrees(p.theta());
			float phi = (float) Math.toDegrees(p.phi());
   		    List<Particle> ecECAL = ev.getECAL((int)p.getProperty("pindex"));
        	for (Particle pec : ecECAL) {
            	int iS = (int) pec.getProperty("sector");
        		if(iS==e_sect) {
        			int   ind = getDet((int) pec.getProperty("layer"));
        			float nrg = (float) pec.getProperty("energy");
            		int    iu =  (int)  pec.getProperty("iu");
        			ecal_rad.add(e_the-the,iS,0,ind);
        			ecal_rad.add(e_phi-phi,iS,1,ind);
        			ecal_rad.add(nrg,iS,2,ind);
        			ecal_rad.add(the,iS,3,ind);
        			ecal_rad.add((float)(iu-iU),iS,4,ind);
        		}
        	}
        }
        return true;
    }

	// SELECT	
		
	public boolean selectElastic(){
		boolean res = false;
		if(prot_part_ind>-1){
			elast_dPhi = prot_phi - e_phi + 180f;
			while(elast_dPhi> 180f)elast_dPhi -= 360f;
			while(elast_dPhi<-180f)elast_dPhi += 360f;
			float tantan = (float) (Math.tan(Ve.theta()/2) * Math.tan(Vprot.theta() ) );
			elast_EB = 0.93827f * (1-tantan)/tantan;

			if(EB<4f &&  prot_the > 60f-30f*prot_mom)res = true;
			if(EB>4f && EB<9f && Math.abs(elast_dPhi)<10f && prot_the > 55f && prot_the > 70f-30f*prot_mom && Math.abs(elast_EB-6.5f)<1.2f )res = true;
			if(EB>9f && Math.abs(elast_dPhi)<10f && prot_the > 70f-30f*prot_mom && Math.abs(elast_EB-10f)<1.2f )res = true;
		}
		return res;
	}
	
	public boolean select_ep(){
		int run = getRunNumber();
		DataGroup ECprot = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECprot"),run);
		if(prot_ecal.size()>0) {			
			prot_mom  = (float) prot_ecal.get(0).p();
            prot_the  = (float) Math.toDegrees(prot_ecal.get(0).theta());
            prot_phi  = (float) Math.toDegrees(prot_ecal.get(0).phi());
            prot_vz   = (float) prot_ecal.get(0).vz();
        	prot_beta = (float) prot_ecal.get(0).getProperty("beta");
         	Vprot     =         prot_ecal.get(0).vector();	
			ep_dPhi   = prot_phi - e_phi + 180f;
			while(ep_dPhi> 180f)ep_dPhi -= 360f;
			while(ep_dPhi<-180f)ep_dPhi += 360f;
			LorentzVector VmissN = new LorentzVector(0,0,0,0);
			VmissN.add(VT);
			VmissN.add(VB);
			VmissN.sub(Ve);
			VmissN.sub(Vprot);
			ep_MM = (float)VmissN.mass2();
			((H2F) ECprot.getData(e_sect-1).get(0)).fill(ep_MM,e_W);
			((H1F) ECprot.getData(e_sect-1+6).get(0)).fill(ep_MM);  	
			pi0_mom=-1f;pi0_the=-1f;pi0_phi=-1f;
			if(ep_MM<0.1) {
//		        List<Particle> protECAL = ev.getECAL((int)prot_ecal.get(0).getProperty("pindex"));
//				System.out.println("PROT: "+protECAL.get(0).getProperty("x")+" "+protECAL.get(0).getProperty("y"));
				pi0_mom=(float)VmissN.p();
				pi0_the=(float)Math.toDegrees(VmissN.theta());
				pi0_phi=(float)Math.toDegrees(VmissN.phi());
				pi0_cx =(float)VmissN.px()/pi0_mom;
				pi0_cy =(float)VmissN.py()/pi0_mom;
				return true;
			}
		}		
		return false;
	}	
	
	public boolean select_epbar(){
		int run = getRunNumber();
		DataGroup ECpbar = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECpbar"),run);
		if(prot_ecal.size()>0 && pbar_ecal.size()>0) {			
			pbar_mom  = (float) pbar_ecal.get(0).p();
            pbar_the  = (float) Math.toDegrees(pbar_ecal.get(0).theta());
            pbar_phi  = (float) Math.toDegrees(pbar_ecal.get(0).phi());
            pbar_vz   = (float) pbar_ecal.get(0).vz();
        	pbar_beta = (float) pbar_ecal.get(0).getProperty("beta");
         	Vpbar     =         pbar_ecal.get(0).vector();	
			epbar_dPhi = pbar_phi - e_phi + 180f;
			while(epbar_dPhi> 180f)epbar_dPhi -= 360f;
			while(epbar_dPhi<-180f)epbar_dPhi += 360f;
			LorentzVector VmissN = new LorentzVector(0,0,0,0);
			VmissN.add(VT);
			VmissN.add(VB);
			VmissN.sub(Ve);
			VmissN.sub(Vprot);
			VmissN.sub(Vpbar);
			epbar_MM = (float)VmissN.mass2();
			((H2F) ECpbar.getData(e_sect-1).get(0)).fill(epbar_MM,e_W);
			((H1F) ECpbar.getData(e_sect-1+6).get(0)).fill(epbar_MM);  	
			if(epbar_MM<4.2) {
				return true;
			}
		}		
		return false;
	}
	
	public boolean select_epip(){
		int run = getRunNumber();
		DataGroup ECpip = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECpip"),run);
		if(pip_ecal.size()>0) {			
			pip_mom  = (float) pip_ecal.get(0).p();
            pip_the  = (float) Math.toDegrees(pip_ecal.get(0).theta());
            pip_phi  = (float) Math.toDegrees(pip_ecal.get(0).phi());
            pip_vz   = (float) pip_ecal.get(0).vz();
        	pip_beta = (float) pip_ecal.get(0).getProperty("beta");
         	Vpip     =         pip_ecal.get(0).vector();	
			epip_dPhi = pip_phi - e_phi + 180f;
			while(epip_dPhi> 180f)epip_dPhi -= 360f;
			while(epip_dPhi<-180f)epip_dPhi += 360f;
			LorentzVector VmissN = new LorentzVector(0,0,0,0);
			VmissN.add(VT);
			VmissN.add(VB);
			VmissN.sub(Ve);
			VmissN.sub(Vpip);
			epip_MM = (float)VmissN.mass2();
			((H2F) ECpip.getData(e_sect-1).get(0)).fill(epip_MM,e_W);    
			((H1F) ECpip.getData(e_sect-1+6).get(0)).fill(epip_MM);
			neut_mom=-1f;neut_the=-1f;neut_phi=-1f;
			if(epip_MM<1.2) {
				neut_mom=(float)VmissN.p();
				neut_the=(float)Math.toDegrees(VmissN.theta());
				neut_phi=(float)Math.toDegrees(VmissN.phi());
				neut_cx =(float)VmissN.px()/neut_mom;
				neut_cy =(float)VmissN.py()/neut_mom;
				return true;
			}
		}
		return false;
	}
	
	public void debug() {
	    IndexGenerator ig = new IndexGenerator();
	    for (Map.Entry<Long,List<Particle>>  entry : ev.part.getMap().entrySet()){
	           long hash = entry.getKey();
	           int pid = ig.getIndex(hash, 0);
	           int sec = ig.getIndex(hash, 1);
	           for (Particle pp : ev.part.getItem(pid,sec)) {	        	   
	               System.out.println(pid+" "+sec+" "+(int)pp.getProperty("layer")
	                                             +" "+(int)pp.getProperty("status")
	                                             +" "+     pp.getProperty("energy"));
	           }
	    }		
	} 

	// FILL
	
	public void fillHists(){
				
		fillECkin();
		fillECelec();
		fillSCelec();
		
		if(goodPIM) fillECpim();	
		
	    if(goodPIP) { // FC pi+
	    	if(select_epip()) { //tagged neutron (e' pi+)
	    		fillECpip();
	    		fillECneut();
	    		fillECphot();	    		
	    	}	
	    }

	    if(goodPROT) { // FC proton
	    	if(select_ep()) { // tagged pizero (e' p)
	    		fillECprot();
	    		fillECpi0();
//	    		System.out.println(" ");
	    	}
	    	if(select_epbar()) {
	    		fillECpbar();
	    	}
	    }
	}
	
	public void fillECkin() {
		int run = getRunNumber();
		DataGroup ECkin = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECkin"),run);				
		((H2F) ECkin.getData(e_sect-1).get(0)).fill(e_the,e_mom);
		((H2F) ECkin.getData(e_sect-1+6).get(0)).fill(e_W,e_mom);
		((H2F) ECkin.getData(e_sect-1+12).get(0)).fill(e_W,e_the);
		if (e_the>6) ((H2F) ECkin.getData(e_sect-1+18).get(0)).fill(e_W,(e_phi>-30?e_phi:360+e_phi)-(e_sect-1)*60);
		ECkin = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECkin"),run);				
		((H1F) ECkin.getData(e_sect-1).get(0)).fill(lU);
		((H1F) ECkin.getData(e_sect-1+6).get(0)).fill(lV);
		((H1F) ECkin.getData(e_sect-1+12).get(0)).fill(lW);
		ECkin = this.getDataGroup().getItem(0,2,getDetectorTabNames().indexOf("ECkin"),run);				
		((H2F) ECkin.getData(0).get(0)).fill(-e_x,e_y);
	}
	
	public void fillECelec() {
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("ECelec");
		IndexGenerator ig = new IndexGenerator();
		
		DataGroup dg0 = this.getDataGroup().getItem(0,0,k,run);				
//		DataGroup dg1 = this.getDataGroup().getItem(0,1,k,run);				
		DataGroup dg2 = this.getDataGroup().getItem(0,2,k,run);				
		DataGroup dg3 = this.getDataGroup().getItem(0,3,k,run);				
		DataGroup dg4 = this.getDataGroup().getItem(0,4,k,run);				
		
		((H2F) dg0.getData(e_sect-1).get(0)).fill(e_mom,e_ecal_esum/1000f/e_mom);
		((H2F) dg0.getData(e_sect-1+12).get(0)).fill(e_ecal_pcsum/1000f,e_ecal_ecsum/1000f);

		for (Map.Entry<Long,Float>  entry : elec_ecal_resid.getMap().entrySet()){
			long hash = entry.getKey();
			int is = ig.getIndex(hash, 0); int ic = ig.getIndex(hash, 1); int il = ig.getIndex(hash, 2);
			DataGroup dg1 = this.getDataGroup().getItem(il,1,k,run);				
			if(ic==3) ((H2F) dg0.getData(e_sect-1+6).get(0)).fill(elec_ecal_resid.getItem(e_sect,3,0),e_ecal_esum/1000f/e_mom);
			if(ic<3) {
//			((H2F)dg1.getData(is-1+ic*6+il*12).get(0)).fill(e_the,entry.getValue());
			((H2F)dg1.getData(is-1+ic*6).get(0)).fill(e_mom,entry.getValue());
//			((H2F)dg1.getData(is-1+ic*6).get(0)).fill(e_cz,entry.getValue());
			}
		}
		counter[e_sect-1][0]++;
		
		for (Map.Entry<Long,Float>  entry : ecal_rad.getMap().entrySet()){
			long hash = entry.getKey();
			int is = ig.getIndex(hash, 0); int ic = ig.getIndex(hash, 1); int il = ig.getIndex(hash, 2);
			if(ic<2) {
		        if(il==0) {
		        	if(ic==0&&ecal_rad.getItem(is,0,il)>-1&&ecal_rad.getItem(is,0,il)<1)   counter[e_sect-1][1]++;
		        	if(ic==1&&ecal_rad.getItem(is,1,il)>-2&&ecal_rad.getItem(is,1,il)<1.5) counter[e_sect-1][2]++;
		        }
				int sgn = getTorusCurrent(run)<0?1:-1;
		        if(ic==1) sgn=1;
				((H2F)dg2.getData(is-1+ic*6+il*12).get(0)).fill(e_the,sgn*entry.getValue()); 
				((H2F)dg4.getData(is-1+ic*6+il*12).get(0)).fill(e_mom,sgn*entry.getValue());
				((H2F)dg3.getData(is-1+ic*6+il*12).get(0)).fill(ecal_rad.getItem(is,2,il)/1e3,sgn*entry.getValue());
			}
		}
		
	}
	
	public void fillSCelec() {
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("SCelec");
		 
		IndexGenerator ig = new IndexGenerator();
		
		for (Map.Entry<Long,Float>  entry : elec_ftof_resid.getMap().entrySet()){
			long hash = entry.getKey();
			int is = ig.getIndex(hash, 0); int ic = ig.getIndex(hash, 1); int il = ig.getIndex(hash, 2);
			DataGroup dg1 = this.getDataGroup().getItem(il,1,k,run);
			((H2F)dg1.getData(is-1+ic*6).get(0)).fill(e_the,entry.getValue());
		}		
	}
    
    public void fillECprot() {
		int run = getRunNumber();   
		DataGroup dg1 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECprot"),run);
		((H2F)dg1.getData(0).get(0)).fill(e_mom,e_the);
		((H2F)dg1.getData(1).get(0)).fill(prot_mom,prot_the);
		((H2F)dg1.getData(2).get(0)).fill(e_vz,prot_vz);
		((H2F)dg1.getData(3).get(0)).fill(prot_phi,prot_vz-e_vz);
		((H2F)dg1.getData(4).get(0)).fill(prot_the,prot_vz-e_vz);
		((H2F)dg1.getData(5).get(0)).fill(prot_vz,prot_vz-e_vz);
		((H2F)dg1.getData(6).get(0)).fill(prot_phi,ep_dPhi);
		((H2F)dg1.getData(7).get(0)).fill(prot_the,ep_dPhi);
		((H2F)dg1.getData(8).get(0)).fill(prot_vz,ep_dPhi);
		((H2F)dg1.getData(9).get(0)).fill(prot_mom,prot_beta);  	
    } 
    
    public void fillECpbar() {
		int run = getRunNumber();   
		DataGroup dg1 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECpbar"),run);
		((H2F)dg1.getData(0).get(0)).fill(e_mom,e_the);
		((H2F)dg1.getData(1).get(0)).fill(pbar_mom,pbar_the);
		((H2F)dg1.getData(2).get(0)).fill(e_vz,pbar_vz);
		((H2F)dg1.getData(3).get(0)).fill(pbar_phi,pbar_vz-e_vz);
		((H2F)dg1.getData(4).get(0)).fill(pbar_the,pbar_vz-e_vz);
		((H2F)dg1.getData(5).get(0)).fill(pbar_vz, pbar_vz-e_vz);
		((H2F)dg1.getData(6).get(0)).fill(pbar_phi,epbar_dPhi);
		((H2F)dg1.getData(7).get(0)).fill(pbar_the,epbar_dPhi);
		((H2F)dg1.getData(8).get(0)).fill(pbar_vz, epbar_dPhi);
		((H2F)dg1.getData(9).get(0)).fill(pbar_mom,pbar_beta);  	
    } 
    
    public void fillECpip() {
		int run = getRunNumber();
		DataGroup dg1 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECpip"),run);
		((H2F)dg1.getData(0).get(0)).fill(e_mom,e_the);
		((H2F)dg1.getData(1).get(0)).fill(pim_mom,pip_the);
		((H2F)dg1.getData(2).get(0)).fill(e_vz,pip_vz);
		((H2F)dg1.getData(3).get(0)).fill(pip_phi,pip_vz-e_vz);
		((H2F)dg1.getData(4).get(0)).fill(pip_the,pip_vz-e_vz);
		((H2F)dg1.getData(5).get(0)).fill(pip_vz,pip_vz-e_vz);
		((H2F)dg1.getData(6).get(0)).fill(pip_phi,epip_dPhi);
		((H2F)dg1.getData(7).get(0)).fill(pip_the,epip_dPhi);
		((H2F)dg1.getData(8).get(0)).fill(pip_vz,epip_dPhi);
		((H2F)dg1.getData(9).get(0)).fill(pip_mom,pip_beta);  	
    }
	
	public void fillECpim() {		
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("ECpim");
		for (Particle p : pim_ecal) {
       	    int    is = (int) p.getProperty("sector"); 
    		int    il = getDet((int) p.getProperty("layer"));
    		float the = (float) Math.toDegrees(p.theta());
    		float phi = (float) Math.toDegrees(p.phi());
    		float nrg = (float) p.getProperty("energy");	
	        Point3D xyz = getResidual(p);
			DataGroup dg1 = this.getDataGroup().getItem(il,1,k,run);	
	        ((H2F)dg1.getData(is-1+   0).get(0)).fill(the,xyz.x());
	        ((H2F)dg1.getData(is-1+   6).get(0)).fill(the,xyz.y());	        
	        ((H2F)dg1.getData(is-1+  12).get(0)).fill(the,xyz.z());	        
		}
			
	}
	
    public void fillECpi0() {	 	    
		int run = getRunNumber();
		DataGroup ECpi0  = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECpi0"),run);
		DataGroup ECprot = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECprot"),run);
		
		boolean good_tagged_fiduc = false;
        ((H2F)ECpi0.getData(0).get(0)).fill(pi0_mom, pi0_the);
        
	    float cxmm = (float) (Math.sin(pi0_the*3.14159f/180f)*Math.cos(pi0_phi*3.141259f/180f));
	    float cymm = (float) (Math.sin(pi0_the*3.14159f/180f)*Math.sin(pi0_phi*3.141259f/180f));  
	    
        if(pi0_mom>0.4) {
        	float nphi = newPhi(pi0_phi);
        	float   cx = (float) (Math.sin(pi0_the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
        	float   cy = (float) (Math.sin(pi0_the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
        	if(neutFiduc(2,cx,cy)) {
    			((H1F)ECprot.getData(e_sect-1+6).get(1)).fill(ep_MM);  	
                ((H2F)ECpi0.getData(8).get(0)).fill(-cx,cy);
                ((H1F)ECpi0.getData(10).get(0)).fill(pi0_mom);
        		good_tagged_fiduc = true;
        	}
        }
        
        if(!makePI0()) return;
        
        for (NeutralMeson nm: pi0_ecal) {
        	
        	if (nm.getMesonKin()) {
        	
        	float       cx = (float) (Math.sin(nm.ecal_pi0_the*3.14159f/180f)*Math.cos(nm.ecal_pi0_phi*3.141259f/180f));
        	float       cy = (float) (Math.sin(nm.ecal_pi0_the*3.14159f/180f)*Math.sin(nm.ecal_pi0_phi*3.141259f/180f));   
        
        	float      dcx = cxmm-cx;
        	float      dcy = cymm-cy; 
        
        	boolean  cxcut = Math.abs(dcx)<0.06;	
        	boolean  cycut = Math.abs(dcy)<0.06;
        
        	if(nm.getPhotonSector(0)!=e_sect) {         		
            	((H2F)ECpi0.getData(1).get(0)).fill(dcx,dcy);
            	if(cxcut) ((H1F)ECpi0.getData(3).get(0)).fill(dcy);
            	if(cycut) ((H1F)ECpi0.getData(2).get(0)).fill(dcx);                
            	((H1F)ECpi0.getData(6).get(0)).fill(nm.ecal_pi0_mass);            	
        	}
        
        	if (cxcut && cycut && good_tagged_fiduc) {        	
        		((H2F)ECpi0.getData(5).get(0)).fill(pi0_mom, nm.ecal_pi0_mom);	    	
        		((H1F)ECpi0.getData(7).get(0)).fill(nm.ecal_pi0_mass);
        		if(nm.ecal_pi0_mom>0.4) {
        			((H2F)ECpi0.getData(4).get(0)).fill(nm.ecal_pi0_opa,nm.VG1.e()*nm.VG2.e());
        			float nphi = newPhi(nm.ecal_pi0_phi);
        			cx = (float) (Math.sin(nm.ecal_pi0_the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
        			cy = (float) (Math.sin(nm.ecal_pi0_the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
        			((H2F)ECpi0.getData(9).get(0)).fill(-cx,cy);
        			((H1F)ECpi0.getData(10).get(1)).fill(pi0_mom);
        		}
        	} 
        	}
        }
    }
    
    public void fillECneut() {
		int run = getRunNumber();
		DataGroup ECpip  = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECpip"),run);
		DataGroup ECneut = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECneut"),run);
		DataGroup ECphot = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECphot"),run);
		
        boolean good_tagged_fiduc = false;
        ((H2F)ECneut.getData(0).get(0)).fill(neut_mom, neut_the);
        
	    float cxmm = (float) (Math.sin(neut_the*3.14159f/180f)*Math.cos(neut_phi*3.141259f/180f));
	    float cymm = (float) (Math.sin(neut_the*3.14159f/180f)*Math.sin(neut_phi*3.141259f/180f));
       
        if(neut_mom>0.4) {
        	float nphi = newPhi(neut_phi);
        	float   cx = (float) (Math.sin(neut_the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
        	float   cy = (float) (Math.sin(neut_the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
        	if(neutFiduc(1,cx,cy)) {
    			((H1F)ECpip.getData(e_sect-1+6).get(1)).fill(epip_MM);
                ((H2F)ECneut.getData(8).get(0)).fill(-cx,cy);
                ((H1F)ECneut.getData(10).get(0)).fill(neut_mom);
        		good_tagged_fiduc = true;
        	}
        }        
		
        if(!goodNEUT) return;

        for (Particle p: neut_ecal) {
        	
    		List<Particle> neutECAL = ev.getECAL((int)p.getProperty("pindex"));            
            ecal_neut_sec  = (int)   neutECAL.get(0).getProperty("sector");
            ecal_neut_the  = (float) Math.toDegrees(p.theta());
            ecal_neut_phi  = (float) Math.toDegrees(p.phi());
            ecal_neut_cx   = (float) (p.px()/p.p());
            ecal_neut_cy   = (float) (p.py()/p.p());
            ecal_neut_beta = (float) p.getProperty("beta");
            
            ecal_neut_esum[0] = 0;
            for (Particle pp : neutECAL) ecal_neut_esum[0] += pp.getProperty("energy"); 
        
            float    mass2 = neut_mom*neut_mom*(1f/(ecal_neut_beta*ecal_neut_beta)-1);
            boolean mass2cut = neut_mom<1.2?mass2>0.45:true;
            
            float       cx = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.cos(ecal_neut_phi*3.141259f/180f));
            float       cy = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.sin(ecal_neut_phi*3.141259f/180f));
   		
            float      dcx = cxmm-cx;
            float      dcy = cymm-cy;
   		
            boolean  cxcut = Math.abs(dcx)<0.1;	
            boolean  cycut = Math.abs(dcy)<0.1;
        
            if(ecal_neut_sec!=e_sect) {        	
            	((H2F)ECneut.getData(1).get(0)).fill(dcx,dcy);
            	if(cxcut&&mass2cut) ((H1F)ECneut.getData(3).get(0)).fill(dcy);
            	if(cycut&&mass2cut) ((H1F)ECneut.getData(2).get(0)).fill(dcx); 
        	
             	((H2F)ECphot.getData(2).get(0)).fill(ecal_neut_esum[0],ecal_neut_the);        
             	((H2F)ECphot.getData(4).get(0)).fill(ecal_neut_the,e_the);        
             	((H2F)ECneut.getData(4).get(0)).fill(neut_mom, ecal_neut_beta);        
            	((H1F)ECneut.getData(6).get(0)).fill(mass2);
            }
		
            if(ecal_neut_sec==e_sect) {
             	((H2F)ECphot.getData(1).get(0)).fill(ecal_neut_esum[0],ecal_neut_the);        
             	((H2F)ECphot.getData(3).get(0)).fill(ecal_neut_the,e_the);        
             	((H2F)ECphot.getData(5).get(0)).fill(ecal_neut_phi,ecal_neut_the);        
             	((H2F)ECphot.getData(6).get(0)).fill(e_phi,ecal_neut_the);        
            }
        
            if (cxcut && cycut && good_tagged_fiduc) {
            	((H2F)ECneut.getData(5).get(0)).fill(neut_mom, ecal_neut_beta);        
               	((H1F)ECneut.getData(7).get(0)).fill(mass2);
               	if(mass2cut) {
            		float nphi = newPhi(ecal_neut_phi);
            		cx = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
            		cy = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
                	((H2F)ECneut.getData(9).get(0)).fill(-cx,cy);
           			((H1F)ECneut.getData(10).get(1)).fill(neut_mom);
            	}
            }            
        }
		
    }
       
    public void fillECphot() {
		int run = getRunNumber();
		DataGroup ECphot = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECphot"),run);
    	int pin = -1;
    	for (Particle p : phot_ecal) {
        	int pindex     = (int)   p.getProperty("pindex"); 
        	if (pindex!=pin) {
        		ecal_phot_sec  = (int)   p.getProperty("sector");
        		ecal_phot_the  = (float) Math.toDegrees(p.theta());
        		ecal_phot_phi  = (float) Math.toDegrees(p.phi());
        		ecal_phot_beta = (float) p.getProperty("beta");   
        		pin = pindex;
        		if(ecal_phot_sec == e_sect) ((H2F)ECphot.getData(0).get(0)).fill(p.e(),ecal_phot_the);
        	}
    	}
    }
    
    public float newPhi(float phi) {
    	float newphi=phi;
    	if(phi<-30) newphi=360+phi;
    	if(newphi<30) return newphi;
    	if(newphi>30&&newphi<90)   return newphi-60;
    	if(newphi>90&&newphi<150)  return newphi-120;
    	if(newphi>150&&newphi<210) return newphi-180;
    	if(newphi>210&&newphi<270) return newphi-240;
    	if(newphi>270&&newphi<330) return newphi-300;
    	return phi;
    }
    
    //PLOT
    
    @Override       
    public void plotHistos(int run) {
    	setRunNumber(run); 
    	ECkinPlot("ECkin");
    	ECelecPlot("ECelec"); 
        SCelecPlot("SCelec");
        ECprotPlot("ECprot");
        ECpbarPlot("ECpbar");
        ECpipPlot("ECpip");
        ECpimPlot("ECpim");
        ECpi0Plot("ECpi0");
        ECneutPlot("ECneut");
        ECphotPlot("ECphot");
        if(!isAnalyzeDone) return;
        showNeutronEff();
        showPi0Eff();
    }
    
    public void ECkinPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6)        plot123(index);    	
    }
    
    public void ECelecPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6)        plot123(index);
        if(getActive123()>5) ECelecPlotFits(index);
    }
    
    public void SCelecPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6)        plot123(index);
        if(getActive123()>5) SCelecPlotFits(index);
    }
	
	public void ECprotPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
		if(getActive123()<2) plot123(index);		
	}
	
	public void ECpbarPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
		if(getActive123()<2) plot123(index);		
	}
	
	public void ECpi0Plot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
		if(getActive123()<2) plot123(index);		
	}
	
	public void ECpipPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
		if(getActive123()<2) plot123(index);	
	}
    
    public void ECpimPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6)        plot123(index);    	
        if(getActive123()>5)  ECpimPlotFits(index);
    }
	
	public void ECneutPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6) plot123(index);   
	}
	
	public void ECphotPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6) plot123(index);    
	}	
	
    @Override
    public void plotEvent(DataEvent de) {
       analyze();
    }
    
    public void analyze() {    
       System.out.println(getDetectorName()+".Analyze() ");
       
       if(!dropSummary) {
    	   analyzeECelec();
    	   analyzeECpim();
    	   analyzeECneut();    	             
       }
       isAnalyzeDone = true;
    }
    
    public void analyzeECelec() {
        getECelecSummary(1,7,0,3,0,3);
        getSCelecSummary(1,7,0,2,0,3);
    }
    
    public void analyzeECpim() {
        getECpimSummary(1,7,0,3,0,3);    	
    }
    
    public void analyzeECneut() {
        showNeutronEff();    	
    }
    	
	public void getRADYLD(int run) {
		DataGroup dg5 = this.getDataGroup().getItem(0,5,0,run);	
		GraphErrors g1 = ((GraphErrors)dg5.getData(0).get(0)); g1.reset();
		GraphErrors g2 = ((GraphErrors)dg5.getData(1).get(0)); g2.reset();
		GraphErrors g3 = ((GraphErrors)dg5.getData(2).get(0)); g3.reset();
		for (int i=0; i<6; i++) g1.addPoint(i,counter[i][0],0f,0f);
		for (int i=0; i<6; i++) g2.addPoint(i,counter[i][1],0f,0f);
		for (int i=0; i<6; i++) g3.addPoint(i,counter[i][2],0f,0f);	
	}
	
	public void ECelecPlotFits(int index) {
		int icmax=3, ilmax=3;
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.clear();
        if(getActive123()==6) {c.divide(6,3); icmax=3;}
        if(getActive123()==7) {c.divide(6,4); ilmax=2; icmax=2;}
        
        int n=0;
        for (int il=0; il<ilmax; il++) {
        	for (int ic=0; ic<icmax; ic++ ) {
        		for (int is=1; is<7; is++) { int isk = is+10*index;
        			if (getActive123()==6 && getActiveLayer()==il) {
        				c.cd(n); H1F h = tl.fitData.getItem(isk,ic,il,getRunNumber()).getHist(); 
        				h.setOptStat((ic<2)?"0":"1100"); if(ic==0) h.setTitle("SECTOR "+is);
        				h.setTitleX("DC - "+det[il]+xyz[ic]+" RESIDUAL (CM)");
        				h.setFillColor(4);
        				c.draw(h);n++;
        				if(ic<2) c.draw(tl.fitData.getItem(isk,ic,il,getRunNumber()).getGraph(),"same");
        			}
        			if (getActive123()==7) {
           		     	c.cd(n); H1F h = tl.fitData.getItem(isk,ic,il+3,getRunNumber()).getHist(); 
            			h.setOptStat("0"); if(il==0&&ic==0) h.setTitle("SECTOR "+is);
            			h.setTitleX("e - #gamma "+(ic==0?" #Delta#theta ":" #Delta#phi ")+det[il]+" (DEG)");
            			h.setFillColor(4);
            			c.draw(h);n++;
            			c.draw(tl.fitData.getItem(isk,ic,il+3,getRunNumber()).getGraph(),"same");
        			}
        		}
        	}
        }

	}
	
	public void ECpimPlotFits(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.clear();
        c.divide(6,3);
        
        int n=0;
        for (int il=0; il<3; il++) {
        	for (int ic=0; ic<3; ic++ ) {
        		for (int is=1; is<7; is++) { int isk = is+10*index;
        			if (getActive123()==6 && getActiveLayer()==il) {
        				c.cd(n); H1F h = tl.fitData.getItem(isk,ic,il,getRunNumber()).getHist(); 
        				h.setOptStat((ic<2)?"0":"1100"); if(ic==0) h.setTitle("SECTOR "+is);
        				h.setTitleX("DC - "+det[il]+xyz[ic]+" RESIDUAL (CM)");
        				h.setFillColor(4);
        				c.draw(h);n++;
        				if(ic<2) c.draw(tl.fitData.getItem(isk,ic,il,getRunNumber()).getGraph(),"same");
        			}
        		}
        	}
        }

	} 
	
	public void SCelecPlotFits(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.clear();
        c.divide(6,3);
        
        int n=0;
        for (int il=0; il<2; il++) {
        	for (int ic=0; ic<3; ic++ ) {
        		for (int is=1; is<7; is++) { int isk = is+10*index;
        			if (getActive123()==6 && getActiveLayer()==il) {
        				c.cd(n); H1F h = tl.fitData.getItem(isk,ic,il,getRunNumber()).getHist(); 
        				h.setOptStat((ic<2)?"0":"1100"); if(ic==0) h.setTitle("SECTOR "+is);
        				h.setTitleX("DC - "+scdet[il]+xyz[ic]+" RESIDUAL (CM)");
        				h.setFillColor(4);
        				c.draw(h);n++;
        				if(ic<2)c.draw(tl.fitData.getItem(isk,ic,il,getRunNumber()).getGraph(),"same");
        			}
        		}
        	}
        }

	} 	
	
    public void getECelecSummary(int is1, int is2, int il1, int il2, int ic1, int ic2) {
    	
    	float p1=0,p2=0,f1=0,f2=0;
        int run=getRunNumber(), k=getDetectorTabNames().indexOf("ECelec");
        
        float plim[][] = {{-0.5f,-2.2f,-0.8f,-2.2f,-0.8f,-2.2f},{0.5f,2f,0.6f,2f,0.6f,2f}};
        float flim[][] = {{-1.5f,-4f,-1.5f,-4f,-1.5f,-4f},{1.5f,4f,1.5f,3f,1.5f,3f}};
        
        cfitEnable = true;
        DataGroup dg2 = this.getDataGroup().getItem(0,2,k,run);	
        DataGroup dg5 = this.getDataGroup().getItem(0,5,k,run);	
                
        for (int il=il1; il<il2; il++) {            
           DataGroup  dg1 = this.getDataGroup().getItem(il,1,k,run);	
     	   GraphErrors g1 = ((GraphErrors)dg5.getData(il  ).get(0)); g1.reset();
     	   GraphErrors g2 = ((GraphErrors)dg5.getData(il  ).get(1)); g2.reset();
     	   GraphErrors g3 = ((GraphErrors)dg5.getData(il+3).get(0)); g3.reset();
     	   GraphErrors g4 = ((GraphErrors)dg5.getData(il+3).get(1)); g4.reset();
     	   g1.setTitle("STRIP WIDTH "+((il==0)?4.5:10)+" CM"); 
     	   g1.setTitleX("SECTOR"); g1.setTitleY("DC-"+det[il]+" RESIDUALS (CM)");
     	   g3.setTitleX("SECTOR"); g3.setTitleY("ELEC - #gamma  "+det[il]+" (DEG)");
     	   for (int ic=ic1; ic<ic2; ic++) {
        	  if(ic<2) {p1=plim[0][ic+il*2]; p2=plim[1][ic+il*2]; f1=flim[0][ic+il*2]; f2=plim[1][ic+il*2];}
    		  for (int is=is1; is<is2; is++) { int isk = is+10*k;        		 
    			  if(!dg1.getData(is-1+ic*6).isEmpty()) {
            	  tl.fitData.add(fitEngine(((H2F)dg1.getData(is-1+ic*6).get(0)).projectionY(),1,-3.5,3.5,-3.5,3.5),isk,ic,il,run);
            	  if(ic<2) tl.fitData.add(fitEngine(((H2F)dg2.getData(is-1+ic*6+il*12).get(0)).projectionY(),3,p1,p2,f1,f2),isk,ic,il+3,run);                  
             	  if(ic==0) g1.addPoint(is,     tl.fitData.getItem(isk,ic,il,run).mean,0,tl.fitData.getItem(isk,ic,il,run).sigma);
             	  if(ic==1) g2.addPoint(is+0.2, tl.fitData.getItem(isk,ic,il,run).mean,0,tl.fitData.getItem(isk,ic,il,run).sigma);
             	  if(ic==0) g3.addPoint(is,     tl.fitData.getItem(isk,ic,il+3,run).mean,0,tl.fitData.getItem(isk,ic,il+3,run).sigma);
             	  if(ic==1) g4.addPoint(is+0.2, tl.fitData.getItem(isk,ic,il+3,run).mean,0,tl.fitData.getItem(isk,ic,il+3,run).sigma);
    			  }
    		  }
           }
        } 
    }
    
    public void getECpimSummary(int is1, int is2, int il1, int il2, int ic1, int ic2) {
    	
        int run=getRunNumber(), k=getDetectorTabNames().indexOf("ECpim");
        
        cfitEnable = true;
        DataGroup dg5 = this.getDataGroup().getItem(0,5,k,run);	
        
        for (int il=il1; il<il2; il++) {            
           DataGroup dg1 = this.getDataGroup().getItem(il,1,k,run);	
      	   GraphErrors g1 = ((GraphErrors)dg5.getData(il  ).get(0)); g1.reset();
      	   GraphErrors g2 = ((GraphErrors)dg5.getData(il  ).get(1)); g2.reset();
     	   g1.setTitle("STRIP WIDTH "+((il==0)?4.5:10)+" CM"); 
     	   g1.setTitleX("SECTOR"); g1.setTitleY("DC-"+det[il]+" RESIDUALS (CM)");
      	   for (int ic=ic1; ic<ic2; ic++) {
        		  for (int is=is1; is<is2; is++) { int isk = is+10*k;
    			  if(!dg1.getData(is-1+ic*6).isEmpty()) {
                	  tl.fitData.add(fitEngine(((H2F)dg1.getData(is-1+ic*6).get(0)).projectionY(),1,-3.5,6.5,-3.5,6.5),isk,ic,il,run);
                 	  if(ic==0) g1.addPoint(is,     tl.fitData.getItem(isk,ic,il,run).mean,0,tl.fitData.getItem(isk,ic,il,run).sigma);
                 	  if(ic==1) g2.addPoint(is+0.2, tl.fitData.getItem(isk,ic,il,run).mean,0,tl.fitData.getItem(isk,ic,il,run).sigma);
    			  }
     		  }
      	   }
        }
   	
    }  
    
    public void getSCelecSummary(int is1, int is2, int il1, int il2, int ic1, int ic2) {
    	
        int run=getRunNumber(), k=getDetectorTabNames().indexOf("SCelec");
        
        cfitEnable = true;
        DataGroup dg5 = this.getDataGroup().getItem(0,5,k,run);	
        
        for (int il=il1; il<il2; il++) {            
           DataGroup  dg1 = this.getDataGroup().getItem(il,1,k,run);	
      	   GraphErrors g1 = ((GraphErrors)dg5.getData(il  ).get(0)); g1.reset();
      	   GraphErrors g2 = ((GraphErrors)dg5.getData(il  ).get(1)); g2.reset();
     	   for (int ic=ic1; ic<ic2; ic++) { 
        		  for (int is=is1; is<is2; is++) { int isk = is+10*k;
    			  if(!dg1.getData(is-1+ic*6).isEmpty()) {
                	  tl.fitData.add(fitEngine(((H2F)dg1.getData(is-1+ic*6).get(0)).projectionY(),3,-4.5,4.5,-4.5,4.5),isk,ic,il,run);
                 	  if(ic==0) g1.addPoint(is,     tl.fitData.getItem(isk,ic,il,run).mean,0,tl.fitData.getItem(isk,ic,il,run).sigma);
                 	  if(ic==1) g2.addPoint(is+0.2, tl.fitData.getItem(isk,ic,il,run).mean,0,tl.fitData.getItem(isk,ic,il,run).sigma);
    			  }
     		  }
      	   }
        }
   	
    }
    
    public boolean neutFiduc(int opt, float cx, float cy) {
    	float off = 0.52f;
    	if(opt==1) return (cx>0.07 && Math.sqrt((cx+off)*(cx+off)+cy*cy)<(0.52+off) && Math.abs(cy)<0.57*(cx-0.07));
    	if(opt==2) return Math.sqrt((cx-0.35)*(cx-0.35)+cy*cy)<0.1;
    	return false;    
    }
    
    public GraphErrors getEff(int index) {
		DataGroup dg = this.getDataGroup().getItem(0,1,index,getRunNumber());    	
    	GraphErrors geff = new GraphErrors();
    	geff.getAttributes().setTitleX("p_mm (GeV)");
    	geff.getAttributes().setTitleY("Efficiency");
    	GraphErrors g1 = ((H1F)dg.getData(10).get(0)).getGraph();
    	GraphErrors g2 = ((H1F)dg.getData(10).get(1)).getGraph();
    	for (int ix=0; ix<g1.getDataSize(0); ix++) {
    		float x1 = (float) g1.getDataX(ix); float y1 = (float) g1.getDataY(ix);
    		if(y1>0) {
        		float   y2 = (float)g2.getDataY(ix);
    			float  eff = y2/y1;
//    			float effe = (float) Math.sqrt((y1-1)*eff*(1-eff))/(y1-1); 
    			float effe = (float) Math.sqrt(eff*(1-eff)/(y1-1)); 
    			geff.addPoint(x1, eff, 0f, effe);
    		}
    	}
    	return geff;
    } 
    	
    public Point3D getResidual(Particle p) {    	
 //       float  dx = ((float)p.getProperty("hx")-(float)p.getProperty("x"));
 //       float  dy = ((float)p.getProperty("hy")-(float)p.getProperty("y"));
 //       float  dz = ((float)p.getProperty("hz")-(float)p.getProperty("z"));
        float  dx = p.hasProperty("tx")?((float)p.getProperty("tx")-(float)p.getProperty("x")):0;
        float  dy = p.hasProperty("ty")?((float)p.getProperty("ty")-(float)p.getProperty("y")):0;
        float  dz = p.hasProperty("tz")?((float)p.getProperty("tz")-(float)p.getProperty("z")):0;
        Point3D xyz = new Point3D(dx,dy,dz);
        xyz.rotateZ(Math.toRadians(-60*(p.getProperty("sector")-1)));
        xyz.rotateY(Math.toRadians(-25)); 	
        return xyz;
    }
	
    public double Vangle(Vector3 v1, Vector3 v2){
        double res = 0;
        double l1 = v1.mag();
        double l2 = v2.mag();
        double prod = v1.dot(v2);
        if( l1 * l2 !=0 && Math.abs(prod)<l1*l2 )res = Math.toDegrees( Math.acos(prod/(l1*l2) ) );
        return res; 
    }
    
    public H1F makeH1(String name, String tag, int nx, double x1, double x2, String tit, String titx) {
    	H1F h1 = new H1F(name+tag,name+tag,nx,x1,x2);
    	    if(tit!="") h1.setTitle(tit);
    	    h1.setTitleX(titx); 
    	    return h1;
    }
    
    public H2F makeH2(String name, String tag, int nx, double x1, double x2, int ny, double y1, double y2, String tit, String titx, String tity) {
    	H2F h2 = new H2F(name+tag,name+tag,nx,x1,x2,ny,y1,y2);
    	    if(tit!="") h2.setTitle(tit);
    	    h2.setTitleX(titx); h2.setTitleY(tity);
    	    return h2;
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
    
    public void plot123(int index) {      	
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveLayer(),getActive123(),index,getRunNumber()));   	
    } 
    
    public void plotUVW(int index) {  
      	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,getActive123()-1,index,getRunNumber()));
    }  
    
    public void showNeutronEff() {
		int index = getDetectorTabNames().indexOf("ECneut");
		EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.cd(11);c.getPad().getAxisX().setRange(0., 2.5); c.getPad().getAxisY().setRange(0., 1.); c.draw(getEff(index)); c.draw(neuteff,"same");
    }
    
    public void showPi0Eff() {
		int index = getDetectorTabNames().indexOf("ECpi0");
		EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.cd(11);c.getPad().getAxisX().setRange(0., 5.5); c.getPad().getAxisY().setRange(0., 1.); c.draw(getEff(index));   	
    }
    
    @Override
    public void timerUpdate() {

    	if(nc>100) {
    		showNeutronEff();
    		showPi0Eff();
    		nc=0;
    	}
    	nc++;
    	return;
    }
    
}
