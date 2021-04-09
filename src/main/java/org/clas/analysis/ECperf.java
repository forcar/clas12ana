package org.clas.analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Map;


import org.clas.tools.Event;

import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.geom.prim.Point3D;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;

import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECperf extends DetectorMonitor {
	
    Event    ev = new Event();
	
	DataGroup dg = null;
	int nc=0;
	
	int trigger_sect;
	
	public boolean goodE,goodP,goodPROT,goodPBAR,goodPIP,goodPIM,goodNEUT,goodPHOT,goodPI0;
	public int nE,nP,nPROT,nPBAR,nPIP,nPIM,nNEUT,nPHOT;
	
	public boolean taggedPI0, taggedETA;
	String[]   det = new String[]{"PCAL ","ECIN ","ECOU "};
	String[]  det5 = new String[]{"P1B","P1A","PCAL ","ECIN ","ECOU "};
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
	public float[] lU= {0,0,0},lV={0,0,0},lW={0,0,0};
	public int[]   iU= {0,0,0},iV={0,0,0},iW={0,0,0};
	public float cZ;
	public int piU,piV,piW;
	
	public LorentzVector VB, VT, Ve, VGS, Vprot, Vpbar, Vpip, Vpim;
	public int   e_part_ind, e_sect, e_FTOF_pad1a, e_FTOF_pad1b, e_HTCC_bin_phi, e_HTCC_bin_theta;
	public float e_mom, e_the, e_phi, ephi, e_vx, e_vy, e_vz, e_cz, e_x, e_y, e_t, e_pat;
	public float e_xB, e_Q2, e_W;
	
	public int htcc_sect;
	public float htcc_x,htcc_y,htcc_npe;

	public int   p_sect;
	public float p_mom, p_the, p_phi, p_vx, p_vy, p_vz, p_cz, p_x, p_y;
	
	public float e_ecal_esum,e_ecal_pcsum,e_ecal_ecsum,e_ecal_rat,e_ecal_sf;
	public float p_ecal_esum,p_ecal_pcsum,p_ecal_ecsum;
        
	public int   prot_part_ind;
	public float prot_mom, prot_the, prot_phi, prot_vx, prot_vy, prot_vz, prot_beta, prot_beta_pcal, prot_beta_ecin;
	public float pbar_mom, pbar_the, pbar_phi, pbar_vx, pbar_vy, pbar_vz, pbar_beta;

	public int   pim_part_ind, pip_part_ind, pip_FTOF_pad1b;
	public float pip_mom, pip_the, pip_phi, pip_vx, pip_vy, pip_vz, pip_beta, pip_x, pip_y, pip_z;
	public float pip_ecal_esum;
	
	public float pim_mom, pim_the, pim_phi, pim_vx, pim_vy, pim_vz, pim_beta;

	public float elast_dPhi, elast_EB;
	public float epip_dPhi, epip_MM, ep_dPhi, epbar_dPhi, ep_MM, epbar_MM;
//	public float pi0_mass, pi0_mom, pi0_e, pi0_the, pi0_phi, pi0_open, pi0_cx, pi0_cy;
	public float nm_mass, nm_mom, nm_e, nm_the, nm_phi, nm_open, nm_cx, nm_cy, nm_ggp;
	
	public float ecal_pi0_mass, ecal_pi0_mom, ecal_pi0_e, ecal_pi0_the, ecal_pi0_phi;
	public float ecal_pi0_opa, ecal_pi0_cx, ecal_pi0_cy, ecal_pi0_X;	
	public float neut_mom,neut_the,neut_phi,neut_cx,neut_cy;
	public float ecal_neut_the,ecal_neut_phi,ecal_neut_beta,ecal_neut_cx,ecal_neut_cy;
	public float ecal_phot_the,ecal_phot_phi,ecal_phot_beta,ecal_phot_nrg;
	public int   ecal_neut_sec, ecal_phot_sec, ecal_pi0_sec;
	public int[] ecal_neut_esum = new int[6];
	public int[] ecal_phot_esum = new int[6];
	public int   newPhiSector;
	
	public IndexedList<Float> elec_ecal_resid = new IndexedList<Float>(3);
	public IndexedList<Float> posi_ecal_resid = new IndexedList<Float>(3);
	public List<Particle>   e_ecal  = new ArrayList<Particle>();
	public List<Particle>   p_ecal  = new ArrayList<Particle>();
	public List<Particle> pim_ecal  = new ArrayList<Particle>();
	public List<Particle> pip_ecal  = new ArrayList<Particle>();
	public List<Particle> prot_ecal = new ArrayList<Particle>();
	public List<Particle> pbar_ecal = new ArrayList<Particle>();
	public List<Particle> phot_ecal = new ArrayList<Particle>();
	public List<Particle> neut_ecal = new ArrayList<Particle>();
	public List<NeutralMeson> nm_ecal = new ArrayList<NeutralMeson>();
//	public IndexedList<List<Float>>  ecal_rad = new IndexedList<List<Float>>(3);
	public IndexedList<List<Particle>> ecphot      = new IndexedList<List<Particle>>(2);	
	public IndexedList<List<Particle>> ecphot_filt = new IndexedList<List<Particle>>(1);	
	public IndexedList<List<Particle>>  ecpim      = new IndexedList<List<Particle>>(2);	
	public IndexedList<List<Particle>>  ecpim_filt = new IndexedList<List<Particle>>(1);	
	public IndexedList<Float> elec_ftof_resid = new IndexedList<Float>(3);
	public IndexedList<Float> posi_ftof_resid = new IndexedList<Float>(3);
	
	public float[][]  counter = new float[6][3];
	
	public IndexedTable rfTable;	
	GraphErrors neuteff = null;
	
	public ECperf(String name) {
        super(name);
        this.setDetectorTabNames("ECkin",
        		                 "ECelec",
        		                 "ECposi",
        		                 "ECtime",
        		                 "SCelec",
        		                 "ECprot",
        		                 "ECpbar",
        		                 "ECpip",
        		                 "ECpim",
        		                 "ECmip",
        		                 "ECpi0",
        		                 "ECeta",
        		                 "ECneut",
        		                 "ECphot",
        		                 "ECtrig");

        this.useCALButtons(true);
        this.use123Buttons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();      
    }
    
    public void localinit() {
    	System.out.println("ECperf.localinit()");
    	config="phot";
    	configEngine("phot");
    	tl.setFitData(Fits);  
//        part.setGeom("2.5");  
//        part.setConfig("pi0");  
//        part.setGoodPhotons(1212);    	
        neuteff = getGraph(outPath+"files/neuteff.vec",50); neuteff.setMarkerColor(1); neuteff.setLineColor(1);
    }  
    
    public void localclear() {
    	System.out.println("ECperf:localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	slider.setValue(0);
    } 
    
    public void dstinit(int run) {
    	System.out.println("ECperf:dstinit("+run+")");
		Mp = 0.93827f;
        Eb = EB = getBeamEnergy(run);
        System.out.println("Eb="+Eb+" run="+runNum);
		rfPeriod = 4.008f;       
		rf_large_integer = 1000;    	
        rfTable = ebcm.getConstants(run,"/calibration/eb/rf/config");
        if (rfTable.hasEntry(1, 1, 1)){
            System.out.println(String.format("RF period from ccdb for run %d: %f",runNum,rfTable.getDoubleValue("clock",1,1,1)));
            rfPeriod = (float) rfTable.getDoubleValue("clock",1,1,1);
        }
        VB = new LorentzVector(0,0,Eb,Eb);
        VT = new LorentzVector(0,0,0,Mp);    
    	System.out.println("ECperf:dstinit complete");
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
    	createECkin(3);
    	createECkin(4);
    	createECkin(5);
    	createECkin(6);
    	createECelec(0);
    	createECelec(1);
    	createECelec(2);
    	createECelec(3);
    	createECelec(4);
    	createECelec(5);
    	createECposi(0);
    	createECtime(0);
    	createECtime(3);
    	createECtime(4);
    	createSCelec(1);
    	createSCelec(5);
    	createECprot(0);
    	createECprot(1);
    	createECpbar(0);
    	createECpbar(1);
    	createECpim(0);    	
    	createECpim(1);
    	createECpim(5);
    	createECpip(0);
    	createECpip(1);
    	createECpip(2);
    	createECpip(3);
    	createECpip(4);
    	createECpip(5);
    	createECneut(0);
    	createECneut(1);
    	createECphot(0);
    	createECphot(1);
    	createECpi0(0);
    	createECpi0(1);
    	createECeta(1);
    	createECtrig(0);
    }
    
    public void processEvent(DataEvent event) {
    	
        if(dropBanks) dropBanks(event);
            	
    	ev.init(event);
    	ev.setHipoEvent(isHipo3Event);
    	ev.setEventNumber(getEventNumber());
//    	ev.requireOneElectron(!event.hasBank("MC::Event"));
    	ev.requireOneElectron(false);
   	    ev.setElecTriggerSector(ev.getElecTriggerSector());
         
        fillECtrig(0);        
    	if(!ev.procEvent(event)) return;
        fillECtrig(1); 
        
    	e_ecal    = ev.getPART(0.1, 11);      nE = e_ecal.size();
    	p_ecal    = ev.getPART(0.1, 12);      nP = p_ecal.size();    	
	    pim_ecal  = ev.getPART(0.5, 212);   nPIM = pim_ecal.size(); 
	    
	    if (nE==1)            fillECtrig(2);
	    if (nE==0)            fillECtrig(3);
	    if (nPIM>0)           fillECtrig(4);
	    if (nPIM>0 && nE==0) {fillECtrig(5); fillECpim();}
//	    if (nPIM>0 && nE>0) {fillECtrig(5); fillECpim();}
	    
	    if (!makeELEC()) return;
        
	    prot_ecal = ev.getPART(0.2,2212);    nPROT = prot_ecal.size();
	    pbar_ecal = ev.getPART(0.2,2213);    nPBAR = pbar_ecal.size();
	    pip_ecal  = ev.getPART(0.2, 211);     nPIP = pip_ecal.size();
	    neut_ecal = ev.getPART(0.043,2112);  nNEUT = neut_ecal.size(); 
	    phot_ecal = ev.getPART(0.001,22);    nPHOT = phot_ecal.size();
	    
		DataGroup ECnm0   = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECpi0"),getRunNumber());
        if (ev.pmc.size()==2) {
	    	double  opa = ev.pmv.get(0).theta(ev.pmv.get(1));
	    	Vector3 ggc = ev.pmv.get(0).cross(ev.pmv.get(1));
	    	nm_ggp = (float) Math.toDegrees(Math.atan2(ggc.y(),ggc.x()));  
	    	if(nm_ggp<0) nm_ggp=-nm_ggp;
	    	nm_ggp=nm_ggp-90;
	    	if(nm_ggp<0) nm_ggp=nm_ggp+180;
	    	

	    	for (Particle p : phot_ecal) {
	    		((H1F)ECnm0.getData(4).get(0)).fill(p.getProperty("beta"));    	
	    	}
        }
	    
	    filterEvent();	
	    fillHists();
    } 
  
    void filterEvent() {
	    goodP    = nP>0;
	    goodPIP  = nPIP>0;
	    goodPIM  = nPIM>0;
	    goodPROT = nPROT>0;
	    goodPBAR = nPBAR>0;
	    goodNEUT = nNEUT>0;
	    goodPHOT = nPHOT>0;	    
	    ecphot.clear(); ecphot_filt.clear(); ecpim_filt.clear();
	    if(goodP) makePOSI();
	    if(goodPHOT) {
	    	ecphot      = getECALClusters(phot_ecal); 
	    	ecphot_filt = filterECALClusters(ecphot,2,22); 
	    }
	    if(goodPIM) {
	    	ecpim      = getECALClusters(pim_ecal);
	    	ecpim_filt = filterECALClusters(ecpim,4,212);
	    }
    } 
    
    public IndexedList<List<Particle>> getECALClusters(List<Particle> list) {
    	
//    	System.out.println(" ");
//    	System.out.println("Event Number "+ev.eventNumber);
        IndexedList<List<Particle>> olist = new IndexedList<List<Particle>>(2);       
    	for (Particle p : list) {
    		int ip = (int)p.getProperty("pindex");
    		for (Particle ec : ev.getECAL(ip)) {
    			int is = (int) ec.getProperty("sector");
    			int il = (int) ec.getProperty("layer");
    			float en = (float) ec.getProperty("energy");
//    			if(ev.part.get(ip).getProperty("ppid")==22) System.out.println(ip+" "+is+" "+il+" "+en);
    			if (!olist.hasItem(is,il)) olist.add(new ArrayList<Particle>(), is,il); 
    			     olist.getItem(is,il).add(ec);
    		}
    	}
    	return olist;
    }
    
    public IndexedList<List<Particle>> filterECALClusters(IndexedList<List<Particle>> list, int n, int pid) {
    	
        IndexedList<List<Particle>> olist = new IndexedList<List<Particle>>(1);       
		IndexGenerator ig = new IndexGenerator();
		
    	for (Map.Entry<Long,List<Particle>>  entry : list.getMap().entrySet()){
			int is = ig.getIndex(entry.getKey(), 0); 
			int ip = (int) entry.getValue().get(0).getProperty("pindex");
			if(entry.getValue().size()>0 && entry.getValue().size()<n && Math.abs(ev.part.get(ip).getProperty("ppid"))==pid) {
				if(!olist.hasItem(is)) olist.add(new ArrayList<Particle>(), is); 
				    Particle p = entry.getValue().get(0);
	    			int il = (int) p.getProperty("layer");
	    			float en = (float) p.getProperty("energy");
//	    			System.out.println("Filter "+ip+" "+is+" "+il+" "+en);
				    olist.getItem(is).add(p);
//				    for (Particle p : olist.getItem(is)) {
//				    	olist.getItem(is).add(p);
//				    }
			}
    	}
    	return olist;    	
    }
    
	public void fillHists(){
				
		fillECkin();
		fillECphot0(); 
		fillECelec(); 
		fillECposi();
		fillECtime();
		fillSCelec();
    	fillTrajXY();
		
	    if(goodPIP) { // FC pi+
	    	if(select_epip()) { //(e' pi+) tagged neutron
	    		fillECpip();
	    		fillECneut();
	    		fillECphot1();	    		
	    	}
	    	fillECpipFocus();
	    }

	    if(goodPROT && !goodPIP && !goodPIM) { 
	    	if(select_ep()) { // (e'p) tagged neutral meson nm
	    		fillECprot();
	    		if(taggedPI0) fillECnm("ECpi0");
	    		if(taggedETA) fillECnm("ECeta");
	    	}
	    	if(select_epbar()) fillECpbar();
	    }
	}
	    
    @Override       
    public void plotHistos(int run) {
    	setRunNumber(run); 
    	
    	ECkinPlot("ECkin");
    	ECelecPlot("ECelec"); 
    	ECposiPlot("ECposi");
        ECtimePlot("ECtime");
        SCelecPlot("SCelec");
        ECprotPlot("ECprot");
        ECpbarPlot("ECpbar");
        ECpipPlot("ECpip");
        ECpimPlot("ECpim");
        ECpi0Plot("ECpi0");
        ECpi0Plot("ECeta");
        ECneutPlot("ECneut"); 
        ECphotPlot("ECphot");
        ECtrigPlot("ECtrig");
        if(!isAnalyzeDone) return;
        showNeutronEff();
//        showPi0Eff();
    }	
	
// HISTOS	   
    public void createECkin(int st) {
    	
    	String tab = "ECkin", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
		F1D f1,f2;
		String[] luvw = {" LU"," LV"," LW"}; String[] uvw = {" U"," V"," W"}; 
		int[]     npmt = {68,62,62,36,36,36,36,36,36};
		
    	switch (st) {        
        case 0: 
        dg = new DataGroup(6,4);
        for(int is=1;is<7;is++){    
	        tag = is+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_1_",tag,60,5.0,12, 60,EB*0.75,EB*0.95,"","#theta (^o)","p (GeV)"),is-1);
	        dg.addDataSet(makeH2(tab+"_2_",tag,60,0.6,1.2,60,EB*0.75,EB*0.95,"","W (GeV)","p (GeV)"),is-1+6);
	        dg.addDataSet(makeH2(tab+"_3_",tag,60,0.6,1.2,60,5,12,   "","W (GeV)","#theta (^o)"),is-1+12);
	        dg.addDataSet(makeH2(tab+"_4_",tag,60,0.6,1.2,60,-20,30, "","W (GeV)","#phi (^o"),is-1+18);
        }
		break;		
        case 1:
        dg = new DataGroup(6,3);
        for(int is=1; is<7; is++) {
	        tag = is+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_1_",tag,50,0,100,30,0.15,0.3,"Sector "+is,"LU (cm)","E/P"),is-1);
	        dg.addDataSet(makeH2(tab+"_2_",tag,50,0,100,30,0.15,0.3," ",         "LV (cm)","E/P"),is-1+6);
	        dg.addDataSet(makeH2(tab+"_3_",tag,50,0,100,30,0.15,0.3," ",         "LW (cm)","E/P"),is-1+12);        	
        }
        break;        
        case 2: 
        f1 = new F1D("H_e_EC_resid0_f+"+run,"[a]",1,69);  f1.setParameter(0,0f); f1.setLineWidth(1);
        f2 = new F1D("H_e_EC_resid1_f+"+run,"[a]",0,EB/4);f2.setParameter(0,0f); f2.setLineWidth(1);
        dg = new DataGroup(6,4);
        for(int is=1;is<7;is++){
	        tag = is+"_"+st+"_"+k+"_"+run;
			dg.addDataSet(makeH2(tab+"_1_",tag,60,0,EB,60,0.12,0.35,"","p (GeV)","E/P"),is-1);
			dg.addDataSet(makeH2(tab+"_2_",tag,68,1,69,50,0.12,0.35,"","PCAL U STRIP","E/P"),is-1+6);
			dg.addDataSet(makeH2(tab+"_3_",tag,68,1,69,50,-5,4,     "","PCAL U STRIP","#chi PID"),is-1+12); dg.addDataSet(f1,is-1+12);
//			dg.addDataSet(makeH2(tab+"_4_",tag,70,0,EB/6,70,0,EB/5, "","PC (GeV)","EC (GeV)"),is-1+18);	
			dg.addDataSet(makeH2(tab+"_4_",tag,68,1,69,50,0,1, "","PCAL U STRIP","PCAL / ETOT"),is-1+18); dg.addDataSet(f2,is-1+18);	
        }
		break;		
        case 3:
        dg = new DataGroup(6,3);
        for(int is=1; is<7; is++) {
	        tag = is+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_1_",tag,50,-10,20,80,0,40,"Sector "+is,"#phi (^o)","LU (cm)"),is-1);
	        dg.addDataSet(makeH2(tab+"_2_",tag,50,-20,25,80,0,40," ",         "#phi (^o)","LV (cm)"),is-1+6);
	        dg.addDataSet(makeH2(tab+"_3_",tag,50,-10,35,80,0,40," ",         "#phi (^o)","LW (cm)"),is-1+12);        	
        }                
    	break;
		case 4:
		dg = new DataGroup(6,6);    	
		String lab4[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};    
		f1 = new F1D("H_e_EC_rad5_f+"+run,"[a]",-20,30); 
		f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1); int in=0;
		for(int i=0;i<3;i++) { //pcal,ecin,ecou
			for(int n=0; n<2; n++) { //thet,phi
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,-20,30,60,n==0?-5:-15,n==0?10:20,"","S"+is+" #phi_e",lab4[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
				}
			}	
		}
		break;
		case 5:
		for(int i=0;i<1;i++) { //pcal
			dg = new DataGroup(6,3); in=0; 
			for(int j=0; j<3; j++) { //lu,lv,lw	
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+j+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,0,100,60,-20,25,"","S"+is+luvw[j],"PCAL-ECIN #Delta"+luvw[j]), in);in++;  
				}
			}
			this.getDataGroup().add(dg,i,st,k,run);
		}
		return;
		case 6:  
		dg = new DataGroup(6,6); in=0;	
		for(int i=0;i<2;i++) { //pcal, ecin
			for(int j=0; j<3; j++) { //lu,lv,lw	
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+j+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,70,0,420,npmt[3*i+j],0,npmt[3*i+j]+1,"","S"+is+luvw[j],det[i]+uvw[j]+" STRIP"+luvw[j]), in);in++;  
				}
			}
		}
    	}    	    	
    	this.getDataGroup().add(dg,0,st,k,run);              	    	
    }
    
    public void createECelec(int st) {
    	
    	String tab = "ECelec", tag = null;
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	 
        F1D     f1,f2;
		
    	switch (st) {
    	
        case 0:   		
		f1 = new F1D("H_e_EC_luvw_f1+"+run,"[a]",0,430); 
		f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
		String[] xtxt0 = {" LU"," LV"," LW"}; 

		for(int i=0;i<3;i++) { //pcal,ecin,ecou
	        dg = new DataGroup(6,6); in=0; float y = i<1?5:15;
			for(int j=0; j<3; j++) { //lu,lv,lw	
			for(int n=0; n<2; n++) { //x,y
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+j+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,100,0,430,40,-y,y,"","S"+is+xtxt0[j],"DC-"+det[i]+" "+xyz[n]), in);
					dg.addDataSet(f1, in); in++;  
				}
			}
			}
		 	this.getDataGroup().add(dg,i,st,k,run);
		}		
		return;	
		
        case 1:   
		float[] x1={0,5}; float[] x2={EB,35};
		f1 = new F1D("H_e_EC_resid_f1+"+run,"[a]",x1[0],x2[0]); 
		f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
		f2 = new F1D("H_e_EC_resid_f2+"+run,"[a]",x1[1],x2[1]); 
		f2.setParameter(0, 0f); f2.setLineColor(1); f2.setLineWidth(1);
		String[] xtxt1 = {" p_e"," #theta_e"}; 

		for(int i=0;i<3;i++) { //pcal,ecin,ecou
	        dg = new DataGroup(6,4); in=0; float y = i<1?5:15;
			for(int j=0; j<2; j++) { //p,the	
			for(int n=0; n<2; n++) { //x,y
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+j+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,x1[j],x2[j],40,-y,y,"","S"+is+xtxt1[j],"DC-"+det[i]+" "+xyz[n]), in);
					dg.addDataSet((j==0?f1:f2), in); in++;  
				}
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
					dg.addDataSet(makeH2(tab+"_1_",tag,60,5,35,100,n==0?-2:-10,n==0?12:27,"","S"+is+" #theta_e",lab2[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
			}
		}	
		}        	
        
    	break;
    	
    	case 3:
    	dg = new DataGroup(6,6);    	
        String lab3[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};    
    	f1 = new F1D("H_e_EC_rad3_f+"+run,"[a]",0.0,0.15); 
    	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
    	for(int i=0;i<3;i++) { //pcal,ecin,ecou
    		for(int n=0; n<2; n++) { //thet,phi
    			for(int is=1;is<7;is++){  //sector  	
    				tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,0,0.15,100,n==0?-2:-10,n==0?12:27,"","S"+is+" "+det[i]+" E",lab3[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
    			}
    		}	
    	}
    	
    	break;
    	
    	case 4:
    	dg = new DataGroup(6,6);    	
        String lab4[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};    
    	f1 = new F1D("H_e_EC_rad4_f+"+run,"[a]",0.0,EB); 
    	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
    	for(int i=0;i<3;i++) { //pcal,ecin,ecou
    		for(int n=0; n<2; n++) { //thet,phi
    			for(int is=1;is<7;is++){  //sector  	
    				tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,0,EB,100,n==0?-2:-10,n==0?12:27,"","S"+is+" p_e",lab4[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
    			}
    		}	
    	}
    	
    	break;
    	
    	case 5:
        dg = new DataGroup(3,2); GraphErrors g = null;
    	f1 = new F1D("H_e_EC_zero_f+"+run,"[a]",0,7); 
    	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
        for (int i=0; i<6; i++) { //pcal,ecin,ecou (dcres) pcal,ecin,ecou (egres) 
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
    
    public void createECpim(int st) {

    	String tab = "ECpim", tag = null;   	
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);

		F1D f1,f2,f3; 
		
    	switch (st) {   	
        case 0: 
        f1 = new F1D("ECposi_1+"+run,"[a]/x",0.4,3.0);      f1.setParameter(0,0.108f);  f1.setLineWidth(1);
        f2 = new F1D("H_pim_EC_resid0_f+"+run,"[a]",1,69);  f2.setParameter(0,0f);      f2.setLineWidth(1);
        f3 = new F1D("H_pim_EC_resid0_f+"+run,"[a]",0,1.6);f3.setParameter(0,0f);      f3.setLineWidth(1);
        dg = new DataGroup(6,4);
        for(int is=1;is<7;is++){
	        tag = is+"_"+st+"_"+k+"_"+run;
			dg.addDataSet(makeH2(tab+"_1_",tag,60,0,10,60,0.0,0.35,"","p (GeV)","E/P"),is-1);                dg.addDataSet(f1,is-1);
			dg.addDataSet(makeH2(tab+"_2_",tag,68,1,69,50,0.5,7,"","PCAL U STRIP","EC / PCAL"),is-1+6);
			dg.addDataSet(makeH2(tab+"_3_",tag,68,1,69,50,-4,4,     "","PCAL U STRIP","#chi PID"),is-1+12); dg.addDataSet(f2,is-1+12);
//			dg.addDataSet(makeH2(tab+"_4_",tag,70,0,EB/6,70,0,EB/5, "","PC (GeV)","EC (GeV)"),is-1+18);	
			dg.addDataSet(makeH2(tab+"_4_",tag,60,0,1.6,50,-5,5, "","Em (GeV)","#chi PID"),is-1+18);       dg.addDataSet(f3,is-1+18);	
        }		
		break;
    	
        case 1:   
		float[] x1={0,5}; float[] x2={EB,35};
		f1 = new F1D("H_pim_EC_resid_f1+"+run,"[a]",x1[0],x2[0]); 
		f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
		f2 = new F1D("H_pim_EC_resid_f2+"+run,"[a]",x1[1],x2[1]); 
		f2.setParameter(0, 0f); f2.setLineColor(1); f2.setLineWidth(1);
		String[] xtxt = {" p_pim"," #theta_pim"}; 

		for(int i=0;i<3;i++) { //pcal,ecin,ecou
	        dg = new DataGroup(6,4); int in=0; float y = i<1?5:15;
			for(int j=0; j<2; j++) { //p,the	
			for(int n=0; n<2; n++) { //x,y
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+j+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,x1[j],x2[j],40,-y,y,"","S"+is+xtxt[j],"DC-"+det[i]+xyz[n]), in);
					dg.addDataSet((j==0?f1:f2), in); in++;  
				}
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
    
    public void createECposi(int st) {
    	
    	String tab = "ECposi", tag = null;
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
    	
        F1D f1, f2, f3; 
        
    	switch (st) {
        case 0: 
        f1 = new F1D("ECposi_1+"+run,"[a]/x",0.4,0.9); f1.setParameter(0,0.108f); f1.setLineWidth(1);
        f2 = new F1D("ECposi_2+"+run,"[a]",1,69);      f2.setParameter(0,0f);     f2.setLineWidth(1);
        f3 = new F1D("ECposi_3+"+run,"[a]",0,0.75);    f3.setParameter(0,0f);     f3.setLineWidth(1);
        dg = new DataGroup(6,4);
        for(int is=1;is<7;is++){
	        tag = is+"_"+st+"_"+k+"_"+run;
			dg.addDataSet(makeH2(tab+"_1_",tag,60,0,3,60,0.12,0.35,"","p (GeV)","E/P"),is-1);                dg.addDataSet(f1,is-1);
			dg.addDataSet(makeH2(tab+"_2_",tag,68,1,69,50,0.12,0.35,"","PCAL U STRIP","E/P"),is-1+6);
			dg.addDataSet(makeH2(tab+"_3_",tag,68,1,69,50,-4,4,      "","PCAL U STRIP","#chi PID"),is-1+12); dg.addDataSet(f2,is-1+12);
//			dg.addDataSet(makeH2(tab+"_4_",tag,70,0,EB/6,70,0,EB/5, "","PC (GeV)","EC (GeV)"),is-1+18);	
			dg.addDataSet(makeH2(tab+"_4_",tag,60,0,0.75,50,-5,5, "","Em (GeV)","#chi PID"),is-1+18);        dg.addDataSet(f3,is-1+18);	
        } 
    	}
    	
    	this.getDataGroup().add(dg,0,st,k,run);     	
    	
    }
    
    public void createECtime(int st) {
    	
    	String tab = "ECtime", tag = null;
    	int run = getRunNumber(), in=0, k=getDetectorTabNames().indexOf(tab);
    	 
        F1D     f1=null;
		
    	switch (st) {
    	
        case 0:
        dg = new DataGroup(2,5);
        for(int i=0; i<2; i++) {
	        tag = in+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_1_",tag,200,-200,200,200,-200,200,det[i],"X (CM)","Y (CM)"),in); in++;         	
        } 
        for(int i=0; i<2; i++) {
	        tag = in+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_2_",tag,200,-200,200,200,-200,200,"","X (CM)","Y (CM)"),in); in++;  	
	        dg.addDataSet(makeH2(tab+"_3_",tag,200,-200,200,200,-200,200,"","X (CM)","Y (CM)"),in); in++;  	
	        dg.addDataSet(makeH2(tab+"_4_",tag,200,-200,200,200,-200,200,"","X (CM)","Y (CM)"),in); in++;  	
	        dg.addDataSet(makeH2(tab+"_5_",tag,200,-200,200,200,-200,200,"","X (CM)","Y (CM)"),in); in++;  	
	    }        
        break;
    	
    	case 2:
    		
        dg = new DataGroup(6,4);    	
        for(int i=0;i<2;i++) { //pcal,ecin
            for(int n=0; n<2; n++) { //phi_e,phi_ecal
            	for(int is=1;is<7;is++){  //sector  	
            		tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
        			dg.addDataSet(makeH2(tab+"_1_",tag,60,-60,60,60,5,20,"","S"+is+" #phi_"+(n==0?"e":det[i])+" (deg)","#theta_"+det[i]+" (deg)"), in);
        			in++; 
            	}
            }	
        }    		
        break;
        
    	case 3:
    		
        dg = new DataGroup(6,3);    	
        for(int i=0;i<3;i++) { //pcal,ecin,ecou
        	for(int n=0; n<1; n++) { 
        		for(int is=1;is<7;is++){  //sector  	
        			tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
    				dg.addDataSet(makeH2(tab+"_1_",tag,60,-6,6,60,0,0.5,"","S"+is+" "+det[i]+" t-tv-path/c",det[i]+" E"), in);
    				in++; 
        		}
        	}	
        }    		
    	break;
    	
    	case 4:
    	dg = new DataGroup(6,6);    	
        String lab4[] = new String[]{"e - #gamma #Delta#theta","e - #gamma #Delta#phi"};    
    	f1 = new F1D("H_e_ECt_4_f+"+run,"[a]",-6,6); 
    	f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
    	for(int i=0;i<3;i++) { //pcal,ecin,ecou
    		for(int n=0; n<2; n++) { //thet,phi
    			for(int is=1;is<7;is++){  //sector  	
    				tag = is+"_"+n+"_"+i+"_"+st+"_"+k+"_"+run;
    				float sca=1; if (i>0&&n==1) sca=2;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,-6,6,60,-5*sca,20*sca,"","S"+is+" t-tv-path/c ",lab4[n]+" "+det[i]), in);
					dg.addDataSet(f1, in); in++; 
    			}
    		}	
    	}
        
    	}
    	
    	this.getDataGroup().add(dg,0,st,k,run);    
    	
    } 
        
    public void createSCelec(int st) {
    	
    	String tab = "SCelec", tag = null;
    	int    run=getRunNumber(), k=getDetectorTabNames().indexOf(tab);
    	
		F1D     f1 = null, f2=null;    
		
    	switch (st) {
    	
        case 1:   
		String[] xtxt = {" p_e"," bar number",}; 

		for(int i=0;i<2;i++) { //p1a,p1b
			float[] x1={0,1}; float[] x2={EB,i==0?24:63};
			f1 = new F1D("H_e_SC_resid_f1+"+run,"[a]",x1[0],x2[0]); 
			f1.setParameter(0, 0f); f1.setLineColor(1); f1.setLineWidth(1);
			f2 = new F1D("H_e_SC_resid_f2+"+run,"[a]",x1[1],x2[1]); 
			f2.setParameter(0, 0f); f2.setLineColor(1); f2.setLineWidth(1);
	        dg = new DataGroup(6,4); int in=0;  
			for(int j=0; j<2; j++) { //p,the	
			for(int n=0; n<2; n++) { //x,y
				for(int is=1;is<7;is++){  //sector  	
					tag = is+"_"+n+"_"+j+"_"+i+"_"+st+"_"+k+"_"+run;
					float y1 = i==0?n==0?15:5:5; int nx=j==0?60:i==0?23:62;
					dg.addDataSet(makeH2(tab+"_1_",tag,nx,x1[j],x2[j],40,-y1,y1,"","S"+is+xtxt[j],"DC-"+scdet[i]+xyz[n]), in);
					dg.addDataSet((j==0?f1:f2), in); in++;  
				}
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
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
        F1D f = new F1D("f_prot"+tag,"1/(1+[a]^2/x^2)^0.5", 0.5,4.1); f.setParameter(0,0.93827); f.setLineColor(1); f.setLineStyle(1);   
    	
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
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,4.1,100,0.3,1.05,     "","p (GeV)","#beta"),n);dg.addDataSet(f,n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,4.1,100,0.3,1.05,     "","p (GeV)","#beta_pcal"),n);dg.addDataSet(f,n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0,4.1,100,0.3,1.05,     "","p (GeV)","#beta_ecin"),n);dg.addDataSet(f,n);n++;
    	}    	
    	this.getDataGroup().add(dg,0,st,k,run);      
   	
    }
    
    public void createECpbar(int st) {
    	
    	String tab = "ECpbar", tag = null;   	
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
    	
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
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
    	
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
    	break;
    	case 2:
        dg = new DataGroup(6,3);  
        for(int is=1;is<7;is++){ 
        	tag = is+"_"+st+"_"+k+"_"+run;        	
           	dg.addDataSet(makeH1(tab+"_0_",tag, 50,0,100," ","PCAL E (MeV)"),is-1); 
           	dg.addDataSet(makeH1(tab+"_1_",tag, 50,0,100," ","ECIN E (MeV)"),is-1+6); 
           	dg.addDataSet(makeH1(tab+"_2_",tag, 50,0,100," ","ECOU E (MeV)"),is-1+12); 
        }
        break;
    	case 3:
        dg = new DataGroup(5,3);
        for(int iv=0; iv<5; iv++) {
	        tag = iv+"_"+st+"_"+k+"_"+run;
	        dg.addDataSet(makeH2(tab+"_1_",tag,90,-400,400,90,-400,400,det5[iv],"X (CM)","Y(CM)"),iv);        	
	        dg.addDataSet(makeH2(tab+"_2_",tag,90,-400,400,90,-400,400,det5[iv],"X (CM)","Y(CM)"),iv+5);        	
	        dg.addDataSet(makeH2(tab+"_3_",tag,90,-400,400,90,-400,400,det5[iv],"X (CM)","Y(CM)"),iv+10);        	
        }  
        break;
    	case 4:
    		dg = new DataGroup(6,3);int inn=0;
			String[]   fdet = new String[]{"FTOF1B ","FTOF1A ","PCAL "};
			for(int i=0; i<3; i++) { //P1A, P1B, PCAL
				for(int is=1;is<7;is++){  //sector 
					tag = is+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,0.5,0.9,40,-50,50,"","S"+is+" p_pi+",fdet[i]+" Y(CM)"),inn); inn++;
				}
			}
        break;
    	case 5:
    		dg = new DataGroup(6,3);inn=0;
			fdet = new String[]{"FTOF1B ","FTOF1A ","PCAL "};
			for(int i=0; i<3; i++) { //P1A, P1B, PCAL
				for(int is=1;is<7;is++){  //sector 
					tag = is+"_"+i+"_"+st+"_"+k+"_"+run;
					dg.addDataSet(makeH2(tab+"_1_",tag,60,25,35,40,-50,50,"","S"+is+" #theta_pi+",fdet[i]+" Y(CM)"),inn); inn++;
				}
			} 
    	}
     		
    	this.getDataGroup().add(dg,0,st,k,run);      
   	
    }
    
    public void createECpi0(int st) {
    	
    	String tab = "ECpi0", tag = null;   	
    	int run = getRunNumber(), n=0, k=getDetectorTabNames().indexOf(tab); 
    	
    	switch (st) {
    	
    	case 0:
    	dg = new DataGroup(5,2);  
    	tag = st+"_"+k+"_"+run;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0.9,1.25,100,0,8,   "photon from tagged pizero","#gamma1 #beta","ECAL energy (GeV)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,100,0.9,1.25,100,0,8,   "photon from tagged pizero","#gamma2 #beta","ECAL energy (GeV)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,100,0.9,1.25,           "photon from tagged pizero","#gamma1 #beta"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,100,0.9,1.25,           "photon from tagged pizero","#gamma2 #beta"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,100,0.9,1.25,           "all photons","#gamma #beta"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,180,0,180,              "N#gamma(GEN)=2","#gamma #gamma #phi (deg)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,180,0,180,              "N#gamma(REC)=2","#gamma #gamma #phi (deg)"),n);n++;
    	break;   
        case 1:		
        dg = new DataGroup(4,3);  
        tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0.0,5.5,50,0,100,   "pizero","p_mm (GeV)","#theta_mm (^o)"),n);n++;
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
    
    public void createECeta(int st) {
    	
    	String tab = "ECeta", tag = null;   	
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab); 
    	
    	switch (st) {
    	   
        case 1:		
        dg = new DataGroup(4,3); int n=0;
        tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0.0,5.5,50,0,100,      "eta","p_mm (GeV)","#theta_mm (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.5,0.5,50,-0.5,0.5,    " ","cx_mm - cx_ecal","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                " ","cx_mm - cx_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                " ","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,45,50,0,4,             " ","Opening Angle (^o)","E1*E2 (GeV^2)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,5.5,50,0,5.5,          " ","p_mm (GeV)","p_ecal (GeV)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,60,0,0.8,                   " ","IVM (GeV)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,60,0,0.8,                   " ","IVM (GeV)"),n);n++;
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
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
        F1D f1 = new F1D("neut","1/(1+[a]^2/x^2)^0.5", 0.4,4.0); 
        f1.setParameter(0,0.93957); f1.setLineColor(0); f1.setLineStyle(1); f1.setLineWidth(2);  
        
    	switch (st) {
    	case 0:
    	dg = new DataGroup(6,4); 
        for(int is=1;is<7;is++){ 
        	tag = is+"_"+st+"_"+k+"_"+run;        	
        	dg.addDataSet(makeH2(tab+"_0_",tag,50,0,2.5,50,0.,100,"S"+is,"p_mm (GeV)","PCAL #pi^+-n (cm)"),is-1); 
        	dg.addDataSet(makeH2(tab+"_1_",tag,50,0,2.5,50,0.,100,"",    "p_mm (GeV)","ECIN #pi^+-n (cm)"),is-1+6); 
        	dg.addDataSet(makeH2(tab+"_2_",tag,50,0,2.5,50,0.,100,"",    "p_mm (GeV)","ECOU #pi^+-n (cm)"),is-1+12); 
        	dg.addDataSet(makeH2(tab+"_3_",tag,50,-0.5,2.0,50,-0.5,0.5,"",    "Mass^2 (GeV^2)","cx_mm - cx_ecal"),is-1+18); 
        }
        break;
        case 1:		
        dg = new DataGroup(4,3); int n=0;
        tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,80,0,4,50,0,40,             "neutron","p_mm (GeV)","#theta_mm (^o)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.5,0.5,50,-0.5,0.5,    " ","cx_mm - cx_ecal","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                "cycut.M2cut ","cx_mm - cx_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.5,0.5,                "cxcut.M2cut ","cy_mm - cy_ecal"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,80,0,4,50,0.3,1.2,  "no electron","p_mm (GeV)","#beta_ecal"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,80,0,4,50,0.3,1.2,            " ","p_mm (GeV)","#beta_ecal"),n);n++;
    	dg.addDataSet(f1, n-1);
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.2,2.0,          "no electron","Mass^2 (GeV^2)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,50,-0.2,2.0,            "cx,cy cut","Mass^2 (GeV^2)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.6,0.0,50,-0.3,0.3,         " ","cx_mm","cy_mm"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,-0.6,0.0,50,-0.3,0.3,     " ","cx_ecal","cy_ecal"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,80,0,4,                   " ","p_mm (GeV)"),n  );n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,80,0,4,                   " ","p_mm (GeV)"),n-1);n++;  
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,80,0,4,                   " ","p_mm (GeV)"),n-2);n++;  
    	((H1F)dg.getData(10).get(1)).setFillColor(4);
    	((H1F)dg.getData(10).get(2)).setFillColor(2);
    	}
    	this.getDataGroup().add(dg,0,st,k,run);      
    	
    }
   	
    public void createECphot(int st) {
    	
    	String tab = "ECphot", tag = null;   	
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
    	
    	switch (st) {
        case 0:       
        for(int i=0; i<3; i++) { //pcal,ecin,ecou
            dg = new DataGroup(6,4);  int n=0;
        	for(int id=0; id<4; id++) {
                for(int is=1; is<7; is++) {
                	tag = is+"_"+id+"+"+i+"_"+st+"_"+k+"_"+run;
                	int col = id==0||id==2?43:44;
                	String txt = det[i]+"S"+is+(id<2?" ELEC-"+(id==0?"n":"#gamma"):" NOELEC-"+(id==2?"n":"#gamma"));
                	dg.addDataSet(makeH1(tab+"_1_",tag,15,0.5,15.5,txt," ",1,0),n); 
                	dg.addDataSet(makeH1(tab+"_2_",tag,15,0.5,15.5,txt," ",col,col),n); n++;                	
                }
            }
        	this.getDataGroup().add(dg,i,st,k,run);    
        } 
		return;	   
        case 1:    	
        dg = new DataGroup(5,3); int n=0;
        tag = st+"_"+k+"_"+run;        	
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,500,60,5,25,   "electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,300,50,5,35,   "electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
       	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,300,50,5,35,"no electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,25, 50,5,25,   "electron","#theta_elec (deg)","#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,25, 50,5,25,"no electron","#theta_elec (deg)","#theta_ecal (deg)"),n);n++;
     	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,300,50,0,35,"no electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
       	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,50,0,300,50,0,35,"no electron","esum_ecal (MeV)",  "#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,6,1,7,               "electron","Multiplicity"),n)                         ;n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,6,1,7,            "no electron","Multiplicity"),n)                         ;n++;
    	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,6,1,7,         "tagged neutron","Multiplicity"),n)                         ;n++;
       	dg.addDataSet(makeH1(tab+"_"+n+"_",tag,6,1,7,   "tagged neutron < 1.2","Multiplicity"),n)                         ;n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,80,-180,180,80,5,22,"electron","#phi_ecal (deg)",  "#theta_ecal (deg)"),n);n++;
    	dg.addDataSet(makeH2(tab+"_"+n+"_",tag,80,-180,180,80,5,22,"electron","#phi_e (deg)"   ,  "#theta_ecal (deg)"),n);n++;
    	}
    	this.getDataGroup().add(dg,0,st,k,run);  
    	
    }
    
    public void createECtrig(int st) {
    	String tab = "ECtrig", tag = null;   	
    	int run = getRunNumber(), k=getDetectorTabNames().indexOf(tab);
    	
    	switch (st) {
    	   
        case 0:
        dg = new DataGroup(3,3); int n=0;
        tag = st+"_"+k+"_"+run;
        dg.addDataSet(makeH1(tab+"_"+n+"_",tag,32,-0.5,31.5,"All Events","Trigger Bits"),n);n++;
        dg.addDataSet(makeH1(tab+"_"+n+"_",tag,32,-0.5,31.5,"REC::event,particle","Trigger Bits"),n);n++;
        dg.addDataSet(makeH1(tab+"_"+n+"_",tag,32,-0.5,31.5,"1 Electron","Trigger Bits"),n);n++;
        dg.addDataSet(makeH1(tab+"_"+n+"_",tag,32,-0.5,31.5,"No Electron","Trigger Bits"),n);n++;
        dg.addDataSet(makeH1(tab+"_"+n+"_",tag,32,-0.5,31.5,"#pi-","Trigger Bits"),n);n++;
        dg.addDataSet(makeH1(tab+"_"+n+"_",tag,32,-0.5,31.5,"#pi- No Electron","Trigger Bits"),n);n++;
    	}
    	this.getDataGroup().add(dg, 0,st,k,run);
    }
	
    // MAKE
    
    public boolean makeELEC(){
    	
    	if(nE!=1) return false; //e_ecal.size()!=1;
    	
    	boolean[] good_fiduc = {false,false,false};
        e_ecal_esum = 0f; e_ecal_pcsum=0; e_ecal_ecsum=0; e_ecal_rat=0;
        
        elec_ecal_resid.clear();
        elec_ftof_resid.clear();  
        
        Particle epart = e_ecal.get(0);
                
        short   status = (short) epart.getProperty("status");
        float   chipid = (float) epart.getProperty("chi2pid");
        boolean   inDC = (status>=2000 && status<3000);
        
        if(!inDC) return false;
        
        e_mom = (float) epart.p();      
        e_vz  = (float) epart.vz();
        
        List<Particle> elecECAL = ev.getECAL((int)epart.getProperty("pindex"));
    	List<Particle> elecHTCC = ev.getHTCC((int)epart.getProperty("pindex"));
    	
    	for (int i=0; i<3; i++) {iU[i]=0; iV[i]=0; iW[i]=0; lU[i]=0; lV[i]=0; lW[i]=0;}

//    	Boolean good_e = elecECAL.size()==1 && elecECAL.get(0).getProperty("sector")==ev.getElecTriggerSector();
//    	if(!good_e) return false;
    	  
    	for (Particle p : elecECAL) {    		
            e_sect   =   (int) p.getProperty("sector");
    		e_x      = (float) p.getProperty("x");
    		e_y      = (float) p.getProperty("y");
            e_cz     = p.hasProperty("cz")?(float) p.getProperty("cz"):0;
            e_t      = (float) p.getProperty("time");
            e_pat    = (float) p.getProperty("path");
    		float en = (float) p.getProperty("energy");   		
    		int  lay = getDet((int) p.getProperty("layer"));    		    		
    		
    		lU[lay] = p.hasProperty("lu")?(float)p.getProperty("lu"):0;
    		lV[lay] = p.hasProperty("lv")?(float)p.getProperty("lv"):0;
    		lW[lay] = p.hasProperty("lw")?(float)p.getProperty("lw"):0; 
    		
    		iU[lay] = p.hasProperty("iu")?(int)p.getProperty("iu"):(int)(lU[lay]*0.891/(lay==0?4.5:10));
    		iV[lay] = p.hasProperty("iv")?(int)p.getProperty("iv"):(int)(lV[lay]*0.891/(lay==0?4.5:10));
    		iW[lay] = p.hasProperty("iw")?(int)p.getProperty("iw"):(int)(lW[lay]*0.891/(lay==0?4.5:10));
    		    		
    		Point3D xyz = getResidual(p);	    

    		elec_ecal_resid.add((float)xyz.x(),e_sect,0,lay);
    		elec_ecal_resid.add((float)xyz.y(),e_sect,1,lay);
//    		elec_ecal_resid.add((float)xyz.z(),e_sect,2,ind);

    		elec_ecal_resid.add((float)iU[lay],e_sect,3,lay);   		
    		elec_ecal_resid.add(chipid,e_sect,4,lay);
    		
    		if(lay>-1) e_ecal_esum  += en;
    		if(lay==0) e_ecal_pcsum  = en;
    		if(lay>0)  e_ecal_ecsum += en;

    		if (lay==0) good_fiduc[0] = iU[0]>2 && iV[0]<62 && iW[0]<62; //PCAL
    	    if (lay==1) good_fiduc[1] = iU[1]>2 && iV[1]<36 && iW[1]<36; //ECIN
    	    if (lay==2) good_fiduc[2] = iU[2]>2 && iV[2]<36 && iW[2]<36; //ECOU	
   	    }
    	
    	e_ecal_rat = e_ecal_pcsum/e_ecal_esum;
    	e_ecal_sf = e_ecal_esum/e_mom/1000f;
    	
    	for (Particle p : elecHTCC) {
            htcc_sect   =   (int) p.getProperty("sector");
            htcc_x      =   (float) p.getProperty("x");
            htcc_y      =   (float) p.getProperty("y");
            htcc_npe    =   (float) p.getProperty("nphe");    		
    	}
          
//    	if (fiduCuts && !((good_fiduc1)||(good_fiduc1&&good_fiduc2)||(good_fiduc1&&good_fiduc2&&good_fiduc3))) return false;
    	
    	if (fiduCuts && !(good_fiduc[0]&&good_fiduc[1]&&good_fiduc[2])) return false;
    	
        if(Math.abs(e_vz+3)<12 && e_mom>0.1){
            e_the  = (float) Math.toDegrees(epart.theta());
            e_phi  = (float) Math.toDegrees(epart.phi()); ephi = (e_phi>-30?e_phi:360+e_phi)-(e_sect-1)*60;
            e_vx   = (float) epart.vx(); 
            e_vy   = (float) epart.vy();
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
    
    public boolean makePOSI(){
    	
    	boolean good_fiduc1 = false, good_fiduc2 = false, good_fiduc3 = false; 
        p_ecal_esum = 0f; p_ecal_pcsum=0; p_ecal_ecsum=0; p_sect=0;
        posi_ecal_resid.clear();
        posi_ftof_resid.clear();  
        
        Particle ppart = p_ecal.get(0);
                
        short   status = (short) ppart.getProperty("status");
        float   chipid = (float) ppart.getProperty("chi2pid");
        boolean   inDC = (status>=2000 && status<3000);
        
        if(!inDC) return false;
        
        p_mom = (float) ppart.p();      
        p_vz  = (float) ppart.vz();
        
        List<Particle> posiECAL = ev.getECAL((int)ppart.getProperty("pindex"));
    	List<Particle> posiFTOF = ev.getFTOF((int)ppart.getProperty("pindex"));
    	
    	for (Particle p : posiFTOF) {
            p_sect    = (int) p.getProperty("sector");   			
	    	int scind = (int) p.getProperty("layer");
		    Point3D xyz = getResidual(p);
	        posi_ftof_resid.add((float)xyz.x(),p_sect,0,scind-1);
	    	posi_ftof_resid.add((float)xyz.y(),p_sect,1,scind-1);
	    	posi_ftof_resid.add((float)xyz.z(),p_sect,2,scind-1); 		 		
    	}
    	  
    	for (Particle p : posiECAL) {    		
            p_sect   =   (int) p.getProperty("sector");
    		p_x      = (float) p.getProperty("x");
    		p_y      = (float) p.getProperty("y");
            p_cz     = p.hasProperty("cz")?(float) p.getProperty("cz"):0;
    		float en = (float) p.getProperty("energy");   		
    		int  ind = getDet((int) p.getProperty("layer"));
    		
    		piU = p.hasProperty("iu")?(int)p.getProperty("iu"):0;
    		piV = p.hasProperty("iv")?(int)p.getProperty("iv"):0;
    		piW = p.hasProperty("iw")?(int)p.getProperty("iw"):0; 
    		
    		Point3D xyz = getResidual(p);	        
    		posi_ecal_resid.add((float)xyz.x(),p_sect,0,ind);
    		posi_ecal_resid.add((float)xyz.y(),p_sect,1,ind);
    		posi_ecal_resid.add((float)xyz.z(),p_sect,2,ind);
    		posi_ecal_resid.add(p.hasProperty("iu")?(float)p.getProperty("iu"):0,p_sect,3,ind);
    		posi_ecal_resid.add(chipid,p_sect,4,ind);
    		
    		if(ind>-1) p_ecal_esum  += en;
    		if(ind==0) p_ecal_pcsum  = en;
    		if(ind>0)  p_ecal_ecsum += en;
    		if (ind==0) good_fiduc1 = piU>2 && piV<62 && piW<62; //PCAL
    	    if (ind==1) good_fiduc2 = piU>2 && piV<36 && piW<36; //ECIN
    	    if (ind==2) good_fiduc3 = piU>2 && piV<36 && piW<36; //ECOU	
   	    }
          
//    	if (fiduCuts && !((good_fiduc1)||(good_fiduc1&&good_fiduc2)||(good_fiduc1&&good_fiduc2&&good_fiduc3))) return false;
    	
    	if (fiduCuts && !(good_fiduc1&&good_fiduc2&&good_fiduc3)) return false;
    	
        if(Math.abs(p_vz+3)<12 && p_mom>0.1){
            p_the  = (float) Math.toDegrees(ppart.theta());
            p_phi  = (float) Math.toDegrees(ppart.phi());
            p_vx   = (float) ppart.vx(); 
            p_vy   = (float) ppart.vy();
            return true;
         }  
        
        return false;

    }

    public boolean makeNM() {
    	
    	if(!goodPHOT) return false;
   	
    	int   G1_part_ind=-1, G2_part_ind=-1, n=0;
    	float G1_mom=0,G2_mom=0;
        
    	nm_ecal.clear();
    	
        for (Particle p : phot_ecal) {
			if(p.p()>0.2 && Math.toDegrees(p.theta())>6 && G1_mom < p.p()){
				G1_mom = (float) p.p();
				G1_part_ind = n;				    			
			}
			if(G1_part_ind>-1 && n!=G1_part_ind && p.p()>0.2 && Math.toDegrees(p.theta())>6 && G2_mom < p.p() && p.p() < G1_mom){
				G2_mom = (float) p.p();
				G2_part_ind = n;
			}
			n++;			
        }
        
        if(!(G1_part_ind>-1 && G2_part_ind>-1)) return false;
        
		NeutralMeson nm = new NeutralMeson(taggedPI0);
		nm.addPhoton(phot_ecal.get(G1_part_ind)); //System.out.println("makeNM "+nm.toString(0));			
		nm.addPhoton(phot_ecal.get(G2_part_ind)); //System.out.println("makeNM "+nm.toString(1));				
		nm_ecal.add(nm);
		
		return true;
       
    }
    
    public class NeutralMeson {
    	
    	public float mass;
    	public float e;
    	public float mom;
    	public float the;
    	public float phi;
    	public float opa;
    	public float X;
    	
    	public  LorentzVector VPI0;
    	public  LorentzVector VG1;
    	public  LorentzVector VG2;
    	
    	private boolean tag;
    	
    	public List<Particle>  plist = new ArrayList<Particle>();
    	
    	public NeutralMeson(boolean val) {
    		this.tag = val;
    	}
    	
    	public void addPhoton(Particle p) {
    		plist.add(p);
    	}
    	
    	public Particle getPhoton(int n) {
    		return plist.get(n);    		
    	}
    	
    	public int getPINDEX(int n) {
    		return (int) getPhoton(n).getProperty("pindex");
    	}
    	
    	public int getPhotonSector(int n) {
    		return (int) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("sector");
    	}
    	
    	public int getPhotonLayer(int n) {
    		return (int) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("layer");
    	}
    	
    	public float getPhotonEnergy(int n) {
    		return (float) ev.getECAL((int)getPhoton(n).getProperty("pindex")).get(0).getProperty("energy");
    	}
    	
    	public float getParticleEnergy(int n) {
    		return (float) ev.part.get(getPINDEX(n)).p();
    	}
    	
    	public float getBeta(int n) {
    		return (float) getPhoton(n).getProperty("beta");    		
    	}
    	
    	public float getEnergy(int n) {
    		return (float) getPhoton(n).p();    		
    	}  
    	
    	public Boolean filter(boolean val) {
    		if (val) return  sameSector();
    		return !sameSector();
    	}
    	
    	public boolean sameSector() {
    		return getPhotonSector(0)==getPhotonSector(1);		
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
			this.mass = (float)VPI0.mass();
			this.e    = (float)VPI0.e();
			this.mom  = (float)VPI0.p();
			this.the  = (float)Math.toDegrees(VPI0.theta());
			this.phi  = (float)Math.toDegrees(VPI0.phi());
			this.opa  = (float)Vangle(VG1.vect(),VG2.vect());
			this.X    = (float)((VG1.e()-VG2.e())/(VG1.e()+VG2.e()));
			return this.mass > 0.08 && filter(tag);			
    	}
    	
    	public String toString(int n) {
    		return getPINDEX(n)+" "+getPhotonSector(n)+" "+getPhotonLayer(n)+" "+getPhotonEnergy(n)+" "+getParticleEnergy(n);
    	}
    	    	    	
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
		DataGroup dg1 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECprot"),run);
		DataGroup dg0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECprot"),run);
		if(nPROT==1) { 			
	        List<Particle> protECAL = ev.getECAL((int)prot_ecal.get(0).getProperty("pindex"));
			prot_mom  = (float) prot_ecal.get(0).p();
            prot_the  = (float) Math.toDegrees(prot_ecal.get(0).theta());
            prot_phi  = (float) Math.toDegrees(prot_ecal.get(0).phi());
            prot_vz   = (float) prot_ecal.get(0).vz();
        	prot_beta = (float) prot_ecal.get(0).getProperty("beta");
        	for (Particle p : protECAL) {
        		int mult = ev.getECALMULT((int)p.getProperty("sector"),(int)p.getProperty("layer")); //count number of clusters in sector,layer
        		prot_beta_pcal = (mult==1 && p.getProperty("layer")==1) ? (float) p.getProperty("beta"):0;
        		prot_beta_ecin = (mult==1 && p.getProperty("layer")==4) ? (float) p.getProperty("beta"):0;
        		((H2F)dg1.getData(10).get(0)).fill(prot_mom,prot_beta_pcal);  	
        		((H2F)dg1.getData(11).get(0)).fill(prot_mom,prot_beta_ecin);  	
        	}
         	Vprot     = prot_ecal.get(0).vector();	
			ep_dPhi   = prot_phi - e_phi + 180f;
			while(ep_dPhi> 180f)ep_dPhi -= 360f;
			while(ep_dPhi<-180f)ep_dPhi += 360f;
			LorentzVector VmissN = new LorentzVector(0,0,0,0);
			VmissN.add(VT);
			VmissN.add(VB);
			VmissN.sub(Ve);
			VmissN.sub(Vprot);
			ep_MM = (float)VmissN.mass2();
			((H2F) dg0.getData(e_sect-1  ).get(0)).fill(ep_MM,e_W);  
			((H1F) dg0.getData(e_sect-1+6).get(0)).fill(ep_MM);      
			nm_mom=-1f;nm_the=-1f;nm_phi=-1f;
			taggedPI0 = ep_MM<0.1; taggedETA = ep_MM>0.24 && ep_MM<0.36;
			if(taggedPI0 || taggedETA) {
//		        List<Particle> protECAL = ev.getECAL((int)prot_ecal.get(0).getProperty("pindex"));
//				System.out.println("PROT: "+protECAL.get(0).getProperty("x")+" "+protECAL.get(0).getProperty("y"));
				nm_mom=(float)VmissN.p();
				nm_the=(float)Math.toDegrees(VmissN.theta());
				nm_phi=(float)Math.toDegrees(VmissN.phi());
				nm_cx =(float)VmissN.px()/nm_mom;
				nm_cy =(float)VmissN.py()/nm_mom;
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
		DataGroup dg0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECpip"),run);
					
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
		((H2F) dg0.getData(e_sect-1).get(0)).fill(epip_MM,e_W);    
		((H1F) dg0.getData(e_sect-1+6).get(0)).fill(epip_MM);
			
		neut_mom=-1f;neut_the=-1f;neut_phi=-1f;
			
		if(epip_MM<1.0 && pip_ecal.size()==1) {
			neut_mom=(float)VmissN.p();
			neut_the=(float)Math.toDegrees(VmissN.theta());
			neut_phi=(float)Math.toDegrees(VmissN.phi());
			return true;
		}
		
		return false;
	}
    
	public void debug() {
	    IndexGenerator ig = new IndexGenerator();
	    for (Map.Entry<Long,List<Particle>>  entry : ev.partmap.getMap().entrySet()){
	           long hash = entry.getKey();
	           int pid = ig.getIndex(hash, 0);
	           int sec = ig.getIndex(hash, 1);
	           for (Particle pp : ev.partmap.getItem(pid,sec)) {	        	   
	               System.out.println(pid+" "+sec+" "+(int)pp.getProperty("layer")
	                                             +" "+(int)pp.getProperty("status")
	                                             +" "+     pp.getProperty("energy")); 
	           }
	    }		
	} 

	// FILL
	
	public void fillECtrig(int i) {
		int run = getRunNumber();
		int k = getDetectorTabNames().indexOf("ECtrig");
		DataGroup dg0 = this.getDataGroup().getItem(0,0,k,run);
		for (int j=0; j<32; j++) {
			if(isTrigBitSet(j)) ((H1F) dg0.getData(i).get(0)).fill(j);
		}
	}
	
	public void fillECkin() {
		int run = getRunNumber();
		int k = getDetectorTabNames().indexOf("ECkin");
		DataGroup dg0 = this.getDataGroup().getItem(0,0,k,run);				
		DataGroup dg1 = this.getDataGroup().getItem(0,1,k,run);						
		DataGroup dg3 = this.getDataGroup().getItem(0,3,k,run);	
		
		((H2F) dg0.getData(e_sect-1   ).get(0)).fill(e_the,e_mom);
		((H2F) dg0.getData(e_sect-1+ 6).get(0)).fill(e_W,e_mom);
		((H2F) dg0.getData(e_sect-1+12).get(0)).fill(e_W,e_the);
		if (e_the>6) ((H2F) dg0.getData(e_sect-1+18).get(0)).fill(e_W,ephi);
		if(lV[0]>19 && lW[0]>19) ((H2F) dg1.getData(e_sect-1   ).get(0)).fill(lU[0],e_ecal_sf);
		if(lU[0]>60)             ((H2F) dg1.getData(e_sect-1+ 6).get(0)).fill(lV[0],e_ecal_sf);
		if(lU[0]>60)             ((H2F) dg1.getData(e_sect-1+12).get(0)).fill(lW[0],e_ecal_sf);
		((H2F) dg3.getData(e_sect-1   ).get(0)).fill(ephi,lU[0]);
		((H2F) dg3.getData(e_sect-1+ 6).get(0)).fill(ephi,lV[0]);
		((H2F) dg3.getData(e_sect-1+12).get(0)).fill(ephi,lW[0]);
		
//		ECkin = this.getDataGroup().getItem(0,2,k,run);				
//		((H2F) ECkin.getData(0).get(0)).fill(-e_x,e_y);
	}
	
	public void fillECelec() {
		int run = getRunNumber();
		int    k = getDetectorTabNames().indexOf("ECelec");//
		int   k1 = getDetectorTabNames().indexOf("ECkin");
		
		IndexGenerator ig = new IndexGenerator();
		
		DataGroup dg2 = this.getDataGroup().getItem(0,2,k,run);				
		DataGroup dg3 = this.getDataGroup().getItem(0,3,k,run);				
		DataGroup dg4 = this.getDataGroup().getItem(0,4,k,run);	
		
		DataGroup dg02 = this.getDataGroup().getItem(0,2,k1,run);				
		DataGroup dg04 = this.getDataGroup().getItem(0,4,k1,run);				
		DataGroup dg06 = this.getDataGroup().getItem(0,6,k1,run);	
		
		// Forward tracking residuals
		
		((H2F) dg02.getData(e_sect-1).get(0)).fill(e_mom,e_ecal_sf);
//		((H2F) dg0.getData(e_sect-1+18).get(0)).fill(e_ecal_pcsum/1000f,e_ecal_ecsum/1000f);
//		((H2F) dg0.getData(e_sect-1+18).get(0)).fill(e_ecal_esum/1000f,elec_ecal_resid.getItem(e_sect,4,0));

		for (Map.Entry<Long,Float>  entry : elec_ecal_resid.getMap().entrySet()){
			long hash = entry.getKey(); float val = entry.getValue();
			int is = ig.getIndex(hash, 0); int ic = ig.getIndex(hash, 1); int il = ig.getIndex(hash, 2);
			DataGroup dg0 = this.getDataGroup().getItem(il,0,k,run);				
			DataGroup dg1 = this.getDataGroup().getItem(il,1,k,run);	
			DataGroup dg05 = this.getDataGroup().getItem(il,5,k1,run);	

			if(ic==0) {((H2F)dg0.getData(is-1  ).get(0)).fill(lU[il],val); ((H2F)dg0.getData(is-1+12).get(0)).fill(lV[il],val); ((H2F)dg0.getData(is-1+24).get(0)).fill(lW[il],val);}
			if(ic==1) {((H2F)dg0.getData(is-1+6).get(0)).fill(lU[il],val); ((H2F)dg0.getData(is-1+18).get(0)).fill(lV[il],val); ((H2F)dg0.getData(is-1+30).get(0)).fill(lW[il],val);}
			
			if(ic==0) {((H2F)dg1.getData(is-1  ).get(0)).fill(e_mom,val); ((H2F)dg1.getData(is-1+12).get(0)).fill(e_the,val);}
			if(ic==1) {((H2F)dg1.getData(is-1+6).get(0)).fill(e_mom,val); ((H2F)dg1.getData(is-1+18).get(0)).fill(e_the,val);}
			
			if(ic==0 && il==0 && lU[0]>0 && lU[1]>0 && e_mom>7 && e_mom<8) {
				((H2F)dg05.getData(is-1  ).get(0)).fill(lU[0],lU[0]-lU[1]); ((H2F)dg05.getData(is-1+6).get(0)).fill(lV[0],lV[0]-412+lV[1]); ((H2F)dg05.getData(is-1+12).get(0)).fill(lW[0],lW[0]-375+lW[1]);	
			}
			
			if(ic==0 && il<2) {
				((H2F)dg06.getData(is-1+18*il).get(0)).fill(lU[il],iU[il]); ((H2F)dg06.getData(is-1+6+18*il).get(0)).fill(lV[il],iV[il]); ((H2F)dg06.getData(is-1+12+18*il).get(0)).fill(lW[il],iW[il]);		
			}
						
			if(e_mom>0) {
			if(ic==3) ((H2F) dg02.getData(e_sect-1+ 6).get(0)).fill(elec_ecal_resid.getItem(e_sect,3,0),e_ecal_sf);
			if(ic==3) ((H2F) dg02.getData(e_sect-1+12).get(0)).fill(elec_ecal_resid.getItem(e_sect,3,0),elec_ecal_resid.getItem(e_sect,4,0));
			if(ic==3) ((H2F) dg02.getData(e_sect-1+18).get(0)).fill(elec_ecal_resid.getItem(e_sect,3,0),e_ecal_rat);
			}	
		}
		
		counter[e_sect-1][0]++;
		
		// Radiative photon residuals
		
    	for (Map.Entry<Long,List<Particle>>  entry : ecphot_filt.getMap().entrySet()){
			int is = ig.getIndex(entry.getKey(), 0);
            	for (Particle ec : entry.getValue()) {				
            		int         il = getDet((int) ec.getProperty("layer"));
            		float ecal_nrg = (float) (ec.getProperty("energy")/1e3);
            		float ecal_the = (float) Math.toDegrees(ec.theta());
            		float ecal_phi = (float) Math.toDegrees(ec.phi());
            		float  ecalphi = (ecal_phi>-30?ecal_phi:360+ecal_phi)-(e_sect-1)*60;
            		float    thdif = (float)(e_the-ecal_the)*Math.signum(ev.tpol);           		
                    float    phdif = (float)(ephi-ecalphi)*Math.signum(ev.spol);
                    
            		if(dropEsect?is!=e_sect:is==e_sect) {                    
            			((H2F)dg2.getData(is-1+   0+il*12).get(0)).fill(e_the,thdif); 
            			((H2F)dg2.getData(is-1+   6+il*12).get(0)).fill(e_the,phdif);
            			((H2F)dg3.getData(is-1+   0+il*12).get(0)).fill(ecal_nrg,thdif);           		
            			((H2F)dg3.getData(is-1+   6+il*12).get(0)).fill(ecal_nrg,phdif);     
            			((H2F)dg4.getData(is-1+   0+il*12).get(0)).fill(e_mom,thdif); 
            			((H2F)dg4.getData(is-1+   6+il*12).get(0)).fill(e_mom,phdif);   
            			((H2F)dg04.getData(is-1+  0+il*12).get(0)).fill(ephi,thdif); 
            			((H2F)dg04.getData(is-1+  6+il*12).get(0)).fill(ephi,phdif); 

            			if(il==0 && thdif>-1 && thdif<1)   counter[e_sect-1][1]++;
            			if(il==0 && phdif>-2 && phdif<1.5) counter[e_sect-1][2]++;
            		}

            	}           	
    	}
				
	}
	
	public void fillECposi() {
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("ECposi");
		
		if(p_sect==0) return;
		
		IndexGenerator ig = new IndexGenerator();
		
		DataGroup dg0 = this.getDataGroup().getItem(0,0,k,run);				
//		DataGroup dg1 = this.getDataGroup().getItem(0,1,k,run);				
//		DataGroup dg2 = this.getDataGroup().getItem(0,2,k,run);				
//		DataGroup dg3 = this.getDataGroup().getItem(0,3,k,run);				
//		DataGroup dg4 = this.getDataGroup().getItem(0,4,k,run);				
		
		// Forward tracking residuals
		
		((H2F) dg0.getData(p_sect-1).get(0)).fill(p_mom,p_ecal_esum/1000f/p_mom);
//		((H2F) dg0.getData(e_sect-1+18).get(0)).fill(e_ecal_pcsum/1000f,e_ecal_ecsum/1000f);
		((H2F) dg0.getData(p_sect-1+18).get(0)).fill(p_ecal_esum/1000f,posi_ecal_resid.getItem(p_sect,4,0));

		for (Map.Entry<Long,Float>  entry : posi_ecal_resid.getMap().entrySet()){
			long hash = entry.getKey();
			int is = ig.getIndex(hash, 0); int ic = ig.getIndex(hash, 1); int il = ig.getIndex(hash, 2);
//			DataGroup dg1 = this.getDataGroup().getItem(il,1,k,run);				
			if(ic==3) ((H2F) dg0.getData(p_sect-1+6 ).get(0)).fill(posi_ecal_resid.getItem(p_sect,3,0),p_ecal_esum/1000f/p_mom);
			if(ic==3) ((H2F) dg0.getData(p_sect-1+12).get(0)).fill(posi_ecal_resid.getItem(p_sect,3,0),posi_ecal_resid.getItem(p_sect,4,0));
			if(ic<3) {
//			((H2F)dg1.getData(is-1+ic*6+il*12).get(0)).fill(e_the,entry.getValue());
//			((H2F)dg1.getData(is-1+ic*6).get(0)).fill(e_mom,entry.getValue());
//			((H2F)dg1.getData(is-1+ic*6).get(0)).fill(e_cz,entry.getValue());
			}
		}
		counter[p_sect-1][0]++;
		
		// Radiative photon residuals
/*		
    	for (Map.Entry<Long,List<Particle>>  entry : ecphot.getMap().entrySet()){
			int is = ig.getIndex(entry.getKey(), 0);
    		if(is==e_sect) {
            	for (Particle ec : entry.getValue()) {				
            		int    il = getDet((int) ec.getProperty("layer"));
            		float nrg =      (float) ec.getProperty("energy");
            		int    iu =        (int) ec.getProperty("iu");
            		float   t =      (float) ec.getProperty("time");
            		float pat =      (float) ec.getProperty("path");
            		float   x =      (float) ec.getProperty("x");
            		float   y =      (float) ec.getProperty("y");
//            		System.out.println("pec iS,ind,the,phi "+is+" "+il+" "+Math.toDegrees(ec.theta())+" "+Math.toDegrees(ec.phi()));
            		float thdif = (float)(e_the-Math.toDegrees(ec.theta()))*ev.tpol;
            		float phdif = (float)(e_phi-Math.toDegrees(ec.phi()))*ev.spol;

            		if(il==0 && thdif>-1 && thdif<1)   counter[e_sect-1][1]++;
            		if(il==0 && phdif>-2 && phdif<1.5) counter[e_sect-1][2]++;
    				((H2F)dg2.getData(is-1+   0+il*12).get(0)).fill(e_the,thdif); 
    				((H2F)dg2.getData(is-1+   6+il*12).get(0)).fill(e_the,phdif);
    				((H2F)dg4.getData(is-1+   0+il*12).get(0)).fill(e_mom,thdif); 
    				((H2F)dg4.getData(is-1+   6+il*12).get(0)).fill(e_mom,phdif);
    				((H2F)dg3.getData(is-1+   0+il*12).get(0)).fill(nrg/1e3,thdif);           		
    				((H2F)dg3.getData(is-1+   6+il*12).get(0)).fill(nrg/1e3,phdif);     
            	}
    		}
    	}		
*/		
	}
	
	public void fillECtime() {
		
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("ECtime");
		
		IndexGenerator ig = new IndexGenerator();
		
		DataGroup dg0 = this.getDataGroup().getItem(0,0,k,run);				
		DataGroup dg3 = this.getDataGroup().getItem(0,3,k,run);				
		DataGroup dg4 = this.getDataGroup().getItem(0,4,k,run);
		//System.out.println("Event "+getEventNumber());
    	for (Map.Entry<Long,List<Particle>>  entry : ecphot_filt.getMap().entrySet()){
			int is = ig.getIndex(entry.getKey(), 0);
    		if(is==e_sect) {
            	for (Particle ec : entry.getValue()) {			
            		int    il = getDet((int) ec.getProperty("layer"));
            		float nrg =      (float) ec.getProperty("energy");
            		int    iu =        (int) ec.getProperty("iu");
            		int    iv =        (int) ec.getProperty("iv");
            		int    iw =        (int) ec.getProperty("iw");
            		float   t =      (float) ec.getProperty("time");
            		float newt=      (float) ec.getProperty("newtime");
            		float pat =      (float) ec.getProperty("path");
            		float   x =      (float) ec.getProperty("x");
            		float   y =      (float) ec.getProperty("y");
            		int cstat =        (int) ec.getProperty("cstat");
            		
            		float thdif = (float)(e_the-Math.toDegrees(ec.theta()))*ev.tpol;
            		float phdif = (float)(e_phi-Math.toDegrees(ec.phi()))*ev.spol;
            		float  tdif = (float)((useUnsharedTime?newt:t)-ev.starttime-pat/29.97);
            		float  udif = iu-iU[il], vdif = iv-iV[il], wdif = iw-iW[il];			
            		
            		if(dbgAnalyzer) {
            		  if(nrg/1e3<0.05) System.out.println(il+" "+t+" "+" "+newt+" "+tdif+" "+" "+udif+" "+vdif+" "+wdif+" "+cstat);
            		}
            		
            		if(Math.abs(thdif)>0.3)     ((H2F)dg3.getData(is-1+0+il*6).get(0)).fill(tdif,nrg/1e3);
            		if(Math.abs(e_mom-6.5)<0.5) ((H2F)dg4.getData(is-1+0+il*12).get(0)).fill(tdif,thdif);
            		if(Math.abs(e_mom-6.5)<0.5) ((H2F)dg4.getData(is-1+6+il*12).get(0)).fill(tdif,phdif);
            		if(il<2 && Math.abs(e_mom-6.5)<0.5) ((H2F)dg0.getData(il  ).get(0)).fill(-x,y);            		
            		if(il<2 && Math.abs(e_mom-6.5)<0.5) ((H2F)dg0.getData(il+2).get(0)).fill(-x,y,tdif);            		
               		if(il<2 && Math.abs(e_mom-6.5)<0.5) ((H2F)dg0.getData(il+4).get(0)).fill(-x,y,nrg/1e3);            		
            	}
			}
		}
//    	System.out.println(" ");
		
	}
	
	public void fillSCelec() {
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("SCelec");
		
        Particle epart = e_ecal.get(0);
    	List<Particle> elecFTOF = ev.getFTOF((int)epart.getProperty("pindex"));
		    	
    	for (Particle p : elecFTOF) {
            int is = (int) p.getProperty("sector");   			
	    	int il = (int) p.getProperty("layer");
	    	int ip = (int) p.getProperty("strip");
		    Point3D xyz = getResidual(p);
		    if(il!=3) {
		    	DataGroup dg1 = this.getDataGroup().getItem(il-1,1,k,run);
		    	((H2F)dg1.getData(is-1  ).get(0)).fill(e_mom,xyz.x()); ((H2F)dg1.getData(is-1+12).get(0)).fill(ip,xyz.x());
		    	((H2F)dg1.getData(is-1+6).get(0)).fill(e_mom,xyz.y()); ((H2F)dg1.getData(is-1+18).get(0)).fill(ip,xyz.y());	
		    }
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
		((H2F)dg1.getData(1).get(0)).fill(pip_mom,pip_the);
		((H2F)dg1.getData(2).get(0)).fill(e_vz,pip_vz);
		((H2F)dg1.getData(3).get(0)).fill(pip_phi,pip_vz-e_vz);
		((H2F)dg1.getData(4).get(0)).fill(pip_the,pip_vz-e_vz);
		((H2F)dg1.getData(5).get(0)).fill(pip_vz,pip_vz-e_vz);
		((H2F)dg1.getData(6).get(0)).fill(pip_phi,epip_dPhi);
		((H2F)dg1.getData(7).get(0)).fill(pip_the,epip_dPhi);
		((H2F)dg1.getData(8).get(0)).fill(pip_vz,epip_dPhi);
		((H2F)dg1.getData(9).get(0)).fill(pip_mom,pip_beta);  
    } 
    
    public void fillTrajXY() {
		fillFTOFTrajXY();
		fillECALTrajXY();    	
    }
    
	public void fillECpim() {		
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("ECpim");
		
		DataGroup dg0 = this.getDataGroup().getItem(0,0,k,run);

		for (Particle p : pim_ecal) {
			float mom = (float) p.p();
	        float the = (float) Math.toDegrees(p.theta());
	        float phi = (float) Math.toDegrees(p.phi());
	        float pid = (float) p.getProperty("chi2pid");
	        List<Particle> pimECAL = ev.getECAL((int)p.getProperty("pindex"));
	        int is=0,il=0,iU=0; float pcsum=0, ecsum=0, esum = 0, epcal = 0;
			for (Particle pp : pimECAL) {
	    		is = (int) pp.getProperty("sector");
	    		il = (int) pp.getProperty("layer");				
				DataGroup dg1 = this.getDataGroup().getItem(getDet(il),1,k,run);				
				if(il==1) {
					iU = pp.hasProperty("iu")?(int)pp.getProperty("iu"):(int)(pp.getProperty("lu")*0.891/4.5); 
					epcal=(float)pp.getProperty("energy")/1000;
					pcsum+=(float) pp.getProperty("energy")/1000;
			    }
				if(il==4||il==7) ecsum+=(float) pp.getProperty("energy")/1000;
				Point3D xyz = getResidual(pp);
				((H2F)dg1.getData(is-1+   0).get(0)).fill(mom,xyz.x());
				((H2F)dg1.getData(is-1+   6).get(0)).fill(mom,xyz.y());	        
				((H2F)dg1.getData(is-1+  12).get(0)).fill(the,xyz.x());
				((H2F)dg1.getData(is-1+  18).get(0)).fill(the,xyz.y());
				esum+=(float) pp.getProperty("energy")/1000;
			}
			
			if(esum>0) {
				float eop = esum/mom;
				((H2F) dg0.getData(is-   1).get(0)).fill(mom,eop);
				((H2F) dg0.getData(is-1+ 6).get(0)).fill(iU,ecsum/pcsum);
				((H2F) dg0.getData(is-1+12).get(0)).fill(iU,pid);
				((H2F) dg0.getData(is-1+18).get(0)).fill(esum,pid);
			}
		}
			
	}
    
    public void fillECpipFocus() {
		int run = getRunNumber();
		DataGroup dg0 = this.getDataGroup().getItem(0,4,getDetectorTabNames().indexOf("ECpip"),run);
		DataGroup dg1 = this.getDataGroup().getItem(0,5,getDetectorTabNames().indexOf("ECpip"),run);
		
		for (Particle p : pip_ecal) {
			pip_mom  = (float) p.p();
	        pip_the  = (float) Math.toDegrees(p.theta());
	        pip_phi  = (float) Math.toDegrees(p.phi());			
	        List<Particle> pipECAL = ev.getECAL((int)p.getProperty("pindex"));
	    	List<Particle> pipFTOF = ev.getFTOF((int)p.getProperty("pindex"));   
	    	for (Particle pp : pipFTOF) {
	    		int is = (int) pp.getProperty("sector");
	    		int la = (int) pp.getProperty("layer");
	    		Point3D xyz = new Point3D(pp.getProperty("x"),pp.getProperty("y"),pp.getProperty("z"));
	    		xyz.rotateZ(Math.toRadians(-60*(is-1)));
	    		xyz.rotateY(Math.toRadians(-25)); 
	    		int laa=la==1?2:1;
	    	    if(la<3 && Math.abs(pip_the-35.2)<2)    ((H2F)dg0.getData(is-1+(laa-1)*6).get(0)).fill(pip_mom,xyz.y());
	    	    if(la<3 && Math.abs(pip_mom-0.62)<0.03) ((H2F)dg1.getData(is-1+(laa-1)*6).get(0)).fill(pip_the,xyz.y());
	    	}
	    	for (Particle pp : pipECAL) {
	    		int is = (int) pp.getProperty("sector");
	    		int la = (int) pp.getProperty("layer");
	    		Point3D xyz = new Point3D(pp.getProperty("x"),pp.getProperty("y"),pp.getProperty("z"));
	    		xyz.rotateZ(Math.toRadians(-60*(is-1)));
	    		xyz.rotateY(Math.toRadians(-25)); 
	    	    if(la==1 && Math.abs(pip_the-35.2)<2)    ((H2F)dg0.getData(is-1+12).get(0)).fill(pip_mom,xyz.y());
	    	    if(la<3  && Math.abs(pip_mom-0.62)<0.03) ((H2F)dg1.getData(is-1+12).get(0)).fill(pip_the,xyz.y());
	    	}	    	
		}    	
    }
    
    public void fillECALTrajXY() {
		int       run = getRunNumber();
		int         k = getDetectorTabNames().indexOf("ECpip");
		DataGroup dg1 = this.getDataGroup().getItem(0,3,k,run);
				
		BitSet bits = new BitSet(3);
		BitSet b000 = new BitSet(3);                                        // PCAL=0 ECIN=0 ECOU=0
		BitSet b001 = new BitSet(3);                           b001.set(2); // PCAL=0 ECIN=0 ECOU=1
		BitSet b010 = new BitSet(3);              b010.set(1);              // PCAL=0 ECIN=1 ECOU=0
		BitSet b011 = new BitSet(3);              b011.set(1); b011.set(2); // PCAL=0 ECIN=1 ECOU=1
		BitSet b100 = new BitSet(3); b100.set(0);                           // PCAL=1 ECIN=0 ECOU=0
		BitSet b101 = new BitSet(3); b101.set(0);              b101.set(2); // PCAL=1 ECIN=0 ECOU=1
		BitSet b110 = new BitSet(3); b110.set(0); b110.set(1);              // PCAL=1 ECIN=1 ECOU=0
		BitSet b111 = new BitSet(3); b111.set(0); b111.set(1); b111.set(2); // PCAL=1 ECIN=1 ECIN=1
		
		Boolean goodPC=false,goodECi=false,goodECo=false;
		
//		System.out.println("New ECAL Event");
    	if(ev.trajDet.containsKey(7)) { 
			for(int tmap : ev.trajDet.get(7)) { //retrieve ECAL trajectory planes
				int pin = ev.trajBank.getInt("pindex",tmap);
				float x = -ev.trajBank.getFloat("x",tmap); float y = ev.trajBank.getFloat("y",tmap); float z = ev.trajBank.getFloat("z",tmap);				
				int lay = getDet(ev.trajBank.getInt("layer",tmap));
				if(ev.part.get(pin).pid()==211 && ev.part.get(pin).p()>0.7) { 
					bits.clear(); 
					for (Particle p : ev.getECAL(pin)) {
						bits.set(getDet((int) p.getProperty("layer"))); //convert 1,4,7 to 0,1,2 and set bitmask
//						System.out.println(bits.toString()+"\t"+lay+" "+(-x)+" "+p.getProperty("x")+" "+y+" "+p.getProperty("y")+" "+z+" "+p.getProperty("z"));
					}
					
					switch (lay) {
					case 0: ((H2F)dg1.getData(lay+2).get(0)).fill(x,y);
					        goodPC = bits.equals(b111) || bits.equals(b100) || bits.equals(b101) || bits.equals(b110);
					        if(goodPC) ((H2F)dg1.getData(lay+7).get(0)).fill(x,y);
					        break;
					case 1: if(goodPC) ((H2F)dg1.getData(lay+2).get(0)).fill(x,y);
					        goodECi = bits.equals(b111) || bits.equals(b011) || bits.equals(b110) || bits.equals(b010);
					        if(goodPC && goodECi) ((H2F)dg1.getData(lay+7).get(0)).fill(x,y);
					        break; 
					case 2: if(goodPC) ((H2F)dg1.getData(lay+2).get(0)).fill(x,y);
					        goodECo = bits.equals(b111) || bits.equals(b001) || bits.equals(b101) || bits.equals(b011);
					        if(goodPC && goodECo) ((H2F)dg1.getData(lay+7).get(0)).fill(x,y);
					}						
				}
			}
		}   	
   	
    }
    
    public void fillFTOFTrajXY() {
		int       run = getRunNumber();
		int         k = getDetectorTabNames().indexOf("ECpip");
		DataGroup dg1 = this.getDataGroup().getItem(0,3,k,run);
				
		BitSet bits = new BitSet(2);
		BitSet b00  = new BitSet(2);                                       // P1A=0 P1B=0 
		BitSet b10  = new BitSet(2);              b10.set(0);              // P1A=1 P1B=0   
		BitSet b01  = new BitSet(2);                          b01.set(1);  // P1A=0 P1B=1  
		BitSet b11  = new BitSet(2);              b11.set(0); b11.set(1);  // P1A=1 P1B=1   

		Boolean goodP1A=false,goodP1B=false;
		
//		System.out.println("New FTOF Event");
    	if(ev.trajDet.containsKey(12)) { 
			for(int tmap : ev.trajDet.get(12)) { //retrieve FTOF trajectory planes
				int pin = ev.trajBank.getInt("pindex",tmap);
				float x = -ev.trajBank.getFloat("x",tmap); float y = ev.trajBank.getFloat("y",tmap); float z = ev.trajBank.getFloat("z",tmap);					
				int lay = ev.trajBank.getInt("layer",tmap)-1;
				if(ev.part.get(pin).pid()==211 && ev.part.get(pin).p()>0.7) { 
					bits.clear(); 
					for (Particle p : ev.getFTOF(pin)) {
						bits.set((int)p.getProperty("layer")-1);
//						System.out.println(bits.toString()+"\t"+lay+" "+(-x)+" "+p.getProperty("x")+" "+y+" "+p.getProperty("y")+" "+z+" "+p.getProperty("z"));
					}
					
					switch (lay) {
					case 1: ((H2F)dg1.getData(0).get(0)).fill(x,y);
					        goodP1B = bits.equals(b01) || bits.equals(b11);
					        if(goodP1B) ((H2F)dg1.getData(5).get(0)).fill(x,y);
					        break;
					case 0: if(goodP1B) ((H2F)dg1.getData(1).get(0)).fill(x,y);
					        goodP1A = bits.equals(b10) || bits.equals(b11) ;
					        if(goodP1B && goodP1A) ((H2F)dg1.getData(6).get(0)).fill(x,y);
					        break; 
					}						
				}
			}
		}   	
   	
    }
    public void fillECnm(String name) {	
    	
		int run = getRunNumber();
		boolean good_tagged_fiduc = false;		
		
		DataGroup ECnm0   = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf(name),run);
		DataGroup ECnm1   = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf(name),run);
		DataGroup ECprot = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECprot"),run);
		
        ((H2F)ECnm1.getData(0).get(0)).fill(nm_mom, nm_the);
        
	    float cxmm = (float) (Math.sin(nm_the*3.14159f/180f)*Math.cos(nm_phi*3.141259f/180f));
	    float cymm = (float) (Math.sin(nm_the*3.14159f/180f)*Math.sin(nm_phi*3.141259f/180f));  
	    
        if(nm_mom>0.2) {
        	float nphi = newPhi(nm_phi);
        	float   cx = (float) (Math.sin(nm_the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
        	float   cy = (float) (Math.sin(nm_the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
        	if(neutFiduc(2,cx,cy)) {
    			((H1F)ECprot.getData(e_sect-1+6).get(1)).fill(ep_MM);  	
                ((H2F)ECnm1.getData(8).get(0)).fill(-cx,cy);
                ((H1F)ECnm1.getData(10).get(0)).fill(nm_mom);
        		good_tagged_fiduc = true;
                ((H1F)ECnm0.getData(5).get(0)).fill(nm_ggp);
        	}
        }
        
        if(!makeNM()) return;
        
        for (NeutralMeson nm: nm_ecal) {
        	
        	if (nm.getMesonKin()) {
        	
        	float       cx = (float) (Math.sin(nm.the*3.14159f/180f)*Math.cos(nm.phi*3.141259f/180f));
        	float       cy = (float) (Math.sin(nm.the*3.14159f/180f)*Math.sin(nm.phi*3.141259f/180f));   
        
        	float      dcx = cxmm-cx;
        	float      dcy = cymm-cy; 
        
        	boolean  cxcut = Math.abs(dcx)<0.06;	
        	boolean  cycut = Math.abs(dcy)<0.06;
        
        	boolean goodSector = nm.getPhotonSector(0)!=e_sect && nm.getPhotonSector(1)!=e_sect;
        	
        	if(goodSector) {         		
            	((H2F)ECnm1.getData(1).get(0)).fill(dcx,dcy);
            	if(cxcut) ((H1F)ECnm1.getData(3).get(0)).fill(dcy);
            	if(cycut) ((H1F)ECnm1.getData(2).get(0)).fill(dcx);                
            	((H1F)ECnm1.getData(6).get(0)).fill(nm.mass);            	
        	}
        
        	if (goodSector && cxcut && cycut && good_tagged_fiduc) {        	
        		((H2F)ECnm1.getData(5).get(0)).fill(nm_mom, nm.mom);	    	
        		((H1F)ECnm1.getData(7).get(0)).fill(nm.mass);
        		if(nm.mom>0.2) {
        			((H2F)ECnm1.getData(4).get(0)).fill(nm.opa,nm.VG1.e()*nm.VG2.e());
        			float nphi = newPhi(nm.phi);
        			cx = (float) (Math.sin(nm.the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
        			cy = (float) (Math.sin(nm.the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
        			((H2F)ECnm1.getData(9).get(0)).fill(-cx,cy);
        			((H1F)ECnm1.getData(10).get(1)).fill(nm_mom);
        			((H2F)ECnm0.getData(0).get(0)).fill(nm.getBeta(0), nm.getEnergy(0));
        			((H2F)ECnm0.getData(1).get(0)).fill(nm.getBeta(1), nm.getEnergy(1));
        			((H1F)ECnm0.getData(2).get(0)).fill(nm.getBeta(0));
        			((H1F)ECnm0.getData(3).get(0)).fill(nm.getBeta(1));            
        			((H1F)ECnm0.getData(6).get(0)).fill(nm_ggp);
        		}
        	} 
        	}
        }
    }
    
    public void fillECneut() {
    	
        int run = getRunNumber(),taggedSector = -1;
        boolean good_tagged_fiduc = false;
        
        DataGroup ECpip   = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECpip"),run);
        DataGroup ECpip2  = this.getDataGroup().getItem(0,2,getDetectorTabNames().indexOf("ECpip"),run);
        DataGroup ECneut0 = this.getDataGroup().getItem(0,0,getDetectorTabNames().indexOf("ECneut"),run);
        DataGroup ECneut  = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECneut"),run);
        DataGroup ECphot  = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECphot"),run);
		
        ((H2F)ECneut.getData(0).get(0)).fill(neut_mom, neut_the);
        
        float cxmm = (float) (Math.sin(neut_the*3.14159f/180f)*Math.cos(neut_phi*3.141259f/180f));
        float cymm = (float) (Math.sin(neut_the*3.14159f/180f)*Math.sin(neut_phi*3.141259f/180f));
       
        if(neut_mom>0.4) {
          float nphi = newPhi(neut_phi); taggedSector = newPhiSector;
          float   cx = (float) (Math.sin(neut_the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
          float   cy = (float) (Math.sin(neut_the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
          if(neutFiduc(1,cx,cy)) {
            if(taggedSector!=e_sect) {
              ((H1F)ECpip.getData(e_sect-1+6).get(1)).fill(epip_MM);
              ((H2F)ECneut.getData(8).get(0)).fill(-cx,cy);
              ((H1F)ECneut.getData(10).get(0)).fill(neut_mom);
              good_tagged_fiduc = true; 
            }
          }
        }  
        
        neut_ecal.addAll(phot_ecal); //add photons to neutron list to get beta>0.9 neutrons
       
        if(neut_ecal.size()==0) return;
        
        boolean h12_neutron_found = false, h11_neutron_found=false;
        Vector3 pvec[] = new Vector3[3]; 
        
        List<Particle> pipECAL = ev.getECAL((int)pip_ecal.get(0).getProperty("pindex"));
        
        for(Particle p : pipECAL) {
          int   pip_sec = (int)   p.getProperty("sector");
          float pip_nrg = (float) p.getProperty("energy");
          int   pip_lay = (int)   p.getProperty("layer");
          float x = (float)p.getProperty("x");float y = (float)p.getProperty("y");float z = (float)p.getProperty("z");
          if(pip_lay==1) {((H1F) ECpip2.getData(pip_sec-1   ).get(0)).fill(pip_nrg); pvec[0]=new Vector3(x,y,z);}
          if(pip_lay==4) {((H1F) ECpip2.getData(pip_sec-1+ 6).get(0)).fill(pip_nrg); pvec[1]=new Vector3(x,y,z);}
          if(pip_lay==7) {((H1F) ECpip2.getData(pip_sec-1+12).get(0)).fill(pip_nrg); pvec[2]=new Vector3(x,y,z);}		  
        }
		
        float[] mult = new float[4];

        for (Particle p: neut_ecal) {
          List<Particle> neutECAL = ev.getECAL((int)p.getProperty("pindex"));            
          ecal_neut_sec  = (int)   neutECAL.get(0).getProperty("sector");
          ecal_neut_beta = (float) neutECAL.get(0).getProperty("beta");            
          ecal_neut_the  = (float) Math.toDegrees(p.theta());
          ecal_neut_phi  = (float) Math.toDegrees(p.phi());

          float      mass2 = neut_mom*neut_mom*(1f/(ecal_neut_beta*ecal_neut_beta)-1);
          boolean mass2cut = neut_mom<1.2?mass2>0.45:true;   
           
          ecal_neut_esum[0] = 0;  

          for (Particle pp : neutECAL) {
            ecal_neut_esum[0] += pp.getProperty("energy");  int layer = (int)pp.getProperty("layer");
            Vector3 nvec = new Vector3(pp.getProperty("x"),pp.getProperty("y"),pp.getProperty("z"));            	
            if (pvec[0]!=null && layer==1) {nvec.sub(pvec[0]);((H2F)ECneut0.getData(ecal_neut_sec-1   ).get(0)).fill(neut_mom,nvec.mag());}
            if (pvec[1]!=null && layer==4) {nvec.sub(pvec[1]);((H2F)ECneut0.getData(ecal_neut_sec-1+ 6).get(0)).fill(neut_mom,nvec.mag());}
            if (pvec[2]!=null && layer==7) {nvec.sub(pvec[2]);((H2F)ECneut0.getData(ecal_neut_sec-1+12).get(0)).fill(neut_mom,nvec.mag());};
          }
            
          float       cx = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.cos(ecal_neut_phi*3.141259f/180f));
          float       cy = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.sin(ecal_neut_phi*3.141259f/180f));  		
          float      dcx = cxmm-cx;
          float      dcy = cymm-cy;
          
          boolean  cxcut = Math.abs(dcx)<0.1;	
          boolean  cycut = Math.abs(dcy)<0.1;
    		
          if(ecal_neut_sec==e_sect) {
            ((H2F)ECphot.getData( 1).get(0)).fill(ecal_neut_esum[0],ecal_neut_the);        
            ((H2F)ECphot.getData( 3).get(0)).fill(e_the,            ecal_neut_the);        
            ((H2F)ECphot.getData(11).get(0)).fill(ecal_neut_phi,    ecal_neut_the);        
            ((H2F)ECphot.getData(12).get(0)).fill(e_phi,            ecal_neut_the);   
            mult[0]++;
          }
            
          if(good_tagged_fiduc && ecal_neut_sec==taggedSector) {        	
            ((H2F)ECneut0.getData(ecal_neut_sec-1+18).get(0)).fill(mass2,dcx);
            ((H2F)ECneut.getData(1).get(0)).fill(dcx,dcy);
            if(cxcut&&mass2cut) ((H1F)ECneut.getData(3).get(0)).fill(dcy);
            if(cycut&&mass2cut) ((H1F)ECneut.getData(2).get(0)).fill(dcx); 
 /*           	
                for (Particle pp : neutECAL) {
                if(pp.getProperty("cstat")>0) {
                System.out.println("Event:"+getEventNumber());
           	    System.out.println("cstat: "+(byte)pp.getProperty("cstat"));
            	System.out.println("ustat: "+(byte)pp.getProperty("ustat")+" vstat: "+(byte)pp.getProperty("vstat")+" wstat: "+(byte)pp.getProperty("wstat"));            			
				System.out.println("beta0: "+ecal_neut_beta+" beta: "+pp.getProperty("beta")+" time: "+pp.getProperty("time")+" newtime: "+pp.getProperty("newtime"));
				System.out.println("sector: "+pp.getProperty("sector")+" layer: "+pp.getProperty("layer")+" u: "+pp.getProperty("iu")+" v: "+pp.getProperty("iv")+" w: "+pp.getProperty("iw"));
                }                
                }
*/
            
            ((H2F)ECphot.getData(2).get(0)).fill(ecal_neut_esum[0],                  ecal_neut_the);        
            if(p.pid()==2112) ((H2F)ECphot.getData(5).get(0)).fill(ecal_neut_esum[0],ecal_neut_the);        
            ((H2F)ECphot.getData(4).get(0)).fill(e_the,                              ecal_neut_the);        
            ((H2F)ECneut.getData(4).get(0)).fill(neut_mom,                           ecal_neut_beta);        
            ((H1F)ECneut.getData(6).get(0)).fill(mass2);
            if(!h12_neutron_found) ((H1F)ECneut.getData(10).get(2)).fill(neut_mom); 
            h12_neutron_found = true;
            mult[1]++;
          }
        
          if (cxcut && cycut && good_tagged_fiduc && ecal_neut_esum[0]>0.01) {
            ((H2F)ECneut.getData(5).get(0)).fill(neut_mom, ecal_neut_beta);        
            ((H1F)ECneut.getData(7).get(0)).fill(mass2);
            if(mass2cut) {
              float nphi = newPhi(ecal_neut_phi); 
              cx = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.cos(nphi*3.141259f/180f));
              cy = (float) (Math.sin(ecal_neut_the*3.14159f/180f)*Math.sin(nphi*3.141259f/180f));
              ((H2F)ECneut.getData(9).get(0)).fill(-cx,cy);
              if(!h11_neutron_found) ((H1F)ECneut.getData(10).get(1)).fill(neut_mom);
              ((H2F)ECphot.getData(6).get(0)).fill(ecal_neut_esum[0],ecal_neut_the);  
              h11_neutron_found = true;
              mult[2]++;
              if(neut_mom<1.2) mult[3]++;
           	}
          }            
        }
     	for (int i=0; i<4&&mult[i]>0; i++) ((H1F)ECphot.getData(i+7).get(0)).fill(mult[i]);
    }
    
    public void fillECphot0() {
		int run = getRunNumber(); 
		int k = getDetectorTabNames().indexOf("ECphot");	

		int [][][] pid = new int[6][3][6];
		pid = ev.getECALPID(e_sect);
		
		for (int il=0; il<3; il++) {
			DataGroup dg = this.getDataGroup().getItem(il,0,k,run);				
			for (int id=0; id<4; id++) { int iid = id<2?4:5;
				for (int is=0; is<6; is++) {
					((H1F) dg.getData(is+id*6).get(0)).fill(pid[iid][il][is]);
					((H1F) dg.getData(is+id*6).get(1)).fill(pid[id][il][is]);
				}
			}
		}
    }
       
    public void fillECphot1() {
    	
		int run = getRunNumber();
		DataGroup dg1 = this.getDataGroup().getItem(0,1,getDetectorTabNames().indexOf("ECphot"),run);
		
    	int pin = -1;
    	
    	for (Particle p : phot_ecal) {
    		
    		List<Particle> photECAL = ev.getECAL((int)p.getProperty("pindex"));            
    		ecal_phot_sec  = (int)   photECAL.get(0).getProperty("sector");
    		ecal_phot_beta = (float) photECAL.get(0).getProperty("beta");   
            ecal_phot_esum[0] = 0;
            
            for (Particle pp : photECAL) ecal_phot_esum[0] += pp.getProperty("energy");
            
    		int pindex     = (int)   p.getProperty("pindex"); 
        	if (pindex!=pin) {
        		ecal_phot_the  = (float) Math.toDegrees(p.theta());
        		ecal_phot_phi  = (float) Math.toDegrees(p.phi());
        		pin = pindex;
        		if(ecal_phot_sec == e_sect) ((H2F)dg1.getData(0).get(0)).fill(ecal_phot_esum[0],ecal_phot_the);
        	}
    	}
    }
    
    public float newPhi(float phi) {
    	float newphi=phi;
    	newPhiSector = -1;
    	if(phi<-30) newphi=360+phi;
    	if(newphi<30)              {newPhiSector = 1; return newphi;}
    	if(newphi>30&&newphi<90)   {newPhiSector = 2; return newphi-60;}
    	if(newphi>90&&newphi<150)  {newPhiSector = 3; return newphi-120;}
    	if(newphi>150&&newphi<210) {newPhiSector = 4; return newphi-180;}
    	if(newphi>210&&newphi<270) {newPhiSector = 5; return newphi-240;}
    	if(newphi>270&&newphi<330) {newPhiSector = 6; return newphi-300;}
    	return phi;
    }
    
    //PLOT
    
    public void ECkinPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<7)        plot123(index);    	
    }
    
    public void ECelecPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<5)        plot123(index);
        if(getActive123()>4) ECelecPlotFits(index);
    }
	
    public void ECposiPlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6)        plot123(index);
        if(getActive123()>5) ECposiPlotFits(index);
    }
    
    public void ECtimePlot(String tabname) {
		int index = getDetectorTabNames().indexOf(tabname);
        if(getActive123()<6) plot123(index);    
        if(getActive123()==0) plotECtimeXY(index);
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
		if(getActive123()<3 || getActive123()>=4)  plot123(index);
		if(getActive123()==3) plotECpipXY(index);
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
	
	public void ECtrigPlot(String tabname) {
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

        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.clear(); 
        DataGroup dg5 = this.getDataGroup().getItem(0,5,index,getRunNumber());	
        
        if (getActive123()==5) {c.divide(3, 2);
        	for (int il=0; il<3; il++) {            	 
        		c.cd(il);   float y = il<1?4:10;  c.getPad(il).getAxisY().setRange(-y, y);
        		c.draw((GraphErrors)dg5.getData(il  ).get(0));        //dcres x
        		c.draw((GraphErrors)dg5.getData(il  ).get(1),"same"); //dcres y
        		c.draw((F1D)dg5.getData(il).get(2),"same"); 
        		c.cd(il+3);       y = il<1?1f:1f; c.getPad(il+3).getAxisY().setRange(-y, y);
        		c.draw((GraphErrors)dg5.getData(il+3).get(0));        //egres th
        		c.draw((GraphErrors)dg5.getData(il+3).get(1),"same"); //egres ph             		
        		c.draw((F1D)dg5.getData(il+3).get(2),"same"); 
        	} 
        return;	
        }
        
       int n=0;
        
       for (int il=0; il<3; il++) {
        	for (int ic=0; ic<2; ic++ ) {
        		for (int is=1; is<7; is++) { int isk = is+10*index;
        			if (getActive123()==6) {c.divide(6, 6);
        				c.cd(n); c.getPad(n).setAutoScale();
        				H1F h = tl.fitData.getItem(isk,ic,il,getRunNumber()).getHist(); 
        				h.setOptStat("0"); if(il==0&&ic==0) h.setTitle("SECTOR "+is);
        				h.setTitleX("DC - "+det[il]+xyz[ic]+" RESIDUAL (CM)");
        				h.setFillColor(4);
        				c.draw(h);n++;
        				c.draw(tl.fitData.getItem(isk,ic,il,getRunNumber()).getGraph(),"same");
        			}
        			if (getActive123()==7) {c.divide(6, 6);
           		     	c.cd(n); c.getPad(n).setAutoScale();
           		     	H1F h = tl.fitData.getItem(isk,ic,il+3,getRunNumber()).getHist(); 
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
	
	public void ECposiPlotFits(int index) {
	
	}
	
	public void ECpimPlotFits(int index) {
		
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
        c.clear(); c.divide(6,6);
        
        int n=0;
        for (int il=0; il<3; il++) {
        	for (int ic=0; ic<2; ic++ ) {
        		for (int is=1; is<7; is++) { int isk = is+10*index;
        			if (getActive123()==6) {
        				c.cd(n); H1F h = tl.fitData.getItem(isk,ic,il,getRunNumber()).getHist(); 
        				h.setOptStat("0"); if(il==0&&ic==0) h.setTitle("SECTOR "+is);
        				h.setTitleX("DC - "+det[il]+xyz[ic]+" RESIDUAL (CM)");
        				h.setFillColor(4);
        				c.draw(h);n++;
        				c.draw(tl.fitData.getItem(isk,ic,il,getRunNumber()).getGraph(),"same");
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
        	for (int ic=0; ic<2; ic++ ) {
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
        
        float plim[][] = {{-0.5f,-2.2f,-0.8f,-2.2f,-0.8f,-2.2f},{0.8f,2f,0.6f,2f,0.6f,2f}};
        float flim[][] = {{-1.5f,-4f,-1.5f,-4f,-1.5f,-4f},{1.5f,4f,1.5f,3f,1.5f,3f}};
        
        cfitEnable = true;
        DataGroup dg2 = this.getDataGroup().getItem(0,2,k,run);	
        DataGroup dg5 = this.getDataGroup().getItem(0,5,k,run);	
                
        for (int il=il1; il<il2; il++) {            
     	   GraphErrors g1 = ((GraphErrors)dg5.getData(il  ).get(0)); g1.reset(); //dcres x
     	   GraphErrors g2 = ((GraphErrors)dg5.getData(il  ).get(1)); g2.reset(); //dcres y
     	   GraphErrors g3 = ((GraphErrors)dg5.getData(il+3).get(0)); g3.reset(); //egres th
     	   GraphErrors g4 = ((GraphErrors)dg5.getData(il+3).get(1)); g4.reset(); //egres ph
     	   g1.setTitle("STRIP WIDTH "+((il==0)?4.5:10)+" CM");
     	   g1.setTitleX("SECTOR"); g1.setTitleY("DC-"+det[il]+" RESIDUALS (CM)");
     	   g3.setTitleX("SECTOR"); g3.setTitleY("ELEC - #gamma  "+det[il]+" (DEG)");
     	   
           DataGroup  dg1 = this.getDataGroup().getItem(il,1,k,run);	
     	   for (int ic=ic1; ic<ic2; ic++) { //dc
        	  if(ic<2) {p1=plim[0][ic+il*2]; p2=plim[1][ic+il*2]; f1=flim[0][ic+il*2]; f2=plim[1][ic+il*2];} //egres fitrange
    		  for (int is=is1; is<is2; is++) { int isk = is+10*k;        		 
    			  if(!dg1.getData(is-1+ic*6).isEmpty()) {
    			  double fl = il<1?3.5:8; //dcres fitrange
            	           tl.fitData.add(fitEngine(((H2F)dg1.getData(is-1+ic*6      ).get(0)).projectionY(),1,-fl,fl,-fl,fl),isk,ic,il,  run);
            	  if(ic<2) tl.fitData.add(fitEngine(((H2F)dg2.getData(is-1+ic*6+il*12).get(0)).projectionY(),3, p1, p2,f1,f2),isk,ic,il+3,run);                  
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
      	   GraphErrors g1 = ((GraphErrors)dg5.getData(il  ).get(0)); g1.reset();
      	   GraphErrors g2 = ((GraphErrors)dg5.getData(il  ).get(1)); g2.reset();
     	   g1.setTitle("STRIP WIDTH "+((il==0)?4.5:10)+" CM"); 
     	   g1.setTitleX("SECTOR"); g1.setTitleY("DC-"+det[il]+" RESIDUALS (CM)");
     	   
           DataGroup dg1 = this.getDataGroup().getItem(il,1,k,run);	
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
      	   GraphErrors g1 = ((GraphErrors)dg5.getData(il  ).get(0)); g1.reset();
      	   GraphErrors g2 = ((GraphErrors)dg5.getData(il  ).get(1)); g2.reset();
      	   
           DataGroup  dg1 = this.getDataGroup().getItem(il,1,k,run);	
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
    
    public GraphErrors getEff(int index, int i2, int icol) {
		DataGroup dg = this.getDataGroup().getItem(0,1,index,getRunNumber());    	
    	GraphErrors geff = new GraphErrors(); geff.setMarkerColor(icol); geff.setLineColor(icol);
    	geff.getAttributes().setTitleX("p_mm (GeV)");
    	geff.getAttributes().setTitleY("Efficiency");
    	GraphErrors g1 = ((H1F)dg.getData(10).get(0)).getGraph();
    	GraphErrors g2 = ((H1F)dg.getData(10).get(i2)).getGraph();
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
    	float hx = p.hasProperty("tx")?(float)p.getProperty("tx"):(float)p.getProperty("hx");
    	float hy = p.hasProperty("ty")?(float)p.getProperty("ty"):(float)p.getProperty("hy");
    	float hz = p.hasProperty("tz")?(float)p.getProperty("tz"):(float)p.getProperty("hz");
        Point3D xyz = new Point3D(hx-(float)p.getProperty("x"),hy-(float)p.getProperty("y"),hz-(float)p.getProperty("z"));
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
    
    public void plotECtimeXY(int index) {        
    	
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));        
        int          run = getRunNumber();
        
        c.clear(); c.divide(3,2);
        
        H2F np = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(0).get(0); 
        H2F ne = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(1).get(0); 
        H2F tp = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(2).get(0); 
        H2F te = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(3).get(0);  
        H2F ep = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(4).get(0); 
        H2F ee = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(5).get(0);  
        H2F rpt = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(6).get(0); 
        H2F ret = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(7).get(0); 
        H2F rpe = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(8).get(0); 
        H2F ree = (H2F) this.getDataGroup().getItem(0,0,index,run).getData(9).get(0); 
                
        c.cd(0); c.getPad(0).getAxisZ().setLog(true); c.draw(np);
        c.cd(3); c.getPad(3).getAxisZ().setLog(true); c.draw(ne);
        
        for(int loop = 0; loop < np.getDataBufferSize(); loop++) {
        	float nep = np.getDataBufferBin(loop);
            if (nep>0) {rpe.setDataBufferBin(loop,ep.getDataBufferBin(loop)/nep);}
        	float nee = ne.getDataBufferBin(loop);
            if (nee>0) {ree.setDataBufferBin(loop,ee.getDataBufferBin(loop)/nee);}
        }            
        
        c.cd(1); c.getPad(1).getAxisZ().setLog(true); c.getPad(1).getAxisZ().setRange(0.01,0.3); c.draw(rpe);        
        c.cd(4); c.getPad(4).getAxisZ().setLog(true); c.getPad(4).getAxisZ().setRange(0.01,0.3); c.draw(ree);         
        
        for(int loop = 0; loop < np.getDataBufferSize(); loop++) {
        	float nep = np.getDataBufferBin(loop);
            if (nep>0) {rpt.setDataBufferBin(loop,tp.getDataBufferBin(loop)/nep);}
        	float nee = ne.getDataBufferBin(loop);
            if (nee>0) {ret.setDataBufferBin(loop,te.getDataBufferBin(loop)/nee);}
        }            
        
        c.cd(2); c.getPad(2).getAxisZ().setLog(false); c.getPad(2).getAxisZ().setRange(-6,6); c.draw(rpt);        
        c.cd(5); c.getPad(5).getAxisZ().setLog(false); c.getPad(5).getAxisZ().setRange(-6,6); c.draw(ret);
           
    }    
    
    public void plotECpipXY(int index) {
      	EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));        
 	    DataGroup dg = getDataGroup().getItem(0,3,index,getRunNumber());
        
        c.clear(); c.divide(5, 3);
        
        H2F h2;
        
        for (int i=0; i<5; i++) {
    	    c.cd(i);    c.getPad().getAxisZ().setLog(getLogZ());                                          h2 = (H2F) dg.getData(i).get(0);    c.draw(h2);
    	    c.cd(i+5);  c.getPad().getAxisZ().setLog(getLogZ());                                          h2 = (H2F) dg.getData(i+5).get(0);  c.draw(h2);
    	    c.cd(i+10); c.getPad().getAxisZ().setRange(0.5,1.5); c.getPad().getAxisZ().setLog(getLogZ()); h2 = (H2F) dg.getData(i+10).get(0); c.draw(h2);       	
        }
	    	    	    
    }
       
    public void getPipEff() {    	
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("ECpip");		
		DataGroup dg1 = this.getDataGroup().getItem(0,3,k,run);
		
        for(int i=0; i<5; i++) {
        	H2F e  = (H2F) dg1.getData(i).get(0);   //efficiency denominator (events)
        	H2F w  = (H2F) dg1.getData(i+5).get(0); //efficiency numerator (weights)
        	for(int loop = 0; loop < e.getDataBufferSize(); loop++) {
        		float ne = e.getDataBufferBin(loop);
        		if (ne>0) {H2F h = (H2F) dg1.getData(i+10).get(0); h.setDataBufferBin(loop,w.getDataBufferBin(loop)/ne);}
        	}
        }		
		
    }  
    
    public void getRadMap() {    	
		int run = getRunNumber();
		int   k = getDetectorTabNames().indexOf("ECkin");		
		DataGroup dg5 = this.getDataGroup().getItem(0,5,k,run);
		
        for(int i=0; i<6; i++) {
        	H2F e  = (H2F) dg5.getData(i+6).get(0);   //denominator (events)
        	H2F w  = (H2F) dg5.getData(i+12).get(0);  //numerator (weights)
        	for(int loop = 0; loop < e.getDataBufferSize(); loop++) {
        		float ne = e.getDataBufferBin(loop);
        		if (ne>0) {H2F h = (H2F) dg5.getData(i).get(0); h.setDataBufferBin(loop,w.getDataBufferBin(loop)/ne);}
        	}
        }		
		
    } 
    
    public void showNeutronEff() {
		int index = getDetectorTabNames().indexOf("ECneut");
		EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		if(getActive123()==1) {
		c.cd(11);c.getPad().getAxisX().setRange(0., 4.0); c.getPad().getAxisY().setRange(0., 1.); 
		c.draw(getEff(index,1,4)); 
		c.draw(getEff(index,2,2),"same"); 
		c.draw(neuteff,"same");
		}
    }
    
    public void showPi0Eff() {
		int index = getDetectorTabNames().indexOf("ECpi0");
		EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.cd(11);c.getPad().getAxisX().setRange(0., 5.5); c.getPad().getAxisY().setRange(0., 1.0); c.draw(getEff(index,1,2));   	
    }
    
    public void showPhotTim() {
		int index = getDetectorTabNames().indexOf("ECtime");		  	
    }
    
/* CALIBRATION FILES */    
    
    @Override
	public void writeFile(String table, int is1, int is2, int il1, int il2, int iv1, int iv2) {
		
    	if(!dumpGraphs) return;
    	
		String path = "/Users/colesmith/CLAS12ANA/ECperf/ccdb/";
		String line = new String();
		
		try { 
			File outputFile = new File(path+table+"_"+detcal[0]);
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("ECmip.writefile("+table+")");

			for (int is=is1; is<is2; is++) {
				for (int il=il1; il<il2; il++ ) {
					for (int iv=iv1; iv<iv2; iv++) {
							switch (table) {
							case "offset": line = getOFFSETS(is,il,iv,getRunNumber()); break;
							}
						    System.out.println(line);
						    outputBw.write(line);
						    outputBw.newLine();
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
    
    public String getOFFSETS(int is, int il, int iv, int run) {
    	return "";
    }

    
    @Override
    public void timerUpdate() {

    	if(nc>100) {
    		showNeutronEff();
 //   		showPi0Eff();
    		getPipEff();
//    		getRadMap();
    		nc=0;
    	}
    	nc++;
    	return;
    }
    
}
