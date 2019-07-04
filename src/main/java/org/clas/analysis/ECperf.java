package org.clas.analysis;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.clas.tools.Event;
import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECperf extends DetectorMonitor {
	
	Event ev = new Event();
	
	public boolean goodELEC,goodPIP,goodNEUT,goodPHOT;
	
	public int Nevts, Nelecs, Ntrigs, runNum;
	boolean[] trigger_bits = new boolean[32];
	public int[] Ntrigs_sect = new int[6];
    public int[] Nelecs_sect = new int[6];
	float n_neut_ecal = 0, n_neut_mm = 0, n_neut_mm_save=0;

	public float EB, Eb, Mp;
	public float RFT, STT;
	public long TriggerWord;
	public float rfPeriod;
	public int rf_large_integer;
	
	public LorentzVector VB, VT, Ve, VGS, Vprot, Vpip, VG1, VG2, VPI0;
	public boolean found_eTraj, found_eECAL, found_eFTOF1a, found_eFTOF1b, found_eLTCC, found_eHTCC;
	public int   e_part_ind, e_sect, e_FTOF_pad1a, e_FTOF_pad1b, e_HTCC_bin_phi, e_HTCC_bin_theta;
	public float e_mom, e_the, e_phi, e_vx, e_vy, e_vz;
	public float e_xB, e_Q2, e_W;
	public float e_HTCC_tX, e_HTCC_tY, e_HTCC_tZ, e_LTCC_tX, e_LTCC_tY, e_FTOF_tX, e_FTOF_tY, e_PCAL_tX, e_PCAL_tY;
	public float e_DCSL1_tX, e_DCSL1_tY, e_DCSL2_tX, e_DCSL2_tY, e_DCSL3_tX, e_DCSL3_tY, e_DCSL4_tX, e_DCSL4_tY, e_DCSL5_tX, e_DCSL5_tY, e_DCSL6_tX, e_DCSL6_tY;
	public float e_PCAL_X, e_PCAL_Y, e_PCAL_Z, e_PCAL_edep, e_PCAL_t, e_EC_ein, e_EC_eout, e_EC_etot, e_PCAL_path, e_PCAL_vt;
	public float e_FTOF1a_X, e_FTOF1a_Y, e_FTOF1a_Z, e_FTOF1a_edep, e_FTOF1a_t, e_FTOF1a_path, e_FTOF1a_vt;
	public float e_FTOF1b_X, e_FTOF1b_Y, e_FTOF1b_Z, e_FTOF1b_edep, e_FTOF1b_t, e_FTOF1b_path, e_FTOF1b_vt;
	public float e_LTCC_X, e_LTCC_Y, e_LTCC_Z, e_LTCC_t, e_LTCC_nphe, e_LTCC_path, e_LTCC_vt;
	public float e_HTCC_X, e_HTCC_Y, e_HTCC_Z, e_HTCC_theta, e_HTCC_phi, e_HTCC_t, e_HTCC_nphe, e_HTCC_path, e_HTCC_vt;
    public float e_ecal_esum;
    
	public int prot_part_ind;
	public float prot_mom, prot_the, prot_phi, prot_vx, prot_vy, prot_vz, prot_beta;

	public int   pim_part_ind, pip_part_ind, pip_FTOF_pad1b;
	public float pip_mom, pip_the, pip_phi, pip_vx, pip_vy, pip_vz, pip_beta, pip_FTOF1b_t, pip_FTOF1b_path, pip_FTOF1b_vt;
	public float pip_FTOF1a_t, pip_FTOF1a_path, pip_FTOF1a_vt;
	public float pip_ecal_esum;
	
	public float pim_mom, pim_FTOF1a_t, pim_FTOF1a_path, pim_FTOF1a_vt;
	public float pim_FTOF1b_t, pim_FTOF1b_path, pim_FTOF1b_vt;
	public float thisTime;

	public int   G1_part_ind, G2_part_ind, G1_pcal_ind, G2_pcal_ind, G1_cal_layers, G2_cal_layers;
	public float G1_mom, G1_the, G1_phi, G2_mom, G2_the, G2_phi;
	public float G1_pcal_X, G1_pcal_Y, G1_pcal_Z, G1_pcal_t, G1_pcal_R, G1_pcal_vt, G2_pcal_X, G2_pcal_Y, G2_pcal_Z, G2_pcal_t, G2_pcal_R, G2_pcal_vt;
	public float G1_ein_t, G1_eout_t, G2_ein_t, G2_eout_t;

	public float elast_dPhi, elast_EB;
	public float epip_dPhi, epip_MM;
	public float pi0_mass, pi0_E, pi0_the, pi0_phi, pi0_open;
	
	public float neut_mom,neut_the,neut_phi;
	public float ecal_neut_the,ecal_neut_phi,ecal_neut_beta;
	public float ecal_phot_the,ecal_phot_phi,ecal_phot_beta;
	public int   ecal_neut_sec, ecal_phot_sec;
	public int[] ecal_neut_esum = new int[6];
	public int[] ecal_phot_esum = new int[6];

	public H2F   H_e_t_f, H_e_p_f, H_e_vz_f, H_e_vt_vz, H_e_vt_p, H_e_vt_t;
	public H2F   H_e_PCAL, H_e_FTOF, H_e_LTCC, H_e_DCSL6, H_e_DCSL5, H_e_DCSL4, H_e_DCSL3, H_e_DCSL2, H_e_DCSL1, H_e_HTCC;
    public H2F   H_e_nphe_HTCC, H_e_bin_theta_HTCC, H_e_bin_phi_HTCC, H_e_theta_HTCC, H_e_phi_HTCC;
    
	public H2F[] H_e_HTCC_cut = new H2F[6]; 
	public H2F[] H_e_t_p=new H2F[6], H_e_vz_t=new H2F[6], H_e_vz_p=new H2F[6];
	public H2F[] H_e_EC_etot_p=new H2F[6], H_e_EC_vt_theta=new H2F[6], H_e_EC_XY=new H2F[6];
	public H2F[] H_e_FTOF_vt_pad1a=new H2F[6], H_e_FTOF_edep_pad1a=new H2F[6], H_e_FTOF_XY_pad1a=new H2F[6];
	public H2F[] H_e_FTOF_vt_pad1b=new H2F[6], H_e_FTOF_edep_pad1b=new H2F[6], H_e_FTOF_XY_pad1b=new H2F[6];
	public H2F[] H_FTOF_pos_beta_mom_pad1a=new H2F[6], H_FTOF_neg_beta_mom_pad1a=new H2F[6], H_FTOF_pos_beta_mom_pad1b=new H2F[6], H_FTOF_neg_beta_mom_pad1b=new H2F[6];
	public H2F[] H_FTOF_pos_mass_mom_pad1a=new H2F[6], H_FTOF_pos_mass_the_pad1a=new H2F[6], H_FTOF_neg_mass_mom_pad1a=new H2F[6], H_FTOF_neg_mass_the_pad1a=new H2F[6];
	public H2F[] H_FTOF_pos_mass_mom_pad1b=new H2F[6], H_FTOF_pos_mass_the_pad1b=new H2F[6], H_FTOF_neg_mass_mom_pad1b=new H2F[6], H_FTOF_neg_mass_the_pad1b=new H2F[6];
	//from tof_monitor.java for timing and gain calibration, but from Dan's comment
	//use leptons/pions (both charges) for p1a, p1b and all particles (both charges) for p2.
	public H2F[] p1a_pad_vt_elec=new H2F[6], p1a_pad_vt_pion=new H2F[6], p1b_pad_vt_elec=new H2F[6], p1b_pad_vt_pion=new H2F[6], p2_pad_vt=new H2F[6];
	public H1F[] p1a_pad_edep_elec=new H1F[6], p1a_pad_edep_pion=new H1F[6], p1b_pad_edep_elec=new H1F[6], p1b_pad_edep_pion=new H1F[6], p2_pad_edep=new H1F[6];

	public H2F[] H_e_LTCC_vt_theta=new H2F[6], H_e_LTCC_nphe_theta=new H2F[6], H_e_LTCC_XY=new H2F[6];
	public H2F[] H_e_HTCC_vt_theta=new H2F[6], H_e_HTCC_nphe_theta=new H2F[6], H_e_HTCC_XY=new H2F[6];
	public H1F[][][] H_e_bin_nphe_HTCC;

	public H2F   H_elast_e_th_p, H_elast_p_th_p, H_elast_vz_vz, H_elast_dvz_phi, H_elast_dvz_theta_all, H_elast_dvz_vz;
	public H2F   H_elast_Dphi_phi, H_elast_Dphi_theta, H_elast_Dphi_vz, H_elast_EB_phi, H_elast_EB_theta, H_elast_EB_vz ;
	public H2F[] H_elast_W_theta=new H2F[6], H_elast_W_Q2=new H2F[6], H_elast_inc_W_theta=new H2F[6], H_elast_dvz_theta=new H2F[6];

	public H2F   H_epip_e_th_p, H_epip_p_th_p, H_epip_vz_vz, H_epip_dvz_phi, H_epip_dvz_theta, H_epip_dvz_vz;
	public H2F   H_epip_Dphi_phi, H_epip_Dphi_theta, H_epip_Dphi_vz, H_epip_beta_p, H_epip_FTOF1b_dt_epad, H_epip_FTOF1b_dt_pippad;
	public H2F[] H_epip_W_theta=new H2F[6], H_epip_inc_W_theta=new H2F[6];
	public H1F[] H_epip_W=new H1F[6];
	
	public H2F H_pi0_G1_XY, H_pi0_G1_TR, H_pi0_G1_vt_evt, H_pi0_G1_layer_E;//4
	public H2F H_pi0_G2_XY, H_pi0_G2_TR, H_pi0_G2_vt_evt, H_pi0_G2_layer_E;//8
	public H2F H_pi0_G1_mom_the, H_pi0_G1_phi_the, H_pi0_G2_mom_the, H_pi0_G2_phi_the;//12
	public H2F H_pi0_open_E, H_pi0_E_the, H_pi0_phi_the;//15
    public H1F H_pi0_mass, H_pi0_G1_layers, H_pi0_G2_layers;//18
    
    public H2F   H_neut_e_th_p, H_neut_dth_dph, H_neut_p_beta,H_neut_p_beta_cut, H_neut_esum_the_elec, H_neut_esum_the;
    public H2F   H_neut_phi_the,H_neut_phie_the;
    public H2F[] H_neut_rad_tail = new H2F[6];
    public H2F[] H_neut_phi_the_eff = new H2F[6];
    public H1F[] H_neut_mass2 = new H1F[6];
    public H1F   H_neut_phi1_phi2,H_neut_th1_th2, H_neut_avg_mom;
    
    public H2F   H_phot_esum_the_elec;
    
	public IndexedTable rfTable;	
	
    public ECperf(String name) {
        super(name);
        this.setDetectorTabNames("ECelec",
        		                 "ECpip",
        		                 "ECpi0",
        		                 "ECneut",
        		                 "ECphot",
        		                 "ECtime");
        
        this.init();
        this.localinit();      
    }
    
    public void localinit() {
    	System.out.println("ECperf.localinit()");
        configEngine("muon"); 
    	tl.setFitData(Fits);   
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
		trigger_bits = new boolean[32];
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
    
    @Override    
    public void createHistos(int run) {  
	    System.out.println("ECperf:createHistos("+run+")");
    	setRunNumber(run);
    	dstinit(run);
    	
		for(int s=0;s<6;s++){    	
			H_e_EC_etot_p[s] = new H2F(String.format("H_e_EC_etot_p_%d",s+1),String.format("H_e_EC_etot_p_%d",s+1),100,0,EB,100,0.0,EB/4);
			H_e_EC_etot_p[s].setTitle(String.format("ECAL vs p S%d",s+1));
			H_e_EC_etot_p[s].setTitleX("p (GeV)");
			H_e_EC_etot_p[s].setTitleY("ECAL (GeV)");
		}
		
		H_epip_e_th_p = new H2F("H_epip_e_th_p","H_epip_e_th_p",100,0,EB,100,0,40);
		H_epip_e_th_p.setTitle("electron #theta vs p");
		H_epip_e_th_p.setTitleX("p (GeV)");
		H_epip_e_th_p.setTitleY("#theta (^o)");
		H_epip_p_th_p = new H2F("H_epip_p_th_p","H_epip_p_th_p",100,0,4,100,0,50);
		H_epip_p_th_p.setTitle("pip #theta vs p");
		H_epip_p_th_p.setTitleX("p (GeV)");
		H_epip_p_th_p.setTitleY("#theta (^o)");
		H_epip_vz_vz = new H2F("H_epip_vz_vz","H_epip_vz_vz",100,-15,15,100,-15,15);
		H_epip_vz_vz.setTitle("epipic vz p vs e");
		H_epip_vz_vz.setTitleX("e vz (cm)");
		H_epip_vz_vz.setTitleY("pip vz (cm)");
		H_epip_dvz_phi = new H2F("H_epip_dvz_phi","H_epip_dvz_phi",100,-180,180,100,-15,15);
		H_epip_dvz_phi.setTitle("e pip #Delta vz vs #phi");
		H_epip_dvz_phi.setTitleX("#phi (^o)");
		H_epip_dvz_phi.setTitleY("#Delta vz (cm)");
		H_epip_dvz_theta = new H2F("H_epip_dvz_theta","H_epip_dvz_theta",100,0,80,100,-15,15);
		H_epip_dvz_theta.setTitle("e pip #Delta vz vs #theta");
		H_epip_dvz_theta.setTitleX("#theta (^o)");
		H_epip_dvz_theta.setTitleY("#Delta vz (cm)");
		H_epip_dvz_vz = new H2F("H_epip_dvz_vz","H_epip_dvz_vz",100,-15,15,100,-15,15);
		H_epip_dvz_vz.setTitle("e pip #Delta vz vs z");
		H_epip_dvz_vz.setTitleX("vz (cm)");
		H_epip_dvz_vz.setTitleY("#Delta vz (cm)");
		H_epip_Dphi_phi = new H2F("H_epip_Dphi_phi","H_epip_Dphi_phi",100,-180,180,100,-180,180);
		H_epip_Dphi_phi.setTitle("#Delta#phi vs #phi");
		H_epip_Dphi_phi.setTitleX("#phi (^o)");
		H_epip_Dphi_phi.setTitleY("#Delta#phi (^o)");
		H_epip_Dphi_theta = new H2F("H_epip_Dphi_theta","H_epip_Dphi_theta",100,0,80,100,-180,180);
		H_epip_Dphi_theta.setTitle("#Delta#phi vs #theta");
		H_epip_Dphi_theta.setTitleX("#theta (^o)");
		H_epip_Dphi_theta.setTitleY("#Delta#phi (^o)");
		H_epip_Dphi_vz = new H2F("H_epip_Dphi_vz","H_epip_Dphi_vz",100,-15,15,100,-180,180);
		H_epip_Dphi_vz.setTitle("#Delta#phi vs vz");
		H_epip_Dphi_vz.setTitleX("vz (cm)");
		H_epip_Dphi_vz.setTitleY("#Delta#phi (^o)");

		H_epip_beta_p = new H2F("H_epip_beta_p","H_epip_beta_p",100,0,EB,100,0.8,1.2);
		H_epip_beta_p.setTitle("pip #beta vs p");
		H_epip_beta_p.setTitleX("p (GeV)");
		H_epip_beta_p.setTitleY("#beta");    
		
		for(int s=0;s<6;s++){
			H_epip_W_theta[s] = new H2F("H_epip_W_theta_"+(s+1),"H_epip_W_theta_"+(s+1),100,0,4,100,0,40);
			H_epip_W_theta[s].setTitle("#theta vs W S"+(s+1));
			H_epip_W_theta[s].setTitleX("W ( GeV)");
			H_epip_W_theta[s].setTitleY("#theta (^o)");
			H_epip_inc_W_theta[s] = new H2F("H_epip_inc_W_theta_"+(s+1),"H_epip_inc_W_theta_"+(s+1),100,0,4,100,0,40);
			H_epip_inc_W_theta[s].setTitle("#theta vs W S"+(s+1));
			H_epip_inc_W_theta[s].setTitleX("W ( GeV)");
			H_epip_inc_W_theta[s].setTitleY("#theta (^o)");
			H_epip_W[s] = new H1F("H_epip_W_"+(s+1),"H_epip_W_"+(s+1),100,0,4);
		}
		
		H_neut_e_th_p = new H2F("H_neut_e_th_p","H_neut_e_th_p",50,0,2.5,50,0,40);
		H_neut_e_th_p.setTitle("neutron #theta vs p");
		H_neut_e_th_p.setTitleX("p_mm (GeV)");
		H_neut_e_th_p.setTitleY("#theta_mm (^o)");
		
		H_neut_dth_dph = new H2F("H_neut_dth_dph","H_neut_dth_dph",50,-20,20,50,-180,180);
		H_neut_dth_dph.setTitle("#delta #phi vs #delta #theta");
		H_neut_dth_dph.setTitleX("#theta_mm - #theta_ecal (^o)");
		H_neut_dth_dph.setTitleY("#phi_mm - #phi_ecal (^o)");
		
		H_neut_p_beta = new H2F("H_neut_p_beta","H_neut_p_beta",50,0,2.5,50,0,1.0);
		H_neut_p_beta.setTitle("#beta_ecal vs p_mm (no electron)");
		H_neut_p_beta.setTitleX("p_mm (GeV)");
		H_neut_p_beta.setTitleY("#beta_ecal");	
		
		H_neut_p_beta_cut = new H2F("H_neut_p_beta_cut","H_neut_p_beta_cut",50,0,2.5,50,0,1.0);
		H_neut_p_beta_cut.setTitle("#beta_ecal vs p_mm");
		H_neut_p_beta_cut.setTitleX("p_mm (GeV)");
		H_neut_p_beta_cut.setTitleY("#beta_ecal");	
		
		H_neut_esum_the = new H2F("H_neut_esum_the","H_neut_esum_the",50,0,300.,50,0,35.);
		H_neut_esum_the.setTitle("#theta_ecal vs esum_ecal (no electron)");
		H_neut_esum_the.setTitleX("esum_ecal (MeV)");
		H_neut_esum_the.setTitleY("#theta_ecal (deg)");		
		
		H_neut_phi_the = new H2F("H_neut_phi_the","H_neut_phi_the",80,-180,180.,80,5,22.);
		H_neut_phi_the.setTitle("#theta_ecal vs #phi_ecal (electron)");
		H_neut_phi_the.setTitleX("#phi_ecal (deg)");
		H_neut_phi_the.setTitleY("#theta_ecal (deg)");	
		
		H_neut_phie_the = new H2F("H_neut_phie_the","H_neut_phie_the",80,-180,180.,80,5,22.);
		H_neut_phie_the.setTitle("#theta_ecal vs #phi_e (electron)");
		H_neut_phie_the.setTitleX("#phi_e (deg)");
		H_neut_phie_the.setTitleY("#theta_ecal (deg)");		
		
		H_neut_esum_the_elec = new H2F("H_neut_esum_the_elec","H_neut_esum_the_elec",50,0,300.,50,0,35.);
		H_neut_esum_the_elec.setTitle("#theta_ecal vs esum_ecal (electron)");
		H_neut_esum_the_elec.setTitleX("esum_ecal (MeV)");
		H_neut_esum_the_elec.setTitleY("#theta_ecal (deg)");	
		
		H_neut_phi1_phi2 = new H1F("H_neut_phi1_phi2","H_neut_phi1_phi2",50,-180.,180.);
		H_neut_phi1_phi2.setTitle("#phi_mm - #phi_ecal");
		H_neut_phi1_phi2.setTitleX("#phi_mm - #phi_ecal (deg)");	
		
		H_neut_th1_th2 = new H1F("H_neut_th1_th2","H_neut_th1_th2",50,-20.,20.);
		H_neut_th1_th2.setTitle("#theta_mm - #theta_ecal");
		H_neut_th1_th2.setTitleX("#theta_mm - #theta_ecal (deg)");
		
		H_neut_avg_mom = new H1F("H_neut_avg_mom","H_neut_avg_mom",50,0.,2.5);
		H_neut_avg_mom.setTitle("p_mm (Mass^2 > 0.3)");
		H_neut_avg_mom.setTitleX("p_mm (GeV)");
		
		String[] tit1 = {"(electron)","(no electron)"};
		String[] tit2 = {"(#theta,#phi cut)","(no electron)"};
		String[] tit3 = {"mm","ecal"};
		for (int is=0; is<2; is++ ) {
			H_neut_rad_tail[is] = new H2F("H_neut_rad_tail_"+is,"H_neut_rad_tail_"+is,50,0.,25.,50,0.,25.);
			H_neut_rad_tail[is].setTitle("#theta_elec vs #theta_ecal "+tit1[is]);
			H_neut_rad_tail[is].setTitleX("#theta_ecal (deg)");
			H_neut_rad_tail[is].setTitleY("#theta_elec (deg)");
			H_neut_mass2[is] = new H1F("H_neut_mass2_"+is,"H_neut_mass2_"+is,60,-0.1,2.0);
			H_neut_mass2[is].setTitle("Mass^2 "+tit2[is]);
			H_neut_mass2[is].setTitleX("Mass^2 (GeV^2)");
			H_neut_phi_the_eff[is] = new H2F("H_neut_phi_the_eff_"+tit3[is],"H_neut_phi_the_eff_"+tit3[is],40,-40.,40.,21,5.,37.);
			H_neut_phi_the_eff[is].setTitle("#theta_"+tit3[is]+" vs phi_"+tit3[is]);
			H_neut_phi_the_eff[is].setTitleX("#phi_"+tit3[is]+" (deg)");
			H_neut_phi_the_eff[is].setTitleY("#theta_"+tit3[is]+" (deg)");
		}
		
		H_phot_esum_the_elec = new H2F("H_phot_esum_the_elec","H_neut_phot_the_elec",50,0,500.,60,5,25.);
		H_phot_esum_the_elec.setTitle("#theta_ecal vs esum_ecal (electron)");
		H_phot_esum_the_elec.setTitleX("esum_ecal (MeV)");
		H_phot_esum_the_elec.setTitleY("#theta_ecal (deg)");		
    }
    
	public void myinit(){
		goodELEC = false; goodPIP  = false; goodNEUT = false; goodPHOT = false;   
		e_EC_etot = 0; e_PCAL_edep=0; e_EC_ein=0; e_EC_eout=0;
    }
	
    public double Vangle(Vector3 v1, Vector3 v2){
        double res = 0;
        double l1 = v1.mag();
        double l2 = v2.mag();
        double prod = v1.dot(v2);
        if( l1 * l2 !=0 && Math.abs(prod)<l1*l2 )res = Math.toDegrees( Math.acos(prod/(l1*l2) ) );
        return res; 
    }
   
    public int getSect(DataBank bank, int partInd){
        for(int k = 0; k < bank.rows(); k++) if(bank.getShort("pindex",k)==partInd)return bank.getInt("sector",k);  
        return -1; 
    }  
    
	public boolean fillConfBank(DataBank confbank){
		int[] tb = new int[32];
		boolean selectTrig = false;
		long TriggerWord = confbank.getLong("trigger",0);		
		for (int i = 31; i >= 0; i--) {trigger_bits[i] = (TriggerWord & (1 << i)) != 0;} 
		for (int i = 31; i >= 0; i--) {tb[i] = ((trigger & (1 << i))!=0)?1:0;}		
		if(trigger_bits[1] || trigger_bits[2] || trigger_bits[3] || trigger_bits[4] || trigger_bits[5] || trigger_bits[6]){
			selectTrig = true;
			Ntrigs++;
			if(trigger_bits[1])Ntrigs_sect[0]++;
			if(trigger_bits[2])Ntrigs_sect[1]++;
			if(trigger_bits[3])Ntrigs_sect[2]++;
			if(trigger_bits[4])Ntrigs_sect[3]++;
			if(trigger_bits[5])Ntrigs_sect[4]++;
			if(trigger_bits[6])Ntrigs_sect[5]++;
		}
		System.out.println(Ntrigs_sect[0]+" "+Ntrigs_sect[1]+" "+Ntrigs_sect[2]+" "+Ntrigs_sect[3]+" "+Ntrigs_sect[4]+" "+Ntrigs_sect[5]);
		System.out.println("Bits "+tb[0]+" "+tb[1]+tb[2]+tb[3]+tb[4]+tb[5]+tb[6]);
		return selectTrig;
	}
	
	public void fillRecBank(DataBank recBank){
	    STT = isHipo3Event ? recBank.getFloat("STTime", 0):recBank.getFloat("startTime", 0);
		RFT = recBank.getFloat("RFTime",0);
	}
	
	public void fillECAL(DataBank bank){
		for(int r=0;r<bank.rows();r++){
			if(bank.getShort("pindex",r)==e_part_ind){
				if(bank.getByte("layer",r)==1){
					found_eECAL = true;
					e_PCAL_X = bank.getFloat("x",r);
					e_PCAL_Y = bank.getFloat("y",r);
					e_PCAL_Z = bank.getFloat("z",r);
					e_PCAL_edep += bank.getFloat("energy",r);
					e_PCAL_t = bank.getFloat("time",r);
					e_PCAL_path = bank.getFloat("path",r);
					e_PCAL_vt = e_PCAL_t - e_PCAL_path/29.98f - STT;
					if(e_sect==0)e_sect = bank.getByte("sector",r);
				}
				if(bank.getByte("layer",r)==4)e_EC_ein += bank.getFloat("energy",r);
				if(bank.getByte("layer",r)==7)e_EC_eout += bank.getFloat("energy",r);
			}
			if(bank.getShort("pindex",r)==G1_part_ind){
				if(bank.getByte("layer",r)==1 && G1_pcal_ind==-1){
					G1_cal_layers++;G1_pcal_ind = r;
					G1_pcal_X = bank.getFloat("x",r);
					G1_pcal_Y = bank.getFloat("y",r);
					G1_pcal_Z = bank.getFloat("z",r);
					G1_pcal_t = bank.getFloat("time",r);
					//G1_pcal_R = G1_pcal_X*G1_pcal_X + G1_pcal_Y*G1_pcal_Y + (G1_pcal_Z-e_vz)*(G1_pcal_Z-e_vz);
					G1_pcal_R = G1_pcal_X*G1_pcal_X + G1_pcal_Y*G1_pcal_Y + G1_pcal_Z*G1_pcal_Z;
					G1_pcal_R = (float)Math.sqrt(G1_pcal_R);
					G1_pcal_vt = G1_pcal_t - G1_pcal_R/29.98f - STT;
				}
				//else if(bank.getByte("layer",r)==1 && G1_pcal_ind!=-1){
				//	System.out.println("error: found a photon with TWO PCAL CLUSTERS");
				//}
				else if(bank.getByte("layer",r)>1)G1_cal_layers++;
			}
			if(bank.getShort("pindex",r)==G2_part_ind){
				if(bank.getByte("layer",r)==1 && G2_pcal_ind==-1){
					G2_cal_layers++;G2_pcal_ind = r;
					G2_pcal_X = bank.getFloat("x",r);
					G2_pcal_Y = bank.getFloat("y",r);
					G2_pcal_Z = bank.getFloat("z",r);
					G2_pcal_t = bank.getFloat("time",r);
					//G2_pcal_R = G2_pcal_X*G2_pcal_X + G2_pcal_Y*G2_pcal_Y + (G2_pcal_Z-e_vz)*(G2_pcal_Z-e_vz);
					G2_pcal_R = G2_pcal_X*G2_pcal_X + G2_pcal_Y*G2_pcal_Y + G2_pcal_Z*G2_pcal_Z;
					G2_pcal_R = (float)Math.sqrt(G2_pcal_R);
					G2_pcal_vt = G2_pcal_t - G2_pcal_R/29.98f - STT;
				}
			}
		}
		e_EC_etot = e_PCAL_edep+e_EC_ein+e_EC_eout;
	}
	
	public void fillTraj(DataBank trajBank){
		for(int r=0;r<trajBank.rows();r++){
			//System.out.println("comparing "+e_part_ind+" and "+trajBank.getInt("pindex",r));
			if(trajBank.getShort("pindex",r)==e_part_ind){
				found_eTraj=true;
				switch(trajBank.getShort("detId",r)) {
					case 0:
						e_HTCC_tX = trajBank.getFloat("x",r);
						e_HTCC_tY = trajBank.getFloat("y",r);
						e_HTCC_tZ = trajBank.getFloat("z",r);
						e_HTCC_phi = (float)Math.toDegrees(Math.atan2(e_HTCC_tY,e_HTCC_tX)) + 30f;
						while(e_HTCC_phi<0)e_HTCC_phi+=60;
						while(e_HTCC_phi>60)e_HTCC_phi-=60;
						float htccR = (float)Math.sqrt(e_HTCC_tX*e_HTCC_tX+e_HTCC_tY*e_HTCC_tY+e_HTCC_tZ*e_HTCC_tZ);
						//System.out.println("htccR="+htccR);
						e_HTCC_theta = (float)Math.toDegrees(Math.acos( e_HTCC_tZ/htccR ));
						e_HTCC_bin_theta = -1;e_HTCC_bin_phi=-1;
						if(e_HTCC_theta>5 && e_HTCC_theta<10)e_HTCC_bin_theta = (int) ((e_HTCC_theta - 5f)/1f);
						if(e_HTCC_theta>10 && e_HTCC_theta<20)e_HTCC_bin_theta = 5 + (int) ((e_HTCC_theta - 10f)/2f);
						if(e_HTCC_theta>20 && e_HTCC_theta<40)e_HTCC_bin_theta = 10 + (int) ((e_HTCC_theta - 20f)/4f);
						if(e_HTCC_phi<20)e_HTCC_bin_phi = (int) (e_HTCC_phi/4f);
						if(e_HTCC_phi>20 && e_HTCC_phi<40)e_HTCC_bin_phi = 5 + (int) ((e_HTCC_phi-20f)/1f);
						if(e_HTCC_phi>40)e_HTCC_bin_phi = 25 + (int) ((e_HTCC_phi-40f)/4f);
					case 12:
						e_DCSL1_tX = trajBank.getFloat("x",r);
						e_DCSL1_tY = trajBank.getFloat("y",r);
					case 18:
						e_DCSL2_tX = trajBank.getFloat("x",r);
						e_DCSL2_tY = trajBank.getFloat("y",r);
					case 24:
						e_DCSL3_tX = trajBank.getFloat("x",r);
						e_DCSL3_tY = trajBank.getFloat("y",r);
					case 30:
						e_DCSL4_tX = trajBank.getFloat("x",r);
						e_DCSL4_tY = trajBank.getFloat("y",r);
					case 36:
						e_DCSL5_tX = trajBank.getFloat("x",r);
						e_DCSL5_tY = trajBank.getFloat("y",r);
					case 42:
						e_DCSL6_tX = trajBank.getFloat("x",r);
						e_DCSL6_tY = trajBank.getFloat("y",r);
					case 43:
						e_LTCC_tX = trajBank.getFloat("x",r);
						e_LTCC_tY = trajBank.getFloat("y",r);
					case 46:
						e_FTOF_tX = trajBank.getFloat("x",r);
						e_FTOF_tY = trajBank.getFloat("y",r);
					case 47:
						e_PCAL_tX = trajBank.getFloat("x",r);
						e_PCAL_tY = trajBank.getFloat("y",r);
				}
			}
		}
	}
	
    public boolean makeELEC(){
    	
        List<Particle> nlist = ev.getParticle(11);
        
        Particle epart = nlist.get(0);
        
        e_mom = (float) epart.p();       
        e_vz  = (float) epart.vz();
        short status = (short) epart.getProperty("status");
        boolean inDC = (status>=2000 && status<3000);
        
        e_ecal_esum = 0f;
    	for (Particle p : nlist) e_ecal_esum += p.getProperty("energy");        	
        
        if(inDC && Math.abs(e_vz+3)<12 && (e_mom>1.5 || runNum<2600 ) ){
        	e_sect = (int)   epart.getProperty("sector");
            e_the  = (float) Math.toDegrees(epart.theta());
            e_phi  = (float) Math.toDegrees(epart.phi());
            e_vx   = (float) epart.vx(); 
            e_vy   = (float) epart.vy();
            Ve     =         epart.vector();
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
    
    public boolean makePIP() {
    	
        List<Particle> nlist = ev.getParticle(211);
        
        Particle pipart = nlist.get(0);
        
        pip_mom  = (float) pipart.p();
        short status = (short) pipart.getProperty("status");
        boolean inDC = (status>=2000 && status<3000);
        
        pip_ecal_esum = 0f;
    	for (Particle p : nlist) pip_ecal_esum += p.getProperty("energy");        	
        if(inDC && (pip_mom>0.5||runNum<2600) ){ 
            pip_the  = (float) Math.toDegrees(pipart.theta());
            pip_phi  = (float) Math.toDegrees(pipart.phi());
            pip_vz   = (float) pipart.vz();
        	pip_beta = (float) pipart.getProperty("beta");
        	Vpip     =         pipart.vector();
            return true;       	
        }
        return false;
    }
    
    public boolean makeNEUT() {
    	
        List<Particle> nlist = ev.getParticle(2112);
        
        for (int is=0; is<6; is++) ecal_neut_esum[is]=0;
        if(nlist.size()==1) {
        	ecal_neut_sec  = (int)   nlist.get(0).getProperty("sector");
        	ecal_neut_the  = (float) Math.toDegrees(nlist.get(0).theta());
        	ecal_neut_phi  = (float) Math.toDegrees(nlist.get(0).phi());
        	ecal_neut_beta = (float) nlist.get(0).getProperty("beta");
        	for (Particle p : nlist) ecal_neut_esum[0] += p.getProperty("energy");        	
        	return true;
        }
    	return false;
    }
    
    public boolean makePHOT() {
    	
        List<Particle> nlist = ev.getParticle(22);
        
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
    
	public void makeOthers(DataBank bank){
		for(int k = 0; k < bank.rows(); k++){
			int pid = bank.getInt("pid", k);
			float px = bank.getFloat("px", k);
			float py = bank.getFloat("py", k);
			float pz = bank.getFloat("pz", k);
			float vx = bank.getFloat("vx", k);
			float vy = bank.getFloat("vy", k);
			float vz = bank.getFloat("vz", k);
			float be = bank.getFloat("beta", k);
            short status = (short) Math.abs(bank.getShort("status", k));
			boolean inDC = (status>=2000 && status<4000);
			float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
			float the = (float)Math.toDegrees(Math.acos(pz/mom));
			//if(pid == 211 && pip_part_ind==-1 && inDC && Math.abs(bank.getFloat("chi2pid", k))<3 ){}
			if(pid == 211 && pip_part_ind==-1 && inDC && (mom>0.5||runNum<2600) ){
				pip_part_ind = k;
				pip_mom = mom;
				pip_the = the;
				pip_phi = (float)Math.toDegrees(Math.atan2(py,px));
				pip_vx = vx;pip_vy = vy;pip_vz = vz;
				pip_beta = be;
                Vpip = new LorentzVector(px,py,pz,Math.sqrt(pip_mom*pip_mom+0.13957f*0.13957f));
			}
			if(pid == -211 && pim_part_ind==-1 && inDC && (mom>0.5||runNum<2600) ){
				pim_part_ind = k;
				pim_mom = mom;
			}
			if(pid == 2212 && prot_part_ind==-1){
				prot_part_ind = k;
				prot_mom = mom;
				prot_the = the;
				prot_phi = (float)Math.toDegrees(Math.atan2(py,px));
				prot_vx = vx;prot_vy = vy;prot_vz = vz;
				prot_beta = be;
                Vprot = new LorentzVector(px,py,pz,Math.sqrt(prot_mom*prot_mom+0.93827f*0.93827f));
			}
			//if( (pid == 22 || (runNum<2600&&bank.getByte("charge",k)==0) ) && (mom>0.6||runNum<2600) && the>6 && G1_mom < mom ){}
			if( bank.getByte("charge",k)==0 && mom>0.4 && the>6 && G1_mom < mom ){
				G1_part_ind = k;
				G1_mom = mom;
				G1_the = the;
				G1_phi = (float)Math.toDegrees(Math.atan2(py,px));
				VG1 = new LorentzVector(px,py,pz,mom);
			}
			//if( G1_part_ind>-1 && k!=G1_part_ind && (mom>0.6||runNum<2600) && the>6 && (pid == 22 || (runNum<2600&&bank.getByte("charge",k)==0)) 
					//&& G2_mom < mom && mom < G1_mom){}
			if( G1_part_ind>-1 && k!=G1_part_ind && mom>0.4 && the>6 && bank.getByte("charge",k)==0 && G2_mom < mom && mom < G1_mom){
				G2_part_ind = k;
				G2_mom = mom;
				G2_the = the;
				G2_phi = (float)Math.toDegrees(Math.atan2(py,px));
				VG2 = new LorentzVector(px,py,pz,mom);
			}
		}
	}
	
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
	
	public boolean select_epip(){
		if(goodELEC&&goodPIP) {
			epip_dPhi = pip_phi - e_phi + 180f;
			while(epip_dPhi> 180f)epip_dPhi -= 360f;
			while(epip_dPhi<-180f)epip_dPhi += 360f;
			LorentzVector VmissN = new LorentzVector(0,0,0,0);
			VmissN.add(VT);
			VmissN.add(VB);
			VmissN.sub(Ve);
			VmissN.sub(Vpip);
			epip_MM = (float)VmissN.mass2();
			neut_mom=-1f;neut_the=-1f;neut_phi=-1f;
			if(epip_MM<1.1) {
				neut_mom=(float)VmissN.p();
				neut_the=(float)Math.toDegrees(VmissN.theta());
				neut_phi=(float)Math.toDegrees(VmissN.phi());
				return goodNEUT;
			}
		}
		return false;
	}
	
	public boolean select_epi0(){
		boolean res = false;
		if( true
		  //&& ( runNum < 2600 || (Math.abs(e_FTOF1b_vt) < 0.5 && Math.abs(e_vz+5)<5) ) 
		  && G1_part_ind>-1 && G2_part_ind>-1 
		  //&& (runNum<2600 || (e_PCAL_t>190 && e_PCAL_t<210 && Math.abs(e_PCAL_t - G1_pcal_t) < 5 && Math.abs(e_PCAL_t - G2_pcal_t) < 5 ))
		  //&& Math.abs(G1_pcal_vt)<3 && Math.abs(G2_pcal_vt)<3 
		  ){
			VPI0 = new LorentzVector(0,0,0,0);
			VPI0.add(VG1);
			VPI0.add(VG2);
			pi0_mass = (float)VPI0.mass();
			pi0_E = (float)VPI0.e();
		    pi0_the = (float)Math.toDegrees(VPI0.theta());
			pi0_phi = (float)Math.toDegrees(VPI0.phi());
			pi0_open = (float)Vangle(VG1.vect(),VG2.vect()) ;
			if(runNum > 2600 && pi0_open > 3 && pi0_open>9 * (1 - pi0_E/4) && pi0_the>8 && pi0_mass>0.05 && pi0_mass<0.5)res = true;
			if(runNum < 2600)res = true;
		}
		return res;
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
	
	
    public void processEvent(DataEvent event) {
    	
    	ev.setHipoEvent(isHipo3Event);
	    ev.procEvent(event);
	   	    
	    if(!ev.countElectronTriggers(false)) return;
	    
	    this.myinit();	    	    
	    goodELEC = makeELEC();
	    goodPIP  = makePIP();
	    goodNEUT = makeNEUT();
	    goodPHOT = goodELEC&&makePHOT();
	    
	    fillHists();
	    
    }  
    
	public void fillHists(){
		
		if(e_sect>0) H_epip_W[e_sect-1].fill(epip_MM);
		
		if(select_epip()) {
			fillHelec();
			fillHepip();
			fillHneut();
			fillHphot();
		}	
	}    
	
	public void fillHelec() {
		H_e_EC_etot_p[0].fill(e_mom,e_ecal_esum/1000f);
	}
	
    public void fillHepip() {
		H_e_EC_etot_p[0].fill(pip_mom,pip_ecal_esum/1000f);
		H_epip_e_th_p.fill(e_mom,e_the);
		H_epip_p_th_p.fill(pip_mom,pip_the);
		H_epip_vz_vz.fill(e_vz,pip_vz);
		H_epip_dvz_phi.fill(pip_phi,pip_vz-e_vz);
		H_epip_dvz_theta.fill(pip_the,pip_vz-e_vz);
		H_epip_dvz_vz.fill(pip_vz,pip_vz-e_vz);
		H_epip_Dphi_phi.fill(pip_phi,epip_dPhi);
		H_epip_Dphi_theta.fill(pip_the,epip_dPhi);
		H_epip_Dphi_vz.fill(pip_vz,epip_dPhi);
		H_epip_beta_p.fill(pip_mom,pip_beta);
		//H_epip_W_theta[s].fill(e_W,e_the);
		H_epip_W_theta[e_sect-1].fill(epip_MM,e_the);
		
		if(found_eFTOF1b && pip_FTOF_pad1b>-1){
			H_epip_FTOF1b_dt_epad.fill(    e_FTOF_pad1b, pip_FTOF1b_vt - e_FTOF1b_vt);
			H_epip_FTOF1b_dt_pippad.fill(pip_FTOF_pad1b, pip_FTOF1b_vt - e_FTOF1b_vt);
		}
    	
    }
    
    public void fillHneut() {
    	
		boolean phicut = Math.abs(neut_phi-ecal_neut_phi)<30;	
	    boolean thecut = Math.abs(neut_the-ecal_neut_the)<8;
	    
	    float mass2 = neut_mom*neut_mom*(1f/(ecal_neut_beta*ecal_neut_beta)-1);
	    	    
		H_neut_e_th_p.fill(neut_mom, neut_the);
    	if(neut_mom>0.4) H_neut_phi_the_eff[0].fill(newPhi(neut_phi),neut_the);
		
		if(ecal_neut_sec!=e_sect) {
			H_neut_dth_dph.fill(neut_the-ecal_neut_the,neut_phi-ecal_neut_phi);
			if(thecut) H_neut_phi1_phi2.fill(neut_phi-ecal_neut_phi);
			if(phicut) H_neut_th1_th2.fill(neut_the-ecal_neut_the); 
			H_neut_esum_the.fill(ecal_neut_esum[0],ecal_neut_the);  
			H_neut_p_beta.fill(neut_mom, ecal_neut_beta);        
        	H_neut_rad_tail[1].fill(ecal_neut_the,e_the);
           	H_neut_mass2[1].fill(mass2);
		}
		
        if(ecal_neut_sec==e_sect) {
        	H_neut_esum_the_elec.fill(ecal_neut_esum[0],ecal_neut_the); 
        	H_neut_rad_tail[0].fill(ecal_neut_the,e_the);
        	H_neut_phi_the.fill(ecal_neut_phi,ecal_neut_the);
           	H_neut_phie_the.fill(e_phi,ecal_neut_the);
        }
        
		if (thecut && phicut) {
			H_neut_p_beta_cut.fill(neut_mom, ecal_neut_beta);	    	
        	H_neut_mass2[0].fill(mass2);
        	if(mass2>0.3&&neut_mom>0.4) {
        		H_neut_phi_the_eff[1].fill(newPhi(ecal_neut_phi),ecal_neut_the);
        		H_neut_avg_mom.fill(neut_mom);
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
    
    public void fillHphot() {
    	for (int is=0; is<6; is++) {
    		if(ecal_phot_esum[is]>0 && (is+1)==e_sect) {
    			H_phot_esum_the_elec.fill(ecal_phot_esum[is],ecal_phot_the); 
    		}
        }    	
    }

    public void elecPlot(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.divide(3,2);
		c.cd(0);c.draw(H_e_EC_etot_p[0]);   	
    }
	
	public void epipPlot(int index) {
		System.out.println("isAnalyzeDone = "+isAnalyzeDone);
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.divide(6,4);
		c.cd(0);c.draw(H_epip_e_th_p);
		c.cd(1);c.draw(H_epip_p_th_p);
		c.cd(2);c.draw(H_epip_vz_vz);
		c.cd(3);c.draw(H_epip_dvz_phi);
		c.cd(4);c.draw(H_epip_dvz_theta);
		c.cd(5);c.draw(H_epip_dvz_vz);
		c.cd(6);c.draw(H_epip_Dphi_phi);
		c.cd(7);c.draw(H_epip_Dphi_theta);
		c.cd(8);c.draw(H_epip_Dphi_vz);
		c.cd(9);c.draw(H_epip_beta_p);
//		c.cd(10);c.draw(H_epip_FTOF1b_dt_epad);
//		c.cd(11);c.draw(H_epip_FTOF1b_dt_pippad);
		for(int s=0;s<6;s++){
			H2F h2 =  H_epip_W_theta[s];
			c.cd(12+s);c.draw(h2);
			c.cd(18+s);c.draw(H_epip_W[s]); c.draw(H_epip_W_theta[s].projectionX(),"same");
//			c.getPad(18+s).getAxisX().setRange(0.25,3.25);
		}	
		c.repaint();
	}
	
	public void neutPlot(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.divide(4,3);
        F1D f1 = new F1D("neut","1/(1+[a]^2/x^2)^0.5", 0.21,2.5); 
        f1.setParameter(0,0.93957); f1.setLineColor(0); f1.setLineStyle(1); f1.setLineWidth(2);  
		
		c.cd(0);c.draw(H_neut_e_th_p);
		c.cd(1);c.draw(H_neut_dth_dph);
		c.cd(2);c.draw(H_neut_phi1_phi2);
		c.cd(3);c.draw(H_neut_th1_th2);		
		c.cd(4);c.draw(H_neut_p_beta);     c.draw(f1,"same");
		c.cd(5);c.draw(H_neut_p_beta_cut); c.draw(f1,"same");
		c.cd(6);c.draw(H_neut_mass2[1]);
		c.cd(7);c.draw(H_neut_mass2[0]);
		c.cd(8);c.draw(H_neut_phi_the_eff[0]);
		c.cd(9);c.draw(H_neut_phi_the_eff[1]);
		c.cd(10);c.draw(H_neut_avg_mom);
		c.repaint();		
	}
	
	public void photPlot(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.divide(4,3);		
		c.cd(0);c.draw(H_phot_esum_the_elec);
		c.cd(1);c.draw(H_neut_esum_the_elec);
		c.cd(2);c.draw(H_neut_esum_the);
		c.cd(3);c.draw(H_neut_rad_tail[0]);
		c.cd(4);c.draw(H_neut_rad_tail[1]);
		c.cd(5);c.draw(H_neut_phi_the);
		c.cd(6);c.draw(H_neut_phie_the);
		c.repaint();
	}
	
    @Override       
    public void plotHistos(int run) {
    	setRunNumber(run);
    	elecPlot(0);
        epipPlot(1);
        neutPlot(3);
        photPlot(4);
        neutronEff();
    }
    
    @Override
    public void plotEvent(DataEvent de) {
       analyze();
    }
    
    public void analyze() {    
        System.out.println("I am in analyze()");
        neutronEff();
        isAnalyzeDone = true;
    }
    
    public void neutronEff() {
    	float eff_rat=0, eff_err=0;
    	n_neut_ecal = 0; n_neut_mm = 0;
    	
		for(int iy=0; iy<18; iy++) { //theta bins
	    	int ix1=0, ix2=0;
			for (int ix=0; ix<H_neut_phi_the_eff[0].getDataSize(0); ix++) { //phi bins
    			if (H_neut_phi_the_eff[1].getData(ix, iy)>0) { //get phi fiducial limits for neut_mm
    				if (ix1==0) ix1=ix; //phi_min
    			    if (ix1>0 ) ix2=ix; //phi_max
    				n_neut_ecal+=(float)H_neut_phi_the_eff[1].getData(ix, iy); //accumulate efficiency numerator
    			}
    		}
    		for (int ix=ix1; ix<ix2+1; ix++) { //loop over neut_mm w/ phi fiducial limits from data
				  n_neut_mm+=(float)H_neut_phi_the_eff[0].getData(ix, iy); //accumulate efficiency denominator
    		}
    	}
    	
    	if (n_neut_mm>0&&n_neut_ecal>0) {
    	  eff_rat = n_neut_ecal/n_neut_mm; 
//    	  eff_err = (float) (eff_rat*Math.sqrt(1/n_neut_ecal+1/n_neut_mm)); //only true for uncorrelated numerator/denominator
    	  eff_err = (float) Math.sqrt((n_neut_mm-1)*eff_rat*(1-eff_rat))/(n_neut_mm-1);
    	}
    
    	
    	if(n_neut_mm>n_neut_mm_save) {
    	    System.out.println("N_ecal = "+n_neut_ecal+" N_mm = "+n_neut_mm+" EFF = "+eff_rat+" +/- "+eff_err);
    	    n_neut_mm_save = n_neut_mm;
    	}    	
    }
    
    @Override
    public void timerUpdate() {
    	neutronEff();
    	return;
    }
    
}
