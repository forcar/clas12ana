package org.clas.analysis;

import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

public class ECperf extends DetectorMonitor {
	
	public int Nevts, Nelecs, Ntrigs, runNum;
	boolean[] trigger_bits;
	public int[] Ntrigs_sect, Nelecs_sect;

	public float RFT, STT;
	
	public LorentzVector VB, VT, Ve, VGS, Vprot, Vpip, VG1, VG2, VPI0;
	public boolean found_eTraj, found_eECAL, found_eFTOF1a, found_eFTOF1b, found_eLTCC, found_eHTCC;
	public int e_part_ind, e_sect, e_FTOF_pad1a, e_FTOF_pad1b, e_HTCC_bin_phi, e_HTCC_bin_theta;
	public float e_mom, e_the, e_phi, e_vx, e_vy, e_vz;
	public float e_xB, e_Q2, e_W;
	public float e_HTCC_tX, e_HTCC_tY, e_HTCC_tZ, e_LTCC_tX, e_LTCC_tY, e_FTOF_tX, e_FTOF_tY, e_PCAL_tX, e_PCAL_tY;
	public float e_DCSL1_tX, e_DCSL1_tY, e_DCSL2_tX, e_DCSL2_tY, e_DCSL3_tX, e_DCSL3_tY, e_DCSL4_tX, e_DCSL4_tY, e_DCSL5_tX, e_DCSL5_tY, e_DCSL6_tX, e_DCSL6_tY;
	public float e_PCAL_X, e_PCAL_Y, e_PCAL_Z, e_PCAL_edep, e_PCAL_t, e_EC_ein, e_EC_eout, e_EC_etot, e_PCAL_path, e_PCAL_vt;
	public float e_FTOF1a_X, e_FTOF1a_Y, e_FTOF1a_Z, e_FTOF1a_edep, e_FTOF1a_t, e_FTOF1a_path, e_FTOF1a_vt;
	public float e_FTOF1b_X, e_FTOF1b_Y, e_FTOF1b_Z, e_FTOF1b_edep, e_FTOF1b_t, e_FTOF1b_path, e_FTOF1b_vt;
	public float e_LTCC_X, e_LTCC_Y, e_LTCC_Z, e_LTCC_t, e_LTCC_nphe, e_LTCC_path, e_LTCC_vt;
	public float e_HTCC_X, e_HTCC_Y, e_HTCC_Z, e_HTCC_theta, e_HTCC_phi, e_HTCC_t, e_HTCC_nphe, e_HTCC_path, e_HTCC_vt;

	public int prot_part_ind;
	public float prot_mom, prot_the, prot_phi, prot_vx, prot_vy, prot_vz, prot_beta;

	public int pim_part_ind, pip_part_ind, pip_FTOF_pad1b;
	public float pip_mom, pip_the, pip_phi, pip_vx, pip_vy, pip_vz, pip_beta, pip_FTOF1b_t, pip_FTOF1b_path, pip_FTOF1b_vt;
	public float pip_FTOF1a_t, pip_FTOF1a_path, pip_FTOF1a_vt;
	public float pim_mom, pim_FTOF1a_t, pim_FTOF1a_path, pim_FTOF1a_vt;
	public float pim_FTOF1b_t, pim_FTOF1b_path, pim_FTOF1b_vt;
	public float thisTime;

	public int G1_part_ind, G2_part_ind, G1_pcal_ind, G2_pcal_ind, G1_cal_layers, G2_cal_layers;
	public float G1_mom, G1_the, G1_phi, G2_mom, G2_the, G2_phi;
	public float G1_pcal_X, G1_pcal_Y, G1_pcal_Z, G1_pcal_t, G1_pcal_R, G1_pcal_vt, G2_pcal_X, G2_pcal_Y, G2_pcal_Z, G2_pcal_t, G2_pcal_R, G2_pcal_vt;
	public float G1_ein_t, G1_eout_t, G2_ein_t, G2_eout_t;

	public float elast_dPhi, elast_EB;
	public float epip_dPhi, epip_MM;
	public float pi0_mass, pi0_E, pi0_the, pi0_phi, pi0_open;

	public H2F H_e_t_f, H_e_p_f, H_e_vz_f, H_e_vt_vz, H_e_vt_p, H_e_vt_t;
	public H2F H_e_PCAL, H_e_FTOF, H_e_LTCC, H_e_DCSL6, H_e_DCSL5, H_e_DCSL4, H_e_DCSL3, H_e_DCSL2, H_e_DCSL1, H_e_HTCC;
    public H2F H_e_nphe_HTCC, H_e_bin_theta_HTCC, H_e_bin_phi_HTCC, H_e_theta_HTCC, H_e_phi_HTCC;
	public H2F[] H_e_HTCC_cut; 
	public H2F[] H_e_t_p, H_e_vz_t, H_e_vz_p;
	public H2F[] H_e_EC_etot_p, H_e_EC_vt_theta, H_e_EC_XY;
	public H2F[] H_e_FTOF_vt_pad1a, H_e_FTOF_edep_pad1a, H_e_FTOF_XY_pad1a;
	public H2F[] H_e_FTOF_vt_pad1b, H_e_FTOF_edep_pad1b, H_e_FTOF_XY_pad1b;
	public H2F[] H_FTOF_pos_beta_mom_pad1a, H_FTOF_neg_beta_mom_pad1a, H_FTOF_pos_beta_mom_pad1b, H_FTOF_neg_beta_mom_pad1b;
	public H2F[] H_FTOF_pos_mass_mom_pad1a, H_FTOF_pos_mass_the_pad1a, H_FTOF_neg_mass_mom_pad1a, H_FTOF_neg_mass_the_pad1a;
	public H2F[] H_FTOF_pos_mass_mom_pad1b, H_FTOF_pos_mass_the_pad1b, H_FTOF_neg_mass_mom_pad1b, H_FTOF_neg_mass_the_pad1b;
	//from tof_monitor.java for timing and gain calibration, but from Dan's comment
	//use leptons/pions (both charges) for p1a, p1b and all particles (both charges) for p2.
	public H2F[] p1a_pad_vt_elec, p1a_pad_vt_pion, p1b_pad_vt_elec, p1b_pad_vt_pion, p2_pad_vt;
	public H1F[] p1a_pad_edep_elec, p1a_pad_edep_pion, p1b_pad_edep_elec, p1b_pad_edep_pion, p2_pad_edep;

	public H2F[] H_e_LTCC_vt_theta, H_e_LTCC_nphe_theta, H_e_LTCC_XY;
	public H2F[] H_e_HTCC_vt_theta, H_e_HTCC_nphe_theta, H_e_HTCC_XY;
	public H1F[][][] H_e_bin_nphe_HTCC;

	public H2F H_elast_e_th_p, H_elast_p_th_p, H_elast_vz_vz, H_elast_dvz_phi, H_elast_dvz_theta_all, H_elast_dvz_vz;
	public H2F H_elast_Dphi_phi, H_elast_Dphi_theta, H_elast_Dphi_vz, H_elast_EB_phi, H_elast_EB_theta, H_elast_EB_vz ;
	public H2F[] H_elast_W_theta, H_elast_W_Q2, H_elast_inc_W_theta, H_elast_dvz_theta;

	public H2F H_epip_e_th_p, H_epip_p_th_p, H_epip_vz_vz, H_epip_dvz_phi, H_epip_dvz_theta, H_epip_dvz_vz;
	public H2F H_epip_Dphi_phi, H_epip_Dphi_theta, H_epip_Dphi_vz, H_epip_beta_p, H_epip_FTOF1b_dt_epad, H_epip_FTOF1b_dt_pippad;
	public H2F[] H_epip_W_theta, H_epip_inc_W_theta;
	
	public H2F H_pi0_G1_XY, H_pi0_G1_TR, H_pi0_G1_vt_evt, H_pi0_G1_layer_E;//4
	public H2F H_pi0_G2_XY, H_pi0_G2_TR, H_pi0_G2_vt_evt, H_pi0_G2_layer_E;//8
	public H2F H_pi0_G1_mom_the, H_pi0_G1_phi_the, H_pi0_G2_mom_the, H_pi0_G2_phi_the;//12
	public H2F H_pi0_open_E, H_pi0_E_the, H_pi0_phi_the;//15
    public H1F H_pi0_mass, H_pi0_G1_layers, H_pi0_G2_layers;//18

	public float EB, Eb, Mp;
	
    public ECperf(String name) {
        super(name);
        this.setDetectorTabNames("ECelec",
        		                 "ECpip",
        		                 "ECpi0",
        		                 "ECneut",
        		                 "ECtime");
        
        this.init();
        this.localinit();  
    
    }
    
    public void localinit() {
    	System.out.println("ECperf.localinit()");
    	tl.setFitData(Fits);   
    	EB = 10.607f; Mp = 0.93827f;
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
    
    @Override    
    public void createHistos(int run) {   
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
    }
    
	public void resetCounters(){
        e_part_ind = -1;found_eTraj=false;found_eECAL=false;found_eFTOF1a=false;found_eFTOF1b=false;found_eLTCC=false;found_eHTCC=false;
        pim_part_ind = -1; pip_part_ind = -1;pip_FTOF_pad1b = -1;prot_part_ind = -1;G1_part_ind = -1;G2_part_ind = -1;G1_pcal_ind = -1;G2_pcal_ind = -1;
        G1_cal_layers = 0;G2_cal_layers = 0;G1_mom = 0;G2_mom = 0;
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
        for(int k = 0; k < bank.rows(); k++){
                if(bank.getShort("pindex",k)==partInd)return bank.getInt("sector",k);
        }   
        return -1; 
    } 
    
	public boolean fillConfBank(DataBank confbank){
		boolean selectTrig = false;
		long TriggerWord = confbank.getLong("trigger",0);
		for (int i = 31; i >= 0; i--) {trigger_bits[i] = (TriggerWord & (1 << i)) != 0;} 
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
		return selectTrig;
	}
	

	
	public void fillRecBank(DataBank recBank){
		STT = recBank.getFloat("STTime",0);
		RFT = recBank.getFloat("RFTime",0);
	}
	
	public void fillECAL(DataBank bank){
		e_EC_etot = 0;e_PCAL_edep=0;e_EC_ein=0;e_EC_eout=0;
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
	
    public int makeElectron(DataBank bank){
        for(int k = 0; k < bank.rows(); k++){
                int pid = bank.getInt("pid", k); 
                byte q = bank.getByte("charge", k); 
                float px = bank.getFloat("px", k); 
                float py = bank.getFloat("py", k); 
                float pz = bank.getFloat("pz", k); 
                int status = bank.getShort("status", k); 
                boolean inDC = (status>=2000 && status<4000);
                e_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
                e_the = (float)Math.toDegrees(Math.acos(pz/e_mom));
                e_vz = bank.getFloat("vz", k); 
                //if( pid == 11 && inDC && e_the>6 && Math.abs(e_vz)<200 ){}
                if( pid == 11 && inDC && Math.abs(e_vz+3)<12 && (e_mom>1.5 || runNum<2600 ) ){
                        e_phi = (float)Math.toDegrees(Math.atan2(py,px));
                        e_vx = bank.getFloat("vx", k); 
                        e_vy = bank.getFloat("vy", k); 
                        Ve = new LorentzVector(px,py,pz,e_mom);
                        VGS = new LorentzVector(0,0,0,0);
                        VGS.add(VB);
                        VGS.sub(Ve);
                        e_Q2 = (float) -VGS.mass2();
                        e_xB = e_Q2/(2f*Mp*(Eb-e_mom));
                        e_W  = (float) Math.sqrt(Mp*Mp + e_Q2*(1f/e_xB-1f) );
                        return k;
                }   
        }
        return -1;
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
			int status = bank.getShort("status", k);
			boolean inDC = (status>=2000 && status<4000);
			float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
			float the = (float)Math.toDegrees(Math.acos(pz/mom));
			//if(pid == 211 && pip_part_ind==-1 && inDC && Math.abs(bank.getFloat("chi2pid", k))<3 ){}
			if(pid == 211 && pip_part_ind==-1 && inDC && (mom>1.5||runNum<2600) ){
				pip_part_ind = k;
				pip_mom = mom;
				pip_the = the;
				pip_phi = (float)Math.toDegrees(Math.atan2(py,px));
				pip_vx = vx;pip_vy = vy;pip_vz = vz;
				pip_beta = be;
                                Vpip = new LorentzVector(px,py,pz,Math.sqrt(pip_mom*pip_mom+0.13957f*0.13957f));
			}
			if(pid == -211 && pim_part_ind==-1 && inDC && (mom>1.5||runNum<2600) ){
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
		if(pip_part_ind>-1){
			epip_dPhi = pip_phi - e_phi + 180f;
			while(epip_dPhi> 180f)epip_dPhi -= 360f;
			while(epip_dPhi<-180f)epip_dPhi += 360f;
			LorentzVector VmissN = new LorentzVector(0,0,0,0);
			VmissN.add(VT);
			VmissN.add(VB);
			VmissN.sub(Ve);
			VmissN.sub(Vpip);
			epip_MM = (float)VmissN.mass2();
			return true;
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
	
    public void processEvent(DataEvent event) {
	resetCounters();
	if(event.hasBank("RUN::config") && fillConfBank(event.getBank("RUN::config")) ){
		if(event.hasBank("REC::Event"))fillRecBank(event.getBank("REC::Event"));
		if(event.hasBank("REC::Particle"))e_part_ind = makeElectron(event.getBank("REC::Particle"));
		if(e_part_ind>-1){
			makeOthers(event.getBank("REC::Particle"));
			e_sect = 0;
			if(event.hasBank("REC::Track"))e_sect = getSect(event.getBank("REC::Track"),e_part_ind);
			if(event.hasBank("REC::Traj"))fillTraj(event.getBank("REC::Traj"));
//			if(event.hasBank("REC::Cherenkov"))fillCerenkov(event.getBank("REC::Cherenkov"));
//			if(event.hasBank("REC::Scintillator"))fillFTOF(event.getBank("REC::Particle"),event.getBank("REC::Scintillator"));
			if(event.hasBank("REC::Calorimeter"))fillECAL(event.getBank("REC::Calorimeter"));
			//if(e_sect>0 && found_eTraj){}
			if(e_sect>0 ){
				Nelecs++;Nelecs_sect[e_sect-1]++;
				FillHists();
			}
		}//if e_part_ind>-1
	}//event.hasBank("RUN::config")
    }//processEvent  
    
	public void FillHists(){
		
		int s = e_sect-1;

		if(select_epip()){
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
			H_epip_W_theta[s].fill(epip_MM,e_the);
			if(found_eFTOF1b && pip_FTOF_pad1b>-1){
				H_epip_FTOF1b_dt_epad.fill(    e_FTOF_pad1b, pip_FTOF1b_vt - e_FTOF1b_vt);
				H_epip_FTOF1b_dt_pippad.fill(pip_FTOF_pad1b, pip_FTOF1b_vt - e_FTOF1b_vt);
			}
		}
	}
	
	public void epipPlot(int index) {
        EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
		c.setSize(3600,2400);
		c.divide(6,4);
		c.setAxisTitleSize(24);
		c.setAxisFontSize(24);
		c.setTitleSize(24);
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
		c.cd(10);c.draw(H_epip_FTOF1b_dt_epad);
		c.cd(11);c.draw(H_epip_FTOF1b_dt_pippad);
		for(int s=0;s<6;s++){
			c.cd(12+s);c.draw(H_epip_W_theta[s]);
			c.cd(18+s);c.draw(H_epip_W_theta[s].projectionX());
			c.getPad(18+s).getAxisX().setRange(0.25,3.25);
		}		
	}
	
    @Override       
    public void plotHistos(int run) {
    	setRunNumber(run);
        epipPlot(1);
    }
    
}
