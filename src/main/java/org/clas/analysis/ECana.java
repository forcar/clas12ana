package org.clas.analysis;

import org.clas.viewer.DetectorMonitor;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.StatNumber;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;


public class ECana  extends DetectorMonitor {

    private final int[] npaddles = new int[]{68,62,62,36,36,36,36,36,36};
    
	int Nevts, Nelecs, Ntrigs, runNum;
	boolean[] trigger_bits;
	public float EB, Ebeam;
	public float RFtime1, RFtime2, startTime;
	public long  TriggerWord;
	public int   trig_part_ind, trig_sect, trig_track_ind;
	public int   e_part_ind, e_sect, e_track_ind, hasLTCC, ngammas;
	public int   pip_part_ind, pip_track_ind, pip_sect, pim_part_ind, pim_track_ind, pim_sect;
	public int   foundCVT, CVTcharge;
	public float e_mom, e_theta, e_phi, e_vx, e_vy, e_vz, e_Ivy, e_Ivz, e_ecal_X, e_ecal_Y, e_ecal_Z, e_ecal_E;
	public float pim_ecal_E;
	public float pip_ecal_E;
	public float e_ecal_TH[] = new float[3];
	public float e_ecal_EL[] = new float[3];
	public float pim_ecal_TH[] = new float[3];
	public float pim_ecal_EL[] = new float[3];
	public float pip_ecal_TH[] = new float[3];
	public float pip_ecal_EL[] = new float[3];
	public float iU[] = new float[3];
	public float iV[] = new float[3];
	public float iW[] = new float[3];	
	public float pim_iU[] = new float[3];
	public float pim_iV[] = new float[3];
	public float pim_iW[] = new float[3];	
	public float pip_iU[] = new float[3];
	public float pip_iV[] = new float[3];
	public float pip_iW[] = new float[3];	
	public float x_ecal[] = new float[3];
	public float y_ecal[] = new float[3];
	public float pim_x_ecal[] = new float[3];
	public float pim_y_ecal[] = new float[3];
	public float pip_x_ecal[] = new float[3];
	public float pip_y_ecal[] = new float[3];
	public float e_track_chi2, e_vert_time, e_vert_time_RF, e_Q2, e_xB, e_W;
	public float e_HTCC, e_LTCC, e_pcal_e, e_etot_e, e_TOF_X, e_TOF_Y, e_TOF_Z, e_HTCC_X, e_HTCC_Y, e_HTCC_Z;
	public float e_DCR1_X, e_DCR1_Y, e_DCR1_Z, e_DCR2_X, e_DCR2_Y, e_DCR2_Z, e_DCR3_X, e_DCR3_Y, e_DCR3_Z;
	public float g1_e, g1_theta, g1_phi, g2_e, g2_theta, g2_phi;
	public float pip_mom, pip_theta, pip_phi, pip_vx, pip_vy, pip_vz, pip_vert_time, pip_beta, pip_track_chi2;
	public float pim_mom, pim_theta, pim_phi, pim_vx, pim_vy, pim_vz, pim_vert_time, pim_beta, pim_track_chi2;
	public float CVT_mom, CVT_theta, CVT_phi, CVT_vz, CVT_chi2, CVT_pathlength;
	public int CVT_ndf;
	public LorentzVector VB, VT, Ve, VG1, VG2, VPI0, VPIP, VPIM;        
    
    public ECana(String name) {
        super(name);
        this.setDetectorTabNames("E/P v P","E/P v ThV","E/P v ThD","EPC/P v ThD","EECi/P v ThD","EECo/P v ThD","E/P v XY","E/P v UVW","PIM v UVW","PIP v UVW");
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init(false);
        EB=12; Ebeam = 10.6f;
      	VB = new LorentzVector(0,0,Ebeam,Ebeam);
	    VT = new LorentzVector(0,0,0,0.93827);
    }

    @Override
    public void createHistos() {
        
        this.setNumberOfEvents(0);
         
        String[] layer = new String[]{"PCAL","ECin","ECout"};
        String[] view  = new String[]{"u","v","w"};
        
        DataGroup dg = new DataGroup(3,2);
        H2F h;
        
        for(int i=1; i<7; i++) {
            h = new H2F("ep.p"+i,50,0.5,10.5, 50, 0., 0.5);
            h.setTitleX("Sector "+i+" Momentum (GeV)");
            h.setTitleY("E / P");
            h.setTitle("");
            dg.addDataSet(h, 1);
            h = new H2F("ep.thv"+i,30,6.,36., 50, 0., 0.5);
            h.setTitleX("Sector "+i+" Vertex Theta (deg)");
            h.setTitleY("E / P");
            h.setTitle("");
            dg.addDataSet(h, 1);
            h = new H2F("ep.thd"+i,48,3.,27., 50, 0., 0.5);
            h.setTitleX("Sector "+i+" Detector Theta (deg)");
            h.setTitleY("E / P");
            h.setTitle("");
            dg.addDataSet(h, 1);
           h = new H2F("ep.th0"+i,48,3.,27., 50, 0., 0.5);
            h.setTitleX("Sector "+i+" PC Theta (deg)");
            h.setTitleY("EPC / P");
            h.setTitle("");
            dg.addDataSet(h, 1);
            h = new H2F("ep.th1"+i,48,3.,27., 50, 0., 0.5);
            h.setTitleX("Sector "+i+" ECi Theta (deg)");
            h.setTitleY("EECi / P");
            h.setTitle("");
            dg.addDataSet(h, 1);
            h = new H2F("ep.th2"+i,48,3.,27., 50, 0., 0.5);
            h.setTitleX("Sector "+i+" ECo Theta (deg)");
            h.setTitleY("EECo / P");
            h.setTitle("");
            dg.addDataSet(h, 1);
        }
        
        String tit1[] = {"PCAL EVENTS","ECIN EVENTS","ECOUT EVENTS"};
        String tit2[] = {"PCAL SF","ECIN SF","ECOUT SF"};
        String tit3[] = {"SF v PCAL","SF v ECIN","SF v ECOUT"};
        
        for (int i=1; i<4; i++) {
            h = new H2F("ep_xyc_e"+i,200, -200., 200., 200, -200.,200);   h.setTitle(tit1[i-1]); dg.addDataSet(h, 1); 
            h = new H2F("ep_xyc_w"+i,200, -200., 200., 200, -200.,200);                          dg.addDataSet(h, 1); 
            h = new H2F("ep_xyc_ww"+i,200, -200., 200., 200, -200.,200);                         dg.addDataSet(h, 1); 
            h = new H2F("ep_xyc_sf"+i,200, -200., 200., 200, -200.,200);  h.setTitle(tit2[i-1]); dg.addDataSet(h, 1);  
            h = new H2F("ep_xyc_sff"+i,200, -200., 200., 200, -200.,200); h.setTitle(tit3[i-1]); dg.addDataSet(h, 1);  
        }
        
        for (int is=1; is<7; is++) {
            h = new H2F("ep_pcal_up_"+is,"ep_pcal_up_"+is,25, 0., 0.5, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL SF");
            h.setTitleY("coordU");    
            dg.addDataSet(h, 1);  
            h = new H2F("ep_pcal_vp_"+is,"ep_pcal_vp_"+is,25, 0., 0.5, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL SF");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1);            
            h = new H2F("ep_pcal_wp_"+is,"ep_pcal_wp_"+is,25, 0., 0.5, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL SF");
            h.setTitleY("coordW");  
            dg.addDataSet(h, 1);  
            
            h = new H2F("ep_ecin_up_"+is,"ep_ecin_up_"+is,25, 0., 0.3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin SF");
            h.setTitleY("coordU");        
            dg.addDataSet(h, 1);  
            h = new H2F("ep_ecin_vp_"+is,"ep_ecin_vp_"+is,25, 0., 0.3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin SF");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1); 
            h = new H2F("ep_ecin_wp_"+is,"ep_ecin_wp_"+is,25, 0., 0.3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin SF");
            h.setTitleY("coordW");  
            dg.addDataSet(h, 1);  
           
            h = new H2F("ep_ecou_up_"+is,"ep_ecou_up_"+is,25, 0., 0.1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou SF");
            h.setTitleY("coordU");        
            dg.addDataSet(h, 1);  
            h = new H2F("ep_ecou_vp_"+is,"ep_ecou_vp_"+is,25, 0., 0.1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou SF");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1);  
            h = new H2F("ep_ecou_wp_"+is,"ep_ecou_wp_"+is,25, 0., 0.1, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou SF");
            h.setTitleY("coordW");
            dg.addDataSet(h, 1); 
        }   
        
        for (int is=1; is<7; is++) {
            h = new H2F("pm_pcal_up_"+is,"pm_pcal_up_"+is,25, 0., 60., 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL PIM (MeV)");
            h.setTitleY("coordU");    
            dg.addDataSet(h, 1);  
            h = new H2F("pm_pcal_vp_"+is,"pm_pcal_vp_"+is,25, 0., 60., 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL PIM (MeV)");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1);            
            h = new H2F("pm_pcal_wp_"+is,"pm_pcal_wp_"+is,25, 0., 60., 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL PIM (MeV)");
            h.setTitleY("coordW");  
            dg.addDataSet(h, 1);  
            
            h = new H2F("pm_ecin_up_"+is,"pm_ecin_up_"+is,25, 0., 60., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin PIM (MeV)");
            h.setTitleY("coordU");        
            dg.addDataSet(h, 1);  
            h = new H2F("pm_ecin_vp_"+is,"pm_ecin_vp_"+is,25, 0., 60., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin PIM (MeV)");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1); 
            h = new H2F("pm_ecin_wp_"+is,"pm_ecin_wp_"+is,25, 0., 60., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin PIM (MeV)");
            h.setTitleY("coordW");  
            dg.addDataSet(h, 1);  
           
            h = new H2F("pm_ecou_up_"+is,"pm_ecou_up_"+is,25, 0., 80., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou PIM (MeV)");
            h.setTitleY("coordU");        
            dg.addDataSet(h, 1);  
            h = new H2F("pm_ecou_vp_"+is,"pm_ecou_vp_"+is,25, 0., 80., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou PIM (MeV)");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1);  
            h = new H2F("pm_ecou_wp_"+is,"pm_ecou_wp_"+is,25, 0., 80., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou PIM (MeV)");
            h.setTitleY("coordW");
            dg.addDataSet(h, 1); 
        }  
        
        for (int is=1; is<7; is++) {
            h = new H2F("pp_pcal_up_"+is,"pp_pcal_up_"+is,25, 0., 60., 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL PIP (MeV)");
            h.setTitleY("coordU");    
            dg.addDataSet(h, 1);  
            h = new H2F("pp_pcal_vp_"+is,"pp_pcal_vp_"+is,25, 0., 60., 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL PIP (MeV)");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1);            
            h = new H2F("pp_pcal_wp_"+is,"pp_pcal_wp_"+is,25, 0., 60., 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL PIP (MeV)");
            h.setTitleY("coordW");  
            dg.addDataSet(h, 1);  
            
            h = new H2F("pp_ecin_up_"+is,"pp_ecin_up_"+is,25, 0., 60., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin PIP (MeV)");
            h.setTitleY("coordU");        
            dg.addDataSet(h, 1);  
            h = new H2F("pp_ecin_vp_"+is,"pp_ecin_vp_"+is,25, 0., 60., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin PIP (MeV)");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1); 
            h = new H2F("pp_ecin_wp_"+is,"pp_ecin_wp_"+is,25, 0., 60., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECin PIP (MeV)");
            h.setTitleY("coordW");  
            dg.addDataSet(h, 1);  
           
            h = new H2F("pp_ecou_up_"+is,"pp_ecou_up_"+is,25, 0., 80., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou PIP (MeV)");
            h.setTitleY("coordU");        
            dg.addDataSet(h, 1);  
            h = new H2F("pp_ecou_vp_"+is,"pp_ecou_vp_"+is,25, 0., 80., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou PIP (MeV)");
            h.setTitleY("coordV");        
            dg.addDataSet(h, 1);  
            h = new H2F("pp_ecou_wp_"+is,"pp_ecou_wp_"+is,25, 0., 80., 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECou PIP (MeV)");
            h.setTitleY("coordW");
            dg.addDataSet(h, 1); 
        }         
        
        this.getDataGroup().add(dg,0,0,0);
        
    }
    
    @Override        
    public void plotHistos() {    
    	
    	   EmbeddedCanvas c = new EmbeddedCanvas();
    	    
    	   DataGroup dg = this.getDataGroup().getItem(0,0,0);
    	   
    	    c = this.getDetectorCanvas().getCanvas("E/P v P");
        c.setGridX(false); c.setGridY(false);
    	    c.divide(3, 2);
    	    
        for(int is=1; is<7; is++) {
            c.cd(is-1);
            c.getPad().getAxisZ().setLog(getLogZ());
            c.getPad().getAxisZ().setRange(0.1*zMin, 60*zMax);
            c.draw(dg.getH2F("ep.p"+is));
        }
        
	    c = this.getDetectorCanvas().getCanvas("E/P v ThV");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 2);
	    
        for(int is=1; is<7; is++) {
            c.cd(is-1);
            c.getPad().getAxisZ().setLog(getLogZ());
            c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax);
            c.draw(dg.getH2F("ep.thv"+is));
        }
        
	    c = this.getDetectorCanvas().getCanvas("E/P v ThD");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 2);
	    
        for(int is=1; is<7; is++) {
            c.cd(is-1);
            c.getPad().getAxisZ().setLog(getLogZ());
            c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax);
            c.draw(dg.getH2F("ep.thd"+is));
        }
        
	    c = this.getDetectorCanvas().getCanvas("EPC/P v ThD");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 2);
	    
        for(int is=1; is<7; is++) {
            c.cd(is-1);
            c.getPad().getAxisZ().setLog(getLogZ());
            c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax);
            c.draw(dg.getH2F("ep.th0"+is));
        }
        
	    c = this.getDetectorCanvas().getCanvas("EECi/P v ThD");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 2);
	    
        for(int is=1; is<7; is++) {
            c.cd(is-1);
            c.getPad().getAxisZ().setLog(getLogZ());
            c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax);
            c.draw(dg.getH2F("ep.th1"+is));
        }
        
	    c = this.getDetectorCanvas().getCanvas("EECo/P v ThD");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 2);
	    
        for(int is=1; is<7; is++) {
            c.cd(is-1);
            c.getPad().getAxisZ().setLog(getLogZ());
            c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax);
            c.draw(dg.getH2F("ep.th2"+is));
        }
        
	    c = this.getDetectorCanvas().getCanvas("E/P v XY");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 3);
        
	    H2F h2 ;
	    
	    c.cd(0); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_xyc_e1"); c.draw(h2);
	    c.cd(1); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_xyc_e2"); c.draw(h2);
	    c.cd(2); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_xyc_e3"); c.draw(h2);
	   
	    c.cd(3); c.getPad().getAxisZ().setRange(0.001*zMin, 0.030*zMax); h2 = dg.getH2F("ep_xyc_sf1"); c.draw(h2);
	    c.cd(4); c.getPad().getAxisZ().setRange(0.001*zMin, 0.020*zMax); h2 = dg.getH2F("ep_xyc_sf2"); c.draw(h2);
	    c.cd(5); c.getPad().getAxisZ().setRange(0.001*zMin, 0.006*zMax); h2 = dg.getH2F("ep_xyc_sf3"); c.draw(h2);
	    
	    c.cd(6); c.getPad().getAxisZ().setRange(0.001*zMin, 0.050*zMax); h2 = dg.getH2F("ep_xyc_sff1"); c.draw(h2);
	    c.cd(7); c.getPad().getAxisZ().setRange(0.001*zMin, 0.050*zMax); h2 = dg.getH2F("ep_xyc_sff2"); c.draw(h2);
	    c.cd(8); c.getPad().getAxisZ().setRange(0.001*zMin, 0.050*zMax); h2 = dg.getH2F("ep_xyc_sff3"); c.draw(h2);
	    
	    c = this.getDetectorCanvas().getCanvas("E/P v UVW");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 3);
	    
	    int s = getActiveSector();
	    c.cd(0); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_pcal_up_"+s); c.draw(h2);
	    c.cd(1); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_pcal_vp_"+s); c.draw(h2);
	    c.cd(2); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_pcal_wp_"+s); c.draw(h2);	    
	    c.cd(3); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_ecin_up_"+s); c.draw(h2);
	    c.cd(4); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_ecin_vp_"+s); c.draw(h2);
	    c.cd(5); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_ecin_wp_"+s); c.draw(h2);	    
	    c.cd(6); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_ecou_up_"+s); c.draw(h2);
	    c.cd(7); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_ecou_vp_"+s); c.draw(h2);
	    c.cd(8); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("ep_ecou_wp_"+s); c.draw(h2);	    
	    
	    c = this.getDetectorCanvas().getCanvas("PIM v UVW");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 3);
	    
	    c.cd(0); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_pcal_up_"+s); c.draw(h2);
	    c.cd(1); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_pcal_vp_"+s); c.draw(h2);
	    c.cd(2); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_pcal_wp_"+s); c.draw(h2);	    
	    c.cd(3); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_ecin_up_"+s); c.draw(h2);
	    c.cd(4); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_ecin_vp_"+s); c.draw(h2);
	    c.cd(5); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_ecin_wp_"+s); c.draw(h2);	    
	    c.cd(6); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_ecou_up_"+s); c.draw(h2);
	    c.cd(7); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_ecou_vp_"+s); c.draw(h2);
	    c.cd(8); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pm_ecou_wp_"+s); c.draw(h2);	 
	    
	    c = this.getDetectorCanvas().getCanvas("PIP v UVW");
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 3);
	    
	    c.cd(0); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_pcal_up_"+s); c.draw(h2);
	    c.cd(1); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_pcal_vp_"+s); c.draw(h2);
	    c.cd(2); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_pcal_wp_"+s); c.draw(h2);	    
	    c.cd(3); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_ecin_up_"+s); c.draw(h2);
	    c.cd(4); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_ecin_vp_"+s); c.draw(h2);
	    c.cd(5); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_ecin_wp_"+s); c.draw(h2);	    
	    c.cd(6); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_ecou_up_"+s); c.draw(h2);
	    c.cd(7); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_ecou_vp_"+s); c.draw(h2);
	    c.cd(8); c.getPad().getAxisZ().setLog(getLogZ());   h2 = dg.getH2F("pp_ecou_wp_"+s); c.draw(h2);	    
                       
    } 
    
    @Override
    public void processEvent(DataEvent event) {
    	
    	    DataGroup dg = this.getDataGroup().getItem(0,0,0);
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[5] && super.eventResetTime_current[5] > 0){
            resetEventListener();
        }
        
		trig_part_ind=-1; e_part_ind=-1; trig_sect=0;
		e_sect=0; e_ecal_E = 0;  e_pcal_e=0; e_etot_e=0;
		trig_track_ind = -1; e_track_ind = -1; pim_track_ind = -1; pip_track_ind = -1;
		pim_sect=0 ; pip_sect=0; pip_ecal_E=0; pim_ecal_E=0; pip_part_ind = -1; pim_part_ind = -1;
		
		for (int i=0; i<3; i++) {
			e_ecal_TH[i]=0; e_ecal_EL[i]=0; x_ecal[i]=1000; y_ecal[i]=1000; iU[i]=0; iV[i]=0; iW[i]=0;
			pim_ecal_TH[i]=0; pim_ecal_EL[i]=0; pim_x_ecal[i]=1000; pim_y_ecal[i]=1000; pim_iU[i]=0; pim_iV[i]=0; pim_iW[i]=0;
			pip_ecal_TH[i]=0; pip_ecal_EL[i]=0; pip_x_ecal[i]=1000; pip_y_ecal[i]=1000; pip_iU[i]=0; pip_iV[i]=0; pip_iW[i]=0;
			
		}
		
		//TRIGGER BIT SECTOR
      	int trigger_sect = getElecTriggerSector(); 
      	
      	//HTCC*PCAL Q<0
      	if(event.hasBank("REC::Particle"))  trig_part_ind = makeTrigElectron(event.getBank("REC::Particle"),event); 
      	
        //GET TB TRACK SECTOR, TRIG_TRACK_IND OF TRIG_PART_IND
		if(event.hasBank("REC::Track") && 
		   event.hasBank("TimeBasedTrkg::TBTracks")) getTrigTBTrack(event.getBank("TimeBasedTrkg::TBTracks"),event.getBank("REC::Track")); 
		
		// IS TRIG_PART_IND ID=11?
		if(event.hasBank("REC::Particle"))     e_part_ind = makeElectron(event.getBank("REC::Particle")); 
		
		// FIND PION FOR COMPARISON
		if(event.hasBank("REC::Particle")) makePiPlusPimPID(event.getBank("REC::Particle"));

		if(e_part_ind==-1)return;	
	
		// GET TRACK INDEX of E_PART_IND
		if(event.hasBank("REC::Track")) fillEBTrack(event.getBank("REC::Track"));	
		
        LorentzVector VGS = new LorentzVector(0,0,0,0);
        VGS.add(VB);
        VGS.sub(Ve);
        e_Q2 = (float) -VGS.mass2();
        e_xB = e_Q2/(2f*0.93827f*(Ebeam-e_mom));
        e_W  = (float) Math.sqrt(0.93827f*0.93827f + e_Q2*(1f/e_xB-1f) );
        
        // GET PCAL, ECin, ECou ENERGY ETC. OF E_PART_IND   
		if(event.hasBank("REC::Calorimeter")&&event.hasBank("ECAL::clusters")) {
			getElecEBECal(event.getBank("REC::Calorimeter"),event.getBank("ECAL::clusters"));
			getPiEBECal(event.getBank("REC::Calorimeter"),event.getBank("ECAL::clusters"));
		}
		
		// GET SECTOR OF TRACK INDEX
		if(event.hasBank("TimeBasedTrkg::TBTracks")) getTBTrack(event.getBank("TimeBasedTrkg::TBTracks"));
		
//        System.out.println(e_Q2+" "+e_W+" "+e_mom+" "+e_ecal_E+" "+trig_track_ind+" "+e_sect);
		
        float  sf = e_ecal_E/e_mom;
        float sf0 = e_ecal_EL[0]/e_mom;
        float sf1 = e_ecal_EL[1]/e_mom;
        float sf2 = e_ecal_EL[2]/e_mom;
        
		if(e_mom>Ebeam*0.02 && sf > 0.02 && trig_track_ind>-1 && e_sect==trig_sect){
			if(e_sect>0&&e_sect<7){            
	              dg.getH2F("ep.p"+e_sect).fill(e_mom,sf);
	              dg.getH2F("ep.thv"+e_sect).fill(e_theta,sf);
	              dg.getH2F("ep.thd"+e_sect).fill(e_ecal_TH[0],sf);
	              dg.getH2F("ep.th0"+e_sect).fill(e_ecal_TH[0],sf0);
	              dg.getH2F("ep.th1"+e_sect).fill(e_ecal_TH[1],sf1);
	              dg.getH2F("ep.th2"+e_sect).fill(e_ecal_TH[2],sf2);	              
                  dg.getH2F("ep_pcal_up_"+e_sect).fill(sf0<0.5?sf0:0., iU[0]);
                  dg.getH2F("ep_pcal_vp_"+e_sect).fill(sf0<0.5?sf0:0., iV[0]);
                  dg.getH2F("ep_pcal_wp_"+e_sect).fill(sf0<0.5?sf0:0., iW[0]);
                  dg.getH2F("ep_ecin_up_"+e_sect).fill(sf1<0.5?sf1:0., iU[1]);
                  dg.getH2F("ep_ecin_vp_"+e_sect).fill(sf1<0.5?sf1:0., iV[1]);
                  dg.getH2F("ep_ecin_wp_"+e_sect).fill(sf1<0.5?sf1:0., iW[1]);
                  dg.getH2F("ep_ecou_up_"+e_sect).fill(sf2<0.5?sf2:0., iU[2]);
                  dg.getH2F("ep_ecou_vp_"+e_sect).fill(sf2<0.5?sf2:0., iV[2]);
                  dg.getH2F("ep_ecou_wp_"+e_sect).fill(sf2<0.5?sf2:0., iW[2]);
            }
			if(pim_sect>0&&pim_sect<7){                          
                dg.getH2F("pm_pcal_up_"+pim_sect).fill(pim_ecal_EL[0], pim_iU[0]);
                dg.getH2F("pm_pcal_vp_"+pim_sect).fill(pim_ecal_EL[0], pim_iV[0]);
                dg.getH2F("pm_pcal_wp_"+pim_sect).fill(pim_ecal_EL[0], pim_iW[0]);
                dg.getH2F("pm_ecin_up_"+pim_sect).fill(pim_ecal_EL[1], pim_iU[1]);
                dg.getH2F("pm_ecin_vp_"+pim_sect).fill(pim_ecal_EL[1], pim_iV[1]);
                dg.getH2F("pm_ecin_wp_"+pim_sect).fill(pim_ecal_EL[1], pim_iW[1]);
                dg.getH2F("pm_ecou_up_"+pim_sect).fill(pim_ecal_EL[2], pim_iU[2]);
                dg.getH2F("pm_ecou_vp_"+pim_sect).fill(pim_ecal_EL[2], pim_iV[2]);
                dg.getH2F("pm_ecou_wp_"+pim_sect).fill(pim_ecal_EL[2], pim_iW[2]);
          }
			if(pip_sect>0&&pip_sect<7){                          
                dg.getH2F("pp_pcal_up_"+pip_sect).fill(pip_ecal_EL[0], pip_iU[0]);
                dg.getH2F("pp_pcal_vp_"+pip_sect).fill(pip_ecal_EL[0], pip_iV[0]);
                dg.getH2F("pp_pcal_wp_"+pip_sect).fill(pip_ecal_EL[0], pip_iW[0]);
                dg.getH2F("pp_ecin_up_"+pip_sect).fill(pip_ecal_EL[1], pip_iU[1]);
                dg.getH2F("pp_ecin_vp_"+pip_sect).fill(pip_ecal_EL[1], pip_iV[1]);
                dg.getH2F("pp_ecin_wp_"+pip_sect).fill(pip_ecal_EL[1], pip_iW[1]);
                dg.getH2F("pp_ecou_up_"+pip_sect).fill(pip_ecal_EL[2], pip_iU[2]);
                dg.getH2F("pp_ecou_vp_"+pip_sect).fill(pip_ecal_EL[2], pip_iV[2]);
                dg.getH2F("pp_ecou_wp_"+pip_sect).fill(pip_ecal_EL[2], pip_iW[2]);
          }
			
			dg.getH2F("ep_xyc_e1").fill(-x_ecal[0], y_ecal[0],1f);
			dg.getH2F("ep_xyc_w1").fill(-x_ecal[0], y_ecal[0],sf0<0.5?sf0:0.);
			dg.getH2F("ep_xyc_e2").fill(-x_ecal[1], y_ecal[1],1f);
			dg.getH2F("ep_xyc_w2").fill(-x_ecal[1], y_ecal[1],sf1<0.5?sf1:0.);
			dg.getH2F("ep_xyc_e3").fill(-x_ecal[2], y_ecal[2],1f);
			dg.getH2F("ep_xyc_w3").fill(-x_ecal[2], y_ecal[2],sf2<0.5?sf2:0.);
			
			dg.getH2F("ep_xyc_ww1").fill(-x_ecal[0], y_ecal[0],sf<0.5?sf:0.);
			dg.getH2F("ep_xyc_ww2").fill(-x_ecal[1], y_ecal[1],sf<0.5?sf:0.);
			dg.getH2F("ep_xyc_ww3").fill(-x_ecal[2], y_ecal[2],sf<0.5?sf:0.);
		}
       
    }
    
	public int makePiPlusPimPID(DataBank bank){
		
		boolean foundelec = false;
		boolean goodbeta = false;
		int npositives = 0;
		int nnegatives = 0;
		float mybetap = 0;
		float mybetan = 0;
		
		for(int k = 0; k < bank.rows(); k++){
			int        pid = bank.getInt("pid", k);
			byte         q = bank.getByte("charge", k);
			float thisbeta = bank.getFloat("beta", k);
			
			if(pid==11) foundelec=true;
			goodbeta  = thisbeta>0;
			
			if(q<0 && goodbeta) nnegatives++;
			if(npositives==0 && q>0 && goodbeta) mybetap=thisbeta;
			if(pid!=11 && nnegatives<2 && q<0 && goodbeta) mybetan=thisbeta;
			if(q>0 && goodbeta)npositives++;
		}
		
		//if(foundelec && nnegatives==2 && npositives==1)System.out.println(foundelec+" , "+nnegatives+" , "+npositives+" , "+mybetap+" , "+mybetan);
		
		if(foundelec && nnegatives==2 && npositives==1 && mybetap>0 && mybetan>0){}
		
		if(foundelec){
			
			for(int k = 0; k < bank.rows(); k++){
				int pid = bank.getInt("pid", k);
				byte q = bank.getByte("charge", k);
				
				if(q>0){
					float px = bank.getFloat("px", k);
					float py = bank.getFloat("py", k);
					float pz = bank.getFloat("pz", k);
					pip_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
					pip_theta = (float)Math.toDegrees(Math.acos(pz/pip_mom));
					pip_phi = (float)Math.toDegrees(Math.atan2(py,px));
					pip_vx = bank.getFloat("vx", k);
					pip_vy = bank.getFloat("vy", k);
					pip_vz = bank.getFloat("vz", k);
					pip_beta = bank.getFloat("beta", k);
					if(pip_mom>0.4 && pip_theta<40 && pip_theta>6 && pip_beta>0){
						VPIP = new LorentzVector(px,py,pz,Math.sqrt(pip_mom*pip_mom+0.139*0.139));
						pip_part_ind = k;
					}
				}
				
				if(q<0 && pid!=11){
					float px = bank.getFloat("px", k);
					float py = bank.getFloat("py", k);
					float pz = bank.getFloat("pz", k);
					pim_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
					pim_theta = (float)Math.toDegrees(Math.acos(pz/pim_mom));
					pim_phi = (float)Math.toDegrees(Math.atan2(py,px));
					pim_vx = bank.getFloat("vx", k);
					pim_vy = bank.getFloat("vy", k);
					pim_vz = bank.getFloat("vz", k);
					pim_beta = bank.getFloat("beta", k);
					//System.out.println(pim_mom+" , "+pim_theta+" , "+pim_beta);
					if(pim_mom>0.4 && pim_theta<40 && pim_theta>6 && pim_beta>0){
						VPIM = new LorentzVector(px,py,pz,Math.sqrt(pim_mom*pim_mom+0.139*0.139));
						pim_part_ind = k;
					}
				}
			}
		}
		//if(pim_part_ind>-1)System.out.println("DEBUG PIMPIP part_ind : "+pim_part_ind+" , "+pip_part_ind);
		return -1;
	}
    
	public int makeElectron(DataBank bank){
		for(int k = 0; k < bank.rows(); k++){
			int pid = bank.getInt("pid", k);
			byte q = bank.getByte("charge", k);
			float px = bank.getFloat("px", k);
			float py = bank.getFloat("py", k);
			float pz = bank.getFloat("pz", k);
			e_mom = (float)Math.sqrt(px*px+py*py+pz*pz);			
			e_theta = (float)Math.toDegrees(Math.acos(pz/e_mom));
			e_vz = bank.getFloat("vz", k);
			
			if( pid == 11 && e_mom>0.02*Ebeam && e_mom<EB && e_theta>6 && Math.abs(e_vz)<200 ){
				e_phi = (float)Math.toDegrees(Math.atan2(py,px));
				e_vx = bank.getFloat("vx", k);
				e_vy = bank.getFloat("vy", k);
				Ve = new LorentzVector(px,py,pz,e_mom);
				return k;
			}
		}
		return -1;
	}
	
	public int makeTrigElectron(DataBank bank, DataEvent event){
		
		for(int k = 0; k < bank.rows(); k++){
			int  pid = bank.getInt("pid", k);
			byte   q = bank.getByte("charge", k);
			float px = bank.getFloat("px", k);
			float py = bank.getFloat("py", k);
			float pz = bank.getFloat("pz", k);
			e_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
			e_theta = (float)Math.toDegrees(Math.acos(pz/e_mom));
			e_vz = bank.getFloat("vz", k);
			
			if( pid == 11 && e_mom>0.02*Ebeam && e_mom<EB && e_theta>6 && Math.abs(e_vz)<20 ){}
			
			if( q<0 && e_theta>6 ) {
				float e_ecal_E=0;
				if(event.hasBank("REC::Calorimeter")) {
					DataBank ECALbank = event.getBank("REC::Calorimeter");
					for(int l = 0; l < ECALbank.rows(); l++) {
//						System.out.println(l+" "+ECALbank.getInt("layer",l)+" "+k+" "+ECALbank.getShort("pindex",l)+" "+ECALbank.getByte("sector",l));
						if(ECALbank.getShort("pindex",l)==k) {
							if(ECALbank.getInt("layer",l)==1) {
								trig_sect = ECALbank.getByte("sector",l);
							}
//							System.out.println(x);
					         e_ecal_E += ECALbank.getFloat("energy",l);
						}
					}
				}
				int HTCCnphe = 0;
				if(event.hasBank("REC::Cherenkov")){
					DataBank HTCCbank = event.getBank("REC::Cherenkov");
					for(int l = 0; l < HTCCbank.rows(); l++){
						if(HTCCbank.getShort("pindex",l)==k && HTCCbank.getInt("detector",l)==15){
							HTCCnphe = HTCCbank.getInt("nphe",l);
						}
					}
				}
//				System.out.println("TrigSect= "+trig_sect+" E= "+e_ecal_E+" "+HTCCnphe+" "+e_mom);
				if( HTCCnphe>1 && e_ecal_E/e_mom > 0.18){}
				if( HTCCnphe>1 && e_ecal_E/e_mom > 0.02){
					e_phi = (float)Math.toDegrees(Math.atan2(py,px));
					e_vx = bank.getFloat("vx", k);
					e_vy = bank.getFloat("vy", k);
					Ve = new LorentzVector(px,py,pz,e_mom);
					return k;
				}
			}
		}
		return -1;
	}	

    public int getDet(int layer) {
	    int[] il = {0,0,0,1,1,1,2,2,2}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
	    return il[layer-1];
	}	
	   
	public void getElecEBECal(DataBank bank, DataBank clust){
		
		for(int k = 0; k < bank.rows(); k++){
			   int det = bank.getInt("layer", k);
			short pind = bank.getShort("pindex",k);
			short bind = bank.getShort("index",k);
               float x = bank.getFloat("x",k);
               float y = bank.getFloat("y",k);
               float z = bank.getFloat("z",k);					         
               float r = (float) Math.sqrt(x*x+y*y+z*z);
               float e = bank.getFloat("energy",k);

			if(pind==e_part_ind){				
				 int ind = getDet(det);				
		         x_ecal[ind] = x; y_ecal[ind] = y;
		         e_ecal_TH[ind] = (float) Math.toDegrees(Math.acos(z/r));
		         e_ecal_EL[ind] += e;
                 e_ecal_E += e;
                 if(det==0) e_sect = bank.getByte("sector",k);
                 
                 iU[ind] = (clust.getInt("coordU", bind)-4)/8+1;
                 iV[ind] = (clust.getInt("coordV", bind)-4)/8+1;
                 iW[ind] = (clust.getInt("coordW", bind)-4)/8+1;        
			}

		}
	}
	public void getPiEBECal(DataBank bank, DataBank clust){
		
		for(int k = 0; k < bank.rows(); k++){
			   int det = bank.getInt("layer", k);
			short pind = bank.getShort("pindex",k);
			short bind = bank.getShort("index",k);
               float x = bank.getFloat("x",k);
               float y = bank.getFloat("y",k);
               float z = bank.getFloat("z",k);					         
               float r = (float) Math.sqrt(x*x+y*y+z*z);
               float e = bank.getFloat("energy",k)*1000;            
			if(pind==pim_part_ind){				
				 int ind = getDet(det);				
		         pim_x_ecal[ind] = x; pim_y_ecal[ind] = y;
		         pim_ecal_TH[ind] = (float) Math.toDegrees(Math.acos(z/r));
		         pim_ecal_EL[ind] += e;
                 pim_ecal_E += e;
                 if(det==0) pim_sect = bank.getByte("sector",k);
                 
                 pim_iU[ind] = (clust.getInt("coordU", bind)-4)/8+1;
                 pim_iV[ind] = (clust.getInt("coordV", bind)-4)/8+1;
                 pim_iW[ind] = (clust.getInt("coordW", bind)-4)/8+1;        
			}
			if(pind==pip_part_ind){				
				 int ind = getDet(det);				
		         pip_x_ecal[ind] = x; pim_y_ecal[ind] = y;
		         pip_ecal_TH[ind] = (float) Math.toDegrees(Math.acos(z/r));
		         pip_ecal_EL[ind] += e;
                 pip_ecal_E += e;
                 if(det==0) pip_sect = bank.getByte("sector",k);
                
                pip_iU[ind] = (clust.getInt("coordU", bind)-4)/8+1;
                pip_iV[ind] = (clust.getInt("coordV", bind)-4)/8+1;
                pip_iW[ind] = (clust.getInt("coordW", bind)-4)/8+1;        
			}

		}
	}	
	public void fillEBTrack(DataBank bank){
		e_track_ind=-1;pip_track_ind=-1;pim_track_ind=-1;
		for(int k = 0; k < bank.rows(); k++){
			short pind = bank.getShort("pindex",k);
			if(pind==e_part_ind){
				e_track_chi2 = 	bank.getFloat("chi2",k);
				e_track_ind = bank.getShort("index",k);
			}
			if(pind==pip_part_ind){
				pip_track_chi2 = bank.getFloat("chi2",k);
				pip_track_ind = bank.getShort("index",k);
				//System.out.println("fillEBTrack found pip track "+pip_track_ind);
			}
			if(pind==pim_part_ind){
				pim_track_chi2 = bank.getFloat("chi2",k);
				pim_track_ind = bank.getShort("index",k);
				//System.out.println("fillEBTrack found pim track "+pim_track_ind);
			}
		}
		//System.out.println("fillEBTrack : "+pim_part_ind+" , "+pip_part_ind+" ; "+pim_track_ind+" , "+pip_track_ind);
	}
		
	public void getTBTrack(DataBank bank){ 
		 if(e_track_ind>-1){
			 e_track_chi2 = bank.getFloat("chi2" , e_track_ind);
			 e_sect = bank.getInt("sector", e_track_ind);
		 }
		 if(pip_track_ind>-1)pip_sect = bank.getInt("sector", pip_track_ind);
		 if(pim_track_ind>-1)pim_sect = bank.getInt("sector", pim_track_ind);
	}
	
    public void getTrigTBTrack(DataBank bank, DataBank recBank){
    	
        for(int k = 0; k < bank.rows(); k++){
        	    if(recBank.getShort("pindex",k)==trig_part_ind) trig_track_ind = recBank.getShort("index",k);
        }
        
        if(trig_track_ind>-1 && trig_sect == bank.getInt("sector", trig_track_ind)){
                e_track_chi2 = bank.getFloat("chi2" , trig_track_ind);
                e_sect = bank.getInt("sector", trig_track_ind);
                e_DCR1_X = bank.getFloat("c1_x" , trig_track_ind);
                e_DCR1_Y = bank.getFloat("c1_y" , trig_track_ind);
                e_DCR1_Z = bank.getFloat("c1_z" , trig_track_ind);
                e_DCR3_X = bank.getFloat("c3_x" , trig_track_ind);
                e_DCR3_Y = bank.getFloat("c3_y" , trig_track_ind);
                e_DCR3_Z = bank.getFloat("c3_z" , trig_track_ind);
                e_DCR2_X = bank.getFloat("t1_x" , trig_track_ind);
                e_DCR2_Y = bank.getFloat("t1_y" , trig_track_ind);
                e_DCR2_Z = bank.getFloat("t1_z" , trig_track_ind);
	            Vector3 DCR1POS = new Vector3(e_DCR2_X,e_DCR2_Y,e_DCR2_Z);
	            Vector3 DCR1DIR = new Vector3(bank.getFloat("t1_px" , trig_track_ind),bank.getFloat("t1_py" , trig_track_ind),bank.getFloat("t1_pz" , trig_track_ind));
	            DCR1POS.rotateZ( -3.141597f*(e_sect-1)/3f );
	            DCR1DIR.rotateZ( -3.141597f*(e_sect-1)/3f );
	            float er1X = (float)DCR1POS.x();
	            float er1Y = (float)DCR1POS.y();
	            float er1Z = (float)DCR1POS.z();
	            float er1dX = (float)DCR1DIR.x();
	            float er1dY = (float)DCR1DIR.y();
	            float er1dZ = (float)DCR1DIR.z();
	            e_Ivy = er1Y + er1dY * (0f-er1X) / er1dX;
	            e_Ivz = er1Z + er1dZ * (0f-er1X) / er1dX;
	            float checkPh1 = (float)Math.toDegrees(Math.atan2( er1dY , er1dX ));
	            float checkTh1 = (float)Math.toDegrees(Math.acos( er1dZ / Math.sqrt( er1dX*er1dX+er1dY*er1dY+er1dZ*er1dZ ) ));
	            float checkPh2 = (float)Math.toDegrees(Math.atan2( e_Ivy-er1Y , -er1X ));
	            float checkTh2 = (float)Math.toDegrees(Math.acos( (e_Ivz-er1Z) / Math.sqrt( er1X*er1X + (e_Ivy-er1Y)*(e_Ivy-er1Y) + (e_Ivz-er1Z)*(e_Ivz-er1Z) )  ));

        }
    }
    
    @Override
    public void resetEventListener() {
        System.out.println("Resetting EC histogram");
        this.createHistos();
        this.plotHistos();
    }

    @Override
    public void timerUpdate() {
    	
        for(int i=1; i<4; i++) {
        H2F e  = this.getDataGroup().getItem(0,0,0).getH2F("ep_xyc_e"+i);
        H2F w  = this.getDataGroup().getItem(0,0,0).getH2F("ep_xyc_w"+i);
        H2F ww = this.getDataGroup().getItem(0,0,0).getH2F("ep_xyc_ww"+i);
        for(int loop = 0; loop < e.getDataBufferSize(); loop++) {
        	    float ne = e.getDataBufferBin(loop);
            if (ne>0) this.getDataGroup().getItem(0,0,0).getH2F("ep_xyc_sf"+i).setDataBufferBin(loop,w.getDataBufferBin(loop)/ne);
            if (ne>0) this.getDataGroup().getItem(0,0,0).getH2F("ep_xyc_sff"+i).setDataBufferBin(loop,ww.getDataBufferBin(loop)/ne);
        }
        }

    }
    
}
