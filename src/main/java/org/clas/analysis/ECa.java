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


public class ECa  extends DetectorMonitor {

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
   
    public ECa(String name) {
        super(name);
        this.setDetectorTabNames("E/P v P",
                                 "E/P v ThV",
                                 "E/P v ThD",
                                 "EPC/P v ThD",
                                 "EECi/P v ThD",
                                 "EECo/P v ThD",
                                 "E/P v XY",
                                 "E/P v UVW",
                                 "PIM v UVW",
                                 "PIP v UVW");
        this.useSectorButtons(true);
        this.useSliderPane(true);
        this.init();
        EB=12; Ebeam = 10.6f;
      	VB = new LorentzVector(0,0,Ebeam,Ebeam);
	    VT = new LorentzVector(0,0,0,0.93827);
    }

    @Override
    public void createHistos(int run) {        
    	    setRunNumber(run);
        this.setNumberOfEvents(0);        
        createEOPHistos(0,50,0.5,10.5,"ep_p",  " Momentum (GeV)",      "E / P");
        createEOPHistos(1,30,  6.,36.,"ep_thv"," VertexTheta (deg)",   "E / P");
        createEOPHistos(2,48,  3.,37.,"ep_thd"," Detector Theta (deg)","E / P");
        createEOPHistos(3,48,  3.,27.,"ep_th0"," PC Theta (deg)",      "EPC / P");
        createEOPHistos(4,48,  3.,27.,"ep_th1"," ECIN Theta (deg)",    "EECi / P");
        createEOPHistos(5,48,  3.,27.,"ep_th2"," ECOU Theta (deg)",    "EECo / P");
        createXYZHistos(6);
        createADCHistos(7,0.5,0.3,0.1,"SF");
        createADCHistos(8,60.,60.,80.,"PIM (MeV)");
        createADCHistos(9,60.,60.,80.,"PIP (MeV)");        
    }
    
    @Override    
    public void plotHistos(int run) {
    	    setRunNumber(run);
    	    plotEOPHistos(0);
    	    plotEOPHistos(1);
    	    plotEOPHistos(2);
    	    plotEOPHistos(3);
    	    plotEOPHistos(4);
    	    plotEOPHistos(5);
    	    plotXYZHistos(6);
    	    plotADCHistos(7);
    	    plotADCHistos(8);
    	    plotADCHistos(9);
    }
    
    public void plotADCHistos(int index) {
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(getActiveSector(),0,index,getRunNumber()));	       	
    }
    
    public void plotEOPHistos(int index) {
  	    drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(0,0,index,getRunNumber()));	       	
    }  
    
    public void plotXYZHistos(int index) {
 	    EmbeddedCanvas c = getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
 	    DataGroup dg = getDataGroup().getItem(0,0,index,getRunNumber());
 	    
        c.setGridX(false); c.setGridY(false);
	    c.divide(3, 3);
        
	    H2F h2 ;
	    
	    c.cd(0); c.getPad().getAxisZ().setLog(getLogZ());                h2 = (H2F) dg.getData(0).get(0); c.draw(h2);
	    c.cd(1); c.getPad().getAxisZ().setLog(getLogZ());                h2 = (H2F) dg.getData(1).get(0); c.draw(h2);
	    c.cd(2); c.getPad().getAxisZ().setLog(getLogZ());                h2 = (H2F) dg.getData(2).get(0); c.draw(h2);
	   
	    c.cd(3); c.getPad().getAxisZ().setRange(0.001*zMin, 0.006*zMax); h2 = (H2F) dg.getData(3).get(0); c.draw(h2);
	    c.cd(4); c.getPad().getAxisZ().setRange(0.001*zMin, 0.004*zMax); h2 = (H2F) dg.getData(4).get(0); c.draw(h2);
	    c.cd(5); c.getPad().getAxisZ().setRange(0.001*zMin, 0.001*zMax); h2 = (H2F) dg.getData(5).get(0); c.draw(h2);
	    
	    c.cd(6); c.getPad().getAxisZ().setRange(0.001*zMin, 0.010*zMax); h2 = (H2F) dg.getData(6).get(0); c.draw(h2);
	    c.cd(7); c.getPad().getAxisZ().setRange(0.001*zMin, 0.010*zMax); h2 = (H2F) dg.getData(7).get(0); c.draw(h2);
	    c.cd(8); c.getPad().getAxisZ().setRange(0.001*zMin, 0.010*zMax); h2 = (H2F) dg.getData(8).get(0); c.draw(h2);    	
    }
    
    public void createEOPHistos(int k, int xb, double x1, double x2, String txt1, String txt2, String txt3) {
    	
	    int run = getRunNumber();
        H2F h;  
        DataGroup dg = new DataGroup(3,2);
    
        for (int is=1; is<7; is++) {
            h = new H2F(txt1+"_"+is+"_"+k+"_"+run, xb, x1, x2, 50, 0., 0.5);
            h.setTitleX("Sector "+is+txt2);
            h.setTitleY(txt3);
            dg.addDataSet(h,is-1);
            this.getDataGroup().add(dg,0,0,k,run);
        }            
    } 
    
    public void createXYZHistos(int k) {
    	    
    	    int run = getRunNumber();
    	    H2F h;
        String tit1[] = {"PCAL EVENTS","ECIN EVENTS","ECOUT EVENTS"};
        String tit2[] = {"PCAL SF","ECIN SF","ECOUT SF"};
        String tit3[] = {"SF v PCAL","SF v ECIN","SF v ECOUT"};
        
        DataGroup dg1 = new DataGroup(3,3);
        DataGroup dg2 = new DataGroup(3,3);
        DataGroup dg3 = new DataGroup(3,3);
        
        for (int i=1; i<4; i++) {
            h = new H2F("ep_xyc_w"+i+"_"+run,  200, -200., 200., 200, -200.,200);                        dg2.addDataSet(h,i-1); 
            h = new H2F("ep_xyc_ww"+i+"_"+run, 200, -200., 200., 200, -200.,200);                        dg3.addDataSet(h,i-1); 
            h = new H2F("ep_xyc_e"+i+"_"+run,  200, -200., 200., 200, -200.,200); h.setTitle(tit1[i-1]); dg1.addDataSet(h,i-1); 
            h = new H2F("ep_xyc_sf"+i+"_"+run, 200, -200., 200., 200, -200.,200); h.setTitle(tit2[i-1]); dg1.addDataSet(h,i+2);  
            h = new H2F("ep_xyc_sff"+i+"_"+run,200, -200., 200., 200, -200.,200); h.setTitle(tit3[i-1]); dg1.addDataSet(h,i+5);  
        }
        this.getDataGroup().add(dg1, 0,0,k,run);
        this.getDataGroup().add(dg2, 0,1,k,run);
        this.getDataGroup().add(dg3, 0,2,k,run);
    }
    
    public void createADCHistos(int k, double x1, double x2, double x3, String txt) {
    	
	    int run = getRunNumber();
	    int nch = 25;
        H2F h;  
    
        for (int is=1; is<7; is++) {
            DataGroup dg = new DataGroup(3,3);
            h = new H2F("adc_pcal_u_"+is+"_"+k+"_"+run,"adc_pcal_u_"+is+"_"+k+"_"+run, nch, 0., x1, 68, 1., 69.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("U"); 
            dg.addDataSet(h,0);  
            h = new H2F("adc_pcal_v_"+is+"_"+k+"_"+run,"adc_pcal_v_"+is+"_"+k+"_"+run, nch, 0., x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("V");        
            dg.addDataSet(h,1);            
            h = new H2F("adc_pcal_w_"+is+"_"+k+"_"+run,"adc_pcal_w_"+is+"_"+k+"_"+run, nch, 0., x1, 62, 1., 63.);
            h.setTitleX("Sector "+is+" PCAL "+txt); h.setTitleY("W");  
            dg.addDataSet(h,2); 
        
            h = new H2F("adc_ecin_u_"+is+"_"+k+"_"+run,"adc_ecin_u_"+is+"_"+k+"_"+run, nch, 0., x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("U");    
            dg.addDataSet(h,3);  
            h = new H2F("adc_ecin_v_"+is+"_"+k+"_"+run,"adc_ecin_v_"+is+"_"+k+"_"+run, nch, 0., x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("V");        
            dg.addDataSet(h,4);            
            h = new H2F("adc_ecin_w_"+is+"_"+k+"_"+run,"adc_ecin_w_"+is+"_"+k+"_"+run, nch, 0., x2, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECIN "+txt); h.setTitleY("W");  
            dg.addDataSet(h,5); 
        
            h = new H2F("adc_ecou_u_"+is+"_"+k+"_"+run,"adc_ecou_u_"+is+"_"+k+"_"+run, nch, 0., x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("U");    
            dg.addDataSet(h,6);  
            h = new H2F("adc_ecou_v_"+is+"_"+k+"_"+run,"adc_ecou_v_"+is+"_"+k+"_"+run, nch, 0., x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("V");        
            dg.addDataSet(h,7);            
            h = new H2F("adc_ecou_w_"+is+"_"+k+"_"+run,"adc_ecou_w_"+is+"_"+k+"_"+run, nch, 0., x3, 36, 1., 37.);
            h.setTitleX("Sector "+is+" ECOU "+txt); h.setTitleY("W");  
            dg.addDataSet(h,8);   
            this.getDataGroup().add(dg,is,0,k,run);
        }            
    } 
        
    @Override
    public void processEvent(DataEvent event) {
    	
    	    int run = getRunNumber();
    	    DataGroup dg = this.getDataGroup().getItem(0,0,0,run);
        
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
        float[] sff = new float[3];
        
        sff[0] = e_ecal_EL[0]/e_mom;
        sff[1] = e_ecal_EL[1]/e_mom;
        sff[2] = e_ecal_EL[2]/e_mom;
        
        H2F h; 
        
        boolean good_e = e_sect>0&&e_sect<7, good_pim = pim_sect>0&&pim_sect<7, good_pip = pip_sect>0&&pip_sect<7 ; 
        
		if(e_mom>Ebeam*0.02 && sf > 0.02 && trig_track_ind>-1 && e_sect==trig_sect){
			if(good_e){
                h = (H2F) this.getDataGroup().getItem(0,0,0,run).getData(e_sect-1).get(0); h.fill(e_mom,sf);
                h = (H2F) this.getDataGroup().getItem(0,0,1,run).getData(e_sect-1).get(0); h.fill(e_theta,sf);
                h = (H2F) this.getDataGroup().getItem(0,0,2,run).getData(e_sect-1).get(0); h.fill(e_ecal_TH[0],sf);
			}
			for (int i=0; i<3; i++) {
				if(good_e){
					h = (H2F) this.getDataGroup().getItem(0,0,3+i,run).getData(e_sect-1).get(0);     h.fill(e_ecal_TH[i],sff[i]);
					h = (H2F) this.getDataGroup().getItem(0,0,6,run).getData(i).get(0); h.fill(-x_ecal[i], y_ecal[i],sff[i]<0.5?1f:0);
					h = (H2F) this.getDataGroup().getItem(0,1,6,run).getData(i).get(0); h.fill(-x_ecal[i], y_ecal[i],sff[i]<0.5?sff[i]:0.);
					h = (H2F) this.getDataGroup().getItem(0,2,6,run).getData(i).get(0); h.fill(-x_ecal[i], y_ecal[i],sf<0.5?sf:0.);
					h = (H2F) this.getDataGroup().getItem(e_sect,0,7,run).getData(3*i+0).get(0); h.fill(sff[i]<0.5?sff[i]:0., iU[i]);
					h = (H2F) this.getDataGroup().getItem(e_sect,0,7,run).getData(3*i+1).get(0); h.fill(sff[i]<0.5?sff[i]:0., iV[i]);
					h = (H2F) this.getDataGroup().getItem(e_sect,0,7,run).getData(3*i+2).get(0); h.fill(sff[i]<0.5?sff[i]:0., iW[i]);				  
				}
				if(good_pim) {
                    h = (H2F) this.getDataGroup().getItem(pim_sect,0,8,run).getData(3*i+0).get(0); h.fill(pim_ecal_EL[i], pim_iU[i]);
                    h = (H2F) this.getDataGroup().getItem(pim_sect,0,8,run).getData(3*i+1).get(0); h.fill(pim_ecal_EL[i], pim_iV[i]);
                    h = (H2F) this.getDataGroup().getItem(pim_sect,0,8,run).getData(3*i+2).get(0); h.fill(pim_ecal_EL[i], pim_iW[i]);				  
				}
				if(good_pip) {
                    h = (H2F) this.getDataGroup().getItem(pip_sect,0,9,run).getData(3*i+0).get(0); h.fill(pip_ecal_EL[i], pip_iU[i]);
                    h = (H2F) this.getDataGroup().getItem(pip_sect,0,9,run).getData(3*i+1).get(0); h.fill(pip_ecal_EL[i], pip_iV[i]);
                    h = (H2F) this.getDataGroup().getItem(pip_sect,0,9,run).getData(3*i+2).get(0); h.fill(pip_ecal_EL[i], pip_iW[i]);						  
				}
			}
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
    public void timerUpdate() {    	
        for(int i=0; i<3; i++) {
        H2F e  = (H2F) this.getDataGroup().getItem(0,0,6,getRunNumber()).getData(i).get(0);
        H2F w  = (H2F) this.getDataGroup().getItem(0,1,6,getRunNumber()).getData(i).get(0);
        H2F ww = (H2F) this.getDataGroup().getItem(0,2,6,getRunNumber()).getData(i).get(0);
        for(int loop = 0; loop < e.getDataBufferSize(); loop++) {
        	    float ne = e.getDataBufferBin(loop);
            if (ne>0) {H2F h = (H2F) this.getDataGroup().getItem(0,0,6,getRunNumber()).getData(i+3).get(0); h.setDataBufferBin(loop,w.getDataBufferBin(loop)/ne);}
            if (ne>0) {H2F h = (H2F) this.getDataGroup().getItem(0,0,6,getRunNumber()).getData(i+6).get(0); h.setDataBufferBin(loop,ww.getDataBufferBin(loop)/ne);}
        }
        }
    }
    
}
