package org.clas.detectors;

import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.Vector3;


public class RECmonitor extends DetectorMonitor {
    

    public RECmonitor(String name) {
        super(name);
        this.setDetectorTabNames("CVT cosmic", "CVT positive tracks", "CVT negative tracks", "DC tracks per event", "DC hits per track", "DC momentum", "DC theta angle");
        this.init(false);
    }
    
    @Override
    public void createHistos() {
        // create histograms
        this.setNumberOfEvents(0);

        H1F summary = new H1F("summary","summary",6,0.5,6.5);
        summary.setTitleX("sector");
        summary.setTitleY("REC");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setDetectorSummary(sum);
        
        H1F cos_trk_mult = new H1F("cosmic_track_multiplicity", "cosmic_track_multiplicity", 11, -0.5, 10.5);
        cos_trk_mult.setTitleX("track multiplicity");
        cos_trk_mult.setTitleY("counts");
        cos_trk_mult.setTitle("cosmic track multiplicity");

        H1F cos_trk_phi = new H1F("cosmic_track_phi", "cosmic_track_phi", 180, 0, 180);
        cos_trk_phi.setTitleX("track phi");
        cos_trk_phi.setTitleY("counts");
        cos_trk_phi.setTitle("cosmic track phi angle");
        
        H1F cos_trk_theta = new H1F("cosmic_trk_theta", "cosmic_trk_theta", 180, 0, 180);
        cos_trk_theta.setTitleX("track theta");
        cos_trk_theta.setTitleY("counts");
        cos_trk_theta.setTitle("cosmic track theta angle");

        H1F cos_trk_norm_chi2 = new H1F("cosmic_track_norm_chi2", "cosmic_track_norm_chi2", 100, 0,5);
        cos_trk_norm_chi2.setTitleX("track norm. chi2");
        cos_trk_norm_chi2.setTitleY("counts");
        cos_trk_norm_chi2.setTitle("cosmic track normalized chi2");
        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        
        H1F pos_trk_mult = new H1F("track_multiplicity_pos", "track_multiplicity_pos", 16, -0.5, 15.5);
        pos_trk_mult.setTitleX("track multiplicity");
        pos_trk_mult.setTitleY("counts");
        pos_trk_mult.setTitle("track multiplicity for q = +1");

        H1F pos_trk_p = new H1F("track_momentum_pos", "track_momentum_pos", 1060, 0.0, 10.6);
        pos_trk_p.setTitleX("track momentum [GeV]");
        pos_trk_p.setTitleY("counts");
        pos_trk_p.setTitle("track momentum for q = +1");
        
        H1F pos_trk_pt = new H1F("track_trans_mom_pos", "track_trans_mom_pos", 1060, 0.0, 10.6);
        pos_trk_pt.setTitleX("track transverse momentum pt [GeV]");
        pos_trk_pt.setTitleY("counts");
        pos_trk_pt.setTitle("track transverse momentum pt for q = +1");
        
        H1F pos_trk_phi0 = new H1F("track_phi0_pos", "track_phi0_pos", 180, 0, 180);
        pos_trk_phi0.setTitleX("track phi0 [deg]");
        pos_trk_phi0.setTitleY("counts");
        pos_trk_phi0.setTitle("track phi0 angle for q = +1");

        H1F pos_trk_z0 = new H1F("track_z0_pos", "track_z0_pos", 200, 0, 200);
        pos_trk_z0.setTitleX("track z0");
        pos_trk_z0.setTitleY("counts");
        pos_trk_z0.setTitle("track z0 for q = +1");

        H1F pos_trk_norm_chi2 = new H1F("track_norm_chi2_pos", "track_norm_chi2_pos", 100, 0,5);
        pos_trk_norm_chi2.setTitleX("track norm. chi2");
        pos_trk_norm_chi2.setTitleY("counts");
        pos_trk_norm_chi2.setTitle("track normalized chi2 for q = +1");
        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        
        H1F neg_trk_mult = new H1F("track_multiplicity_neg", "track_multiplicity_neg", 15, 0.5, 15.5);
        neg_trk_mult.setTitleX("track multiplicity");
        neg_trk_mult.setTitleY("counts");
        neg_trk_mult.setTitle("track multiplicity for q = -1");

        H1F neg_trk_p = new H1F("track_momentum_neg", "track_momentum_neg", 1060, 0.0, 10.6);
        neg_trk_p.setTitleX("track momentum [GeV]");
        neg_trk_p.setTitleY("counts");
        neg_trk_p.setTitle("track momentum for q = -1");
        
        H1F neg_trk_pt = new H1F("track_trans_mom_neg", "track_trans_mom_neg", 1060, 0.0, 10.6);
        neg_trk_pt.setTitleX("track transverse momentum pt [GeV]");
        neg_trk_pt.setTitleY("counts");
        neg_trk_pt.setTitle("track transverse momentum pt for q = -1");
        
        H1F neg_trk_phi0 = new H1F("track_phi0_neg", "track_phi0_neg", 180, 0, 180);
        neg_trk_phi0.setTitleX("track phi0 [deg]");
        neg_trk_phi0.setTitleY("counts");
        neg_trk_phi0.setTitle("track phi0 angle for q = -1");

        H1F neg_trk_z0 = new H1F("track_z0_neg", "track_z0_neg", 200, 0, 200);
        neg_trk_z0.setTitleX("track z0");
        neg_trk_z0.setTitleY("counts");
        neg_trk_z0.setTitle("track z0 for q = -1");

        H1F neg_trk_norm_chi2 = new H1F("track_norm_chi2_neg", "track_norm_chi2_neg", 100, 0,5);
        neg_trk_norm_chi2.setTitleX("track norm. chi2");
        neg_trk_norm_chi2.setTitleY("counts");
        neg_trk_norm_chi2.setTitle("track normalized chi2 for q = -1");

        
        DataGroup dg = new DataGroup(4,4);
        dg.addDataSet(cos_trk_mult, 0);
        dg.addDataSet(cos_trk_phi, 1);
        dg.addDataSet(cos_trk_theta, 2);
        dg.addDataSet(cos_trk_norm_chi2, 3);
        dg.addDataSet(pos_trk_mult, 4);
        dg.addDataSet(pos_trk_p, 5);
        dg.addDataSet(pos_trk_pt, 6);
        dg.addDataSet(pos_trk_phi0, 7);
        dg.addDataSet(pos_trk_z0, 8);
        dg.addDataSet(pos_trk_norm_chi2, 9);
        dg.addDataSet(neg_trk_mult, 10);
        dg.addDataSet(neg_trk_p, 11);
        dg.addDataSet(neg_trk_pt, 12);
        dg.addDataSet(neg_trk_phi0, 13);
        dg.addDataSet(neg_trk_z0, 14);
        dg.addDataSet(neg_trk_norm_chi2, 15);
        this.getDataGroup().add(dg,0,0,0);
        
        
        for(int sector=1; sector <= 6; sector++) {
            H1F cd_trk_ev = new H1F("cd_trk_ev_sec" + sector, "number of tracks per events  sector " + sector, 20, 0.5, 20.5);
            cd_trk_ev.setTitleX("number of tracks per event");
            cd_trk_ev.setTitleY("counts");
            cd_trk_ev.setTitle("sector "+sector);
            H1F cd_hits_per_trk = new H1F("cd_hits_per_trk_sec" + sector, "number of hits per track  sector " + sector, 150, 0.5, 150.5);
            cd_hits_per_trk.setTitleX("number of hits per track");
            cd_hits_per_trk.setTitleY("counts");
            cd_hits_per_trk.setTitle("sector "+sector);
            H1F cd_momentum = new H1F("cd_momentum_sec" + sector, "momentum distribution  sector " + sector, 106, 0.0, 10.6);
            cd_momentum.setTitleX("momentum");
            cd_momentum.setTitleY("counts");
            cd_momentum.setTitle("sector "+sector);
            H1F cd_theta = new H1F("cd_theta_sec" + sector, "theta distribution  sector " + sector, 126, 0.0, 3.15);
            cd_theta.setTitleX("theta");
            cd_theta.setTitleY("counts");
            cd_theta.setTitle("sector "+sector);
  
            dg.addDataSet(cd_trk_ev, 0);
            dg.addDataSet(cd_hits_per_trk, 1);
            dg.addDataSet(cd_momentum, 2);
            dg.addDataSet(cd_theta, 3);
            this.getDataGroup().add(dg, sector,0,0);
        }
            
        
    }
        
    @Override
    public void plotHistos() {
        // initialize canvas and plot histograms
        
        this.getDetectorCanvas().getCanvas("CVT cosmic").divide(2, 2);
        this.getDetectorCanvas().getCanvas("CVT cosmic").setGridX(false);
        this.getDetectorCanvas().getCanvas("CVT cosmic").setGridY(false);
        this.getDetectorCanvas().getCanvas("CVT cosmic").cd(0);
        this.getDetectorCanvas().getCanvas("CVT cosmic").draw(this.getDataGroup().getItem(0,0,0).getH1F("cosmic_track_multiplicity"));
        this.getDetectorCanvas().getCanvas("CVT cosmic").cd(1);
        this.getDetectorCanvas().getCanvas("CVT cosmic").draw(this.getDataGroup().getItem(0,0,0).getH1F("cosmic_track_phi"));
        this.getDetectorCanvas().getCanvas("CVT cosmic").cd(2);
        this.getDetectorCanvas().getCanvas("CVT cosmic").draw(this.getDataGroup().getItem(0,0,0).getH1F("cosmic_trk_theta"));
        this.getDetectorCanvas().getCanvas("CVT cosmic").cd(3);
        this.getDetectorCanvas().getCanvas("CVT cosmic").draw(this.getDataGroup().getItem(0,0,0).getH1F("cosmic_track_norm_chi2"));
        this.getDetectorCanvas().getCanvas("CVT cosmic").update();
        
        this.getDetectorCanvas().getCanvas("CVT positive tracks").divide(3, 2);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").setGridX(false);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").setGridY(false);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").cd(0);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_multiplicity_pos"));
        this.getDetectorCanvas().getCanvas("CVT positive tracks").cd(1);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_momentum_pos"));
        this.getDetectorCanvas().getCanvas("CVT positive tracks").cd(2);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_trans_mom_pos"));
        this.getDetectorCanvas().getCanvas("CVT positive tracks").cd(3);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_phi0_pos"));
        this.getDetectorCanvas().getCanvas("CVT positive tracks").cd(4);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_z0_pos"));
        this.getDetectorCanvas().getCanvas("CVT positive tracks").cd(5);
        this.getDetectorCanvas().getCanvas("CVT positive tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_norm_chi2_pos"));
        this.getDetectorCanvas().getCanvas("CVT positive tracks").update();
        
        this.getDetectorCanvas().getCanvas("CVT negative tracks").divide(3, 2);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").setGridX(false);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").setGridY(false);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").cd(0);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_multiplicity_neg"));
        this.getDetectorCanvas().getCanvas("CVT negative tracks").cd(1);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_momentum_neg"));
        this.getDetectorCanvas().getCanvas("CVT negative tracks").cd(2);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_trans_mom_neg"));
        this.getDetectorCanvas().getCanvas("CVT negative tracks").cd(3);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_phi0_neg"));
        this.getDetectorCanvas().getCanvas("CVT negative tracks").cd(4);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_z0_neg"));
        this.getDetectorCanvas().getCanvas("CVT negative tracks").cd(5);
        this.getDetectorCanvas().getCanvas("CVT negative tracks").draw(this.getDataGroup().getItem(0,0,0).getH1F("track_norm_chi2_neg"));
        this.getDetectorCanvas().getCanvas("CVT negative tracks").update();

        
        this.getDetectorCanvas().getCanvas("DC tracks per event").divide(2, 3);
        this.getDetectorCanvas().getCanvas("DC tracks per event").setGridX(false);
        this.getDetectorCanvas().getCanvas("DC tracks per event").setGridY(false);
        this.getDetectorCanvas().getCanvas("DC hits per track").divide(2, 3);
        this.getDetectorCanvas().getCanvas("DC hits per track").setGridX(false);
        this.getDetectorCanvas().getCanvas("DC hits per track").setGridY(false);
        this.getDetectorCanvas().getCanvas("DC momentum").divide(2, 3);
        this.getDetectorCanvas().getCanvas("DC momentum").setGridX(false);
        this.getDetectorCanvas().getCanvas("DC momentum").setGridY(false);
        this.getDetectorCanvas().getCanvas("DC theta angle").divide(2, 3);
        this.getDetectorCanvas().getCanvas("DC theta angle").setGridX(false);
        this.getDetectorCanvas().getCanvas("DC theta angle").setGridY(false);
        
        for(int sector=1; sector <=6; sector++) {
            this.getDetectorCanvas().getCanvas("DC tracks per event").getPad(sector-1).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("DC tracks per event").cd(sector-1);
            this.getDetectorCanvas().getCanvas("DC tracks per event").draw(this.getDataGroup().getItem(sector,0,0).getH1F("cd_trk_ev_sec"+sector));
            this.getDetectorCanvas().getCanvas("DC hits per track").getPad(sector-1).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("DC hits per track").cd(sector-1);
            this.getDetectorCanvas().getCanvas("DC hits per track").draw(this.getDataGroup().getItem(sector,0,0).getH1F("cd_hits_per_trk_sec" + sector));
            this.getDetectorCanvas().getCanvas("DC momentum").getPad(sector-1).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("DC momentum").cd(sector-1);
            this.getDetectorCanvas().getCanvas("DC momentum").draw(this.getDataGroup().getItem(sector,0,0).getH1F("cd_momentum_sec" + sector));
            this.getDetectorCanvas().getCanvas("DC theta angle").getPad(sector-1).getAxisZ().setLog(getLogZ());
            this.getDetectorCanvas().getCanvas("DC theta angle").cd(sector-1);
            this.getDetectorCanvas().getCanvas("DC theta angle").draw(this.getDataGroup().getItem(sector,0,0).getH1F("cd_theta_sec" + sector));
        }
    
        this.getDetectorCanvas().getCanvas("DC tracks per event").update();
        this.getDetectorCanvas().getCanvas("DC hits per track").update();
        this.getDetectorCanvas().getCanvas("DC momentum").update();
        this.getDetectorCanvas().getCanvas("DC theta angle").update();
        
        
    }

    @Override
    public void processEvent(DataEvent event) {
        
        if (this.getNumberOfEvents() >= super.eventResetTime_current[14] && super.eventResetTime_current[14] > 0){
            resetEventListener();
        }
        
		//if (!testTriggerMask()) return;
        
        // process event info and save into data group
                
        if (event.hasBank("CVTRec::Cosmics")==true){
            DataBank bankSVTCosmics = (DataBank) event.getBank("CVTRec::Cosmics");
            int nTracks_cos = bankSVTCosmics.rows();
            for (int row = 0; row < nTracks_cos; ++row) {
                float  svtCosmicTrack_phi = bankSVTCosmics.getFloat("phi", row);
                float  svtCosmicTrack_theta = bankSVTCosmics.getFloat("theta", row);
                float  svtCosmicTrack_kfChi2 = bankSVTCosmics.getFloat("chi2", row);
                int    svtCosmicTrack_kfNdf = bankSVTCosmics.getInt("ndf", row);
                
                this.getDataGroup().getItem(0,0,0).getH1F("cosmic_track_multiplicity").fill(nTracks_cos);
                if(nTracks_cos > 0){
                this.getDataGroup().getItem(0,0,0).getH1F("cosmic_track_phi").fill(svtCosmicTrack_phi);
                this.getDataGroup().getItem(0,0,0).getH1F("cosmic_trk_theta").fill(svtCosmicTrack_theta);
                if(svtCosmicTrack_kfNdf > 0)this.getDataGroup().getItem(0,0,0).getH1F("cosmic_track_norm_chi2").fill(svtCosmicTrack_kfChi2/svtCosmicTrack_kfNdf);
                }
                
            }
        } 
        
        else if (event.hasBank("CVTRec::Tracks")) {
            DataBank bankSVTTracks = event.getBank("CVTRec::Tracks");
            int nTracks = bankSVTTracks.rows();
            for (int row = 0; row < nTracks; ++row) {
                
                int    svtTrack_q = bankSVTTracks.getInt("q", row);
                float  svtTrack_p = bankSVTTracks.getFloat("p", row);
                float  svtTrack_pt = bankSVTTracks.getFloat("pt", row);
                float  svtTrack_phi0 = bankSVTTracks.getFloat("phi0", row);
                float  svtTrack_z0 = bankSVTTracks.getFloat("z0", row);
                float  svtTrack_circlefit_chi2_per_ndf = bankSVTTracks.getFloat("circlefit_chi2_per_ndf", row);

                if(svtTrack_q == +1){
                    this.getDataGroup().getItem(0,0,0).getH1F("track_multiplicity_pos").fill(nTracks);
                    if(nTracks > 0){
                        this.getDataGroup().getItem(0,0,0).getH1F("track_momentum_pos").fill(svtTrack_p);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_trans_mom_pos").fill(svtTrack_pt);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_phi0_pos").fill(svtTrack_phi0);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_z0_pos").fill(svtTrack_z0);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_norm_chi2_pos").fill(svtTrack_circlefit_chi2_per_ndf);
                    }
                }
                    
                if(svtTrack_q == -1){
                    this.getDataGroup().getItem(0,0,0).getH1F("track_multiplicity_neg").fill(nTracks);
                    if(nTracks > 0){
                        this.getDataGroup().getItem(0,0,0).getH1F("track_momentum_neg").fill(svtTrack_p);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_trans_mom_neg").fill(svtTrack_pt);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_phi0_neg").fill(svtTrack_phi0);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_z0_neg").fill(svtTrack_z0);
                        this.getDataGroup().getItem(0,0,0).getH1F("track_norm_chi2_neg").fill(svtTrack_circlefit_chi2_per_ndf);
                    }
                }
                
            }
        }
        
        else if (event.hasBank("TimeBasedTrkg::TBTracks")) {
            DataBank bankTBTracks = event.getBank("TimeBasedTrkg::TBTracks");
            int rows = bankTBTracks.rows();
            int sector;
            int ndf;
            for(int i = 0; i < rows; i++){
                sector = bankTBTracks.getByte("sector",i);
                ndf = bankTBTracks.getShort("ndf",i);

                Vector3 momentum = new Vector3();
		momentum.setXYZ(bankTBTracks.getFloat("p0_x", i), bankTBTracks.getFloat("p0_y", i),bankTBTracks.getFloat("p0_z", i));
                
                this.getDataGroup().getItem(sector,0,0).getH1F("cd_momentum_sec" + sector).fill(momentum.mag()); 
                this.getDataGroup().getItem(sector,0,0).getH1F("cd_theta_sec" + sector).fill(momentum.theta());   
                this.getDataGroup().getItem(sector,0,0).getH1F("cd_trk_ev_sec" + sector).fill(rows); 
                this.getDataGroup().getItem(sector,0,0).getH1F("cd_hits_per_trk_sec" + sector).fill(ndf+5);
                
            }
                
        } 
                        
    }
                  

    @Override
    public void timerUpdate() {
        
    }

}
