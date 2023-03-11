package org.clas.tools;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

//import org.jlab.service.ec.*;
import org.clas.service.ec.*;
import org.jlab.io.base.DataEvent;

public class EngineControl implements ActionListener {
	
    JTextField              pcalTF = new JTextField(4);  
    JTextField              ecinTF = new JTextField(4);  
    JTextField              ecouTF = new JTextField(4);  
    JTextField              wlogTF = new JTextField(4); 
    JLabel                  cfigLB = new JLabel();
    JComboBox            configCMB = null;
    JCheckBox              debugCB = null;
    public JCheckBox         engCB = null;
    JCheckBox             repeatCB = null;
    ButtonGroup                bG1 = null;
	
	public String config="phot",split="",spthr="",touch="",configField="",mcpart="pi0",asacc=" ",asa=" ";
	public String variation="rga_fall2018_bg",geomVariation="rga_fall2018",pass="pass1";
	public int pcS,eciS,ecoS,pcP,eciP,ecoP;
	public int PCTrackingPlane,ECTrackingPlane;
	public float pcT,eciT,ecoT;
	public double wlogPar=3.0;
	public boolean debug=false,doEng=false,repeatEv=false,isMC=false,dbgECEngine=false ;
	public boolean useFADCTime, useFTpcal, useUnsharedEnergy, useTWCorrections;
	public boolean useDTCorrections, usePass2Timing, usePass2Energy, useCalibPass2, outputECHITS;
	public boolean useASA1, useASA2, useASA3, useASA4, useASA5, useCCPC, useCCEC, useCC;
	
	public ECEngine engine = null;
	
    public List<ECStrip>     strips = new ArrayList<ECStrip>();
    public List<ECPeak>       peaks = new ArrayList<ECPeak>(); 
    public List<ECCluster> clusters = new ArrayList<ECCluster>();
	
    int[][] sthrMuon = {{15,15,15},{15,15,15},{15,15,15}}; //15,20,20
    int[][] sthrPhot = {{10,10,10},{9,9,9},{8,8,8}};
    int[][] sthrElec = {{10,10,10},{10,10,10},{10,10,10}};
    int[][] sthrTest = {{10,10,10},{9,9,9},{13,13,13}};
    int[][] sthrZero = {{1,1,1},{1,1,1},{1,1,1}};
    
    int[][] pthrMuon = {{15,15,15},{20,20,20},{20,20,20}};
    int[][] pthrPhot = {{18,18,18},{20,20,20},{15,15,15}};
    int[][] pthrElec = {{30,30,30},{30,30,30},{30,30,30}};
    int[][] pthrTest = {{18,18,18},{22,22,22},{15,15,15}};
    int[][] pthrZero = {{1,1,1},{1,1,1},{1,1,1}};
    
    double[] cerrMuon = {5.5,10.,10.};
    double[] cerrPhot = {7,15,20}; 
    double[] cerrElec = {10.,10.,10.};
    double[] cerrTest = {4.5,11.,13.};
	
	public EngineControl() {
		engine = new ECEngine();		
	}
	
	public void initEngine() {
	    System.out.println("*** EngineControl.initEngine:Initializing ecEngine ***");
		engine.setIsSingleThreaded(true);
	    engine.setIsMC(isMC);     
	    engine.setVariation(variation);
	    engine.setGeomVariation(geomVariation);	    
	    engine.init();
	    setEngineConfig(ECCommon.config);
	    System.out.println("isMC: "+isMC+" "+mcpart);
	    System.out.println("Configuration: "+config); 
	    System.out.println("Variation: "+variation);
	    System.out.println("GeomVariation: "+geomVariation);
	    System.out.println("SingleThreaded:"+ECCommon.isSingleThreaded);
	    System.out.println("EngineControl.initEngine complete\n");
	}	
	
	public void configEngine() {
		pcS  = getStripThr(config,0,1);
	    eciS = getStripThr(config,1,1);
	    ecoS = getStripThr(config,2,1);
	    pcP  = getPeakThr(config,0,1);
	    eciP = getPeakThr(config,1,1);
	    ecoP = getPeakThr(config,2,1);
	    pcT  = getClusterErr(config,0);
	    eciT = getClusterErr(config,1);
	    ecoT = getClusterErr(config,2); 
	    engine.setStripThresholds(pcS,eciS,ecoS); 
	    engine.setPeakThresholds( pcP,eciP,ecoP); 
	    engine.setClusterCuts(pcT,eciT,ecoT); 
	    split   = "Split"+ECCommon.splitMethod;
	    spthr   = "SpThr"+ECCommon.splitThresh[0]+ECCommon.splitThresh[1]+ECCommon.splitThresh[2]+pass;
	    touch   = "touchID"+ECCommon.touchID;
	    cfigLB.setText(getConfigField()); 		
	}
	
	public void updateConfig(String val) {
        
        switch (val) {
        case   "photon": config="phot"; mcpart=config; break;
        case "electron": config="elec"; break;
        case     "muon": config="muon"; break;
        case   "pizero": config="pi0";  mcpart=config;  break;
        case     "test": config="test"; break;
        case     "None": config="none"; break;
        case   "Split0": split="split0";   engine.setSplitMethod(0); break;
        case   "Split1": split="split1";   engine.setSplitMethod(1); break;
        case   "Split2": split="split2";   engine.setSplitMethod(2); break;
        case   "Split3": split="split3";   engine.setSplitMethod(3); break;
        case   "Split4": split="split4";   engine.setSplitMethod(4); break;
        case "SpThr332": spthr="spthr332"; engine.setSplitThresh(3,3,2); break;
        case "SpThr333": spthr="spthr333"; engine.setSplitThresh(3,3,3); break;
        case "touchID1": touch="touchid1"; engine.setTouchID(1); break;
        case "touchID2": touch="touchid2"; engine.setTouchID(2); break;
        case      "DEF": asacc=" ";        resetASACC(); setUseDEF(true);  break;   
        case     "CCEC": asacc="+ccec";    resetASACC(); setUseCCEC(true); break;
        case   "CCPCEC": asacc="+ccpcec";  resetASACC(); setUseCCPC(true); setUseCCEC(true); break;
        case     "ASA1": asacc="+asa1";    resetASACC(); setUseASA1(true); break;
        case     "ASA2": asacc="+asa2";    resetASACC(); setUseASA2(true); break;
        case     "ASA3": asacc="+asa3";    resetASACC(); setUseASA3(true); break;
        case     "ASA4": asacc="+asa4";    resetASACC(); setUseASA4(true); break;
        case     "ASA5": asacc="+asa5";    resetASACC(); setUseASA5(true); break;
        case  "ASACC 1": asacc="+asa1cc";  resetASACC(); setUseCCPC(true); setUseCCEC(true); setUseASA1(true); break;   
        case  "ASACC 2": asacc="+asa2cc";  resetASACC(); setUseCCPC(true); setUseCCEC(true); setUseASA2(true); break;   
        case  "ASACC 3": asacc="+asa3cc";  resetASACC(); setUseCCPC(true); setUseCCEC(true); setUseASA3(true); break;   
        case  "ASACC 4": asacc="+asa4cc";  resetASACC(); setUseCCPC(true); setUseCCEC(true); setUseASA4(true); break;   
        case  "ASACC 5": asacc="+asa4cc";  resetASACC(); setUseCCPC(true); setUseCCEC(true); setUseASA5(true); break;   
        case   "rga_bg": variation="rga_fall2018_bg"; engine.setVariation(variation); break;
        case  "default": variation="default";         engine.setVariation(variation);
        }		
	}
	
    public void setEngineConfig(String val) {
    	configCMB.setSelectedItem(val);
    	config = val;
    }
	
	public void configDisplay() {
        pcT     = ECCommon.clusterSize[0];
        eciT    = ECCommon.clusterSize[1];
        ecoT    = ECCommon.clusterSize[2];
        wlogPar = ECCommon.logParam;
        
        split   = "Split"+ECCommon.splitMethod;
        spthr   = "SpThr"+ECCommon.splitThresh[0]+ECCommon.splitThresh[1]+ECCommon.splitThresh[2];
        touch   = "touchID"+ECCommon.touchID;
        
        switch (Integer.parseInt(bG1.getSelection().getActionCommand())) {
        case 0:    
        pcalTF.setText(Integer.toString(pcS));
        ecinTF.setText(Integer.toString(eciS));              
        ecouTF.setText(Integer.toString(ecoS)); break;
        case 1:
        pcalTF.setText(Integer.toString(pcP));
        ecinTF.setText(Integer.toString(eciP));              
        ecouTF.setText(Integer.toString(ecoP)); break;
        case 2:
        pcalTF.setText(Double.toString(pcT));
        ecinTF.setText(Double.toString(eciT));              
        ecouTF.setText(Double.toString(ecoT)); 
        }
        
        wlogTF.setText(Double.toString(wlogPar));
        cfigLB.setText(getConfigField());		
	}
	
    public int getStripThr(String config, int idet, int layer) {
        switch (config) {
        case     "pi0": return sthrPhot[idet][layer-1] ;  
        case    "phot": return sthrPhot[idet][layer-1] ; 
        case    "muon": return sthrMuon[idet][layer-1] ;  
        case    "elec": return sthrElec[idet][layer-1] ;
        case    "test": return sthrTest[idet][layer-1] ;
        case    "none": return sthrZero[idet][layer-1] ;
        }
        return 0;
     }
    
    public int getPeakThr(String config, int idet, int layer) {
        switch (config) {
        case     "pi0": return pthrPhot[idet][layer-1] ;  
        case    "phot": return pthrPhot[idet][layer-1] ;  
        case    "muon": return pthrMuon[idet][layer-1] ; 
        case    "elec": return pthrElec[idet][layer-1] ;
        case    "test": return pthrTest[idet][layer-1] ;
        case    "none": return pthrZero[idet][layer-1] ;
        }
        return 0;
     }
    
    public float getClusterErr(String config, int idet) {
        switch (config) {
        case     "pi0": return (float) cerrPhot[idet] ;  
        case    "phot": return (float) cerrPhot[idet] ;  
        case    "muon": return (float) cerrMuon[idet] ; 
        case    "elec": return (float) cerrElec[idet] ;
        case    "test": return (float) cerrTest[idet] ;
        case    "none": return (float) cerrMuon[idet] ;
        }
        return 0;
     }		
	 
	public String getConfigField() {
		return mcpart+"    "+config+"+"+split+"+"+spthr+"+"+touch+asacc+"    "+variation+"+"+geomVariation+"+"+pass;
	}
	
    public String getConfig(int val) {
    	if(val==11)  return "elec";
    	if(val==211) return "muon";
    	if(val==22)  return "phot";
    	return config;    	
    }
    
    public void setMC(Boolean val) {
    	engine.setIsMC(val);
    	isMC = val;
    }
     
    public void setVariation(String val) {
    	engine.setVariation(val);
    	variation = val;
    }
    
    public void setGeomVariation(String val) {
    	engine.setGeomVariation(val);
    	geomVariation = val;
    }
	
    public void setUseUnsharedTime(Boolean val) {
    	engine.setUseUnsharedTime(val);
    }
    
    public void setUseUnsharedEnergy(Boolean val) {
    	engine.setUseUnsharedEnergy(val);
    	useUnsharedEnergy = val;
    }
    
    public void setUseFADCTime(Boolean val) {
    	engine.setUseFADCTime(val);
    	useFADCTime = val;
    }
    
    public void setUsePass2Timing(Boolean val) {
    	engine.setUsePass2Timing(val);
    	pass=val?"pass2":"pass1"; 
    	usePass2Timing = val;
    }
    
    public void setUsePass2Energy(Boolean val) {
    	engine.setUsePass2Energy(val);
    	pass=val?"pass2":"pass1"; 
    	usePass2Energy = val;
    } 
    
    public void setUseCalibPass2(Boolean val) {
    	engine.setUseCalibPass2(val);
    	useCalibPass2 = val;
    }
    
    public void outputECHITS(Boolean val) {
    	engine.outputECHITS(val);
    	outputECHITS = val;
    }
    
    public void resetASACC() {
    	setUseASA1(false); setUseASA2(false); setUseASA3(false); setUseASA4(false); setUseASA5(false); 
    	setUseCCEC(false); setUseCCPC(false); setUseDEF(false);
    }
    
    public void setUseDEF(Boolean val) {
    	engine.setUseDEF(val);
    }
    
    public void setUseASA1(Boolean val) {
    	engine.setUseASA1(val);
    	useASA1 = val;
    } 
    
    public void setUseASA2(Boolean val) {
    	engine.setUseASA2(val);
    	useASA2 = val;
    } 
    
    public void setUseASA3(Boolean val) {
    	engine.setUseASA3(val);
    	useASA3 = val;
    } 
    
    public void setUseASA4(Boolean val) {
    	engine.setUseASA4(val);
    	useASA4 = val;
    } 
    
    public void setUseASA5(Boolean val) {
    	engine.setUseASA5(val);
    	useASA5 = val;
    } 
    
    public void setUseCCPC(Boolean val) {
    	engine.setUseCCPC(val);
    	useCCPC = val;
    }
    
    public void setUseCCEC(Boolean val) {
    	engine.setUseCCEC(val);
    	useCCEC = val;
    }
    
    public void setUseCC(Boolean val) {
    	engine.setUseCCPC(val);
    	engine.setUseCCEC(val);
    	useCC = val;
    }  
    
    public void setUseTWcorr(Boolean val) {
    	engine.setTWCorrections(val);
    	useTWCorrections = val;
    }
    
    public void setUseDTcorr(Boolean val) {
    	engine.setDTCorrections(val);
    	useDTCorrections = val;
    }
    
    public void setUseFTpcal(Boolean val) {
    	engine.setUseFTpcal(val);
    	useFTpcal = val;
    } 
    
    public void setPCTrackingPlane(int val) {
    	engine.setPCTrackingPlane(val);
    	PCTrackingPlane = val;
    }
    
    public void setECTrackingPlane(int val) {
    	engine.setECTrackingPlane(val);
    	ECTrackingPlane = val;
    } 
    
    public void setStripThreshold(String val) {
    	String[] x  = val.split(",",3); 
    	engine.setStripThresholds(Integer.parseInt(x[0]),Integer.parseInt(x[1]),Integer.parseInt(x[2]));
    } 
    
    public void setPeakThreshold(String val) {
    	String[] x  = val.split(",",3); 
    	engine.setPeakThresholds(Integer.parseInt(x[0]),Integer.parseInt(x[1]),Integer.parseInt(x[2]));
    } 
    
    public int[] getStripThresholds() {
    	return engine.getStripThresholds();
    }
    
    public int[] getPeakThresholds() {
    	return engine.getPeakThresholds();
    }
    
    public void setLogParam(float val) {
    	engine.setLogParam(val);
    	wlogPar= val;
    }  
    
    public void setDbgECEngine(Boolean val) {
    	engine.setDebug(val);
    	dbgECEngine = val;
    }  
    
    public void update(String s) {       
        updateConfig(s);
        configEngine();
        configDisplay();     	
    }
    
	public JPanel getECEnginePane() {
    
		JPanel buttonPane = new JPanel();
		buttonPane.setLayout(new FlowLayout());      
		buttonPane.add(new JLabel("Config:"));  
		
		configCMB = new JComboBox();
		DefaultComboBoxModel model = (DefaultComboBoxModel) configCMB.getModel();
		model.addElement("photon");
		model.addElement("electron");
		model.addElement("muon");
		model.addElement("pizero");
		model.addElement("test");
		model.addElement("None");
		model.addElement("Split0");
		model.addElement("Split1");
		model.addElement("Split2");
		model.addElement("Split3");
		model.addElement("Split4");
		model.addElement("SpThr332");
		model.addElement("SpThr333");
		model.addElement("touchID1");
		model.addElement("touchID2");
		model.addElement("DEF");
		model.addElement("ASA1");
		model.addElement("ASA2");
		model.addElement("ASA3");
		model.addElement("ASA4");
		model.addElement("ASA5");
		model.addElement("CCEC");
		model.addElement("CCPCEC");
		model.addElement("ASACC 0");
		model.addElement("ASACC 1");
		model.addElement("ASACC 2");
		model.addElement("ASACC 3");
		model.addElement("rga_bg");
		model.addElement("default");
		configCMB.setModel(model);
		configCMB.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
            update((String)configCMB.getSelectedItem());           
        }
		});

		buttonPane.add(configCMB);
           
		bG1 = new ButtonGroup();
		JRadioButton b0 = new JRadioButton("Strp"); buttonPane.add(b0); b0.setActionCommand("0"); b0.addActionListener(this);
		JRadioButton b1 = new JRadioButton("Peak"); buttonPane.add(b1); b1.setActionCommand("1"); b1.addActionListener(this); 
		JRadioButton b2 = new JRadioButton("Size"); buttonPane.add(b2); b2.setActionCommand("2"); b2.addActionListener(this); 
		bG1.add(b0); bG1.add(b1); bG1.add(b2); b2.setSelected(true); 
    
		buttonPane.add(new JLabel("PC:"));       
		pcalTF.setActionCommand("PC");  pcalTF.addActionListener(this); pcalTF.setText(Double.toString(pcT));  
		buttonPane.add(pcalTF); 
    
		buttonPane.add(new JLabel("ECi:"));
		ecinTF.setActionCommand("ECI"); ecinTF.addActionListener(this); ecinTF.setText(Double.toString(eciT));  
		buttonPane.add(ecinTF); 
    
		buttonPane.add(new JLabel("ECo:"));
		ecouTF.setActionCommand("ECO"); ecouTF.addActionListener(this); ecouTF.setText(Double.toString(ecoT));  
		buttonPane.add(ecouTF); 
    
		buttonPane.add(new JLabel("WLOG:"));
		wlogTF.setActionCommand("WLOG"); 
		wlogTF.addActionListener(this); 
		wlogTF.setText(Double.toString(wlogPar));  
		buttonPane.add(wlogTF); 
    
		debugCB = new JCheckBox("Debug");
		debugCB.addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
            if(e.getStateChange() == ItemEvent.SELECTED) {
                debug = true;
                engine.setDebugSplit(true);
            } else {
                debug = false;
                engine.setDebugSplit(false);
            };
        }
		}); 
		debugCB.setSelected(false);
		buttonPane.add(debugCB);
    
		engCB = new JCheckBox("ECEngine");
		engCB.addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
            if(e.getStateChange() == ItemEvent.SELECTED) {
                doEng = true;
            } else {
                doEng = false;
            };
            initEngine();
        }
		});           
		engCB.setSelected(false);
		buttonPane.add(engCB);
    
		repeatCB = new JCheckBox("Repeat");
		repeatCB.addItemListener(new ItemListener() {
        public void itemStateChanged(ItemEvent e) {
            if(e.getStateChange() == ItemEvent.SELECTED) {
                repeatEv = true;
            } else {
                repeatEv = false;
            };
        }
		});         
		repeatCB.setSelected(false);
		buttonPane.add(repeatCB);

		buttonPane.add(cfigLB);
    
		return buttonPane;
    
	}
	
	public void processSelect(int val, ActionEvent e) {
		switch(val) {
		   case 0:
		   if(e.getActionCommand().compareTo("PC")==0) {pcS  = Integer.valueOf(pcalTF.getText());engine.setStripThresholds(pcS,eciS,ecoS);}
		   if(e.getActionCommand().compareTo("ECI")==0){eciS = Integer.valueOf(ecinTF.getText());engine.setStripThresholds(pcS,eciS,ecoS);}
		   if(e.getActionCommand().compareTo("ECO")==0){ecoS = Integer.valueOf(ecouTF.getText());engine.setStripThresholds(pcS,eciS,ecoS);} 
		   pcalTF.setText(Integer.toString(pcS)); ecinTF.setText(Integer.toString(eciS));ecouTF.setText(Integer.toString(ecoS));
		   break;
		   case 1:
		   if(e.getActionCommand().compareTo("PC")==0) {pcP  = Integer.valueOf(pcalTF.getText());engine.setPeakThresholds(pcP,eciP,ecoP);}
		   if(e.getActionCommand().compareTo("ECI")==0){eciP = Integer.valueOf(ecinTF.getText());engine.setPeakThresholds(pcP,eciP,ecoP);}
		   if(e.getActionCommand().compareTo("ECO")==0){ecoP = Integer.valueOf(ecouTF.getText());engine.setPeakThresholds(pcP,eciP,ecoP);}	   
		   pcalTF.setText(Integer.toString(pcP)); ecinTF.setText(Integer.toString(eciP));ecouTF.setText(Integer.toString(ecoP));
		   break;
		   case 2:
		   if(e.getActionCommand().compareTo("PC")==0) {pcT  = Float.valueOf(pcalTF.getText()); engine.setClusterCuts(pcT,eciT,ecoT);}
		   if(e.getActionCommand().compareTo("ECI")==0){eciT = Float.valueOf(ecinTF.getText()); engine.setClusterCuts(pcT,eciT,ecoT);}
		   if(e.getActionCommand().compareTo("ECO")==0){ecoT = Float.valueOf(ecouTF.getText()); engine.setClusterCuts(pcT,eciT,ecoT);}		   
		   pcalTF.setText(Float.toString(pcT)); ecinTF.setText(Float.toString(eciT));ecouTF.setText(Float.toString(ecoT));
	   }	   
	}
	
	public void processDataEvent(DataEvent de) {
		engine.processDataEvent(de);
		strips   = engine.getStrips();
		peaks    = engine.getPeaks();
		clusters = engine.getClusters();
	}
	   
	public void actionPerformed(ActionEvent e) {
		processSelect(Integer.parseInt(bG1.getSelection().getActionCommand()),e);
	    if(e.getActionCommand().compareTo("WLOG")==0) {wlogPar = Double.valueOf(wlogTF.getText()); engine.setLogParam(wlogPar);}
	}  
}

