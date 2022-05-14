package org.clas.tools;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import org.clas.service.ec.ECCommon;
import org.clas.service.ec.ECEngine;

public class EngineControl implements ActionListener {
	
    JTextField              pcalTF = new JTextField(4);  
    JTextField              ecinTF = new JTextField(4);  
    JTextField              ecouTF = new JTextField(4);  
    JTextField              wlogTF = new JTextField(4); 
    JLabel                  cfigLB = new JLabel();
    JComboBox            configCMB = null;
    JCheckBox              debugCB = null;
    JCheckBox                engCB = null;
    JCheckBox               wlogCB = null;
    JCheckBox             repeatCB = null;
    ButtonGroup                bG1 = null;
	
	public String config,split,spthr,touch,configField;
	public int pcS,eciS,ecoS,pcP,eciP,ecoP;
	public float pcT,eciT,ecoT;
	public double wlogPar;
	public boolean debug,doEng,repeatEv,isMC;
	public ECEngine engine = null;
	
	public EngineControl(ECEngine val) {
		engine = val;
	}
	
	public void initEngineThresh() {
		
	}
	
	public void initEngine() {
		
	}
	 
	public String getConfigField() {
		return config+"+"+split+"+"+spthr+"+"+touch;
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
		model.addElement("SpThr332");
		model.addElement("SpThr333");
		model.addElement("touchID1");
		model.addElement("touchID2");
		configCMB.setModel(model);
		
		configCMB.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
            String s = (String) configCMB.getSelectedItem();
            
            switch (s) {
            case   "photon": config="phot"; break;
            case "electron": config="elec"; break;
            case     "muon": config="muon"; break;
            case   "pizero": config="pi0";  break;
            case     "test": config="test"; break;
            case     "None": config="none"; break;
            case   "Split0": split="split0";   engine.setSplitMethod(0); break;
            case   "Split1": split="split1";   engine.setSplitMethod(1); break;
            case   "Split2": split="split2";   engine.setSplitMethod(2); break;
            case "SpThr332": spthr="spthr332"; engine.setSplitThresh(3,3,2); break;
            case "SpThr333": spthr="spthr333"; engine.setSplitThresh(3,3,3); break;
            case "touchID1": touch="touchid1"; engine.setTouchID(1); break;
            case "touchID2": touch="touchid2"; engine.setTouchID(2); break;
            }
            
            initEngineThresh();
            
            pcT     = ECCommon.clusterSize[0];
            eciT    = ECCommon.clusterSize[1];
            ecoT    = ECCommon.clusterSize[2];
            wlogPar = ECCommon.logParam;
            
            split   = "Split"+ECCommon.splitMethod;
            spthr   = "SpThr"+ECCommon.splitThresh[0]+ECCommon.splitThresh[0]+ECCommon.splitThresh[0];
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
            wlogCB.setText(Double.toString(wlogPar));
            cfigLB.setText(config+"+"+split+"+"+spthr+"+"+touch);
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
		wlogCB.setActionCommand("WLOG"); wlogCB.addActionListener(this); wlogCB.setText(Double.toString(wlogPar));  
		buttonPane.add(wlogCB); 
    
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
	   
	public void actionPerformed(ActionEvent e) {
		processSelect(Integer.parseInt(bG1.getSelection().getActionCommand()),e);
	    if(e.getActionCommand().compareTo("WLOG")==0) {wlogPar = Double.valueOf(wlogTF.getText()); engine.setLogParam(wlogPar);}
	}  
}

