package org.clas.detectors;

import java.awt.BorderLayout;
import java.awt.Font;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextPane;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;


public class Acronyms  extends JTabbedPane {        
    JPanel acronymsPanel =  new JPanel(new BorderLayout());
    
    public Acronyms() {
        this.add("Acronyms", this.acronymsPanel);
        JTextPane acronymsDefinitions = new JTextPane();
        acronymsDefinitions.setText("Detectors:\n\n"
                + "CND:         Central Neutron Detector \n"
                + "CTOF:        Central Time of Flight \n"
                + "DC:            Drift Chambers \n"
                + "ECAL:         Electromagnetic Calorimeter \n"
                + "FMT:           Forward Micromegas Tracker \n"
                + "FTOF:         Forward Time of Flight \n"
                + "FTTRK:       Forward Tagger Tracker \n"
                + "FTHODO:   Forward Tagger Hodoscope \n"
                + "FTCAL:       Forward Tagger Calorimeter \n"
                + "HTTC:         High Threshold Cherenkov Counter \n"
                + "LTTC:         Low Threshold Cherenkov Counter \n"
                + "BMT:          Barrel Micromegas Tracker \n"
                + "FMT:          Forward Micromegas Tracker \n"
                + "PCal:          Preshower Calorimeter  \n"
                + "RICH:         Ring-Imaging Cherenkov Detector \n"
                + "BST:           Barrel Silicon Vertex Tracker \n\n"
                + "Components:\n\n"
                + "ADC:     Analog to Digital Converter \n"
                + "TDC:     Time to Digital Converter \n"
                + "PMT:     Photo Multiplier Tube\n"
                + "SiPM:    Silicon Photomultiplier\n");
        this.acronymsPanel.add(acronymsDefinitions,BorderLayout.CENTER);
        acronymsDefinitions.setFont(new Font("Avenir",Font.PLAIN,16));
        acronymsDefinitions.setEditable(false);
    }
      
}
