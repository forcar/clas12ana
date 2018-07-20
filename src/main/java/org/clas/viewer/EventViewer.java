package org.clas.viewer;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import javax.swing.ImageIcon;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextPane;
import javax.swing.KeyStroke;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;

import org.clas.analysis.ECa;
import org.clas.analysis.ECmip;
import org.clas.analysis.ECpi0;
import org.clas.analysis.ECt;
import org.jlab.detector.decode.CLASDecoder;
import org.jlab.detector.decode.CodaEventDecoder;
import org.jlab.detector.decode.DetectorEventDecoder;
import org.jlab.detector.view.DetectorListener;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.evio.EvioDataEvent;
import org.jlab.io.evio.EvioSource;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.task.DataSourceProcessorPane;
import org.jlab.io.task.IDataEventListener;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.elog.LogEntry; 
        
/**
 *
 * @author lcsmith
 */

public class EventViewer implements IDataEventListener, DetectorListener, ActionListener, ChangeListener {
    
    JTabbedPane                  tabbedpane = null;
    JPanel                        mainPanel = null;
    JMenuBar                        menuBar = null;
    DataSourceProcessorPane   processorPane = null;
    
    CodaEventDecoder               decoder = new CodaEventDecoder();
    CLASDecoder                clasDecoder = new CLASDecoder();
    DetectorEventDecoder   detectorDecoder = new DetectorEventDecoder();
       
    private int   canvasUpdateTime = 2000;
    private int analysisUpdateTime = 100;
    private int          runNumber = 0;
    private int        eventNumber = 0;
    private int      ccdbRunNumber = 0;
    String                 workDir = null;
   
    double PERIOD = 0;
    int     PHASE = 0;
    int    CYCLES = 0;
    
    public String outPath = "/Users/lcsmith/CLAS12ANA";
    
    DetectorMonitor[] monitors= {
//    		new ECa("ECa"),
//    		new ECt("ECt")
   		new ECmip("ECmip")
//    		new ECpi0("ECpi0")
    };
    
    Map<String,DetectorMonitor> Monitors = new LinkedHashMap<String,DetectorMonitor>();
        
//    DetectorMonitor[] monitors = null;
    
//    {
//    		new ECa("ECa")
//    		new ECt("ECt")
//    		new ECmip("ECmip")
//    		new ECpi0("ECpi0")
 //   }  ; 
        
    public EventViewer() {    	
        		
        menuBar = new JMenuBar();
        JMenuItem menuItem;
        
        JMenu file = new JMenu("File");
        file.getAccessibleContext().setAccessibleDescription("File options");
        menuItem = new JMenuItem("Load files...");
        menuItem.getAccessibleContext().setAccessibleDescription("Load files");
        menuItem.addActionListener(this);
        file.add(menuItem);        
        file.addSeparator();
        menuItem = new JMenuItem("Open histograms file");
        menuItem.getAccessibleContext().setAccessibleDescription("Open histograms file");
        menuItem.addActionListener(this);
        file.add(menuItem);
        menuItem = new JMenuItem("Save histograms to file");
        menuItem.getAccessibleContext().setAccessibleDescription("Save histograms to file");
        menuItem.addActionListener(this);
        file.add(menuItem);
        menuItem = new JMenuItem("Print histograms as png");
        menuItem.getAccessibleContext().setAccessibleDescription("Print histograms as png");
        menuItem.addActionListener(this);
        file.add(menuItem);       
        menuBar.add(file);
        
        JMenu settings = new JMenu("Settings");
        settings.getAccessibleContext().setAccessibleDescription("Choose monitoring parameters");
        menuItem = new JMenuItem("Set GUI update interval");
        menuItem.getAccessibleContext().setAccessibleDescription("Set GUI update interval");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuItem = new JMenuItem("Set global z-axis log scale");
        menuItem.getAccessibleContext().setAccessibleDescription("Set global z-axis log scale");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuItem = new JMenuItem("Set global z-axis lin scale");
        menuItem.getAccessibleContext().setAccessibleDescription("Set global z-axis lin scale");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuItem = new JMenuItem("Set run number");
        menuItem.getAccessibleContext().setAccessibleDescription("Set run number");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuBar.add(settings);
        
        JMenu reset = new JMenu("Reset");
        reset.getAccessibleContext().setAccessibleDescription("Reset histograms");
        
        JMenuItem menuItemdefault = new JMenuItem("Default for all");
        menuItemdefault.getAccessibleContext().setAccessibleDescription("Default for all");
        menuItemdefault.addActionListener(this);
        reset.add(menuItemdefault);
        
        JMenuItem menuItemdisable = new JMenuItem("Disable histogram reset");
        menuItemdisable.getAccessibleContext().setAccessibleDescription("Disable histogram reset");
        menuItemdisable.addActionListener(this);
        reset.add(menuItemdisable);
        
        menuBar.add(reset);
        
        String TriggerDef[] = { "Electron",
        		        "Electron S1","Electron S2","Electron S3","Electron S4","Electron S5","Electron S6",
        		        "HTCC S1","HTCC S2","HTCC S3","HTCC S4","HTCC S5","HTCC S6",
        		        "PCAL S1","PCAL S2","PCAL S3","PCAL S4","PCAL S5","PCAL S6",
        		        "ECAL S1","ECAL S2","ECAL S3","ECAL S4","ECAL S5","ECAL S6",
        		        "HT.PC","HT.EC","PC.EC","FTOF.PC","Unused","Unused",
        		        "1K Pulser"};
        		             
        JMenu trigBitsBeam = new JMenu("TriggerBits");
        trigBitsBeam.getAccessibleContext().setAccessibleDescription("Test Trigger Bits");
        
        for (int i=0; i<32; i++) {
        	
            JCheckBoxMenuItem bb = new JCheckBoxMenuItem(TriggerDef[i]);  
            final Integer bit = new Integer(i);
            bb.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
                	
                    if(e.getStateChange() == ItemEvent.SELECTED) {
                        for(int k=0; k<monitors.length; k++) {
                      	monitors[k].setTriggerMask(bit);
                        }
                    } else {
                        for(int k=0; k<monitors.length; k++) {
                     	monitors[k].clearTriggerMask(bit);
                        }
                    };
                }
            });         
            trigBitsBeam.add(bb); 
        	        	
        }

        menuBar.add(trigBitsBeam);
        
        String MonitorList[] = new String[monitors.length];
        for (int i=0; i<monitors.length; i++) MonitorList[i]=monitors[i].getDetectorName();
        
        JMenu monitorMenu = new JMenu("Monitors");
        monitorMenu.getAccessibleContext().setAccessibleDescription("Select active monitors");
        
        for (int i=0; i<MonitorList.length; i++) {
        	
            JCheckBoxMenuItem bb = new JCheckBoxMenuItem(MonitorList[i]);  
            bb.setActionCommand(MonitorList[i]);
                      
            bb.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {                	
                    if(e.getStateChange() == ItemEvent.SELECTED) {
                    	System.out.println(bb.getActionCommand());
//                    	String name = monitorlist[i].getDetectorName();
//                    	if(!monitor.containsKey(name)) monitor.put(monitorlist[i].getDetectorName(), monitorlist[i]);
                    } else {
//                        if monitors.containsKey(key)
                    };
                }
            });         
            bb.doClick();
            monitorMenu.add(bb);         	        	
        }
        
        menuBar.add(monitorMenu);

        // create main panel
        mainPanel = new JPanel();	
        mainPanel.setLayout(new BorderLayout());
        
      	tabbedpane 	= new JTabbedPane();

        processorPane = new DataSourceProcessorPane();
        processorPane.setUpdateRate(analysisUpdateTime);

        mainPanel.add(tabbedpane);
        mainPanel.add(processorPane,BorderLayout.PAGE_END);
           
        GStyle.getAxisAttributesX().setTitleFontSize(18);
        GStyle.getAxisAttributesX().setLabelFontSize(14);
        GStyle.getAxisAttributesY().setTitleFontSize(18);
        GStyle.getAxisAttributesY().setLabelFontSize(14);

        for(int k=0; k<this.monitors.length; k++) {
                this.tabbedpane.add(this.monitors[k].getDetectorPanel(), this.monitors[k].getDetectorName());
        	        this.monitors[k].getDetectorView().getView().addDetectorListener(this);                        
        }
        
        this.tabbedpane.addChangeListener(this);
               
        this.processorPane.addEventListener(this);
        
        this.setCanvasUpdate(canvasUpdateTime);
        
    }
      
    public void actionPerformed(ActionEvent e) {
        System.out.println(e.getActionCommand());
        if(e.getActionCommand() == "Load files...") {
            this.readFiles();
        }           
        if(e.getActionCommand()=="Set GUI update interval") {
            this.chooseUpdateInterval();
        }
        if(e.getActionCommand()=="Set global z-axis log scale") {
        	   for(int k=0; k<this.monitors.length; k++) this.monitors[k].setLogZ(true);
        }
        if(e.getActionCommand()=="Set global z-axis lin scale") {
           for(int k=0; k<this.monitors.length; k++) this.monitors[k].setLogZ(false);
        }
        if(e.getActionCommand()=="Set run number") {
           setRunNumber(e.getActionCommand());
        }

        if(e.getActionCommand()=="Open histograms file") {
            String fileName = null;
            JFileChooser fc = new JFileChooser();
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            File workingDirectory = new File(System.getProperty("user.dir"));  
            fc.setCurrentDirectory(workingDirectory);
            int option = fc.showOpenDialog(null);
            if (option == JFileChooser.APPROVE_OPTION) {
                fileName = fc.getSelectedFile().getAbsolutePath();            
            }
            if(fileName != null) this.loadHistosFromFile(fileName);
        }        
        if(e.getActionCommand()=="Print histograms as png") {
            this.printHistosToFile();
        }
        
        if(e.getActionCommand()=="Save histograms to file") {
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "CLAS12Ana_run_" + this.runNumber + "_" + df.format(new Date()) + ".hipo";
            JFileChooser fc = new JFileChooser();
            File workingDirectory = new File(System.getProperty("user.dir"));   
            fc.setCurrentDirectory(workingDirectory);
            File file = new File(fileName);
            fc.setSelectedFile(file);
            int returnValue = fc.showSaveDialog(null);
            if (returnValue == JFileChooser.APPROVE_OPTION) {
               fileName = fc.getSelectedFile().getAbsolutePath();            
            }
            this.saveHistosToFile(fileName);
        }
            
        if (e.getActionCommand()=="Default for all"){
            for (int k=0;k<monitors.length;k++){
                this.monitors[k].eventResetTime_current[k] = this.monitors[k].eventResetTime_default[k];
            }
        }
         
         if (e.getActionCommand()=="Disable histogram reset"){
            for (int k=0;k<monitors.length;k++){
                this.monitors[k].eventResetTime_current[k] = 0;
            }
        }
        
        if ( e.getActionCommand().substring(0, 5).equals("Reset")){
            resetHistograms(e.getActionCommand());
        }
      
    }

    public void chooseUpdateInterval() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "GUI update interval (ms)",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "1000");
        if(s!=null){
            int time = 1000;
            try { 
                time= Integer.parseInt(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            if(time>0) {
                this.setCanvasUpdate(time);
            }
            else {
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
        }
    }
        
    private JLabel getImage(String path,double scale) {
        JLabel label = null;
        Image image = null;
        try {
            URL url = new URL(path);
            image = ImageIO.read(url);
        } catch (IOException e) {
        	e.printStackTrace();
                System.out.println("Picture upload from " + path + " failed");
        }
        ImageIcon imageIcon = new ImageIcon(image);
        double width  = imageIcon.getIconWidth()*scale;
        double height = imageIcon.getIconHeight()*scale;
        imageIcon = new ImageIcon(image.getScaledInstance((int) width,(int) height, Image.SCALE_SMOOTH));
        label = new JLabel(imageIcon);
        return label;
    }
    
    public JPanel  getPanel(){
        return mainPanel;
    }
    
    public void setTriggerPhaseConstants(int run) {
		IndexedTable   jitter = this.monitors[0].engine.getConstantsManager().getConstants(run, "/calibration/ec/time_jitter");        
        PERIOD = jitter.getDoubleValue("period",0,0,0);
        PHASE  = jitter.getIntValue("phase",0,0,0); 
        CYCLES = jitter.getIntValue("cycles",0,0,0);
    }
    
    public long getTriggerWord(DataEvent event) {    	
 	    DataBank bank = event.getBank("RUN::config");	        
        return bank.getLong("trigger", 0);
    }
  
    public int getTriggerPhase(DataEvent event) {    	
 	    DataBank bank = event.getBank("RUN::config");	        
        long timestamp = bank.getLong("timestamp",0);    
        if (CYCLES==0) return 0;
        return (int) (PERIOD*((timestamp+PHASE)%CYCLES)); // TI derived phase correction due to TDC and FADC clock differences 
    }
    
    private int getRunNumber(DataEvent event) {
        int rNum = this.runNumber;
        DataBank bank = event.getBank("RUN::config");
        if(bank!=null) {
            rNum      = bank.getInt("run", 0);
        }
        return rNum;
    }
/*    
    @Override
    public void dataEventAction(DataEvent event) {
    	
        HipoDataEvent hipo = null;
        
	    if(event!=null ){

            if(event instanceof EvioDataEvent){
             	hipo = (HipoDataEvent) clasDecoder.getDataEvent(event);
                DataBank   header = clasDecoder.createHeaderBank(hipo, this.ccdbRunNumber, 0, (float) 0, (float) 0);
                hipo.appendBanks(header);
            } 
            else {
                hipo = (HipoDataEvent) event;    
            }
            
            if(this.runNumber != this.getRunNumber(hipo)) {
                this.runNumber = this.getRunNumber(hipo);
                System.out.println("Setting run number to: " +this.runNumber);
                resetEventListener();
            }
            
            for(int k=0; k<this.monitors.length; k++) {
                this.monitors[k].setTriggerPhase(getTriggerPhase(hipo));
                this.monitors[k].setTriggerWord(getTriggerWord(hipo));   
                this.monitors[k].dataEventAction(hipo);
            }      
	    }
    }
*/    
    @Override
    public void dataEventAction(DataEvent event) {
    	
       // EvioDataEvent decodedEvent = deco.DecodeEvent(event, decoder, table);
        //decodedEvent.show();
        		
        HipoDataEvent hipo = null;
        
        if(event!=null ){
            //event.show();

            if(event instanceof EvioDataEvent){
             	hipo = (HipoDataEvent) clasDecoder.getDataEvent(event);
                DataBank   header = clasDecoder.createHeaderBank(hipo, 0, 0, (float) 0, (float) 0);
                hipo.appendBanks(header);
            } 
            else {
                hipo = (HipoDataEvent) event;    
            }
            
            int rNum = this.runNumber;
            int eNum = this.eventNumber;
            if(event.hasBank("RUN::config")) {
                DataBank bank = event.getBank("RUN::config");
                 rNum      = bank.getInt("run", 0);
                 eNum      = bank.getInt("event", 0);
            }
            
            if(rNum!=0 && this.runNumber != rNum) {
                this.runNumber = rNum;
                for(int k=0; k<this.monitors.length; k++) {
                    this.monitors[k].setRunNumber(this.runNumber);
                }
                for(int k=0; k<this.monitors.length; k++) {
                    this.monitors[k].createHistos(this.runNumber);
                    this.monitors[k].initGStyle();
                    this.monitors[k].plotHistos(this.runNumber);
                    this.monitors[k].arBtn.setSelected(true);         
                    if (this.monitors[k].sectorButtons) {this.monitors[k].bS2.doClick();}
                }
            } 
           
            this.eventNumber = eNum;
            
            for(int k=0; k<this.monitors.length; k++) {
                this.monitors[k].setEventNumber(this.eventNumber);
            }
            
            setTriggerPhaseConstants(this.runNumber);
            
            for(int k=0; k<this.monitors.length; k++) {
                this.monitors[k].setTriggerPhase(getTriggerPhase(hipo));
                this.monitors[k].setTriggerWord(getTriggerWord(hipo));        	    
                this.monitors[k].dataEventAction(hipo);
            }      
	}
    }
    
    private void readFiles() {
        EvioSource     evioReader = new EvioSource();
        HipoDataSource hipoReader = new HipoDataSource();
        JFileChooser fc = new JFileChooser();
        fc.setDialogTitle("Choose input files directory...");
        fc.setMultiSelectionEnabled(true);
        fc.setAcceptAllFileFilterUsed(false);
        File workingDirectory = new File(this.workDir);
        fc.setCurrentDirectory(workingDirectory);
        int returnValue = fc.showOpenDialog(null);
        if (returnValue == JFileChooser.APPROVE_OPTION) {
            int nf = 0;
            for (File fd : fc.getSelectedFiles()) {
                if (fd.isFile()) {
                    if (fd.getName().contains(".evio") || fd.getName().contains(".hipo")) {
                        Integer current = 0;
                        Integer nevents = 0;
                        DataEvent event = null;
                        if(fd.getName().contains(".hipo")) {
                            hipoReader.open(fd);
                            current = hipoReader.getCurrentIndex();
                            nevents = hipoReader.getSize();                            
                        }
                        else if(fd.getName().contains(".evio")) {
                            evioReader.open(fd);
                            current = evioReader.getCurrentIndex();
                            nevents = evioReader.getSize();
                        }

                        System.out.println("\nFILE: " + nf + " " + fd.getName() + " N.EVENTS: " + nevents.toString() + "  CURRENT : " + current.toString());                        

                        for (int k = 0; k < nevents; k++) {
                            if(fd.getName().contains(".hipo")) {
                                if (hipoReader.hasEvent()) {
                                    event = hipoReader.getNextEvent();                          
                                }
                            }
                            else if(fd.getName().contains(".evio")) {
                                if (evioReader.hasEvent()) {
                                    event = evioReader.getNextEvent();
                                }
                            }
                            if(event != null) {
                                this.dataEventAction(event);
                                if(k % 10000 == 0) System.out.println("Read " + k + " events");
                            }
                        }
                        for(int k=0; k<this.monitors.length; k++) {
                            this.monitors[k].analyze();
                            this.monitors[k].fillSummary();
                            this.monitors[k].initGStyle();
                            this.monitors[k].plotHistos(this.getRunNumber(event));
                        }
                        nf++;
                    }
                }
            }
//            this.updateTable();
            System.out.println("Task completed");
        }
    }    
       
    @Override
    public void resetEventListener() {
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].resetEventListener();
            this.monitors[k].timerUpdate();
        }      
    }
   
    public void printHistosToFile() {
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String data = outPath + "/output" + "/clas12ana_" + this.runNumber + "_" + df.format(new Date());        
        File theDir = new File(data);
        // if the directory does not exist, create it
        if (!theDir.exists()) {
            boolean result = false;
            try{
                theDir.mkdir();
                result = true;
            } 
            catch(SecurityException se){
                //handle it
            }        
            if(result) {    
            System.out.println("Created directory: " + data);
            }
        }
        
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].printCanvas(data);
        }
        
        System.out.println("Histogram pngs succesfully saved in: " + data);
    }
    
    public int getFileRunNumber(String file) {
    	   String[] tokens = file.split("_");
    	   return Integer.parseInt(tokens[2]);
    }
    
    public void loadHistosFromFile(String fileName) {

        System.out.println("Opening file: " + fileName);
        runNumber = getFileRunNumber(fileName);    
        TDirectory dir = new TDirectory();
        dir.readFile(fileName);
        System.out.println(dir.getDirectoryList());
        dir.cd();
        dir.pwd();
        
        for(int k=0; k<this.monitors.length; k++) {        	   
        	this.monitors[k].setRunNumber(runNumber);
            this.monitors[k].createHistos(this.monitors[k].getRunNumber());           
            this.monitors[k].initGStyle();
            this.monitors[k].plotHistos(this.monitors[k].getRunNumber());
            this.monitors[k].arBtn.setSelected(true);         
            if (this.monitors[k].sectorButtons) {this.monitors[k].bS2.doClick();}
            this.monitors[k].readDataGroup(dir);
        }
    }
    
    public void saveHistosToFile(String fileName) {
        TDirectory dir = new TDirectory();
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].writeDataGroup(dir);
        }
        System.out.println("Saving histograms to file " + fileName);
        dir.writeFile(fileName);
    }
        
    public void setCanvasUpdate(int time) {
        System.out.println("Setting " + time + " ms update interval");
        this.canvasUpdateTime = time;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setCanvasUpdate(time);
        }
    }

    public void stateChanged(ChangeEvent e) {
        this.timerUpdate();
    }
    
    @Override
    public void timerUpdate() {
        if(this.runNumber==0) return;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].timerUpdate();
        }
   }

    public static void main(String[] args){
        JFrame frame = new JFrame("CLAS12Ana");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        EventViewer viewer = new EventViewer();
        frame.add(viewer.mainPanel);
        frame.setJMenuBar(viewer.menuBar);
        frame.setSize(1400, 800);
        frame.setVisible(true);
    }
    
    private void setRunNumber(String actionCommand) {
    
        System.out.println("Set run number for CCDB access");
        String  RUN_number = (String) JOptionPane.showInputDialog(null, "Set run number to ", " ", JOptionPane.PLAIN_MESSAGE, null, null, "2284");
        
        if (RUN_number != null) { 
            int cur_runNumber= this.runNumber;
            try {
                cur_runNumber = Integer.parseInt(RUN_number);
            } 
            catch (
                NumberFormatException f) {JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            if (cur_runNumber > 0){ 
                this.ccdbRunNumber = cur_runNumber;
                clasDecoder.setRunNumber(cur_runNumber,true);
            } 
            else {JOptionPane.showMessageDialog(null, "Value must be a positive integer!");}   
        }
        
    }

    private void resetHistograms(String actionCommand) {
        
        if (actionCommand=="Reset ECAL histograms"){
            System.out.println("Reset ECAL histograms");
         	int resetOption = JOptionPane.showConfirmDialog(null, "Do you want to automaticaly reset ECAL plots ?", " ", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
                if (resetOption == JOptionPane.YES_OPTION) {
                    String  resetTiming = (String) JOptionPane.showInputDialog(null, "Update every (number of events)", " ", JOptionPane.PLAIN_MESSAGE, null, null, "10000");
                    if (resetTiming != null) {    
                        int time = this.monitors[5].eventResetTime_default[5];
                        try {time = Integer.parseInt(resetTiming);} 
                        catch (NumberFormatException f) {JOptionPane.showMessageDialog(null, "Value must be a positive integer!");}
                        if (time > 0) {this.monitors[5].eventResetTime_current[5] = time;} 
                        else {JOptionPane.showMessageDialog(null, "Value must be a positive integer!");}   
                    }
                }else if (resetOption == JOptionPane.NO_OPTION){
			     this.monitors[5].eventResetTime_current[5] = 0;
                }	
        } 
        
    }

	@Override
	public void processShape(DetectorShape2D arg0) {
		// TODO Auto-generated method stub
		
	}
}