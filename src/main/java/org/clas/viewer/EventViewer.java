package org.clas.viewer;

import java.awt.BorderLayout;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.swing.ButtonGroup;
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
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileSystemView;

import org.clas.analysis.ECa;
import org.clas.analysis.ECelas;
import org.clas.analysis.ECmc;
import org.clas.analysis.ECmip;
import org.clas.analysis.ECperf;
import org.clas.analysis.ECpi0;
import org.clas.analysis.ECsf;
import org.clas.analysis.ECt;
import org.clas.analysis.ECcalib;

import org.jlab.detector.decode.CLASDecoder4;
import org.jlab.detector.decode.CodaEventDecoder;
import org.jlab.detector.decode.DetectorEventDecoder;
import org.jlab.detector.view.DetectorListener;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.TDirectory;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.evio.EvioDataEvent;
import org.jlab.io.evio.EvioSource;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo3.Hipo3DataSource;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.io.task.DataSourceProcessorPane;
import org.jlab.io.task.IDataEventListener;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.system.ClasUtilsFile;
        
/*
 * @author lcsmith
 * Adapted from KPP-plots framework developed by R. DeVita, V. Ziegler, G. Gavalian
 */

public class EventViewer implements IDataEventListener, DetectorListener, ActionListener, ItemListener, ChangeListener {
    
    JTabbedPane                  tabbedpane = null;
    JPanel                        mainPanel = null;
    JMenuBar                        menuBar = null;    
    DataSourceProcessorPane   processorPane = null;
    
    JCheckBoxMenuItem co0,co1,co2,co3,co4,co4b,co5,co6,co7 = null;   
    JCheckBoxMenuItem       cf,cf0,cf1,cf2,cf3,cf4,cf5 = null;   
    JCheckBoxMenuItem                              ctr = null;  
    JRadioButtonMenuItem               ct0,ct1,ct2,ct3 = null;  
    JRadioButtonMenuItem ctr0,ctr1,ctr2,ctr3,ctr4,ctr5 = null;  
    
    CodaEventDecoder               decoder = new CodaEventDecoder();
    CLASDecoder4               clasDecoder = new CLASDecoder4();
    DetectorEventDecoder   detectorDecoder = new DetectorEventDecoder();
    private SchemaFactory    schemaFactory = new SchemaFactory();
       
    private int   canvasUpdateTime = 2000;
    private int           TVOffset = 0;
    private float         logParam = 0;
    private int analysisUpdateEvnt = 100;
    private int          runNumber = 0;
    private int        eventNumber = 0;
    private int      ccdbRunNumber = 0;
   
    double PERIOD = 0;
    int     PHASE = 0;
    int    CYCLES = 0;
    
    public String outPath = "/Users/cole/CLAS12ANA/";
    public String workDir = outPath;
    
    public Boolean   clearHist = true;
    public Boolean    autoSave = false;
    public Boolean   dropBanks = false;
    public Boolean dropSummary = false;
    public Boolean  dumpGraphs = false;
    public Boolean   dumpFiles = false;
    public Boolean defaultGain = false;
    public Boolean    fiduCuts = false;
    public Boolean   dropEsect = false;
    public Boolean  cfitEnable = false;
    public Boolean  sfitEnable = false;
    public Boolean  dfitEnable = false;
    public Boolean gdfitEnable = false;
    public Boolean  yLogEnable = false;
    public Boolean  zLogEnable = true;
    public Boolean  fitVerbose = false;
    public String       TLname = "UVW";
    public Boolean      TLflag = false;
    public Boolean       clear = true; 
    public Integer       TRpid = 11;
    DetectorMonitor[] monitors = null;
    public JFileChooser     fc = null; 
    
    List<Integer>     runList  = new ArrayList<Integer>();
    
    Map<String,DetectorMonitor> Monitors = new LinkedHashMap<String,DetectorMonitor>();
    
    int    selectedTabIndex = 0;
    String selectedTabName  = " ";
        
    public static void main(String[] args){

        JFrame frame = new JFrame("CLAS12Ana");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        EventViewer viewer = new EventViewer(args);
        frame.add(viewer.mainPanel);
        frame.setJMenuBar(viewer.menuBar);
        frame.setSize(1400, 800);
        frame.setVisible(true);
    }
    
    public EventViewer(String[] args) {  
    	createMonitors(args);
    	createMenuBar();
    	createPanels();
   	    String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        schemaFactory.initFromDirectory(dir);
    }
    
    public void createMonitors(String[] args) {
    	
        workDir = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/CLAS12ANA/";
     	monitors = new DetectorMonitor[args.length==0 ? 1:args.length];
        int n = 0;
    	if (args.length != 0) {
        	for(String s : args) { 
        	   switch (s) {
      	         case     "ECa": monitors[n++]=new ECa(s);    break;  
    	         case    "ECsf": monitors[n++]=new ECsf(s);   break; 
        	     case     "ECt": monitors[n++]=new ECt(s);    break;
        	     case   "ECmip": monitors[n++]=new ECmip(s);  break; 
        	     case "ECcalib": monitors[n++]=new ECcalib(s);break; 
        	     case   "ECpi0": monitors[n++]=new ECpi0(s);  break;
        	     case  "ECperf": monitors[n++]=new ECperf(s); break;
        	     case    "ECmc": monitors[n++]=new ECmc(s);   break;
        	     case  "ECelas": monitors[n++]=new ECelas(s); 
        	   }
        	}
    	} else {
//   		monitors[n] = new ECperf("ECperf"); 
    		monitors[n] = new ECmc("ECmc");
//    		monitors[n] = new ECt("ECt"); 
//      		monitors[n] = new ECsf("ECsf"); 
//    		monitors[n] = new ECcalib("ECcalib"); 
//    		monitors[n] = new ECmip("ECmip"); 
//    		monitors[n] = new ECpi0("ECpi0");

        }
    }
    
    public void createMenuBar() {
        		
        menuBar = new JMenuBar();
        
        JMenu menu;
        JMenuItem menuItem;
               
        menu     = new JMenu("File");
        menuItem = new JMenuItem("Load Run");                menuItem.addActionListener(this); menu.add(menuItem);
        menuItem = new JMenuItem("Analyze Runs");            menuItem.addActionListener(this); menu.add(menuItem);
        menuItem = new JMenuItem("Analyze Histos");          menuItem.addActionListener(this); menu.add(menuItem); menu.addSeparator();
        menuItem = new JMenuItem("Open histograms file");    menuItem.addActionListener(this); menu.add(menuItem);
        menuItem = new JMenuItem("Save histograms to file"); menuItem.addActionListener(this); menu.add(menuItem);
        menuItem = new JMenuItem("Save timelines to file");  menuItem.addActionListener(this); menu.add(menuItem);
        menuItem = new JMenuItem("Print histograms as png"); menuItem.addActionListener(this); menu.add(menuItem);
        menuBar.add(menu);

        menu     = new JMenu("Options");  
        co0  = new JCheckBoxMenuItem("ClearHist");     co0.addItemListener(this);       menu.add(co0);  co0.doClick();
        co1  = new JCheckBoxMenuItem("AutoSave");      co1.addItemListener(this);       menu.add(co1);
        co2  = new JCheckBoxMenuItem("DropBanks");     co2.addItemListener(this);       menu.add(co2);
        co3  = new JCheckBoxMenuItem("DropSummary");   co3.addItemListener(this);       menu.add(co3);
        co4  = new JCheckBoxMenuItem("DumpGraphs");    co4.addItemListener(this);       menu.add(co4);
        co4b = new JCheckBoxMenuItem("DumpFiles");    co4b.addItemListener(this);       menu.add(co4b);
        co5  = new JCheckBoxMenuItem("DefaultGains");  co5.addItemListener(this);       menu.add(co5);
        co6  = new JCheckBoxMenuItem("FiduCuts");      co6.addItemListener(this);       menu.add(co6);
        co7  = new JCheckBoxMenuItem("DropEsect");     co7.addItemListener(this);       menu.add(co7);
        menuBar.add(menu);
        
        menu     = new JMenu("Fitting");
        cf  = new JCheckBoxMenuItem("Verbose");        cf.addItemListener(this);       menu.add(cf);
        cf0 = new JCheckBoxMenuItem("Calibration");   cf0.addItemListener(this);       menu.add(cf0);  
        cf1 = new JCheckBoxMenuItem("Residual");      cf1.addItemListener(this);       menu.add(cf1);  
        cf2 = new JCheckBoxMenuItem("TMF");           cf2.addItemListener(this);       menu.add(cf2);  
        cf3 = new JCheckBoxMenuItem("GTMF");          cf3.addItemListener(this);       menu.add(cf3);
        menuBar.add(menu);
              
        menu     = new JMenu("Settings");       
        menuItem = new JMenuItem("Set logParam");                menuItem.addActionListener(this); menu.add(menuItem);   
        menuItem = new JMenuItem("Set GUI update interval");     menuItem.addActionListener(this); menu.add(menuItem);
        cf4      = new JCheckBoxMenuItem("log Y");                    cf4.addItemListener(this);   menu.add(cf4);  
        cf5      = new JCheckBoxMenuItem("log Z");                    cf5.addItemListener(this);   menu.add(cf5);  cf5.doClick();
        menuItem = new JMenuItem("Set run number");              menuItem.addActionListener(this); menu.add(menuItem);
        menuBar.add(menu);
        
        menu     = new JMenu("Reset");
        menuItem = new JMenuItem("Default for all");         menuItem.addActionListener(this); menu.add(menuItem);
        menuItem = new JMenuItem("Disable histogram reset"); menuItem.addActionListener(this); menu.add(menuItem);
        menuBar.add(menu);
        
        menu     = new JMenu("Timelines");
        ButtonGroup group = new ButtonGroup();
        ct0 = new JRadioButtonMenuItem("UVW");        ct0.addItemListener(this);       group.add(ct0); menu.add(ct0); 
        ct1 = new JRadioButtonMenuItem("FADC Slot");  ct1.addItemListener(this);       group.add(ct1); menu.add(ct1);
        ct2 = new JRadioButtonMenuItem("HV Slot");    ct2.addItemListener(this);       group.add(ct2); menu.add(ct2);
        ct3 = new JRadioButtonMenuItem("Sectors");    ct3.addItemListener(this);       group.add(ct3); menu.add(ct3);
        menuBar.add(menu);

        String TriggerDef[] = { "Electron OR",
		        "e Sector 1","e Sector 2","e Sector 3","e Sector 4","e Sector 5","e Sector 6",
		        "Muons S1+ S4-","Muons S2+ S5-","Muons S3+ S6-",
		        "Muons S4+ S1-","Muons S5+ S2-","Muons S6+ S3-","","","","","","","","","","","","","","","","","","",
		        "1K Pulser"};   
        
        menu = new JMenu("TriggerBits");
        
        for (int i=0; i<32; i++) {
        	
            ctr = new JCheckBoxMenuItem(TriggerDef[i]);  
            final Integer bit = new Integer(i);
            
            ctr.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {               	
                    if(e.getStateChange() == ItemEvent.SELECTED) {
                        for(int k=0; k<monitors.length; k++) monitors[k].setTriggerMask(bit);
                    } else {
                        for(int k=0; k<monitors.length; k++) monitors[k].clearTriggerMask(bit);
                    };
                }
            });         
            menu.add(ctr);         	        	 
        }
        
        menuBar.add(menu);
        
        menu     = new JMenu("Trigger PID");
        group = new ButtonGroup();
        ctr0 = new JRadioButtonMenuItem("Electron"); ctr0.addItemListener(this);       group.add(ctr0); menu.add(ctr0); ctr0.doClick();
        ctr1 = new JRadioButtonMenuItem("Pion");     ctr1.addItemListener(this);       group.add(ctr1); menu.add(ctr1);
        ctr2 = new JRadioButtonMenuItem("Photon");   ctr2.addItemListener(this);       group.add(ctr2); menu.add(ctr2);
        ctr3 = new JRadioButtonMenuItem("Neutron");  ctr3.addItemListener(this);       group.add(ctr3); menu.add(ctr3);
        ctr4 = new JRadioButtonMenuItem("Muon");     ctr4.addItemListener(this);       group.add(ctr4); menu.add(ctr4);
        ctr5 = new JRadioButtonMenuItem("PC Muon");  ctr5.addItemListener(this);       group.add(ctr5); menu.add(ctr5);
        
        menuBar.add(menu);
    
    }
    
    public void createPanels() {

        mainPanel = new JPanel();	
        mainPanel.setLayout(new BorderLayout());
        
      	tabbedpane 	= new JTabbedPane();

        processorPane = new DataSourceProcessorPane();
        processorPane.setUpdateRate(analysisUpdateEvnt);

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
/*        
        this.tabbedpane.addChangeListener(new ChangeListener() {   
        	public void stateChanged(ChangeEvent e) {
        	System.out.println("Change Event");
        	if (e.getSource() instanceof JTabbedPane) {
        		JTabbedPane pane = (JTabbedPane) e.getSource();
        		selectedTabIndex = pane.getSelectedIndex();
        		selectedTabName  = (String) pane.getTitleAt(selectedTabIndex);
        		System.out.println(selectedTabIndex + " " + selectedTabName);
        	}
        	}
        });
*/        
        this.tabbedpane.addChangeListener(this);	               
        this.processorPane.addEventListener(this);        
        this.setCanvasUpdate(canvasUpdateTime);
        
    }
    
    @Override
    public void stateChanged(ChangeEvent e) {
        this.timerUpdate();
    }
    
	@Override
	public void itemStateChanged(ItemEvent e) {
		Object source = e.getItemSelectable();	
		
		if (source==co0)    clearHist = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==co1)     autoSave = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==co2)    dropBanks = (e.getStateChange() == ItemEvent.SELECTED)?true:false;	
		if (source==co3)  dropSummary = (e.getStateChange() == ItemEvent.SELECTED)?true:false;	
		if (source==co4)   dumpGraphs = (e.getStateChange() == ItemEvent.SELECTED)?true:false;	
		if (source==co4b)   dumpFiles = (e.getStateChange() == ItemEvent.SELECTED)?true:false;	
		if (source==co5)  defaultGain = (e.getStateChange() == ItemEvent.SELECTED)?true:false;	
		if (source==co6)     fiduCuts = (e.getStateChange() == ItemEvent.SELECTED)?true:false;	
		if (source==co7)    dropEsect = (e.getStateChange() == ItemEvent.SELECTED)?true:false;	
		if (source==cf )   fitVerbose = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==cf0)   cfitEnable = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==cf1)   sfitEnable = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==cf2)   dfitEnable = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==cf3)  gdfitEnable = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==cf4)   yLogEnable = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==cf5)   zLogEnable = (e.getStateChange() == ItemEvent.SELECTED)?true:false;
		if (source==ct0)       TLname = ct0.getText();
		if (source==ct1)       TLname = ct1.getText();
		if (source==ct2)       TLname = ct2.getText();
		if (source==ct3)       TLflag = (e.getStateChange() == ItemEvent.SELECTED)?true:false; 
		if (source==ctr0)       TRpid = (e.getStateChange() == ItemEvent.SELECTED)?   11:11; 
		if (source==ctr1)       TRpid = (e.getStateChange() == ItemEvent.SELECTED)?  211:11; 
		if (source==ctr2)       TRpid = (e.getStateChange() == ItemEvent.SELECTED)?   22:11; 
		if (source==ctr3)       TRpid = (e.getStateChange() == ItemEvent.SELECTED)? 2112:11; 
		if (source==ctr4)       TRpid = (e.getStateChange() == ItemEvent.SELECTED)?    0:11; 
		if (source==ctr5)       TRpid = (e.getStateChange() == ItemEvent.SELECTED)?   -1:11; 
		
		for(int k=0; k<this.monitors.length; k++) {this.monitors[k].dropBanks   = dropBanks; 
		                                           this.monitors[k].dropSummary = dropSummary; 
		                                           this.monitors[k].dumpGraphs  = dumpGraphs;
		                                           this.monitors[k].dumpFiles   = dumpFiles;
		                                           this.monitors[k].defaultGain = defaultGain;
		                                           this.monitors[k].fiduCuts    = fiduCuts;
		                                           this.monitors[k].dropEsect   = dropEsect;
		                                           this.monitors[k].autoSave    = autoSave;
                                                   this.monitors[k].cfitEnable  = cfitEnable;
                                                   this.monitors[k].sfitEnable  = sfitEnable;
                                                   this.monitors[k].dfitEnable  = dfitEnable;
                                                   this.monitors[k].gdfitEnable = gdfitEnable;
                                                   this.monitors[k].setLogY(yLogEnable);                                                  
                                                   this.monitors[k].setLogZ(zLogEnable);
                                                   this.monitors[k].fitVerbose  = fitVerbose;
                                                   this.monitors[k].TRpid       = TRpid;
                                                   this.monitors[k].initTimeLine(TLname);
                                                   this.monitors[k].setTLflag(TLflag);}
    }  
	
    public void actionPerformed(ActionEvent e) {
        System.out.println(e.getActionCommand());
        switch (e.getActionCommand()) {
          case("Load Run"):                    this.loadHistoFromRunIndex(); break;
//          case("Load Summary"):                this.readHistosFromSummary(); break;
          case("Analyze Runs"):                this.readFiles(); break;
          case("Analyze Histos"):              this.readHistos(); break;
          case("Set GUI update interval"):     this.chooseUpdateInterval(); break;
          case("Set logParam"):                this.chooseLogParam(); break;
          case("Set run number"):              this.setRunNumber(e.getActionCommand()); break;
          case("Open histograms file"):        this.histoChooser(); break;
          case("Save histograms to file"):     this.histoSaver(); break;
          case("Save timelines to file"):      this.timelineSaver(); break;
          case("Print histograms as png"):     this.printHistosToFile(); break;
          case("Default for all"):             for (int k=0;k<monitors.length;k++) this.monitors[k].eventResetTime_current[k] = this.monitors[k].eventResetTime_default[k]; break;       
          case("Disable histogram reset"):     for (int k=0;k<monitors.length;k++) this.monitors[k].eventResetTime_current[k] = 0; break;
          case("Reset"):                       for (int k=0;k<monitors.length;k++) this.monitors[k].eventResetTime_current[k] = 0;
        }
        if (e.getActionCommand().substring(0, 5).equals("Reset")) resetHistograms(e.getActionCommand());
    }
    
    public void histoChooser() {
        JFileChooser fc = new JFileChooser(new File(workDir));
        fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) this.loadHistosFromFile(fc.getSelectedFile().getAbsolutePath());   	
    }
    
    public void histoSaver() {
//      DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String fileName = "CLAS12Ana_run_" + this.runNumber + "_" + monitors[0].getDetectorName() + ".hipo";
        if(autoSave) {this.saveHistosToFile(workDir+fileName);return;}
        JFileChooser fc = new JFileChooser(new File(workDir));
        fc.setSelectedFile(new File(fileName));
        if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION) this.saveHistosToFile(fc.getSelectedFile().getAbsolutePath());   	
    }
    
    public void timelineSaver() {
        for(int k=0; k<this.monitors.length; k++) this.monitors[k].saveTimelines();	
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
            int val = 1000;
            try { 
                val= Integer.parseInt(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            if(val>0) {
                this.setCanvasUpdate(val);
            }
            else {
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
        }
    }
    
    public void chooseLogParam() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "ECEngine logParam",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "0");
        if(s!=null){
            float val = 0;
            try { 
                val= Float.parseFloat(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            this.setLogParam(val);
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
		IndexedTable   jitter = this.monitors[0].cm.getConstants(run, "/calibration/ec/time_jitter");        
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
        DataBank bank = event.getBank("RUN::config"); 
        if(this.ccdbRunNumber >0) return this.ccdbRunNumber;
        return (bank!=null) ? bank.getInt("run",0):this.runNumber;
    }
    
    private int getEventNumber(DataEvent event) {
        DataBank bank = event.getBank("RUN::config");
        return (bank!=null) ? bank.getInt("event", 0): this.eventNumber;
    }
    
    private boolean isGoodRun(int run) {
    	if( runList.isEmpty()) {runList.add(this.runNumber); return true;}
    	if(!runList.isEmpty()) {
    		if(runList.contains(run)) return false;
    	}
    	runList.add(run);
    	return true;
    }
    
    @Override
    public void dataEventAction(DataEvent event) {

	    if(event!=null ) processEvent(filterEvent(decodeEvent(event)));

    }
    
    private boolean processEvent(DataEvent event) {
    	
    	if(event==null) return false; 
      
        this.eventNumber = getEventNumber(event);
        
        setTriggerPhaseConstants(this.runNumber);

        for(int k=0; k<this.monitors.length; k++) {
        	this.monitors[k].setEventNumber(this.eventNumber);
        	this.monitors[k].setTriggerPhase(getTriggerPhase(event));
            this.monitors[k].setTriggerWord(getTriggerWord(event));        	    
            this.monitors[k].dataEventAction(event);
        }  
        
        return true;        
    }
    
    private DataEvent filterEvent(DataEvent event) {
    	
        int rNum = 0; 
        
        rNum = getRunNumber(event);
       
        if(rNum!=0 && this.runNumber!=rNum && clear) {  
        	System.out.println("EventViewer: Processing Run "+rNum);
        	this.runNumber = rNum;
            if(!clearHist) clear=false;
        	return isGoodRun(rNum)? initRun(rNum,event):null;
        }
        
        return event;
    }
    
    private DataEvent decodeEvent(DataEvent event) {
    	
    	DataEvent hipo = null;
        if(event instanceof EvioDataEvent){
         	Event    dump = clasDecoder.getDataEvent(event);    
            Bank   header = clasDecoder.createHeaderBank(this.ccdbRunNumber, getEventNumber(event), (float) 0, (float) 0);
            Bank  trigger = clasDecoder.createTriggerBank();
            if(header!=null)  dump.write(header);
            if(trigger!=null) dump.write(trigger);
            hipo = new HipoDataEvent(dump,schemaFactory); 
        }   
        else {            	
        	hipo = event; 
        } 
        
        return hipo;
    }
    
    private DataEvent initRun(int runno, DataEvent event) {
    	
        for(int k=0; k<this.monitors.length; k++) {
        	if(autoSave && this.runNumber!=0) this.monitors[k].saveHistosToFile();
            this.runNumber = runno; 
            this.monitors[k].isHipo3Event = processorPane.isHipo3Event;
        	this.monitors[k].setRunNumber(this.runNumber); 
           	this.monitors[k].localclear();
           	this.monitors[k].initCCDB(this.runNumber);
           	this.monitors[k].initEBCCDB(this.runNumber);
        	this.monitors[k].createHistos(this.runNumber);
            this.monitors[k].initGStyle();
            this.monitors[k].plotHistos(this.runNumber);
            this.monitors[k].arBtn.setSelected(true);         
            if(this.monitors[k].sectorButtons) this.monitors[k].bS2.doClick();
        } 
        
        return event;
    }
    
    private void readFiles() {
        HipoDataSource   hipoReader = new HipoDataSource();
        JFileChooser fc = new JFileChooser(new File(workDir));
        fc.setDialogTitle("Choose input files directory...");
        fc.setMultiSelectionEnabled(true);
        fc.setAcceptAllFileFilterUsed(true);
        if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
            int nf = 0;
            for (File fd : fc.getSelectedFiles()) {
                if (fd.isFile()) {
                    if (fd.getName().contains(".hipo")) {
                        Integer current = 0;
                        Integer nevents = 0;
                        DataEvent event = null;
                        hipoReader.open(fd);
                        current = hipoReader.getCurrentIndex();
                        nevents = hipoReader.getSize();
                        System.out.println("\nFILE: " + nf + " " + fd.getName() + " N.EVENTS: " + nevents.toString() + "  CURRENT : " + current.toString());                        
                        for (int k = 0; k < nevents; k++) {
                        	event = hipoReader.getNextEvent();                
                            if(event != null) {
                                this.dataEventAction(event);
                                if(k % 10000 == 0) System.out.println("Read " + k + " events");
                            }
                        }
                        hipoReader.close();
                        nf++;
                    }
                }
            }
            System.out.println("Task completed");
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
    
    public String[] getFileTokens(String file) {
   	    return (file.split(Pattern.quote("."))[0].split("_"));
    }
    
    public int getFileRunNumber(String file) {
    	return Integer.parseInt(getFileTokens(file)[2]);
    }
    
    public boolean isCalibrationFile(String file) {
    	return getFileTokens(file).length>5;
    }
    
    public String getFileCalibrationTag(String file) {
    	return getFileTokens(file)[4];
    }
    
    public int getFileCalibrationIndex(String file) {
        return Integer.parseInt(getFileTokens(file)[5]);
    }
  
    @Override
    public void resetEventListener() {
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].resetEventListener();
            this.monitors[k].timerUpdate();
        }      
    } 
    
    private List<File> selectHistos() {
        fc = new JFileChooser(new File(workDir+monitors[0].getDetectorName()));
        fc.setDialogTitle("Choose input histos directory...");
        fc.setMultiSelectionEnabled(true);
        fc.setAcceptAllFileFilterUsed(false);  
        List<File> files = new ArrayList<File>();
        if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
        	if(fc.getSelectedFile().getName().contains("summary")) {
        		String fname = fc.getSelectedFile().getAbsolutePath();
        		try 
        		{
        	    BufferedReader reader = new BufferedReader(new FileReader(fname));
        	    String line;
        	    while ((line = reader.readLine()) != null)
        	    {        	    	
        	      files.add(new File(line));
        	    }
        		reader.close();
        	}
        	  catch (Exception e)
        	  {
        	    System.err.format("Exception occurred trying to read '%s'.", fname);
        	    e.printStackTrace();        	   
        	  }
        	} else {
                for (File fd : fc.getSelectedFiles()) {
                    if (fd.isFile()) {
                        files.add(new File(fd.getAbsolutePath()));
                    }                
                } 
        	}
        } 
        return files;
    }
    
    private void readHistos() {
    	String fname = null;
        for(int k=0; k<this.monitors.length; k++) this.monitors[k].localclear();

        for (File fd : selectHistos()) {
            if (fd.isFile()) {
                fname=fd.getAbsolutePath();
                loadHistosFromFile(fname);
            }                
        }  
        if(isCalibrationFile(fname)) this.monitors[0].writeFile(getFileCalibrationTag(fname),1,7,0,3,0,3);
    }
             
    public void loadHistoFromRunIndex() {
    	File[] f = fc.getSelectedFiles();
    	monitors[0].dropSummary=false;
    	loadHistosFromFile(f[monitors[0].getRunIndex()].getAbsolutePath());
    }
    
    public void loadHistosFromFile(String fileName) {
        System.out.println("EventViwer.loadHistosFromFile("+fileName+")");
        
        runNumber = getFileRunNumber(fileName); 
        
        TDirectory dir = new TDirectory(); dir.readFile(fileName); dir.cd(); dir.pwd();
        
        for(int k=0; k<this.monitors.length; k++) {  
            if(isCalibrationFile(fileName)) this.monitors[k].detcal[getFileCalibrationIndex(fileName)]=runNumber;
         	this.monitors[k].setRunNumber(runNumber);
            this.monitors[k].createHistos(runNumber);  
            this.monitors[k].initCCDB(runNumber);
            this.monitors[k].initEBCCDB(runNumber);
            this.monitors[k].initGStyle();
            this.monitors[k].readDataGroup(runNumber,dir);
            if (this.monitors[k].sectorButtons) {this.monitors[k].bS2.doClick();}
            this.monitors[k].arBtn.setSelected(true);   
        }
        return;
    }
    
    public void saveHistosToFile(String fileName) {
        for(int k=0; k<this.monitors.length; k++) {
        	TDirectory dir = new TDirectory();
            this.monitors[k].writeDataGroup(dir);
            dir.writeFile(fileName);
            System.out.println("EventViewer.saveHistosToFile("+fileName+")");
        }
    } 
    
    public void setCanvasUpdate(int time) {
        System.out.println("EventViewer.setCanvasUpdate("+time+")");
        this.canvasUpdateTime = time;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setCanvasUpdate(time);
        }
    }
    
    public void setTVOffset(int time) {
        System.out.println("Setting Tvertex offset " + time + " ns");
        this.TVOffset = time;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setTVOffset(time);
        }
    }
    
    public void setLogParam(float val) {
        System.out.println("EventViwer.setLogParam("+val+")");
        this.logParam = val;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setLogParam(val);
        }
    }
    
    @Override
    public void timerUpdate() {
        if(this.runNumber==0) return;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].timerUpdate();
        }
   }
    
    private void setRunNumber(String actionCommand) {
    
        System.out.println("EventViewer.setRunNumber("+actionCommand+")");
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

    	System.out.println("EventViewer.resetHistograms("+actionCommand+")");
        if (actionCommand=="Reset ECAL histograms"){
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