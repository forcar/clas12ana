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
import java.util.TimerTask;
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
import org.clas.analysis.ECmc1;
import org.clas.analysis.ECmc2;
import org.clas.analysis.ECmcn;
import org.clas.analysis.ECmip;
import org.clas.analysis.ECmon;
import org.clas.analysis.ECperf;
import org.clas.analysis.ECpi0;
import org.clas.analysis.ECstatus;
import org.clas.analysis.ECsf;
import org.clas.analysis.ECt;
import org.clas.tools.DataSourceProcessorPane;
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
import org.jlab.io.base.DataEventType;
import org.jlab.io.evio.EvioDataEvent;
import org.jlab.io.evio.EvioSource;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo3.Hipo3DataSource;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.io.task.DataSourceProcessor;
//import org.jlab.io.task.DataSourceProcessorPane;
import org.jlab.io.task.IDataEventListener;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.system.ClasUtilsFile;
        
/*
 * @author lcsmith
 * Adapted from KPP-plots framework developed by R. DeVita, V. Ziegler, G. Gavalian
 */

public class EventViewer implements IDataEventListener, DetectorListener, ActionListener, ItemListener, ChangeListener {
    
    JPanel                        mainPanel = null;
    JMenuBar                        menuBar = null;    
    JTabbedPane                  tabbedpane = null;    
    JCheckBoxMenuItem  co0,co1,co2,co3,co4,co4b,co5,co6,co7,co8;   
    JCheckBoxMenuItem   cf,cf0,cf1,cf2,cf3,cf4,cf5,cf6a,cf6b,cf6c,cf6d,cf6e,cf6f,cf6g,cf6h,cf6i,cf6j;
    JCheckBoxMenuItem   cf7,cf8,cf9,cf10a,cf10b,cf10c,cf11,cf12,cf13,cf14,cf15,cf16,cf17,cf18;   
    JCheckBoxMenuItem                                   ctr;    
    JRadioButtonMenuItem                    ct0,ct1,ct2,ct3;  
    JRadioButtonMenuItem ctr0,ctr1,ctr2,ctr3,ctr4,ctr5,ctr6,ctr7; 
    
    DataSourceProcessorPane  processorPane = null;
    CodaEventDecoder               decoder = new CodaEventDecoder();
    CLASDecoder4               clasDecoder = new CLASDecoder4();
    DetectorEventDecoder   detectorDecoder = new DetectorEventDecoder();
    private SchemaFactory    schemaFactory = new SchemaFactory();
       
    private int   canvasUpdateTime = 2000;
    private int           TVOffset = 0;
    private float         logParam = 0;

    private int analysisUpdateEvnt = 100;
    private int          runNumber = 1;
    private int        eventNumber = 0;
    private int      ccdbRunNumber = 0;
   
    double PERIOD = 0;
    int     PHASE = 0;
    int    CYCLES = 0;
    
    public String outPath = "/Users/cole/CLAS12ANA/";
    public String workDir = outPath;
    
    public Boolean      clearHist = false;
    public Boolean       autoSave = false;
    public Boolean      dropBanks = false;
    public Boolean    dropSummary = false;
    public Boolean     dumpGraphs = false;
    public Boolean      dumpFiles = false;
    public Boolean    defaultGain = false;
    public Boolean       fiduCuts = false;
    public Boolean      dropEsect = false;
    public Boolean      onlyEsect = false;
    public Boolean       FTOFveto = true;
    public Boolean     cfitEnable = false;
    public Boolean     sfitEnable = false;
    public Boolean     dfitEnable = false;
    public Boolean    gdfitEnable = false;
    public Boolean    trfitEnable = false;
    public Boolean     yLogEnable = false;
    public Boolean     zLogEnable = true;
    public Boolean    dbgECEngine = false;
    public Boolean    dbgAnalyzer = false;
    public Boolean   unsharedTime = true;
    public Boolean unsharedEnergy = true;
    public Boolean        normPix = false;
    public Boolean        normAtt = false;
    public Boolean         normCz = false;
    public Boolean         SFcorr = false;
    public Boolean          HiRes = false;
    public Boolean         TWcorr = true;
    public Boolean         DTcorr = true;
    public Boolean         FTpcal = true;
    public Boolean     fitVerbose = false;
    public String          TLname = "UVW";
    public Boolean         TLflag = false;
    public Boolean          clear = true; 
    public Integer          TRpid = 11;
    public Boolean      useATDATA = false;
    public Boolean    StatEachRun = false;
    public Boolean    useFADCTime = false;
    public Boolean usePass2Timing = true;
    public Boolean usePass2Energy = true;
    public Boolean  useCalibPass2 = true;
    public Boolean         useGPP = false;
    public Boolean   outputECHITS = false; 
    public Boolean         useASA = false;
    public Boolean          useCC = false;
    
    public JFileChooser      fc = null; 
    
    List<Integer>     runList  = new ArrayList<Integer>();
    
    DetectorMonitor[]           monitors = null;
    Map<String,DetectorMonitor> Monitors = new LinkedHashMap<String,DetectorMonitor>();
    
    int    selectedTabIndex = 0;
    String selectedTabName  = " ";
    
    private DataSourceProcessor  dataProcessor = null;    
    private java.util.Timer      processTimer  = null;
    private int                  eventDelay    = 0;
   
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
    	System.out.println("*** WELCOME TO CLAS12ANA ***\n");
    	dataProcessor = new DataSourceProcessor();    
    	createMonitors(args);
    	initMonitors();
    	createMenuBar();
    	createPanels();   	    
        schemaFactory.initFromDirectory(ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4"));
        System.out.println("EventViewer init complete \n");
    }
    
    public void createMonitors(String[] args) {
        workDir = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/CLAS12ANA/";
    	System.out.println("EventViewer.createMonitors: workDir="+workDir);
     	monitors = new DetectorMonitor[args.length==0 ? 1:args.length];
        int n = 0;
    	if (args.length != 0) {
        	for(String s : args) { 
        	   switch (s) {
      	         case      "ECa": monitors[n++]=new ECa(s);    break;  
    	         case     "ECsf": monitors[n++]=new ECsf(s);   break; 
        	     case      "ECt": monitors[n++]=new ECt(s);    break;
        	     case    "ECmip": monitors[n++]=new ECmip(s);  break; 
        	     case	 "ECmon": monitors[n++]=new ECmon(s);  break;
        	     case  "ECcalib": monitors[n++]=new ECcalib(s);break; 
        	     case    "ECpi0": monitors[n++]=new ECpi0(s);  break;
        	     case   "ECperf": monitors[n++]=new ECperf(s); break;
        	     case     "ECmc": monitors[n++]=new ECmc(s);   break;
        	     case    "ECmcn": monitors[n++]=new ECmcn(s);  break;
        	     case    "ECmc1": monitors[n++]=new ECmc1(s);  break;
        	     case    "ECmc2": monitors[n++]=new ECmc2(s);  break;
        	     case "ECstatus": monitors[n++]=new ECstatus(s,"ECAL"); break;
        	     case   "ECelas": monitors[n++]=new ECelas(s); 
        	   }
        	}
    	} else {
// 		    monitors[n] = new ECperf("ECperf"); 
//    		monitors[n] = new ECelas("ECelas");
//    		monitors[n] = new ECa("ECa");
//    		monitors[n] = new ECmc2("ECmc2");
//    		monitors[n] = new ECstatus("ECstatus","ECAL");
//    		monitors[n] = new ECmc("ECmc");
//    		monitors[n] = new ECmc1("ECmc1");
//    		monitors[n] = new ECmc2("ECmc2");
//    		monitors[n] = new ECmcn("ECmcn");
//     	    monitors[n] = new ECt("ECt"); 
//          monitors[n] = new ECperf("ECperf");
//      		monitors[n] = new ECsf("ECsf"); 
//    		monitors[n] = new ECcalib("ECcalib"); 
    		monitors[n] = new ECmon("ECmon"); 
//    		monitors[n] = new ECmip("ECmip"); 
//    		monitors[n] = new ECpi0("ECpi0"); 

        }
    }
    
    public void createMenuBar() {
    	System.out.println("EventViewer.createMenuBar");
        		
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
        co0  = new JCheckBoxMenuItem("ClearHist");     co0.addItemListener(this);       menu.add(co0);  
        co1  = new JCheckBoxMenuItem("AutoSave");      co1.addItemListener(this);       menu.add(co1);
        co3  = new JCheckBoxMenuItem("DropSummary");   co3.addItemListener(this);       menu.add(co3);
        co4  = new JCheckBoxMenuItem("DumpGraphs");    co4.addItemListener(this);       menu.add(co4);
        co4b = new JCheckBoxMenuItem("DumpFiles");    co4b.addItemListener(this);       menu.add(co4b);
        co5  = new JCheckBoxMenuItem("DefaultGains");  co5.addItemListener(this);       menu.add(co5);
        co6  = new JCheckBoxMenuItem("FiduCuts");      co6.addItemListener(this);       menu.add(co6);
        co7  = new JCheckBoxMenuItem("DropEsect");     co7.addItemListener(this);       menu.add(co7);
        co8  = new JCheckBoxMenuItem("OnlyEsect");     co8.addItemListener(this);       menu.add(co8);
        menuBar.add(menu);
        
        if(!monitors[0].getDetectorName().equals("ECstatus")) co0.doClick(); //must use .equals here
        
        menu     = new JMenu("Fitting");
        cf  = new JCheckBoxMenuItem("Verbose");        cf.addItemListener(this);       menu.add(cf);
        cf0 = new JCheckBoxMenuItem("Calibration");   cf0.addItemListener(this);       menu.add(cf0);  
        cf1 = new JCheckBoxMenuItem("Residual");      cf1.addItemListener(this);       menu.add(cf1);  
        cf2 = new JCheckBoxMenuItem("TMF");           cf2.addItemListener(this);       menu.add(cf2);  
        cf3 = new JCheckBoxMenuItem("GTMF");          cf3.addItemListener(this);       menu.add(cf3);
        cf18= new JCheckBoxMenuItem("TRES");         cf18.addItemListener(this);       menu.add(cf18);
        menuBar.add(menu);
              
        menu     = new JMenu("Settings");       
        menuItem = new JMenuItem("Set GUI update interval");     menuItem.addActionListener(this); menu.add(menuItem);
        cf4      = new JCheckBoxMenuItem("log Y");                      cf4.addItemListener(this); menu.add(cf4);  
        cf5      = new JCheckBoxMenuItem("log Z");                      cf5.addItemListener(this); menu.add(cf5); cf5.doClick();
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
        
        menu     = new JMenu("Debug");
        cf7      = new JCheckBoxMenuItem("ECEngine"); cf7.addItemListener(this);   menu.add(cf7);  
        cf8      = new JCheckBoxMenuItem("Analyzer"); cf8.addItemListener(this);   menu.add(cf8); 
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
        ctr6 = new JRadioButtonMenuItem("Proton");   ctr6.addItemListener(this);       group.add(ctr6); menu.add(ctr6);
        ctr7 = new JRadioButtonMenuItem("ElecPion"); ctr7.addItemListener(this);       group.add(ctr7); menu.add(ctr7);      
        menuBar.add(menu);
       
        menu   	= new JMenu("ECstatus");
        cf9  = new JCheckBoxMenuItem("ATDATA");         cf9.addItemListener(this); menu.add(cf9);
        cf17 = new JCheckBoxMenuItem("StatEachRun");   cf17.addItemListener(this); menu.add(cf17);
        menuItem = new JMenuItem("Set MaxEvents"); menuItem.addActionListener(this); menu.add(menuItem);
        menuBar.add(menu);
        
        menu   	= new JMenu("ECcalib");
        cf10a = new JCheckBoxMenuItem("NormPix");cf10a.addItemListener(this); menu.add(cf10a);
        cf10b = new JCheckBoxMenuItem("NormAtt");cf10b.addItemListener(this); menu.add(cf10b);
        cf10c = new JCheckBoxMenuItem("NormCz") ;cf10c.addItemListener(this); menu.add(cf10c);
        menuBar.add(menu); 
        
        menu   	= new JMenu("ECsf");
        cf11 = new JCheckBoxMenuItem("SFcorr"); cf11.addItemListener(this); menu.add(cf11);
        cf12 = new JCheckBoxMenuItem("HiRes");  cf12.addItemListener(this); menu.add(cf12);
        menuBar.add(menu);
        
        menu   	= new JMenu("ECt");
        cf13 = new JCheckBoxMenuItem("TWcorr"); cf13.addItemListener(this); menu.add(cf13); cf13.doClick();
        menuBar.add(menu); 
        
        menu   	= new JMenu("ECpi0");
        cf16 = new JCheckBoxMenuItem("FTOFveto"); cf16.addItemListener(this); menu.add(cf16); cf16.doClick();
        menuBar.add(menu); 
        
        menu     = new JMenu("ECEngine");
        co2      = new JCheckBoxMenuItem("Enable");              co2.addItemListener(this);  menu.add(co2);;
        menuItem = new JMenuItem("Set logParam");                menuItem.addActionListener(this); menu.add(menuItem);   
        cf6a     = new JCheckBoxMenuItem("RejectSharedTime");    cf6a.addItemListener(this); menu.add(cf6a); cf6a.doClick(); 
        cf6b     = new JCheckBoxMenuItem("RejectSharedEnergy");  cf6b.addItemListener(this); menu.add(cf6b); cf6b.doClick();
        cf6c     = new JCheckBoxMenuItem("Use FADC time");       cf6c.addItemListener(this); menu.add(cf6c);
        cf14     = new JCheckBoxMenuItem("RepairMissingDT");     cf14.addItemListener(this); menu.add(cf14); cf14.doClick();
        cf15     = new JCheckBoxMenuItem("PCAL FTime");          cf15.addItemListener(this); menu.add(cf15); cf15.doClick();
        cf6d     = new JCheckBoxMenuItem("PASS 2 Timing");       cf6d.addItemListener(this); menu.add(cf6d); cf6d.doClick();
        cf6f     = new JCheckBoxMenuItem("PASS 2 Energy");       cf6f.addItemListener(this); menu.add(cf6f); cf6f.doClick();
        cf6e     = new JCheckBoxMenuItem("calibpass2");          cf6e.addItemListener(this); menu.add(cf6e); cf6e.doClick();
        cf6g     = new JCheckBoxMenuItem("ECHITS");              cf6g.addItemListener(this); menu.add(cf6g);
        cf6h     = new JCheckBoxMenuItem("ASA");                 cf6h.addItemListener(this); menu.add(cf6h);
        cf6i     = new JCheckBoxMenuItem("CC");                  cf6i.addItemListener(this); menu.add(cf6i);
        cf6j     = new JCheckBoxMenuItem("GPP");                 cf6j.addItemListener(this); menu.add(cf6j);
        menuItem = new JMenuItem("Set PC Z plane");              menuItem.addActionListener(this); menu.add(menuItem);   
        menuItem = new JMenuItem("Set EC Z plane");              menuItem.addActionListener(this); menu.add(menuItem);   
        menuItem = new JMenuItem("Set Hit Thresh");              menuItem.addActionListener(this); menu.add(menuItem);   
        menuItem = new JMenuItem("Set Peak Thresh");             menuItem.addActionListener(this); menu.add(menuItem);   
        menuBar.add(menu);
         
    }
    
    public void createPanels() {
    	System.out.println("EventViewer.createPanels()");

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
        this.dataProcessor.addEventListener(this);
        this.setCanvasUpdate(canvasUpdateTime);
        
    }
    
    @Override
    public void stateChanged(ChangeEvent e) {
        this.timerUpdate();
    }
    
    public Boolean sc(ItemEvent e) {return (e.getStateChange() == ItemEvent.SELECTED)?true:false;}
  
	@Override
	public void itemStateChanged(ItemEvent e) {
		Object s = e.getItemSelectable();	
		if (s==co0) {clearHist      = sc(e);} 
		if (s==co1) {autoSave       = sc(e); monitors[0].autoSave    = autoSave;}
		if (s==co2) {dropBanks      = sc(e); monitors[0].setDropBanks(dropBanks);}
		if (s==co3) {dropSummary    = sc(e); monitors[0].dropSummary = dropSummary;}
		if (s==co4) {dumpGraphs     = sc(e); monitors[0].dumpGraphs  = dumpGraphs;}
		if (s==co4b){dumpFiles      = sc(e); monitors[0].dumpFiles   = dumpFiles;}
		if (s==co6) {fiduCuts       = sc(e); monitors[0].fiduCuts    = fiduCuts;}
		if (s==co7) {dropEsect      = sc(e); monitors[0].dropEsect   = dropEsect;}
		if (s==co7) {onlyEsect      = sc(e); monitors[0].onlyEsect   = onlyEsect;}
		if (s==cf)  {fitVerbose     = sc(e); monitors[0].fitVerbose  = fitVerbose;}
		if (s==cf0) {cfitEnable     = sc(e); monitors[0].cfitEnable  = cfitEnable;}
		if (s==cf1) {sfitEnable     = sc(e); monitors[0].sfitEnable  = sfitEnable;}
		if (s==cf2) {dfitEnable     = sc(e); monitors[0].dfitEnable  = dfitEnable;}
		if (s==cf3) {gdfitEnable    = sc(e); monitors[0].gdfitEnable = gdfitEnable;}
		if (s==cf18){trfitEnable    = sc(e); monitors[0].trfitEnable = trfitEnable;}
		if (s==cf4) {yLogEnable     = sc(e); monitors[0].setLogY(yLogEnable);} 
		if (s==cf5) {zLogEnable     = sc(e); monitors[0].setLogZ(zLogEnable);}
		if (s==co5) {defaultGain    = sc(e); monitors[0].eng.setDefaultGain(defaultGain);}		
		if (s==cf6a){unsharedTime   = sc(e); monitors[0].eng.setUseUnsharedTime(unsharedTime);} 
		if (s==cf6b){unsharedEnergy = sc(e); monitors[0].eng.setUseUnsharedEnergy(unsharedEnergy);}
		if (s==cf6c){useFADCTime    = sc(e); monitors[0].eng.setUseFADCTime(useFADCTime);}
		if (s==cf6d){usePass2Timing = sc(e); monitors[0].eng.setUsePass2Timing(usePass2Timing);}
		if (s==cf6f){usePass2Energy = sc(e); monitors[0].eng.setUsePass2Energy(usePass2Energy);}
		if (s==cf6e){useCalibPass2  = sc(e); monitors[0].eng.setUseCalibPass2(useCalibPass2);}
		if (s==cf6g){outputECHITS   = sc(e); monitors[0].eng.outputECHITS(outputECHITS);}
		if (s==cf6h){useASA         = sc(e); monitors[0].eng.setUseASA1(useASA);}
		if (s==cf6i){useCC          = sc(e); monitors[0].eng.setUseCC(useCC);}
		if (s==cf6j){useGPP         = sc(e);{monitors[0].eng.setUseGPP(useGPP);monitors[0].useGPP = useGPP;}}
		if (s==cf7) {dbgECEngine    = sc(e); monitors[0].eng.setDbgECEngine(dbgECEngine);}
		if (s==cf8) {dbgAnalyzer    = sc(e); monitors[0].dbgAnalyzer = dbgAnalyzer;}
		if (s==cf9) {useATDATA      = sc(e); monitors[0].useATDATA   = useATDATA;}
		if (s==cf10a){normPix       = sc(e); monitors[0].normPix     = normPix;}
		if (s==cf10b){normAtt       = sc(e); monitors[0].setNormAtt(normAtt);}
		if (s==cf10c){normCz        = sc(e); monitors[0].normCz      = normCz;}
		if (s==cf11){SFcorr         = sc(e); monitors[0].SFcorr      = SFcorr;}
		if (s==cf12){HiRes          = sc(e); monitors[0].HiRes       = HiRes;}
		if (s==cf13){TWcorr         = sc(e); monitors[0].eng.setUseTWcorr(TWcorr);}
		if (s==cf14){DTcorr         = sc(e); monitors[0].eng.setUseDTcorr(DTcorr);}
		if (s==cf15){FTpcal         = sc(e); monitors[0].eng.setUseFTpcal(FTpcal);}
		if (s==cf16){FTOFveto       = sc(e); monitors[0].FTOFveto    = FTOFveto;}
		if (s==cf17){StatEachRun    = sc(e); monitors[0].StatEachRun = StatEachRun;}
		if (s==ct3) {TLflag         = sc(e); monitors[0].setTLflag(TLflag);}
		
		if (s==ct0)  {TLname = ct0.getText(); monitors[0].initTimeLine(TLname);}
		if (s==ct1)  {TLname = ct1.getText(); monitors[0].initTimeLine(TLname);}
		if (s==ct2)  {TLname = ct2.getText(); monitors[0].initTimeLine(TLname);}
		
		if (s==ctr0) {TRpid = sc(e)?   11:11; monitors[0].TRpid = TRpid;}
		if (s==ctr1) {TRpid = sc(e)?  211:11; monitors[0].TRpid = TRpid;}
		if (s==ctr2) {TRpid = sc(e)?   22:11; monitors[0].TRpid = TRpid;}
		if (s==ctr3) {TRpid = sc(e)? 2112:11; monitors[0].TRpid = TRpid;}
		if (s==ctr4) {TRpid = sc(e)?    0:11; monitors[0].TRpid = TRpid;}
		if (s==ctr5) {TRpid = sc(e)?   -1:11; monitors[0].TRpid = TRpid;}
		if (s==ctr6) {TRpid = sc(e)? 2212:11; monitors[0].TRpid = TRpid;}		
		if (s==ctr7) {TRpid = sc(e)?11211:11; monitors[0].TRpid = TRpid;}		
    }  
	
    public void actionPerformed(ActionEvent e) {
        System.out.println(e.getActionCommand());
        switch (e.getActionCommand()) {
          case("Load Run"):                    this.loadHistoFromRunIndex(); break;
//          case("Load Summary"):                this.readHistosFromSummary(); break;
          case("Analyze Runs"):                this.readFiles(); break;
          case("Analyze Histos"):              monitors[0].analyzeHistos = true; this.readHistos(); break;
          case("Set GUI update interval"):     this.setUpdateInterval(); break; 
          case("Set MaxEvents"):               this.setMaxEvents(); break;
          case("Set logParam"):                this.setLogParam(); break;
          case("Set PC Z plane"):              this.setPCZplane(); break;
          case("Set EC Z plane"):              this.setECZplane(); break;
          case("Set Hit Thresh"):              this.setHitThresh(); break;
          case("Set Peak Thresh"):             this.setPeakThresh(); break;
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
    
    public void initMonitors() {
    	System.out.println("EventViewer.initMonitors");
		for(int k=0; k<this.monitors.length; k++) {
			this.monitors[k].dropBanks   = dropBanks; 
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
			this.monitors[k].setTLflag(TLflag);
			this.monitors[k].setDbgAnalyzer(dbgAnalyzer);
    	} 
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

    public void setUpdateInterval() {
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
    
    public void setMaxEvents() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "Number of Events Analyzed",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "10002");
        if(s!=null){
            int val = 10002;
            try { 
                val= Integer.parseInt(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            if(val>0) {
            	this.setMaxEvents(val);
            }
            else {
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
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
    
    public void setLogParam() {
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
    
    public void setPCZplane() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "ECEngine pcTrackingPlane",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "0");
        if(s!=null){
            int val = 0;
            try { 
                val= Integer.parseInt(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            this.setPCTrackingPlane(val);
        }
    } 
    
    public void setECZplane() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "ECEngine ecTrackingPlane",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "0");
        if(s!=null){
            int val = 0;
            try { 
                val= Integer.parseInt(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            this.setECTrackingPlane(val);
        }
    } 
    
    public void setHitThresh() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "ECEngine stripThresholds",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "0");
        if(s!=null){
            String val = "0,0,0";
            try { 
                val= s;
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            this.setStripThreshold(val);
        }
    } 
    
    public void setPeakThresh() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "ECEngine peakThresholds",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "0");
        if(s!=null){
            String val = "0,0,0";
            try { 
                val= s;
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            this.setPeakThreshold(val);
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
    
    public long getTriggerWord(DataEvent event) {    	
 	    DataBank bank = event.getBank("RUN::config");          
        return bank.getLong("trigger", 0);
    }
  
    public void setTriggerPhaseConstants(int run) {
		IndexedTable it = this.monitors[0].cm.getConstants(run, "/calibration/ec/time_jitter");        
        PERIOD = it.getDoubleValue("period",0,0,0);
        PHASE  = it.getIntValue("phase",0,0,0); 
        CYCLES = it.getIntValue("cycles",0,0,0);
    }
    
    public int getTriggerPhase(DataEvent event) {    	
 	    DataBank bank = event.getBank("RUN::config");	        
        long timestamp = bank.getLong("timestamp",0);    
        if (CYCLES==0) return 0;
        return (int) (timestamp>0 ? (PERIOD*((timestamp+PHASE)%CYCLES)):0); // TI derived phase correction due to TDC and FADC clock differences 
    }
    
    private int getRunNumber(DataEvent event) {
        if(this.ccdbRunNumber>0) return this.ccdbRunNumber;
        return event.hasBank("RUN::config") ? event.getBank("RUN::config").getInt("run",0):this.runNumber;
    }
    
    private int getEventNumber(DataEvent event) {
        return event.hasBank("RUN::config") ? event.getBank("RUN::config").getInt("event",0):this.eventNumber;
    }
    
    private int getTotalEvents() {
    	return processorPane.getNevents();
    }
    
    private boolean rejectEvent(DataEvent event) {
    	if(event==null) return true;
    	if(event.hasBank("RUN::scaler") || event.hasBank("RAW::scaler") || event.hasBank("HEL::flip")) return true;
    	return false;
    }
    
    private boolean isFirstRun(int run) {
    	return runList.isEmpty();
    }
    
    private boolean isNewRun(int run) {
    	if(isFirstRun(run)) {runList.add(run); return true;}
    	if(runList.contains(run)) return false;
    	runList.add(run);
    	return true;
    }
    
    @Override
    public void dataEventAction(DataEvent event) {   	
	   if(!rejectEvent(event))  processEvent(filterEvent(decodeEvent(event)));
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
    
    private DataEvent filterEvent(DataEvent event) {

    	int rNum = getRunNumber(event);
    	
        if(rNum!=0 && isNewRun(rNum) && clear) {//clear is initialized true (except ECscaler) and reset true by readFiles()
        	System.out.println("\nEventViewer: Processing Run "+rNum+" Event: "+getEventNumber(event));
        	this.runNumber = rNum;
            if(!clearHist) clear=false; //bypass initRun after first run is analyzed
        	return initRun(rNum,event);
        }        
        return event;
    }
    
    private void processEvent(DataEvent event) {
    	
        setTriggerPhaseConstants(this.runNumber);

        for(int k=0; k<this.monitors.length; k++) {
        	this.monitors[k].setEventNumber(getEventNumber(event));
        	this.monitors[k].setTriggerPhase(getTriggerPhase(event));
            this.monitors[k].setTriggerWord(getTriggerWord(event));        	    
            this.monitors[k].dataEventAction(event);
        }  
    }
    
    private DataEvent initRun(int runno, DataEvent event) {
    	System.out.println("\nEventViewer.initRun("+runno+")");
        for(int k=0; k<this.monitors.length; k++) {
        	if(autoSave && this.runNumber!=0) this.monitors[k].saveHistosToFile();
            this.runNumber = runno; 
            this.monitors[k].isHipo3Event = processorPane.isHipo3Event;
        	this.monitors[k].setRunNumber(this.runNumber); 
            this.monitors[k].setTotalEvents(getTotalEvents());
           	this.monitors[k].localclear();
           	this.monitors[k].initCCDB(this.runNumber);
           	this.monitors[k].initEBCCDB(this.runNumber);
        	this.monitors[k].createHistos(this.runNumber);
            this.monitors[k].initGStyle(18);
            this.monitors[k].plotHistos(this.runNumber);
            this.monitors[k].arBtn.setSelected(true);    
            if(this.monitors[k].sectorButtons) this.monitors[k].bS2.doClick();
        }         
        return event;
    }
    
    private void readFiles() {
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
                        HipoDataSource   source = new HipoDataSource();
                        source.open(fd);
                        this.dataProcessor.setSource(source);
                        current = source.getCurrentIndex();
                        nevents = source.getSize();
                        System.out.println("\nFILE: " + nf + " " + fd.getName() + " N.EVENTS: " + nevents.toString() + "  CURRENT : " + current.toString());
                        boolean hasFinished = false;
                        while (hasFinished==false) {
                        	for (int i=1 ; i<=50 ; i++) {
                        		boolean status = dataProcessor.processNextEvent(eventDelay,DataEventType.EVENT_ACCUMULATE);
                        		if(status==false && hasFinished==false){
                        			hasFinished = true;
                        			System.out.println("[DataProcessingPane] ----> task is done...");
                        		}
                        	}
                        }
//                        this.startProcessorTimer();
//                        int k=0;
//                        while (source.hasEvent()) {
//                        	event = source.getNextEvent();                                      
//                            this.dataEventAction(event); k++;
//                            if(k % 10000 == 0) System.out.println("Read " + k + " events");
//                        }
                        source.close(); clear=true;
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
    
//  Analyze Histos
    
    private void readHistos() {
    	String fname = null;
        for(int k=0; k<this.monitors.length; k++) this.monitors[k].localclear();

        for (File fd : selectHistos()) {
            if (fd.isFile()) {
                fname=fd.getAbsolutePath();
                loadHistosFromFile(fname);
            }                
        }  
        this.monitors[0].writeScript(monitors[0].getDetectorName());
        if(isCalibrationFile(fname)) {
        	this.monitors[0].writeFile(getFileCalibrationTag(fname),1,7,0,3,0,3);
        }
        monitors[0].analyzeHistos = false;
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
        
    public void loadHistosFromFile(String fileName) {
        System.out.println("EventViewer.loadHistosFromFile("+fileName+")");
        
        runNumber = getFileRunNumber(fileName); 
        
        TDirectory dir = new TDirectory(); dir.readFile(fileName); dir.cd(); dir.pwd();
        processorPane.setHistoName(fileName);

        for(int k=0; k<this.monitors.length; k++) {  
            if(isCalibrationFile(fileName)) this.monitors[k].detcal[getFileCalibrationIndex(fileName)]=runNumber;
         	this.monitors[k].setRunNumber(runNumber);
            this.monitors[k].createHistos(runNumber);  
            this.monitors[k].initGStyle(18);
            this.monitors[k].initCCDB(runNumber);
            this.monitors[k].initEBCCDB(runNumber);
            this.monitors[k].readDataGroup(runNumber,dir);
            if (this.monitors[k].sectorButtons) {this.monitors[k].bS2.doClick();}
            this.monitors[k].arBtn.setSelected(true);   
        }
        return;
    }
             
    public void loadHistoFromRunIndex() {
    	File[] f = fc.getSelectedFiles();
    	monitors[0].dropSummary=false;
    	loadHistosFromFile(f[monitors[0].getRunIndex()].getAbsolutePath());
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
        System.out.println("EventViewer.setTVOffset("+time+")");
        this.TVOffset = time;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setTVOffset(time);
        }
    }
    
    public void setLogParam(float val) {
        System.out.println("EventViewer.setLogParam("+val+")");
        this.logParam = val;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].eng.setLogParam(val);
        }
    }
    
    public void setPCTrackingPlane(int val) {
        System.out.println("EventViewer.setPCTrackingPlane("+val+")");
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].eng.setPCTrackingPlane(val);
        }
    }
    
    public void setECTrackingPlane(int val) {
        System.out.println("EventViewer.setECTrackingPlane("+val+")");
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].eng.setECTrackingPlane(val);
        }
    }
    
    public void setStripThreshold(String val) {
        System.out.println("EventViewer.setStripThresholds("+val+")");
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].eng.setStripThreshold(val);
        }
    }
    
    public void setPeakThreshold(String val) {
        System.out.println("EventViewer.setPeakThresholds("+val+")");
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].eng.setPeakThreshold(val);
        }
    }
    
    public void setMaxEvents(int val) {
        System.out.println("EventViewer.setMaxEvents("+val+")");
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setMaxEvents(val);
        }
    } 
    
    @Override
    public void timerUpdate() {
        if(this.runNumber==0) return;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].timerUpdate();
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
    
    private void startProcessorTimer(){
        class CrunchifyReminder extends TimerTask {
            boolean hasFinished = false;
            public void run() {
                //dataProcessor.processNextEvent(0, DataEventType.EVENT_START);
                /*if(hasFinished==true){
                    dataProcessor.processNextEvent(0, DataEventType.EVENT_STOP);
                    return;
                }*/
                //System.out.println("running");
                for (int i=1 ; i<=50 ; i++) {
                    boolean status = dataProcessor.processNextEvent(eventDelay,DataEventType.EVENT_ACCUMULATE);
                    if(status==false&&hasFinished==false){
                        hasFinished = true;
                        System.out.println("[DataProcessingPane] ----> task is done...");
                    }
                }
//                statusLabel.setText(dataProcessor.getStatusString());
            }
        }
        processTimer = new java.util.Timer();
        processTimer.schedule(new CrunchifyReminder(),1,1);        
    }    

	@Override
	public void processShape(DetectorShape2D arg0) {
		// TODO Auto-generated method stub
		
	}

}