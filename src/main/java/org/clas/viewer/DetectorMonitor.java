package org.clas.viewer;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JSplitPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileSystemView;

import org.clas.tools.FitData;
import org.jlab.detector.base.DetectorOccupancy;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.RangeSlider;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.io.task.IDataEventListener;
import org.jlab.rec.eb.EBCCDBConstants;
import org.jlab.service.eb.EBEngine;
import org.jlab.service.eb.EventBuilder;
import org.jlab.service.ec.ECEngine;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;


public class DetectorMonitor implements ActionListener {    
    
    private final String                detectorName;
    private ArrayList<String>       detectorTabNames = new ArrayList();
    private IndexedList<DataGroup>      detectorData = new IndexedList<DataGroup>(4);
    public  List<Integer>                    runlist = new ArrayList<Integer>();
    public  int                       runIndexSlider = 0;
    public  JSlider                           slider = null;
    private DataGroup                detectorSummary = null;
    private DetectorOccupancy      detectorOccupancy = new DetectorOccupancy();
    private JPanel                     detectorPanel = null;
    private JPanel                       actionPanel = null;
    private JPanel                    controlsPanel0 = null;
    private JPanel                    controlsPanel1 = null;
    private JPanel                    controlsPanel2 = null;
    private JPanel                     runIndexPanel = null;
    private EmbeddedCanvasTabbed      detectorCanvas = null;
    private DetectorPane2D              detectorView = null;
    private ButtonGroup                          bT0 = null;
    private ButtonGroup                          bG0 = null;
    private ButtonGroup                          bG1 = null;
    private ButtonGroup                          bG2 = null;
    private ButtonGroup                          bG3 = null;
    private int                       numberOfEvents = 0;
    public  Boolean                    sectorButtons = false;
    private Boolean                       sliderPane = false;
    private int                     detectorActivePC = 1;
    private int                 detectorActiveSector = 1;
    private int                   detectorActiveView = 0;
    private int                  detectorActiveLayer = 0;
    private int                    detectorActivePID = 0;
    private Boolean                     detectorLogZ = true;
    private Boolean                             isTB = false;
    private Boolean                            usePC = false;
    private Boolean                           usePID = false;
    private Boolean                        stopBlink = true;
    Timer timer = null;
    
    public JRadioButton bEL,bPI,bPH,bP,bC,bS1,bS2,bS3,bS4,bS5,bS6,bpcal,becin,becou,bu,bv,bw;
    private JCheckBox tbBtn;
    public  JCheckBox arBtn;
    
    public int        bitsec = 0;
    public long      trigger = 0;
    public int  triggerPhase = 0;
    public int        trigFD = 0;
    public int        trigCD = 0;
    
    public boolean   testTrigger = false;
    public boolean TriggerBeam[] = new boolean[32];
    public int       TriggerMask = 0;
    
    private int   runNumber = 0;
    private int eventNumber = 0;
    private int     viewRun = 0;  
    
    public int eventResetTime_current[]=new int[19];
    public int eventResetTime_default[]=new int[19];    
    
    public double lMin=0, lMax=500, zMin=0.1, zMax=1.0, zMinLab, zMaxLab;
    public double slideMin=1, slideMax=500;
    public  Boolean doAutoRange = false;
    
    public String[] layer = new String[]{"pcal","ecin","ecou"};
    public String[]  view = new String[]{"u","v","w"};    
    
    public ECEngine  engine = new ECEngine();  
	public EBEngine     ebe = new EBEngine("CLAS12ANA");

    public String variation = "default";
    public String      geom = "2.5";
    public String    config = "muon";   
	public EventBuilder eb = null;
    
    int[][] sthrMuon = {{15,15,15},{20,20,20},{20,20,20}};
    int[][] sthrPhot = {{10,10,10},{9,9,9},{8,8,8}};
    int[][] sthrElec = {{10,10,10},{10,10,10},{10,10,10}};
    int[][] sthrZero = {{1,1,1},{1,1,1},{1,1,1}};
    
    int[][] pthrMuon = {{15,15,15},{20,20,20},{20,20,20}};
    int[][] pthrPhot = {{18,18,18},{20,20,20},{15,15,15}};
    int[][] pthrElec = {{30,30,30},{30,30,30},{30,30,30}};
    int[][] pthrZero = {{1,1,1},{1,1,1},{1,1,1}};
        
    double[] cerrMuon = {5.5,10.,10.};
    double[] cerrPhot = {7,15.,20.};
    double[] cerrElec = {10.,10.,10.};  
    
    public String  outPath = "/home/lcsmith/CLAS12ANA/";
    public String  pawPath = outPath+"paw/";
    public String  jawPath = outPath+"jaw/";
    public String  vecPath = null;
    
    public Boolean isAnalyzeDone = false;
    public Boolean      autoSave = false;
    public Boolean     dropBanks = false;
    public Boolean   dropSummary = false; 
    public Boolean    dumpGraphs = false; 
    public Boolean     fitEnable = false; 
    public Boolean    fitVerbose = false; 
    
    public String                 TLname = null;    
    public Map<String,Integer> TimeSlice = new HashMap<String,Integer>();  
    public List<Integer>       BlinkRuns = new ArrayList<Integer>();
    int ntimer = 0;
    
    public DetectorMonitor(String name){
        initGStyle();
        detectorName   = name;
        detectorPanel  = new JPanel();
        detectorCanvas = new EmbeddedCanvasTabbed();
        detectorView   = new DetectorPane2D();
        
        numberOfEvents = 0;        
        eventResetTime_current[0]=0;
        eventResetTime_current[1]=0;     
        eventResetTime_default[0]=10000000;  
        eventResetTime_default[1]=10000000;       
        for (int i=0; i<2; i++){
            eventResetTime_current[i]=eventResetTime_default[i];
        }
        outPath = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/CLAS12ANA/";
        System.out.println(detectorName+" outPath = "+outPath);
        getEnv();
        pawPath = outPath+"paw/";
        vecPath = pawPath+detectorName+"/";  
        jawPath = outPath+"jaw/";
        TimeSlice.put("UVW", 3);
        TimeSlice.put("FADC Slot", 16);
        TimeSlice.put("HV Slot", 24);
    }
    
    public void getEnv() {        
        String ostype = System.getProperty("os.name"); 
        System.out.println("DetectorMonitor.getEnv(): os.name = "+ostype);
        
        if (ostype!=null&&ostype.startsWith("Mac")) {
            outPath = "/Users/cole/CLAS12ANA/";
        } else {
            outPath = "/home/lcsmith/CLAS12ANA/";            	
        }

    }
    
    public void init() {
        initPanel();    	
    }
    
    public void initTimeLine(String name) {
    	TLname = name;
    	if(getRunNumber()!=0) plotHistos(getRunNumber());
    }
    
    public void createTimeLineHistos() {
    	
    }
    
    public void localinit() {
    	
    }
    
    public void localclear() {
    	
    }
    
    public void initGStyle() {
        GStyle.getAxisAttributesX().setTitleFontSize(18);
        GStyle.getAxisAttributesX().setLabelFontSize(14);
        GStyle.getAxisAttributesY().setTitleFontSize(18);
        GStyle.getAxisAttributesY().setLabelFontSize(18);
        GStyle.getAxisAttributesZ().setLabelFontSize(14); 
        GStyle.getAxisAttributesX().setAxisGrid(false);
        GStyle.getAxisAttributesY().setAxisGrid(false);
        GStyle.getAxisAttributesX().setLabelFontName("Avenir");
        GStyle.getAxisAttributesY().setLabelFontName("Avenir");
        GStyle.getAxisAttributesZ().setLabelFontName("Avenir");
        GStyle.getAxisAttributesX().setTitleFontName("Avenir");
        GStyle.getAxisAttributesY().setTitleFontName("Avenir");
        GStyle.getAxisAttributesZ().setTitleFontName("Avenir");
        GStyle.getAxisAttributesZ().setAxisAutoScale(true);    	
        GStyle.getGraphErrorsAttributes().setMarkerStyle(1);
        GStyle.getGraphErrorsAttributes().setMarkerColor(2);
        GStyle.getGraphErrorsAttributes().setMarkerSize(2);
        GStyle.getGraphErrorsAttributes().setLineColor(2);
        GStyle.getGraphErrorsAttributes().setLineWidth(1);
        GStyle.getGraphErrorsAttributes().setFillStyle(1);   
    }
    
    public void initPanel() {
        getDetectorPanel().setLayout(new BorderLayout());
        actionPanel = new JPanel(new FlowLayout());
        runIndexPanel = new JPanel(new FlowLayout());
        controlsPanel0 = new JPanel(new GridBagLayout());
        this.controlsPanel1 = new JPanel();
        this.controlsPanel1.setBorder(BorderFactory.createTitledBorder("Display"));		
        this.controlsPanel2 = new JPanel();
        this.controlsPanel2.setBorder(BorderFactory.createTitledBorder("Timeline"));        
        controlsPanel1.add(packActionPanel());
        controlsPanel2.add(packRunIndexPanel());
        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL; c.weightx = 0.5;		
        c.gridx=0 ; c.gridy=0 ; controlsPanel0.add(controlsPanel1,c);
        c.gridx=0 ; c.gridy=1 ; controlsPanel0.add(controlsPanel2,c);
        getDetectorPanel().add(getDetectorCanvas(),BorderLayout.CENTER);           
        getDetectorPanel().add(controlsPanel0,BorderLayout.SOUTH); 
//        getDetectorPanel().add(packRunIndexPanel(),BorderLayout.SOUTH); 
    }
    
    public void configEngine(String config) {
    	engine.isSingleThreaded=true;
        engine.setVariation(variation);
        engine.init();
        engine.isMC = false;
       
        engine.setStripThresholds(getStripThr(config, 0, 1),
                                  getStripThr(config, 1, 1),
                                  getStripThr(config, 2, 1));  
        engine.setPeakThresholds(getPeakThr(config, 0, 1),
                                 getPeakThr(config, 1, 1),
                                 getPeakThr(config, 2, 1));  
        engine.setClusterCuts(getClusterErr(config,0),
                              getClusterErr(config,1),
                              getClusterErr(config,2));   
    }
    
    public void configEventBuilder() {
    	ebe.init();
        eb = new EventBuilder(new EBCCDBConstants(10,ebe.getConstantsManager()));    	    	
    }
    
    public int getStripThr(String config, int idet, int layer) {
        switch (config) {
        case     "pi0": return sthrPhot[idet][layer-1] ;  
        case    "phot": return sthrPhot[idet][layer-1] ; 
        case    "muon": return sthrMuon[idet][layer-1] ;  
        case    "elec": return sthrElec[idet][layer-1] ;
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
        case    "none": return (float) cerrMuon[idet] ;
        }
        return 0;
     } 
    
    public void dropBanks(DataEvent event) {
        if(event.hasBank("ECAL::hits")) {
            event.removeBank("ECAL::hits");        
            event.removeBank("ECAL::peaks");        
            event.removeBank("ECAL::clusters");        
            event.removeBank("ECAL::calib");
            event.removeBank("ECAL::moments");
         } 
        if(event.hasBank("ECAL::adc")) engine.processDataEvent(event);     	
    }
    
    public void analyze() {
        // analyze detector data at the end of data processing
    }
    
    public void fillSummary() {
       
    }
    
    public void createHistos(int run) {
        // initialize canvas and create histograms
    }
    
    public void createSummary() {
        // initialize canvas and create histograms
    } 
    
    public void dataEventAction(DataEvent event) {
        if (!testTriggerMask()) return;
        setNumberOfEvents(getNumberOfEvents()+1);
        switch (event.getType())  {
        case EVENT_START:      processEvent(event); break;
        case EVENT_SINGLE:    {processEvent(event); plotEvent(event);break;}
        case EVENT_ACCUMULATE: processEvent(event); break;
        case EVENT_STOP:       analyze(); if(autoSave) saveHistosToFile();
	    }
    }

    public void drawDetector() {
    
    }
    
    public void setTriggerPhase(int phase) {
    	   this.triggerPhase = phase;
    }
    
    public int getTriggerPhase() {
    	    return this.triggerPhase;
    }
    
    public void setTriggerWord(long trig) {
    	   trigger = trig;
    }
    
    public void setTestTrigger(boolean test) {
    	   testTrigger = test;
    }
    
    public int     getFDTrigger()            {return (int)(trigger)&0x000000000ffffffff;}
    public int     getCDTrigger()            {return (int)(trigger>>32)&0x00000000ffffffff;}
    public boolean isGoodFD()                {return  getFDTrigger()>0;}    
    public boolean isTrigBitSet(int bit)     {int mask=0; mask |= 1<<bit; return isTrigMaskSet(mask);}
    public boolean isTrigMaskSet(int mask)   {return (getFDTrigger()&mask)!=0;}
    public boolean isGoodECALTrigger(int is) {return (testTrigger)? is==getECALTriggerSector():true;}    
    public int           getElecTrigger()    {return getFDTrigger()&0x1;}
    public int     getElecTriggerSector()    {return (int) (isGoodFD() ? Math.log10((getFDTrigger()&0x7e)>>1)/0.301+1:0);} 
    public int     getECALTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>19)/0.301+1:0);}       
    public int     getPCALTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>13)/0.301+1:0);}       
    public int     getHTCCTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>7)/0.301+1:0);} 
    
    public int    getTriggerMask()        {return TriggerMask;}
    public void   setTriggerMask(int bit) {TriggerMask|=(1<<bit);}  
    public void clearTriggerMask(int bit) {TriggerMask&=~(1<<bit);}  
    public boolean testTriggerMask()      {return TriggerMask!=0 ? isTrigMaskSet(TriggerMask):true;}
    public boolean isGoodTrigger(int bit) {return TriggerBeam[bit] ? isTrigBitSet(bit):true;}
   
    public EmbeddedCanvasTabbed getDetectorCanvas() {
        return detectorCanvas;
    }
    
    public ArrayList<String> getDetectorTabNames() {
        return detectorTabNames;
    }
    
    public IndexedList<DataGroup>  getDataGroup(){
        return detectorData;
    }
    
    public void fillTimeLineGraph() {
    	
    }
    
    public void plotTimeLine(int index) {
    	
    }

    public String getDetectorName() {
        return detectorName;
    }
    
    public DetectorOccupancy getDetectorOccupancy() {
        return detectorOccupancy;
    }
    
    public JPanel getDetectorPanel() {
        return detectorPanel;
    }
    
    public DataGroup getDetectorSummary() {
        return detectorSummary;
    }
    
    public DetectorPane2D getDetectorView() {
        return detectorView;
    }
    
    public void useSectorButtons(boolean flag) {
    	sectorButtons = flag;
    }
    
    public void useSliderPane(boolean flag) {
	    sliderPane = flag;
    }    
    
    public void usePCCheckBox(boolean flag) {
        usePC = flag;
    }
    
    public void usePIDCheckBox(boolean flag) {
        usePID = flag;
    }
    
    public int getActivePC() {
    	    return detectorActivePC;
    }
    
    public int getActiveSector() {
	    return detectorActiveSector;
    }
    
    public int getActiveLayer() {
	    return detectorActiveLayer;
    }
    
    public int getActiveView() {
	    return detectorActiveView;
    }
    
    public int getNumberOfEvents() {
        return numberOfEvents;
    }
    
    public void setLogZ(boolean flag) {
	    detectorLogZ = flag;
	    plotHistos(getRunNumber());
    }
    
    public Boolean getLogZ() {
	    return detectorLogZ;
    }
    
    public JPanel packActionPanel() {
        if (sectorButtons) actionPanel.add(getButtonPane()); 
        if    (sliderPane) actionPanel.add(getSliderPane());
//                           actionPanel.add(getRunSliderPane());
    	return actionPanel;
    }
    
    public JPanel packRunIndexPanel() {
        runIndexPanel.add(getRunSliderPane());	
        return runIndexPanel;
    }
    
    public JPanel getButtonPane() {
        JPanel buttonPane = new JPanel();
                
        if(usePC) {
        bP = new JRadioButton("P"); buttonPane.add(bP); bP.setActionCommand("1"); bP.addActionListener(this);
        bC = new JRadioButton("C"); buttonPane.add(bC); bC.setActionCommand("0"); bC.addActionListener(this); 
        bG0 = new ButtonGroup(); bG0.add(bP); bG0.add(bC);
        bP.setSelected(true);
        }
        
        if(usePID) {
        bEL = new JRadioButton("e-"); buttonPane.add(bEL); bEL.setActionCommand("1"); bP.addActionListener(this);
        bPI = new JRadioButton("pi"); buttonPane.add(bPI); bPI.setActionCommand("2"); bPI.addActionListener(this); 
        bPH = new JRadioButton("ph"); buttonPane.add(bPH); bPH.setActionCommand("2"); bPH.addActionListener(this); 
        bT0 = new ButtonGroup(); bT0.add(bEL); bT0.add(bPI); bT0.add(bPH);
        bEL.setSelected(true);
        }   
        
        bS1 = new JRadioButton("S1"); buttonPane.add(bS1); bS1.setActionCommand("1"); bS1.addActionListener(this);
        bS2 = new JRadioButton("S2"); buttonPane.add(bS2); bS2.setActionCommand("2"); bS2.addActionListener(this); 
        bS3 = new JRadioButton("S3"); buttonPane.add(bS3); bS3.setActionCommand("3"); bS3.addActionListener(this); 
        bS4 = new JRadioButton("S4"); buttonPane.add(bS4); bS4.setActionCommand("4"); bS4.addActionListener(this); 
        bS5 = new JRadioButton("S5"); buttonPane.add(bS5); bS5.setActionCommand("5"); bS5.addActionListener(this);  
        bS6 = new JRadioButton("S6"); buttonPane.add(bS6); bS6.setActionCommand("6"); bS6.addActionListener(this); 
	    bG1 = new ButtonGroup(); bG1.add(bS1);bG1.add(bS2);bG1.add(bS3);bG1.add(bS4);bG1.add(bS5);bG1.add(bS6);
        bS2.setSelected(true);        
        bpcal = new JRadioButton("PC");  buttonPane.add(bpcal); bpcal.setActionCommand("0"); bpcal.addActionListener(this);
        becin = new JRadioButton("ECi"); buttonPane.add(becin); becin.setActionCommand("1"); becin.addActionListener(this); 
        becou = new JRadioButton("ECo"); buttonPane.add(becou); becou.setActionCommand("2"); becou.addActionListener(this); 
        bG2 = new ButtonGroup(); bG2.add(bpcal); bG2.add(becin); bG2.add(becou);
        bpcal.setSelected(true);
        bu = new JRadioButton("U"); buttonPane.add(bu); bu.setActionCommand("0"); bu.addActionListener(this);
        bv = new JRadioButton("V"); buttonPane.add(bv); bv.setActionCommand("1"); bv.addActionListener(this); 
        bw = new JRadioButton("W"); buttonPane.add(bw); bw.setActionCommand("2"); bw.addActionListener(this); 
        bG3 = new ButtonGroup(); bG3.add(bu); bG3.add(bv); bG3.add(bw);
        bu.setSelected(true);                
        return buttonPane;
    } 
    
    public JPanel getRunSliderPane() {
        JPanel sliderPane = new JPanel();
                   slider = new JSlider(JSlider.HORIZONTAL, 0, 450,0); 
        JLabel      label = new JLabel("" + String.format("%d", 0));
        sliderPane.add(new JLabel("Run Index,Run",JLabel.CENTER));
        sliderPane.add(slider);
        sliderPane.add(label);
        slider.setPreferredSize(new Dimension(1200,10));
        slider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                JSlider slider = (JSlider) e.getSource(); 
                runIndexSlider = (slider.getValue()>=runlist.size())?runlist.size()-1:slider.getValue();
                int run = runlist.get(runIndexSlider);
                label.setText(String.valueOf(""+String.format("%d", runIndexSlider)+" "+String.format("%d", run)));
                plotHistos(run);
            }
        });     
        JButton button = new JButton("Blink Run");
        button.addActionListener(this);
        sliderPane.add(button);
        return sliderPane;
    }
    
    public JPanel getSliderPane() {
        JPanel sliderPane = new JPanel();
        JLabel xLabel = new JLabel("Z-Range:");
        RangeSlider rslider = new RangeSlider();
        rslider.setMinimum((int) slideMin);
        rslider.setMaximum((int) slideMax);
        rslider.setValue((int)slideMin);
        rslider.setUpperValue((int)slideMax);            
        zMin =     rslider.getValue();
        zMax = 0.1*rslider.getUpperValue();
        zMinLab = Math.pow(2, zMin/10); zMaxLab = Math.pow(10, zMax/10);
        JLabel rangeSliderValue1 = new JLabel("" + String.format("%4.0f", zMinLab));
        JLabel rangeSliderValue2 = new JLabel("" + String.format("%4.0f", zMaxLab));
        sliderPane.add(xLabel);
        sliderPane.add(rangeSliderValue1);
        sliderPane.add(rslider);
        sliderPane.add(rangeSliderValue2);           
        rslider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                RangeSlider lslider = (RangeSlider) e.getSource();   
                lMin = lslider.getValue(); lMax = lslider.getUpperValue();
                zMax = 0.1*lMax; zMin = Math.min(lMin,zMax); 
                zMaxLab = Math.pow(10, zMax/10); zMinLab = Math.pow(2, zMin/10); 
                rangeSliderValue1.setText(String.valueOf("" + String.format("%4.0f", zMinLab)));
                rangeSliderValue2.setText(String.valueOf("" + String.format("%4.0f", zMaxLab)));
                plotHistos(getRunNumber());
            }
        });  
        arBtn = new JCheckBox("AutoRange");
        arBtn.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                doAutoRange = (e.getStateChange() == ItemEvent.SELECTED) ? true:false;
                plotHistos(getRunNumber());
            }
        });         
        arBtn.setSelected(false);         
        sliderPane.add(arBtn);        
        return sliderPane;
    }  
    
    public void actionPerformed(ActionEvent e) {
    	if(e.getActionCommand().compareTo("Blink Run")==0)  BlinkRunAction();
    	if(bG0!=null) detectorActivePC     = Integer.parseInt(bG0.getSelection().getActionCommand()); 
        if(bT0!=null) detectorActivePID    = Integer.parseInt(bT0.getSelection().getActionCommand()); 
        detectorActiveSector = Integer.parseInt(bG1.getSelection().getActionCommand());
        detectorActiveLayer  = Integer.parseInt(bG2.getSelection().getActionCommand());
        detectorActiveView   = Integer.parseInt(bG3.getSelection().getActionCommand());
        plotHistos(getRunNumber());
    } 
    
    public void BlinkRunAction() {
    	BlinkRuns.add(slider.getValue());
    	if(BlinkRuns.size()==2) {stopBlink=false; doBlink();}
    	if(BlinkRuns.size()==3) {stopBlink=true;  timer.cancel(); timer=null; BlinkRuns.clear();}
    }
    
    public void doBlink() {
     	TimerTask task;   
    	task = new TimerTask() {
    		@Override
    	    public synchronized void run() {     
    	        if(!stopBlink) {
    	           ntimer=(ntimer==0)?1:0;   	          
    	           slider.setValue(BlinkRuns.get(ntimer));
    	        }
    		}
    	};
    	timer = new Timer(); timer.schedule(task,1000,500);
    }
    
    public void processEvent(DataEvent event) {
        // process event
    }
    
    public void plotEvent(DataEvent event) {
        // process event
    }
    
    public void plotHistos(int run) {

    }
        
    public void saveHistosToFile() {
    	String fileName = "CLAS12Ana_run_" + getRunNumber() + "_" + getDetectorName() + ".hipo";
        TDirectory dir = new TDirectory();
        writeDataGroup(dir);
        dir.writeFile(outPath+fileName);
        System.out.println("Saving histograms to file " + outPath+fileName); 
    }
    
    public void printCanvas(String dir) {
        // print canvas to files
        int run = getViewRun();
        if(run==0) run = getRunNumber();
        for(int tab=0; tab<detectorTabNames.size(); tab++) {
            String fileName = dir + "/" + detectorName + "_" + run + "_canvas" + tab + ".png";
            System.out.println(fileName);
            detectorCanvas.getCanvas(detectorTabNames.get(tab)).save(fileName);
        }
    }   
    
    public void resetEventListener() {
        System.out.println("Resetting " +  getDetectorName() + " histogram for run "+ getRunNumber());
        createHistos(getRunNumber());
        plotHistos(getRunNumber());
    }
    
    public void setCanvasUpdate(int time) {
        for(int tab=0; tab<detectorTabNames.size(); tab++) {
            detectorCanvas.getCanvas(detectorTabNames.get(tab)).initTimer(time);
        }
    }
    
    public void setDetectorCanvas(EmbeddedCanvasTabbed canvas) {
        detectorCanvas = canvas;
    }
    
    public void setDetectorTabNames(String... names) {
        for(String name : names) {
            detectorTabNames.add(name);
        }
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed(names);
        setDetectorCanvas(canvas);
    }
 
    public void setDetectorSummary(DataGroup group) {
        detectorSummary = group;
    }
    
    public void setNumberOfEvents(int num) {
       numberOfEvents = num;
    }
    
    public void setEventNumber(int num) {
       eventNumber = num;
    }

    public void setRunNumber(int num) {
       runNumber = num;
    }
    
    public int getEventNumber() {
        return eventNumber;
    }

    public int getRunNumber() {
        return runNumber;
    }
    
    public float getTorusPolarity(int run) {
    	if (run>=3029&&run<=3065) return -1.00f;
    	if (run>=3072&&run<=3087) return -0.75f;
    	if (run>=3097&&run<=3105) return  0.75f;
    	if (run>=3131&&run<=3293) return  1.00f;
    	if (run>=3304&&run<=3817) return -1.00f;
    	if (run>=3819&&run<=3834) return  0.75f;
    	if (run>=3839&&run<=3987) return  1.00f;
    	if (run>=3995&&run<=4326) return -1.00f;  
    	return 0f;
    }

    public int getViewRun() {
        return viewRun;
    }   
     
    public void timerUpdate() { 
    }
        
    public void dumpGraph(String filename, GraphErrors graph) {
    	PrintWriter writer = null;
		try 
		{
			writer = new PrintWriter(filename);
		} 
		catch (FileNotFoundException e) 
		{			// TODO Auto-generated catch block
			e.printStackTrace();
		}				
    	for (int i=0; i<graph.getDataSize(0); i++) writer.println(String.format("%1$.3f %2$.6f %3$.8f",graph.getDataX(i),graph.getDataY(i),graph.getDataEY(i))); 
    	writer.close();    	
    }
    
    
    //Generalization of H2F method rebinY for list of non-equal ngroups. 
    public H2F rebinY(H2F h, int... ngroups) {
	    List<Integer> ig = new ArrayList<Integer>();
        for (int i: ngroups) ig.add(i);
        int nbinsY = (ngroups.length==1)?h.getYAxis().getNBins() / ig.get(0):ngroups.length;
        H2F hrebinY = new H2F(h.getName(), h.getTitle(), h.getXAxis().getNBins(), h.getXAxis().min(), h.getXAxis().max(),nbinsY, 1, nbinsY+1); 
        hrebinY.getAttributes().setTitleX(h.getTitleX());
        for (int ibx = 0; ibx < h.getXAxis().getNBins(); ibx++) {
           int ngp=0;
           for (int iby = 0; iby < nbinsY; iby++) {
              double height = 0.0;
              int ng = (ngroups.length==1)?ig.get(0):ig.get(iby);
              for (int igroup = 0; igroup < ng; igroup++) height += h.getBinContent(ibx, ngp + igroup);
              hrebinY.setBinContent(ibx, iby, height);
              ngp+=ng;
           }
        }
        return hrebinY;
    }
    
    
    //Performs a vertical merge of a list of H2F
    public H2F CombineH2F(H2F...hlist) {    	
    	int   xbins=0,ybins=0;
    	double xmin=0,xmax=0;
    	for (H2F h: hlist) {xbins=h.getXAxis().getNBins(); ybins+=h.getYAxis().getNBins();xmin=h.getXAxis().min();xmax=h.getXAxis().max();}
    	H2F hnew = new H2F("dum","dum",xbins,xmin,xmax,ybins,1.,ybins+1); 
        int nj=-1;
    	for (H2F h: hlist) {
    		for (int j=0; j<h.getYAxis().getNBins(); j++) {
    			nj++;
    			for (int i=0; i<xbins; i++) hnew.setBinContent(i, nj, h.getBinContent(i, j));
    		}
    	}    	
    	return hnew;
    }

    public FitData fitEngine(H1F h,int ff, double pmin, double pmax, double fmin, double fmax) {
       FitData fd = new FitData(h.getGraph()); 
       fd.setInt((int)h.getIntegral());
       fd.setHist(h);
       fd.graph.getAttributes().setTitleX(h.getTitleX()); 
       fd.hist.getAttributes().setTitleX(h.getTitleX()); 
       fd.initFit(ff,pmin,pmax,fmin,fmax); 
       fd.fitGraph("",fitEnable,fitVerbose); 
       return fd;
    }
        
    public HashMap<Integer,ArrayList<Integer>> mapByIndex(DataBank bank) {
        HashMap<Integer,ArrayList<Integer>> map=new HashMap<Integer,ArrayList<Integer>>();
        for (int ii=0; ii<bank.rows(); ii++) {
            final int index = bank.getInt("pindex", ii);
            if (!map.containsKey(index)) map.put(index,new ArrayList<Integer>());
            map.get(index).add(ii);
        }
        return map;
    }
 
    public void readDataGroup(int run, TDirectory dir) {
        String folder = getDetectorName() + "/";
        System.out.println("Reading from: " + folder);
        DataGroup sum = getDetectorSummary();
        IndexGenerator ig = new IndexGenerator();
        if (sum!=null) {
        int nrows = sum.getRows();
        int ncols = sum.getColumns();
        int nds   = nrows*ncols;
        DataGroup newSum = new DataGroup(ncols,nrows);
        for(int i = 0; i < nds; i++){
            List<IDataSet> dsList = sum.getData(i);
            for(IDataSet ds : dsList){
                System.out.println("\t --> " + ds.getName());
                newSum.addDataSet(dir.getObject(folder, ds.getName()),i);
            }
        }            
        setDetectorSummary(newSum);
        
        }
        
        Map<Long, DataGroup> map = this.getDataGroup().getMap();
        for( Map.Entry<Long, DataGroup> entry : map.entrySet()) {
            Long key = entry.getKey();
            if (run==ig.getIndex(key, 3)) {
            DataGroup group = entry.getValue();
            int nrows = group.getRows();
            int ncols = group.getColumns();
            int nds   = nrows*ncols;
            DataGroup newGroup = new DataGroup(ncols,nrows);
            for(int i = 0; i < nds; i++){
                List<IDataSet> dsList = group.getData(i);
                for(IDataSet ds : dsList){
                    System.out.println("\t --> " + ds.getName());
                    newGroup.addDataSet(dir.getObject(folder, ds.getName()),i);    
                }
            }
            map.replace(key, newGroup);
            }
        }
        this.analyze();
        this.plotHistos(getRunNumber());
    }
    
    public void writeDataGroup(TDirectory dir) {
        String folder = "/" + getDetectorName();
        dir.mkdir(folder);
        dir.cd(folder);
        DataGroup sum = getDetectorSummary();
        if (sum!=null) {
        int nrows = sum.getRows();
        int ncols = sum.getColumns();
        int nds   = nrows*ncols;
        for(int i = 0; i < nds; i++){
            List<IDataSet> dsList = sum.getData(i);
            for(IDataSet ds : dsList){
                System.out.println("\t --> " + ds.getName());
                dir.addDataSet(ds);
            }
        }      
        }
        
        Map<Long, DataGroup> map = getDataGroup().getMap();
        for( Map.Entry<Long, DataGroup> entry : map.entrySet()) {
            DataGroup group = entry.getValue();
            int nrows = group.getRows();
            int ncols = group.getColumns();
            int nds   = nrows*ncols;
            for(int i = 0; i < nds; i++){
                List<IDataSet> dsList = group.getData(i);
                for(IDataSet ds : dsList){
                    System.out.println("\t --> " + ds.getName());
                    dir.addDataSet(ds);
                }
            }
        }
    }
    
    public void drawGroup(EmbeddedCanvas c, DataGroup group) {
        int nrows = group.getRows();
        int ncols = group.getColumns();
        c.divide(ncols, nrows); 	    
        int nds = nrows * ncols;
        for (int i = 0; i < nds; i++) {
            List<IDataSet> dsList = group.getData(i);
            //System.out.println(" pad = " + i + " size = " + dsList.size());
            c.cd(i);  String opt = " ";
            c.getPad().getAxisZ().setLog(getLogZ());
            if(!doAutoRange) c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax);
            if( doAutoRange) c.getPad().getAxisZ().setAutoScale(true);
            for (IDataSet ds : dsList) {
//                System.out.println("\t --> " + ds.getName());
            	   c.draw(ds,opt); opt="same";
            }
        } 	       	
    } 
    
    /**
     * @param fromBank the bank containing the index variable
     * @param idxVarName the name of the index variable
     * @return map with keys being the index in toBank and values the indices in fromBank
     */
     public static Map<Integer,List<Integer>> loadMapByIndex( 
             HipoDataBank fromBank,
             String idxVarName) {
         Map<Integer,List<Integer>> map=new HashMap<Integer,List<Integer>>();
         if (fromBank!=null) {
             for (int iFrom=0; iFrom<fromBank.rows(); iFrom++) {
                 final int iTo = fromBank.getInt(idxVarName,iFrom);
                 if (!map.containsKey(iTo)) map.put(iTo,new ArrayList<Integer>()); 
                 map.get(iTo).add(iFrom);
             }
         }
         return map;
     }
    
}
