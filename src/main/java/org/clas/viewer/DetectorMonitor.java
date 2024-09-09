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
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Timer;
import java.util.TimerTask;
import java.util.TreeMap;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileSystemView;

import org.clas.tools.FitData;
import org.clas.tools.ParallelSliceFitter;
import org.clas.tools.FTHashCollection;
import org.clas.tools.DataGroupManager;
import org.clas.tools.TimeLine;
import org.jlab.detector.base.DetectorOccupancy;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.DataVector;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.Axis;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.RangeSlider;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.eb.EBCCDBConstants;

import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

import org.clas.tools.EmbeddedCanvasTabbed;
import org.clas.tools.EngineControl;

public class DetectorMonitor implements ActionListener {    
    
    private final String                detectorName;
    private ArrayList<String>       detectorTabNames = new ArrayList();
    private IndexedList<DataGroup>      detectorData = new IndexedList<DataGroup>(4);
    public  DataGroupManager                     dgm = new DataGroupManager();
    public  List<Integer>                    runlist = new ArrayList<Integer>();
    public  int                       runIndexSlider = 0;
    public  JSlider                        runslider = null;
    public  JSlider                          yslider = null;
    private DataGroup                detectorSummary = null;
    private DetectorOccupancy      detectorOccupancy = new DetectorOccupancy();
    private JPanel                     detectorPanel = null;
    private JPanel                       actionPanel = null;
    private JPanel                    controlsPanel0 = null;
    private JPanel                    controlsPanel1 = null;
    private JPanel                    controlsPanel2 = null;
    private JPanel                    controlsPanel3 = null;
    private JPanel                     runIndexPanel = null;
    private JPanel                     ECEnginePanel = null;
    private EmbeddedCanvasTabbed      detectorCanvas = null;
    private DetectorPane2D              detectorView = null;
    private ButtonGroup                          bT0 = null;
    private ButtonGroup                          bG0 = null;
    private ButtonGroup                          bG1 = null;
    private ButtonGroup                          bG2 = null;
    private ButtonGroup                          bG3 = null;
    private ButtonGroup                          bG4 = null;
    private ButtonGroup                          bG5 = null;
    private ButtonGroup                          bRO = null;
    private int                       numberOfEvents = 0;
    private int                          totalEvents = 0;
    public  Boolean                    sectorButtons = false;
    private Boolean                      YsliderPane = false;
    private Boolean                      ZsliderPane = false;
    private Boolean                     ECEnginePane = false;
    private int                     detectorActivePC = 1;
    private int                 detectorActiveSector = 1;
    private int                   detectorActiveView = 0;
    private int                  detectorActiveLayer = 0;
    private int                    detectorActivePID = 0;
    private int                    detectorActive123 = 1;
    private int                   detectorActiveSCAL = 0;
    private int                   detectorActiveRDIF = 0;
    private Boolean                     detectorLogY = false;
    private Boolean                     detectorLogZ = true;
    private Boolean                             isTB = false;
    private Boolean                            usePC = false;
    private Boolean                           usePID = false;
    private Boolean                           useUVW = false;
    private Boolean                           use123 = false;
    private Boolean                        useSCALER = false;
    private Boolean                           useCAL = false;
    private Boolean                           useSEC = false;
    private Boolean                          useRDIF = false;
    private Boolean                        useEBCCDB = false;
    private Boolean                        stopBlink = true;
    private Boolean                          verbose = false;
    Timer                                      timer = null;
   
    public JRadioButton bEL,bPI,bPH,bP,bC,bT,b0,b1,bS0,bS1,bS2,bS3,bS4,bS5,bS6,bS7,bpcal,becin,becou,bu,bv,bw;
    public  JCheckBox arBtn,arBtn2;
    
    public int        bitsec = 0;
    public long      trigger = 0;
    public int  triggerPhase = 0;
    public int        trigFD = 0;
    public int        trigCD = 0;
    public int         TRpid = 11;
    public int         MCpid = 11;
    
    public boolean   testTrigger = false;
    public boolean TriggerBeam[] = new boolean[32];
    public int       TriggerMask = 0;
    
    private int   runNumber = 0;
    public  int    runIndex = 0;
    public  int      yIndex = 0;
    public  int   maxevents = 0;
    public  int   MaxEvents = 10002;

    private int eventNumber = 0;
    private int     viewRun = 0;  
    
    public int eventResetTime_current[]=new int[19];
    public int eventResetTime_default[]=new int[19];    
    
    public double lMin=0, lMax=500, zMin=0.1, zMax=1.0, zMinLab, zMaxLab;
    public int izMaxLab;
    public double slideMin=1, slideMax=500;
    public Boolean doAutoRange = false;
    public Boolean dNorm = false;
    
    public String[] layer = new String[]{"pcal","ecin","ecou"};
    public String[]  view = new String[]{"u","v","w"};    
    
    public String variation = "default";
    public String      geom = "2.5";
    public String    config = "phot";  //When re-running ECEngine from cooked data this should always be "phot"
	
	public Boolean  isHipo3Event = true;
    
    String   ltcc[] = {"L","R"};
    String   htcc[] = {"L","R"};
    String   ftof[] = {"PANEL1A_L","PANEL1A_R","PANEL1B_L","PANEL1B_R","PANEL2_L","PANEL2_R"};
    String   ctof[] = {"U","D"};
    String    cnd[] = {"Inner","Middle","Outer"};
    String   band[] = {"1L","1R","2L","2R","3L","3R","4L","4R","5L","5R"};
    String   ecal[] = {"U","V","W","UI","VI","WI","UO","VO","WO"};
    int     nltcc[] = {18,18};
    int     nhtcc[] = {4,4};
    int     nftof[] = {23,23,62,62,5,5};
    int     nctof[] = {48,48};
    int      ncnd[] = {2,2,2};
    int     nband[] = {24,24,24,24,24,24,24,24,20,20};
    int     necal[] = {68,62,62,36,36,36,36,36,36};  
    
    public TreeMap<String,String[]> layMap = new TreeMap<String,String[]>();
    public TreeMap<String,int[]>   nlayMap = new TreeMap<String,int[]>();  
    
    private int[] npmt = {68,62,62,36,36,36,36,36,36};    
    
    public String  outPath = "/home/lcsmith/CLAS12ANA/";
    public String  pawPath = outPath+"paw/";
    public String  jawPath = outPath+"jaw/";
    public String  filPath = null;
    public String  tabPath = null;
    public String   tlPath = null;
    public String  vecPath = null;
    
    public Boolean      isAnalyzeDone = false;
    public Boolean           autoSave = false;
    public Boolean          dropBanks = false;
    public Boolean        dropSummary = false; 
    public Boolean         dumpGraphs = false; 
    public Boolean          dumpFiles = false; 
    public Boolean        defaultGain = false; 
    public Boolean          dropEsect = false;
    public Boolean          onlyEsect = false;
    public Boolean           FTOFveto = false;
    public Boolean           fiduCuts = false;
    public Boolean         cfitEnable = false; 
    public Boolean         dfitEnable = false; 
    public Boolean        gdfitEnable = false; 
    public Boolean         sfitEnable = false; 
    public Boolean        trfitEnable = false; 
    public Boolean         fitVerbose = false; 
    public Boolean      isEngineReady = false;
    public Boolean isTimeLineFitsDone = false;
    public Boolean        histosExist = false;
    public Boolean          dgmActive = false;
    public Boolean          useATDATA = false;
    public Boolean             useGPP = false;
    public Boolean        StatEachRun = false;
    public Boolean            normPix = false;
    public Boolean            normAtt = false;
    public Boolean             normCz = false;
    public Boolean             SFcorr = false;
    public Boolean             TWcorr = false;
    public Boolean              HiRes = false;
    public Boolean        dbgAnalyzer = false;
    public Boolean      analyzeHistos = false;
    
    public IndexedList<FitData>            Fits = new IndexedList<FitData>(4);
    public IndexedList<GraphErrors>  FitSummary = new IndexedList<GraphErrors>(4);
    public IndexedList<GraphErrors>       glist = new IndexedList<GraphErrors>(1);
    
    public TimeLine                   tl = new TimeLine();
    
    public String                 TLname = null;
    public Boolean                TLflag = null;
    public int                     TLmax = 1000;
    public Map<String,Integer> TimeSlice = new HashMap<String,Integer>();  
    public List<Integer>       BlinkRuns = new ArrayList<Integer>();
    
    public ConstantsManager           cm = new ConstantsManager();
    public ConstantsManager         ebcm = new ConstantsManager();
    public EBCCDBConstants        ebccdb = null;
    public int[]                  detcal = {0,0,0};
    public float                TVOffset = 0;
    public float                logParam = 3f;
    public int           PCTrackingPlane = 9;
    public int           ECTrackingPlane = 0;
    public String                   root = " ";
    public String                 osType = " ";
    int                           ntimer = 0;
    public int                   normrun = 0;
    public int                   normrng = 0;
    
	Map<String,Object[]> map = new HashMap<String,Object[]>();
    Object[] can = {false, false, 0, 0, 0, 0};
    
    public EngineControl eng = new EngineControl();
    
    public FTHashCollection rtt = null;
    
    String[]  ccdbTables = new String[]{
            "/calibration/ec/attenuation", 
            "/calibration/ec/atten", 
            "/calibration/ec/gain", 
            "/calibration/ec/timing",
            "/calibration/ec/ftiming",
            "/calibration/ec/ftime",
            "/calibration/ec/dtime",
            "/calibration/ec/fveff",
            "/calibration/ec/dveff",
            "/calibration/ec/time_jitter",
            "/calibration/ec/fadc_offset",
            "/calibration/ec/fadc_global_offset",
            "/calibration/ec/tdc_global_offset",
            "/calibration/ec/global_gain_shift",
            "/calibration/ec/torus_gain_shift",
            "/calibration/ec/global_time_walk",
            "/calibration/ec/effective_velocity",
            "/calibration/ec/tmf_offset",
            "/calibration/ec/tmf_window",
            "/calibration/ec/fthr",
            "/calibration/ec/deff",
            "/calibration/ec/ftres",
            "/calibration/ec/dtres",
            "/daq/fadc/ec",
            "/daq/tt/ec"
    };    
    
    public DetectorMonitor(String name){
    	initGStyle(14);
        detectorName   = name;
        detectorPanel  = new JPanel();
        detectorCanvas = new EmbeddedCanvasTabbed();
        detectorView   = new DetectorPane2D();
        root           = detectorName+".DetectorMonitor.";
        numberOfEvents = 0;        
        eventResetTime_current[0]=0;
        eventResetTime_current[1]=0;     
        eventResetTime_default[0]=10000000;  
        eventResetTime_default[1]=10000000;       
        for (int i=0; i<2; i++){
            eventResetTime_current[i]=eventResetTime_default[i];
        }
        outPath = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/CLAS12ANA/";
        osType  = System.getProperty("os.name");
        pawPath = outPath+"paw/";
        vecPath = pawPath+detectorName+"/";  
        jawPath = outPath+"jaw/";
        tabPath = outPath+detectorName+"/tables/";
        tlPath  = outPath+detectorName+"/timelines/";
        filPath = outPath+detectorName+"/";
        TimeSlice.put("UVW", 3);
        TimeSlice.put("FADC Slot", 16);
        TimeSlice.put("HV Slot", 24);
        System.out.println(root+"outPath = "+outPath);
        System.out.println(root+"osType  = "+osType);
        System.out.println(root+"ebcm.init"); ebcm.init(EBCCDBConstants.getAllTableNames());
        System.out.println(root+"cm.init");     cm.init(Arrays.asList(ccdbTables));
    }
    
    public void init() {
    	System.out.println(root+"init");
        initPanel();    	
    }
    
    public void initTimeLine(String name) {
    	TLname = name;
    	if(getRunNumber()!=0) {plotHistos(getRunNumber()) ; plotScalers(getRunNumber());}
    }
    
    public void initCCDB(int runNumber) {    	
    	System.out.println(root+"initCCDB("+runNumber+")");
    }
        
    public void initEBCCDB(int runNumber) {
    	if(!useEBCCDB) return;
    	System.out.println(root+"initEBCCDB("+runNumber+")");
    	ebccdb = new EBCCDBConstants(runNumber,ebcm);
    }
    
    public void initEPICS() {
        layMap.put("LTCC",ltcc); nlayMap.put("LTCC", nltcc);
        layMap.put("HTCC",htcc); nlayMap.put("HTCC", nhtcc);
        layMap.put("FTOF",ftof); nlayMap.put("FTOF", nftof);
        layMap.put("CTOF",ctof); nlayMap.put("CTOF", nctof);
        layMap.put("BAND",band); nlayMap.put("BAND", nband);
        layMap.put("CND",cnd);   nlayMap.put("CND",  ncnd);
        layMap.put("ECAL",ecal); nlayMap.put("ECAL", necal);    	
    }
    
    public void getReverseTT(ConstantsManager ccdb, int run, String table) {
        System.out.println("monitor.getReverseTT()"); 
        IndexedTable tt = ccdb.getConstants(run, table);
        rtt = new FTHashCollection<int[]>(4);
        for(int ic=1; ic<74; ic++) {
            for (int sl=3; sl<21; sl++) {
                int chmax=16;
                if (sl==6||sl==16) chmax=128;
                for (int ch=0; ch<chmax; ch++){
                    if (tt.hasEntry(ic,sl,ch)) {
                        int[] dum = {ic,sl,ch}; rtt.add(dum,tt.getIntValue("sector",    ic,sl,ch),
                                                            tt.getIntValue("layer",     ic,sl,ch),
                                                            tt.getIntValue("component", ic,sl,ch),
                                                            tt.getIntValue("order",     ic,sl,ch));
                    };
                }
            }
        }
    } 

    public void createTimeLineHistos() {    	
    }
    
    public void localinit(String variation) {    	
    }
    
    public void localclear() {    	
    }
    
    public void initGStyle(int fontsize) {
        GStyle.getAxisAttributesX().setTitleFontSize(fontsize);
        GStyle.getAxisAttributesX().setLabelFontSize(fontsize);
        GStyle.getAxisAttributesY().setTitleFontSize(fontsize);
        GStyle.getAxisAttributesY().setLabelFontSize(fontsize);
        GStyle.getAxisAttributesZ().setLabelFontSize(fontsize); 
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
        actionPanel   = new JPanel(new FlowLayout());
        ECEnginePanel = new JPanel(new FlowLayout());
        runIndexPanel = new JPanel(new FlowLayout());
        controlsPanel0 = new JPanel(new GridBagLayout());
        this.controlsPanel1 = new JPanel();
        this.controlsPanel1.setBorder(BorderFactory.createTitledBorder("Display"));		
        this.controlsPanel2 = new JPanel();
        this.controlsPanel2.setBorder(BorderFactory.createTitledBorder("Timeline"));        
        controlsPanel1.add(packActionPanel());
        controlsPanel2.add(packRunIndexPanel());
        if(ECEnginePane) {
        	this.controlsPanel3 = new JPanel();
        	this.controlsPanel3.setBorder(BorderFactory.createTitledBorder("ECEngine"));        
        	controlsPanel3.add(packECEnginePanel());
        }
        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL; c.weightx = 0.5;		
        c.gridx=0 ; c.gridy=0 ; controlsPanel0.add(controlsPanel1,c);
        c.gridx=0 ; c.gridy=1 ; controlsPanel0.add(controlsPanel2,c);
        if(ECEnginePane) {c.gridx=0 ; c.gridy=2 ; controlsPanel0.add(controlsPanel3,c);}
        getDetectorPanel().add(getDetectorCanvas(),BorderLayout.CENTER);           
        getDetectorPanel().add(controlsPanel0,BorderLayout.SOUTH); 
//        getDetectorPanel().add(packRunIndexPanel(),BorderLayout.SOUTH); 
    }
    
    public void setDbgAnalyzer(Boolean val) {
    	dbgAnalyzer = val;
    }
    
    public String getConfig(int val) {
    	return eng.getConfig(val);   	
    }
    
    public void setDropBanks(Boolean val) {
    	dropBanks = val;
    	if(dropBanks) eng.engCB.doClick();
    }
    
    public void dropBanks(DataEvent de) {
    	
    	if(!isEngineReady) {isEngineReady = true;}
    	
        if(!isHipo3Event&&de.hasBank("ECAL::clusters")) de.removeBanks("ECAL::hits",
        		                                                       "ECAL::peaks",
        		                                                       "ECAL::clusters",
        		                                                       "ECAL::calib",
        		                                                       "ECAL::moments");
        
//        if(!isHipo3Event&&de.hasBank("REC::Event"))       de.removeBanks("REC::Event");
//        if(!isHipo3Event&&de.hasBank("REC::Particle"))    de.removeBanks("REC::Particle");
//        if(!isHipo3Event&&de.hasBank("REC::Calorimeter")) de.removeBanks("REC::Calorimeter");

        if( isHipo3Event&&de.hasBank("ECAL::clusters")) de.removeBank("ECAL::clusters");
        if( isHipo3Event&&de.hasBank("ECAL::hits"))     de.removeBank("ECAL::hits");
        if( isHipo3Event&&de.hasBank("ECAL::peaks"))    de.removeBank("ECAL::peaks");
        if( isHipo3Event&&de.hasBank("ECAL::calib"))    de.removeBank("ECAL::calib");
        if( isHipo3Event&&de.hasBank("ECAL::moments"))  de.removeBank("ECAL::moments");
        
        if(de.hasBank("ECAL::adc")) eng.processDataEvent(de);     	
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
        setNumberOfEvents(getNumberOfEvents()+1);
        switch (event.getType())  {
        case EVENT_START:      processEvent(event); break;
        case EVENT_SINGLE:    {processEvent(event); plotEvent(event);break;}
        case EVENT_ACCUMULATE: processEvent(event); break;
        case EVENT_STOP:       processEVENT_STOP(event);
	    }
    }
    
    public void processEVENT_STOP(DataEvent de) {
        System.out.println(root+"processEVENT_STOP");    	
    	doSTOPEvent(de);
    	analyze(); 
        plotHistos(getRunNumber()); 
        if(autoSave) saveHistosToFile();
        System.out.println(root+"processEVENT_STOP finished");    	
    }
    
    public void doSTOPEvent(DataEvent de) {  	
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
//    public int     getElecTriggerSector()    {return (int) (isGoodFD() ? Math.log10((getFDTrigger()&0x7e)>>1)/0.301+1:0);} 
    public int     getECALTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>19)/0.301+1:0);}       
    public int     getPCALTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>13)/0.301+1:0);}       
    public int     getHTCCTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>7)/0.301+1:0);} 
    
    public int    getTriggerMask()        {return TriggerMask;}
    public void   setTriggerMask(int bit) {TriggerMask|=(1<<bit);}  
    public void clearTriggerMask(int bit) {TriggerMask&=~(1<<bit);}  
    public boolean testTriggerMask()      {return TriggerMask!=0 ? isTrigMaskSet(TriggerMask):true;}
    public boolean isGoodTrigger(int bit) {return TriggerBeam[bit] ? isTrigBitSet(bit):true;}
   
    public int getElecTriggerSector(Boolean shift) { //shift: true=outbending e- false=inbending e-
    	int[] tb = new int[32];
    	int tbsum=0, ts=0;
    	for (int i = 31; i >= 0; i--) {tb[i] = ((trigger & (1 << i))!=0)?1:0;}
    	for (int i=shift?8:1; i<(shift?14:7); i++) {
    		tbsum+=tb[i]; if(tb[i]>0) ts=shift?i-7:i;
    	}
    	return (tbsum==0||tbsum>1) ? 0:ts;
    }
    public int getElecTriggerSector(int f) { //shift: 0,1,2
    	int shift=1+f*7;
    	int[] tb = new int[32];
    	int tbsum=0, ts=0;
    	for (int i = 31; i >= 0; i--) {tb[i] = ((trigger & (1 << i))!=0)?1:0;}
    	for (int i=shift; i<shift+6; i++) {
    		tbsum+=tb[i]; if(tb[i]>0) ts=i-f*7;
    	}
    	return (tbsum==0||tbsum>1) ? 0:ts;
    }
       
    public EmbeddedCanvasTabbed getDetectorCanvas() {
         return detectorCanvas;
    }
    
    public ArrayList<String> getDetectorTabNames() {
        return detectorTabNames;
    }
    
    public IndexedList<DataGroup>  getDataGroup(){
    	if(dgmActive) return dgm.detectorData;
        return detectorData;
    }
    
    public void fillTimeLineGraph() {
    	
    }
    
    public void plotTimeLine(String tab) {
    	
    }
     
    public int getDet(int layer) {
        int[] il = {0,0,0,1,1,1,2,2,2}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
    }
    
    public int getLay(int layer) {
        int[] il = {1,2,3,1,2,3,1,2,3}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
        return il[layer-1];
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
    
    public void useEBCCDB(boolean flag) {
    	useEBCCDB = flag;
    }
    
    public void use123Buttons(boolean flag) {
    	use123 = flag;
    }
    
    public void useSCALERButtons(boolean flag) {
    	useSCALER = flag;
    }
    
    public void useSECButtons(boolean flag) {
    	useSEC = flag;
    }
    
    public void useCALUVWSECButtons(boolean flag) {
    	useCAL = flag;    	
    	useUVW = flag;
    	useSEC = flag;
    }
    
    public void useUVWButtons(boolean flag) {
    	useUVW = flag;
    }
    
    public void useRDIFButtons(boolean flag) {
    	useRDIF = flag;
    } 
    
    public void useCALButtons(boolean flag) {
    	useCAL = flag;
    }
    
    public void useZSliderPane(boolean flag) {
	    ZsliderPane = flag;
    } 
    
    public void useYSliderPane(boolean flag) {
	    YsliderPane = flag;
    } 
    
    public void useECEnginePane(boolean flag) {
    	ECEnginePane = flag;
    }
    
    public void useSliderPane(boolean flag) {
	     useZSliderPane(flag);
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
    
    public int getActiveLayer() { //PCAL,ECIN,ECOU
	    return detectorActiveLayer;
    }
    
    public int getActiveView() { //UVW
	    return detectorActiveView;
    }
    
    public int getActive123() {
	    return detectorActive123;
    }
    
    public int getActiveSCAL() {
	    return detectorActiveSCAL;
    }  
    
    public int getActiveRDIF() {
    	return detectorActiveRDIF;
    }
    
    public int getNumberOfEvents() {
        return numberOfEvents;
    }
    
    public int getTotalEvents() {
    	return totalEvents;
    }
    
    public void setLogY(boolean flag) {
	    detectorLogY = flag;
	    dgm.detectorLogY = flag;
	    if(histosExist) {plotHistos(getRunNumber()); plotScalers(getRunNumber());}
    }   
    
    public void setLogZ(boolean flag) {
	    detectorLogZ = flag;
	    dgm.detectorLogZ = flag;
	    if(histosExist) {plotHistos(getRunNumber()); plotScalers(getRunNumber());}
    }
    
    public void setNormAtt(boolean flag) {
	    normAtt = flag;
	    if(histosExist) {plotHistos(getRunNumber()); plotScalers(getRunNumber());}
    } 
    
    public Boolean getLogY() {
	    return detectorLogY;
    }
    
    public Boolean getLogZ() {
	    return detectorLogZ;
    }
    
    public JPanel packActionPanel() {
        if (use123||useUVW||useSEC) actionPanel.add(getButtonPane());      
        if    (YsliderPane) actionPanel.add(getYSliderPane());
        if    (ZsliderPane) actionPanel.add(getZSliderPane());
    	return actionPanel;
    }
    
    public JPanel packECEnginePanel() {
    	ECEnginePanel.add(eng.getECEnginePane());
    	return ECEnginePanel;    	
    }
    
    public JPanel packRunIndexPanel() {
        runIndexPanel.add(getRunSliderPane());	
        return runIndexPanel;
    }
    
    public JPanel getButtonPane() {
        JPanel buttonPane = new JPanel();
        
        if(useRDIF) {
        b0 = new JRadioButton("RDIF"); buttonPane.add(b0); b0.setActionCommand("1"); b0.addActionListener(this); 
        b1 = new JRadioButton("MEAN"); buttonPane.add(b1); b1.setActionCommand("0"); b1.addActionListener(this); 
        bRO = new ButtonGroup(); bRO.add(b0); bRO.add(b1); 
        b1.setSelected(true);        	
        }
                
        if(usePC) {
        bP = new JRadioButton("P"); buttonPane.add(bP); bP.setActionCommand("2"); bP.addActionListener(this);
        bC = new JRadioButton("C"); buttonPane.add(bC); bC.setActionCommand("1"); bC.addActionListener(this); 
        bT = new JRadioButton("T"); buttonPane.add(bT); bT.setActionCommand("0"); bT.addActionListener(this);
        bG0 = new ButtonGroup(); bG0.add(bP); bG0.add(bC); bG0.add(bT); 
        bT.setSelected(true); bT.doClick();
        }
        
        if(usePID) {
        bEL = new JRadioButton("e-"); buttonPane.add(bEL); bEL.setActionCommand("1"); bP.addActionListener(this);
        bPI = new JRadioButton("pi"); buttonPane.add(bPI); bPI.setActionCommand("2"); bPI.addActionListener(this); 
        bPH = new JRadioButton("ph"); buttonPane.add(bPH); bPH.setActionCommand("3"); bPH.addActionListener(this); 
        bT0 = new ButtonGroup(); bT0.add(bEL); bT0.add(bPI); bT0.add(bPH);
        bEL.setSelected(true);
        }

        if(useSCALER) {
        bG5 = new ButtonGroup();
        bS0 = new JRadioButton("C"); buttonPane.add(bS0); bS0.setActionCommand("0"); bS0.addActionListener(this);
        bS1 = new JRadioButton("D"); buttonPane.add(bS1); bS1.setActionCommand("2"); bS1.addActionListener(this); 
        bS2 = new JRadioButton("T"); buttonPane.add(bS2); bS2.setActionCommand("4"); bS2.addActionListener(this);
        bG5.add(bS0);bG5.add(bS1);bG5.add(bS2);
        bS0.setSelected(true);        		
        }
       
        if(use123) {
        bG4 = new ButtonGroup();
        bS0 = new JRadioButton("0"); buttonPane.add(bS0); bS0.setActionCommand("0"); bS0.addActionListener(this);
        bS1 = new JRadioButton("1"); buttonPane.add(bS1); bS1.setActionCommand("1"); bS1.addActionListener(this); 
        bS2 = new JRadioButton("2"); buttonPane.add(bS2); bS2.setActionCommand("2"); bS2.addActionListener(this); 
        bS3 = new JRadioButton("3"); buttonPane.add(bS3); bS3.setActionCommand("3"); bS3.addActionListener(this); 
        bS4 = new JRadioButton("4"); buttonPane.add(bS4); bS4.setActionCommand("4"); bS4.addActionListener(this);  
        bS5 = new JRadioButton("5"); buttonPane.add(bS5); bS5.setActionCommand("5"); bS5.addActionListener(this); 
        bS6 = new JRadioButton("6"); buttonPane.add(bS6); bS6.setActionCommand("6"); bS6.addActionListener(this); 
        bS7 = new JRadioButton("7"); buttonPane.add(bS7); bS7.setActionCommand("7"); bS7.addActionListener(this); 
   	    bG4.add(bS0);bG4.add(bS1);bG4.add(bS2);bG4.add(bS3);
   	    bG4.add(bS4);bG4.add(bS5);bG4.add(bS6);bG4.add(bS7);
        bS0.setSelected(true);   bS0.doClick();     	
        }
        
        if(useCAL) {
        bpcal = new JRadioButton("PC");  buttonPane.add(bpcal); bpcal.setActionCommand("0"); bpcal.addActionListener(this);
        becin = new JRadioButton("ECi"); buttonPane.add(becin); becin.setActionCommand("1"); becin.addActionListener(this); 
        becou = new JRadioButton("ECo"); buttonPane.add(becou); becou.setActionCommand("2"); becou.addActionListener(this); 
        bG2 = new ButtonGroup(); bG2.add(bpcal); bG2.add(becin); bG2.add(becou);
        bpcal.setSelected(true);
        }
        
        if(useUVW) {
        bu = new JRadioButton("U"); buttonPane.add(bu); bu.setActionCommand("0"); bu.addActionListener(this);
        bv = new JRadioButton("V"); buttonPane.add(bv); bv.setActionCommand("1"); bv.addActionListener(this); 
        bw = new JRadioButton("W"); buttonPane.add(bw); bw.setActionCommand("2"); bw.addActionListener(this); 
        bG3 = new ButtonGroup(); bG3.add(bu); bG3.add(bv); bG3.add(bw);
        bu.setSelected(true); 
        }
        
        if(useSEC) {
        bS1 = new JRadioButton("S1"); buttonPane.add(bS1); bS1.setActionCommand("1"); bS1.addActionListener(this);
        bS2 = new JRadioButton("S2"); buttonPane.add(bS2); bS2.setActionCommand("2"); bS2.addActionListener(this); 
        bS3 = new JRadioButton("S3"); buttonPane.add(bS3); bS3.setActionCommand("3"); bS3.addActionListener(this); 
        bS4 = new JRadioButton("S4"); buttonPane.add(bS4); bS4.setActionCommand("4"); bS4.addActionListener(this); 
        bS5 = new JRadioButton("S5"); buttonPane.add(bS5); bS5.setActionCommand("5"); bS5.addActionListener(this);  
        bS6 = new JRadioButton("S6"); buttonPane.add(bS6); bS6.setActionCommand("6"); bS6.addActionListener(this); 
	    bG1 = new ButtonGroup(); bG1.add(bS1);bG1.add(bS2);bG1.add(bS3);bG1.add(bS4);bG1.add(bS5);bG1.add(bS6);
        bS2.setSelected(true);  
        }   
        
        return buttonPane;
    } 
    
    
    public JPanel getRunSliderPane() {
        JPanel sliderPane = new JPanel();
        JLabel      label = new JLabel("" + String.format("%d", 0));       
        runslider         = new JSlider(JSlider.HORIZONTAL, 0, TLmax, 1); 
        runslider.setPreferredSize(new Dimension(TLmax,10));

        sliderPane.add(new JLabel("Run Index,Run",JLabel.CENTER));
        sliderPane.add(runslider);
        sliderPane.add(label);       
        runslider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                JSlider slider = (JSlider) e.getSource(); 
                if(runlist.size()!=0) {
                	runIndexSlider = (slider.getValue()>runlist.size()-1)?runlist.size()-1:slider.getValue();
                	int run = runlist.get(runIndexSlider);
                    label.setText(String.valueOf(""+String.format("%d", runIndexSlider)+" "+String.format("%d", run)));
                    plotScalers(run);
                    plotHistos(run);
                } else {
                	runIndexSlider = slider.getValue();
                	label.setText(String.valueOf(""+String.format("%d", runIndexSlider)));
                }
            }
        });  
        
        JButton button1 = new JButton("Blink Run");
        button1.addActionListener(this);
        sliderPane.add(button1);
        
        JButton button2 = new JButton("Norm Run");
        button2.addActionListener(this);
        sliderPane.add(button2);
        
        arBtn2 = new JCheckBox("Normalize");
        arBtn2.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                dNorm = (e.getStateChange() == ItemEvent.SELECTED) ? true:false;                
                plotScalers(getRunNumber());
            }
        });         
        arBtn2.setSelected(false);         
        sliderPane.add(arBtn2); 
        
        return sliderPane;
    }
    
    public JPanel getYSliderPane() {
        JPanel sliderPane = new JPanel();
        JLabel     xLabel = new JLabel("Y-Range:");    	
                  yslider = new JSlider(JSlider.HORIZONTAL, 1, 68, 1); 
        JLabel      label = new JLabel("" + String.format("%d", 1));
        sliderPane.add(xLabel);
        sliderPane.add(yslider);
        sliderPane.add(label);        
        yslider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                JSlider slider = (JSlider) e.getSource(); 
                int yIndexSlider = slider.getValue();               
                label.setText(String.valueOf(""+String.format("%d", yIndexSlider)));
                setYIndex(yIndexSlider);
                plotTimeLine("TimeLine");
            }
        });     

        return sliderPane;        
    }
    
    public JPanel getZSliderPane() {
        JPanel sliderPane = new JPanel();
        JLabel     xLabel = new JLabel("Z-Range:");
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
                dgm.zMin=zMin ; dgm.zMax=zMax;
                zMaxLab = Math.pow(10, zMax/10); zMinLab = Math.pow(2, zMin/10); 
                rangeSliderValue1.setText(String.valueOf("" + String.format("%4.0f", zMinLab)));
                rangeSliderValue2.setText(String.valueOf("" + String.format("%4.0f", zMaxLab)));
                if(doAutoRange) izMaxLab = (int) zMaxLab;
                plotScalers(getRunNumber());
                plotHistos(getRunNumber());
            }
        });  
        arBtn = new JCheckBox("AutoRange");
        arBtn.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                doAutoRange = (e.getStateChange() == ItemEvent.SELECTED) ? true:false;
                dgm.doAutoRange = doAutoRange;
                plotScalers(getRunNumber());
                plotHistos(getRunNumber());
            }
        });         
        arBtn.setSelected(false);         
        sliderPane.add(arBtn);        
        return sliderPane;
    } 

    public void actionPerformed(ActionEvent e) {
    	if(e.getActionCommand().compareTo("Blink Run")==0)  BlinkRunAction();
    	if(e.getActionCommand().compareTo("Norm Run")==0)   NormRunAction();
    	if(bG0!=null) detectorActivePC     = Integer.parseInt(bG0.getSelection().getActionCommand()); 
        if(bT0!=null) detectorActivePID    = Integer.parseInt(bT0.getSelection().getActionCommand()); 
        if(bG1!=null) detectorActiveSector = Integer.parseInt(bG1.getSelection().getActionCommand());
        if(bG2!=null) detectorActiveLayer  = Integer.parseInt(bG2.getSelection().getActionCommand());
        if(bG3!=null) detectorActiveView   = Integer.parseInt(bG3.getSelection().getActionCommand());
        if(bG4!=null) detectorActive123    = Integer.parseInt(bG4.getSelection().getActionCommand());
        if(bG5!=null) detectorActiveSCAL   = Integer.parseInt(bG5.getSelection().getActionCommand());
        if(bRO!=null) detectorActiveRDIF   = Integer.parseInt(bRO.getSelection().getActionCommand());
        plotScalers(getRunNumber());
        plotHistos(getRunNumber());
    } 
    
    public void NormRunAction() { 
    	normrun=runIndexSlider;
    	if(StatEachRun) {
    		int imax=Math.min(getNormrng(),maxevents);
    		for(int i=0; i<imax; i++) {normrng=1; NormRunFunction(); normrun++;}
    		normrun=1; normrng=1;
    	} else {
        	normrng = getNormrng();
        	NormRunFunction();
    	}
    }
    
    public int getNormrng() {    
    	int rng = (int) zMaxLab;
    	return normrun+rng-1>=maxevents ? maxevents-normrun:rng;
    }
        
    public void NormRunFunction() {
    	
    }
    
    public void BlinkRunAction() {
    	BlinkRuns.add(runslider.getValue());
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
    	           runslider.setValue(BlinkRuns.get(ntimer));
    	        }
    		}
    	};
    	timer = new Timer(); timer.schedule(task,1000,500);
    }
    
    public void processEvent(DataEvent event) {
        
    }
    
    public void plotEvent(DataEvent event) {
       
    }
    
    public void plotHistos(int run) {

    }    
    
    public void plotScalers(int run) {

    }
        
    public void saveHistosToFile() {
    	String fileName = "CLAS12Ana_run_" + getRunNumber() + "_" + getDetectorName() + ".hipo";
        TDirectory dir = new TDirectory();
        writeDataGroup(dir);
        dir.writeFile(outPath+fileName);
        System.out.println(root+"saveHistosToFile():" + outPath+fileName); 
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
        System.out.println(root+"resetEventListener():" +  getDetectorName() + " histogram for run "+ getRunNumber());
        createHistos(getRunNumber());
        plotHistos(getRunNumber());
    }
    
    public void setCanvasUpdate(int time) {
        for(int tab=0; tab<detectorTabNames.size(); tab++) {
           getDetectorCanvas().getCanvas(detectorTabNames.get(tab)).initTimer(time);
        }
    }
    
    public void setTVOffset(int time) {
    	TVOffset = time;
    }
    
    void setDetectorCanvas(EmbeddedCanvasTabbed canvas) {
        detectorCanvas = canvas;
    }
    
    public void setDetectorTabNames(String... names) {
        for(String name : names) detectorTabNames.add(name); 
        if(dgmActive) {dgm.setDetectorTabNames(names); setDetectorCanvas(dgm.getDetectorCanvas()); return;}
        setDetectorCanvas(new EmbeddedCanvasTabbed(names));      
    }
 
    public void setDetectorSummary(DataGroup group) {
        detectorSummary = group;
    }
    
    public void setNumberOfEvents(int num) {
       numberOfEvents = num;
    }
    
    public void setTotalEvents(int num) {
        totalEvents = num;
    }    
    
    public void setMaxEvents(int num) {
        MaxEvents = num;
    }
         
    public void setEventNumber(int num) {
       eventNumber = num;
    }

    public void setRunNumber(int num) {
       runNumber = num;
    }
    
    public void setYIndex(int num) {
    	yIndex = num;
    }
    
    public int getEventNumber() {
        return eventNumber;
    }

    public int getRunNumber() {
        return runNumber;
    }
    
    public int getRunIndex() {
    	return runIndexSlider;
    }
    
    public int getYIndex() {
    	return yIndex;
    }
    
    public float getBeamEnergy(int run) { 
    	if (run<=100)   return  7.54626f; //rgk-f18
    	if (run<=2360)  return 10.731f;
    	if (run<=2597)  return  2.2219f;
    	if (run<=3002)  return 10.731f;
    	if (run<=3120)  return  6.423f;
    	if (run<=3818)  return 10.594f;
    	if (run<=3861)  return  6.423f;
    	if (run<=4326)  return 10.594f;
    	if (run<=5699)  return 10.6041f;
    	if (run<=5875)  return  7.54626f;
    	if (run<=6000)  return  6.53536f;
    	if (run<=6399)  return 10.5986f;
    	if (run<=6783)  return 10.1998f;
    	if (run<=11283) return 10.4096f;
    	if (run<=11315) return  4.17179f;
    	if (run<=11391) return 10.2129f;
    	if (run<=11599) return 10.3894f;
    	if (run<=11656) return  2.14418f;
    	if (run<=12282) return 10.3894f;
    	if (run<=12447) return  2.1864f;
    	if (run<=12616) return 10.1967f;
    	if (run<=12716) return 10.3394f;
    	if (run<=12955) return 10.4057f;    
    	if (run<=15490) return  5.986f;
//    	if (run<=15727) return  2.07052f;    	
    	if (run<=15727) return  2.06252f;    	
      	if (run<=16079) return  2.21205f;
       	if (run<=18437) return 10.5473f;
       	if (run<=19131) return 10.5322f;
       	if (run<=19659) return 6.39463f;
       	if (run<=19892) return 8.47757f;
        if (run<=99999) return 10.54726f;
    	return 0.0f;
    }
    
    public int getTorusPolarity(int run) {
    	if (run>=3029&&run<=3065) return getTorusColor("-1.0");
    	if (run>=3072&&run<=3087) return getTorusColor("-0.75");
    	if (run>=3097&&run<=3105) return getTorusColor("+0.75");
    	if (run>=3131&&run<=3293) return getTorusColor("+1.00");
    	if (run>=3304&&run<=3817) return getTorusColor("-1.00");
    	if (run>=3819&&run<=3834) return getTorusColor("+0.75");
    	if (run>=3839&&run<=3987) return getTorusColor("+1.00");
    	if (run>=3995&&run<=4326) return getTorusColor("-1.00");  
    	if (run>=4624&&run<=5419) return getTorusColor("-1.00");  
    	if (run>=5420&&run<=5995) return getTorusColor("+1.00");  
    	if (run>=5996&&run<=6000) return getTorusColor("+0.50"); 
    	if (run>=6001&&run<=6141) return 1;
    	if (run>=6142&&run<=6783) return getTorusColor("-1.00");  
    	if (run>=11285&&run<=11285) return getTorusColor("+1.00");  
    	if (run>=11286)           return getTorusColor("-1.00");  
    	return getTorusColor("0");
    }
    
    public float getSolenoidCurrent(int run) {
    	if (run<100) return 1.00f;    	
    	if (run>=2379&&run<=2432) return -1.00f;
    	if (run>=2433&&run<=2449) return -0.60f;
    	if (run>=2450&&run<=2488) return +0.60f;
    	if (run>=2489&&run<=2552) return +1.00f;
    	if (run>=2556&&run<=2571) return +0.60f;
    	if (run>=2572&&run<=2597) return -0.60f;
    	if (run>=15000)           return -1.0f;
    	return 0.00f;
    }
    
    public float getTorusCurrent(int run) {
    	if (run<100) return 1.00f;
    	if (run>=2379&&run<=2394)   return +1.00f;
     	if (run>=2394&&run<=2500)   return +0.60f;
    	if (run>=2526&&run<=2635)   return -0.60f;
    	if (run>=2636&&run<=3065)   return -1.00f;
    	if (run>=3066&&run<=3095)   return -0.75f;
    	if (run>=3096&&run<=3107)   return +0.75f;
    	if (run>=3129&&run<=3293)   return +1.00f;
    	if (run>=3304&&run<=3817)   return -1.00f;
    	if (run>=3819&&run<=3834)   return +0.75f;
    	if (run>=3839&&run<=3987)   return +1.00f;
    	if (run>=3995&&run<=4326)   return -1.00f;  
    	if (run>=4624&&run<=5419)   return -1.00f;  
    	if (run>=5420&&run<=5995)   return +1.00f;  
    	if (run>=5996&&run<=6000)   return +0.50f;  
    	if (run>=6001&&run<=6141)   return  0f;
    	if (run>=6142&&run<=6783)   return -1.00f;  
    	if (run>=11285&&run<=11285) return +1.00f; 
    	if (run>=15000&&run<=15490) return -1.00f;
    	if (run>=15491&&run<=15732) return +0.50f;  
    	if (run>=15733&&run<=18417) return -1.0f;
    	if (run>=18419&&run<=19131) return +1.0f;
    	if (run>19131)              return +1.0f;
    	return 0.00f;
    }
    
    public int shiftTrigBits(int run) { //primary e- trigger bits were not in 0-7
    	if (run>=16043 && run<=16078) return 2; //RG-C 2 GeV
    	if (run>=18419 && getTorusCurrent(run)>0) return 1; //RG-D outbending
    	return 0;
    }
    
    public String getRunGroup(int run) {
    	if (run>=5674&&run<=6000) return "rgk";
    	if (run>=6132&&run<=6604) return "rgb-s19";
    	if (run>=6604)            return "rga-s19";
    	if (run>=11000)           return "rgb";
    	if (run>=15000)           return "rgm";
    	if (run>=16043)           return "rgc";
    	if (run>=18127)           return "rgd";
    	if (run>=19209)           return "rgk";
    	return "rga";
    }
    
    public int getRGIndex(String val) {
    	switch (val) {
    	case "rgk": return 6;
    	case "rgb": return 3;
    	case "rgb-s19": return 3;
    	case "rga": return 1;
    	case "rga-s19": return 1;
    	case "rgm": return 2;
    	}
    	return 0;    	
    }
    
    public int getTorusColor(String val) {
    	switch (val) {
    	case "-1.00": return 1;  
    	case "-0.75": return 6;  
    	case "+0.50": return 4;  
    	case "+0.75": return 3;  
    	case "+1.00": return 2; 
    	}
    	return 1;
    }

    public int getViewRun() {
        return viewRun;
    }   
     
    public void timerUpdate() { 
    }
    
    public void writeFile(String table, int is1, int is2, int il1, int il2, int iv1, int iv2) {   	
    }
    
    public void writeScript(String table) {   	
    }
    
// CANVAS HELPERS
    
    public EmbeddedCanvas getCanvas(String tabname) {
    	int index = getDetectorTabNames().indexOf(tabname);
    	return getDetectorCanvas().getCanvas(getDetectorTabNames().get(index));
    }
        
// HISTO HELPERS    
    
    public H1F makeH1(String name, String tag, int nx, double x1, double x2, String tit, String titx, int ... color) {
	    H1F h1 = new H1F(name+tag,name+tag,nx,x1,x2);
	    if(tit!="") h1.setTitle(tit);
	    h1.setTitleX(titx);
	    int n=0;
        for (int col : color) {
    	    if(n==0) h1.setLineColor(col); 
    	    if(n==1) h1.setFillColor(col); h1.setOptStat("1000000");
    	    if(n==2 && col==1) h1.setOptStat("1000000");
    	    if(n==2 && col==2) h1.setOptStat("1000100");
    	    n++;
        }
	    return h1;
    }
    
    public H1F makeH1(String name, int nx, double x1, double x2, String tit, String titx, String tity, int ... color) {
	    H1F h1 = new H1F(name,nx,x1,x2);
	    if(tit!="") h1.setTitle(tit);
	    h1.setTitleX(titx);  h1.setTitleY(tity);
	    int n=0;	    
        for (int col : color) {
    	    if(n==0) h1.setLineColor(col); 
    	    if(n==1) h1.setFillColor(col);
    	    if(n==2 && col==1) h1.setOptStat("1000000");
    	    if(n==2 && col==2) h1.setOptStat("1000100");
   	        n++;
        }
	    return h1;
    }    
    
    
    public H1F makeH1(String name, int nx, double x1, double x2, String tit, String titx, int ... color) {
	    H1F h1 = new H1F(name,nx,x1,x2);
	    if(tit!="") h1.setTitle(tit);
	    h1.setTitleX(titx);
	    int n=0;
        for (int col : color) {
    	    if(n==0) h1.setLineColor(col); 
    	    if(n==1) h1.setFillColor(col);
    	    if(n==2 && col==1) h1.setOptStat("1000000");
    	    if(n==2 && col==2) h1.setOptStat("1000100");
    	    n++;
        }
	    return h1;
    }

    public H2F makeH2(String name, String tag, int nx, double x1, double x2, int ny, double y1, double y2, String tit, String titx, String tity) {
	    H2F h2 = new H2F(name+tag,name+tag,nx,x1,x2,ny,y1,y2);
	    if(tit!="") h2.setTitle(tit);
	    h2.setTitleX(titx); h2.setTitleY(tity);
	    return h2;
    }
    
    public H2F makeH2(String name, int nx, double x1, double x2, int ny, double y1, double y2, String tit, String titx, String tity) {
	    H2F h2 = new H2F(name,nx,x1,x2,ny,y1,y2);
	    if(tit!="") h2.setTitle(tit);
	    h2.setTitleX(titx); h2.setTitleY(tity);
	    return h2;
    }
    
    public GraphErrors makeGraph(String name, int ... options) {
    	GraphErrors graph = new GraphErrors();
    	graph.setName(name); 
    	int n=0;
    	for (int opt : options) {
    		if(n==0) graph.setMarkerColor(opt); 
    		if(n==1) graph.setMarkerSize(opt); 
    		if(n==2) graph.setMarkerStyle(opt);
    		n++;
    	}    	
    	return graph;
    }
    
    public H1F subset(H1F hin, float x1, float x2, int col) {
    	H1F hout = hin.histClone("dup"); hout.reset();
    	for (int i=0; i<hin.getDataSize(0); i++) {
    	  if(hin.getDataX(i)>x1 && hin.getDataX(i)<x2) hout.setBinContent(i, hin.getDataY(i));
    	}
    	hout.setFillColor(col); hout.setOptStat("1000000");
    	return hout;
    }

    public H1F projectionX(H2F h2, float ymin, float ymax ) {
        String name = "X Projection";
        Axis xAxis = h2.getXAxis(), yAxis = h2.getYAxis();
        H1F projX = new H1F(name, xAxis.getNBins(), xAxis.min(), xAxis.max());
        for (int x = 0; x < xAxis.getNBins(); x++) {
            double height = 0.0;
            for (int y = yAxis.getBin(ymin); y < yAxis.getBin(ymax); y++) {
                height += h2.getBinContent(x, y);
            }
            projX.setBinContent(x, height);
        }        
        return projX;
    }
    
    public H1F projectionY(H2F h2, float xmin, float xmax ) {
        String name = "Y Projection";
        Axis xAxis = h2.getXAxis(), yAxis = h2.getYAxis();
        H1F projY = new H1F(name, yAxis.getNBins(), yAxis.min(), yAxis.max());
        for (int y = 0; y < yAxis.getNBins(); y++) {
            double height = 0.0;
            for (int x = xAxis.getBin(xmin); x < xAxis.getBin(xmax); x++) {
                height += h2.getBinContent(x, y);
            }
            projY.setBinContent(y, height);
        }
        
        return projY;
    }

    
// GRAPH HELPERS   
        
    public void GraphPlot(GraphErrors graph, EmbeddedCanvas c, int zone, float xmin, float xmax, float ymin, float ymax, int mcol, int msiz, int msty, String xtit, String ytit, String opt) {
    	c.cd(zone); c.getPad(zone).getAxisX().setRange(xmin, xmax); c.getPad(zone).getAxisY().setRange(ymin, ymax); 
    	graph.setMarkerColor(mcol); graph.setMarkerSize(msiz); graph.setMarkerStyle(msty); graph.setLineColor(mcol);
    	if(!(xtit=="")) graph.setTitleX(xtit); 
    	if(!(ytit=="")) graph.setTitleY(ytit); 
    	c.draw(graph,opt);
    }
    
    public GraphErrors getFitSlices(H2F h, String xy, int col, float... mnmx) {
    	
    	double[] dum = new double[6];
    	float min=0, max=0;
    	if(h.getDataBufferSize()==0) return new GraphErrors("dum",dum,dum,dum,dum);
    	
    	ParallelSliceFitter fitter;    	
    	fitter = new ParallelSliceFitter(h); fitter.setBackgroundOrder(0); 
    	if(mnmx.length==2) fitter.setRange(mnmx[0], mnmx[1]);
    	if(xy=="x") fitter.fitSlicesX(); 
    	if(xy=="y") fitter.fitSlicesY(); 
    	
    	GraphErrors graph = new GraphErrors("graph",fitter.getMeanSlices().getVectorX().getArray(),
    			                                    fitter.getMeanSlices().getVectorY().getArray(),
    			                                    new double[h.getDataSize(0)],
    			                                    fitter.getSigmaSlices().getVectorY().getArray());
    	graph.setTitleX(h.getTitleX()); graph.setTitleY(h.getTitleY()); graph.setMarkerColor(col);
    	return graph;
    }
    
    
    public GraphErrors scaleGraph(GraphErrors g, double sca) {
        int np = g.getDataSize(0);
    	double[] x = new double[np]; double[] xe = new double[np]; 
    	double[] y = new double[np]; double[] ye = new double[np];
    	for (int i=0; i<np; i++) {x[i]=g.getDataX(i); y[i]=g.getDataY(i)/sca; xe[i]=g.getDataEX(i); ye[i]=g.getDataEY(i)/sca;}     	
    	return new GraphErrors("SCALED",x,y,xe,ye);
    }  
    
    public GraphErrors sliceToGraph(GraphErrors gin, int il, int iv) {
    	int np = npmt[3*il+iv]; 
    	double[] x = new double[np]; double[] xe = new double[np]; 
    	double[] y = new double[np]; double[] ye = new double[np]; 
    	double[] dx =  gin.getVectorX().getArray(); double[] dy = gin.getVectorY().getArray();
    	for (int i=0; i<np; i++) x[i]=i+1;
    	int n=0; for(double ddy : dy) {y[(int)(dx[n]-1)] = ddy; n++;}
        return new GraphErrors("TMF",x,y,xe,ye);    	
    } 
    
    public DataVector getGoodMean(DataVector v1, DataVector v2, double err) {
        DataVector v3= new DataVector();
        int np = v1.size();
        for (int i=0; i<np; i++) {
        	v3.add(Math.abs(v1.getValue(i)-v2.getValue(i))>err ? v2.getValue(i):v1.getValue(i));
        }
        return v3;
    }
    
    public GraphErrors findGoodMean(GraphErrors g1, GraphErrors g2, double err) {

    	double[] dx1 = g1.getVectorX().getArray();
        double[] dy1 = g1.getVectorY().getArray();
        double[] dy2 = g2.getVectorY().getArray();
        int np=dx1.length;
    	double[] x = new double[np]; double[] xe = new double[np]; 
    	double[] y = new double[np]; double[] ye = new double[np]; 
        for (int i=0; i<np; i++) {
        	x[i]=i+1;
        	y[i]=Math.abs(dy1[i]-dy2[i])>err ? dy2[i]:dy1[i];
        }
    	return new GraphErrors("GM",x,y,xe,ye);
    }
        
    public void dumpGraph(String filename, GraphErrors graph) {
    	PrintWriter writer = null;
		try 
		{
			writer = new PrintWriter(filename);
		} 
		catch (FileNotFoundException e) 
		{
			e.printStackTrace();
		}				
    	for (int i=0; i<graph.getDataSize(0); i++) writer.println(String.format("%1$.3f %2$.6f %3$.8f",graph.getDataX(i),graph.getDataY(i),graph.getDataEY(i))); 
    	writer.close();    	
    }
    
    public GraphErrors getGraph(String filename) {   
        
    	GraphErrors g = new GraphErrors(); 
        try{
            FileReader       file = new FileReader(filename);
            BufferedReader     br = new BufferedReader(file);
            String s; int n = 0 ;
            while ((s = br.readLine()) != null) {              
              String[] col = s.trim().split("\\s+"); 
              float i = Float.parseFloat(col[0]); float j = Float.parseFloat(col[1]);float k = Float.parseFloat(col[2]);
              g.addPoint(i,j,0,k);
              n++;
            }    
            br.close();
            file.close();
         }  
         
         catch(FileNotFoundException ex) {
            ex.printStackTrace();            
         }     
         catch(IOException ex) {
             ex.printStackTrace();
         }
         return g;

    } 
    
    //TIMELINE HELPERS
    
    public void setTLflag(Boolean flag) {
    	this.TLflag = flag;
    	if (getRunNumber()>0) plotHistos(getRunNumber());
    }
    
    public void drawTimeLine(EmbeddedCanvas c, int is, int i, float norm, String title) {
    	
    	GraphErrors g2 = null;
    	
		List<GraphErrors> gglist = getGraph(((H2F)tl.Timeline.getItem(i,0)),((H2F)tl.Timeline.getItem(i,1)),is-1);

    	DataLine line = new DataLine(-0.5,norm,runIndex,norm); line.setLineColor(3); line.setLineWidth(2);
        
		for (int ii=1; ii<gglist.size(); ii++) {    
    		gglist.get(ii).setTitleX("Run Index"); gglist.get(ii).setTitleY(title); 
			c.draw(gglist.get(ii),(ii==1)?" ":"same"); c.draw(line);
		}
		
		g2 = new GraphErrors(); g2.setMarkerSize(5); g2.setMarkerColor(4); g2.setLineColor(2);
		g2.addPoint(runIndexSlider,gglist.get(0).getDataY(runIndexSlider),0,0); c.draw(g2,"same");    	
    }
    
    
    public List<GraphErrors> getGraph(H2F h2a, H2F h2b, int ybin) {
	    int[] col = {1,1,2,5,6,7,8};
	    H1F h1a = h2a.getSlicesY().get(ybin); H1F h1b = h2b.getSlicesY().get(ybin); 
	    glist.clear();
	    GraphErrors g = new GraphErrors() ; g.setLineColor(col[0]); g.setMarkerColor(col[0]); g.setMarkerSize(3); glist.add(g,0);
	    List<GraphErrors> gglist = new ArrayList<GraphErrors>();
	    for (int i=0; i<runlist.size(); i++) {
	    	int it = getTorusPolarity(runlist.get(i)); int im = getRGIndex(getRunGroup(runlist.get(i)));
	    	glist.getItem(0).addPoint(h1a.getDataX(i)-0.5, h1a.getDataY(i),0, h1b.getDataY(i)); 
	    	if (!glist.hasItem(it)) {g = new GraphErrors() ; g.setLineColor(col[it]); g.setMarkerStyle(im); g.setMarkerColor(col[it]); g.setMarkerSize(3); glist.add(g,it);} 
	    	glist.getItem(it).addPoint(h1a.getDataX(i)-0.5, h1a.getDataY(i),0, h1b.getDataY(i));
	    }   
	    for (int i=0; i<6; i++) if(glist.hasItem(i)) gglist.add(glist.getItem(i));	    
	    return gglist;
    }
    
    public void saveTimelines() {
    	
    }
    
    public void saveTimeLine(int i, int il, int ip, String fname, String tag) {
		 TDirectory dir = new TDirectory();
		 GraphErrors g[]= new GraphErrors[6];
		 for (int is=0; is<6; is++) g[is] = new GraphErrors("sec"+(is+1));		 
		 H2F h2a = (H2F) tl.Timeline.getItem(i,0); H2F h2b = (H2F) tl.Timeline.getItem(i,1);
		 for (int ir=0; ir<runlist.size(); ir++) {
			 dir.mkdir("/"+runlist.get(ir)); dir.cd("/"+runlist.get(ir));
			 for (int is=0; is<6; is++) {
				 g[is].addPoint(runlist.get(ir),h2a.getSlicesY().get(is).getDataY(ir),0,h2b.getSlicesY().get(is).getDataY(ir));				  			 
		         FitData fd = tl.fitData.getItem(is+1,il,ip,runlist.get(ir));
		         H1F h1 = fd.getHist(); h1.setTitle("hsec"+(is+1)); h1.setName(tag+" Sector "+(is+1));
		         F1D f1 = new F1D("fit:"+h1.getName()); f1 = (F1D) fd.graph.getFunction(); f1.setName("fit:"+h1.getName());
		         dir.addDataSet(h1); dir.addDataSet(f1);
			 }
			 
		 }
		 dir.mkdir("/timelines");  dir.cd("/timelines");
		 for (int is=0; is<6; is++) dir.addDataSet(g[is]);
    	 System.out.println("Saving timeline to "+tlPath+fname+".hipo");
		 dir.writeFile(tlPath+fname+".hipo");		 
    }
    
    //Generalization of H2F rebinY method for list of non-equal ngroups. 
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
    
    
    //Performs a y-axis concatenation of a list of H2F
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
    
    //FIT HELPERS
    
    public FitData fitEngine(GraphErrors g, int ff, int fmin) {
        FitData fd = new FitData(g);        
        fd.doFit = g.getDataSize(0)==0 ? false:true;
        double fmax = fd.doFit ? g.getDataX(g.getDataSize(0)-1)*1.05 : 0;
        fd.initFit(ff,0,0,fmin,fmax); 
        fd.fitGraph("",cfitEnable,fitVerbose); 
        return fd;
     }
    
    public FitData fitEngine(GraphErrors g, int ff, int fmin, int fmax) {
        FitData fd = new FitData(g);        
        fd.initFit(ff,0,1,fmin,fmax); 
        fd.doFit = g.getDataSize(0)==0 ? false:true; 	
        fd.fitGraph("",cfitEnable,fitVerbose); 
        return fd;
     } 
    
    public FitData fitEngine(GraphErrors g, int ff, double pmin, double pmax, double fmin, double fmax ) {
        FitData fd = new FitData(g);        
        fd.initFit(ff,pmin,pmax,fmin,fmax); 
        fd.doFit = g.getDataSize(0)==0 ? false:true; 	
        fd.fitGraph("",cfitEnable,fitVerbose); 
        return fd;
     }

    public FitData fitEngine(H1F h, int ff, double hmax) {
       FitData fd = new FitData(h.getGraph()); 
       fd.setInt((int)h.getIntegral());
       fd.setHist(h);
       fd.graph.getAttributes().setTitleX(h.getTitleX()); 
       fd.hist.getAttributes().setTitleX(h.getTitleX()); 
       fd.doMeanHist(hmax); //use Mean of blue MIP histo
       fd.initFit(ff,hmax); 
       fd.fitGraph("",cfitEnable,fitVerbose); 
       return fd;
    }
    
    public FitData fitEngine(H1F h, int ff, double pmin, double pmax, double fmin, double fmax) {
       FitData fd = new FitData(h.getGraph()); 
       fd.setInt((int)h.getIntegral());
       String tit = h.getTitle();
       fd.setHist(h);
       fd.graph.getAttributes().setTitleX(h.getTitleX()); 
       fd.hist.getAttributes().setTitleX(h.getTitleX());
       fd.hist.setTitle(tit);
       fd.graph.setTitle(tit);
       fd.initFit(ff,pmin,pmax,fmin,fmax);
       fd.fitGraph("",cfitEnable,fitVerbose); 
       return fd;
    }
    
    public FitData fitEngine(H1F h, double pmin, double pmax, double fmin, double fmax) {
        FitData fd = new FitData(h.getGraph()); 
        fd.setHist(h);
        fd.simpleFit(pmin,pmax,fmin,fmax);
        return fd;
    }
    
    public FitData fitEngine(H1F h, int ff, double pmin, double pmax, double fmin, double fmax, double sig1, double sig2) {
        FitData fd = new FitData(h.getGraph()); 
        fd.setInt((int)h.getIntegral());
        fd.setHist(h);
        fd.graph.getAttributes().setTitleX(h.getTitleX()); 
        fd.hist.getAttributes().setTitleX(h.getTitleX()); 
        fd.setSigmas(sig1,sig2);
        fd.initFit(ff,pmin,pmax,fmin,fmax); 
        fd.fitGraph("",cfitEnable,fitVerbose); 
        return fd;
     }
        
    public HashMap<Integer,ArrayList<Integer>> mapByIndex(DataBank bank, String val) {
        HashMap<Integer,ArrayList<Integer>> map=new HashMap<Integer,ArrayList<Integer>>();
        for (int ii=0; ii<bank.rows(); ii++) {
            final int index = bank.getInt(val, ii);
            if (!map.containsKey(index)) map.put(index,new ArrayList<Integer>());
            map.get(index).add(ii);
        }
        return map;
    }
    
    //DATAGROUP HELPERS
    
    //Used to associate canvas configuration with histogram name
	public void cc(String name, boolean linlogy, boolean linlogz, float ymin, float ymax, float zmin, float zmax) {
		Object[] obj = {linlogy,linlogz,ymin,ymax,zmin,zmax};
		map.put(name,obj);  
	}    
    
    //Implementation of canvas configuration map
	public Boolean config(String name, EmbeddedCanvas canvas) {
		if(!map.containsKey(name)) return false;
		can = map.get(name); 
		canvas.getPad().getAxisY().setLog((Boolean)can[0]);
		canvas.getPad().getAxisZ().setLog((Boolean)can[1]);
		float ymin=(float)can[2], ymax=(float)can[3], zmin=(float)can[4], zmax=(float)can[5]; 
		if(ymin!=ymax) canvas.getPad().getAxisY().setRange(ymin,ymax);
		if(zmin!=zmax) canvas.getPad().getAxisZ().setRange(zmin,zmax);		
		return true;
	}    
 
    public void readDataGroup(int run, TDirectory dir) {
        String folder = getDetectorName() + "/";
        System.out.println("Reading from: " + folder);
        IndexGenerator ig = new IndexGenerator();

        if (getDetectorSummary()!=null) { //seems to be always true?        
        	DataGroup sum = getDetectorSummary();
        	int nrows = sum.getRows();
        	int ncols = sum.getColumns();
        	int nds   = nrows*ncols;
        	DataGroup newSum = new DataGroup(ncols,nrows);
        	for(int i = 0; i < nds; i++){
        		List<IDataSet> dsList = sum.getData(i);
        		for(IDataSet ds : dsList){
        			if(verbose)System.out.println("\t --> " + ds.getName());
        			newSum.addDataSet(dir.getObject(folder, ds.getName()),i);
        		}
        	}            
        	setDetectorSummary(newSum);
        }
        
        Map<Long, DataGroup> map = this.getDataGroup().getMap(); //All defined histos are here
        
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
                    if(verbose)System.out.println("\t --> " + ds.getName());                 
                    if(dir.getObject(folder, ds.getName())!=null) newGroup.addDataSet(dir.getObject(folder, ds.getName()),i);    
                    if(dir.getObject(folder, ds.getName())==null) newGroup.addDataSet(ds,i);    
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
                if(verbose)System.out.println("\t --> " + ds.getName());
                dir.addDataSet(ds);
            }
        }      
        }
        
        Map<Long, DataGroup> map = getDataGroup().getMap();  //All defined histos are here
        
        for( Map.Entry<Long, DataGroup> entry : map.entrySet()) {
            DataGroup group = entry.getValue(); 
            int nrows = group.getRows();
            int ncols = group.getColumns();
            int nds   = nrows*ncols;
            for(int i = 0; i < nds; i++){
                List<IDataSet> dsList = group.getData(i);
                for(IDataSet ds : dsList){
                    if(verbose)System.out.println("\t --> " + ds.getName());
                    dir.addDataSet(ds);
                }
            }
        }
    }
    
	public void drawGroup(String tabname, int is, int st, int run) {
    	int index = getDetectorTabNames().indexOf(tabname);
    	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),getDataGroup().getItem(is,st,index,run));   			
	}
    
    public void drawGroup(EmbeddedCanvas c, DataGroup group) {
    	if(group==null) return;
        int nrows = group.getRows();
        int ncols = group.getColumns();
        c.clear(); c.divide(ncols, nrows); 	    
        int nds = nrows * ncols;
        for (int i = 0; i < nds; i++) {
            boolean maxtest = false;
            List<IDataSet> dsList = group.getData(i);
            c.cd(i);  String opt = " "; 
            if(dsList.size()>2) maxtest = true; //prevents AutoScale from triggering on empty or low yield histogram
            if(dsList.size()>0) {
            	IDataSet ds0 = dsList.get(0);            	
            	if(!config(ds0.getName(),c)) { //override canvas auto configuration if ds0 is tagged
            		if (ds0 instanceof H1F) {
            			c.getPad().getAxisY().setAutoScale(true);
            			c.getPad().getAxisY().setLog(getLogY()); 
            		}
            		if (ds0 instanceof H2F) {
            			c.getPad().getAxisZ().setLog(getLogZ());
            			if(!doAutoRange)  c.getPad().getAxisZ().setRange(0.1*zMin, 200*zMax); 
            			if( doAutoRange) {c.getPad().getAxisX().setAutoScale(true);
		                                  c.getPad().getAxisY().setAutoScale(true);
		                                  c.getPad().getAxisZ().setAutoScale(true);}
            		}
            		if (ds0 instanceof GraphErrors) {
            			if( doAutoRange) {c.getPad().getAxisX().setAutoScale(true);
                                          c.getPad().getAxisY().setAutoScale(true);
                                          c.getPad().getAxisZ().setAutoScale(true);}            			
            		}
            	}
            	for (IDataSet ds : dsList) if(maxtest?ds.getMax()>0.5:true) {c.draw(ds,opt); opt="same";}
            }
        } 	       	
    }
    
    /**
     * @param fromBank the bank containing the index variable
     * @param idxVarName the name of the index variable
     * @return map with keys being the index in toBank and values the indices in fromBank
     */
     public static Map<Integer,List<Integer>> loadMapByIndex( 
             DataBank fromBank,
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
