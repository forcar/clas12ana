package org.clas.viewer;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSplitPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.jlab.detector.base.DetectorOccupancy;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.RangeSlider;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.task.IDataEventListener;
import org.jlab.utils.groups.IndexedList;


public class DetectorMonitor implements IDataEventListener, ActionListener {    
    
    private final String           detectorName;
    private ArrayList<String>      detectorTabNames  = new ArrayList();
    private IndexedList<DataGroup> detectorData      = new IndexedList<DataGroup>(4);
    private DataGroup              detectorSummary   = null;
    private DetectorOccupancy      detectorOccupancy = new DetectorOccupancy();
    private JPanel                 detectorPanel     = null;
    private JPanel                   actionPanel     = null;
    private EmbeddedCanvasTabbed   detectorCanvas    = null;
    private DetectorPane2D         detectorView      = null;
    private ButtonGroup                          bG1 = null;
    private ButtonGroup                          bG2 = null;
    private ButtonGroup                          bG3 = null;
    private int                       numberOfEvents = 0;
    private Boolean                    sectorButtons = false;
    private Boolean                       sliderPane = false;
    private int                 detectorActiveSector = 1;
    private int                   detectorActiveView = 0;
    private int                  detectorActiveLayer = 0;
    private Boolean                     detectorLogZ = true;
    private Boolean                             isTB = false;
    
    private JRadioButton bS1,bS2,bS3,bS4,bS5,bS6,bpcal,becin,becou,bu,bv,bw;
    private JCheckBox tbBtn;
    
    public int        bitsec = 0;
    public long      trigger = 0;
    public long triggerPhase = 0;
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
    
    public double zMin=0.1, zMax=1.0, zMinLab, zMaxLab;
    public double slideMin=1, slideMax=500;
    public  Boolean doSlider = false;
    
    public String[] layer = new String[]{"pcal","ecin","ecou"};
    public String[]  view = new String[]{"u","v","w"};    
    
    public DetectorMonitor(String name){
        GStyle.getAxisAttributesX().setTitleFontSize(14);
        GStyle.getAxisAttributesX().setLabelFontSize(14);
        GStyle.getAxisAttributesY().setTitleFontSize(14);
        GStyle.getAxisAttributesY().setLabelFontSize(14);
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
//        GStyle.getAxisAttributesZ().setLog(true);
        
        this.detectorName = name;
        this.detectorPanel  = new JPanel();
        this.detectorCanvas = new EmbeddedCanvasTabbed();
        this.detectorView   = new DetectorPane2D();
        this.numberOfEvents = 0;
        
        eventResetTime_current[0]=0;
        eventResetTime_current[1]=0;
     
        eventResetTime_default[0]=10000000;  
        eventResetTime_default[1]=10000000;  
     
        for (int i=0; i<2; i++){
            eventResetTime_current[i]=eventResetTime_default[i];
        }
    }
    
    public void init() {
        getDetectorPanel().setLayout(new BorderLayout());
        actionPanel = new JPanel();
        actionPanel.setLayout(new FlowLayout());
        getDetectorPanel().add(getDetectorCanvas(),BorderLayout.CENTER);           
        getDetectorPanel().add(packActionPanel(),BorderLayout.PAGE_END); 
        createHistos(0);
        plotHistos(0); 
        if (sectorButtons) {bS2.doClick();}
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
    
    @Override
    public void dataEventAction(DataEvent event) {
        if (!testTriggerMask()) return;
        this.setNumberOfEvents(this.getNumberOfEvents()+1);
        if (event.getType() == DataEventType.EVENT_START) {
//            resetEventListener();
            processEvent(event);
	} else if (event.getType() == DataEventType.EVENT_SINGLE) {
            processEvent(event);
            plotEvent(event);
	} else if (event.getType() == DataEventType.EVENT_ACCUMULATE) {
            processEvent(event);
	} else if (event.getType() == DataEventType.EVENT_STOP) {
            analyze();
	}
    }

    public void drawDetector() {
    
    }
    
    public void setTriggerPhase(long phase) {
    	   this.triggerPhase = phase;
    }
    
    public long getTriggerPhase() {
    	    return this.triggerPhase;
    }
    
    public void setTriggerWord(long trig) {
    	   this.trigger = trig;
    }
    
    public void setTestTrigger(boolean test) {
    	   this.testTrigger = test;
    }
    
    public int     getFDTrigger()            {return (int)(this.trigger)&0x000000000ffffffff;}
    public int     getCDTrigger()            {return (int)(this.trigger>>32)&0x00000000ffffffff;}
    public boolean isGoodFD()                {return  getFDTrigger()>0;}    
    public boolean isTrigBitSet(int bit)     {int mask=0; mask |= 1<<bit; return isTrigMaskSet(mask);}
    public boolean isTrigMaskSet(int mask)   {return (getFDTrigger()&mask)!=0;}
    public boolean isGoodECALTrigger(int is) {return (testTrigger)? is==getECALTriggerSector():true;}    
    public int           getElecTrigger()    {return getFDTrigger()&0x1;}
    public int     getElecTriggerSector()    {return (int) (isGoodFD() ? Math.log10((getFDTrigger()&0x7e)>>1)/0.301+1:0);} 
    public int     getECALTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>19)/0.301+1:0);}       
    public int     getPCALTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>13)/0.301+1:0);}       
    public int     getHTCCTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>7)/0.301+1:0);} 
    
    public int    getTriggerMask()        {return this.TriggerMask;}
    public void   setTriggerMask(int bit) {this.TriggerMask|=(1<<bit);}  
    public void clearTriggerMask(int bit) {this.TriggerMask&=~(1<<bit);}  
    public boolean testTriggerMask()      {return this.TriggerMask!=0 ? isTrigMaskSet(this.TriggerMask):true;}
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
    	    this.sectorButtons = flag;
    }
    
    public void useSliderPane(boolean flag) {
    	    this.sliderPane = flag;
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
	    this.detectorLogZ = flag;
    }
    
    public Boolean getLogZ() {
	    return this.detectorLogZ;
    }
    
    public JPanel packActionPanel() {
        if (sectorButtons) actionPanel.add(getButtonPane()); 
        if    (sliderPane) actionPanel.add(getSliderPane());
    	    return actionPanel;
    }
    
    public JPanel getButtonPane() {
        JPanel buttonPane = new JPanel();
        bS1 = new JRadioButton("Sector 1"); buttonPane.add(bS1); bS1.setActionCommand("1"); bS1.addActionListener(this);
        bS2 = new JRadioButton("Sector 2"); buttonPane.add(bS2); bS2.setActionCommand("2"); bS2.addActionListener(this); 
        bS3 = new JRadioButton("Sector 3"); buttonPane.add(bS3); bS3.setActionCommand("3"); bS3.addActionListener(this); 
        bS4 = new JRadioButton("Sector 4"); buttonPane.add(bS4); bS4.setActionCommand("4"); bS4.addActionListener(this); 
        bS5 = new JRadioButton("Sector 5"); buttonPane.add(bS5); bS5.setActionCommand("5"); bS5.addActionListener(this);  
        bS6 = new JRadioButton("Sector 6"); buttonPane.add(bS6); bS6.setActionCommand("6"); bS6.addActionListener(this); 
	    bG1 = new ButtonGroup(); bG1.add(bS1);bG1.add(bS2);bG1.add(bS3);bG1.add(bS4);bG1.add(bS5);bG1.add(bS6);
        bS2.setSelected(true);        
        bpcal = new JRadioButton("PCAL"); buttonPane.add(bpcal); bpcal.setActionCommand("0"); bpcal.addActionListener(this);
        becin = new JRadioButton("ECin"); buttonPane.add(becin); becin.setActionCommand("1"); becin.addActionListener(this); 
        becou = new JRadioButton("ECou"); buttonPane.add(becou); becou.setActionCommand("2"); becou.addActionListener(this); 
        bG2 = new ButtonGroup(); bG2.add(bpcal); bG2.add(becin); bG2.add(becou);
        bpcal.setSelected(true);
        bu = new JRadioButton("U"); buttonPane.add(bu); bu.setActionCommand("0"); bu.addActionListener(this);
        bv = new JRadioButton("V"); buttonPane.add(bv); bv.setActionCommand("1"); bv.addActionListener(this); 
        bw = new JRadioButton("W"); buttonPane.add(bw); bw.setActionCommand("2"); bw.addActionListener(this); 
        bG3 = new ButtonGroup(); bG3.add(bu); bG3.add(bv); bG3.add(bw);
        bu.setSelected(true);                
        return buttonPane;
    } 
    
    public JPanel getSliderPane() {
    	    JPanel sliderPane = new JPanel();
        JLabel xLabel = new JLabel("Z-Range:");
        RangeSlider slider = new RangeSlider();
        slider.setMinimum((int) slideMin);
        slider.setMaximum((int) slideMax);
        slider.setValue((int)slideMin);
        slider.setUpperValue((int)slideMax);            
        zMin =     slider.getValue();
        zMax = 0.1*slider.getUpperValue();
        zMinLab = Math.pow(2, zMin/10); zMaxLab = Math.pow(10, zMax/10);
        JLabel rangeSliderValue1 = new JLabel("" + String.format("%4.0f", zMinLab));
        JLabel rangeSliderValue2 = new JLabel("" + String.format("%4.0f", zMaxLab));
        sliderPane.add(xLabel);
        sliderPane.add(rangeSliderValue1);
        sliderPane.add(slider);
        sliderPane.add(rangeSliderValue2);           
        slider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                RangeSlider slider = (RangeSlider) e.getSource();
                zMin =     slider.getValue(); 
                zMax = 0.1*slider.getUpperValue();
                zMinLab = Math.pow(2, zMin/10); zMaxLab = Math.pow(10, zMax/10);
                rangeSliderValue1.setText(String.valueOf("" + String.format("%4.0f", zMinLab)));
                rangeSliderValue2.setText(String.valueOf("" + String.format("%4.0f", zMaxLab)));
                plotHistos(getRunNumber());
            }
        });  
        JCheckBox sliderBtn = new JCheckBox("Enable");
        sliderBtn.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                doSlider = (e.getStateChange() == ItemEvent.SELECTED) ? true:false;
                plotHistos(getRunNumber());
            }
        });         
        sliderBtn.setSelected(false);         
        sliderPane.add(sliderBtn);        
        return sliderPane;
    }  
    
    public void actionPerformed(ActionEvent e) {
        // TODO Auto-generated method stub
        this.detectorActiveSector = Integer.parseInt(bG1.getSelection().getActionCommand());
        this.detectorActiveLayer  = Integer.parseInt(bG2.getSelection().getActionCommand());
        this.detectorActiveView   = Integer.parseInt(bG3.getSelection().getActionCommand());
        plotHistos(getRunNumber());
    } 
    
    public void processEvent(DataEvent event) {
        // process event
    }
    
    public void plotEvent(DataEvent event) {
        // process event
    }
    
    public void plotHistos(int run) {

    }
    
    public void printCanvas(String dir) {
        // print canvas to files
        int run = this.getViewRun();
        if(run==0) run = this.getRunNumber();
        for(int tab=0; tab<this.detectorTabNames.size(); tab++) {
            String fileName = dir + "/" + this.detectorName + "_" + run + "_canvas" + tab + ".png";
            System.out.println(fileName);
            this.detectorCanvas.getCanvas(this.detectorTabNames.get(tab)).save(fileName);
        }
    }   
    
    @Override
    public void resetEventListener() {
        System.out.println("Resetting " + this.getDetectorName() + " histogram for run "+this.getRunNumber());
        this.createHistos(this.getRunNumber());
        this.plotHistos(this.getRunNumber());
    }
    
    public void setCanvasUpdate(int time) {
        for(int tab=0; tab<this.detectorTabNames.size(); tab++) {
            this.detectorCanvas.getCanvas(this.detectorTabNames.get(tab)).initTimer(time);
        }
    }
    
    public void setDetectorCanvas(EmbeddedCanvasTabbed canvas) {
        this.detectorCanvas = canvas;
    }
    
    public void setDetectorTabNames(String... names) {
        for(String name : names) {
            this.detectorTabNames.add(name);
        }
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed(names);
        this.setDetectorCanvas(canvas);
    }
 
    public void setDetectorSummary(DataGroup group) {
        this.detectorSummary = group;
    }
    
    public void setNumberOfEvents(int numberOfEvents) {
        this.numberOfEvents = numberOfEvents;
    }
    
    public void setEventNumber(int eventNumber) {
        this.eventNumber = eventNumber;
    }

    public void setRunNumber(int runNumber) {
        this.runNumber = runNumber;
    }
    
    public int getEventNumber() {
        return eventNumber;
    }

    public int getRunNumber() {
        return runNumber;
    }

    public int getViewRun() {
        return viewRun;
    }   
    
    @Override
    public void timerUpdate() {      
    }
 
    public void readDataGroup(TDirectory dir) {
        String folder = this.getDetectorName() + "/";
        System.out.println("Reading from: " + folder);
        DataGroup sum = this.getDetectorSummary();
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
        this.setDetectorSummary(newSum);
        
        }
        
        Map<Long, DataGroup> map = this.getDataGroup().getMap();
        for( Map.Entry<Long, DataGroup> entry : map.entrySet()) {
            Long key = entry.getKey();
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
        this.plotHistos(this.getRunNumber());
    }
    
    public void writeDataGroup(TDirectory dir) {
        String folder = "/" + this.getDetectorName();
        dir.mkdir(folder);
        dir.cd(folder);
        DataGroup sum = this.getDetectorSummary();
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
        
        Map<Long, DataGroup> map = this.getDataGroup().getMap();
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
            if( doSlider) c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax);
            if(!doSlider) c.getPad().getAxisZ().setAutoScale(true);
            for (IDataSet ds : dsList) {
//                System.out.println("\t --> " + ds.getName());
            	   c.draw(ds,opt); opt="same";
            }
        } 	       	
    }
    
}
