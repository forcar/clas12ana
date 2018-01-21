package org.clas.viewer;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSplitPane;
import org.jlab.detector.base.DetectorOccupancy;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.task.IDataEventListener;
import org.jlab.utils.groups.IndexedList;


public class DetectorMonitor implements IDataEventListener, ActionListener {    
    
    private final String           detectorName;
    private ArrayList<String>      detectorTabNames  = new ArrayList();
    private IndexedList<DataGroup> detectorData      = new IndexedList<DataGroup>(3);
    private DataGroup              detectorSummary   = null;
    private DetectorOccupancy      detectorOccupancy = new DetectorOccupancy();
    private JPanel                 detectorPanel     = null;
    private EmbeddedCanvasTabbed   detectorCanvas    = null;
    private DetectorPane2D         detectorView      = null;
    private ButtonGroup            bG1               = null;
    private int                    numberOfEvents;
    private Boolean                sectorButtons     = false;
    private int                 detectorActiveSector = 1;
    private Boolean                     detectorLogZ = true;
    private Boolean                             isTB = false;
    private JRadioButton bS1,bS2,bS3,bS4,bS5,bS6;
    private JCheckBox        tbBtn;
    
    public int bitsec = 0;
    public long trigger = 0;
    public long triggerPhase = 0;
    public long DC_max_occ = 15;
    public int trigFD = 0;
    public int trigCD = 0;
    
    public boolean testTrigger = false;
    public boolean TriggerBeam[] = new boolean[32];
    public int TriggerMask = 0;
    
    public int eventResetTime_current[]=new int[19];
    public int eventResetTime_default[]=new int[19];    
    
    public DetectorMonitor(String name){
        GStyle.getAxisAttributesX().setTitleFontSize(18);
        GStyle.getAxisAttributesX().setLabelFontSize(14);
        GStyle.getAxisAttributesY().setTitleFontSize(18);
        GStyle.getAxisAttributesY().setLabelFontSize(14);
        this.detectorName = name;
        this.detectorPanel  = new JPanel();
        this.detectorCanvas = new EmbeddedCanvasTabbed();
        this.detectorView   = new DetectorPane2D();
        this.numberOfEvents = 0;
        
        eventResetTime_current[0]=0;
        eventResetTime_current[1]=0;
        eventResetTime_current[2]=0;
        eventResetTime_current[3]=0;
        eventResetTime_current[4]=0;
        eventResetTime_current[5]=0;
        eventResetTime_current[6]=0;
        eventResetTime_current[7]=0;
        eventResetTime_current[8]=0;
        eventResetTime_current[9]=0;
        eventResetTime_current[10]=0;
        eventResetTime_current[11]=0;
        eventResetTime_current[12]=0;
        eventResetTime_current[13]=0;
        eventResetTime_current[14]=0;
        eventResetTime_current[15]=0;
        eventResetTime_current[16]=0;
        eventResetTime_current[17]=0;
        eventResetTime_current[18]=0;
        
        eventResetTime_default[0]=10000000;   // BMT
        eventResetTime_default[1]=10000000;   // BST
        eventResetTime_default[2]=10000000;   // CND
        eventResetTime_default[3]=10000000;   // CTOF
        eventResetTime_default[4]=10000000;   // DC
        eventResetTime_default[5]=10000000;   // ECAL
        eventResetTime_default[6]=10000000;   // FMT
        eventResetTime_default[7]=10000000;   // FTCAL
        eventResetTime_default[8]=10000000;   // FTHODO
        eventResetTime_default[9]=10000000;   // FTOF
        eventResetTime_default[10]=10000000;  // FTTRK
        eventResetTime_default[11]=10000000;  // HTTC
        eventResetTime_default[12]=10000000;  // LTTC
        eventResetTime_default[13]=10000000;  // RICH
        eventResetTime_default[14]=10000000;  // RECONN
        eventResetTime_default[15]=10000000;  // RF
        eventResetTime_default[16]=10000000;  // HEL
        eventResetTime_default[17]=10000000;  // Faraday Cup
        eventResetTime_default[18]=10000000;  // Trigger
     
        for (int i=0; i<19; i++){
            eventResetTime_current[i]=eventResetTime_default[i];
        }
    }

    
    public void analyze() {
        // analyze detector data at the end of data processing
    }

    public void createHistos() {
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
/*    
    public boolean isGoodCDTrigger()       {return (testTrigger)? isGoodCD():true;}  
    public boolean isGoodHTCCTrigger()     {return (testTrigger)? isGoodHTCC():true;}
    public boolean isGoodFTOFTrigger()     {return (testTrigger)? isGoodFTOF():true;}
    public boolean isGoodBSTTrigger()      {return (testTrigger)? isGoodBST():true;}
    public boolean isGoodCTOFTrigger()     {return (testTrigger)? isGoodCTOF():true;}
    public boolean isGoodCNDTrigger()      {return (testTrigger)? isGoodCND():true;}
    public boolean isGoodBMTTrigger()      {return (testTrigger)? isGoodBMT():true;}
    public boolean isGoodFD()              {return  this.trigFD>=256&&this.trigFD<=8196;}    
    public boolean isGoodCD()              {return  isGoodBST()||isGoodCTOF()||isGoodCND()||isGoodBMT();}
    public boolean isGoodHTCC()            {return  this.trigFD==1;}    
    public boolean isGoodFTOF()            {return  isGoodFD();}   
    public boolean isGoodBST()             {return  this.trigCD==256;}    
    public boolean isGoodCTOF()            {return  this.trigCD==512;}
    public boolean isGoodCND()             {return  this.trigCD==1024;}   
    public boolean isGoodBMT()             {return  this.trigCD==2048;}       
    public int     getFDTrigger()          {return (this.trigger>>16)&0x0000ffff;}    
    public int     getCDTrigger()          {return this.trigger&0x00000fff;} 
*/
    
    public int     getFDTrigger()            {return (int)(this.trigger)&0x000000000ffffffff;}
    public int     getCDTrigger()            {return (int)(this.trigger>>32)&0x00000000ffffffff;}
    public boolean isGoodFD()                {return  getFDTrigger()>0;}    
    public boolean isTrigBitSet(int bit)     {int mask=0; mask |= 1<<bit; return isTrigMaskSet(mask);}
    public boolean isTrigMaskSet(int mask)   {return (getFDTrigger()&mask)!=0;}
    public boolean isGoodECALTrigger(int is) {return (testTrigger)? is==getECALTriggerSector():true;}    
    public int           getElecTrigger()    {return getFDTrigger()&0x1;}
    public int     getElecTriggerSector()    {return (int) (isGoodFD() ? Math.log10(getFDTrigger()>>1)/0.301+1:0);} 
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
    
    public int getActiveSector() {
    	    return detectorActiveSector;
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

    public void init(boolean flagDetectorView) {
        // initialize monitoring application
        // detector view is shown if flag is true
        getDetectorPanel().setLayout(new BorderLayout());
        drawDetector();
        JSplitPane   splitPane = new JSplitPane();
        splitPane.setLeftComponent(getDetectorView());
        splitPane.setRightComponent(getDetectorCanvas());
        if(flagDetectorView) {
            getDetectorPanel().add(splitPane,BorderLayout.CENTER);  
        }
        else {
            getDetectorPanel().add(getDetectorCanvas(),BorderLayout.CENTER); 
            if (sectorButtons) getDetectorPanel().add(getButtonPane(),BorderLayout.PAGE_END);  
        }
        createHistos();
        plotHistos(); 
        if (sectorButtons) bS2.doClick();
    }
    
    public JPanel getButtonPane() {
    	    bG1 = new ButtonGroup();
        JPanel buttonPane = new JPanel();
        bS1 = new JRadioButton("Sector 1"); buttonPane.add(bS1); bS1.setActionCommand("1"); bS1.addActionListener(this);
        bS2 = new JRadioButton("Sector 2"); buttonPane.add(bS2); bS2.setActionCommand("2"); bS2.addActionListener(this); 
        bS3 = new JRadioButton("Sector 3"); buttonPane.add(bS3); bS3.setActionCommand("3"); bS3.addActionListener(this); 
        bS4 = new JRadioButton("Sector 4"); buttonPane.add(bS4); bS4.setActionCommand("4"); bS4.addActionListener(this); 
        bS5 = new JRadioButton("Sector 5"); buttonPane.add(bS5); bS5.setActionCommand("5"); bS5.addActionListener(this);  
        bS6 = new JRadioButton("Sector 6"); buttonPane.add(bS6); bS6.setActionCommand("6"); bS6.addActionListener(this); 
        bG1.add(bS1);bG1.add(bS2);bG1.add(bS3);bG1.add(bS4);bG1.add(bS5);bG1.add(bS6);
        //tbBtn = new JCheckBox("TrigBit");
        //tbBtn.addItemListener(new ItemListener() {
        //    public void itemStateChanged(ItemEvent e) {
        //    	  isTB = e.getStateChange()==ItemEvent.SELECTED;
        //    }
        //});         
        //tbBtn.setSelected(false);        
        //buttonPane.add(tbBtn);       
        return buttonPane;
    } 
    
    public void actionPerformed(ActionEvent e) {
        // TODO Auto-generated method stub
        this.detectorActiveSector   = Integer.parseInt(bG1.getSelection().getActionCommand());
        plotHistos();
    } 
    
    public void processEvent(DataEvent event) {
        // process event
    }
    
    public void plotEvent(DataEvent event) {
        // process event
    }
    
    public void plotHistos() {

    }
    
    public void printCanvas(String dir) {
        // print canvas to files
        for(int tab=0; tab<this.detectorTabNames.size(); tab++) {
            String fileName = dir + "/" + this.detectorName + "_canvas" + tab + ".png";
            System.out.println(fileName);
            this.detectorCanvas.getCanvas(this.detectorTabNames.get(tab)).save(fileName);
        }
    }
    
    @Override
    public void resetEventListener() {
        System.out.println("Resetting " + this.getDetectorName() + " histogram");
        this.createHistos();
        this.plotHistos();
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

    @Override
    public void timerUpdate() {
        
    }
 
    public void readDataGroup(TDirectory dir) {
        String folder = this.getDetectorName() + "/";
        System.out.println("Reading from: " + folder);
        DataGroup sum = this.getDetectorSummary();
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
        Map<Long, DataGroup> map = this.getDataGroup().getMap();
        for( Map.Entry<Long, DataGroup> entry : map.entrySet()) {
            Long key = entry.getKey();
            DataGroup group = entry.getValue();
            nrows = group.getRows();
            ncols = group.getColumns();
            nds   = nrows*ncols;
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
        this.plotHistos();
    }
    
    public void writeDataGroup(TDirectory dir) {
        String folder = "/" + this.getDetectorName();
        dir.mkdir(folder);
        dir.cd(folder);
        DataGroup sum = this.getDetectorSummary();
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
        Map<Long, DataGroup> map = this.getDataGroup().getMap();
        for( Map.Entry<Long, DataGroup> entry : map.entrySet()) {
            DataGroup group = entry.getValue();
            nrows = group.getRows();
            ncols = group.getColumns();
            nds   = nrows*ncols;
            for(int i = 0; i < nds; i++){
                List<IDataSet> dsList = group.getData(i);
                for(IDataSet ds : dsList){
                    System.out.println("\t --> " + ds.getName());
                    dir.addDataSet(ds);
                }
            }
        }
    }
        
}
