package org.clas.tools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.utils.groups.IndexedList;

public class DataGroupManager {
	
	public IndexedList<DataGroup>       detectorData = new IndexedList<>(4);
	private EmbeddedCanvasTabbed      detectorCanvas = null;
	public ArrayList<String>        detectorTabNames = new ArrayList();
	Map<String,Object[]>                         map = new HashMap<>();
	public Map<String,int[]>                    imap = new HashMap<>();
	public Map<String,String>                   hmap = new HashMap<>();
	
	DataGroup dg = null;
	String   tag = null;
	int  ilist[] = null, n=0;
	public double zMin=0, zMax=0;
    public Boolean  doAutoRange = false;
	public Boolean detectorLogY = false;
	public Boolean detectorLogZ = false;
	
	public DataGroupManager() {
		
	}
	
	public void add(String name, int row, int col, int is, int st, int run) {
		dg = new DataGroup(row,col); 
		int k = detectorTabNames.indexOf(name);
		ilist = new int[4]; ilist[0]=is; ilist[1]=st; ilist[2]=k; ilist[3]=run;
		tag = "_"+is+"_"+st+"_"+k+"_"+run;
		System.out.println("dgm.add("+name+tag+")");
		detectorData.add(dg,ilist);
		System.out.println("Done");
		n=0;
	}
	
	public void setDataGroup(String name, int is, int st, int run) {
		dg = detectorData.getItem(is,st,detectorTabNames.indexOf(name),run);
	}
	
	public void setDetectorCanvas(EmbeddedCanvasTabbed canvas) {
		detectorCanvas = canvas;
	}
	
    public void setDetectorTabNames(String... names) {
        setDetectorCanvas(new EmbeddedCanvasTabbed(names));    	      
        for(String name : names) detectorTabNames.add(name);     
    }
    
    public void setLogY(boolean val) {
	    detectorLogY = val;
    }
    
    public void setLogZ(boolean val) {
	    detectorLogZ = val;
    }
    
    public Boolean getLogY() {
	    return detectorLogY;
    }
    
    public Boolean getLogZ() {
	    return detectorLogZ;
    }
    
    public IndexedList<DataGroup>  getDataGroup(){
        return detectorData;
    }
    
    public ArrayList<String> getDetectorTabNames() {
        return detectorTabNames;
    }
    
    public EmbeddedCanvasTabbed getDetectorCanvas() {
        return detectorCanvas;
    }
    
    public int[] get_ilist() {
    	return ilist;
    }
    
    public String getName(String name) {
		imap.put(name,get_ilist());
		hmap.put(name,name+tag);		
		return hmap.get(name);
    }
    
    public void printMap(String name, int i) {
    	int[] dum = imap.get(name); System.out.println(i+" "+name+" "+dum[0]+" "+dum[1]+" "+dum[2]+" "+dum[3]);
    	System.out.println(detectorData.hasItem(imap.get(name))+" "+ hmap.get(name));     	
    }
    
    public boolean isH1(String name) {
    	return (detectorData.getItem(imap.get(name)).getData(hmap.get(name)) instanceof H1F);
   }
    
    public boolean isH2(String name) {
    	return (detectorData.getItem(imap.get(name)).getData(hmap.get(name)) instanceof H2F);
   }
    
   
    
// HISTO HELPERS  
    
    public void fill(String name, double ... val) {    	
    	if(!imap.containsKey(name))               {System.out.println(name+" missing from imap"); return;} 
    	if(!detectorData.hasItem(imap.get(name))) {System.out.println(imap.get(name)+" missing from detectorData"); return;} 
    	if(!hmap.containsKey(name))               {System.out.println(name+" missing from hmap"); return;}     	
    	if(val.length==1)               detectorData.getItem(imap.get(name)).getH1F(hmap.get(name)).fill(val[0]);
    	if(val.length==2 && isH1(name)) detectorData.getItem(imap.get(name)).getH1F(hmap.get(name)).fill(val[0],val[1]);
        if(val.length==2 && isH2(name)) detectorData.getItem(imap.get(name)).getH2F(hmap.get(name)).fill(val[0],val[1]);
    	if(val.length==3)               detectorData.getItem(imap.get(name)).getH2F(hmap.get(name)).fill(val[0],val[1],val[2]);
    	return;
    }
    
    public H1F getH1F(String name) {
    	return detectorData.getItem(imap.get(name)).getH1F(hmap.get(name));
    }
    
    public H2F getH2F(String name) {
    	return detectorData.getItem(imap.get(name)).getH2F(hmap.get(name));
    }
    
    public GraphErrors getGraph(String name) {
    	return detectorData.getItem(imap.get(name)).getGraph(hmap.get(name));    	
    }
    
    public void geteff(String eff, String numer, String denom) {
    	if(isH1(eff)) getH1F(eff).reset();
    	if(isH2(eff)) getH2F(eff).reset();
    	if(isH1(eff)) getH1F(eff).add(H1F.divide(getH1F(numer), getH1F(denom)));
    	if(isH2(eff)) getH2F(eff).add(H2F.divide(getH2F(numer), getH2F(denom)));
    }
    
    public void makeH1(String name) { //use this to force zone based display of H1 created with overlap tag f=-2
    	dg.addDataSet(getH1F(name),n++);
    }
    
    public void makeH1(String name, int nx, double x1, double x2, int f, String tit, String titx, int ... color) {
	    H1F h1 = new H1F(getName(name),nx,x1,x2);
	    if(tit!="") h1.setTitle(tit);
	    h1.setTitleX(titx);
	    int nn=0;
        for (int col : color) {
    	    if(nn==0) h1.setLineColor(col); 
    	    if(nn==1) h1.setFillColor(col);
    	    if(nn==2 && col==1) h1.setOptStat("1000000");
    	    if(nn==2 && col==2) h1.setOptStat("1000100");
    	    nn++;
        }
        dg.addDataSet(h1, (f==-2)? n-1: n++);
	    if(f>=0)  dg.addDataSet(makeF1D("f"+n+tag,x1,x2,f), n-1); 	   
	    System.out.println("makeH1("+hmap.get(name)+" "+n+")");
    }
    
    public void makeH1(String name, int nx, double x1, double x2, int f, String tit, String titx, String tity, int ... color) {
	    H1F h1 = new H1F(getName(name),nx,x1,x2);
	    if(tit!="") h1.setTitle(tit);
	    h1.setTitleX(titx);  h1.setTitleY(tity);
	    int nn=0;	    
        for (int col : color) {
    	    if(nn==0) h1.setLineColor(col); 
    	    if(nn==1) h1.setFillColor(col);
    	    if(nn==2 && col==1) h1.setOptStat("1000000");
    	    if(nn==2 && col==2) h1.setOptStat("1000100");
   	        nn++;
        }
        dg.addDataSet(h1, (f==-2)? n-1: n++);
	    if(f>=0)  dg.addDataSet(makeF1D("f"+n+tag,x1,x2,f), n-1);
	    System.out.println("makeH1("+hmap.get(name)+" "+n+")");
    }    
    
    public void makeH2(String name, int nx, double x1, double x2, int ny, double y1, double y2, int f, String tit, String titx, String tity) {
	    H2F h2 = new H2F(getName(name),nx,x1,x2,ny,y1,y2);
	    if(tit!="") h2.setTitle(tit);
	    h2.setTitleX(titx); h2.setTitleY(tity);
	    dg.addDataSet(h2, n++);
	    if(f!=-1) dg.addDataSet(makeF1D("f"+n+tag,x1,x2,f), n-1);
	    System.out.println("makeH2("+hmap.get(name)+" "+n+")");
    }
    
    public F1D makeF1D(String name, double x1, double x2, double val) {
    	F1D f1 = new F1D("p1","[a]",x1,x2); f1.setParameter(0,val);
    	return f1;
    }
    
    public void makeGraph(String name, int f, int ... options) {
    	GraphErrors graph = new GraphErrors();
    	graph.setName(getName(name)); 
    	int nn=0;
    	for (int opt : options) {
    		if(nn==0) graph.setMarkerColor(opt); 
    		if(nn==1) graph.setMarkerSize(opt); 
    		if(nn==2) graph.setMarkerStyle(opt);
    		nn++;
    	}    	
    	dg.addDataSet(graph, (f==-2)? n-1: n++);
    }
    
    public void addDataSet(IDataSet ds, int f) {
    	dg.addDataSet(ds, (f==-2)? n-1: n++);
    }
    
    //DATAGROUP HELPERS	    

	public void cc(String name, boolean linlogy, boolean linlogz, float ymin, float ymax, float zmin, float zmax) {
		Object[] obj = {linlogy,linlogz,ymin,ymax,zmin,zmax};
		map.put(hmap.get(name),obj);  
	}
	
	public Boolean config(String name, EmbeddedCanvas canvas) {
		if(!map.containsKey(name)) return false;                  
		Object[] can = map.get(name);            
		canvas.getPad().getAxisY().setLog((Boolean)can[0]);               
		canvas.getPad().getAxisZ().setLog((Boolean)can[1]);            
		float ymin=(float)can[2], ymax=(float)can[3], zmin=(float)can[4], zmax=(float)can[5];            
		if(ymin!=ymax) canvas.getPad().getAxisY().setRange(ymin,ymax);                
		if(zmin!=zmax) canvas.getPad().getAxisZ().setRange(zmin,zmax);                
		return true;
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
            List<IDataSet> dsList = group.getData(i);           
            c.cd(i);  String opt = " ";                  
            if(dsList.size()>0) {	               
          	    IDataSet ds0 = dsList.get(0);                   
            	if(!config(ds0.getName(),c)) {                      
            		if (ds0 instanceof H1F) {
            			c.getPad().getAxisY().setAutoScale(true);
            			c.getPad().getAxisY().setLog(getLogY());
            		}                
            		if (ds0 instanceof H2F) {
            			c.getPad().getAxisZ().setLog(getLogZ());
            			if(!doAutoRange) c.getPad().getAxisZ().setRange(0.1*zMin, 20*zMax); 
            			if( doAutoRange) c.getPad().getAxisZ().setAutoScale(true);
            		}            		
            	}
            	for (IDataSet ds : dsList) {c.draw(ds,opt); opt="same";}
            }
        } 	       	
    } 
	
}
