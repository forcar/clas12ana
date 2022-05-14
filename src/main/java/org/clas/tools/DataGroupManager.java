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
	public Map<String,Boolean>                  tmap = new HashMap<>();
	public Map<String,H1F[]>                   r1map = new HashMap<>();
	public Map<String,H2F[]>                   r2map = new HashMap<>();
	
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
		tag = "-"+is+"-"+st+"-"+k+"-"+run;
		System.out.println("dgm.add("+name+tag+")");
		detectorData.add(dg,ilist); 
		n=0;
	}
	
	public DataGroup setDataGroup(String name, int is, int st, int run) {
		return detectorData.getItem(is,st,detectorTabNames.indexOf(name),run);
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
    
    public IndexedList<DataGroup> getDataGroup(){
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
    
    public int getRun() {
    	return ilist[3];
    }
    
    public String getName(String name) {
		imap.put(name,get_ilist());
		hmap.put(name,name+tag);
		tmap.put(name,true);
		return hmap.get(name);
    }
    
    public void setTest(String name,boolean val) {
    	tmap.put(name,val);    	
    }
    
    public void printMap(String name, int i) {
    	int[] dum = imap.get(name); System.out.println(i+" "+name+" "+dum[0]+" "+dum[1]+" "+dum[2]+" "+dum[3]);
    	System.out.println(getDataGroup().hasItem(imap.get(name))+" "+ hmap.get(name));     	
    }
    
    public boolean isH1(String name) {
    	return (getDataGroup().getItem(imap.get(name)).getData(hmap.get(name)) instanceof H1F);
    }
    
    public boolean isH2(String name) {
    	return (getDataGroup().getItem(imap.get(name)).getData(hmap.get(name)) instanceof H2F);
    }
    
    public boolean isH1R(String name) {
    	return (r1map.containsKey(name));
    }
    
    public boolean isH2R(String name) {
    	return (r2map.containsKey(name));
    }
    
    public boolean isGraph(String name) {
    	return (getDataGroup().getItem(imap.get(name)).getData(hmap.get(name)) instanceof GraphErrors);
    }
    
   
    
// HISTO HELPERS  
    
    public void fill(String name, double ... val) {    	
        if(!goodName(name)) return; 
        if(!tmap.get(name)) return;
        if(val.length==3 && isH1R(name)) {if(val[2]==1) {getH1FR(name)[0].fill(val[0],val[1]);} getH1FR(name)[1].fill(val[0]);        return;}
        if(val.length==2 && isH1R(name)) {getH1FR(name)[0].fill(val[0],val[1]);        getH1FR(name)[1].fill(val[0]);        return;}
        if(val.length==3 && isH2R(name)) {getH2FR(name)[0].fill(val[0],val[1],val[2]); getH2FR(name)[1].fill(val[0],val[1]); return;} 
    	if(val.length==1)                getH1F(name).fill(val[0]);
    	if(val.length==2 && isH1(name))  getH1F(name).fill(val[0],val[1]);
        if(val.length==2 && isH2(name))  getH2F(name).fill(val[0],val[1]);
    	if(val.length==3)                getH2F(name).fill(val[0],val[1],val[2]);
    	return;
    }
    
    public void draw(String name, EmbeddedCanvas c, int is, int st, int i) {      
        int index = getDetectorTabNames().indexOf(name);
        drawGroupItem(getDataGroup().getItem(is,st,index,getRun()),c,i);        
    }
    
    public boolean goodName(String name) {
    	if(!imap.containsKey(name))                 {System.out.println(name+" missing from imap"); return false;} 
    	if(!getDataGroup().hasItem(imap.get(name))) {System.out.println(imap.get(name)+" missing from detectorData"); return false;} 
    	if(!hmap.containsKey(name))                 {System.out.println(name+" missing from hmap"); return false;} 
    	return true;
    }
    
    public H1F[] getH1FR(String name) {
    	return r1map.get(name);
    }
    
    public H2F[] getH2FR(String name) {
    	return r2map.get(name);
    }
    
    public H1F getH1F(String name) {
    	return getDataGroup().getItem(imap.get(name)).getH1F(hmap.get(name));
    }
    
    public H2F getH2F(String name) {
    	return getDataGroup().getItem(imap.get(name)).getH2F(hmap.get(name));
    }
    
    public GraphErrors getGraph(String name) {
    	return getDataGroup().getItem(imap.get(name)).getGraph(hmap.get(name));    	
    }
    
    public void geteff(String eff, String ... val) {
    	if(isH1R(eff))   {getH1F(eff).reset(); getH1F(eff).add(H1F.divide(getH1FR(eff)[0],getH1FR(eff)[1])); return;}
    	if(isH2R(eff))   {getH2F(eff).reset(); getH2F(eff).add(H2F.divide(getH2FR(eff)[0],getH2FR(eff)[1])); return;}
    	if(isH1(eff))    {getH1F(eff).reset(); getH1F(eff).add(H1F.divide(getH1F(val[0]), getH1F(val[1])));}
    	if(isH2(eff))    {getH2F(eff).reset(); getH2F(eff).add(H2F.divide(getH2F(val[0]), getH2F(val[1])));}
    	if(isGraph(eff)) {getGraph(eff).reset(); getGraph(eff).copy(H1F.divide(getH1F(val[0]), getH1F(val[1])).getGraph());}
    }
    
    public class H2FF extends H2F {
   	 
   	 public H2FF(String str1, String str2, int nbx, double x1, double x2, int nby, double y1, double y2) {
   		 super("H2FF");
   		 setName(str1); setTitle(str2);
   		 set(nbx,x1,x2,nby,y1,y2);
   	 }
   	 
    }
    
    public void makeH1R(String name, int nx, double x1, double x2, int f, String tit, String titx, String tity, Object ... color) {
    	makeH1(name, nx, x1, x2, f, tit, titx, tity, color);
    	H1F[] hr = {new H1F(name,name,nx,x1,x2),new H1F(name,name,nx,x1,x2)};
    	r1map.put(name,hr);
    }
    
    public void makeH2R(String name, int nx, double x1, double x2, int ny, double y1, double y2, int f, String tit, String titx, String tity) {
    	makeH2(name, nx, x1, x2, ny, y1, y2, f, tit, titx, tity);
    	H2F[] hr = {new H2F(name,name,nx,x1,x2,ny,y1,y2),new H2F(name,name,nx,x1,x2,ny,y1,y2)};
    	r2map.put(name,hr);
    }
    
    public void makeH1(String name) { //use this to force zone based display of H1 created with overlap tag f=-2
    	dg.addDataSet(getH1F(name),n++);
    }
    
    public void makeH1(String name, int nx, double x1, double x2, int f, String tit, String titx, Object ... color) {
	    H1F h1 = new H1F(getName(name),getName(name),nx,x1,x2);
	    if(tit!="") h1.setTitle(tit);
	    h1.setTitleX(titx);
	    int nn=0;
        for (Object col : color) {
    	    if(nn==0) h1.setLineColor((Integer)col); 
    	    if(nn==1) h1.setFillColor((Integer)col);
    	    if(nn==2) h1.setOptStat((String)col);
    	    nn++;
        }
        dg.addDataSet(h1, (f==-2)? n-1: n++);
	    if(f>=0)  makeF1D("f"+getName(name),"[a]",x1,x2,f);	   
	    System.out.println("makeH1("+hmap.get(name)+" "+n+")");
    }
    
    public void makeH1(String name, int nx, double x1, double x2, int f, String tit, String titx, String tity, Object ... color) {
	    H1F h1 = new H1F(getName(name),getName(name),nx,x1,x2);
	    if(tit!="") h1.setTitle(tit);
	    h1.setTitleX(titx);  h1.setTitleY(tity);
	    int nn=0;	    
        for (Object col : color) {
    	    if(nn==0) h1.setLineColor((Integer)col); 
    	    if(nn==1) h1.setFillColor((Integer)col);
    	    if(nn==2) h1.setOptStat((String)col);
   	        nn++;
        }
        dg.addDataSet(h1, (f==-2)? n-1: n++);
	    if(f>=0)  makeF1D("f"+getName(name),"[a]",x1,x2,f);
	    System.out.println("makeH1("+hmap.get(name)+" "+n+")");
    }    
    
    public void makeH2(String name, int nx, double x1, double x2, int ny, double y1, double y2, int f, String tit, String titx, String tity) {
	    H2F h2 = new H2F(getName(name),getName(name),nx,x1,x2,ny,y1,y2);
	    if(tit!="") h2.setTitle(tit);
	    h2.setTitleX(titx); h2.setTitleY(tity);
	    dg.addDataSet(h2, n++);
	    if(f!=-1) makeF1D("f"+getName(name),"[a]",x1,x2,f);
	    System.out.println("makeH2("+hmap.get(name)+" "+n+")");
    }
    
    public void makeF1D(String name, String f, double x1, double x2, double val) {
    	F1D f1 = new F1D(name,f,x1,x2); f1.setParameter(0,val);
    	dg.addDataSet(f1, n-1);
    }
    
    public void makeGraph(String name, int f, String tit, String titx, String tity, int ... options) {
    	GraphErrors graph = new GraphErrors();
    	graph.setName(getName(name)); 
    	graph.setTitle(tit); graph.setTitleX(titx); graph.setTitleY(tity);
    	int nn=0;
    	for (int opt : options) {
    		if(nn==0) {graph.setMarkerColor(opt); graph.setLineColor(opt);}
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

    // custom configuration: associates canvas based drawing options with IDataSet object name
	public void cc(String name, boolean linlogy, boolean linlogz, float ymin, float ymax, float zmin, float zmax) {
		Object[] obj = {linlogy,linlogz,ymin,ymax,zmin,zmax};
		map.put(hmap.get(name),obj);  
	}
	
	//custom configuration: unpacks IDataSet object map and performs canvas operations at draw time
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
    	drawGroup(getDetectorCanvas().getCanvas(getDetectorTabNames().get(index)),
    			  getDataGroup().getItem(is,st,index,run));   			
	}
    
    public void drawGroup(EmbeddedCanvas c, DataGroup group) {
    	if(group==null) return;
        int nrows = group.getRows(), ncols = group.getColumns();
        c.clear(); c.divide(ncols, nrows); 	    
        int nds = nrows * ncols;
        for (int i = 0; i < nds; i++) {
          c.cd(i); drawGroupItem(group,c,i);                 
        } 	       	
    }
    
    public void drawGroupItem(DataGroup group, EmbeddedCanvas c, int i) {
    	boolean maxtest = false; String opt = " "; 
    	List<IDataSet> dsList = group.getData(i); 
        if(dsList.size()>1) maxtest = true; //prevents AutoScale from triggering on empty or low yield histogram
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
        		if (ds0 instanceof GraphErrors) {
        			maxtest = false;
        			if( doAutoRange) {c.getPad().getAxisX().setAutoScale(true);
                                      c.getPad().getAxisY().setAutoScale(true);
                                      c.getPad().getAxisZ().setAutoScale(true);}            			
        		}        		
        	}
        	for (IDataSet ds : dsList) if(maxtest?ds.getMax()>0.5:true) {c.draw(ds,opt); opt="same";}
        }   	
    }
	
}
