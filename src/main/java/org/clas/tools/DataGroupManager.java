package org.clas.tools;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.utils.groups.IndexedList;

public class DataGroupManager {
	
	private IndexedList<DataGroup>  detectorData = new IndexedList<DataGroup>(4);
	Map<String,Object[]> map = new HashMap<String,Object[]>();
	Boolean  doAutoRange = false;
	Boolean detectorLogY = false;
	Boolean detectorLogZ = false;
	Object[] can = {false, false, 0, 0, 0, 0};
	public float zMin=0, zMax=0;

	public DataGroupManager() {
		
	}

	public void cc(String name, boolean linlogy, boolean linlogz, float ymin, float ymax, float zmin, float zmax) {
		Object[] obj = {linlogy,linlogz,ymin,ymax,zmin,zmax};
		map.put(name,obj);  
	}
	
	public Boolean config(String name, EmbeddedCanvas canvas) {
		if(!map.containsKey(name)) return false;
		can = map.get(name); 
		canvas.getPad().getAxisY().setLog((Boolean)can[0]);
		canvas.getPad().getAxisZ().setLog((Boolean)can[1]);
		float ymin=(float)can[1], ymax=(float)can[2], zmin=(float)can[3], zmax=(float)can[4]; 
		if(ymin!=ymax) canvas.getPad().getAxisY().setRange(ymin,ymax);
		if(zmin!=zmax) canvas.getPad().getAxisZ().setRange(zmin,zmax);
		return true;
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
