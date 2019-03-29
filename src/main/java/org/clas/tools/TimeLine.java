package org.clas.tools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.data.TDirectory;
import org.jlab.utils.groups.IndexedList;

public class TimeLine {
	
	 public IndexedList<IDataSet> Timeline = new IndexedList<IDataSet>(2);
	 public IndexedList<FitData>   fitData = null;
	 private Map<Integer,Integer>   NYbins = new HashMap<Integer,Integer>();  
	 private TDirectory                dir = null;
	 
	 static void Timeline() {
		 
	 }
	 
	 public void defineTimeLine(String title, int k, int xbins, String ytit, int ybins, float ymin, float ymax) {
		 createTimeLineHisto(k,title,ytit,xbins,ybins,ymin,ymax);		 
	 }
	 
	 public void createTimeLineHisto(int k, String tit, String ytit, int nx, int ny, float ymin, float ymax) {
		 H2F h;
		 h =  new H2F(tit, tit, nx, 0, nx, ny, ymin, ymax);
		 h.setTitleX("Run Index") ; h.setTitleY(ytit); Timeline.add(h,k,0); //mean
		 h =  new H2F(tit, tit, nx, 0, nx, ny, ymin, ymax);
		 h.setTitleX("Run Index") ; h.setTitleY(ytit); Timeline.add(h,k,1); //error
	 }
	 
	 public void initTimelineFile() {
		 dir = new TDirectory();
		 dir.mkdir("/timelines");
	 }
	 
	 public void setFitData(IndexedList<FitData> fitData) {
		 this.fitData = fitData;
	 }
	 
	 public void setNYbins(int key, int value) {
		 NYbins.putIfAbsent(key, value);
	 }
	 
	 public Integer getNYbins(int key) {
		 return NYbins.get(key);	 
	 }
	 	 	
}
