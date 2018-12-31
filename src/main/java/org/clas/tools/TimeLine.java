package org.clas.tools;

import org.jlab.groot.data.H2F;
import org.jlab.groot.data.IDataSet;
import org.jlab.utils.groups.IndexedList;

public class TimeLine {
	
	 IndexedList<IDataSet>       Timeline = new IndexedList<IDataSet>(2);
	 
	 static void Timeline() {
		 
	 }
	 
	 public void defineTimeLine(String title, int k, int xbins, String ytit, int ybins, float ymin, float ymax) {
		 createTimeLineHisto(k,title,ytit,xbins,ybins,ymin,ymax);		 
	 }
	 
	 public void createTimeLineHisto(int k, String tit, String ytit, int nx, int ny, float ymin, float ymax) {
		 H2F h;
		 h =  new H2F(tit, tit, nx, 0, nx, ny, 0, ny);
		 h.setTitleX("Run Index") ; h.setTitleY(ytit); Timeline.add(h,k,0); //mean
		 h =  new H2F(tit, tit, nx, 0, nx, ny, ymin, ymax);
		 h.setTitleX("Run Index") ; h.setTitleY(ytit); Timeline.add(h,k,1); //error
	 }

	
}
