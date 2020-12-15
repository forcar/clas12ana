package org.clas.tools;

import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;

public class DataGroupManager extends DataGroup {
	
	public EmbeddedCanvas canvas = null;
	
	public DataGroupManager(String name) {
		super(name);
	}
	
	public DataGroupManager(int ncols, int nrows) {
		super(ncols,nrows);
	}
	
	public DataGroupManager(String name, int ncols, int nrows, EmbeddedCanvas canvas) {
		super(name, ncols,nrows);
		this.canvas = canvas;
	}
	
	public void cc(int n, boolean linlogy, boolean linlogz, float ymin, float ymax, float zmin, float zmax) {
		canvas.cd(n);
		canvas.getPad().getAxisY().setLog(linlogy);
		canvas.getPad().getAxisZ().setLog(linlogz);
		if(ymin!=ymax) canvas.getPad().getAxisY().setRange(ymin,ymax);
		if(zmin!=zmax) canvas.getPad().getAxisZ().setRange(zmin,zmax);
	}
	
	public EmbeddedCanvas getCanvas() {
		return this.canvas;
	}
	

}
