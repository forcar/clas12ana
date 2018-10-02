package org.clas.tools;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;

public class FitData {
	
	public GraphErrors graph = null;
	public H1F          hist = null;

	public double xmin;
	public double xmax;
	public double amp;
	public double mean;
	public double meane;
	public double sigma;
	public int integral;
	public int intmin=30;
	public int fitcol=4;
	public String optstat="1100";
	public Boolean doFit = false;

	public FitData(GraphErrors graph, double xmin, double xmax) {
	    setGraph(graph);
	    this.graph.setFunction(new F1D("f","[amp]*gaus(x,[mean],[sigma])", xmin, xmax));            
	    this.graph.getFunction().setLineWidth(1);
	};

	public void setGraph(GraphErrors graph) {
	    this.graph = graph;
	}
	
	public void setHist(H1F hist) {
		this.hist = hist;
		this.hist.setLineWidth(1);
		this.hist.setOptStat("100");
	}
	
	public H1F getHist() {
		return this.hist;
	}
	
	public double getMean() {
		return this.hist.getMean();
	}

	public void setInt(int integral) {
	    this.integral = integral;
	    fitcol = 4;
	    if(integral>intmin) {
	        doFit = true;
	        fitcol = 1;
	    }
	}

	public void setIntMin(int intmin) {
	    this.intmin = intmin;
	}

	public void initFit(double xmin, double xmax) {
	    this.xmin = xmin;
	    this.xmax = xmax;
	    amp   = getMaxYIDataSet(graph,xmin,xmax);
	    mean  = getMeanIDataSet(graph,xmin,xmax); 
	    sigma = getRMSIDataSet(graph,xmin,xmax);
	    graph.getFunction().setParameter(0, amp);
	    graph.getFunction().setParameter(1, mean);
	    graph.getFunction().setParameter(2, sigma);
	    graph.getFunction().setRange(mean-2.5*sigma, mean+1.7*sigma);
	 }

	public void fitGraph(String opt) {
	    if (doFit) DataFitter.fit(graph.getFunction(), graph, "Q");
	    amp   = graph.getFunction().getParameter(0);
	    mean  = graph.getFunction().parameter(1).value();
	    meane = graph.getFunction().parameter(1).error();
	    sigma = graph.getFunction().getParameter(2);   
	    graph.getFunction().setLineColor(fitcol);
	    graph.getFunction().setOptStat(opt=="Q"?"0":optstat);
	}

	public void plotGraph(EmbeddedCanvas c, int col) {
	    graph.getFunction().setLineColor(col); 
	    c.draw(graph);
	}

	public GraphErrors getGraph() {
	    return graph;
	}

	private double getMaxYIDataSet(IDataSet data, double min, double max) {
	    double max1 = 0;
	    double xMax = 0;
	    for (int i = 0; i < data.getDataSize(0); i++) {
	        double x = data.getDataX(i);
	        double y = data.getDataY(i);
	        if (x > min && x < max && y != 0) {
	            if (y > max1) {
	                max1 = y;
	                xMax = x;
	            }
	        }
	    }
	    return max1;
	}

	private double getMeanIDataSet(IDataSet data, double min, double max) {
	    int nsamples = 0;
	    double sum = 0;
	    double nEntries = 0;
	    for (int i = 0; i < data.getDataSize(0); i++) {
	        double x = data.getDataX(i);
	        double y = data.getDataY(i);
	        if (x > min && x < max && y != 0) {
	            nsamples++;
	            sum += x * y;
	            nEntries += y;
	        }
	    }
	    return sum / (double) nEntries;
	}

	private double getRMSIDataSet(IDataSet data, double min, double max) {
	    int nsamples = 0;
	    double mean = getMeanIDataSet(data, min, max);
	    double sum = 0;
	    double nEntries = 0;

	    for (int i = 0; i < data.getDataSize(0); i++) {
	        double x = data.getDataX(i);
	        double y = data.getDataY(i);
	        if (x > min && x < max && y != 0) {
	            nsamples++;
	            sum += Math.pow(x - mean, 2) * y;
	            nEntries += y;
	        }
	    }
	    return Math.sqrt(sum / (double) nEntries);
	}         

        
}


