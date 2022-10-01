package org.clas.tools;

import javax.swing.JFrame;

import org.clas.analysis.ECt;
import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.IDataSet;
//import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;

public class FitData {
	
	public GraphErrors graph = null;
	public H1F          hist = null;

	public double pmin;
	public double pmax;
	public double fmin;
	public double fmax;
	public double amp;
	public double mean;
	public double meane;
	public double sigma;
	public double sigmae;
	public double p0,p1,p2,p3,p4,p0e,p1e,p2e,p3e,p4e;
	public double sig1=2.5;
	public double sig2=1.7;
	public int func;
	public int integral;
	public int intmin=30;
	public int fitcol=4;
	public String f_optstat="1110";
	public String g_optstat="1100";
	public String h_optstat="100";
	public Boolean doFit = false;
	
	String predefFunctionsF1D[] = {"[amp]*gaus(x,[mean],[sigma])",
			                       "[amp]*gaus(x,[mean],[sigma])+[p0]",
			                       "[amp]*gaus(x,[mean],[sigma])+[p0]+x*[p1]",
			                       "[amp]*gaus(x,[mean],[sigma])+[p0]+x*[p1]+x*x*[p2]",
			                       "[amp]*gaus(x,[mean],[sigma])+[p0]+x*[p1]+x*x*[p2]",
			                       "[p0]", "[p0]+[p1]*x", "[p0]+[p1]*x+[p2]*x*x","[p0]+[p1]*x+[p2]*x*x+[p3]*x*x*x",
			                       "[a]*exp(x,[b])",
			                       "[a]+[b]*cos(x*[c])",
			                       "[a]+[b]*cos(x*[d])+[c]*cos(2*x*[e])",
			                       "1/((1-[p])+[p]/x)",
                                   "[p0]+exp(-(x-[p1])/[p2])+1-exp(-([p3]-x)/[p4])",
                                   "[sf1]*(1+[sf3]/x+[sf4]/x/x)",
                                   "[sf1]*(1+[sf3]*0.1/x+[sf4]*0.001/x/x)",
                                   "[p0]+[p1]/x +[p2]/x^0.5"};
	
	public FitData(GraphErrors graph) {
	    setGraph(graph);
	}

	public void setGraph(GraphErrors graph) {
	    this.graph = graph;
	    this.graph.getAttributes().setOptStat(g_optstat);
	}
	
	public void setHist(H1F hist) {
		this.hist = hist;
		this.hist.setLineWidth(1);
		this.hist.setOptStat(h_optstat);
		this.hist.setTitle("");
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
	
	public void setSigmas(double sig1, double sig2) {
		this.sig1 = sig1;
		this.sig2 = sig2;
	}
	
	public void simpleFit(double pmin, double pmax, double fmin, double fmax) {
	    hist.setFunction(new F1D("f",predefFunctionsF1D[0], fmin, fmax)); 
	    hist.getFunction().setLineWidth(1);
	    amp   = getMaxYIDataSet(graph,pmin,pmax,true);
	    mean  = getMaxYIDataSet(graph,pmin,pmax,false);
	    sigma = getRMSIDataSet(graph,pmin,pmax);
	    hist.getFunction().setParameter(0, 30); hist.getFunction().setParameter(1, 1); hist.getFunction().setParameter(2, 0.1);
//	    hist.getFunction().setRange(mean-2*sigma, mean+2*sigma);
	    DataFitter.fit(hist.getFunction(), hist, "Q");
	}

	public void initFit(int func, double pmin, double pmax, double fmin, double fmax) {
		this.func = func;
	    this.pmin = pmin; this.fmin=fmin;
	    this.pmax = pmax; this.fmax=fmax;	    
	    graph.setFunction(new F1D("f",predefFunctionsF1D[func], fmin, fmax)); 
	    graph.getFunction().setLineWidth(1); 
	    if(func<5) {
	      amp   = getMaxYIDataSet(graph,pmin,pmax,true);
	      mean  = getMaxYIDataSet(graph,pmin,pmax,false);
//	      mean  = getMeanIDataSet(graph,pmin,pmax); 
	      sigma = getRMSIDataSet(graph,pmin,pmax);
	      initFunc(0,amp);
	      initFunc(1,mean);
	      initFunc(2,sigma);
		  if(func==3) {initFunc(3,-5); initFunc(4,0.3); initFunc(5, -1e-3);}	    
		  if(func==0) graph.getFunction().setRange(mean-sig1*sigma, mean+sig2*sigma);
		  if(func==3) graph.getFunction().setRange(fmin,fmax);
	    }
	    if (func==6)  {initFunc(0,20.0); initFunc(1,0.057) ;                  graph.getFunction().setRange(fmin, fmax);g_optstat="1110";}
	    if (func==7)  {initFunc(0,0.23); initFunc(1,0.56) ; initFunc(2,-0.3); graph.getFunction().setRange(fmin, fmax);}
	    if (func==14) {initFunc(0,0.25); initFunc(1,-0.018,-0.040,-0.016); initFunc(2,0.0006,0.0005,0.0007); graph.getFunction().setRange(fmin, fmax);}
	    if (func==15) {initFunc(0,0.27); initFunc(1,-0.146);initFunc(2,0.117);graph.getFunction().setRange(fmin, fmax);g_optstat="";f_optstat="1110";}
	    if (func==16) {initFunc(0,0.5);  initFunc(1,0.001); initFunc(2,100)  ;graph.getFunction().setRange(fmin, fmax);g_optstat="";f_optstat="1110";}
	    if (func==13) {graph.getFunction().setRange(fmin, fmax);g_optstat="";f_optstat="111110";}
	}
	
	public void initFunc(int par, double val) {
        graph.getFunction().setParameter(par, val);
	}
	
	public void initFunc(int par, double val, double min, double max) {
        graph.getFunction().setParameter(par, val);
        graph.getFunction().setParLimits(par, min, max);
	}

	public void fitGraph(String opt, Boolean fitEnable, Boolean fitVerbose) {
	    if (doFit&&fitEnable) DataFitter.fit(graph.getFunction(), graph, "Q");
	    if (func>4) {
	    	for (int i=0; i<graph.getFunction().getNPars(); i++) {
	    		if(i==0) {p0 = graph.getFunction().parameter(0).value(); p0e  = graph.getFunction().parameter(0).error();}
	    		if(i==1) {p1 = graph.getFunction().parameter(1).value(); p1e  = graph.getFunction().parameter(1).error();}
	    		if(i==2) {p2 = graph.getFunction().parameter(2).value(); p2e  = graph.getFunction().parameter(2).error();}
	    		if(i==3) {p3 = graph.getFunction().parameter(3).value(); p3e  = graph.getFunction().parameter(3).error();}
	    		if(i==4) {p4 = graph.getFunction().parameter(4).value(); p4e  = graph.getFunction().parameter(4).error();}
	    	}
	    }	                  
	    if (func<5) {
	       amp   = graph.getFunction().getParameter(0);
	       mean  = graph.getFunction().parameter(1).value();
	       meane = graph.getFunction().parameter(1).error();
	       sigma = graph.getFunction().parameter(2).value();   
	      sigmae = graph.getFunction().parameter(2).error();  
	    }
	    graph.getFunction().setLineColor(fitcol);
	    graph.getAttributes().setOptStat(g_optstat);
	    graph.getFunction().setOptStat(opt=="Q"?"0":f_optstat);
//	    graph.getFunction().setOptStat("1110");
	}

	public void plotGraph(EmbeddedCanvas c, int col) {
	    graph.getFunction().setLineColor(col); 
	    c.draw(graph);
	}

	public GraphErrors getGraph() {
	    return graph;
	}

	private double getMaxYIDataSet(IDataSet data, double min, double max, boolean opt) {
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
	    return opt?max1:xMax;
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
	            sum += (x-mean)*(x-mean)* y;
	            nEntries += y;
	        }
	    }
	    return Math.sqrt(sum / (double) nEntries);
	}

    
    public static void fitTest(int is, int il, int iv, int nmax) {
    	int fnum=13, fmin=10;
    	GraphErrors g;
    	DetectorMonitor mon = new DetectorMonitor("ECt");
    	mon.cfitEnable = true; mon.fitVerbose = true;
        JFrame frame = new JFrame("fitTest");
        frame.setSize(2100,1700);
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        if (il==0) canvas.divide(9,8);
        if (il>0)  canvas.divide(6,6);
        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
        for (int i=0; i<nmax; i++) {
            String nam ="s"+is+"l"+il+"v"+iv+"p"+i;
        	g = mon.getGraph("/Users/colesmith/CLAS12ANA/ECt/twplots/"+nam+".vec");
        	canvas.cd(i); canvas.getPad(i).getAxisY().setRange(-3,6);
        	if(g.getDataSize(0)==0) continue;
        	FitData fd = new FitData(g);
            fd.initFit(fnum,0,1,fmin,g.getDataX(g.getDataSize(0)-1)*1.05);
            switch (il) {
            case 0: fd.initFunc(0,-0.5); fd.initFunc(1,20,18,22);  fd.initFunc(2,9,7,11);   fd.initFunc(3,170,160,180); fd.initFunc(4,20,15,25); break;
            case 1: fd.initFunc(0,-0.5); fd.initFunc(1,13.9,5,22); fd.initFunc(2,5.4,5,11); fd.initFunc(3,125,100,180); fd.initFunc(4,15,14,25); break;
            case 2: fd.initFunc(0,-0.5); fd.initFunc(1,17,10,18);  fd.initFunc(2,5,3,11);   fd.initFunc(3,145,110,180); fd.initFunc(4,15,14,20);
            }
            fd.doFit = true; 
            fd.fitGraph("",true,true);fd.getGraph().setTitle(nam);
            canvas.draw(fd.getGraph());
//        	canvas.draw(mon.fitEngine(g,fnum,fmin).getGraph());
        }
    }
	
    public static void main(String[] args) {
    	fitTest(6,0,1,36);   	
    }
        
}


