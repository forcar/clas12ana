package org.clas.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;

import org.clas.viewer.DetectorMonitor;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.group.DataGroup;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;


public class EChv extends DetectorMonitor {
	
    DataGroup dg = new DataGroup(6,3);    
    IndexedList<List<Float>> mean = new IndexedList<List<Float>>(2);
    float[] ngain = {9.8f,8.8f,8.8f}; 
    static String gainPath = "/Users/colesmith/CLAS12ANA/ECcalib/ccdb/";
    static String  snpPath = "/Users/colesmith/clas12/HV/";
    public static final int[]     DHV = {20,150};
    public static final float[] HVEXP = {12,11};

    H1F h;
    int detid;
	
    public EChv(String name) {
    	super(name);
    }
    
    public List<String> getList(String path) { 
    	
    	List<String> out = new ArrayList();
    	String strCurrentLine;
       
    	try (BufferedReader br = new BufferedReader(new FileReader(path))) {
    	while ((strCurrentLine = br.readLine()) != null) out.add(strCurrentLine);
    	} catch (IOException e) {
    		e.printStackTrace();    		
    	}    	
    	return out;
    }
    
//  *** runGenerateSNP() *** 
    
    public static void runGenerateSNP(int run, int detid, String snpName) {
    	EChv reader = new EChv("EChv");
    	String gainFile = gainPath+"hvgain_"+run;
    	String  snpFile = snpPath+(detid==0?"PCAL_HV/":"ECAL_HV/")+snpName;
    	IndexedList<List<Float>> gain = parseGain(reader.getList(gainFile)); 
    	List<String>  snp = reader.getList(snpFile); 
    	writeSNP(run, detid, gain, snp);
    }
    
	public static void writeSNP(int run, int detid, IndexedList<List<Float>> gain, List<String> snp) {
		
		String line = new String();
		
		try { 
			File outputFile = new File(gainPath+(detid==0?"PCAL_HV_":"ECAL_HV_")+"hvgain_"+run+".snp");
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
	    	int n=0;
	    	for (String item : snp) {
	    		if(n<11) {outputBw.write(item); outputBw.newLine();}
	    		if(n>10) {
	    			String parts1[] = item.split(" ");
	    			String parts2[] = parts1[0].split(":");
	    			String parts3[] = parts2[0].split("_");
	    			float hvold = Float.parseFloat(parts1[2]);			
	    			
	    			if(parts2[1].equals("vset")) {
	    				int id=detid;
	    				int is=decodeSNP(parts3[4]);
	    				int il=decodeSNP(parts3[5]);
	    				int ip=Integer.parseInt(parts3[6].substring(1));
	    				
	    				float gv = gain.getItem(is,il,ip).get(0);
	    				float ge = gain.getItem(is,il,ip).get(1);
	    				float cn = gain.getItem(is,il,ip).get(2);
	    				float gerr = (gv>0) ? ge/gv:0;
	                    if (gv<0.5||gv>1.7||(gerr>0.25 && cn<6)||cn<2) gv=1f;
	                    double ratio=Math.pow(gv, 1./HVEXP[detid]);
	                    double hvnew = (ratio>0.5) ? hvold/ratio : hvold;
	                    if((hvnew-hvold)>100) hvnew=hvold;
	                    if (gerr>0.25 && cn>5) System.out.println(is+" "+il+" "+ip+" "+hvold+" "+cn+" "+gerr+" "+gv);
	    				line = parts1[0]+" "+parts1[1]+" "+hvnew;
	    				outputBw.write(line); outputBw.newLine();
	    		    }
	    		}
	    		n++;
	    	}

			outputBw.close();
			outputFw.close();
		}
		catch(IOException ex) {
			System.out.println("Error writing file '" );                   
			ex.printStackTrace();
		}

	} 
	
    public static int decodeSNP(String input) {
    	switch (input) {
    	case  "U": return 1;
    	case  "V": return 2;
    	case  "W": return 3;
    	case "UI": return 4;
    	case "VI": return 5;
    	case "WI": return 6;
    	case "UO": return 7;
    	case "VO": return 8;
    	case "WO": return 9;
    	case "SEC1": return 1;
    	case "SEC2": return 2;
    	case "SEC3": return 3;
    	case "SEC4": return 4;
    	case "SEC5": return 5;
    	case "SEC6": return 6;
    	}
    	return 0;
    }

    
    public static IndexedList<List<Float>> parseGain(List<String> list) {
    	IndexedList<List<Float>> gainList = new IndexedList<List<Float>>(3);
    	for (String item : list) {
			String parts[] = item.split(" ");
			int   is = Integer.parseInt(parts[0]);
			int   il = Integer.parseInt(parts[1]);
			int   ip = Integer.parseInt(parts[2]);
			float  g = Float.parseFloat(parts[3]);
			float ge = Float.parseFloat(parts[4]);
			float cn = Float.parseFloat(parts[5]);
			
    		if (!gainList.hasItem(is,il,ip)) gainList.add(new ArrayList<Float>(), is,il,ip);
    		gainList.getItem(is,il,ip).add(g); 
    		gainList.getItem(is,il,ip).add(ge); 
    		gainList.getItem(is,il,ip).add(cn); 
    	}
    	return gainList;
    } 
    
//  *** runCompareFiles() ***    
		
    public static void runCompareFiles(int id) {
    	EChv reader = new EChv("EChv");
		if(id==0) reader.compareFiles(0,"PCAL_HV-2017_11_23-09_15_28.snp","PCAL_HV-2022_06_21-21_55_52.snp"); 
        if(id==1) reader.compareFiles(1,"ECAL_HV-2017_11_23-09_15_28.snp","ECAL_HV-2022_06_27-20_38_08.snp"); 
        if(id==2) reader.compareFiles(2,"ECAL_HV-2017_11_23-09_15_28.snp","ECAL_HV-2022_06_27-20_38_08.snp"); 
		
    }
        
    public void compareFiles(int detid, String... item) {
		String dir = "/Users/colesmith/clas12/HV/";
		this.detid = detid; String det = detid==0?"PCAL":"ECAL";
		int n=0;
		for (String snp: item) {
			setGStyle(n);
			createHVHistos(detid,"(VOLTS)");  this.getDataGroup().add(dg,detid,0,0,0);
//			fillHistos(parseList(getList(dir+det+"_HV/"+getList(dir+"dump_"+det).get(it))),n);
			fillHistos(parseList(getList(dir+det+"_HV/"+snp)),n);
			getMean(n); n++;
		}
		if(item.length==2) setTitles();
		plotHistos();
    }
    
    public IndexedList<List<Integer>> parseList(List<String> list) {
    	
    	IndexedList<List<Integer>> hvlist = new IndexedList<List<Integer>>(4);
    	int n=0;
    	for (String item : list) {
    		if(n>10) {
    			String parts1[] = item.split(" ");
    			String parts2[] = parts1[0].split(":");
    			String parts3[] = parts2[0].split("_");
    			int val=(int)Float.parseFloat(parts1[2]);			
    			
    			if(parts2[1].equals("vset")) {
    				int id=detid;
    				int is=decode(parts3[4]);
    				int il=decode(parts3[5]);
    				int ip=Integer.parseInt(parts3[6].substring(1))-1;
    				if(is==4 && val>2500 && val<=2550) System.out.println(parts2[0]+" "+val);
    				int iil = (il>2)?il-3:il;
    				if(id<2 || (id>1&&il>2)) {
        			if (!hvlist.hasItem(id,is,iil,ip)) hvlist.add(new ArrayList<Integer>(),id,is,iil,ip);
                         hvlist.getItem(id,is,iil,ip).add(val);	
    				}
    		    }
    		}
    		n++;
    	}     	
    	return hvlist;
    }
    
    public void setGStyle(int icol) {
    	GStyle.getH1FAttributes().setOptStat("100");
    	GStyle.getAxisAttributesX().setTitleFontSize(18);    	
    	switch (icol) {
    	case 0: GStyle.getH1FAttributes().setFillColor(0);
    	        GStyle.getH1FAttributes().setLineWidth(5); break;
    	case 1: GStyle.getH1FAttributes().setFillColor(27);
                GStyle.getH1FAttributes().setLineWidth(2);
   	    }
    }
	
    public void createHVHistos(int id, String txt) {
    	
       for (int is=1; is<7; is++) {
    	   switch (id) {
    	   case 0:         
           h = new H1F("hv_pcal_u_"+is, 10, 750,1050);   
           h.setTitleX("Sector "+is+" PCAL U "+txt); 
           dg.addDataSet(h,is-1);  
           h = new H1F("hv_pcal_v_"+is, 10, 750, 1050); 
           h.setTitleX("Sector "+is+" PCAL V "+txt);        
           dg.addDataSet(h,is-1+6);            
           h = new H1F("hv_pcal_w_"+is, 10, 750, 1050);  
           h.setTitleX("Sector "+is+" PCAL W "+txt); 
           dg.addDataSet(h,is-1+12); 
           break;
    	   case 1:        
           h = new H1F("hv_ecin_u_"+is, 20, 1400, 3000);  
           h.setTitleX("Sector "+is+" ECIN U "+txt);    
           dg.addDataSet(h,is-1);  
           h = new H1F("hv_ecin_v_"+is, 20, 1400, 3000); 
           h.setTitleX("Sector "+is+" ECIN V "+txt);        
           dg.addDataSet(h,is-1+6);            
           h = new H1F("hv_ecin_w_"+is, 20, 1400, 3000);  
           h.setTitleX("Sector "+is+" ECIN W "+txt);  
           dg.addDataSet(h,is-1+12); 
           break;
    	   case 2:       
           h = new H1F("hv_ecou_u_"+is, 20, 1400, 3000);  
           h.setTitleX("Sector "+is+" ECOU U "+txt);     
           dg.addDataSet(h,is-1);  
           h = new H1F("hv_ecou_v_"+is, 20, 1400, 3000);  
           h.setTitleX("Sector "+is+" ECOU V "+txt);      
           dg.addDataSet(h,is-1+6);           
           h = new H1F("hv_ecou_v_"+is, 20, 1400, 3000); 
           h.setTitleX("Sector "+is+" ECOU W "+txt);  
           dg.addDataSet(h,is-1+12);   
    	   }
      }            
    } 
 
    public int decode(String input) {
    	switch (input) {
    	case  "U": return 0;
    	case  "V": return 1;
    	case  "W": return 2;
    	case "UI": return 0;
    	case "VI": return 1;
    	case "WI": return 2;
    	case "UO": return 3;
    	case "VO": return 4;
    	case "WO": return 5;
    	case "SEC1": return 0;
    	case "SEC2": return 1;
    	case "SEC3": return 2;
    	case "SEC4": return 3;
    	case "SEC5": return 4;
    	case "SEC6": return 5;
    	}
    	return 0;
    }
    
    public boolean reject(int idet, int is, int il, int ip) {
    	if(idet==0) return false;
    	if(is==5 && il==1 && ip==15) return true;
    	if(is==5 && il==1 && ip==16) return true;
    	if(is==6 && il==2 && ip==13) return true;
    	if(is==6 && il==0 && ip==4)  return true;
    	if(is==6 && il==0 && ip==6)  return true;
    	return false;
    }
    
    public void fillHistos(IndexedList<List<Integer>> hvlist, int n) {
    	IndexGenerator ig = new IndexGenerator();
	    for (Map.Entry<Long,List<Integer>>  entry : hvlist.getMap().entrySet()){
	           int id = ig.getIndex(entry.getKey(), 0);   
	           int is = ig.getIndex(entry.getKey(), 1);   
	           int il = ig.getIndex(entry.getKey(), 2);   
	           int ip = ig.getIndex(entry.getKey(), 3);  
	           int hv = entry.getValue().get(0);
	           if(!reject(id,is,il,ip)) ((H1F) this.getDataGroup().getItem(id,0,0,0).getData(is+il*6).get(n)).fill(hv);	           
	    }    	
    }
    
    public void plotHistos() {
        JFrame frame = new JFrame("ECAL HV");
        frame.setSize(2000,800);
        EmbeddedCanvas canvas = new EmbeddedCanvas();
        drawGroup(canvas,getDataGroup().getItem(detid,0,0,0));	       	
        frame.add(canvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);        
    }
    
    public void getMean(int n) {   	
    	for (int il=0; il<3; il++) {
    	    for (int is=0; is<6; is++) {
    	    	float val = (float) ((H1F) this.getDataGroup().getItem(detid,0,0,0).getData(is+il*6).get(n)).getMean();  			 	
    	    	if(!mean.hasItem(is,il)) mean.add(new ArrayList<Float>(),is,il);
                    mean.getItem(is,il).add(val); 
    		}
    	}
    }
    
    public void setTitles() {   	
    	for (int il=0; il<3; il++) {
    	    for (int is=0; is<6; is++) {    	    	
    	    	float v1 = mean.getItem(is,il).get(0); float v2 = mean.getItem(is,il).get(1);
    	    	float delV = v2-v1;  float delG = ngain[detid]*delV/v1;
    	    	String tit = "#Delta V = "+String.format("%.0f",delV)+"   #Delta G/G = "+String.format("%.2f",delG);
    	    	((H1F) this.getDataGroup().getItem(detid,0,0,0).getData(is+il*6).get(0)).getAttributes().setTitle(tit);  	    	
    	    }
    	}
    }

	public static void main(String[] args) {
		runCompareFiles(2);	
// RGD		
//		runGenerateSNP(18312,0,"PCAL_HV-2023_09_23-09_56_53.snp");  
//		runGenerateSNP(18312,1,"ECAL_HV-2023_09_25-07_29_35.snp");  
// RGC		
//		runGenerateSNP(16066,0,"PCAL_HV-2022_06_17-18_24_13.snp");
//		runGenerateSNP(16066,1,"ECAL_HV-2022_06_17-18_23_29.snp");
//		runGenerateSNP(16149,1,"ECAL_HV_hvgain_16066.snp");		
	}

}

