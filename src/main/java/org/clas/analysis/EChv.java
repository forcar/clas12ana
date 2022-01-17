package org.clas.analysis;

import java.io.BufferedReader;
import java.io.FileReader;
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

	H1F h;
	int detid;
	
	public EChv(String name) {
        super(name);		
	}
	
	public void setGStyle(int icol) {
    	GStyle.getH1FAttributes().setOptStat("1110");
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
           h = new H1F("hv_pcal_u_"+is, 10, 750,1000);   
           h.setTitleX("Sector "+is+" PCAL U "+txt); 
           dg.addDataSet(h,is-1);  
           h = new H1F("hv_pcal_v_"+is, 10, 750, 1000); 
           h.setTitleX("Sector "+is+" PCAL V "+txt);        
           dg.addDataSet(h,is-1+6);            
           h = new H1F("hv_pcal_w_"+is, 10, 750, 1000);  
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
    
    public void fillHistos(IndexedList<List<Integer>> hvlist, int n) {
    	IndexGenerator ig = new IndexGenerator();
	    for (Map.Entry<Long,List<Integer>>  entry : hvlist.getMap().entrySet()){
	           int id = ig.getIndex(entry.getKey(), 0);   
	           int is = ig.getIndex(entry.getKey(), 1);   
	           int il = ig.getIndex(entry.getKey(), 2);   
	           int ip = ig.getIndex(entry.getKey(), 3);  
	           int hv = entry.getValue().get(0);
	           if(hv>100) ((H1F) this.getDataGroup().getItem(id,0,0,0).getData(is+il*6).get(n)).fill(hv);	           
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
    
    public void processFile(int detid, int... item) {
		String dir = "/Users/colesmith/clas12/HV/";
		this.detid = detid; String det = detid==0?"PCAL":"ECAL";
		int n=0;
		for (int it: item) {
			setGStyle(n);
			createHVHistos(detid,"(VOLTS)");  this.getDataGroup().add(dg,detid,0,0,0);
			fillHistos(parseList(getList(dir+det+"_HV/"+getList(dir+"dump_"+det).get(it))),n);
			n++;
		}
		plotHistos();
    }
    
	public static void main(String[] args) {		
		EChv reader = new EChv("EChv");
		//0,6,0:PCAL_HV-2017_11_23-09_16_01.snp PCAL_HV-2020_10_09-10_43_25.snp
		//1,6,0:ECAL_HV-2020_10_09-10_43_57.snp ECAL_HV-2017_11_23-09_15_28.snp
		//2,6,0:ECAL_HV-2020_10_09-10_43_57.snp ECAL_HV-2017_11_23-09_15_28.snp
		int det=0, snp=0;
		if(args.length==0) reader.processFile(det,snp);
		if(args.length==1) reader.processFile(Integer.parseInt(args[0]),snp);
		if(args.length==2) reader.processFile(Integer.parseInt(args[0]),Integer.parseInt(args[1]));
		if(args.length==3) reader.processFile(Integer.parseInt(args[0]),Integer.parseInt(args[1]),Integer.parseInt(args[2]));	
	}

}

